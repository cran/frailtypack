module natural_effects
  use Autres_fonctions!,only:bgos,Cholesky_Factorisation,splinebasisIndiv,bb
  use tailles
  use var_mediation
  implicit none
  !private
  double precision,dimension(:),allocatable,save,private::b_ind,bgamma
  !Beware of out-of-scope modification of variable and multithreading !
  !One thread may set var = 1, do some computation and then use var !
  !But in th mean time another thread changed the value of var !
  !integer,save,private::zs,zt
  !integer,save::posrec !records position for writting in unformatted binary file
  !integer,private,save::indic_int
  !double precision,save,private::xs,xt
  !double precision,private,save::time_ev
  double precision,save,private::bs,bt
  integer,private::integ_type
  double precision,dimension(:,:),allocatable,save,private::vectgamma
  double precision::sig_s,sig_t,sig_st,sig_omega,sig_u,eta_ind,alpha_ind
  double precision,dimension(:),allocatable,private::vls2,vss2,vlt2,vst2,tS2,tT2
  double precision,dimension(30),parameter,private::lag_nodes=(/0.04740718d0,&
                0.24992392d0,0.61483345d0,1.14319583d0,1.83645455d0,&
                2.69652187d0,3.72581451d0,4.92729377d0,6.30451559d0,&
                7.86169329d0,9.60377599d0,11.53654660d0,13.66674469d0,&
                16.00222119d0,18.55213484d0,21.32720432d0,24.34003576d0,&
                27.60555480d0,31.14158670d0,34.96965201d0,39.11608495d0,&
                43.61365291d0,48.50398616d0,53.84138541d0,59.69912186d0,&
                66.18061779d0,73.44123860d0,81.73681051d0,91.55646652d0,&
                104.15752443d0/)
  double precision,dimension(30),parameter,private::lag_weights=(/1.160441D-01,&
                2.208511D-01,2.413998D-01,1.946368D-01,1.237284D-01,&
                6.367878D-02,2.686048D-02,9.338071D-03,2.680697D-03,&
                6.351291D-04,1.239075D-04,1.982879D-05,2.589351D-06,&
                2.740943D-07,2.332831D-08,1.580746D-09,8.427479D-11,&
                3.485161D-12,1.099018D-13,2.588313D-15,4.437838D-17,&
                5.365918D-19,4.393947D-21,2.311410D-23,7.274588D-26,&
                1.239150D-28,9.832375D-32,2.842324D-35,1.878608D-39,&
                8.745980D-45/)
  contains
  subroutine compute_rt(b,covar,knots,np,nparams,zi,times,res,boot,integ)
      !!$ use omp_lib
      implicit none
      integer,dimension(10)::nparams
      integer::np
      integer::nsplines2
      integer::boot
      double precision,dimension(ntimes)::times
      double precision,dimension(ntimes,4,1+boot*nboot),intent(out)::res
      double precision,dimension(2+nparams(7))::knots
      integer::i,j
      integer::integ
      double precision::xx1,xx2,xx3
      double precision,dimension(np)::b
      double precision,dimension(np,np)::covar
      double precision,dimension(nparams(1)+4)::zi
      !integer,dimension(4)::mtaille
      !double precision,dimension(mtaille(1)),intent(in)::lamOut
      !double precision,dimension(mtaille(3)),intent(in)::suOut
      !double precision,dimension(mtaille(2)),intent(in)::lam2Out
      !double precision,dimension(mtaille(4)),intent(in)::su2Out
      !double precision,dimension(mtaille(1)),intent(in)::x1Out
      !double precision,dimension(mtaille(2)),intent(in)::x2Out
      !nparams = numbers of parameters related to
      !nparams(1) = nz1+2
      !nparams(2) = nz2 +2
      !nparams(3) = nva
      !nparams(4) = nombre frailty (4 : var omega + 3 param signa, or 3 : var omega + diagonal sigma (no covariance))
      !nparams(5) = indicator zeta
      !nparams(6) = indicator alpha
      !nparams(7) = parameter gamma(s) = nknots (interior knots)
      !nparams(8) = parameter gamma(s) = splines_ord
      !nparams(9) = mtaille(1)
      !nparams(10)= mtailles(2)
      !np = sum(nparams)

      !!!!
      nsplines2 = nparams(7)+nparams(8)
      allocate(vls2(nparams(9)),vss2(nparams(9)),vlt2(nparams(10)),vst2(nparams(10)))
      allocate(tS2(nparams(9)),tT2(nparams(9)))
      allocate(bgamma(nsplines2))
      allocate(b_ind(np))
      allocate(basissplines(nparams(7)+nparams(8)))
      allocate(vectgamma(1001,2))

      integ_type=integ
      

      call random_params(np,nparams,b,covar,zi,knots,.FALSE.)
      !write(*,*) "Computation for R(t)."
      do i=1,ntimes
        call montecarlo(times(i),1,1,nmc,xx1)
        call montecarlo(times(i),0,1,nmc,xx2)
        call montecarlo(times(i),0,0,nmc,xx3)
        res(i,2,1) = xx1
        res(i,3,1) = xx2
        res(i,4,1) = xx3
      enddo
      !write(*,*) "End of computation for R(t)."
      if(boot.eq.1) then
        !write(*,*) "End of computation for R(t)."
        !write(*,*) "Now computing boostraped samples"
        do i=1,nboot
          call random_params(np,nparams,b,covar,zi,knots,.TRUE.)
          res(:,1,i+1) = times
          do j=1,ntimes
            call montecarlo(times(j),1,1,nmc,xx1)
            call montecarlo(times(j),0,1,nmc,xx2)
            call montecarlo(times(j),0,0,nmc,xx3)
            res(j,2,i+1) = xx1
            res(j,3,i+1) = xx2
            res(j,4,i+1) = xx3
          enddo
        end do
        !write(*,*) "End of boostrap computing"
      end if
      deallocate(basissplines,vectgamma,b_ind,bgamma)
      deallocate(tS2,tT2,vls2,vss2,vlt2,vst2)
      return
  end subroutine
  subroutine montecarlo(t,zs,zt,nsim,result)
       !$ use OMP_LIB
       implicit none
       double precision::t
       integer::zs,zt
       integer::nsim
       integer::i
       double precision,dimension(2,2)::sigma
       double precision,external::h
       double precision,dimension(:,:),allocatable::x1
       double precision::result
       double precision,dimension(nsim)::res

       sigma(1,1)=sig_omega+sig_u+sig_s+zs*sig_s
       sigma(2,2)=(eta_ind**2.d0)*sig_omega+(alpha_ind**2.d0)*sig_u+sig_t
       sigma(1,2)=eta_ind*sig_omega+alpha_ind*sig_u+zs*zt*sig_st
       sigma(2,1)=sigma(1,2)

       call Cholesky_Factorisation2(sigma,2)
       allocate(x1(nsim,2))
       do i=1,nsim
         call bgos(1.0d0,2,x1(i,1),x1(i,2),0.d0)
       enddo
       x1=matmul(x1,transpose(sigma))
       result=0.0d0
       !$omp parallel do default(none) private(i) shared(nsim,x1,zs,zt,t,res)
       do i=1,nsim
         call one_dim(t,x1(i,1),x1(i,2),zs,zt,res(i))
       enddo
       !$omp end parallel do
       result=sum(res)/nsim
       deallocate(x1)
       return
  end subroutine
  subroutine one_dim(time,xs,xt,zs,zt,res)
    implicit none
    double precision::r1,r2
    integer::i,zs,zt
    double precision::time,xs,xt,s,res
    r1=0.d0
    r2=0.d0
    if(integ_type.eq.1)then
      do i=0,199
        s = time*(2*i+1)/400.d0
        r1 = r1 + (time/200.0d0) * survival_t(time,s,xs,xt,zs,zt,1)
      end do
    else if (integ_type.eq.2)then
      do i=1,30
        r1 = r1 + survival_t(time,lag_nodes(i),xs,xt,zs,zt,1)*&
                  dexp(lag_nodes(i))*lag_weights(i)
      end do
    else
      !call qagi(survival_t, time_ev, inf, epsabs, epsrel, r1, abserr, neval, ier)
    end if
    r2 = survival_t(time,0.0d0,xs,xt,zs,zt,2)
    res=r1 + r2
    return
  end subroutine
  double precision function survival_t(time,s,xs,xt,zs,zt,indic_int)
    IMPLICIT NONE
    double precision::time,s,xs,xt
    double precision::preds,predt
    double precision::h0tt,h0ts,h0ss,l0ss,h0st
    double precision::gamma
    integer::zs,zt
    integer::i,j,indic_int 
    survival_t = 0.d0 
    if(s.lt.0.0d0) then
      survival_t = 0.d0
      goto 180
    endif

    preds=xs+bs*zs
    predt=xt+bt*zt

    !h0ss & l0ss
    if(s.le.0.d0) then
      h0ss = 0.d0
      l0ss = 0.d0
    else
      if(s.le.tS2(1)) then
        h0ss = -dlog(vss2(1))
        l0ss = vls2(1)
      else if (s.ge.tS2(size(tS2))) then
        h0ss = -dlog(vss2(size(tS2)))
        l0ss = vls2(size(tS2))
      else
        i=0
        do while (tS2(i+1).lt.s)
          i = i+1
        end do
        h0ss = -dlog(vss2(i))
        l0ss = vls2(i)
      end if
    end if
    !h0tt
    if(time.le.0.d0) then
      h0tt = 0.0d0
    else
      if(time.le.tT2(1)) then
        h0tt = -dlog(vst2(1))
      else if (time.ge.tT2(size(tT2))) then
        h0tt = -dlog(vst2(size(tT2)))
      else
        i=0
        do while (tT2(i+1).lt.time)
          i = i+1
        end do
        h0tt = -dlog(vst2(i))
      end if
    end if

    if(indic_int.eq.1) then
      !call calcul_risque(nz1,nz2,b_ind,s,dummy,h0ts,2)
      !h0ts = -dlog(h0ts)
      if(s.gt.time)then
        survival_t = 0.0d0
        goto 180
      end if
      !h0ts
      if(s.le.0.d0)then
        h0ts = 0.0d0
      else
        if(s.le.tT2(1)) then
          h0ts = -dlog(vst2(1))
        else if (s.ge.tT2(size(tT2))) then
          h0ts = -dlog(vst2(size(tT2)))
        else
          i=0
          do while (tT2(i+1).lt.s)
            i = i+1
          end do
          h0ts = -dlog(vst2(i))
        end if
      end if
      if((s.lt.boundaryknotsurro(1)).or.(s.gt.boundaryknotsurro(2)))then
        gamma  = 0.0d0
      else
        j=1
        do while(vectgamma(j,1).lt.s)
          j = j + 1
        end do
        gamma = vectgamma(j,2)
      end if
      survival_t=exp(-(h0ts*dexp(predt)+dexp(predt+gamma)*(h0tt-h0ts)))*&
                  l0ss*dexp(preds)*dexp(-h0ss*dexp(preds))
    endif
    if(indic_int.eq.2)then
      !call calcul_risque(nz1,nz2,b_ind,time_ev,dummy,h0st,1)
      !h0st
      if(time.le.0)then
        h0st = 0.0d0
      else
        if(time.le.tS2(1)) then
          h0st = -dlog(vss2(1))
        else if (time.ge.tS2(size(tS2))) then
          h0st = -dlog(vss2(size(tS2)))
        else
          i=0
          do while (tS2(i+1).lt.time)
            i = i+1
          end do
          h0st = -dlog(vss2(i))
        end if
      end if
      survival_t=exp(-h0tt*exp(predt))*exp(-h0st*exp(preds))
    endif
    180 continue
    return
  end function
  subroutine random_params(np,nparams,b,covar,zi,knots,ran)
    implicit none
    integer::np
    integer,dimension(10)::nparams
    double precision,dimension(np),intent(in)::b
    double precision,dimension(np,np),intent(in)::covar
    double precision,dimension(2+nparams(7))::knots
    double precision,dimension(np,np)::cdummy
    !double precision,dimension(np,np)::c2
    double precision,dimension(np)::x1
    double precision::dummy,s
    double precision,dimension(2,2)::temp_mat
    double precision,dimension(nparams(1)+4)::zi
    logical::ran
    integer::i,j,nsplines2
    nsplines2 = nparams(7)+nparams(8)
    if(.NOT.ran)then
      b_ind = b
      bgamma=b_ind((np-nsplines2+1):np)
      vectgamma = 0.0d0
      basissplines=0.0d0
      bs=b_ind(np-nsplines2-1)
      bt=b_ind(np-nsplines2)
      if(nparams(5).eq.1) then
        eta_ind = b_ind(nparams(1)+nparams(2)+nparams(5))
      else
        eta_ind = 1.0d0
      end if
      if(nparams(6).eq.1)then
        alpha_ind = b_ind(nparams(1)+nparams(2)+nparams(5)+5)
      else
        alpha_ind=1.0d0
      end if
      sig_omega = b_ind(nparams(1)+nparams(2)+nparams(5)+1)**2
      sig_s =b_ind(nparams(1)+nparams(2)+nparams(5)+2)
      sig_t  =b_ind(nparams(1)+nparams(2)+nparams(5)+3)
      sig_st =b_ind(nparams(1)+nparams(2)+nparams(5)+4)

      sig_u = b_ind(nparams(1)+nparams(2)+nparams(5)+nparams(6)+5)
      temp_mat(1,:) = (/sig_s,0.d0/)
      temp_mat(2,:) = (/sig_st, sig_t/)
      temp_mat = matmul(temp_mat,transpose(temp_mat))

      sig_s = temp_mat(1,1)
      sig_t = temp_mat(2,2)
      sig_st =temp_mat(1,2)
      do i=0,1000
        s = (i/1000.d0)*(knots(2)-knots(1)) + knots(1)
        call splinebasisIndiv(splines_ord-1,nknots+2*splines_ord,nknots,nknots&
                              +splines_ord,s,knots(3:(2+nparams(7))),knots(1:2),&
                               basissplines)
                               vectgamma(i+1,1) = s
        do j=1,(nknots+splines_ord)
          vectgamma(i+1,2) = vectgamma(i+1,2) + bgamma(j)*basissplines(j)
        end do
      end do
      call baseline_rs(nparams(1)-2,nparams(2)-2,b_ind,nparams(9),&
                       nparams(10),zi,tS2,vls2,vss2,tT2,vlt2,vst2)
    end if
    if(ran) then
      cdummy=covar
      call Cholesky_Factorisation2(cdummy,np)
      do i=1,np
        call bgos(1.0d0,2,x1(i),dummy,0.d0)
      enddo
      x1=matmul(x1,transpose(cdummy))
      do i=1,np
        b_ind(i) = b(i) + x1(i)
      end do
      bgamma=b_ind((np-nsplines2+1):np)
      vectgamma = 0.0d0
      basissplines=0.0d0
      bs=b_ind(np-nsplines2-1)
      bt=b_ind(np-nsplines2)
      if(nparams(5).eq.1) then
        eta_ind = b_ind(nparams(1)+nparams(2)+nparams(5))
      else
        eta_ind = 1.0d0
      end if
      if(nparams(6).eq.1)then
        alpha_ind = b_ind(nparams(1)+nparams(2)+nparams(5)+5)
      else
        alpha_ind=1.0d0
      end if
      sig_omega = b_ind(nparams(1)+nparams(2)+nparams(5)+1)**2
      sig_s =b_ind(nparams(1)+nparams(2)+nparams(5)+2)
      sig_t  =b_ind(nparams(1)+nparams(2)+nparams(5)+3)
      sig_st =b_ind(nparams(1)+nparams(2)+nparams(5)+4)

      sig_u = b_ind(nparams(1)+nparams(2)+nparams(5)+nparams(6)+5)
      temp_mat(1,:) = (/sig_s,0.d0/)
      temp_mat(2,:) = (/sig_st, sig_t/)
      temp_mat = matmul(temp_mat,transpose(temp_mat))
      sig_s = temp_mat(1,1)
      sig_t = temp_mat(2,2)
      sig_st =temp_mat(1,2)
      !!
      do i=0,1000
        s = (i/1000.0d0)*(knots(2)-knots(1)) + knots(1)
        call splinebasisIndiv(splines_ord-1,nknots+2*splines_ord,nknots,nknots&
                              +splines_ord,s,knots(3:(2+nparams(7))),knots(1:2),&
                               basissplines)
        vectgamma(i+1,1) = s
        do j=1,(nknots+splines_ord)
          vectgamma(i+1,2) = vectgamma(i+1,2) + bgamma(j)*basissplines(j)
        end do
      end do
      call baseline_rs(nparams(1)-2,nparams(2)-2,b_ind,nparams(9),&
                       nparams(10),zi,tS2,vls2,vss2,tT2,vlt2,vst2)
    end if
  end subroutine
  subroutine baseline_rs(nz1,nz2,b,mt1,mt2,zi,tS3,vls3,vss3,tT3,vlt3,vst3)
    Implicit none
    integer,intent(in):: nz1,nz2,mt1,mt2
    double precision ,dimension(npmax),intent(in):: b
    integer::i,j,n,k,l
    double precision::x1,x2,h,su,lam
    double precision ,dimension(-2:nz1)::the2
    double precision ,dimension(-2:nz2)::the1T
    double precision,dimension(-2:(nz1+3))::zi
    double precision,dimension(mt1)::tS3
    double precision,dimension(mt2)::tT3
    double precision,dimension(mt1)::vls3,vss3
    double precision,dimension(mt2)::vlt3,vst3
    n  = nz1+2
    k = 0
    do i=1,(nz1+2)
        k = k + 1
        l = 0
        do j=1,(nz1+2)
            l = l + 1
        end do
    end do
    k = 0
    do i=(nz1+2)+1,(nz1+2)+nz2+2
        k = k + 1
        l = 0
        do j=(nz1+2)+1,(nz1+2)+nz2+2
            l = l + 1
        end do
    end do
    do i=1,nz1+2
        the1T(i-3)=b(i)**2
    end do
    do i=1,nz2+2
        j = nz1+2+i
        the2(i-3)=b(j)**2
    end do
    h = (zi(n)-zi(1))*0.01d0
    x1 = zi(1)
        do i=1,mt1
            if(i .ne.1)then
                x1 = x1 + h
            end if
            call calcul_rs(x1,the1T,nz1+2,zi,su,lam)
            tS3(i)=x1
            vls3(i)=lam
            vss3(i)=su
        end do
    x2 = zi(1)
    do i=1,mt2
        if(i.ne.1)then
            x2 = x2 + h
        endif
        call calcul_rs(x2,the2,nz2+2,zi,su,lam)
        tT3(i)=x2
        vlt3(i)=lam
        vst3(i)=su
    end do
    return
  end subroutine baseline_rs
  subroutine calcul_rs(x,the,n,zi,su,lam)
    IMPLICIT NONE
    integer,intent(in)::n
    double precision,intent(in)::x
    double precision,intent(out)::lam,su
    double precision,dimension(-2:(n-2)),intent(in)::the
    double precision,dimension(-2:(n+1))::zi
    integer::j,k,i
    double precision::ht,ht2,h2,som,htm,h2t,h3,h2n,hn, &
    im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2, &
    h,gl,hh
    j=0
    gl=0.d0
    som = 0.d0
    do k = 2,n-1
        if ((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
            j = k-1
            if (j.gt.1)then
                do i=2,j
                som = som+the(i-4)
                end do
            endif
            ht = x-zi(j)
            htm= x-zi(j-1)
            h2t= x-zi(j+2)
            ht2 = zi(j+1)-x
            ht3 = zi(j+3)-x
            hht = x-zi(j-2)
            h = zi(j+1)-zi(j)
            hh= zi(j+1)-zi(j-1)
            h2= zi(j+2)-zi(j)
            h3= zi(j+3)-zi(j)
            h4= zi(j+4)-zi(j)
            h3m= zi(j+3)-zi(j-1)
            h2n=zi(j+2)-zi(j-1)
            hn= zi(j+1)-zi(j-2)
            hh3 = zi(j+1)-zi(j-3)
            hh2 = zi(j+2)-zi(j-2)
            mm3 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
            mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
            *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
            h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3 = (0.25d0*(x-zi(j-3))*mm3)+(0.25d0*hh2*mm2) &
            +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
            im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
            im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
            im  = ht*mm*0.25d0
            gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)+(the(j)*im)
            lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)+(the(j)*mm)
        endif
    end do
    if(x.ge.zi(n))then
        som = 0.d0
        do i=1,n
            som = som+the(i-3)
        end do
        gl = som
    endif
    su  = dexp(-gl)
    return
  end subroutine calcul_rs
end module
module natural_effects_longi 
  implicit none 
  integer::integ_type
  double precision,dimension(30),parameter,private::lag_nodes=(/0.04740718d0,&
      0.24992392d0,0.61483345d0,1.14319583d0,1.83645455d0,&
      2.69652187d0,3.72581451d0,4.92729377d0,6.30451559d0,&
      7.86169329d0,9.60377599d0,11.53654660d0,13.66674469d0,&
      16.00222119d0,18.55213484d0,21.32720432d0,24.34003576d0,&
      27.60555480d0,31.14158670d0,34.96965201d0,39.11608495d0,&
      43.61365291d0,48.50398616d0,53.84138541d0,59.69912186d0,&
      66.18061779d0,73.44123860d0,81.73681051d0,91.55646652d0,&
      104.15752443d0/)
  double precision,dimension(30),parameter,private::lag_weights=(/1.160441D-01,&
      2.208511D-01,2.413998D-01,1.946368D-01,1.237284D-01,&
      6.367878D-02,2.686048D-02,9.338071D-03,2.680697D-03,&
      6.351291D-04,1.239075D-04,1.982879D-05,2.589351D-06,&
      2.740943D-07,2.332831D-08,1.580746D-09,8.427479D-11,&
      3.485161D-12,1.099018D-13,2.588313D-15,4.437838D-17,&
      5.365918D-19,4.393947D-21,2.311410D-23,7.274588D-26,&
      1.239150D-28,9.832375D-32,2.842324D-35,1.878608D-39,&
      8.745980D-45/)
  contains 
  double precision function uniran()
      !
      !     Random number generator(RCARRY), adapted from F. James
      !     "A Review of Random Number Generators"
      !      Comp. Phys. Comm. 60(1990), pp. 329-344.

      double precision,save::carry
      double precision,dimension(24),save::seeds
      double precision,parameter::one=1
      double precision,parameter::twom24 = ONE/16777216
      integer,save::i,j
      data i, j, carry / 24, 10, 0.0 /
      data seeds / &
      0.8804418, 0.2694365, 0.0367681, 0.4068699, 0.4554052, 0.2880635, &
      0.1463408, 0.2390333, 0.6407298, 0.1755283, 0.7132940, 0.4913043, &
      0.2979918, 0.1396858, 0.3589528, 0.5254809, 0.9857749, 0.4612127, &
      0.2196441, 0.7848351, 0.4096100, 0.9807353, 0.2689915, 0.5140357/

      uniran = seeds(i) - seeds(j) - carry

      if (uniran .lt. 0) then
          uniran = uniran + 1
          carry = twom24
      else
          carry = 0
      end if

      seeds(I) = uniran
      I = 24 - MOD( 25-I, 24 )
      J = 24 - MOD( 25-J, 24 )

  end function uniran
  SUBROUTINE BGOS(SX,ID,X1,X2,RO)

      !C     ID=1:U(0,SX); ID DIFF DE 1 :N(0,SX)
          !use var_surrogate, only: random_generator
      
          implicit none
          double precision ::RO,SX
          integer ::ID
          double precision ::F,V1,V2,S,DLS,RO2
          double precision ::X1,X2!,UNIRAN
          integer::random_generator
          random_generator=2
      !C     !write(*,*)'dans bgos'


      5    CONTINUE

      !C     !write(*,*)'avant rand :'

      !C     X1=RAND()
      !C     X2=RAND()
        ! scl 27/03/2018: remplacement de uniran() par random_number(), pour pouvoir gerer le seed
          if(random_generator==2)then ! on generer avec uniran(mais gestion du seed pas garanti)
              X1=UNIRAN()
              X2=UNIRAN()
          else !on generer avec RANDOM_NUMBER(avec gestion du seed garanti)
              !CALL RANDOM_NUMBER(X1)
              !CALL RANDOM_NUMBER(X2)
              X1=UNIRAN()
              X2=UNIRAN()
          endif

          IF(ID.NE.1) GO TO 10
          F=2.d0*dSQRT(3.d0)
          X1=(X1-0.5)*F
          X2=(X2-0.5)*F
          GO TO 20
      10    CONTINUE
          V1=2.d0*X1-1
          V2=2.d0*X2-1
          S=V1*V1+V2*V2
          IF(S.GE.1.) GO TO 5
          DLS=dSQRT(-2.d0*dLOG(S)/S)
          X1=V1*DLS
          X2=V2*DLS
      20    CONTINUE
          RO2=RO*RO
          IF(ABS(RO).GT.1.E-10) X2=(X1+X2*dSQRT(1.d0/RO2-1.d0))*RO
          X1=X1*SX
          X2=X2*SX

      !C      !write(*,*) 'X1 ',X1,' X2 ',X2
      !C OK, X1 et X2 sont créés

      !C      !write(*,*)'fin bgos'

          RETURN
  END subroutine bgos
  subroutine Cholesky_Factorisation(vc)
      ! VC: matrice de variance-covariance (symetrique) a factoriser
  
      implicit none
      integer :: jj,j,k,ier,maxmes !maxmes= nombre de dimension ou encore dimension de X
      double precision::eps
      double precision,dimension(:,:),intent(inout)::vc
      double precision,dimension(:),allocatable::vi
  
      !=============debut de la fonction=============================
  
      maxmes=size(vc,2)
      allocate(vi(maxmes*(maxmes+1)/2))
      jj=0
      Vi=0.d0
      do j=1,maxmes
          do k=j,maxmes
             jj=j+k*(k-1)/2
             Vi(jj)=VC(j,k)
          end do
      end do
  
      EPS=10.d-10
      CALL DMFSD(Vi,maxmes,eps,ier)! fonction qui fait la factorisation de cholesky
      VC=0.d0
      if (ier.eq.-1) then
          !print*,"Probleme dans la transformation de cholesky pour la generation multinormale"
          ! stop
      else ! on retourne un vecteur de 0 car pas possible de transformer
          do j=1,maxmes
              do k=1,j
                  VC(j,k)=Vi(k+j*(j-1)/2)
              end do
          end do
      end if
  
  end subroutine Cholesky_Factorisation
  subroutine dmfsd(a,n,eps,ier)
          !
          !   FACTORISATION DE CHOLESKY D'UNE MATRICE SDP
          !   MATRICE = TRANSPOSEE(T)*T
          !   ENTREE : TABLEAU A CONTENANT LA PARTIE SUPERIEURE STOCKEE COLONNE
          !            PAR COLONNE DE LA METRICE A FACTORISER
          !   SORTIE : A CONTIENT LA PARTIE SUPPERIEURE DE LA MATRICE triangulaire T
          !
          !   SUBROUTINE APPELE PAR DSINV
          !
          !   N : DIM. MATRICE
          !   EPS : SEUIL DE TOLERANCE
          !   IER = 0 PAS D'ERREUR
          !   IER = -1 ERREUR
          !   IER = K COMPRIS ENTRE 1 ET N, WARNING, LE CALCUL CONTINUE
          !
              implicit none
          
              integer,intent(in)::n
              integer,intent(out)::ier
              double precision,intent(in)::eps
              double precision,dimension(n*(n+1)/2),intent(inout)::A
              double precision :: dpiv,dsum,tol
              integer::i,k,l,kpiv,ind,lend,lanf,lind
          
          !
          !   TEST ON WRONG INPUT PARAMETER N
          !
              dpiv=0.d0
              if (n-1.lt.0) goto 12
              if (n-1.ge.0) ier=0
          !
          !   INITIALIZE DIAGONAL-LOOP
          !
              kpiv=0
              do k=1,n
                  kpiv=kpiv+k
                  ind=kpiv
                  lend=k-1
          !
          !   CALCULATE TOLERANCE
          !
                  tol=dabs(eps*sngl(A(kpiv)))
          !
          !   START FACTORIZATION-LOOP OVER K-TH ROW
          !
                  do i=k,n
                      dsum=0.d0
                      if (lend.lt.0) goto 2
                      if (lend.eq.0) goto 4
                      if (lend.gt.0) goto 2
          !
          !   START INNER LOOP
          !
          2           do l=1,lend
                      lanf=kpiv-l
                      lind=ind-l
                      dsum=dsum+A(lanf)*A(lind)
                      end do
          
          !
          !   END OF INNEF LOOP
          !
          !   TRANSFORM ELEMENT A(IND)
          !
          4           dsum=A(ind)-dsum
                      if (i-k.ne.0) goto 10
                      if (i-k.eq.0) goto 5
          !   TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE
          !
          
          
          5           if (sngl(dsum)-tol.le.0) goto 6
                      if (sngl(dsum)-tol.gt.0) goto 9
          6           if (dsum.le.0) goto 12
                      if (dsum.gt.0) goto 7
          7           if (ier.le.0) goto 8
                      if (ier.gt.0) goto 9
          8           ier=k-1
          !
          !   COMPUTE PIVOT ELEMENT
          !
          9           dpiv=dsqrt(dsum)
                      A(kpiv)=dpiv
                      dpiv=1.D0/dpiv
                      goto 11
          !
          !   CALCULATE TERMS IN ROW
          !
          10          A(ind)=dsum*dpiv
          11          ind=ind+i
                  end do
              end do
          
          !
          !   END OF DIAGONAL-LOOP
          !
              if(ier.eq.-1) then
                  !print*,'Erreur dans le calcul de la cholesky, subroutine dmfsd: ier1=',ier
              end if
          
              return
          12    ier=-1
              !print*,'Erreur dans le calcul de la cholesky, subroutine dmfsd: ier1=',ier
        return
  end subroutine dmfsd
  subroutine compute_rt(b,covar,nparams,zi,times,sizes,meta,res)
      !!$ use omp_lib
      implicit none
      logical::meta
      integer,dimension(9)::nparams
      integer,dimension(4),intent(in)::sizes 
      integer::boot
      integer::nmc
      integer::nboot
      integer::ntimes
      double precision,dimension(sizes(3))::times
      integer::i,k
      double precision::xx1,xx2,xx3
      double precision,dimension(nparams(1))::b
      double precision,dimension(nparams(1),nparams(1))::covar
      double precision,dimension((nparams(2)+4))::zi
      double precision,dimension(nparams(3))::tT,vst,vlt
      double precision,dimension(8)::params 
      double precision,dimension(2,2)::chol1,chol2
      double precision,dimension(sizes(3),4,1+sizes(1)*sizes(4)),intent(out)::res

      !structure sizes
      !(1) = boot (0/1)
      !(2) = nmc
      !(3) = ntimes 
      !(4) = nboot 
      boot = sizes(1) 
      nmc = sizes(2)
      ntimes = sizes(3)
      nboot  = sizes(4)
      !structure nparams 
      !(1) = np 
      !(2) = nz
      !(3) = mt
      !(4) = nvarm 
      !(5) = nvart
      !(6) = rank_trtM
      !(7) = rank_trtT
      !!!!
      call random_params(nparams,b,covar,zi,params,chol1,chol2,tT,vlt,vst,.FALSE.,meta)
      do i=1,ntimes
        call montecarlo(times(i),1,1,nmc,chol1,chol2,nparams,params,tT,vst,meta,xx1)
        call montecarlo(times(i),0,1,nmc,chol1,chol2,nparams,params,tT,vst,meta,xx2)
        call montecarlo(times(i),0,0,nmc,chol1,chol2,nparams,params,tT,vst,meta,xx3)
        res(i,1,1) = times(i)
        res(i,2,1) = xx1
        res(i,3,1) = xx2
        res(i,4,1) = xx3
      enddo
      if(boot.eq.1) then
          do k=1,nboot
              call random_params(nparams,b,covar,zi,params,chol1,chol2,tT,vlt,vst,.TRUE.,meta)
              res(:,1,i+1) = times
              do i=1,ntimes
                  call montecarlo(times(i),1,1,nmc,chol1,chol2,nparams,params,tT,vst,meta,xx1)
                  call montecarlo(times(i),0,1,nmc,chol1,chol2,nparams,params,tT,vst,meta,xx2)
                  call montecarlo(times(i),0,0,nmc,chol1,chol2,nparams,params,tT,vst,meta,xx3)
                  res(i,1,k+1) = times(i)
                  res(i,2,k+1) = xx1
                  res(i,3,k+1) = xx2
                  res(i,4,k+1) = xx3
              enddo
          end do
        end if
      return
  end subroutine
  subroutine montecarlo(t,zm,zt,nsim,chol1,chol2,nparams,params,tT,vst,meta,result)
      !$ use OMP_LIB
      implicit none
      logical::meta
      double precision::t
      integer,dimension(9)::nparams
      integer::zm,zt
      integer::nsim
      integer::i
      double precision,dimension(2,2)::chol1,chol2
      double precision,dimension(:,:),allocatable::w,nu
      double precision::result
      double precision,dimension(nsim)::res
      double precision,dimension(8)::params
      double precision,dimension(nparams(3))::tT,vst
      allocate(w(nsim,2),nu(nsim,2))
      if(meta) then 
          do i=1,nsim
              call bgos(1.0d0,2,w(i,1),w(i,2),0.d0)
              call bgos(1.0d0,2,nu(i,1),nu(i,2),0.d0)
          enddo
          w=matmul(w,chol1)
          nu=matmul(nu,chol2)
      else
          do i=1,nsim
              call bgos(1.0d0,2,w(i,1),w(i,2),0.d0)
              nu(i,1) = 0.d0 
              nu(i,2) = 0.d0 
          enddo
          w=matmul(w,chol1)
      end if 
      result=0.0d0
      !$omp parallel do default(none) private(i) shared(nsim,t,w,nu,zm,zt, &
      !$omp nparams,params,tT,vst,res)
      do i=1,nsim
        call survival_t(t,w(i,1),w(i,2),nu(i,1),nu(i,2),zm,zt,nparams,params,tT,vst,res(i))
      enddo
      !$omp end parallel do
      result=sum(res)/nsim
      deallocate(w,nu)
      return
  end subroutine
  subroutine survival_t(time,w0,w1,nu_m,nu_t,zm,zt,nparams,params,tT,vst,res)
      IMPLICIT NONE
      integer,dimension(9)::nparams
      double precision::time
      double precision::w0,w1,nu_m,nu_t
      double precision,dimension(8)::params
      integer::zm,zt
      integer::i
      double precision,intent(out)::res
      double precision,dimension(nparams(3))::tT,vst
      double precision::x1,x2 
      integer::N 
      res=0.0d0
      if(time.lt.0.0d0) then
        res = 0.d0
        goto 180
      endif
      if(integ_type.eq.1) then 
          N=30
          do i=0,N-1 
              x1 = (i+1)/(N*time)
              x2 = i/(N*time) 
              res = res + (hazard_t(x1,w0,w1,nu_m,nu_t,zm,zt,nparams,params,tT,vst)+&
                           hazard_t(x2,w0,w1,nu_m,nu_t,zm,zt,nparams,params,tT,vst))/2.0d0 
          !s = time*(2*i+1)/400.d0 
          !res = res + (time/200.0d0) * hazard_t(s,w0,w1,nu_m,nu_t,zm,zt,nparams,params,tT,vst)
          end do
          res=(time/N)*res
      else 
          do i=1,30
              if(lag_nodes(i).le.time) then 
                  res = res + hazard_t(lag_nodes(i),w0,w1,nu_m,nu_t,zm,zt,nparams,params,tT,vst)*&
                          dexp(lag_nodes(i))*lag_weights(i)
              end if
          end do
      end if 
      180 continue
      res=dexp(-res)
      return
  end subroutine
  double precision function hazard_t(time,w0,w1,nu_m,nu_t,zm,zt,nparams,params,tT,vst)
      IMPLICIT NONE
      integer,dimension(9)::nparams
      double precision::time
      double precision::predm,predt
      double precision::h0tt
      double precision::w0,w1,nu_m,nu_t
      double precision,dimension(8)::params
      double precision::beta0,beta1,bzm,bzt,bzmt,link,covM,covT 
      double precision,dimension(nparams(3))::tT,vst
      integer::zm,zt
      integer::i

      if(time.lt.0.0d0) then
          hazard_t = 0.d0
          goto 180
      endif
      !(1) = beta0
      !(2) = beta1
      !(3) = bzt
      !(4) = bzm 
      !(5) = link
      !(6) = covM
      !(7) = covT
      beta0 = params(1)
      beta1 = params(2)
      bzt = params(3)
      bzm = params(4)
      bzmt = params(5)
      link = params(6)
      covM = params(7)
      covT = params(8)
      predm=beta0+w0+(beta1+w1)*time + (bzm+nu_m)*zm + bzmt*time*zm + covM
      predt=(bzt+nu_t)*zt + covT
      !h0tt
      if(time.le.0.d0) then
          h0tt = 0.0d0
      else
      if(time.le.tT(1)) then
          h0tt = -dlog(vst(1))
      else if (time.ge.tT(size(tT))) then
          h0tt = -dlog(vst(size(tT)))
      else
          i=0
          do while (tT(i+1).lt.time)
          i = i+1
          end do
          h0tt = -dlog(vst(i))
      end if
      end if

      hazard_t=h0tt*dexp(predt+link*predm)
      !print*,'hazard_t=',hazard_t,time,tT(1),tT(size(tT))
      180 continue
      return
  end function
  subroutine random_params(nparams,b,covar,zi,params,chol1,chol2,tT,vlt,vst,ran,meta)
      implicit none
      logical::meta
      integer::np,mt
      integer::rank_trtM,rank_trtT,rank_trtMt,rank_time
      integer,dimension(9)::nparams 
      double precision,dimension(nparams(1)),intent(in)::b 
      double precision,dimension(nparams(1))::b_ind,x1
      double precision,dimension(nparams(1))::b2
      double precision,dimension(nparams(1),nparams(1)),intent(in)::covar 
      double precision,dimension(nparams(1),nparams(1))::cdummy
      double precision,dimension(nparams(1)+4)::zi
      double precision,dimension(8)::params 
      double precision,dimension(nparams(3))::tT,vlt,vst 
      double precision ,dimension(-2:nparams(2))::the
      double precision,dimension(2,2)::chol1,chol2
      double precision::dummy
      logical::ran
      integer::nz,nvarm,nvart
      integer::i
      np = nparams(1)
      nz = nparams(2)
      mt = nparams(3)
      nvarm = nparams(4)
      nvart = nparams(5)
      rank_trtM = nparams(6)
      rank_trtT = nparams(7)
      rank_trtMt = nparams(8)
      rank_time = nparams(9)
      !structure nparams 
      !(1) = np 
      !(2) = nz
      !(3) = mt
      !(4) = nvarm 
      !(5) = nvart

      !structure params
      !(1) = beta0
      !(2) = beta1
      !(3) = bzt
      !(4) = bzm 
      !(5) = link
      !(6) = covM
      !(7) = covT

      !beta0 = params(1)
      !beta1 = params(2)
      !bzt = params(3)
      !bzm = params(4)
      !bzmt = params(5)
      !link = params(6)
      !covM = params(7)
      !covT = params(8)
      if(.NOT.ran)then
          b2 = b 
          do i=1,nz+2
              the(i-3)=b2(i)**2
          end do 
          call baseline_rs(nparams(2),nparams(3),zi,the,tT,vlt,vst)
          chol1(1,:) = (/b2(nz+5),b2(nz+6)/)
          chol1(2,:) = (/0.d0,b2(nz+7)/)
          params(1) = b2(nz+7+nvart+1)
          params(2) = b2(nz+7+nvart+rank_time)
          params(3) = b2(nz+7+rank_trtT)
          params(4) = b2(nz+7+nvart+rank_trtM)
          if(rank_trtMt.ne.0) then 
            params(5) = b2(nz+7+nvart+rank_trtMt)
          else 
            params(5) = 0.0d0 
          end if 
          params(6) = b2(nz+3)
          params(7) = 0.d0 !covM
          params(8) = 0.d0 !covT
          chol2 = 0.d0 
          if(meta) then 
              chol2(1,:) = (/b2(np-2),b2(np-1)/)
              chol2(2,:) = (/0.d0,b2(np)/)
          end if
      else 
          cdummy=covar
          call Cholesky_Factorisation(cdummy)
          do i=1,np
              call bgos(1.0d0,2,x1(i),dummy,0.d0)
          enddo
          x1=matmul(x1,transpose(cdummy))
          do i=1,np
              b_ind(i) = b(i) + x1(i)
          end do
          do i=1,nz+2
            !the(i-3)=b_ind(i)**2
            the(i-3)=b(i)**2
          end do 
        call baseline_rs(nparams(2),nparams(3),zi,the,tT,vlt,vst)
        chol1(1,:) = (/b_ind(nz+5),b_ind(nz+6)/)
        chol1(2,:) = (/0.d0,b_ind(nz+7)/)
        params(1) = b_ind(nz+7+nvart+1)
        params(2) = b_ind(nz+7+nvart+rank_time)
        params(3) = b_ind(nz+7+rank_trtT)
        params(4) = b_ind(nz+7+nvart+rank_trtM)
        if(rank_trtMt.ne.0) then 
          params(5) = b_ind(nz+7+nvart+rank_trtMt)
        else 
          params(5) = 0.0d0 
        end if 
        params(6) = b_ind(nz+3)
        params(7) = 0.d0 !covM
        params(8) = 0.d0 !covT
        chol2 = 0.d0 
        if(meta) then 
            chol2(1,:) = (/b_ind(np-2),b_ind(np-1)/)
            chol2(2,:) = (/0.d0,b_ind(np)/)
        end if
      end if
  end subroutine 
  subroutine baseline_rs(nz,mt,zi,theT,tT,vlt,vst)  
      Implicit none
      integer,intent(in)::nz
      integer::i,n,k,mt
      double precision::x1,h,su,lam
      double precision ,dimension(-2:nz)::theT
      double precision,dimension((nz+6))::zi
      double precision,dimension(mt)::tT,vlt,vst 
      n  = nz+2
      k = 0
      h = (zi(n)-zi(1))*0.01d0
      x1 = zi(1)
          do i=1,mt
              if(i .ne.1)then
                  x1 = x1 + h
              end if
              call calcul_rs(x1,theT,nz+2,zi,su,lam)
              tT(i)=x1
              vlt(i)=lam
              vst(i)=su
          end do
      return
  end subroutine baseline_rs
  subroutine calcul_rs(x,the,n,zi,su,lam)
      IMPLICIT NONE
      integer,intent(in)::n
      double precision,intent(in)::x
      double precision,intent(out)::lam,su
      double precision,dimension(-2:(n-2)),intent(in)::the
      double precision,dimension((n+4))::zi
      integer::j,k,i
      double precision::ht,ht2,h2,som,htm,h2t,h3,h2n,hn, &
      im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2, &
      h,gl,hh
      j=0
      gl=0.d0
      som = 0.d0
      do k = 2,n-1
          if ((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
              j = k-1
              if (j.gt.1)then
                  do i=2,j
                  som = som+the(i-4)
                  end do
              endif
              ht = x-zi(j)
              htm= x-zi(j-1)
              h2t= x-zi(j+2)
              ht2 = zi(j+1)-x
              ht3 = zi(j+3)-x
              hht = x-zi(j-2)
              h = zi(j+1)-zi(j)
              hh= zi(j+1)-zi(j-1)
              h2= zi(j+2)-zi(j)
              h3= zi(j+3)-zi(j)
              h4= zi(j+4)-zi(j)
              h3m= zi(j+3)-zi(j-1)
              h2n=zi(j+2)-zi(j-1)
              hn= zi(j+1)-zi(j-2)
              hh3 = zi(j+1)-zi(j-3)
              hh2 = zi(j+2)-zi(j-2)
              mm3 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
              mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
              *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
              mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
              h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
              mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
              im3 = (0.25d0*(x-zi(j-3))*mm3)+(0.25d0*hh2*mm2) &
              +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
              im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
              im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
              im  = ht*mm*0.25d0
              gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)+(the(j)*im)
              lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)+(the(j)*mm)
          endif
      end do
      if(x.ge.zi(n))then
          som = 0.d0
          do i=1,n
              som = som+the(i-3)
          end do
          gl = som
      endif
      su  = dexp(-gl)
      return
  end subroutine calcul_rs
end module 