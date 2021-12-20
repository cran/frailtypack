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
  double precision,dimension(:),allocatable,private::vls,vss,vlt,vst,tS,tT
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
      integer::nsplines
      integer::boot
      double precision,dimension(ntimes)::times
      double precision,dimension(ntimes,4,1+boot*nboot),intent(out)::res
      double precision,dimension(2,2)::sigma_ind
      double precision,dimension(2+nparams(7))::knots
      integer::i,j
      integer::integ
      double precision::s
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
      allocate(vls(nparams(9)),vss(nparams(9)),vlt(nparams(10)),vst(nparams(10)))
      allocate(tS(nparams(9)),tT(nparams(9)))
      allocate(bgamma(nsplines))
      allocate(b_ind(np))
      allocate(basissplines(nparams(7)+nparams(8)))
      allocate(vectgamma(1001,2))

      integ_type=integ
      nsplines = nparams(7)+nparams(8)

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
      deallocate(tS,tT,vls,vss,vlt,vst)
      return
  end subroutine
  subroutine montecarlo(t,zs,zt,nsim,result)
       !$ use OMP_LIB
       implicit none
       double precision::t
       integer::zs,zt
       integer::nsim
       integer::i
       double precision,dimension(2,2)::sigma,sigma2
       double precision,external::h
       double precision,dimension(:,:),allocatable::x1
       double precision::result
       double precision,dimension(nsim)::res

       sigma(1,1)=sig_omega+sig_u+sig_s+zs*sig_s
       sigma(2,2)=(eta_ind**2.d0)*sig_omega+(alpha_ind**2.d0)*sig_u+sig_t
       sigma(1,2)=eta_ind*sig_omega+alpha_ind*sig_u+zs*zt*sig_st
       sigma(2,1)=sigma(1,2)

       call Cholesky_Factorisation(sigma)
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
    double precision::dummy
    integer::zs,zt
    integer::i,j,indic_int
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
      if(s.le.tS(1)) then
        h0ss = -dlog(vss(1))
        l0ss = vls(1)
      else if (s.ge.tS(size(tS))) then
        h0ss = -dlog(vss(size(tS)))
        l0ss = vls(size(tS))
      else
        i=0
        do while (tS(i+1).lt.s)
          i = i+1
        end do
        h0ss = -dlog(vss(i))
        l0ss = vls(i)
      end if
    end if
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
        if(s.le.tT(1)) then
          h0ts = -dlog(vst(1))
        else if (s.ge.tT(size(tT))) then
          h0ts = -dlog(vst(size(tT)))
        else
          i=0
          do while (tT(i+1).lt.s)
            i = i+1
          end do
          h0ts = -dlog(vst(i))
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
        if(time.le.tS(1)) then
          h0st = -dlog(vss(1))
        else if (time.ge.tS(size(tS))) then
          h0st = -dlog(vss(size(tS)))
        else
          i=0
          do while (tS(i+1).lt.time)
            i = i+1
          end do
          h0st = -dlog(vss(i))
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
    integer::i,j,nsplines
    nsplines = nparams(7)+nparams(8)
    if(.NOT.ran)then
      b_ind = b
      bgamma=b_ind((np-nsplines+1):np)
      vectgamma = 0.0d0
      basissplines=0.0d0
      bs=b_ind(np-nsplines-1)
      bt=b_ind(np-nsplines)
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
                       nparams(10),zi,tS,vls,vss,tT,vlt,vst)
    end if
    if(ran) then
      cdummy=covar
      call Cholesky_Factorisation(cdummy)
      do i=1,np
        call bgos(1.0d0,2,x1(i),dummy,0.d0)
      enddo
      x1=matmul(x1,transpose(cdummy))
      do i=1,np
        b_ind(i) = b(i) + x1(i)
      end do
      bgamma=b_ind((np-nsplines+1):np)
      vectgamma = 0.0d0
      basissplines=0.0d0
      bs=b_ind(np-nsplines-1)
      bt=b_ind(np-nsplines)
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
                       nparams(10),zi,tS,vls,vss,tT,vlt,vst)
    end if
  end subroutine
  subroutine baseline_rs(nz1,nz2,b,mt1,mt2,zi,tS,vls,vss,tT,vlt,vst)
    Implicit none
    integer,intent(in):: nz1,nz2,mt1,mt2
    double precision ,dimension(npmax),intent(in):: b
    integer::i,j,n,k,l,jj
    double precision::x1,x2,h,su,lam
    double precision ,dimension(-2:nz1)::the2
    double precision ,dimension(-2:nz2)::the1T
    double precision,dimension(-2:(nz1+3))::zi
    double precision,dimension(mt1)::tS
    double precision,dimension(mt2)::tT
    double precision,dimension(mt1)::vls,vss
    double precision,dimension(mt2)::vlt,vst
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
            tS(i)=x1
            vls(i)=lam
            vss(i)=su
        end do
    x2 = zi(1)
    do i=1,mt2
        if(i.ne.1)then
            x2 = x2 + h
        endif
        call calcul_rs(x2,the2,nz2+2,zi,su,lam)
        tT(i)=x2
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
    double precision,dimension(-2:(n+1))::zi
    integer::j,k,i
    double precision::ht,ht2,h2,som,pm,htm,h2t,h3,h2n,hn, &
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
