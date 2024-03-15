!========================          FUNCPAJ_SPLINES         ====================
double precision function funcpajlongisplines2(b,np,id,thi,jd,thj,k0)
  !use donnees, only:MC1,MC2,MC3,MC4,MC5,MC6,MC7,MC8,MC9,MC10,MC11,MC12,MC13,&
  !MC14,MC15,MC16,MC17,MC18,MC19,MC20,MC21,MC22,MC23,MC24,MC25,&
  !x5,w5,x7,w7,x9,w9,x12,w12,x15,w15,x3,w3!,x2,w2,x3,w3
  use donnees_indiv
  use lois_normales
  use tailles
  use comon
  use Autres_fonctions, only:init_random_seed, pos_proc_domaine, bgos, uniran,&
                             Cholesky_Factorisation! Monte-carlo random generation
  use var_surrogate, only: mediation,Chol! a_deja_simul,nbre_sim,Chol,Vect_sim_MC,graine,aleatoire,&
                           !mediation!,frailt_base,nb_procs
  use var_mediation,only:nuzm,nuzt,treat_ind,center_ind,&
                         nparammed,ncenters,nmcmed,niter,&
                         method_int,type_ma,matmc,betazm,&
                         betaztime,betazt,Cmult,matmc_frail,nmcfrail
  !use ParametresPourParallelisation
  use residusM
  use optim
  !$ use OMP_LIB
  IMPLICIT NONE
  ! *** NOUVELLLE DECLARATION F90 :
  integer,intent(in)::id,jd,np
  double precision,dimension(np),intent(in)::b
  double precision,dimension(2),intent(in)::k0
  double precision,intent(in)::thi,thj
  ! for the numerical integral hrmsym
  integer :: restar,nf2,jj,ier
  double precision:: epsabs,epsrel
  !double precision,dimension(2):: result, abserr2
  !double precision,dimension(1000) :: work
  external :: vraistot,vraistot_splines,vraistot_weib
  double precision,dimension(nea) :: xea
  !integer ::neval,ifail
  integer::n,i,j,k,vj,ig,choix,l!,it
  integer,dimension(ngmax)::cpt
  double precision::pe1,pe2,som1,som2,res,vet2,h1!,vet,sum
  double precision :: eps_s
  double precision,dimension(-2:npmax):: the1,the2
  double precision,dimension(np)::bh
  double precision,dimension(ngmax)::res2,res1dc,res2dc &
  ,res3dc,integrale1,integrale2,integrale3,integrale4
  !AD: for death,change dimension
  double precision,dimension(ndatemax)::dut1
  double precision,dimension(ndatemaxdc)::dut2
  !AD:end
  double precision,dimension(0:ndatemax)::ut1
  double precision,dimension(0:ndatemaxdc)::ut2
  double precision::int,eps
  double precision,parameter::pi=3.141592653589793d0
  double precision,dimension(:,:),allocatable :: mat_sigma,varcov_marg_inv
  double precision,dimension(:),allocatable :: matv
  double precision,dimension(nva3,nva3) :: element
  !double precision,dimension(:,:),allocatable :: mat_sigmaB, varcov_marg_invB ! add TwoPart
  !double precision,dimension(:),allocatable :: matvB ! add TwoPart
  !double precision,dimension(nvaB,nvaB) :: elementB ! add TwoPart
  !double precision,dimension(3):: resultatInt ! add Monte-carlo
  double precision::func8J,func9J,func10J,func11J, funcTP4J,funcG
  external::func8J,func9J,func10J,func11J, funcTP4J,funcG ! add Monte-carlo
  !integer::vcdiag ! add Monte-carlo
  !double precision,dimension(nb1)::mu_mc
  !double precision,dimension(nb1,nb1)::vcjm
  double precision,dimension(nodes_number,nb1)::fraili
  !double precision::SX,xMC ! for random generation
  !double precision,dimension(25000)::MC ! for Monte-carlo pre-generated points
  !integer::m 
  !mediation~meta-analysis
  integer::ll,itcenter!,itcenter1,itcenter2
  !integer::nre_center=3
  integer::kk
  !double precision,dimension(nmcmed,2)::matnu
  double precision,dimension(2,2)::matrecenter!,sigmaMC,invmatrecenter
  double precision,dimension(ncenters)::rescenter,nindcenter
  double precision::resG1,resG2,resG3,resG4
  !double precision,dimension(:,:),allocatable::matmc !matrix for MC integration
  double precision::resCcurr,res_indiv,res4re
  !double precision,dimension(nquad)::xherm,wherm
  !double precision::resherm,resherm1,resherm2,detmatrecenter,densrecenter
  !double precision::varnu_m,varnu_t,dummy
  double precision::minint,maxint,minres2dc,maxres2dc
  integer::ncenters2! sumtrtid
  double precision,dimension(2)::vectnu
  !double precision,dimension(5)::rescenter2,rescenter2p1,rescenter2p2,rescenter2p3
  !double precision::rescenterp1,rescenterp2,rescenterp3
  !integer::sumnmes
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !double precision,dimension(:,:),allocatable::bb1,bb2,grandb
  !double precision :: sigmav,range
  double precision,dimension(:,:),allocatable::muOMP!,mu1
  double precision,dimension(:,:),allocatable::Z1OMP!,Z2OMP
  integer::numpatOMP,nmescurOMP,it_curOMP,nmescurrOMP,nmescurr1OMP!,nmesOMP,nmescur2OMP
  !integer,parameter ::nf=1
  !double precision,dimension(:),allocatable,save::b1
  double precision,dimension(:),allocatable::ycurrentOMP!,current_meanOMP!,ycurrent2,
  !double precision,dimension(:),allocatable,save :: xeacurrent,part
  !double precision,dimension(:,:),allocatable:: x2,x22,z22,z11,x2cur,z1cur
  !integer,dimension(:),allocatable,save:: nii,nii2,nmes_o,nmes_o2
  !integer,save:: it_rec
  !double precision,dimension(ndatemaxdc)::intsurv
  double precision:: ut2curOMP!,frailpolOMP,frailpol2OMP,frailpol3OMP,frailpol4OMP
  !double precision,dimension(:),allocatable::res1cur,res3cur,res2cur
  !double precision,dimension(:,:),allocatable::Z1B, muB,XB,mu1B,x2Bcur,z1Bcur ! add TwoPart
  !double precision,dimension(:),allocatable:: Bcurrent, current_meanRaw ! add TwoPart
  !integer::nmescurB, it_curB, interceptBin !add TwoPart
  !double precision::fixed_Binary
  !integer::GLMloglink0,MTP0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real::tstart,tstop  

  !$ integer::nthreads 
  method_int=1
  
  !$ nthreads = omp_get_max_threads()

  !call cpu_time(tstart)
  niter=niter+1
  minint=0.0d0
  maxint=0.0d0
  minres2dc=0.0d0
  maxres2dc=0.0d0
  ncenters2 = 5
  !double precision,dimension(1,1)::densrecenter
  npp = np
  eps_s = 1.d-7
  choix=0
  ig=0
  k=0
  vj=0
  n=0
  j=0
  ut1=0.d0
  ut2=0.d0
  dut2=0.d0
  dut1=0.d0
  do i=1,np
    bh(i)=b(i)
  end do
  fraili=0.d0
  if (id.ne.0) bh(id)=bh(id)+thi
  if (jd.ne.0) bh(jd)=bh(jd)+thj
  b1 = bh
  n = (np-nva-effet-indic_ALPHA-1-nb_re-netadc-netar-nparammed)/(effet+1) !nst        !to znaczy ze dzielimy lliczbe wezlow na 2
  the1 = 0.d0
  the2 = 0.d0
  if(typeJoint.ne.2) then
    do i=1,n
      the1(i-3)=(bh(i))*(bh(i))
      j = n+i
      if (nst.eq.2) then
        the2(i-3)=(bh(j))*(bh(j))
      endif
    end do
  else
    do i=1,n
      the2(i-3)=(bh(i))*(bh(i))
    end do
  end if
  if(typeJoint.eq.3) then
    sigmav = bh(np-nva-nb_re-1-netadc - netar-1)
    alpha = bh(np-nva-nb_re-1-netadc - netar)
  else
    sigmav = -1.d0
  end if
  if(nea.ge.1) then
    etaydc(1:netadc) = bh(np-nva-nb_re-netadc-nparammed:np-nva-nb_re-nparammed-1)
    etayr(1:netar) =bh(np-nva-nb_re-netadc-netar-nparammed:np-nva-nb_re-netadc-nparammed-1)
    sigmae = bh(np-nva-nb_re-nparammed)*bh(np-nva-nb_re-nparammed)

    Ut = 0.d0
    Utt = 0.d0
    do j=1,nb1
      do k=1,j
        if(j.eq.k) then
          Ut(j,k)=sqrt(bh(np-nva-nb_re-nparammed+k+j*(j-1)/2)**2.d0)
          Utt(k,j)=sqrt(bh(np-nva-nb_re-nparammed+k+j*(j-1)/2)**2.d0)
        else
          Ut(j,k)=bh(np-nva-nb_re-nparammed+k+j*(j-1)/2)
          Utt(k,j)=bh(np-nva-nb_re-nparammed+k+j*(j-1)/2)
        end if
      end do
    end do
    if(typeJoint.eq.3) then
      Ut(nea,nea) = sqrt(sigmav**2.d0)
      Utt(nea,nea) = sqrt(sigmav**2.d0)
    end if
  end if
  !----------  calcul de ut1(ti) et ut2(ti) ---------------------------
  !    attention the(1)  sont en nz=1
  !        donc en ti on a the(i)
  if(effet.eq.1) then
    dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))
    ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))
  end if
  dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))
  ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))
  ut1(0) = 0.d0
  ut2(0) = 0.d0
  !//// NEW AMADOU vvv :
  !--- strate1
  som1 = 0.d0
  vj = 0
  if(effet.eq.1) then
    do i=2,ndate-1
      do k = 2,n-2
        if (((date(i)).ge.(zi(k-1))).and.(date(i).lt.zi(k)))then
          j = k-1
          if ((j.gt.1).and.(j.gt.vj))then
            som1 = som1 + the1(j-4)
            vj  = j
          endif
        endif
      end do
      ut1(i) = som1 +(the1(j-3)*im3(i))+(the1(j-2)*im2(i)) &
      +(the1(j-1)*im1(i))+(the1(j)*im(i))
      dut1(i) = (the1(j-3)*mm3(i))+(the1(j-2)*mm2(i)) &
      +(the1(j-1)*mm1(i))+(the1(j)*mm(i))
    end do
  end if
  !--- strate2
  vj = 0

  som2 = 0.d0
  do i=2,ndatedc-1
    do k = 2,n-2
      if (((datedc(i)).ge.(zi(k-1))).and.(datedc(i).lt.zi(k)))then
        j = k-1
        if ((j.gt.1).and.(j.gt.vj))then
          som2 = som2 + the2(j-4)
          vj  = j
        endif
      endif
    end do
    ut2(i) = som2 +(the2(j-3)*im3dc(i))+(the2(j-2)*im2dc(i)) &
    +(the2(j-1)*im1dc(i))+(the2(j)*imdc(i))
    dut2(i) = (the2(j-3)*mm3dc(i))+(the2(j-2)*mm2dc(i)) &
    +(the2(j-1)*mm1dc(i))+(the2(j)*mmdc(i))
  end do
  i = n-2
  h1 = (zi(i)-zi(i-1))
  if(effet.eq.1) then
    ut1(ndate)=som1+the1(i-4)+the1(i-3)+the1(i-2)+the1(i-1)
    dut1(ndate) = (4.d0*the1(i-1)/h1)
  end if
  ut2(ndatedc)=som2+the2(i-4)+the2(i-3)+the2(i-2)+the2(i-1)!am the1(i-4)
  dut2(ndatedc) = (4.d0*the2(i-1)/h1)

  if(mediation) then
    if(type_ma.eq.3) then
      !b = (/.....,a1,a2,a3,betazm,betazt,link)
      if(link.eq.4) then
        matrecenter(1,1) = bh(np)
        matrecenter(1,2) = 0.d0
        matrecenter(2,2) = bh(np)
        matrecenter(2,1) = bh(np-1)
        !varnu_t=bh(np-3)*2.0d0
        !betazm=bh(np-1)
        !betazm=bh(np-2)
        !betazt=bh(np-1)
        !betaztime=bh(np)
      else
        matrecenter(1,1) = bh(np-2)
        matrecenter(1,2) = 0.d0
        matrecenter(2,2) = bh(np)
        matrecenter(2,1) = bh(np-1)
        !varnu_t=bh(np-3)*2.0d0
        !betazm=bh(np-1)

        !betazm=bh(np-2)
        !betazt=bh(np-1)
        !betaztime=bh(np)
      end if
      !ma_link=bh(np)
      !betazt=0.d0
      !ma_link=0.0d0
      !Afin d'éviter des pb de matrice non definie matrecenter
      !correspond à la cholesky (supérieure) de la matrice de covariance
      ! i.e cov = LL^t ou L=matrecenter
      !Ainsi, var(nuzm) =  bh(np-4)**2+bh(np-3)**2
      !       var(nuzt) =  bh(np-2)**2
      !       cov(nuzm,nuzt) = bh(np-3)*bh(np-2)
      !matnu = matmul(matmc,matrecenter)
      res=0.0d0

      resG3 = 0.d0
      resG4 =  0.0d0

      do ll=1,ncenters
        rescenter(ll) = 0.d0
        resG1 = 0.0d0
        resG2 = 0.0d0
        nindcenter(ll) = 0
        do k=1,ng
          if(center_ind(k).eq.ll) then
            nindcenter(ll) = nindcenter(ll) +1
          end if
        end do

        !maxint=0.d0
        !maxres2dc=0.0d0
        if(method_int.eq.1) then
          !call omp_set_num_threads()

          !!private(nuzm,nuzt,itcenter,&
          !!$omp&  resCcurr,k,res_indiv,ig,res1,res2,res3,res1dc,res2dc,&
          !!$omp&   res3dc,res4re,cpt,integrale1,integrale3,integrale4,aux1,aux2,&
          !!$omp&   vet2,RisqCumul,Rdc_res,it,it_rec,epsabs,epsrel,restar,nf2,nmes_o,&
          !!$omp&   sum_mat,mat_sigma,ycurrent,auxig,numpat,nmescur,nmescurr,nmescurr1,&
          !!$omp&   x2,x2cur,z1cur,current_mean,res1cur,res2cur,res3cur,Z1,l,kk,&
          !!$omp&   varcov_marg,matv,j,jj,ier,eps,element,mu,xea,ut2cur,choix,int),&


          !$omp parallel do default(private), shared(nmcmed,ng,matmc,matrecenter,&
          !$omp&     center_ind,ll,treat_ind,bh,np,npp,nva,nva2,nva3,nvaB,nparammed,vedc,&
          !$omp&     cdc,dut2,nt1dc,ut2,nea,nb1,Chol,nmesy,nmesrec,nmesrec1,&
          !$omp&     yy,typeJoint,s_cag_id,sigmae,s_cag,ziy,vey,nby,Ut,Utt,&
          !$omp&     TwoPart,nodes_number,methodGH,method_GH,nf2,genz,&
          !$omp&     epsabs,epsrel,restar,nmesB,matmc_frail,nmcfrail,nthreads,&
          !$omp&     invBi_cholDet,invBi_chol,b1,mediation,t0dc,t1dc,&
          !$omp&     betazt,betazm,betaztime,effet,link,maxmesy,Cmult), reduction(+:resG1,resG2,resG3,resG4)
          do itcenter=1,nmcmed
            !call dblepr('omp iter 1',-1,1.0d0,1)
            !allocate(nmes_o(ng))
            allocate(Z1OMP(maxmesy,nb1))
            allocate(muOMP(maxmesy,1),ycurrentOMP(maxmesy))
            vectnu = matmul(matrecenter,matmc(itcenter,:))
            !nuzm=matnu(itcenter,1)
            !nuzt=matnu(itcenter,2)

            nuzm=vectnu(1)
            nuzt=vectnu(2)

            do k=1,ng
              if(center_ind(k).eq.ll) then
                !calcul vraisemblance individuelle sujet k
                res_indiv=0.0d0
                ig=k
                res1(k) = 0.d0
                res2(k) = 0.d0
                res3(k) = 0.d0
                res1dc(k) = 0.d0
                res2dc(k) = 0.d0
                res3dc(k) = 0.d0
                res4re = 0.0d0
                cpt(k) = 0
                integrale1(k) = 0.d0
                integrale2(k) = 0.d0
                integrale3(k) = 0.d0
                integrale4(k) = 0.0d0
                aux1(k)=0.d0
                aux2(k)=0.d0

                !ccccccccccccccccccccccccccccccccccccccccc
                ! pour le deces
                !ccccccccccccccccccccccccccccccccccccccccc
                !do k=1,ng

                if(nva2.gt.0)then
                  vet2 = 0.d0
                  do j=1,nva2
                    vet2 = vet2 + bh(np-nva3-nva2-nvaB-nparammed+j)*dble(vedc(k,j))
                  end do
                  vet2 = dexp(vet2)
                else
                  vet2=1.d0
                endif
                if(cdc(k).eq.1)then
                  res2dc(k) = dlog(dut2(nt1dc(k))*vet2)
                  if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
                    funcpajLongisplines2=-1.d9
                    !goto 123
                  end if
                endif

                RisqCumul(k) = ut2(nt1dc(k))*vet2
                Rdc_res(k) = ut2(nt1dc(k))!*vet2
              end if
            end do
            !**************INTEGRALES ****************************
            if(nea.ge.1) then
              it = 0
              it_rec = 1
              epsabs = 1.d-100
              epsrel = 1.d-100
              restar = 0
              nf2 = nf
              nmes_o=0
              integrale4 = 0.d0
              sum_mat= 0.d0
              if(nb1.eq.2) then
                Chol=0.d0
                Chol(1,1)=bh(np-nva-nb_re-nparammed+1)
                Chol(2,1)=bh(np-nva-nb_re-nparammed+2)
                !Chol(1,2)=bh(np-nva-nb_re+2)
                Chol(2,2)=bh(np-nva-nb_re-nparammed+3)
              end if
              do ig=1,ng
                !if(center_ind(ig).eq.ll) then
                  ycurrentOMP  = 0.d0
                  auxig=ig
                  choix = 4
                  numpatOMP = ig
                  nmescurOMP =nmesy(ig)
                  nmescurrOMP = nmesrec(ig)
                  nmescurr1OMP = nmesrec1(ig)
                  it_curOMP = it
                  allocate(mat_sigma(nmescurOMP,nmescurOMP))
                  x2 = 0.d0
                  x2cur = 0.d0
                  z1cur = 0.d0
                  current_mean = 0.d0
                  mat_sigma = 0.d0
                  if(nmescurOMP.gt.0) then
                    do i= 1,nmescurOMP
                      ycurrentOMP(i) = yy(it+i)
                      mat_sigma(i,i) = sigmae!**2.d0 !sigma is the variance ?!
                      if(s_cag_id.eq.1)then
                        if(ycurrentOMP(i).gt.s_cag) then
                          nmes_o(ig) = nmes_o(ig)+1
                        end if
                      else
                        nmes_o(ig) = nmescurOMP
                      end if
                    end do
                  else
                    nmes_o(ig) = nmescurOMP
                  end if
                  res1cur = 0.d0
                  res2cur = 0.d0
                  res3cur = 0.d0
                  if(typeJoint.eq.3) then
                    do i= 0,nmescurrOMP-1
                      res1cur(i+1) = ut1(nt1(it_rec+i))
                      res3cur(i+1)  = ut1(nt0(it_rec+i))
                      res2cur(i+1) = dut1(nt1(it_rec+i))
                    end do
                  end if
                  ! creation de Zi
                  Z1OMP=0.d0
                  l=0
                  if(nmescurOMP.gt.0) then
                    do kk=1,nby
                      l=l+1
                      do i=1,nmescurOMP
                        Z1OMP(i,l)=dble(ziy(it+i,kk))
                      end do
                    end do
                  else
                    do i=1,nmescurOMP
                      Z1OMP(i,1)=0.d0
                    end do
                  end if
                  l = 0
                  X2 = 0.d0
                  if(nmescurOMP.gt.0) then
                    do kk=1,nva3
                      l = l + 1
                      do j=1,nmescurOMP
                        X2(j,l) = dble(vey(it+j,kk))
                      end do
                    end do
                  end if
                  if(nmescurOMP.gt.0) then
                    varcov_marg((it+1):(it+nmescurOMP),1:nmescurOMP) =Matmul( MATMUL(ziy((it+1):(it+nmescurOMP),1:nby), &
                    MATMUL(Ut(1:nby,1:nby),Utt(1:nby,1:nby))),transpose(ziy((it+1):(it+nmescurOMP),1:nby)))+ &
                    mat_sigma
                  end if
                  if(nmescurOMP.gt.0) then
                    allocate(matv(nmescurOMP*(nmescurOMP+1)/2),varcov_marg_inv(nmescurOMP,nmescurOMP))
                    matv = 0.d0
                    do j=1,nmescurOMP
                      do kk=j,nmescurOMP
                        jj=j+kk*(kk-1)/2
                        matv(jj)=varcov_marg(it+j,kk)
                      end do
                    end do
                    ier = 0
                    eps = 1.d-10
                    call dsinvj(matv,nmescurOMP,eps,ier)
                    varcov_marg_inv=0.d0
                    do j=1,nmescurOMP
                      do kk=1,nmescurOMP
                        if (kk.ge.j) then
                          varcov_marg_inv(j,kk)=matv(j+kk*(kk-1)/2)
                        else
                          varcov_marg_inv(j,kk)=matv(kk+j*(j-1)/2)
                        end if
                      end do
                    end do
                    element=0.d0
                    element =  Matmul(Matmul(Transpose(vey(it+1:it+nmescurOMP,1:nva3)), &
                    varcov_marg_inv(1:nmescurOMP,1:nmescurOMP)), vey(it+1:it+nmescurOMP,1:nva3))
                    do j=1,nva3
                      do kk=1,nva3
                        sum_mat(j,kk) = sum_mat(j,kk) +  element(j,kk)
                      end do
                    end do
                    deallocate(matv,varcov_marg_inv)
                    muOMP = 0.d0
                    if(TwoPart.eq.0) then
                      muOMP(1:nmescurOMP,1) = matmul(X2(1:nmescurOMP,1:(nva3)),bh((np-nva3-nparammed+1):(np-nparammed)))&
                      +(nuzm)*treat_ind(ig) ! +
                      if((link.eq.4).or.(link.eq.2))then
                        do j=1,nmescurOMP
                          !muOMP(j,1) = muOMP(j,1) + X2(j,3)*betaztime*treat_ind(ig)
                        end do
                      end if
                      !mu(1:nmescur,1) = matmul(X2(1:nmescur,1:(nva3)),bh((np-nva-nparammed+1):(np-nparammed)))
                    else if(TwoPart.eq.1) then
                      muOMP(1:nmescurOMP,1) = matmul(X2(1:nmescurOMP,1:(nva3)),bh((np-nva3-nvaB-nparammed+1):(np-nvaB-nparammed)))
                    end if
                    xea = 0.d0
                  else
                    muOMP = 0.d0
                    do kk=1,nb1
                      Z1OMP(:,kk)=0.d0
                    end do
                  end if
                  ut2curOMP = ut2(nt1dc(ig))
                  choix = 3
                  if(center_ind(ig).eq.ll) then
                    if(method_GH.le.1) then
                      if(nmesy(numpatOMP).gt.0) then
                        allocate(mu1(nmesy(numpatOMP),1))
                      else
                        allocate(mu1(1,1))
                      end if
                      if(typeJoint.eq.2.and.nb1.eq.1) then
                      else if(typeJoint.eq.2.and.nb1.eq.2) then
                        !if(itcenter.eq.1) then
                          !nodes_number=1
                          !call gauherJ22OMP(int,choix,nodes_number,numpatOMP,&
                          !                  nmescurOMP,ut2curOMP,muOMP,ycurrentOMP,it_curOMP,Z1OMP)
                          call dummyMCOMP(int,numpatOMP,&
                          nmescurOMP,ut2curOMP,muOMP,ycurrentOMP,it_curOMP,Z1OMP)
                        !end if
                      end if
                      deallocate(mu1)
                      integrale4(ig) =int !result(1) !
                    end if
                  end if
                  it_rec = it_rec + nmescurOMP
                  it = it + nmescurOMP
                  if(integrale4(ig).gt.1.E+30) then
                    integrale4(ig) = 1.E+30
                  end if
                  deallocate(mat_sigma)
                !end if
              end do
            else
              sigmae = 1.d0
            end if
            resCcurr=1.0d0
            do k=1,ng
              if(center_ind(k).eq.ll) then
                if(nb1.ge.1) then
                  nmescurOMP = nmesy(k)
                else
                  nmescurOMP = 0
                end if
                if((effet.eq.0.or.link.eq.2).or.link.eq.4) then
                  res2(k) = 0.d0
                end if
                if(itcenter.eq.1) then
                  resG1 = resG1 + res2dc(k)-0.5d0*dlog(2.d0*pi)*&
                              nmes_o(k)-0.5d0*dlog(sigmae)*nmes_o(k)
                  resG4 = resG4+ dlog(integrale4(k))
                  resG3 = resG3 +  res2dc(k)-0.5d0*dlog(2.d0*pi)*&
                              nmes_o(k)-0.5d0*dlog(sigmae)*nmes_o(k)
                end if
                if (integrale4(k).le.0.d0) then
                  res_indiv= ((nuzt)*treat_ind(k))*cdc(k)&
                              -112.0d0
                else
                  res_indiv=((nuzt)*treat_ind(k))*cdc(k)&
                             + dlog(integrale4(k))
                endif
                resCcurr= resCcurr*dexp(res_indiv)!*(10.0d0**Cmult)
              end if
            end do 
            resG2 = resG2+resCcurr
            deallocate(Z1OMP,muOMP,ycurrentOMP)
          end do
          !$omp end parallel do
          rescenter(ll) = dlog(resG2/nmcmed) +  resG1 -nindcenter(ll)*log(10.0d0)*Cmult
        end if
        res=res+rescenter(ll)
      end do
      !---------- calcul de la penalisation -------------------
      pe1 = 0.d0
      pe2 = 0.d0
      do i=1,n-3
          if(effet.eq.1) then
            pe1 = pe1+(the1(i-3)*the1(i-3)*m3m3(i))+(the1(i-2) &
            *the1(i-2)*m2m2(i))+(the1(i-1)*the1(i-1)*m1m1(i))+( &
            the1(i)*the1(i)*mmm(i))+(2.d0*the1(i-3)*the1(i-2)* &
            m3m2(i))+(2.d0*the1(i-3)*the1(i-1)*m3m1(i))+(2.d0* &
            the1(i-3)*the1(i)*m3m(i))+(2.d0*the1(i-2)*the1(i-1)* &
            m2m1(i))+(2.d0*the1(i-2)*the1(i)*m2m(i))+(2.d0*the1(i-1) &
            *the1(i)*m1m(i))
          else
            pe1 = 0.d0
          end if
          if(nst.eq.1)then
            pe2=0.d0
          else
            pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2) &
            *the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+( &
            the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)* &
            m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0* &
            the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)* &
            m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1) &
            *the2(i)*m1m(i))
          endif
      end do
      pe = k0(1)*pe1 + k0(2)*pe2
      resnonpen = res
      res = res - pe
      if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
          funcpajLongisplines2=-1.d9
          if(typeJoint.eq.3)    Rrec = 0.d0
          if(typeJoint.eq.3)    Nrec = 0.d0
          Rdc = 0.d0
          Ndc = 0.d0
          goto 123
      else
          funcpajLongisplines2 = res
          do k=1,ng
            if(typeJoint.eq.3)Rrec(k)=res1(k)
            if(typeJoint.eq.3)Nrec(k)=nig(k)
            Rdc(k)=RisqCumul(k)
            Ndc(k)=cdc(k)
          end do
      end if
    end if
  end if
  !Ad:
  123     continue 
  !call cpu_time(tstop)
  return
end function funcpajlongisplines2
