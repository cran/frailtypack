!==========================  DISTANCE   =================================

    subroutine distancessplines(nz1,nz2,b,effet,mt,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)

    
use tailles,only:ndatemax,npmax,NSUJETMAX
    use comon,only:date,zi,t0,t1,c,nt0,nt1,nsujet,nva,ndate, &
    nst,I_hess,H_hess,Hspl_hess,hess

    implicit none

    integer::nz1,nz2,i,j,n,np,k,l,effet,mt
    double precision::x1,x2,h,su,bsup,binf,lam,lbinf, &
    lbsup
    double precision,dimension(npmax,npmax)::hes1,hes2
    double precision,dimension(-2:npmax)::the1,the2
    double precision,dimension(npmax)::b
    double precision,dimension(mt,3)::lamOut,suOut,lam2Out,su2Out
    double precision,dimension(mt)::x1Out,x2Out


    n  = nz1+2
    if(nst.eq.2)then 
        np  = nz1+2+nz2+2+effet+nva
    else
        np  = nz1+2+effet+nva
    endif

    do i=1,nz1+2
        do j=1,nz1+2
            hes1(i,j)=hess(i,j)
        end do
    end do 

    if(nst.eq.2)then
        k = 0
        do i=nz1+3,nz1+2+nz2+2
            k = k + 1 
            l = 0
            do j=nz1+3,nz1+2+nz2+2
                l = l + 1
                hes2(k,l)=hess(i,j)
            end do
        end do
    endif

    do i=1,nz1+2
        the1(i-3)=(b(i))*(b(i))
    end do
    if(nst.eq.2)then
        do i=1,nz2+2
            j = nz1+2+i
            the2(i-3)=(b(j))*(b(j))
        end do
    endif

    h = (zi(n)-zi(1))/(mt-1) ! Al modif : *0.01d0 ! attention depend de mt 
    x1 = zi(1)
    x2 = zi(1)

    do i=1,mt
        
        if(i.ne.1)then
            x1 = x1 + h
        end if
        call cosps(x1,the1,nz1+2,hes1,zi,binf,su,bsup,lbinf,lam,lbsup)
        if(bsup.lt.0.d0)then
            bsup = 0.d0
        endif
        if(binf.gt.1.d0)then
            binf = 1.d0
        endif
        if(lbinf.lt.0.d0)then
            lbinf = 0.d0
        endif

!   Replaced by next sentences and add new ones JRG January 05
!
!
        x1Out(i)=x1
        lamOut(i,1)=lam
        lamOut(i,2)=lbinf
        lamOut(i,3)=lbsup
        suOut(i,1)=su

        suOut(i,2)=binf
        suOut(i,3)=bsup


        if(nst.eq.2)then
            if(i.ne.1)then
                x2 = x2 + h
            endif
            call cosps(x2,the2,nz2+2,hes2,zi,binf,su,bsup,lbinf,lam,lbsup)
            if(bsup.lt.0.d0)then
                bsup = 0.d0
            endif
            if(binf.gt.1.d0)then
                binf = 1.d0
            endif
            if(lbinf.lt.0.d0)then
                lbinf = 0.d0
            endif

            x2Out(i)=x2
            lam2Out(i,1)=lam
            lam2Out(i,2)=lbinf
            lam2Out(i,3)=lbsup
            su2Out(i,1)=su
            su2Out(i,2)=binf
            su2Out(i,3)=bsup
        endif
    end do

    return

    end subroutine distancessplines

!=========================================================================================================
!=========================================================================================================
!=========================================================================================================
!=========================================================================================================

    subroutine distancecpm(b,m,mt,xR1,moyLamR1,xSu1,moysuR1,xR2,moyLamR2,xSu2,moysuR2)

    use tailles
    use comon,only:cens,vvv,t0,t1,c,nsujet,nva,nva1,nva2,nst,nbrecu,nbintervR,ttt,date
    use optim

    implicit none

    integer::m,i,j,ier,t,k,gg,jj,uu,ii,ncur,l,mt
    double precision::sx,som,som11R,som21R,som11RW,som21RW, &
    glRW,x
    double precision,dimension(m)::b
    double precision,dimension(m)::bgen
    double precision,dimension(1000,m)::u,v
    double precision,dimension(1000)::lamR,suR,glR    
    double precision,dimension((m*(m+1)/2))::vv
    double precision,dimension(nbintervR)::tempsR

    double precision::lamR25,suR25,lamR975,suR975, &
    LamRW,suRW
!AD: sorties
    double precision::ep
    double precision,dimension(mt)::xR1,xR2
    double precision,dimension(mt,3)::moyLamR1,moyLamR2
    double precision,dimension(100)::xSu1,xSu2
    double precision,dimension(100,3)::moysuR1,moysuR2


    sx=1.d0 ! ecart-type ou variance des réalisations gaussiennes générées
    uu=0
    do i=1,nbintervR
        tempsR(i)=(ttt(i-1)+ttt(i))/(2.d0)
    end do

    do i=1,m*(m+1)/2
        vv(i)=vvv(i)
    end do

!AD: ep=10.d-10
    ep=10.d-10
    call dmfsdj(vv,m,ep,ier)

!      ns=m-nva-3! ns est le nombre de paramètres intervenat dans la fonctionde
! risque ici 3 car on estime p param : alpha,eta,theta
!     Pour chaque temps de 0 à la censure du décès par pas cens/100
!      x=-cens/100 !permet de commencer à 0 dans la boucle

!      n=0 ! permet de commencer à 0 dans la boucle


    ncur = 1

!!!!!!!!    HAZ     !!!!!!!!!!
!!!!!!!!
    do t=1,nbintervR
        lamR=0.d0
        glR=0.d0

!         n=n+1 ! compteur sur nième temps
        x=tempsR(t)
! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Pour chaque paramètre estimé
            do i=1,nbintervR
                som=0.d0
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    call bgos(SX,0,u(k,j),v(k,j),0.d0)
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do
                bgen(i)=(b(i)+som)**2
            end do              ! en sortie on récupère le nouveau vecteur b


!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR
!     fct de survie           : su  = dexp(-gl)

            lamR(k)=0.d0
            lamRW=0.d0

            do gg=1,nbintervR
                if ((x.ge.(ttt(gg-1))).and.(x.lt.(ttt(gg)))) then
                    lamR(k)=bgen(gg)
                endif
            end do

            do l=1,nbintervR
                if ((x.ge.(ttt(l-1))).and.(x.lt.(ttt(l)))) then
                    lamRW=b(l)**2
                endif
            end do
            

        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles

        call percentile(lamR,lamR25,lamR975)


        do ii=1,nbintervR
            if (x.gt.ttt(ii-1).and.x.lt.ttt(ii))then
                uu=ii
            endif
        end do
        if(t.ne.1) then 
            xR1(ncur) = real(ttt(uu-1))
        else
            xR1=date(1)
        end if
        xR1(ncur+1) = real(x)
        xR1(ncur+2) = real(ttt(uu))

        moyLamR1(ncur,1) = real(LamRW)
        moyLamR1(ncur+1,1) = moyLamR1(ncur,1)
        moyLamR1(ncur+2,1) = moyLamR1(ncur,1)

        moyLamR1(ncur,2) = real(lamR25)
        moyLamR1(ncur+1,2) = moyLamR1(ncur,2)
        moyLamR1(ncur+2,2) = moyLamR1(ncur,2)


        moyLamR1(ncur,3) = real(lamR975)
        moyLamR1(ncur+1,3) = moyLamR1(ncur,3)
        moyLamR1(ncur+2,3) = moyLamR1(ncur,3)


        ncur=ncur+3

    end do


!!!! SURV 1

    bgen=0.d0
    x=ttt(0) !date(1)

    do t=1,100
        if(t .ne.1) then
            x=x+(ttt(nbintervR)-ttt(0))/99 !100 changement 04/07/2013 by Al
        end if

        do k=1,1000
            do i=1,nbintervR
                som=0.d0
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    call bgos(SX,0,u(k,j),v(k,j),0.d0)
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do
                bgen(i)=(b(i)+som)**2

            end do
!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR
!     fct de survie           : su  = dexp(-gl)
            som11R=0.d0
            som21R=0.d0
            glR(k)=0.d0
            suR(k)=0.d0

            do gg=1,nbintervR
                if ((x.ge.(ttt(gg-1))).and.(x.lt.(ttt(gg)))) then

                    som11R=bgen(gg)*(x-ttt(gg-1))

                    if (gg.ge.2)then
                        do jj=1,gg-1
                            som21R=som21R+bgen(jj)*(ttt(jj)-ttt(jj-1))
                        end do
                    endif

                    glR(k)=(som11R+som21R)
                    suR(k)=dexp(-glR(k))
                endif
            end do
            
        end do

        som11RW=0.d0
        som21RW=0.d0
        glRW=0.d0
        suRW=0.d0
        do l=1,nbintervR
            if ((x.ge.(ttt(l-1))).and.(x.lt.(ttt(l)))) then

                som11RW=(b(l)**2)*(x-ttt(l-1))

                if (l.ge.2)then
                    do jj=1,l-1
                    som21RW=som21RW+(b(jj)**2)*(ttt(jj)-ttt(jj-1))
                    end do
                endif

                glRW=(som11RW+som21RW)
                suRW=dexp(-glRW)
            endif
        end do

        suR25=0.d0
        suR975=0.d0
        call percentile(suR,suR25,suR975)
        if(t.ne.1) then
            xSu1(t) = real(x)
        else
            xSu1(t) = ttt(0) !date(1)
        end if


        moysuR1(t,1) = suRW

        moysuR1(t,2) = suR25

        moysuR1(t,3) = suR975

        if(moysuR1(t,1).lt.0.d0)then
            moysuR1(t,1) = 0.d0
        endif
        if(moysuR1(t,1).gt.1.d0)then
            moysuR1(t,1) = 1.d0
        endif

        if(moysuR1(t,2).lt.0.d0)then
            moysuR1(t,2) = 0.d0
        endif
        if(moysuR1(t,2).gt.1.d0)then
            moysuR1(t,2) = 1.d0
        endif
        if(moysuR1(t,3).lt.0.d0)then
            moysuR1(t,3) = 0.d0
        endif
        if(moysuR1(t,3).gt.1.d0)then
            moysuR1(t,3) = 1.d0
        endif
    end do
!!!!!!!   FIN SURV     !!!!!!!!!!

!!!! HAZ 2
!!!!


    if (nst == 2) then

        ncur = 1

        do t=1,nbintervR
            lamR=0.d0
            glR=0.d0
!         n=n+1 ! compteur sur nième temps
            x=tempsR(t)
! On simule 1000 réalisations gaussienne par paramètres
            do k=1,1000
!     Pour chaque paramètre estimé
                do i=nbintervR+1,nbintervR*nst
                    som=0.d0
                    do j=1,i         ! cela correspond au produit trp(L)%*%U
                        call bgos(SX,0,u(k,j),v(k,j),0.d0)
                        som=som+vv(i*(i-1)/2+j)*u(k,j)
                    end do

                    bgen(i-nbintervR)=(b(i)+som)**2

                end do              ! en sortie on récupère le nouveau vecteur b
                      
!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR
!     fct de survie           : su  = dexp(-gl)

                lamR(k)=0.d0
                lamRW=0.d0

                do gg=1,nbintervR
                    if ((x.ge.(ttt(gg-1))).and.(x.lt.(ttt(gg)))) then
                        lamR(k)=bgen(gg)
                    endif
                end do

                do l=1,nbintervR
                    if ((x.ge.(ttt(l-1))).and.(x.lt.(ttt(l)))) then
                        lamRW=b(l+nbintervR)**2
                    endif
                end do


            end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles

            call percentile(lamR,lamR25,lamR975)


            do ii=1,nbintervR
                if (x.gt.ttt(ii-1).and.x.lt.ttt(ii))then
                    uu=ii
                endif
            end do
            if(t.ne.1) then 
                xR2(ncur) = real(ttt(uu-1))
            else
                xR2=date(1)
            end if
            xR2(ncur+1) = real(x)
            xR2(ncur+2) = real(ttt(uu))

            moyLamR2(ncur,1) = real(LamRW)
            moyLamR2(ncur+1,1) = moyLamR2(ncur,1)
            moyLamR2(ncur+2,1) = moyLamR2(ncur,1)

            moyLamR2(ncur,2) = real(lamR25)
            moyLamR2(ncur+1,2) = moyLamR2(ncur,2)
            moyLamR2(ncur+2,2) = moyLamR2(ncur,2)

            moyLamR2(ncur,3) = real(lamR975)
            moyLamR2(ncur+1,3) = moyLamR2(ncur,3)
            moyLamR2(ncur+2,3) = moyLamR2(ncur,3)
            ncur=ncur+3
        end do

        bgen=0.d0

        x=ttt(0) !date(1)
        do t=1,100
            if(t .ne. 1) then
                x=x+(ttt(nbintervR)-ttt(0))/99 !100 changement 04/07/2013 by Al
            end if

            do k=1,1000
                do i=nbintervR+1,nbintervR*nst
                    som=0.d0
                    do j=1,i         ! cela correspond au produit trp(L)%*%U
                        call bgos(SX,0,u(k,j),v(k,j),0.d0)
                        som=som+vv(i*(i-1)/2+j)*u(k,j)
                    end do
                    if(bgen(i-nbintervR).gt.(1.d0)) then
                        bgen(i-nbintervR)=1.d0
                    else
                        bgen(i-nbintervR)=(b(i)+som)**2
                    endif
                end do

                som11R=0.d0
                som21R=0.d0
                glR(k)=0.d0
                suR(k)=0.d0

                do gg=1,nbintervR
                    if ((x.ge.(ttt(gg-1))).and.(x.lt.(ttt(gg)))) then

                        som11R=bgen(gg)*(x-ttt(gg-1))

                        if (gg.ge.2)then
                            do jj=1,gg-1
                                som21R=som21R+bgen(jj)*(ttt(jj)-ttt(jj-1))
                            end do
                        endif

                        glR(k)=(som11R+som21R)
                        suR(k)=dexp(-glR(k))
                    endif
                end do
            end do

            som11RW=0.d0
            som21RW=0.d0
            glRW=0.d0
            suRW=0.d0

            do l=1,nbintervR
                if ((x.ge.(ttt(l-1))).and.(x.lt.(ttt(l)))) then

                    som11RW=(b(l+nbintervR)**2)*(x-ttt(l-1))

                    if (l.ge.2)then
                        do jj=1,l-1
                            som21RW=som21RW+(b(jj+nbintervR)**2)*(ttt(jj)-ttt(jj-1))
                        end do
                    endif
                    glRW=(som11RW+som21RW)
                    suRW=dexp(-glRW)
                endif
            end do
            suR25=0.d0
            suR975=0.d0
            call percentile(suR,suR25,suR975)
            if(t.ne.1) then
                xSu2(t) = real(x)
            else
                xSu2(t) = ttt(0) !date(1)
            end if
            moysuR2(t,1) = suRW
            moysuR2(t,2) = suR25
            moysuR2(t,3) = suR975

            if(moysuR2(t,1).lt.0.d0)then
                moysuR2(t,1) = 0.d0
            endif
            if(moysuR2(t,1).gt.1.d0)then
                moysuR2(t,1) = 1.d0
            endif

            if(moysuR2(t,2).lt.0.d0)then
                moysuR2(t,2) = 0.d0
            endif
            if(moysuR2(t,2).gt.1.d0)then
                moysuR2(t,2) = 1.d0
            endif
            if(moysuR2(t,3).lt.0.d0)then
                moysuR2(t,3) = 0.d0
            endif
            if(moysuR2(t,3).gt.1.d0)then
                moysuR2(t,3) = 1.d0
            endif
        end do
    end if

    end subroutine distancecpm

!=========================================================================================================
!=========================================================================================================
!=========================================================================================================
!=========================================================================================================

    subroutine distanceweib(b,m,mt,xR1,moyLamR1,xSu1,moysuR1,xR2,moyLamR2,xSu2,moysuR2)
    
    use tailles
    use comon,only:cens,vvv,t0,t1,c,nsujet &
    ,nva,nva1,nva2,nst,typeof2,etaR,etaD,betaR,betaD,date,mint
    use optim
    
    implicit none
    
    integer,intent(in)::m,mt
    integer::i,j,ier,t,k,ns
    double precision::lamR25,lamDC25,suR25,suDC25,lamR975,lamDC975,suR975,suDC975, &
    LamRW,LamDCW,suRW,suDCW
        double precision,dimension(m*(m+1)/2)::vv
    double precision,dimension(m)::b,bgen
    double precision,dimension(1000,m)::u
    double precision::sx,som,x,zz,zy 
    double precision,dimension(1000):: lamR,lamDC,suR,suDC,glDC,glR
! theorique - estimés après Max de Vrais
      double precision::glDCW,glRW

    double precision,dimension(mt)::xR1,xR2    
    double precision,dimension(mt,3)::moyLamR1,moyLamR2
    double precision,dimension(100,3)::moysuR1,moysuR2
    double precision,dimension(100)::xSu1,xSu2
    double precision::ep




    lamR25=0.d0
    lamDC25=0.d0
    suR25=0.d0
    suDC25=0.d0
    lamR975=0.d0
    lamDC975=0.d0
    suR975=0.d0
    suDC975=0.d0
    LamRW=0.d0
    LamDCW=0.d0
    suRW=0.d0
    suDCW=0.d0
    ier=0
    t=0
    k=0
    ns=0
    sx=1.d0 ! ecart-type ou variance des réalisations gaussiennes générées
    !typeof2 = 1 shared weib
    if (typeof2 == 1) then
        ns=m-nva-1
    end if
    !typeof2 = 2  weib
     if (typeof2 == 2) then
        ns=m-nva-2
    end if
    
    do i=1,ns*(ns+1)/2
        vv(i)=vvv(i)
    end do

    ep=10.d-10
    call dmfsdj(vv,m,ep,ier)

    do k=1,1000!0
        do j=1,ns
            call bgos(SX,0,zz,zy,0.d0)
            u(k,j)=zz
        end do
    end do


!     Pour chaque temps de 0 à la censure du décès par pas cens/100

 !     commencer à 0 dans la boucle
    x=mint !date(1)
!       shape et scale

    betaR = b(1)**2
    etaR =  b(2)**2
    if (nst == 2) then
        betaD = b(3)**2
        etaD = b(4)**2
    else
        etaD = 0.d0
        betaD = 0.d0
    end if

    do t=1,mt

        lamR=0.d0
        lamDC=0.d0
!         n=n+1 ! compteur sur nième temps
        if(t.ne.1) then
            x=x+(cens-mint)/(mt-1)
        end if

! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Pour chaque paramètre estimé
            do i=1,ns
                som=0.d0
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do
                bgen(i)=(b(i)+som)**2
            end do              ! en sortie on récupère le nouveau vecteur b



!            Moybgen(k)=Moybgen(k)+bgen(i)/1000
!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR
!     fct de survie           : su  = dexp(-gl)

            lamR(k)=(bgen(1)*(x**(bgen(1)-1.d0)))/(bgen(2)**bgen(1))
            lamRW=((b(1)**2)*(x**((b(1)**2)-1.d0)))/((b(2)**2)**(b(1)**2))

            if(nst == 2) then
                lamDC(k)=(bgen(3)*(x**(bgen(3)-1.d0)))/(bgen(4)**bgen(3))
                lamDCW=((b(3)**2)*(x**((b(3)**2)-1.d0)))/((b(4)**2)**(b(3)**2))
            end if
        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles

        call percentile(lamR,lamR25,lamR975)
        if(nst == 2) then
            call percentile(lamDC,lamDC25,lamDC975)
        end if
!----- strate 1
        if(t == 1) then
            xR1(t)=mint !date(1)
        else
            xR1(t) = real(x)
        end if

        moyLamR1(t,1) = real(LamRW)
        moyLamR1(t,2) = real(lamR25)
        moyLamR1(t,3) = real(lamR975)

        if(nst == 2) then
!----- strate 2
            xR2(t) = xR1(t)!real(x)
            moyLamR2(t,1) = real(lamDCW)
            moyLamR2(t,2) = real(lamDC25)
            moyLamR2(t,3) = real(lamDC975)
        else
            xR2(t) = 0.d0
            moyLamR2(t,1) =  0.d0
            moyLamR2(t,2) =  0.d0
            moyLamR2(t,3) =  0.d0
        end if

    end do

    x=mint !date(1)

    do t=1,100

        glR=0.d0
        suR=0.d0
        glDC=0.d0
        suDC=0.d0
!         n=n+1 ! compteur sur nième temps
        if(t.ne.1) then
            x=x+(cens-mint)/99 !100
        end if
! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Pour chaque paramètre estimé
            do i=1,ns
                som=0.d0
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do
                bgen(i)=(b(i)+som)**2
            end do              ! en sortie on récupère le nouveau vecteur b



!            Moybgen(k)=Moybgen(k)+bgen(i)/1000
!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR
!     fct de survie           : su  = dexp(-gl)

            glR(k)  =  (x/bgen(2))**bgen(1)
            suR(k)  = dexp(-glR(k))
            glRW  =  (x/(b(2)**2))**(b(1)**2)
            suRW  = dexp(-glRW)

            if(nst == 2) then
                glDC(k)  =  (x/bgen(4))**bgen(3)
                suDC(k)  = dexp(-glDC(k))
                glDCW  =  (x/(b(4)**2))**(b(3)**2)
                suDCW  = dexp(-glDCW)
            end if

        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles

        call percentile(suR,suR25,suR975)
        if(nst == 2) then
            call percentile(suDC,suDC25,suDC975)
        end if

        if(t==1) then
            xSu1(t) = mint !date(1)
        else
            xSu1(t) = real(x)
        end if
        moysuR1(t,1) = real(suRW)
        moysuR1(t,2) = real(suR25)
        moysuR1(t,3) = real(suR975)

        if(moysuR1(t,1).lt.0.d0)then
            moysuR1(t,1) = 0.d0
        endif
        if(moysuR1(t,1).gt.1.d0)then
            moysuR1(t,1) = 1.d0
        endif

        if(moysuR1(t,2).lt.0.d0)then
            moysuR1(t,2) = 0.d0
        endif
        if(moysuR1(t,2).gt.1.d0)then
            moysuR1(t,2) = 1.d0
        endif
        if(moysuR1(t,3).lt.0.d0)then
            moysuR1(t,3) = 0.d0
        endif
        if(moysuR1(t,3).gt.1.d0)then
            moysuR1(t,3) = 1.d0
        endif

        if(nst == 2) then
            xSu2(t) = xSu1(t)!real(x)
            moysuR2(t,1) = real(suDCW)
            moysuR2(t,2) = real(suDC25)
            moysuR2(t,3) = real(suDC975)
            if(moysuR2(t,1).lt.0.d0)then
                moysuR2(t,1) = 0.d0
            endif
            if(moysuR2(t,1).gt.1.d0)then
                moysuR2(t,1) = 1.d0
            endif

            if(moysuR2(t,2).lt.0.d0)then
                moysuR2(t,2) = 0.d0
            endif
            if(moysuR2(t,2).gt.1.d0)then
                moysuR2(t,2) = 1.d0
            endif    
            if(moysuR2(t,3).lt.0.d0)then
                moysuR2(t,3) = 0.d0
            endif
            if(moysuR2(t,3).gt.1.d0)then
                moysuR2(t,3) = 1.d0
            endif
        else
            xSu2(t) = 0.d0
            moysuR2(t,1) = 0.d0
            moysuR2(t,2) = 0.d0
            moysuR2(t,3) = 0.d0
        end if
    end do


    end subroutine distanceweib

!==========================  DISTANCE   =================================

    subroutine distanceJsplines(nz1,nz2,b,mt1,mt2,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)
    
    use tailles
    use comon,only:date,datedc,zi,t0,t1,t0dc,t1dc,c,cdc &
    ,nt0,nt1,nt0dc,nt1dc,nsujet,nva,nva1,nva2,ndate,ndatedc,nst &
    ,PEN_deri,I_hess,H_hess,Hspl_hess,hess
    
    Implicit none  
    
    
    integer,intent(in):: nz1,nz2,mt1,mt2
    double precision ,dimension(npmax),intent(in):: b
    integer::i,j,n,k,l      
    double precision::x1,x2,h,su,bsup,binf,lam,lbinf,lbsup
    double precision ,dimension(npmax,npmax)::hes1,hes2
    double precision ,dimension(-2:npmax):: the1,the2     
    double precision,dimension(mt1)::x1Out
    double precision,dimension(mt2)::x2Out
    double precision,dimension(mt1,3):: lamOut,suOut
    double precision,dimension(mt2,3):: lam2Out,su2Out
        
    n  = nz1+2
    do i=1,nz1+2
        do j=1,nz1+2
            hes1(i,j)=h_Hess(i,j)
        end do
    end do 

    if(nst.eq.2)then  
        k = 0
        do i=nz1+3,nz1+2+nz2+2
            k = k + 1 
            l = 0
            do j=nz1+3,nz1+2+nz2+2
                l = l + 1
                hes2(k,l)=H_hess(i,j)
            end do
        end do   
    endif

    do i=1,nz1+2
        the1(i-3)=(b(i))*(b(i))
    end do
    if(nst.eq.2)then  
        do i=1,nz2+2
            j = nz1+2+i
            the2(i-3)=(b(j))*(b(j))
        end do
    endif
    h = (zi(n)-zi(1))*0.01d0
    x1 = zi(1)
    x2 = zi(1)     
  
! Recurrent          
    do i=1,mt1
        if(i .ne.1)then
            x1 = x1 + h 
        end if
        call cospJ(x1,the1,nz1+2,hes1,zi,binf,su,bsup,lbinf,lam,lbsup)
                
        if(bsup.lt.0.d0)then
            bsup = 0.d0
        endif
        if(binf.gt.1.d0)then
            binf = 1.d0 
        endif
        if(lbinf.lt.0.d0)then
            lbinf = 0.d0
        endif 
!!!   Replaced by next sentences and add new ones JRG January 05

        x1Out(i)=x1
        lamOut(i,1)=lam
        lamOut(i,2)=lbinf
        lamOut(i,3)=lbsup 
        suOut(i,1)=su
        suOut(i,2)=bsup
        suOut(i,3)=binf
    end do   

! Death  
    if(nst.eq.2)then        
        do i=1,mt2
!!!   Replaced by next sentences and add new ones JRG January 05    
            if(i.ne.1)then
                x2 = x2 + h 
            endif    
            call cospJ(x2,the2,nz2+2,hes2,zi,binf,su,bsup,lbinf,lam,lbsup)
            if(bsup.lt.0.d0)then
                bsup = 0.d0
            endif
            if(binf.gt.1.d0)then
                binf = 1.d0 
            endif
            if(lbinf.lt.0.d0)then
                lbinf = 0.d0
            endif 

            x2Out(i)=x2
            lam2Out(i,1)=lam
            lam2Out(i,2)=lbinf
            lam2Out(i,3)=lbsup 
            su2Out(i,1)=su
            su2Out(i,2)=binf
            su2Out(i,3)=bsup 
        end do  
    endif
              
    return
        
    end subroutine distanceJsplines

!====================================================================

    subroutine distancejcpm(b,m,mt1,mt2,x1R,moyLamR,xSu1,moysuR,x2DC,moyLamDC,xSu2,moysuDC)

    use tailles
    use comon,only:cens,vvv,t0,t1,t0dc,t1dc,c,cdc,nsujet, &
    nva,nva1,nva2,nst,nbintervR,nbintervDC,ttt,tttDC,date
    use optim

    implicit none

    integer::m,i,j,ier,t,k,gg,jj,uu,ii,ncur,mt1,mt2,l
    double precision::sx,som,som11R,som21R,som11RW,som21RW, &
    glRW,x
    double precision,dimension(m)::b
    double precision,dimension(m)::bgen
    double precision,dimension(1000,m)::u,v
    double precision,dimension(1000)::lamR,lamDC,suR,glR    
    double precision,dimension((m*(m+1)/2))::vv
    double precision,dimension(m)::tempsR,tempsDC

    double precision::lamR25,lamDC25,suR25,lamR975,lamDC975,suR975, &
    LamRW,LamDCW,suRW
!AD: sorties
    double precision::ep
    double precision,dimension(mt1)::x1R    !0:nbintervR*3
    double precision,dimension(mt1,3)::moyLamR
    double precision,dimension(100,3)::moysuR    
    double precision,dimension(mt2)::x2DC !0:nbintervDC*3
    double precision,dimension(mt2,3)::moyLamDC
    double precision,dimension(100,3)::moysuDC
    double precision,dimension(100)::xSu1,xSu2


    uu=0
    bgen=0.d0
    sx=1.d0 ! ecart-type ou variance des réalisations gaussiennes générées

    do i=1,nbintervR
        tempsR(i)=(ttt(i-1)+ttt(i))/(2.d0)
    end do

    do i=1,nbintervDC
        tempsDC(i)=(tttdc(i-1)+tttdc(i))/(2.d0)
    end do

    do i=1,m*(m+1)/2
        vv(i)=vvv(i)
    end do
!AD: ep=10.d-10
    ep=10.d-10
    call dmfsdj(vv,m,ep,ier)

!      ns=m-nva-3! ns est le nombre de paramètres intervenat dans la fonctionde
! risque ici 3 car on estime p param : alpha,eta,theta
!     Pour chaque temps de 0 à la censure du décès par pas cens/100
!      x=-cens/100 !permet de commencer à 0 dans la boucle

!      n=0 ! permet de commencer à 0 dans la boucle
    ncur = 1
    do t=1,nbintervR
        lamR=0.d0
!         n=n+1 ! compteur sur nième temps
        x=tempsR(t)
! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Pour chaque paramètre estimé
            do i=1,nbintervR!m
                som=0.d0
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    call bgos(SX,0,u(k,j),v(k,j),0.d0)
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do
                bgen(i)=(b(i)+som)**2
            end do              ! en sortie on récupère le nouveau vecteur b

!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR
!     fct de survie           : su  = dexp(-gl)

            lamR(k)=0.d0
            lamRW=0.d0

            do gg=1,nbintervR
                if ((x.ge.(ttt(gg-1))).and.(x.lt.(ttt(gg)))) then
                    lamR(k)=bgen(gg)
                endif
            end do

            do gg=1,nbintervR
                if ((x.ge.(ttt(gg-1))).and.(x.lt.(ttt(gg)))) then
                    lamRW=b(gg)**2
                endif
            end do

        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles

        call percentile(lamR,lamR25,lamR975)

        do ii=1,nbintervR
            if (x.gt.ttt(ii-1).and.x.lt.ttt(ii))then
                uu=ii
            endif
        end do
        if(t.ne.1) then 
            x1R(ncur) = real(ttt(uu-1))
        else
            x1R=date(1)
        end if
        x1R(ncur+1) = real(x)
        x1R(ncur+2) = real(ttt(uu))

        moyLamR(ncur,1) = real(LamRW)
        moyLamR(ncur+1,1) = moyLamR(ncur,1)
        moyLamR(ncur+2,1) = moyLamR(ncur,1)

        moyLamR(ncur,2) = real(lamR25)
        moyLamR(ncur+1,2) = moyLamR(ncur,2)
        moyLamR(ncur+2,2) = moyLamR(ncur,2)


        moyLamR(ncur,3) = real(lamR975)
        moyLamR(ncur+1,3) = moyLamR(ncur,3)
        moyLamR(ncur+2,3) = moyLamR(ncur,3)

        ncur=ncur+3
    end do


    ncur = 1
    do t=1,nbintervDC
        x=tempsDC(t)
        lamDC=0.d0

!         n=n+1 ! compteur sur nième temps
! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Pour chaque paramètre estimé
            do i=1,nbintervDC!m

                som=0.d0

                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    call bgos(SX,0,u(k,j),v(k,j),0.d0)
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do

                bgen(i)=(b(nbintervR + i)+som)**2
            end do              ! en sortie on récupère le nouveau vecteur b


!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR
!     fct de survie           : su  = dexp(-gl)

            lamDC(k)=0.d0
            lamDCW=0.d0

            do gg=1,nbintervDC
                if((x.ge.(tttdc(gg-1))).and.(x.lt.(tttdc(gg))))then
                    lamDC(k)=bgen(gg)
                endif
            end do

            do gg=1,nbintervDC
                if((x.ge.(tttdc(gg-1))).and.(x.lt.(tttdc(gg))))then
                    lamDCW=b(gg+nbintervR)**2
                endif
            end do
        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles

        call percentile(lamDC,lamDC25,lamDC975)


        do ii=1,nbintervDC
            if (x.gt.tttdc(ii-1).and.x.lt.tttdc(ii))then
                uu=ii
            endif
        end do
        if(t.ne.1) then 
            x2DC(ncur) = real(tttdc(uu-1))
        else
            x2DC(ncur) = date(1)
        end if
        x2DC(ncur+1) = real(x)
        x2DC(ncur+2) = real(tttdc(uu))

        moyLamDC(ncur,1) = real(lamDCW)
        moyLamDC(ncur+1,1) = moyLamDC(ncur,1)
        moyLamDC(ncur+2,1) = moyLamDC(ncur,1)

        moyLamDC(ncur,2) = real(lamDC25)
        moyLamDC(ncur+1,2) = moyLamDC(ncur,2)
        moyLamDC(ncur+2,2) = moyLamDC(ncur,2)

        moyLamDC(ncur,3) = real(lamDC975)
        moyLamDC(ncur+1,3) = moyLamDC(ncur,3)
        moyLamDC(ncur+2,3) = moyLamDC(ncur,3)

        ncur=ncur+3
    end do

!--------------- Fontion de survie recurrent ----------------
    bgen=0.d0

    x=date(1)

    do t=1,100

!         n=n+1 ! compteur sur nième temps
        if(t.ne.1) then
            x=x+ttt(nbintervR)/100
        end if
! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Pour chaque paramètre estimé
            do i=1,nbintervR!m
                som=0.d0
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    call bgos(SX,0,u(k,j),v(k,j),0.d0)
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do
                bgen(i)=(b(i)+som)**2
            end do              ! en sortie on récupère le nouveau vecteur b

!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR
!     fct de survie           : su  = dexp(-gl)

            som11R=0.d0
            som21R=0.d0
            glR(k)=0.d0
            suR(k)=0.d0

            do gg=1,nbintervR
                if ((x.ge.(ttt(gg-1))).and.(x.lt.(ttt(gg)))) then

                    som11R=bgen(gg)*(x-ttt(gg-1))

                    if (gg.ge.2)then
                        do jj=1,gg-1
                            som21R=som21R+bgen(jj)*(ttt(jj)-ttt(jj-1))
                        end do
                    endif

                    glR(k)=(som11R+som21R)
                    suR(k)=dexp(-glR(k))
                endif
            end do

        end do

        som11RW=0.d0
        som21RW=0.d0
        glRW=0.d0
        suRW=0.d0
        do l=1,nbintervR
            if ((x.ge.(ttt(l-1))).and.(x.lt.(ttt(l)))) then

                som11RW=(b(l)**2)*(x-ttt(l-1))

                if (l.ge.2)then
                    do jj=1,l-1
                        som21RW=som21RW+(b(jj)**2)*(ttt(jj)-ttt(jj-1))
                    end do
                endif

                glRW=(som11RW+som21RW)
                suRW=dexp(-glRW)
            endif
        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles
        suR25=0.d0
        suR975=0.d0
        call percentile(suR,suR25,suR975)
        if(t.ne.1) then
            xSu1(t) = real(x)
        else
            xSu1(t) = date(1)
        end if
        moysuR(t,1) = suRW

        moysuR(t,2) = suR25

        moysuR(t,3) = suR975

        if(moysuR(t,1).lt.0.d0)then
            moysuR(t,1) = 0.d0
        endif
        if(moysuR(t,1).gt.1.d0)then
            moysuR(t,1) = 1.d0
        endif

        if(moysuR(t,2).lt.0.d0)then
            moysuR(t,2) = 0.d0
        endif
        if(moysuR(t,2).gt.1.d0)then
            moysuR(t,2) = 1.d0
        endif
        if(moysuR(t,3).lt.0.d0)then
            moysuR(t,3) = 0.d0
        endif
        if(moysuR(t,3).gt.1.d0)then
            moysuR(t,3) = 1.d0
        endif

    end do

!--------------- Fontion de survie ----------------
    bgen=0.d0

    x=date(1)
    do t=1,100

!         n=n+1 ! compteur sur nième temps
        if(t.ne.1) then
            x=x+tttdc(nbintervDC)/100
        end if

! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Pour chaque paramètre estimé
            do i=1,nbintervDC!m
                som=0.d0
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    call bgos(SX,0,u(k,j),v(k,j),0.d0)
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do
                bgen(i)=(b(i+nbintervR)+som)**2
            end do              ! en sortie on récupère le nouveau vecteur b

!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR
!     fct de survie           : su  = dexp(-gl)

            som11R=0.d0
            som21R=0.d0
            glR(k)=0.d0
            suR(k)=0.d0


            do gg=1,nbintervDC
                if ((x.ge.(tttdc(gg-1))).and.(x.lt.(tttdc(gg)))) then

                    som11R=bgen(gg)*(x-tttdc(gg-1))

                    if (gg.ge.2)then
                        do jj=1,gg-1
                            som21R=som21R+bgen(jj)*(tttdc(jj)-tttdc(jj-1))
                        end do
                    endif

                    glR(k)=(som11R+som21R)
                    suR(k)=dexp(-glR(k))
                endif
            end do


!            moysuR0=moysuR0+suR(k)/1000
            
        end do
    
        som11RW=0.d0
        som21RW=0.d0
        glRW=0.d0
        suRW=0.d0
        do l=1,nbintervDC
            if ((x.ge.(tttdc(l-1))).and.(x.lt.(tttdc(l)))) then
            
                som11RW=(b(l+nbintervR)**2)*(x-tttdc(l-1))
                
                if (l.ge.2)then
                    do jj=1,l-1
                        som21RW=som21RW+(b(jj+nbintervR)**2)*(tttdc(jj)-tttdc(jj-1))
                    end do
                endif

                glRW=(som11RW+som21RW)
                suRW=dexp(-glRW)
            endif
        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles
        suR25=0.d0
        suR975=0.d0
        call percentile(suR,suR25,suR975)
        if(t.ne.1) then
            xSu2(t) = real(x)
        else
            xSu2(t) = date(1)
        end if

        moysuDC(t,1) = suRW

        moysuDC(t,2) = suR25

        moysuDC(t,3) = suR975

        if(moysuDC(t,1).lt.0.d0)then
            moysuDC(t,1) = 0.d0
        endif
        if(moysuDC(t,1).gt.1.d0)then
            moysuDC(t,1) = 1.d0
        endif

        if(moysuDC(t,2).lt.0.d0)then
            moysuDC(t,2) = 0.d0
        endif
        if(moysuDC(t,2).gt.1.d0)then
            moysuDC(t,2) = 1.d0
        endif
        if(moysuDC(t,3).lt.0.d0)then
            moysuDC(t,3) = 0.d0
        endif
        if(moysuDC(t,3).gt.1.d0)then
            moysuDC(t,3) = 1.d0
        endif

    end do

    end subroutine distancejcpm


    subroutine distanceJweib(b,m,mt1,x1R,moyLamR,xSu1,moysuR,x2DC,moyLamDC,xSu2,moysuDC)
    
    use tailles
    use comon,only:cens,vvv,t0,t1,t0dc,t1dc,c,cdc,nsujet &
    ,nva,nva1,nva2,nst,etaR,etaD,betaR,betaD,date,datedc
    use optim

    implicit none

    integer::m,i,j,ier,t,k,ns,mt1
    double precision::lamR25,lamDC25,suR25,suDC25,lamR975,lamDC975,suR975,suDC975, &
    LamRW,LamDCW,suRW,suDCW
        double precision,dimension(m*(m+1)/2)::vv
    double precision,dimension(m)::b,bgen
    double precision,dimension(1000,m)::u
    double precision::sx,som,x,zz,zy
    double precision,dimension(1000):: lamR,lamDC,suR,suDC,glDC,glR

! theorique - estimés après Max de Vrais
      double precision::glDCW,glRW

    double precision,dimension(mt1)::x1R,x2DC
    double precision,dimension(mt1,3)::moyLamR,moyLamDC

    double precision,dimension(100,3)::moysuR,moysuDC    
    double precision,dimension(100)::xSu1,xSu2
    double precision::ep


    lamR25=0.d0
    lamDC25=0.d0
    suR25=0.d0
    suDC25=0.d0
    lamR975=0.d0
    lamDC975=0.d0
    suR975=0.d0
    suDC975=0.d0
    LamRW=0.d0
    LamDCW=0.d0
    suRW=0.d0
    suDCW=0.d0
    sx=1.d0 ! ecart-type ou variance des réalisations gaussiennes générées
    ns=m-nva-2

    do i=1,ns*(ns+1)/2
        vv(i)=vvv(i)
    end do

    ep=10.d-10
    call dmfsdj(vv,m,ep,ier)

    do k=1,1000
        do j=1,ns
            call bgos(SX,0,zz,zy,0.d0)
            u(k,j)=zz
        end do
    end do

    betaR = b(1)**2
    etaR =  b(2)**2

    if (nst == 2) then
        betaD = b(3)**2
        etaD = b(4)**2
    else
        etaD = 0.d0
        betaD = 0.d0
    end if


!============================ Surv strat 1
!------------------- Pour les recurents
!     Pour chaque temps de 0 à la censure du décès par pas cens/100
!    x=-cens/100 !permet de commencer à 0 dans la boucle

    x=date(1)
    do t=1,100

        glR=0.d0
        suR=0.d0
!         n=n+1 ! compteur sur nième temps
        if (t.ne.1) then
            x=x+cens/100
        end if

! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR
!     fct de survie           : su  = dexp(-gl)
            do i=1,ns
                som=0.d0
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do

                bgen(i)=(b(i)+som)**2
            end do

            lamR(k)=(bgen(1)*(x**(bgen(1)-1.d0)))/(bgen(2)**bgen(1))
            lamRW=((b(1)**2)*(x**((b(1)**2)-1.d0)))/((b(2)**2)**(b(1)**2))

            if(nst == 2) then
                lamDC(k)=(bgen(3)*(x**(bgen(3)-1.d0)))/(bgen(4)**bgen(3))
                lamDCW=((b(3)**2)*(x**((b(3)**2)-1.d0)))/((b(4)**2)**(b(3)**2))
            end if

            glR(k)  =  (x/bgen(2))**bgen(1)
            suR(k)  = dexp(-glR(k))
            glRW  =  (x/(b(2)**2))**(b(1)**2)
            suRW  = dexp(-glRW)

            if(nst == 2) then
                glDC(k)  =  (x/bgen(4))**bgen(3)
                suDC(k)  = dexp(-glDC(k))
                glDCW  =  (x/(b(4)**2))**(b(3)**2)
                suDCW  = dexp(-glDCW)
            end if
        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles
        call percentile(lamR,lamR25,lamR975)
        if(nst == 2) then
            call percentile(lamDC,lamDC25,lamDC975)
        end if

        call percentile(suR,suR25,suR975)
        if(nst == 2) then
            call percentile(suDC,suDC25,suDC975)
        end if

        if(t == 1) then
            x1R(t) = date(1)
        else
            x1R(t) = real(x)
        end if
        moyLamR(t,1) = real(LamRW)
        moyLamR(t,2) = real(lamR25)
        moyLamR(t,3) = real(lamR975)

        if(nst == 2) then
            x2DC(t) = x1R(t)
            moyLamDC(t,1) = real(lamDCW)
            moyLamDC(t,2) = real(lamDC25)
            moyLamDC(t,3) = real(lamDC975)
        end if

        if(t == 1) then
            xSu1(t) = date(1)
        else
            xSu1(t) = real(x)
        end if
        
        moysuR(t,1) = real(suRW)
        moysuR(t,2) = real(suR25)
        moysuR(t,3) = real(suR975)
        if(nst == 2) then
            xSu2(t) = xSu1(t)
            moysuDC(t,1) = real(suDCW)
            moysuDC(t,2) = real(suDC25)
            moysuDC(t,3) = real(suDC975)
        end if

    end do
    if (nst == 1) then
        moyLamDC = 0.d0
        moysuDC = 0.d0
    end if
    end subroutine distanceJweib



!======================================== Joint mulive   ============================
    subroutine distancej_cpm(b,m,mt1,mt2,mt3,x1R,moyLamR,xSu1,moysuR,x2DC,moyLamDC,xSu2,moysuDC,x3M,moyLamM,xSu3,moysuM)
    
    use taillesmultiv
    use comonmultiv,only:cens,vvv,t0,t1,t0dc,t1dc,t0meta,t1meta,cmeta,c,cdc,nsujet,nsujetmeta, &
    nva,nva1,nva2,nva3,nst,nbintervR,nbintervDC,nbintervM,ttt,tttDC,tttmeta,date
    use optim
    
    implicit none
    
    integer::m,i,j,ier,t,k,gg,jj,uu,ii,ncur,mt1,mt2,mt3,l
    double precision::sx,som,som11R,som21R,som11RW,som21RW, &
    glRW,x
    double precision,dimension(m)::b
    double precision,dimension(m)::bgen
    double precision,dimension(1000,m)::u,v
    double precision,dimension(1000)::lamR,lamDC,suR,glR    
    double precision,dimension((m*(m+1)/2))::vv
    double precision,dimension(m)::tempsR,tempsDC,tempsM
    
    double precision::lamR25,lamDC25,suR25,lamR975,lamDC975,suR975, &
    LamRW,LamDCW,suRW
!AD: sorties
    double precision::ep
    double precision,dimension(mt1)::x1R    !0:nbintervR*3
    double precision,dimension(mt1,3)::moyLamR
    double precision,dimension(100,3)::moysuR
    double precision,dimension(mt2)::x2DC !0:nbintervDC*3
    double precision,dimension(mt2,3)::moyLamDC
    double precision,dimension(100,3)::moysuDC
    double precision,dimension(mt3)::x3M    !0:nbintervR*3
    double precision,dimension(mt3,3)::moyLamM
    double precision,dimension(100,3)::moysuM
    double precision,dimension(100)::xSu1,xSu2,xSu3
    
    
    uu=0
    bgen=0.d0
    sx=1.d0 ! ecart-type ou variance des réalisations gaussiennes générées
    
    do i=1,nbintervR
        tempsR(i)=(ttt(i-1)+ttt(i))/(2.d0)
    end do
    
    do i=1,nbintervDC
        tempsDC(i)=(tttdc(i-1)+tttdc(i))/(2.d0)
    end do
    
    
    do i=1,nbintervM
        tempsM(i)=(tttmeta(i-1)+tttmeta(i))/(2.d0)
    end do
        
    do i=1,m*(m+1)/2
        vv(i)=vvv(i)
    end do
!AD: ep=10.d-10 
       ep=10.d-10
    call dmfsdj(vv,m,ep,ier)

!      ns=m-nva-3! ns est le nombre de paramètres intervenat dans la fonctionde
! risque ici 3 car on estime p param : alpha,eta,theta
!     Pour chaque temps de 0 à la censure du décès par pas cens/100
!      x=-cens/100 !permet de commencer à 0 dans la boucle

!      n=0 ! permet de commencer à 0 dans la boucle
!==========================================================================================
!=============== hazard loco
!==========================================================================================
    ncur = 1
    do t=1,nbintervR
        lamR=0.d0
!         n=n+1 ! compteur sur nième temps
        x=tempsR(t)
! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Pour chaque paramètre estimé                         
            do i=1,nbintervR!m
                som=0.d0
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    call bgos(SX,0,u(k,j),v(k,j),0.d0)  
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do
                bgen(i)=(b(i)+som)**2
            end do              ! en sortie on récupère le nouveau vecteur b
                      
!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR        
!     fct de survie           : su  = dexp(-gl)    

            lamR(k)=0.d0
            lamRW=0.d0

            do gg=1,nbintervR
                if ((x.ge.(ttt(gg-1))).and.(x.lt.(ttt(gg)))) then
                    lamR(k)=bgen(gg)
                endif
            end do

            do gg=1,nbintervR
                if ((x.ge.(ttt(gg-1))).and.(x.lt.(ttt(gg)))) then    
                    lamRW=b(gg)**2
                endif
            end do
                    
        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles

        call percentile(lamR,lamR25,lamR975)

        do ii=1,nbintervR
            if (x.gt.ttt(ii-1).and.x.lt.ttt(ii))then
                uu=ii
            endif
        end do
        if(t.ne.1) then 
            x1R(ncur) = real(ttt(uu-1))    
        else
            x1R=date(1)
        end if        
        x1R(ncur+1) = real(x)
        x1R(ncur+2) = real(ttt(uu))
        
        moyLamR(ncur,1) = real(LamRW)
        moyLamR(ncur+1,1) = moyLamR(ncur,1)
        moyLamR(ncur+2,1) = moyLamR(ncur,1)
        
        moyLamR(ncur,2) = real(lamR25)
        moyLamR(ncur+1,2) = moyLamR(ncur,2)
        moyLamR(ncur+2,2) = moyLamR(ncur,2)
        
        
        moyLamR(ncur,3) = real(lamR975)
        moyLamR(ncur+1,3) = moyLamR(ncur,3)
        moyLamR(ncur+2,3) = moyLamR(ncur,3)
        
        ncur=ncur+3
    end do

!==========================================================================================
!=============== hazard dc
!==========================================================================================
    ncur = 1
    do t=1,nbintervDC    
        x=tempsDC(t)
        lamDC=0.d0

!         n=n+1 ! compteur sur nième temps
! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Pour chaque paramètre estimé                         
            do i=1,nbintervDC!m
            
                som=0.d0
                
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    call bgos(SX,0,u(k,j),v(k,j),0.d0)  
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do
                
                bgen(i)=(b(nbintervR + i)+som)**2
            end do              ! en sortie on récupère le nouveau vecteur b
    
    
!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR        
!     fct de survie           : su  = dexp(-gl)    
    
            lamDC(k)=0.d0
            lamDCW=0.d0
    
            do gg=1,nbintervDC
                if((x.ge.(tttdc(gg-1))).and.(x.lt.(tttdc(gg))))then
                    lamDC(k)=bgen(gg)
                endif
            end do            

            do gg=1,nbintervDC
                if((x.ge.(tttdc(gg-1))).and.(x.lt.(tttdc(gg))))then
                    lamDCW=b(gg+nbintervR)**2
                endif
            end do                 
        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles

        call percentile(lamDC,lamDC25,lamDC975)
        

        
        do ii=1,nbintervDC
            if (x.gt.tttdc(ii-1).and.x.lt.tttdc(ii))then
                uu=ii
            endif
        end do
        if(t.ne.1) then 
            x2DC(ncur) = real(tttdc(uu-1))
        else
            x2DC(ncur) = date(1)
        end if
        x2DC(ncur+1) = real(x)
        x2DC(ncur+2) = real(tttdc(uu))

        moyLamDC(ncur,1) = real(lamDCW)
        moyLamDC(ncur+1,1) = moyLamDC(ncur,1)        
        moyLamDC(ncur+2,1) = moyLamDC(ncur,1)        
        
        moyLamDC(ncur,2) = real(lamDC25)
        moyLamDC(ncur+1,2) = moyLamDC(ncur,2)
        moyLamDC(ncur+2,2) = moyLamDC(ncur,2)
                        
        moyLamDC(ncur,3) = real(lamDC975)
        moyLamDC(ncur+1,3) = moyLamDC(ncur,3)
        moyLamDC(ncur+2,3) = moyLamDC(ncur,3)
    
        ncur=ncur+3
    end do

!==========================================================================================
!=============== hazard Meta
!==========================================================================================
    ncur = 1
    do t=1,nbintervM    
        x=tempsM(t)
        lamDC=0.d0

!         n=n+1 ! compteur sur nième temps
! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Pour chaque paramètre estimé                         
            do i=1,nbintervM!m
            
                som=0.d0
                
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    call bgos(SX,0,u(k,j),v(k,j),0.d0)  
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do
                
                bgen(i)=(b(nbintervR + nbintervDC + i)+som)**2
            end do              ! en sortie on récupère le nouveau vecteur b
    
    
!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR        
!     fct de survie           : su  = dexp(-gl)    
    
            lamDC(k)=0.d0
            lamDCW=0.d0
    
            do gg=1,nbintervM
                if((x.ge.(tttmeta(gg-1))).and.(x.lt.(tttmeta(gg))))then
                    lamDC(k)=bgen(gg)
                endif
            end do            

            do gg=1,nbintervM
                if((x.ge.(tttmeta(gg-1))).and.(x.lt.(tttmeta(gg))))then
                    lamDCW=b(gg+nbintervR+nbintervDC)**2
                endif
            end do                 
        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles

        call percentile(lamDC,lamDC25,lamDC975)
        

        
        do ii=1,nbintervM
            if (x.gt.tttmeta(ii-1).and.x.lt.tttmeta(ii))then
                uu=ii
            endif
        end do
        if(t.ne.1) then 
            x3M(ncur) = real(tttmeta(uu-1))
        else
            x3M(ncur) = date(1)
        end if
        x3M(ncur+1) = real(x)
        x3M(ncur+2) = real(tttmeta(uu))

        moyLamM(ncur,1) = real(lamDCW)
        moyLamM(ncur+1,1) = moyLamM(ncur,1)        
        moyLamM(ncur+2,1) = moyLamM(ncur,1)        
        
        moyLamM(ncur,2) = real(lamDC25)
        moyLamM(ncur+1,2) = moyLamM(ncur,2)
        moyLamM(ncur+2,2) = moyLamM(ncur,2)
                        
        moyLamM(ncur,3) = real(lamDC975)
        moyLamM(ncur+1,3) = moyLamM(ncur,3)
        moyLamM(ncur+2,3) = moyLamM(ncur,3)
    
        ncur=ncur+3
    end do
    
    
!==========================================================================================
!=============== Survie loco
!==========================================================================================    
    
!--------------- Fontion de survie recurrent ----------------
    bgen=0.d0
    
    x=date(1)

    do t=1,100

!         n=n+1 ! compteur sur nième temps
        if(t.ne.1) then
            x=x+ttt(nbintervR)/100
        end if
! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Pour chaque paramètre estimé                         
            do i=1,nbintervR!m
                som=0.d0
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    call bgos(SX,0,u(k,j),v(k,j),0.d0)  
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do
                bgen(i)=(b(i)+som)**2
            end do              ! en sortie on récupère le nouveau vecteur b
                      
!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR        
!     fct de survie           : su  = dexp(-gl)    

            som11R=0.d0
            som21R=0.d0    
            glR(k)=0.d0
            suR(k)=0.d0

            do gg=1,nbintervR
                if ((x.ge.(ttt(gg-1))).and.(x.lt.(ttt(gg)))) then
                
                    som11R=bgen(gg)*(x-ttt(gg-1))
                    
                    if (gg.ge.2)then
                        do jj=1,gg-1
                            som21R=som21R+bgen(jj)*(ttt(jj)-ttt(jj-1))
                        end do
                    endif
                    
                    glR(k)=(som11R+som21R)
                    suR(k)=dexp(-glR(k))
                endif
            end do
       
        end do

        som11RW=0.d0
        som21RW=0.d0
        glRW=0.d0
        suRW=0.d0
        do l=1,nbintervR
            if ((x.ge.(ttt(l-1))).and.(x.lt.(ttt(l)))) then
            
                som11RW=(b(l)**2)*(x-ttt(l-1))
                
                if (l.ge.2)then
                    do jj=1,l-1
                        som21RW=som21RW+(b(jj)**2)*(ttt(jj)-ttt(jj-1))
                    end do
                endif
                
                glRW=(som11RW+som21RW)
                suRW=dexp(-glRW)
            endif
        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles
        suR25=0.d0
        suR975=0.d0
        call percentile(suR,suR25,suR975)
        if(t.ne.1) then
            xSu1(t) = real(x)
        else
            xSu1(t) = date(1)
        end if
        moysuR(t,1) = suRW

        moysuR(t,2) = suR25

        moysuR(t,3) = suR975
        
        if(moysuR(t,1).lt.0.d0)then
            moysuR(t,1) = 0.d0
        endif
        if(moysuR(t,1).gt.1.d0)then
            moysuR(t,1) = 1.d0 
        endif
        
        if(moysuR(t,2).lt.0.d0)then
            moysuR(t,2) = 0.d0
        endif
        if(moysuR(t,2).gt.1.d0)then
            moysuR(t,2) = 1.d0 
        endif    
        if(moysuR(t,3).lt.0.d0)then
            moysuR(t,3) = 0.d0
        endif
        if(moysuR(t,3).gt.1.d0)then
            moysuR(t,3) = 1.d0 
        endif        

    end do
!==========================================================================================
!=============== Survie dc
!==========================================================================================    
!--------------- Fontion de survie ----------------
    bgen=0.d0

    x=date(1)
    do t=1,100

!         n=n+1 ! compteur sur nième temps
        if(t.ne.1) then
            x=x+tttdc(nbintervDC)/100
        end if
        
! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Pour chaque paramètre estimé                         
            do i=1,nbintervDC!m
                som=0.d0
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    call bgos(SX,0,u(k,j),v(k,j),0.d0)  
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do
                bgen(i)=(b(i+nbintervR)+som)**2
            end do              ! en sortie on récupère le nouveau vecteur b
                      
!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR        
!     fct de survie           : su  = dexp(-gl)    

            som11R=0.d0
            som21R=0.d0    
            glR(k)=0.d0
            suR(k)=0.d0



            do gg=1,nbintervDC
                if ((x.ge.(tttdc(gg-1))).and.(x.lt.(tttdc(gg)))) then
                
                    som11R=bgen(gg)*(x-tttdc(gg-1))
                    
                    if (gg.ge.2)then
                        do jj=1,gg-1
                            som21R=som21R+bgen(jj)*(tttdc(jj)-tttdc(jj-1))
                        end do
                    endif
                    
                    glR(k)=(som11R+som21R)
                    suR(k)=dexp(-glR(k))
                endif
            end do


                       
        
!            moysuR0=moysuR0+suR(k)/1000
            
        end do
    
        som11RW=0.d0
        som21RW=0.d0
        glRW=0.d0
        suRW=0.d0
        do l=1,nbintervDC
            if ((x.ge.(tttdc(l-1))).and.(x.lt.(tttdc(l)))) then
            
                som11RW=(b(l+nbintervR)**2)*(x-tttdc(l-1))
                
                if (l.ge.2)then
                    do jj=1,l-1
                        som21RW=som21RW+(b(jj+nbintervR)**2)*(tttdc(jj)-tttdc(jj-1))
                    end do
                endif
                
                glRW=(som11RW+som21RW)
                suRW=dexp(-glRW)
            endif
        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles
        suR25=0.d0
        suR975=0.d0
        call percentile(suR,suR25,suR975)
        if(t.ne.1) then
            xSu2(t) = real(x)
        else
            xSu2(t) = date(1)
        end if

        moysuDC(t,1) = suRW

        moysuDC(t,2) = suR25

        moysuDC(t,3) = suR975
        
        if(moysuDC(t,1).lt.0.d0)then
            moysuDC(t,1) = 0.d0
        endif
        if(moysuDC(t,1).gt.1.d0)then
            moysuDC(t,1) = 1.d0 
        endif
        
        if(moysuDC(t,2).lt.0.d0)then
            moysuDC(t,2) = 0.d0
        endif
        if(moysuDC(t,2).gt.1.d0)then
            moysuDC(t,2) = 1.d0 
        endif    
        if(moysuDC(t,3).lt.0.d0)then
            moysuDC(t,3) = 0.d0
        endif
        if(moysuDC(t,3).gt.1.d0)then
            moysuDC(t,3) = 1.d0 
        endif    

    end do
!==========================================================================================
!=============== Survie loco
!==========================================================================================    
!--------------- Fontion de survie ----------------
    bgen=0.d0

    x=date(1)
    do t=1,100

!         n=n+1 ! compteur sur nième temps
        if(t.ne.1) then
            x=x+tttmeta(nbintervM)/100
        end if
        
! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Pour chaque paramètre estimé                         
            do i=1,nbintervM!m
                som=0.d0
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                    call bgos(SX,0,u(k,j),v(k,j),0.d0)  
                    som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do
                bgen(i)=(b(i+nbintervR+nbintervDC)+som)**2
            end do              ! en sortie on récupère le nouveau vecteur b
                      
!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR        
!     fct de survie           : su  = dexp(-gl)    

            som11R=0.d0
            som21R=0.d0    
            glR(k)=0.d0
            suR(k)=0.d0



            do gg=1,nbintervM
                if ((x.ge.(tttmeta(gg-1))).and.(x.lt.(tttmeta(gg)))) then
                
                    som11R=bgen(gg)*(x-tttmeta(gg-1))
                    
                    if (gg.ge.2)then
                        do jj=1,gg-1
                            som21R=som21R+bgen(jj)*(tttmeta(jj)-tttmeta(jj-1))
                        end do
                    endif
                    
                    glR(k)=(som11R+som21R)
                    suR(k)=dexp(-glR(k))
                endif
            end do


                       
        
!            moysuR0=moysuR0+suR(k)/1000
            
        end do
    
        som11RW=0.d0
        som21RW=0.d0
        glRW=0.d0
        suRW=0.d0
        do l=1,nbintervM
            if ((x.ge.(tttmeta(l-1))).and.(x.lt.(tttmeta(l)))) then
            
                som11RW=(b(l+nbintervR+nbintervDC)**2)*(x-tttmeta(l-1))
                
                if (l.ge.2)then
                    do jj=1,l-1
                        som21RW=som21RW+(b(jj+nbintervR+nbintervDC)**2)*(tttmeta(jj)-tttmeta(jj-1))
                    end do
                endif
                
                glRW=(som11RW+som21RW)
                suRW=dexp(-glRW)
            endif
        end do

! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles
        suR25=0.d0
        suR975=0.d0
        call percentile(suR,suR25,suR975)
        if(t.ne.1) then
            xSu3(t) = real(x)
        else
            xSu3(t) = date(1)
        end if

        moysuM(t,1) = suRW

        moysuM(t,2) = suR25

        moysuM(t,3) = suR975
        
        if(moysuM(t,1).lt.0.d0)then
            moysuM(t,1) = 0.d0
        endif
        if(moysuM(t,1).gt.1.d0)then
            moysuM(t,1) = 1.d0 
        endif
        
        if(moysuM(t,2).lt.0.d0)then
            moysuM(t,2) = 0.d0
        endif
        if(moysuM(t,2).gt.1.d0)then
            moysuM(t,2) = 1.d0 
        endif    
        if(moysuM(t,3).lt.0.d0)then
            moysuM(t,3) = 0.d0
        endif
        if(moysuM(t,3).gt.1.d0)then
            moysuM(t,3) = 1.d0 
        endif    

    end do

    end subroutine distancej_cpm


    subroutine distanceJ_weib(b,m,mt1,x1R,moyLamR,xSu1,moysuR,x2DC,moyLamDC,xSu2,moysuDC,x3M,moyLamM,xSu3,moysuM)
    
    use taillesmultiv
    use comonmultiv,only:cens,vvv,t0,t1,t0dc,t1dc,t0meta,t1meta,c,cmeta,cdc,nsujet,nsujetmeta &
    ,nva,nva1,nva2,nva3,nst,etaR,etaD,etaM,betaR,betaD,betaM,date,datedc,datemeta
    use optim
    
    implicit none
    
    integer::m,i,j,ier,t,k,ns,mt1
    double precision::lamR25,lamM25,lamDC25,suR25,suM25,suDC25,lamR975,lamM975,lamDC975,suR975,suM975,suDC975, &
    LamRW,LamMW,LamDCW,suRW,suMW,suDCW 
        double precision,dimension(m*(m+1)/2)::vv
    double precision,dimension(m)::b,bgen
    double precision,dimension(1000,m)::u
    double precision::sx,som,x,zz,zy 
    double precision,dimension(1000):: lamR,lamDC,lamM,suR,suDC,suM,glDC,glR,glM
    
! theorique - estimés après Max de Vrais
      double precision::glDCW,glRW,glMW

    double precision,dimension(mt1)::x1R,x2DC,x3M
    double precision,dimension(mt1,3)::moyLamR,moyLamDC,moyLamM
    double precision,dimension(100,3)::moysuR,moysuDC,moysuM    
    double precision,dimension(100)::xSu1,xSu2,xSu3
    double precision::ep
    
    
    lamR25=0.d0
    lamM25=0.d0
    lamDC25=0.d0
    suR25=0.d0
    suM25=0.d0    
    suDC25=0.d0
    lamR975=0.d0
    lamM975=0.d0    
    lamDC975=0.d0
    suR975=0.d0
    suM975=0.d0
    suDC975=0.d0
    LamRW=0.d0
    LamMW=0.d0    
    LamDCW=0.d0
    suRW=0.d0
    suMW=0.d0
    suDCW=0.d0
    
    sx=1.d0 ! ecart-type ou variance des réalisations gaussiennes générées
    ns=m-nva-2

    do i=1,ns*(ns+1)/2
        vv(i)=vvv(i)
    end do

    ep=10.d-10
    call dmfsdj(vv,m,ep,ier)

    do k=1,1000
        do j=1,ns   
            call bgos(SX,0,zz,zy,0.d0)
            u(k,j)=zz
        end do
    end do

    betaR = b(1)**2
    etaR =  b(2)**2
    

    if (nst == 3) then
        betaD = b(3)**2
        etaD = b(4)**2
        betaM = b(5)**2
        etaM = b(6)**2        
    else
        etaD = 0.d0
        betaD = 0.d0
    end if


!============================ Surv strat 1    
!------------------- Pour les recurents
!     Pour chaque temps de 0 à la censure du décès par pas cens/100
!    x=-cens/100 !permet de commencer à 0 dans la boucle

    x=date(1)
    do t=1,100

        glR=0.d0
        suR=0.d0  
!         n=n+1 ! compteur sur nième temps
        if (t.ne.1) then
            x=x+cens/100
        end if

! On simule 1000 réalisations gaussienne par paramètres
        do k=1,1000
!     Rappel des notations:
!     Ce sont les paramètres des fcts de type weibull
!     betaR= b(1)**2
!     etaR= b(2)**2
!     betaD= b(3)**2
!     etaD= b(4)**2
!     fct de risque ou hazard : lam=(betaR*(x**(betaR-1.d0)))/(etaR**betaR)
!     fct de risque cumulée   : gl  =  (x/etaR)**betaR
!     fct de survie           : su  = dexp(-gl)
            do i=1,ns
                som=0.d0
                do j=1,i         ! cela correspond au produit trp(L)%*%U
                som=som+vv(i*(i-1)/2+j)*u(k,j)
                end do

                bgen(i)=(b(i)+som)**2
            end do

            lamR(k)=(bgen(1)*(x**(bgen(1)-1.d0)))/(bgen(2)**bgen(1))
            lamRW=((b(1)**2)*(x**((b(1)**2)-1.d0)))/((b(2)**2)**(b(1)**2))

            if(nst == 3) then
                lamDC(k)=(bgen(3)*(x**(bgen(3)-1.d0)))/(bgen(4)**bgen(3))
                lamDCW=((b(3)**2)*(x**((b(3)**2)-1.d0)))/((b(4)**2)**(b(3)**2))
                
                lamM(k)=(bgen(5)*(x**(bgen(5)-1.d0)))/(bgen(6)**bgen(5))
                lamMW=((b(5)**2)*(x**((b(5)**2)-1.d0)))/((b(6)**2)**(b(5)**2))                
            end if

            glR(k)  =  (x/bgen(2))**bgen(1)
            suR(k)  = dexp(-glR(k))
            glRW  =  (x/(b(2)**2))**(b(1)**2)
            suRW  = dexp(-glRW)
            
            if(nst == 3) then
                glDC(k)  =  (x/bgen(4))**bgen(3)
                suDC(k)  = dexp(-glDC(k))
                glDCW  =  (x/(b(4)**2))**(b(3)**2)
                suDCW  = dexp(-glDCW)    

                glM(k)  =  (x/bgen(6))**bgen(5)
                suM(k)  = dexp(-glM(k))
                glMW  =  (x/(b(6)**2))**(b(5)**2)
                suMW  = dexp(-glMW)                    
            end if
        end do
    
! Classer les différent vecteur et en sortir les 2.5 et 97.5 percentiles
        call percentile(lamR,lamR25,lamR975)
        if(nst == 3) then
            call percentile(lamDC,lamDC25,lamDC975)
            call percentile(lamM,lamM25,lamM975)            
        end if

        call percentile(suR,suR25,suR975)
        if(nst == 3) then
            call percentile(suDC,suDC25,suDC975)
            call percentile(suM,suM25,suM975)            
        end if
        
        if(t == 1) then
            x1R(t) = date(1)
        else
            x1R(t) = real(x)
        end if
        moyLamR(t,1) = real(LamRW)
        moyLamR(t,2) = real(lamR25)
        moyLamR(t,3) = real(lamR975)
        
        if(nst == 3) then
            x2DC(t) = x1R(t)
            moyLamDC(t,1) = real(lamDCW)
            moyLamDC(t,2) = real(lamDC25)
            moyLamDC(t,3) = real(lamDC975)
            
            x3M(t) = x1R(t)
            moyLamM(t,1) = real(lamMW)
            moyLamM(t,2) = real(lamM25)
            moyLamM(t,3) = real(lamM975)            
            
        end if

        if(t == 1) then
            xSu1(t) = date(1)
        else
            xSu1(t) = real(x)
        end if
        
        moysuR(t,1) = real(suRW)
        moysuR(t,2) = real(suR25)
        moysuR(t,3) = real(suR975)
            
        if(nst == 3) then
            xSu2(t) = xSu1(t) 
            moysuDC(t,1) = real(suDCW)
            moysuDC(t,2) = real(suDC25)
            moysuDC(t,3) = real(suDC975)
            
            xSu3(t) = xSu1(t) 
            moysuM(t,1) = real(suMW)
            moysuM(t,2) = real(suM25)
            moysuM(t,3) = real(suM975)            
            
        end if
            
    end do
    if (nst == 1) then
        moyLamDC = 0.d0
        moysuDC = 0.d0
        moysuM = 0.d0
    end if
    end subroutine distanceJ_weib
    


