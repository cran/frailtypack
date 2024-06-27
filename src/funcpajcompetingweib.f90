! This file contains 

! funcpaMultivSplines (function): computes the log likelihood for a vector of
! parameters b given the data stored in the module "comonmultiv"


!========================  funcpajcompetingWeib  ====================
! Log likelihood for one recurrent event and two terminal events.
    double precision function funcpajcompetingWeib(b,np,id,thi,jd,thj,k0)

    !use taillesjcompeting
    !use comonjcompeting
    !use residusMjcompeting

    use taillesmultiv
    use comonmultiv
    use residusMmultiv
    use parametersmultiv,only:constraintfrailty
    implicit none

! Variable Definitions
    integer,intent(in)::id,jd,np 
      ! id, jd = indexes at which to approximate first, second derivatives
      ! np = total number parameters, items in vectors b, bh
    double precision,dimension(np),intent(in)::b
      ! b = parameter vector
    double precision,dimension(4)::k0
        ! k0 = smoothing parameters
    double precision,intent(in)::thi,thj
        ! thi, thj = small differences for approximation of first, second derivatives
    integer::n,i,j,k,l,vj,ig,choix
    integer,dimension(ngmax)::cpt!, cptmeta
      ! cpt, cptmax = number of recurrent events for each person
    double precision::res,vet,vet2,vet4,frail,frail2,weight,weight2
        ! vet = temporary value for relative hazard for an individual
        ! frail = used in GH quadrature for integration
        ! weight = used in GH quadrature for integration
    double precision,dimension(np)::bh
        ! bh = local copy of parameter vector
    double precision,dimension(ng)::res2,res1dc,res2dc &
    ,res3dc,res2dc2,integrale3!,res2meta
    double precision,parameter:: pi=3.141592653589793d0

    double precision::ss

    double precision,dimension(6)::restest
    !INTEGER::ii,jj,npg ! delete npg?
    !double precision,dimension(30)::x,w

    !open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
    !write(1,*) "funcpajcompetingWeib.f90:: Setup."
    !close(1)
    funcpajcompetingWeib = 0.0d0 
    !write(*,*) "Computing log likelihood..."
    kkapa=k0
    choix=0
    ig=0
    k=0
    vj=0
    n=0
    j=0
    l = 0

    do i=1,np
        bh(i)=b(i)
    end do

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj


    n = (np-nva-effet-indic_ALPHA)/nst ! delete this?

    ! One Random Effect Model
    if(typeJoint.eq.0)then
        ! b = c(betaR, etaR, betaD, etaD, betaD2, etaD2,theta,alpha1, alpha2, coef...)
        betaR= bh(1)**2.0d0 !shape
        etaR= bh(2)**2.0d0 !scale
        betaD= bh(3)**2.0d0
        etaD= bh(4)**2.0d0
        betaD2= bh(5)**2.0d0
        etaD2= bh(6)**2.0d0 
        if(constraintfrailty.eq.0) then 
            theta = bh(7)**2.0d0
        else 
            theta = dexp(bh(7))**2.0d0 ! log_e random effect sd ---> varaince
        end if 
        alpha1=bh(8) ! random effect link to terminal event process 1
        alpha2=bh(9) ! random effect link to terminal event process 2
        rho = 0.d0
        theta2 = 1.d0
    endif
    if(typeJoint.eq.1)then
        ! Two Random Effect Model
        ! b = c(betaR, etaR, betaD, etaD, betaD2, etaD2, theta, theta2, rho, alpha1, alpha2, coef...)
        betaR= bh(1)**2.0d0 !shape
        etaR= bh(2)**2.0d0 !scale
        betaD= bh(3)**2.0d0
        etaD= bh(4)**2.0d0
        betaD2= bh(5)**2.0d0
        etaD2= bh(6)**2.0d0
        theta = dexp(bh(7))**2.0d0! log_e random effect sd ---> varaince
        theta2 = dexp(bh(8))**2.0d0! log_e random effect sd ---> varaince
        rho = (2.d0*dexp(bh(9))/(1.d0+dexp(bh(9)))) - 1.d0 ! transform correlation from unbounded to (-1,1) using inverse logit
        alpha1=bh(10) ! random effect link to terminal event process 1
        alpha2=bh(11) ! random effect link to terminal event process 2
    endif

    if (abs(rho).ge.1.d0) then
       !debug: open(1, file = '../package_tests/multiv_model_progress.dat',position="append")  
       !debug: write(1,*) "funcpajcompetingWeib.f90:: Correlation out of bounds."
       !debug: close(1)
        funcpajcompetingWeib =-1.d9
        do k=1,ng
            Rrec(k)=0.d0
            Nrec(k)=0
            Rdc(k)=0.d0
            Ndc(k)=0
            !Rrec2(k)=0.d0
            !Nrec2(k)=0
        end do
        goto 123
    endif

    ! Create empty vectors for likelihood components

    do k=1,ng
        res1(k) = 0.d0
        res2(k) = 0.d0
        res3(k) = 0.d0

        res1dc(k) = 0.d0
        res2dc(k) = 0.d0
        res3dc(k) = 0.d0

        res2dc2(k) = 0.d0

        integrale3(k) = 1.d0
        aux1(k)=0.d0
        aux2(k)=0.d0
    end do


!==========================================================================
! Likelihood Contribution for Recurrent Event 1
!==========================================================================
    cpt = 0
    restest = 0.0d0
    do i=1,nsujet
        ! (1) count the number of events per individual
        cpt(g(i))=cpt(g(i))+c(i)

        ! (2) compute the relative hazard for each event
        if(nva1.gt.0)then
            vet = 0.d0
            do j=1,nva1
                vet = vet + bh(np-nva1-nva2-(nva3*event2_ind)-(nva4*terminal2_ind)+j)*dble(ve(i,j))
            end do
            vet = dexp(vet)
        else
            vet=1.d0
        endif

        ! (3a) Compute hazard for non-censored events
        if((c(i).eq.1))then

            
            res2(g(i)) = res2(g(i))+(betaR-1.d0)*dlog(t1(i))+&
            dlog(betaR)-betaR*dlog(etaR)+dlog(vet)

            restest(1) = restest(1) + res2(g(i))


        endif

        if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge.1.d30)) then
            funcpajcompetingWeib=-1.d9
            goto 123
        end if

        ! (3b) Compute cumulative hazard for all events (censored and non-censored)
        ! this part goes into the integral

        res1(g(i)) = res1(g(i))+((t1(i)/etaR)**betaR)*vet 
        !res1(g(i)) = res1(g(i))+vet

         if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
            funcpajcompetingWeib=-1.d9
            goto 123
        end if

        ! (3c) Modification for cumulative hazard with left truncation (for calendar time model 
        ! as opposed to gap time model)
        res3(g(i)) = res3(g(i))+((t0(i)/etaR)**betaR)*vet
        !res3(g(i)) = res3(g(i))+vet
        if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
            funcpajcompetingWeib=-1.d9
            goto 123
        end if
    end do

!==========================================================================
! Likelihood Contribution for Terminal Event 1
!==========================================================================

    do k=1,ng
        if(nva2.gt.0)then
            vet2 = 0.d0
            do j=1,nva2
                vet2 =vet2 + bh(np-nva2-(nva3*event2_ind)-(nva4*terminal2_ind)+j)*dble(vedc(k,j))
            end do
            vet2 = dexp(vet2)
        else
            vet2=1.d0
        endif
        if(cdc(k).eq.1)then

            res2dc(k) = (betaD-1.d0)*dlog(t1dc(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet2)

            if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge.1.d30)) then
                funcpajcompetingWeib=-1.d9
                goto 123
            end if
        endif

        ! cumulative hazard, goes into the integral

        aux1(k)=((t1dc(k)/etaD)**betaD)*vet2
        !aux1(k)=vet2

        if ((aux1(k).ne.aux1(k)).or.(abs(aux1(k)).ge. 1.d30)) then
            funcpajcompetingWeib=-1.d9
            goto 123
        end if
    end do

!==========================================================================
! Likelihood Contribution for Terminal Event 2
!==========================================================================

    do k=1,ng
        if(nva4.gt.0)then
            vet4 = 0.d0
            do j=1,nva4
                vet4 =vet4 + bh(np-(nva4*terminal2_ind)+j)*dble(vedc2(k,j))
            end do
            vet4 = dexp(vet4)
        else
            vet4=1.d0
        endif
        if(cdc2(k).eq.1)then

            res2dc2(k) = (betaD2-1.d0)*dlog(t1dc(k))+dlog(betaD2)-betaD2*dlog(etaD2)+dlog(vet4)
            !res2dc2(k) = dlog(vet4)
            if ((res2dc2(k).ne.res2dc2(k)).or.(abs(res2dc2(k)).ge. 1.d30)) then
                funcpajcompetingWeib=-1.d9
                goto 123
            end if
        endif

        ! cumulative hazard, goes into the integral
        aux2(k)=((t1dc(k)/etaD2)**betaD2)*vet4
        !aux2(k)=vet4
        
        if ((aux2(k).ne.aux2(k)).or.(abs(aux2(k)).ge. 1.d30)) then
            funcpajcompetingWeib=-1.d9
            goto 123
        end if
    end do



    !==========================================================================
    ! Integrals
    !==========================================================================
    ! compute integral for each person
    ! the integral component will depend on the parameterization of the random effects
    ! typeJoint = 0 for 1 random effect
    ! typeJoint = 1 for 2 correlated random effects

    if(typeJoint.eq.0)then
        do k=1, ng
            ss=0.d0
            do j=1, ghPoints
                    weight = ghWeights(j)
                    frail = ghNodes(j)
                    ss = ss + &
                    (weight* dexp( & ! GQ Weights
                    frail * cpt(k) - dexp(frail)*(res1(k)-res3(k)) & ! recurrent event 1
                    + alpha1 * frail * cdc(k) - dexp(frail*alpha1) * aux1(k) & ! terminal event 1
                    + alpha2 * frail * cdc2(k) - dexp(frail*alpha2) * aux2(k) & ! terminal event 2
                    - (frail**2.d0)/(2.d0*theta)))! frailty distribution (normal) 
            end do
            if((ss.eq.0.0d0)) ss = 1.0d-10 ! integral should not be too small when theta is large
            integrale3(k) = ss
            if ((integrale3(k).ne.integrale3(k)).or.(abs(integrale3(k)).ge.1.d30)) then
                funcpajcompetingWeib=-1.d9
                goto 123
            end if
        end do
    else ! typeJoint = 1
        do k=1, ng
            ss=0.d0
            do j=1, ghPoints ! randdom effect for terminal 1
                do l = 1, ghPoints ! random effect for terminal 2
                    weight = ghWeights(j)
                    weight2 = ghWeights(l)
                    frail = ghNodes(j)
                    frail2 = ghNodes(l)
                    ss = ss + &
                    (weight * weight2 * dexp( & ! GQ Weights
                    (frail + frail2) * cpt(k) - dexp(frail + frail2)*(res1(k)-res3(k)) & ! recurrent event 1
                    + alpha1 * frail * cdc(k) - dexp(frail*alpha1) * aux1(k) & ! terminal event 1
                    + alpha2 * frail2 * cdc2(k) - dexp(frail2*alpha2) * aux2(k) & ! terminal event 2
                    - (((frail**2.d0)/theta) + ((frail2**2.d0)/theta2) - ((2.0d0*frail*frail2*rho)/dsqrt(theta)/dsqrt(theta2)))&
                    /(2.d0*(1.d0-(rho**2.d0)))))! frailty distribution (multivariate normal)
                end do
            end do
            if((ss.eq.0.0d0)) ss = 1.0d-10 ! integral should not be too small when theta is large
            integrale3(k) = ss
            if ((integrale3(k).ne.integrale3(k)).or.(abs(integrale3(k)).ge.1.d30)) then
                funcpajcompetingWeib=-1.d9
                goto 123
            end if
        end do
        !open(1, file = '../package_tests/multiv_model_integrals.dat',position="append")  
        !write(1,*) integrale3
        !close(1)
    end if

    !==========================================================================
    ! Combine likelihood
    !==========================================================================


    res = 0.d0
    if(typeJoint.eq.0)then
        do k=1,ng
            res = res + res2(k)+res2dc(k)+res2dc2(k)+dlog(integrale3(k))- &
            dlog(2.d0*pi)/2.d0 - dlog(dsqrt(theta))
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                funcpajcompetingWeib =-1.d9
                Rrec(k)=0.d0
                Nrec(k)=0
                Rdc(k)=0.d0
                Ndc(k)=0
                goto 123
            else
                funcpajcompetingWeib = res
                Rrec(k)=res1(k)
                Nrec(k)=nig(k)
                Rdc(k)=aux1(k)
                Ndc(k)=cdc(k)
            end if
        end do
    end if
    if(typeJoint.eq.1)then
        do k=1,ng
            res = res + res2(k)+res2dc(k)+res2dc2(k)+dlog(integrale3(k))- &
            dlog(dsqrt(theta)) - dlog(dsqrt(theta2)) - dlog(dsqrt(1.0d0-(rho**2.0d0))) - dlog(2.d0*pi)
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                funcpajcompetingWeib =-1.d9
                Rrec(k)=0.d0
                Nrec(k)=0
                Rdc(k)=0.d0
                Ndc(k)=0
                goto 123
            else
                funcpajcompetingWeib = res
                Rrec(k)=res1(k)
                Nrec(k)=nig(k)
                Rdc(k)=aux1(k)
                Ndc(k)=cdc(k)
            end if
        end do
    end if
    123     continue
    return

    

    end function funcpajcompetingWeib

!=================================================================================================
