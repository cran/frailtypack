    
!========================    FUNCPARES RESIDUS MATRINGALE DENSITE A POSTERIORI       ====================

!!!!
!!!! Calcul Residus shared log-normal
!!!!
    double precision function funcpasres(uu,np,id,thi,jd,thj)

    use comon
    use residusM

    implicit none

    integer,intent(in)::id,jd,np
    double precision,dimension(np)::bh
    double precision,dimension(np),intent(in)::uu
    double precision,intent(in)::thi,thj
    double precision::frail

    bh=uu
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    frail = bh(1)*bh(1)

    funcpasres = dexp(nig(indg)*frail - &
    dexp(frail)*cumulhaz(indg) - (frail**2.d0)/(2.d0*sig2))
    
    return

    end function funcpasres


!!!!
!!!! Calcul Residus joint
!!!!
! la differenciation dans le calcul entre joint cluster
! et classique se fait dans le funcpa

    double precision function funcpajres(uu,np,id,thi,jd,thj)

    use comon
        use residusM
    
    implicit none

    integer,intent(in)::id,jd,np
    double precision,intent(in)::thi,thj
    double precision,dimension(np)::uu,bh
    double precision::frail1,res
    double precision,parameter::pi=3.141592653589793d0

    bh=uu
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj
    
    frail1=bh(1)*bh(1)

!    funcpajres = frail1**(Ndc(indg) + Nrec(indg) + 1.d0/theta - 1.d0 + &
!    alpha * (Nrec(indg) + Ndc(indg))) * dexp(-frail1*(1.d0/theta + &
!    Rrec(indg))) * dexp(-(frail1**alpha)*Rdc(indg))

    res = frail1**(Nrec(indg) + 1.d0/theta - 1.d0 + &
    alpha * Ndc(indg)) * dexp(-frail1*(1.d0/theta + &
    Rrec(indg))) * dexp(-(frail1**alpha)*Rdc(indg))

!    res = dexp(-(frail1**alpha)*Rdc(indg)) * frail1**Nrec(indg) &
!    * frail1**(1.d0/theta - 1.d0) * frail1**(alpha * Ndc(indg)) &
!    * dexp(-frail1*(1.d0/theta + Rrec(indg)))

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then ! limite Ã  e+300, ce n'est pas une log-vraisemblance
        funcpajres=-1.d9
        goto 222
    end if

    funcpajres = res

222    continue

    return
    
    end function funcpajres

!!!!
!!!! Calcul Residus joint log-normal
!!!!

    double precision function funcpajres_log(uu,np,id,thi,jd,thj)

    use comon
        use residusM
    
    implicit none

    integer,intent(in)::id,jd,np
    double precision,intent(in)::thi,thj
    double precision,dimension(np)::uu,bh
    double precision::frail,res
    double precision,parameter::pi=3.141592653589793d0

    bh=uu
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj
    
    frail=bh(1)*bh(1)

    res = dexp((Nrec(indg)+alpha*Ndc(indg))*frail- &
    dexp(frail)*Rrec(indg) - dexp(alpha*frail)*Rdc(indg)- &
    (frail**2.d0)/(2.d0*sig2))

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpajres_log=-1.d9
        goto 222
    end if

    funcpajres_log = res

222    continue

    return
    
    end function funcpajres_log

!!!!
!!!! Calcul Residus nested
!!!!

    double precision function funcpanres(uu,np,id,thi,jd,thj)
    
    use comon,only:alpha,eta
        use residusM
    use commun

    
    implicit none

    integer,intent(in)::id,jd,np
    double precision,intent(in)::thi,thj    
    double precision,dimension(np),intent(in)::uu
    integer::j
    double precision,dimension(np)::bh
    double precision::frail1,prod1,prod2,prod3,res
    double precision,dimension(np-1)::frail2

    bh=uu
    
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj    
    
    frail1=bh(1)*bh(1)
    
    do j=1,n_ssgbygrp(indg)
        frail2(j)=bh(j+1)*bh(j+1)
    end do

    prod1 = 1.d0
    prod2 = 1.d0
    prod3 = 1.d0
    
    do j=1,n_ssgbygrp(indg)
        prod1 = prod1 * (frail2(j)**mij(indg,j)) * dexp(-frail1 * frail2(j) * cumulhaz1(indg,j))
        prod2 = prod2 * frail2(j)**((1.d0/eta) - 1.d0) * dexp(-frail2(j)/eta)
        prod3 = prod3 * dexp(-frail1 * frail2(j) * cumulhaz0(indg,j))
    end do

    res = frail1**(mid(indg)+1.d0/alpha - 1.d0) * prod1 * prod3 * dexp(-frail1/alpha) * prod2
    
    if ((res.ne.res).or.(abs(res).ge. 1.d300)) then
        funcpanres=-1.d9
        goto 333
    end if

    funcpanres = res

333    continue

    return
    
    end function funcpanres    
    
!!!!
!!!! Calcul Residus additive
!!!!

    double precision function funcpaares(uu,np,id,thi,jd,thj)
    
    use comon
        use residusM
    use additiv,only:mid,Xbeta,ve2,ut1,ut2    

    implicit none

    integer,intent(in)::id,jd,np
    double precision,intent(in)::thi,thj
    double precision,dimension(np),intent(in)::uu
    integer::k,ip
    double precision,dimension(np)::bh
    double precision::frail1,frail2,som1,som2
    double precision,parameter::pi=3.141592653589793d0
    double precision::result
    double precision,dimension(1,2)::apres
    double precision,dimension(2,1)::avant    
    double precision,dimension(1,1)::res    
    double precision::vet


    bh=uu
!    kapa=k0(1)
        
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj    
    
    frail1=bh(1)
    frail2=bh(2)
    
    apres(1,1)=frail1
    apres(1,2)=frail2    
    
    avant(1,1)=frail1
    avant(2,1)=frail2
    
    res = (-1.d0/2) * matmul(apres,matmul(invsigma,avant))    
    result = res(1,1)

    som1 = 0.d0
    do k=1,nsujet
        if(g(k) == indg)then
            if(c(k).eq.1)then
                som1 = som1 + frail1 + frail2 * ve2(k,1) + dlog(som_Xbeta(indg))
            end if
        end if
    end do
    
    som2 = 0.d0

    do k=1,nsujet
        if(nva.gt.0 .and.g(k).eq.indg)then
            vet = 0.d0 
            do ip=1,nva
                vet = vet + b_temp(np-nva+ip)*ve(k,ip)
            end do
            vet = dexp(vet)
        else
            vet=1.d0
        endif
        if(typeof==0) then
            if(g(k) == indg)then    
                if(stra(k).eq.1)then
                    som2 = som2  - 1.d0 * ut1(nt1(k)) * dexp(frail1 + frail2 * ve2(k,1) + dlog(vet))
                end if
                if(stra(k).eq.2)then
                    som2 = som2 - 1.d0 * ut2(nt1(k)) * dexp(frail1 + frail2 * ve2(k,1) + dlog(vet))
                end if
            end if
        else
            if(g(k) == indg)then    
                som2 = som2  - 1.d0 * cumulhaz(g(k)) * dexp(frail1 + frail2 * ve2(k,1) + dlog(vet))
            end if

        end if
    end do

    funcpaares = 1.d0/(2.d0 * pi * dsqrt(detSigma))* dexp(som1 + som2 + result)
    
    return
    
    end function funcpaares    




!========================    FUNCPARES MULTIVE RESIDUS MATRINGALE DENSITE A POSTERIORI       ====================


    double precision function funcpamultires(uu,np,id,thi,jd,thj)
    
    use comonmultiv,only:eta,theta,alpha,alpha1,alpha2
    use residusMmultiv
    
    IMPLICIT NONE
    
    integer,intent(in)::id,jd,np
    double precision,dimension(np)::bh
    double precision,dimension(np),intent(in)::uu
!    double precision,dimension(3),intent(in)::k0
    double precision,intent(in)::thi,thj
    double precision::frail1,frail2
    double precision,parameter::pi=3.141592653589793d0
!    double precision,dimension(3)::kappa_tmp

!    kappa_tmp = k0
    bh=uu

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj    

    
    frail1=bh(1)
    frail2=bh(2) 


!---------- calcul de la penalisation -------------------

    funcpamultires= frail1*(Ndc(indg)*alpha1+Nrec(indg)) &
    +frail2*(Ndc(indg)*alpha2+Nrec2(indg)) &
    -dexp(frail1)*Rrec(indg)-dexp(frail2)*Rrec2(indg) &
    -dexp(frail1*alpha1+frail2*alpha2)*Rdc(indg) &
    +(2.d0*((2.d0*dexp(alpha)/(dexp(alpha)+1.d0))-1.d0) &
    *frail1*frail2/sqrt(theta*eta)  &
    -(frail1**2.d0)/theta -(frail2**2.d0)/eta) &
    /(2.d0*(1.d0-((2.d0*dexp(alpha)/(dexp(alpha)+1.d0))-1.d0)**2.d0)) 

    return
    
    end function funcpamultires


    
      
