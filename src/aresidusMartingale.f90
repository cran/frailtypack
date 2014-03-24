!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Shared
!=============================================================================
    subroutine ResidusMartingale(b,np,namesfuncres,Resmartingale,frailtypred,frailtyvar,frailtysd)

    use residusM
    use optimres
    use comon

    implicit none
    
    integer::np
    double precision,external::namesfuncres
    double precision,dimension(np),intent(in)::b
    double precision,dimension(ng),intent(out)::Resmartingale
    double precision,dimension(ng),intent(out)::frailtypred,frailtysd,frailtyvar

    
    vecuiRes=0.d0
    moyuiR=0.d0
    varuiR=0.d0
    cares=0.d0
    cbres=0.d0
    ddres=0.d0

! la prediction des effets aleatoires n'est pas la meme pour gamma ou log-normal
    if (logNormal.eq.0) then
        do indg=1,ng
            post_esp(indg)=(nig(indg)+1/(b(np-nva)*b(np-nva)))/(cumulhaz(indg)+1/(b(np-nva)*b(np-nva)))
            post_SD(indg)=dsqrt((nig(indg)+1/(b(np-nva)*b(np-nva)))/((cumulhaz(indg)+1/(b(np-nva)*b(np-nva)))**2))

            Resmartingale(indg)=nig(indg)-(post_esp(indg))*cumulhaz(indg)

            frailtypred(indg) = post_esp(indg)
            frailtysd(indg) = post_SD(indg)
            frailtyvar(indg) = frailtysd(indg)**2
        end do
    else
        do indg=1,ng
            vuu=0.9d0
            call marq98res(vuu,1,nires,vres,rlres,ierres,istopres,cares,cbres,ddres,namesfuncres)

            if (istopres.eq.1) then
                Resmartingale(indg)=nig(indg)-((vuu(1)*vuu(1)))*cumulhaz(indg)
                frailtypred(indg) = vuu(1)*vuu(1)
                frailtyvar(indg) = ((2.d0*vuu(1))**2)*vres(1)
                frailtysd(indg) = dsqrt(frailtyvar(indg))
            else
                ! non convergence ou erreur de calcul de la fonction a maximiser
                Resmartingale(indg) = 0.d0
                frailtypred(indg) = 0.d0
                frailtyvar(indg) = 0.d0
                frailtysd(indg) = 0.d0
            endif
        end do
    endif

    end subroutine ResidusMartingale


!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Joint
!=============================================================================

    subroutine ResidusMartingalej(b,np,namesfuncres,Resmartingale,Resmartingaledc,&
    frailtypred,frailtyvar)

    use residusM
    use optimres
    use comon

    implicit none
    
    integer::np
    double precision,external::namesfuncres
    double precision,dimension(np),intent(in)::b
    double precision,dimension(np)::bint
    double precision,dimension(ng),intent(out)::Resmartingale,Resmartingaledc
    double precision,dimension(ng),intent(out)::frailtypred,frailtyvar
    
    bint=b
    ResidusRec=0.d0
    Residusdc=0.d0
    vecuiRes=0.d0
    moyuiR=0.d0
    
    do indg=1,ng

        vuu=0.9d0

        call marq98res(vuu,1,nires,vres,rlres,ierres,istopres,cares,cbres,ddres,namesfuncres)
        !print*,"6:",vres
        !print*,"----------",indg
        if (istopres.eq.1) then
            ResidusRec(indg)=Nrec(indg)-((vuu(1)*vuu(1)))*Rrec(indg)
            Residusdc(indg)=Ndc(indg)-((vuu(1)*vuu(1))**alpha)*Rdc(indg)
            vecuiRes(indg) = vuu(1)*vuu(1)

            Resmartingale(indg) = ResidusRec(indg)
            Resmartingaledc(indg) = Residusdc(indg)

            frailtypred(indg) = vecuiRes(indg)

            frailtyvar(indg) = ((2.d0*vuu(1))**2)*vres(1)
        else
            ! non convergence ou erreur de calcul de la fonction a maximiser
            Resmartingale(indg) = 0.d0
            Resmartingaledc(indg) = 0.d0
            frailtypred(indg) = 0.d0
            frailtyvar(indg) = 0.d0
        endif

    end do
    
    end subroutine ResidusMartingalej

!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Nested
!=============================================================================
    
    subroutine ResidusMartingalen(namesfuncres,Resmartingale,frailtypred,maxng,frailtypredg,&
    frailtyvar,frailtyvarg,frailtysd,frailtysdg)

    use residusM
    use optimres
    use comon,only:alpha,eta
    use commun

    implicit none
    
    integer::i,j,maxng
    double precision,external::namesfuncres
    double precision,dimension(ngexact),intent(out)::Resmartingale
    double precision,dimension(ngexact),intent(out)::frailtypred,frailtysd,frailtyvar
    double precision,dimension(ngexact,maxng),intent(out)::frailtypredg,frailtysdg,frailtyvarg
    double precision,dimension(:),allocatable::vuuu
    double precision,dimension(:,:),allocatable::H_hess0

    cares=0.d0
    cbres=0.d0
    ddres=0.d0
    Resmartingale = mid

    do indg=1,ngexact

        allocate(H_hess0(n_ssgbygrp(indg)+1,n_ssgbygrp(indg)+1))
        
        allocate(vuuu(n_ssgbygrp(indg)+1),vres((n_ssgbygrp(indg)+1)*((n_ssgbygrp(indg)+1)+3)/2))

        vuuu=0.9d0

        call marq98res(vuuu,(n_ssgbygrp(indg)+1),nires,vres,rlres,ierres,istopres,cares,cbres,ddres,namesfuncres)
        !print*,indg,istopres

        do i=1,n_ssgbygrp(indg)+1
            do j=i,n_ssgbygrp(indg)+1
                H_hess0(i,j)=vres((j-1)*j/2+i)
            end do
        end do
        do i=1,(n_ssgbygrp(indg)+1)
            do j=1,i-1
                H_hess0(i,j) = H_hess0(j,i)
            end do
        end do

        if (istopres.eq.1) then

            do i=1,n_ssgbygrp(indg)
                Resmartingale(indg) = Resmartingale(indg) - ((vuuu(1)*vuuu(1+i))**2)*cumulhaz1(indg,i)
                frailtypredg(indg,i) = vuuu(1+i)**2
            end do

            frailtypred(indg) = vuuu(1)**2

            frailtysd(indg) = dsqrt(((2.d0*vuuu(1))**2)*H_hess0(1,1)) ! correction de la variance le 2 est dans le carre
            frailtyvar(indg) = ((2.d0*vuuu(1))**2)*H_hess0(1,1)

            do i=1,n_ssgbygrp(indg)
                frailtysdg(indg,i) = dsqrt(((2.d0*vuuu(1+i))**2)*H_hess0(1+i,1+i))
                frailtyvarg(indg,i) = ((2.d0*vuuu(1+i))**2)*H_hess0(1+i,1+i)
            end do

        else
            Resmartingale(indg) = 0.d0
            frailtypredg(indg,:) = 0.d0
            frailtysdg(indg,:) = 0.d0
            frailtyvarg(indg,:) = 0.d0
            frailtysd(indg) = 0.d0
            frailtyvar(indg) = 0.d0
        end if

        deallocate(vuuu,vres,H_hess0)!,I_hess,H_hess)

    end do

    end subroutine ResidusMartingalen
    


!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Additive
!=============================================================================
    
    subroutine ResidusMartingalea(b,np,namesfuncres,Resmartingale,frailtypred,frailtyvar,frailtysd,&
    frailtypred2,frailtyvar2,frailtysd2,frailtycov)

    use parameters
    use residusM,only:indg,cumulhaz
    use optimres
    use comon,only:alpha,eta,nst,nig,nsujet,g,stra,nt1,nva,ve,typeof!,H_hess
    use additiv,only:ve2,ngexact,ut1,ut2,mid

    implicit none
    
    integer::np,k,ip,i,j,ier,istop,ni
    double precision::vet,ca,cb,dd,rl
    double precision,dimension(np),intent(in)::b
    double precision,external::namesfuncres    
    double precision,dimension(ngexact),intent(out)::Resmartingale
    double precision,dimension(ngexact),intent(out)::frailtypred,frailtysd,frailtyvar,frailtycov
    double precision,dimension(ngexact),intent(out)::frailtypred2,frailtysd2,frailtyvar2
    double precision,dimension(2,2)::H_hess0
    double precision,dimension(2)::vu
    double precision,dimension(2*(2+3)/2)::v

    vet=0.d0
    H_hess0=0.d0
    
    Resmartingale = mid

    do indg = 1,ngexact

        vu=0.0d0
        v=0.d0
        call marq98res(vu,2,ni,v,rl,ier,istop,ca,cb,dd,namesfuncres)

         do i=1,2
             do j=i,2
                 H_hess0(i,j)=v((j-1)*j/2+i)
             end do
         end do
         H_hess0(2,1) = H_hess0(1,2)
        
        

        do k=1,nsujet
            if(nva.gt.0 .and.g(k).eq.indg)then
                vet = 0.d0 
                do ip = 1,nva
                    vet = vet + b(np-nva +ip)*ve(k,ip)
                end do

                vet = dexp(vet)
            else
                vet=1.d0
            endif
            if(typeof==0) then
                if(g(k) == indg)then    
                    if(stra(k).eq.1)then
                        Resmartingale(indg) = Resmartingale(indg) - ut1(nt1(k)) * &
                            dexp(vu(1) + vu(2) * ve2(k,1) + dlog(vet))
                    end if
                    if(stra(k).eq.2)then
                        Resmartingale(indg) = Resmartingale(indg) - ut2(nt1(k)) * dexp(vu(1) &
                        + vu(2) * ve2(k,1) + dlog(vet))
                    end if
                end if
            else
                if(g(k) == indg)then    
                    Resmartingale(indg) = Resmartingale(indg) - cumulhaz(g(k)) * &
                    dexp(vu(1) + vu(2) * ve2(k,1) + dlog(vet))
                end if
            end if
        end do
    
        frailtypred(indg) = vu(1)
        frailtypred2(indg) = vu(2)

        if(istop==1) then
            frailtyvar(indg) = H_hess0(1,1)
            frailtysd(indg) = dsqrt(H_hess0(1,1))
            
            frailtyvar2(indg) = H_hess0(2,2)
            frailtysd2(indg) = dsqrt(H_hess0(2,2))
            
            frailtycov(indg) = H_hess0(1,2)
            
        else
            frailtysd(indg) = 0.d0
            frailtyvar(indg) = 0.d0 
            frailtysd2(indg) = 0.d0
            frailtyvar2(indg) = 0.d0 
            frailtycov(indg) = 0.d0
        end if
    end do

    end subroutine ResidusMartingalea
    



!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Joint multive
!=============================================================================
        
    subroutine Residus_Martingale_multive(b,np,names_func_res,Res_martingale,Res_martingaledc,Res_martingale2,&
    frailtypred,frailtypred2,frailtyvar,frailtyvar2,frailtyCorr)
    

    use residusMmultiv
!    use optim
    use optimres
    use comonmultiv

    implicit none
    
    integer::np
    double precision,external::names_func_res
    double precision,dimension(np),intent(in)::b
    double precision,dimension(np)::bint
    double precision,dimension(ng),intent(out)::Res_martingale,Res_martingaledc,Res_martingale2
    double precision,dimension(ng),intent(out)::frailtypred,frailtypred2,frailtyvar,frailtyvar2,frailtyCorr    
!    double precision,dimension(3)::k0
    double precision::ca,cb,dd,rl
    integer::ni,ier
    
    bint=b
    ResidusRec=0.d0
    Residusdc=0.d0
    ResidusRec2=0.d0!Residusmeta
    moyuiR=0.d0
!     open(15,file='NNNNN.txt',status='old')
!     open(16,file='NNNNNres.txt',status='old')    
    do indg=1,ng
!        print*,'theta',theta,eta,alpha1,alpha2,alpha
!        k0=0.d0
        ca=0.d0
        cb=0.d0
        dd=0.d0
        ni=0
        vuu=0.1d0
        
!        write(15,*)indg,Nrec(indg),Rrec(indg),Nrec2(indg),Rrec2(indg),Ndc(indg),Rdc(indg)

        
!        call marq98j(k0,vuu,2,ni,vres,rl,ier,istopres,1,ca,cb,dd,names_func_res)
        call marq98res(vuu,2,ni,vres,rl,ier,istopres,ca,cb,dd,names_func_res)
!        write(*,*)'groupe',indg,'ier',ier,'istopres',istopres
        
        ResidusRec(indg)=Nrec(indg)-dexp(vuu(1))**Rrec(indg)
        Residusdc(indg)=Ndc(indg)-dexp(vuu(1)*alpha1+vuu(2)*alpha2)*Rdc(indg)
        ResidusRec2(indg)=Nrec2(indg)-dexp(vuu(2))*Rrec2(indg)
        
        Res_martingale(indg) = ResidusRec(indg)
        Res_martingaledc(indg) = Residusdc(indg)
        Res_martingale2(indg) = ResidusRec2(indg)    
            
!        write(16,*)indg,Res_martingale(indg),Res_martingaledc(indg),Res_martingale2(indg),vuu(1),vuu(2)
        
        frailtypred(indg) = vuu(1)
        frailtypred2(indg) = vuu(2)        
        

        frailtyvar(indg) = vres(1)
        frailtyvar2(indg) = vres(3)    
        frailtyCorr(indg) = vres(2)/dsqrt(vres(1)*vres(2))
    end do    
!     close(15)
!     close(16)    
    end subroutine Residus_Martingale_multive    



