!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Shared
!=============================================================================
	subroutine Residus_Martingale(b,np,names_func_res,Res_martingale,frailtypred,frailtyvar,frailtysd)

	use residusM
	use optim
	use comon

	implicit none
	
	integer::np
	double precision,external::names_func_res
	double precision,dimension(np),intent(in)::b	
	double precision,dimension(ng),intent(out)::Res_martingale
	double precision,dimension(ng),intent(out)::frailtypred,frailtysd,frailtyvar

	
	vecuiRes=0.d0
	moyuiR=0.d0
	varuiR=0.d0
	cares=0.d0
	cbres=0.d0
	ddres=0.d0
	
	
	do indg=1,ng 
		post_esp(indg)=(nig(indg)+1/(b(np-nva)*b(np-nva)))/(cumulhaz(indg)+1/(b(np-nva)*b(np-nva)))
		
		post_SD(indg)=dsqrt((nig(indg)+1/(b(np-nva)*b(np-nva)))/((cumulhaz(indg)+1/(b(np-nva)*b(np-nva)))**2))
		
		Res_martingale(indg)=nig(indg)-(post_esp(indg))*cumulhaz(indg)
		
		frailtypred(indg) = post_esp(indg)
		
		frailtysd(indg) = post_SD(indg)
		
		frailtyvar(indg) = frailtysd(indg)**2
	end do

	
	end subroutine Residus_Martingale
	

!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Joint
!=============================================================================
		
	subroutine Residus_Martingalej(b,np,names_func_res,Res_martingale,Res_martingaledc,&
	frailtypred,frailtyvar)

	use residusM
	use optim2
	use comon

	implicit none
	
	integer::np
	double precision,external::names_func_res
	double precision,dimension(np),intent(in)::b
	double precision,dimension(np)::bint
	double precision,dimension(ng),intent(out)::Res_martingale,Res_martingaledc
	double precision,dimension(ng),intent(out)::frailtypred,frailtyvar	
	
	bint=b
	ResidusRec=0.d0
	Residusdc=0.d0
	vecuiRes=0.d0
	moyuiR=0.d0
	
	do indg=1,ng

		vuu=0.9d0
		call marq98(vuu,1,nires,vres,rlres,ierres,istopres,cares,cbres,ddres,names_func_res)
		ResidusRec(indg)=Nrec(indg)-((vuu(1)*vuu(1)))*Rrec(indg)
		Residusdc(indg)=Ndc(indg)-((vuu(1)*vuu(1))**alpha)*Rdc(indg)
		vecuiRes(indg) = vuu(1)*vuu(1)
		
		Res_martingale(indg) = ResidusRec(indg)
		Res_martingaledc(indg) = Residusdc(indg)

		frailtypred(indg) = vecuiRes(indg)

		frailtyvar(indg) = ((2.d0*vuu(1))**2)*vres(1)
		
	end do	
	
	end subroutine Residus_Martingalej	

!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Nested
!=============================================================================
	
	subroutine Residus_Martingalen(names_func_res,Res_martingale,frailtypred,maxng,frailtypredg,&
	frailtyvar,frailtyvarg,frailtysd,frailtysdg)

	use residusM
	use optim2
	use comon,only:alpha,eta,H_hess,I_hess
	use commun

	implicit none
	
	integer::i,j,maxng
	double precision,external::names_func_res	
	double precision,dimension(ngexact),intent(out)::Res_martingale
	double precision,dimension(ngexact),intent(out)::frailtypred,frailtysd,frailtyvar
	double precision,dimension(ngexact,maxng),intent(out)::frailtypredg,frailtysdg,frailtyvarg
	double precision,dimension(:),allocatable::vuuu
	double precision,dimension(:,:),allocatable::H_hess0

	cares=0.d0
	cbres=0.d0
	ddres=0.d0
	
	Res_martingale = mid 

	do indg=1,ngexact 
		allocate(H_hess0(n_ssgbygrp(indg)+1,n_ssgbygrp(indg)+1))
		
		allocate(vuuu(n_ssgbygrp(indg)+1),vres((n_ssgbygrp(indg)+1)*((n_ssgbygrp(indg)+1)+3)/2),&
		I_hess((n_ssgbygrp(indg)+1),(n_ssgbygrp(indg)+1)),H_hess((n_ssgbygrp(indg)+1),(n_ssgbygrp(indg)+1)))
		
		vuuu=0.9d0
		
		call marq98(vuuu,(n_ssgbygrp(indg)+1),nires,vres,rlres,ierres,istopres,cares,cbres,ddres,names_func_res)

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
		
		do i=1,n_ssgbygrp(indg)
			Res_martingale(indg) = Res_martingale(indg) - ((vuuu(1)*vuuu(1+i))**2)*cumulhaz1(indg,i)
			frailtypredg(indg,i) = vuuu(1+i)**2
		end do
	
		frailtypred(indg) = vuuu(1)**2

		if(istopres==1) then
			
			frailtysd(indg) = dsqrt((2.d0*vuuu(1)**2)*H_hess0(1,1))
			frailtyvar(indg) = 2.d0*(vuuu(1)**2)*H_hess0(1,1)

			do i=1,n_ssgbygrp(indg)
				frailtysdg(indg,i) = dsqrt((2.d0*vuuu(1+i)**2)*H_hess0(1+i,1+i))
				frailtyvarg(indg,i) = 2.d0*(vuuu(1+i)**2)*H_hess0(1+i,1+i)
			end do
		else
			frailtysdg(indg,:) = 0.d0
			frailtyvarg(indg,:) = 0.d0 
			frailtysd(indg) = 0.d0
			frailtyvar(indg) = 0.d0 
		end if
		deallocate(vuuu,vres,I_hess,H_hess,H_hess0)
	end do
	
	end subroutine Residus_Martingalen
	


!=============================================================================
!                       CALCUL DES RESIDUS de MARTINGALES Additive
!=============================================================================
	
	subroutine Residus_Martingalea(b,np,names_func_res,Res_martingale,frailtypred,frailtyvar,frailtysd,&
	frailtypred2,frailtyvar2,frailtysd2,frailtycov)

	use parameters
	use residusM
	use optim2
	use optim
	use comon,only:alpha,eta,nst,nig,nsujet,g,stra,nt1,nva,ve,H_hess
	use additiv,only:ve2,ngexact,ut1,ut2,mid

	implicit none
	
	integer::np,k,ip,i,j
	double precision::vet
	double precision,dimension(np),intent(in)::b
	double precision,external::names_func_res	
	double precision,dimension(ngexact),intent(out)::Res_martingale
	double precision,dimension(ngexact),intent(out)::frailtypred,frailtysd,frailtyvar,frailtycov
	double precision,dimension(ngexact),intent(out)::frailtypred2,frailtysd2,frailtyvar2
	double precision,dimension(2,2)::H_hess0
	
	cares=0.d0
	cbres=0.d0
	ddres=0.d0
	vet=0.d0
	H_hess0=0.d0
	
	Res_martingale = mid


	do indg = 1,ngexact

		vuu=0.9d0
		vres=0.d0
		effetres=0
		call marq98(vuu,2,nires,vres,rlres,ierres,istopres,cares,cbres,ddres,names_func_res)
 		do i=1,2
 			do j=i,2
 				H_hess0(i,j)=vres((j-1)*j/2+i)
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
			if(g(k) == indg)then	
				if(stra(k).eq.1)then
					Res_martingale(indg) = Res_martingale(indg) - ut1(nt1(k)) * dexp(vuu(1) + vuu(2) * ve2(k,1) + dlog(vet))
				end if
				if(stra(k).eq.2)then
					Res_martingale(indg) = Res_martingale(indg) - ut2(nt1(k)) * dexp(vuu(1) + vuu(2) * ve2(k,1) + dlog(vet))
				end if
			end if
		end do
	
		frailtypred(indg) = vuu(1)
		frailtypred2(indg) = vuu(2)

		if(istopres==1) then
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
	
	end subroutine Residus_Martingalea
	