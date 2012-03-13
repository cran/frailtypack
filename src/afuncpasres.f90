	
	module commun
	implicit none
	integer,save::ngexact,nssgexact
	integer,dimension(:,:),allocatable,save::ssg
	integer,dimension(:),allocatable,save:: mid 
	integer,dimension(:,:),allocatable,save::mij 
	integer,save::nbpara
	double precision,dimension(:,:),allocatable,save::aux1,aux2
	end module commun 

!========================    FUNCPARES RESIDUS MATRINGALE DENSITE A POSTERIORI       ====================

!!!!
!!!! Calcul Residus shared
!!!!
	double precision function funcpasres(uu,np,id,thi,jd,thj,k0)     
 
	use comon
        use residusM
	
	implicit none

	integer,intent(in)::id,jd,np
	double precision,dimension(np)::bh
	double precision,dimension(np),intent(in)::uu
	double precision,intent(in)::thi,thj
	double precision,dimension(2),intent(in)::k0
	double precision::frail1,kapp0

	bh=uu
	kapp0=k0(1)
	if (id.ne.0) bh(id)=bh(id)+thi
	if (jd.ne.0) bh(jd)=bh(jd)+thj    

	frail1=bh(1)*bh(1)

	funcpasres = frail1**(nig(indg) + 1.d0/theta - 1.d0) &
	*dexp(-frail1*(1.d0/theta + cumulhaz(indg)))
	
	return
		
	end function funcpasres
	
	
	

!!!!
!!!! Calcul Residus joint
!!!!

	double precision function funcpajres(uu,np,id,thi,jd,thj,k0)

	use comon
        use residusM
	
	implicit none

	integer,intent(in)::id,jd,np
	double precision,intent(in)::thi,thj
	double precision,dimension(np)::uu,bh
	double precision,dimension(2),intent(in)::k0
	double precision::frail1,kapp0
	double precision,parameter::pi=3.141592653589793d0

	bh=uu
	kapp0=k0(1)
	if (id.ne.0) bh(id)=bh(id)+thi
	if (jd.ne.0) bh(jd)=bh(jd)+thj    
	
	frail1=bh(1)*bh(1)

!---------- calcul de la penalisation -------------------

	funcpajres = frail1**(Ndc(indg) + Nrec(indg) + 1.d0/theta - 1.d0 + &
	alpha * (Nrec(indg) + Ndc(indg))) * dexp(-frail1*(1/theta + &
	Rrec(indg))) * dexp(-(frail1**alpha)*Rdc(indg))

	return
	
	end function funcpajres
	  


!!!!
!!!! Calcul Residus nested
!!!!

	double precision function funcpanres(uu,np,id,thi,jd,thj,k0)
	
	use comon,only:alpha,eta
        use residusM
	use commun

	
	implicit none

	integer,intent(in)::id,jd,np
	double precision,intent(in)::thi,thj
	double precision,dimension(2),intent(in)::k0	
	double precision,dimension(np),intent(in)::uu
	integer::j
	double precision,dimension(np)::bh
	double precision::frail1,prod1,prod2
	double precision,dimension(np-1)::frail2
	double precision,parameter::pi=3.141592653589793d0
	double precision::kapa

	bh=uu
	kapa=k0(1)
	
	
	if (id.ne.0) bh(id)=bh(id)+thi
	if (jd.ne.0) bh(jd)=bh(jd)+thj    
	
	frail1=bh(1)*bh(1)
	
	do j=1,n_ssgbygrp(indg)
		frail2(j)=bh(j+1)*bh(j+1)
	end do

	prod1 = 1.d0
	prod2 = 1.d0
	

	do j=1,n_ssgbygrp(indg)
		prod1 = prod1 * (frail2(j)**mij(indg,j)) * dexp(-frail1 * frail2(j) * cumulhaz1(indg,j))
		prod2 = prod2 * frail2(j)**((1.d0/eta) - 1) * dexp(-frail2(j)/eta)
	end do

	funcpanres = frail1**(mid(indg)+1.d0/alpha - 1) * prod1 * dexp(-frail1/alpha) * prod2
	
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
!	kapa=k0(1)
		
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
	  