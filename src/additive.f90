	

	
	module propred	
		integer,dimension(:),allocatable,save::filtre,filtre2      
		integer,dimension(:),allocatable,save::stracross!pour crossvalidation   
		double precision,dimension(:),allocatable,save::aux
		double precision,dimension(:,:),allocatable,save::y
		double precision,dimension(:),allocatable,save::v		
		double precision,dimension(:,:),allocatable,save::I1_hess,H1_hess
		double precision,dimension(:,:),allocatable,save::I2_hess,H2_hess
		double precision,dimension(:,:),allocatable,save::HI1,HI2
		double precision,dimension(:,:),allocatable,save::HIH,IH,HI	
	end module propred
	
!
!
!#########################################################################################################
!
!
	subroutine additive(ns0,ng0,nst0, &
	nz0,xmin10,xmin20,tt00,tt10,ic0,groupe0,nva0, &
	str0,vax0,interaction,ag0,noVar,maxiter0,irep10, &
	correl0,np,b,&
	coef,varcoef,varcoef2,rhoEnd,covEnd,varcovEnd, &
	varSigma2,varTau2, &
	ni,res,traceLCV,k0,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out,ier,ddl)

	use parameters		
	use tailles
	use propred
	use comon
	use additiv
	use optim
	
	implicit none	
     
	integer::noVar,interaction
	integer::groupe,j,k,nz,n,np,cpt,ii,iii,iiii,ver, &
	cptstr1,cptstr2,i,ic,ni,ier,istop,l,str, &
        nb_echec,nb_echecor,auxng, &
	irep1,nvacross,nstcross,effetcross
!----------------- ajout      
	integer::ns0,ng0,nst0,correl0,maxiter0,nva0,ag0,nz0,irep10
	double precision,dimension(ns0)::tt00,tt10
	integer,dimension(ns0)::groupe0,ic0,str0
	double precision, dimension(ns0,nva0)::vax0
	double precision::xmin10,xmin20
	double precision,dimension(np)::b
!------------------
	double precision::tt0,tt1
	double precision::h,xmin1,xmin2,res,min,max,maxt, &
	lrs,trace,trace1,trace2,auxi,ax,bx,cx,tol,ddl, &
	fa,fb,fc,goldenadd,estimvadd,f1,f2,f3,varcov    
	double precision,dimension(2):: auxkappa,k0 
	real,dimension(:),allocatable::vax
	double precision::rhoEnd,covEnd,varcovEnd
	double precision,dimension(nva0)::coef,varcoef,varcoef2	
	double precision,dimension(2)::varSigma2,varTau2
	double precision,dimension(99)::x1Out,x2Out
	double precision,dimension(99,3)::lamOut,suOut,lam2Out,su2Out
!AD: add traceLCV	
	double precision,intent(out)::traceLCV
!AD: add for new marq
	double precision::ca,cb,dd,funcpaadd
	external::funcpaadd
	ca=0.d0
	cb=0.d0
	dd=0.d0	
!AD:end	
	epsa=1.d-4
	epsb=1.d-4
	epsd=1.d-3
	model=2
	auxng=0
	maxiter=maxiter0
!----------------
	str=0
!---------------- 
	allocate(invd(2,2))
	 
!=== initialisations 
        lrs=0.d0    
        nb_echec=0 
        nb_echecor=0

!=== fin initialisations 
    
	xmin2=xmin20
	xmin1=xmin10

!-----------------------------------------------------------------------------------------
	
	nsujet=ns0
	nsujetmax=nsujet
	allocate(t0(nsujetmax),t1(nsujetmax),c(nsujetmax),nt0(nsujetmax),nt1(nsujetmax), &
	stra(nsujetmax),g(nsujetmax),stracross(nsujetmax),aux(2*nsujetmax))

	ndatemax=2*nsujet
	allocate(date(ndatemax),mm3(ndatemax),mm2(ndatemax),mm1(ndatemax),mm(ndatemax), &
	im3(ndatemax),im2(ndatemax),im1(ndatemax),im(ndatemax),dut1(ndatemax),dut2(ndatemax), &
	ut1(0:ndatemax),ut2(0:ndatemax))
!-----------------------------------------------------------------------------------------
	ng=ng0
	ngmax=ng
	allocate(nig(ngmax),mid(ngmax))	
!-----------------------------------------------------------------------------------------
	nst=nst0
	correl=correl0
	  
	correlini=correl
!-----------------------------------------------------------------------------------------
	ver=nva0
	nvarmax=ver
	allocate(ve(nsujetmax,nvarmax),ve2(nsujetmax,nvarmax),betaaux(nvarmax))
	allocate(vax(nvarmax))
	
	allocate(filtre(nva0),filtre2(nva0))
	if (noVar.eq.1) then 
		do i=1,nva0
			filtre(i)=0
		enddo  
		nva=0  
	else
		do i=1,nva0
			filtre(i)=1
			filtre2(i)=0
		enddo  
		filtre2(interaction)=1
		nva=nva0  
	end if  
	
			
	ag=ag0
	

	nz=nz0
	
	effet = 1
	nig = 0
	g=0

    
!------------  lecture fichier -----------------------

	maxt = 0.d0
	cpt = 0
	k = 0
	cptstr1 = 0
	cptstr2 = 0

	do i = 1,nsujet    
		if(nst.eq.2)then
			tt0=tt00(i)
			tt1=tt10(i)
			groupe=groupe0(i)
			ic=ic0(i)
			str=str0(i)
			do j=1,nva
				vax(j)=vax0(i,j)
			end do
		else
			tt0=tt00(i)
			tt1=tt10(i)
			groupe=groupe0(i)
			ic=ic0(i)
			do j=1,nva
				vax(j)=vax0(i,j)
			end do			
		endif
		k = k +1
		if(k.eq.1)then
			auxng=groupe
			ngexact=1
			g(k)=1
		else
			g(k)=groupe
		endif
!------------------   observation c=1
		if(ic.eq.1)then
			cpt = cpt + 1
			c(k)=1
			if(str.eq.1.and.nst.eq.2)then
				stra(k) = 1
				cptstr1 = cptstr1 + 1
			endif
			if(str.eq.2.and.nst.eq.2)then
				stra(k) = 2
				cptstr2 = cptstr2 + 1
			endif

			if(nst.eq.1)then
				stra(k) = 1
				cptstr1 = cptstr1 + 1
			endif
			t0(k) = tt0
			t1(k) = tt1

			if(auxng.ne.groupe)then !chgt de groupes
				ngexact=ngexact+1
				auxng=groupe

				if(k.ne.1)then
					g(k)=g(k-1)+1
				endif
				goto 100             
			endif
  
100  continue

			nig(g(k)) = nig(g(k))+1
			iii = 0
			iiii = 0
			do ii = 1,ver
				if(filtre(ii).eq.1)then
					iii = iii + 1
					ve(i,iii) = dble(vax(ii))
					ve2(i,iii) = ve(i,iii) 
!====================================================================
				endif
				if(filtre2(ii).eq.1)then
					iiii = iiii + 1
					ve(k,iiii) = dble(vax(ii))
					ve2(k,iiii) = ve(k,iiii)  
				endif
!====================================================================
			end do   
		else 
!------------------   censure a droite  c=0
			if(ic.eq.0)then
				c(k) = 0 
				if(str.eq.1.and.nst.eq.2)then
					stra(k) = 1
					cptstr1 = cptstr1 + 1
				endif
				if(str.eq.2.and.nst.eq.2)then
					stra(k) = 2
					cptstr2 = cptstr2 + 1
				endif
				if(nst.eq.1)then
					stra(k) = 1
					cptstr1 = cptstr1 + 1
				endif
				iii = 0
				iiii = 0
				do ii = 1,ver
					if(filtre(ii).eq.1)then
						iii = iii + 1
						ve(k,iii) = dble(vax(ii))
						ve2(k,iii) =  ve(k,iii) 
!==================recodage en +-1/2 ===============================
					endif
					if(filtre2(ii).eq.1)then
						iiii = iiii + 1
						ve(k,iiii) = dble(vax(ii))
						ve2(k,iiii) =  ve(k,iiii) 
					endif
				end do 

				t0(k) =  tt0
				t1(k) = tt1
				if(auxng.ne.groupe)then !chgt de groupes
					ngexact=ngexact+1
					auxng=groupe
					if(k.ne.1)then
						g(k)=g(k-1)+1
					endif
					goto 101
				endif   
101   continue
				nig(g(k)) = nig(g(k))+1
			endif
		endif
		if (maxt.lt.t1(k))then
			maxt = t1(k)
		endif
	end do 
	
	nz1=nz
	nz2=nz
	if(nz.gt.20)then
		nz = 20
	endif
	if(nz.lt.4)then
		nz = 4
	endif

!***************************************************
!--------------- zi- ----------------------------------

!      construire vecteur zi (des noeuds)

	min = 1.d-10
	max = maxt

	do i = 1,2*nsujet
		do k = 1,nsujet
			if((t0(k).ge.min))then
				if(t0(k).lt.max)then
					max = t0(k)
				endif
			endif
			if((t1(k).ge.min))then
				if(t1(k).lt.max)then
					max = t1(k)
				endif
			endif
		end do   
		aux(i) = max
		min = max + 1.d-12
		max = maxt
	end do

	date(1) = aux(1)
	k = 1
	do i=2,2*nsujet
		if(aux(i).gt.aux(i-1))then
			k = k+1
			date(k) = aux(i)
		endif 
	end do 
	
	nzmax=nz+3
	allocate(zi(-2:nzmax))
	
	ndate = k

	zi(-2) = date(1)
	zi(-1) = date(1)
	zi(0) = date(1)
	zi(1) = date(1)
	h = (date(ndate)-date(1))/dble(nz-1)
	do i=2,nz-1
		zi(i) =zi(i-1) + h   
	end do

	zi(nz) = date(ndate)
	zi(nz+1)=zi(nz)
	zi(nz+2)=zi(nz)
	zi(nz+3)=zi(nz)

!---------- affectation nt0,nt1----------------------------

	indictronq=0
	do i=1,nsujet 
		if(t0(i).eq.0.d0)then
			nt0(i) = 0
		endif
		if(t0(i).ne.0.d0)then
			indictronq=1
		endif
		do j=1,ndate
			if(date(j).eq.t0(i))then
				nt0(i)=j
			endif
			if(date(j).eq.t1(i))then
				nt1(i)=j
			endif
		end do
	end do  

!---------- affectation des vecteurs de splines -----------------
	n  = nz+2
	
	call vecspliadd(n,ndate)
	
	
	allocate(m3m3(nzmax),m2m2(nzmax),m1m1(nzmax),mmm(nzmax),m3m2(nzmax), &
	m3m1(nzmax),m3m(nzmax),m2m1(nzmax),m2m(nzmax),m1m(nzmax))
		
	call vecpenadd(n)  
	    
	np = nst*n + nva + correl + 2*effet 

!-----------------------
	npmax=np

	allocate(I_hess(npmax,npmax),H_hess(npmax,npmax),Hspl_hess(npmax,npmax) &
	,hess(npmax,npmax),PEN_deri(npmax,1))
	allocate(y(npmax,npmax),v((npmax*(npmax+3)/2)),&
	I1_hess(npmax,npmax),H1_hess(npmax,npmax),I2_hess(npmax,npmax),H2_hess(npmax,npmax), &
	HI1(npmax,npmax),HI2(npmax,npmax),HIH(npmax,npmax),IH(npmax,npmax),HI(npmax,npmax))
!---------------------	
! term de correlation entre frailty + 2 termes de var de la meme frailty   
	nbpara =np
            
!------- initialisation des parametres                  
	do i=1,np
		b(i)= 1.d-1!5.d-1!
	end do
! traitement, intercept et pente initialis�s plus loin


!***********************************************************
!************** NEW : cross validation  ***********************
!!!!! sur une seule strate, sans var expli , sans frailties ****
!*****************************************************************
	nvacross=nva !pour la recherche du parametre de lissage sans var expli
	nva=0
	effetcross=effet
	effet=0
	nstcross=nst
	nst=1
	correl=0

	do l=1,nsujet  
		stracross(l)=stra(l)
	end do
	do l=1,nsujet  
		stra(l)=1
	end do

	
	irep1=irep10

	if(irep1.eq.0.and.nst.ge.2)then
! ne se produit jamais maintenant (1 seule strate pour CV): changement de juin 2009
		stop
	endif
	
	xmin1=xmin10
		
	if(xmin1.le.0.d0)then
		xmin1 = 0.d0
	endif  

	
!*************************************************
!	nvat=nva
	nva=0


!on  travaille d'abord sur une seule strate

	if(irep1.eq.1)then   !pas recherche du parametre de lissage
		xmin1 = dsqrt(xmin1)

		auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)

		if (ni.ge.250) then

			do i=1,nz+2
				b(i)=1.d-1
			end do     
			xmin1 = sqrt(10.d0)*xmin1
			auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)

			if (ni.lt.250) then

			else
				do i=1,nz+2
					b(i)=1.d-1
				end do     
				xmin1 = sqrt(10.d0)*xmin1
				auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)
 
			endif
		else

		endif
!----------------------------------------------------
	else                   !recherche du parametre de lissage

		if(xmin1.le.0.d0)then
			xmin1 = 1.d0
		endif  

		xmin1 = dsqrt(xmin1)
		auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)

		if(ddl.gt.-2.5d0)then
			xmin1 = dsqrt(xmin1)
			auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)
	
			if(ddl.gt.-2.5d0)then
				xmin1 = dsqrt(xmin1)
				auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)


				if(ddl.gt.-2.5d0)then
					xmin1 = dsqrt(xmin1)
					auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)

					if(ddl.gt.-2.5d0)then
						xmin1 = dsqrt(xmin1)
						auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)

						if(ddl.gt.-2.5d0)then
							xmin1 = dsqrt(xmin1)
						endif   
					endif   
				endif   
			endif
		endif 

		if (ni.ge.250) then
			do i=1,nz+2
				b(i)=1.d-1
			end do     
			xmin1 = sqrt(10.d0)*xmin1
			auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)

			if (ni.ge.250) then
				do i=1,nz+2
					b(i)=1.d-1
				end do     
				xmin1 = sqrt(10.d0)*xmin1
			endif
		endif 
		ax = xmin1
		bx = xmin1*dsqrt(1.5d0)  

		call mnbrakadd(ax,bx,cx,fa,fb,fc,b,n)            
		tol = 0.001d0
		res = goldenadd(ax,bx,cx,tol,xmin1,n,b,y,ddl)

		effet=0
		correl=0

		auxkappa(1)=xmin1*xmin1
		auxkappa(2)=0.d0

		call marq98J(auxkappa,b,n,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaadd)
	

	endif  

	
	nva=nvacross ! pour la recherche des parametres de regression
	nst=nstcross ! avec stratification si n�cessaire
	effet=effetcross ! avec effet initial
	do l=1,nsujet  
		stra(l)=stracross(l) !r�tablissement stratification
	end do
	
!********************************************************************
	if(nst.eq.2)then
		xmin2=xmin20
	endif


	k0(1) = xmin1*xmin1
	k0(2) = 0.d0
	
	if(nst.eq.2)then
		k0(2) = xmin2
	endif

!=============================- fin cross validation
	
!===== initialisation des parametres de regression/pas effets aleatopires
!	write(*,*)'====================================='
!	write(*,*)'== avec var explicatives============='
!	write(*,*)'==================================='
	
	effet=0

	np = nst*n + nva

	b(nst*n+1)=-0.15d0           !initialisation traitementc


!	write(*,*)'===avant marq98==========',n,np,nva,nst,effet
	call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaadd)
!AD:	
!	if (istop.eq.4) goto 1000
!AD:		
!===== recherche de l'ensemble des parametres
	
!	write(*,*)'====================================='
!	write(*,*)'== ensemble des parametres ======='
!	write(*,*)'====================================='
	
	
	effet=1
	correl=correlini
	do i=1,nva
		b(np-i+2+correl+1)=b(np-i+1)
	end do
	np=nbpara 
	
	if(correl.eq.1)then
		b(np-nva-2)=0.1d0!-0.69d0!-1.61d0!-0.64d0!initialisation cov avec contrainte 
	endif
	
	b(np-nva-1)=0.5d0!0.5477d0      !      initialisation intercept
	b(np-nva)=0.15d0!0.5477d0  !0.13d0!  !   0.15d0!   initialisation pente
	
	call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaadd)
!AD:	
!	if (istop.eq.4) goto 1000
!AD:		
	j=(np-nva)*(np-nva+1)/2
	
	trace=0
	trace1=0
	trace2=0
	
	if(indic_sousv.eq.0)then 
! que lorsque la matrice totale est inversible, ie sans echec inversion
! strate1 : 
		do i=1,nz1+2
			do j=1,nz1+2
				H1_hess(i,j)=H_hess(i,j)
				I1_hess(i,j)=I_hess(i,j)
			end do
		end do 
		call multiJ(H1_hess,I1_hess,nz1+2,nz1+2,nz1+2,HI1)
		do i =1,nz1+2
			trace1=trace1+HI1(i,i)
		end do

! strate2 :
		if(nst.eq.2)then
			do i=1,nz2+2
				k=nz1+2+i
				do j=1,nz2+2
					l=nz1+2+j                   
					H2_hess(i,j)=H_hess(k,l)
					I2_hess(i,j)=I_hess(k,l)
				end do
			end do
			call multiJ(H2_hess,I2_hess,nz2+2,nz2+2,nz2+2,HI2)
			do i =1,nz2+2
				trace2=trace2+HI2(i,i)
			end do
		endif ! pour la deuxieme strate
	endif
! fin calcul trace
	
	if(indic_sousv.eq.0)then
		call multiJ(I_hess,H_hess,np,np,np,IH)
		call multiJ(H_hess,IH,np,np,np,HIH)
	endif
	if(indic_sousv.eq.1)then
		call multiJ(I_hess,H_hess,sousm,sousm,sousm,IH)
		call multiJ(H_hess,IH,sousm,sousm,sousm,HIH)
	endif
	
! Salida Juan Aug'07
	
	f1 = (b(np-nva))*(2.d0* dexp(b(np-nva-2))/(1.d0+dexp(b(np-nva-2)))-1.d0)
	f2 = (b(np-nva-1))*(2.d0* dexp(b(np-nva-2))/(1.d0+dexp(b(np-nva-2)))-1.d0)
	f3= 2.d0*(b(np-nva-1))*(b(np-nva))*(dexp(b(np-nva-2))/(1.d0+dexp(b(np-nva-2)))- &
	(dexp(b(np-nva-2))/(1.d0+dexp(b(np-nva-2))))**2)

	varcov=f1*f1*H_hess(np-nva-1,np-nva-1)+f2*f2*H_hess(np-nva,np-nva)+ &
	f3*f3*H_hess(np-nva-2,np-nva-2)+2.d0*f1*f3*H_hess(np-nva-1,np-nva-2)+ &
	2.d0*f2*f3*H_hess(np-nva,np-nva-2)+2.d0*f1*f2*H_hess(np-nva-1,np-nva)
     
	
	if (correl.eq.1)then
		rho=cov/(dsqrt(b(np-nva-1)*b(np-nva-1)*b(np-nva)*b(np-nva)))
	else
		rho=-1
	end if

	
	rhoEnd=rho
	covEnd=cov
	varcovEnd=varcov


	if(indic_sousv.eq.0)then
		varSigma2(1)=((2.d0*b(np-nva-1))**2)*H_hess(np-nva-1,np-nva-1)
		varSigma2(2)=((2.d0*b(np-nva-1))**2)*HIH(np-nva-1,np-nva-1)
		
		varTau2(1)=((2.d0*b(np-nva))**2)*H_hess(np-nva,np-nva)
		varTau2(2)=((2.d0*b(np-nva))**2)*HIH(np-nva,np-nva)
	end if
	
	if(indic_sousv.eq.1)then
		varSigma2(1)=((2.d0*b(np-nva-1))**2)*H_hess(sousm-nva-1,sousm-nva-1)
		varSigma2(2)=((2.d0*b(np-nva-1))**2)*HIH(sousm-nva-1,sousm-nva-1)
		
		varTau2(1)=((2.d0*b(np-nva))**2)*H_hess(sousm-nva,sousm-nva)
		varTau2(2)=((2.d0*b(np-nva))**2)*HIH(sousm-nva,sousm-nva)
	end if  

! Covariates

	do i=1,nva
	
		coef(i)=b(np-nva+i)
		
		if(indic_sousv.eq.0)then
			varcoef(i)=H_hess(np-nva+i,np-nva+i) 
			varcoef2(i)=HIH(np-nva+i,np-nva+i)
		end if     
			
		if(indic_sousv.eq.1)then
			varcoef(i)=H_hess(sousm-nva+i,sousm-nva+i)
			varcoef2(i)=HIH(sousm-nva+i,sousm-nva+i)
		end if 
		
	end do

	
	if(effet.eq.1.and.ier.eq.-1)then
		v((np-nva)*(np-nva+1)/2)=10.d10
	endif
		
	j=(np-nva)*(np-nva+1)/2
	
! --------------  Lambda and survival and cumulative hazard estimates 
	
!	call distance(nz1,nz2,b,fich3,fich4,fich3b,fich4b,fich3c,fich4c,effet)

	call distanceadd(nz1,nz2,b,effet,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)
	
!AD:add LCV
!     calcul de la trace, pour le LCV (likelihood cross validation)
	traceLCV=0.d0
	call multiJ(H_hess,I_hess,np,np,np,HI)
	
	do i =1,np
		traceLCV = traceLCV + HI(i,i)
	end do 
	
	traceLCV = (traceLCV - resnonpen) / nsujet
!AD:end     
!1000    continue
	deallocate(invd)
	deallocate(t0,t1,c,nt0,nt1,stra,g,stracross,aux)
	deallocate(date,mm3,mm2,mm1,mm,im3,im2,im1,im,dut1,dut2,ut1,ut2)
	deallocate(nig,mid)	
	deallocate(ve,ve2,betaaux)
	deallocate(vax)
	deallocate(filtre,filtre2)		
	deallocate(zi)	
	deallocate(m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m)
	deallocate(I_hess,H_hess,Hspl_hess,hess,PEN_deri)
	deallocate(y,v,I1_hess,H1_hess,I2_hess,H2_hess,HI1,HI2,HIH,IH,HI)	
	
	end subroutine additive

	
	
	
	!========================== VECSPLI =====================
	subroutine vecspliadd(n,ndate) 
	
	use tailles
	use comon,only:zi,date,mm3,mm2,mm1,mm,im3,im2,im1,im
	
	implicit none
	
	integer::n,ndate,i,j,k
	double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2
	double precision::h3,h4,h3m,h2n,hn,hh3,hh2
	
	
!----------  calcul de u(ti) ---------------------------
	j=0
	do i=1,ndate-1
		do k = 2,n-2
			if ((date(i).ge.zi(k-1)).and.(date(i).lt.zi(k)))then
				j = k-1
			endif
		end do 
		ht = date(i)-zi(j)
		htm= date(i)-zi(j-1)
		h2t= date(i)-zi(j+2)
		ht2 = zi(j+1)-date(i)
		ht3 = zi(j+3)-date(i)
		hht = date(i)-zi(j-2)
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
		mm3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
		mm2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
		*ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
		mm1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
		h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
		mm(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
		im3(i) = (0.25d0*(date(i)-zi(j-3))*mm3(i))+(0.25d0*hh2 &
		*mm2(i))+(0.25d0*h3m*mm1(i))+(0.25d0*h4*mm(i))
		im2(i) = (0.25d0*hht*mm2(i))+(h3m*mm1(i)*0.25d0)+(h4*mm(i)*0.25d0)
		im1(i) = (htm*mm1(i)*0.25d0)+(h4*mm(i)*0.25d0)
		im(i)  = ht*mm(i)*0.25d0
	end do
	
	end subroutine vecspliadd
	
!========================== VECPEN ==============================
	subroutine vecpenadd(n) 
	
	use tailles
	use comon,only:zi,date,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m
	
	implicit none
	
	integer::n,i
	double precision::h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2, &
	a3,a2,b2,c2,a1,b1,c1,a0,x3,x2,x
	
	do i=1,n-3
		h = zi(i+1)-zi(i)
		hh= zi(i+1)-zi(i-1)
		h2= zi(i+2)-zi(i)
		h3= zi(i+3)-zi(i)
		h4= zi(i+4)-zi(i)
		h3m= zi(i+3)-zi(i-1)
		h2n=zi(i+2)-zi(i-1)
		hn= zi(i+1)-zi(i-2)
		hh3 = zi(i+1)-zi(i-3)
		hh2 = zi(i+2)-zi(i-2)
		a3 = h*hh*hn*hh3
		a2 = hh2*hh*h*hn
		b2 = hh2*h2n*hh*h
		c2 = hh2*h2*h*h2n
		a1 = h3m*h2n*hh*h
		b1 = h3m*h2*h*h2n
		c1 = h3m*h3*h2*h
		a0 = h4*h3*h2*h
		x3 = zi(i+1)*zi(i+1)*zi(i+1)-zi(i)*zi(i)*zi(i)
		x2 = zi(i+1)*zi(i+1)-zi(i)*zi(i)
		x  = zi(i+1)-zi(i)
		m3m3(i) = (192.d0*h/(hh*hn*hh3*hh*hn*hh3))
		m2m2(i) = 64.d0*(((3.d0*x3-(3.d0*x2*(2.d0*zi(i+1)+zi(i-2) &
		))+x*(4.d0*zi(i+1)*zi(i+1)+zi(i-2)*zi(i-2)+4.d0*zi(i+1) &
		*zi(i-2)))/(a2*a2)))
		m2m2(i) = m2m2(i) + 64.d0*(((3.d0*x3-(3.d0*x2*(zi(i+2) &
		+zi(i-1)+zi(i+1)))+x*(zi(i+2)*zi(i+2)+zi(i-1)*zi(i-1) &
		+zi(i+1)*zi(i+1)+2.d0*zi(i+2)*zi(i-1)+2.d0*zi(i+2) &
		*zi(i+1)+2.d0*zi(i-1)*zi(i+1)))/(b2*b2)))
		m2m2(i) = m2m2(i) +64.d0*((3.d0*x3-(3.d0*x2*(2.d0*zi(i+2) &
		+zi(i)))+x*(4.d0*zi(i+2)*zi(i+2)+zi(i)*zi(i)+4.d0*zi(i+2) &
		*zi(i)))/(c2*c2))
		m2m2(i) = m2m2(i) +128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+2) &
		+zi(i-1)+3.d0*zi(i+1)+zi(i-2)))+x*(2.d0*zi(i+1)*zi(i+2) &
		+2.d0*zi(i+1)*zi(i-1)+2.d0*zi(i+1)*zi(i+1)+zi(i-2)*zi(i+2) &
		+zi(i-2)*zi(i-1)+zi(i-2)*zi(i+1)))/(a2*b2))
		m2m2(i) = m2m2(i) + 128.d0*((3.d0*x3-(1.5d0* &
		x2*(2.d0*zi(i+2)+zi(i)+2.d0*zi(i+1)+zi(i-2)))+x* &
		(4.d0*zi(i+1)*zi(i+2)+2.d0*zi(i+1)*zi(i)+2.d0*zi(i-2) &
		*zi(i+2)+zi(i-2)*zi(i)))/(a2*c2))
		m2m2(i) = m2m2(i) + 128.d0*((3.d0*x3-(1.5d0*x2 &
		*(3.d0*zi(i+2)+zi(i)+zi(i-1)+zi(i+1)))+x*(zi(i+2)*zi(i)+ &
		2.d0*zi(i-1)*zi(i+2)+zi(i)*zi(i-1)+2.d0*zi(i+1)*zi(i+2) &
		+zi(i+1)*zi(i)+2.d0*zi(i+2)*zi(i+2)))/(b2*c2))
		m1m1(i) = 64.d0*((3.d0*x3-(3.d0*x2*(2.d0*zi(i-1)+zi(i+1))) &
		+x*(4.d0*zi(i-1)*zi(i-1)+zi(i+1)*zi(i+1)+4.d0*zi(i-1) &
		*zi(i+1)))/(a1*a1))
		m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i-1)+zi(i)   &  
		+zi(i+2)))+x*(zi(i-1)*zi(i-1)+zi(i)*zi(i)+zi(i+2)* &
		zi(i+2)+2.d0*zi(i-1)*zi(i)+2.d0*zi(i-1)*zi(i+2)+2.d0* &
		zi(i)*zi(i+2)))/(b1*b1))
		m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i+3) &
		+2.d0*zi(i)))+x*(zi(i+3)*zi(i+3)+4.d0*zi(i)*zi(i) &
		+4.d0*zi(i+3)*zi(i)))/(c1*c1)) 
		m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(3.d0 &
		*zi(i-1)+zi(i)+zi(i+2)+zi(i+1)))+x*(2.d0*zi(i-1)*zi(i-1) &
		+2.d0*zi(i-1)*zi(i)+2.d0*zi(i-1)*zi(i+2)+zi(i+1)*zi(i-1) &
		+zi(i+1)*zi(i)+zi(i+1)*zi(i+2)))/(a1*b1))
		m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+3)+ &
		2.d0*zi(i)+2.d0*zi(i-1)+zi(i+1)))+x*(2.d0*zi(i-1)*zi(i+3) &
		+4.d0*zi(i-1)*zi(i)+zi(i+1)*zi(i+3)+2.d0*zi(i+1)*zi(i))) &
		/(a1*c1))    
		m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+3)+3.d0 &
		*zi(i)+zi(i-1)+zi(i+2)))+x*(zi(i-1)*zi(i+3)+2.d0*zi(i-1)  &   
		*zi(i)+zi(i+3)*zi(i)+2.d0*zi(i)*zi(i)+zi(i+2)*zi(i+3) &
		+2.d0*zi(i+2)*zi(i)))/(b1*c1))
		mmm(i) = (192.d0*h/(h4*h3*h2*h4*h3*h2))
		m3m2(i) = 192.d0*(((-x3+(0.5d0*x2*(5.d0*zi(i+1)+zi(i-2) &
		))-x*(2.d0*zi(i+1)*zi(i+1)+zi(i+1)*zi(i-2)))/(a3*a2)) &
		+((-x3+(0.5d0*x2*(4.d0*zi(i+1)+zi(i-1)+zi(i+2)))-x* &
		(zi(i+1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*b2)) &
		+((-x3+(0.5d0*x2*(3.d0*zi(i+1)+2.d0*zi(i+2)+zi(i)))-x* &
		(2.d0*zi(i+1)*zi(i+2)+zi(i+1)*zi(i)))/(a3*c2)))
		m3m1(i) = 192.d0*(((x3-(0.5d0*x2*(4.d0*zi(i+1)+2.d0*zi(i-1) &
		))+x*(2.d0*zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*a1)) &
		+((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+2)+zi(i-1)+zi(i))) &
		+x*(zi(i+1)*zi(i-1)+zi(i+1)*zi(i)+zi(i+1)*zi(i+2)))/(b1*a3)) &
		+((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+3)+2.d0*zi(i)))+x*(zi(i+1) &
		*zi(i+3)+2.d0*zi(i+1)*zi(i)))/(c1*a3)) )
		m3m(i) = 576.d0*((-(x3/3.d0)+(0.5d0*x2*(zi(i+1)+zi(i))) &
		-x*zi(i+1)*zi(i))/(a3*a0))
		m2m1(i) = 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0*zi(i-1)+3.d0* &
		zi(i+1)+zi(i-2)))-x*(4.d0*zi(i+1)*zi(i-1)+2.d0*zi(i+1) &
		*zi(i+1)+2.d0*zi(i-2)*zi(i-1)+zi(i-2)*zi(i+1)))/(a2*a1))
		m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i-1)+ &
		zi(i)+zi(i+2)+2.d0*zi(i+1)+zi(i-2)))-x*(2.d0*zi(i+1)*zi(i-1) &
		+2.d0*zi(i+1)*zi(i)+2.d0*zi(i+1)*zi(i+2)+zi(i-2)*zi(i-1)+ &
		zi(i-2)*zi(i)+zi(i-2)*zi(i+2)))/(a2*b1))
		m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3)+2.d0 &
		*zi(i)+2.d0*zi(i+1)+zi(i-2)))-x*(2.d0*zi(i+1)*zi(i+3)+4.d0 &
		*zi(i+1)*zi(i)+zi(i-2)*zi(i+3)+2.d0*zi(i-2)*zi(i)))/(a2*c1))
		m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2* &
		(3.d0*zi(i-1)+2.d0*zi(i+1)+zi(i+2)))-x*(2.d0*zi(i+2)*zi(i-1) &
		+zi(i+2)*zi(i+1)+2.d0*zi(i-1)*zi(i-1)+3.d0 &
		*zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(b2*a1))
		m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0 &
		*zi(i-1)+zi(i)+2.d0*zi(i+2)+zi(i+1)))-x*(zi(i+2)*zi(i-1) &
		+zi(i+2)*zi(i)+zi(i+2)*zi(i+2)+zi(i-1)*zi(i-1)+zi(i-1) &
		*zi(i)+zi(i-1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i) &
		+zi(i+1)*zi(i+2)))/(b2*b1))
		m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3) &
		+2.d0*zi(i)+zi(i+2)+zi(i-1)+zi(i+1)))-x*(zi(i+2)*zi(i+3) &
		+2.d0*zi(i+2)*zi(i)+zi(i-1)*zi(i+3)+2.d0*zi(i-1)*zi(i) &
		+zi(i+1)*zi(i+3)+2.d0*zi(i+1)*zi(i)))/(b2*c1))
		m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0*zi(i-1) &
		+zi(i+1)+2.d0*zi(i+2)+zi(i)))-x*(4.d0*zi(i+2)*zi(i-1)+2.d0* &
		zi(i+2)*zi(i+1)+2.d0*zi(i)*zi(i-1)+zi(i)*zi(i+1)))/(c2*a1))
		m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i-1) &
		+2.d0*zi(i)+3.d0*zi(i+2)))-x*(2.d0*zi(i+2)*zi(i-1)+2.d0 &
		*zi(i+2)*zi(i)+2.d0*zi(i+2)*zi(i+2)+zi(i)*zi(i-1)+zi(i) &
		*zi(i)+zi(i)*zi(i+2)))/(c2*b1))
		m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3) &
		+3.d0*zi(i)+2.d0*zi(i+2)))-x*(2.d0*zi(i+2)*zi(i+3)+4.d0 &
		*zi(i+2)*zi(i)+zi(i)*zi(i+3)+2.d0*zi(i)*zi(i)))/(c2*c1))
		m2m(i) = 192.d0*(((x3-(0.5d0*x2*(3.d0*zi(i)+2.d0*zi(i+1) &
		+zi(i-2)))+x*(2.d0*zi(i+1)*zi(i)+zi(i-2)*zi(i)))/(a2*a0)) &
		+((x3-(0.5d0*x2*(3.d0*zi(i)+zi(i+2)+zi(i-1)+zi(i+1))) &
		+x*(zi(i+2)*zi(i)+zi(i-1)*zi(i)+zi(i+1)*zi(i)))/(b2*a0)) &
		+((x3-(0.5d0*x2*(4.d0*zi(i)+2.d0*zi(i+2)))+x*(2.d0*zi(i+2) &
		*zi(i)+zi(i)*zi(i)))/(c2*a0)) )
		m1m(i) = 192.d0*(((-x3+(0.5d0*x2*(3.d0*zi(i)+2.d0*zi(i-1) &
		+zi(i+1)))-x*(2.d0*zi(i-1)*zi(i)+zi(i+1)*zi(i)))/(a1*a0)) &
		+((-x3+(0.5d0*x2*(4.d0*zi(i)+zi(i-1)+zi(i+2))) &
		-x*(zi(i-1)*zi(i)+zi(i)*zi(i)+zi(i+2)*zi(i)))/(b1*a0)) &
		+((-x3+(0.5d0*x2*(5.d0*zi(i)+zi(i+3)))-x*(zi(i+3)*zi(i) &
		+2.d0*zi(i)*zi(i)))/(c1*a0)) )
	end do
	
	end subroutine vecpenadd


!==========================  DISTANCE   =================================
	
!	subroutine distance(nz1,nz2,b,fic1,fic2,fic1b,fic2b,fic1c,fic2c,effet)

	subroutine distanceadd(nz1,nz2,b,effet,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)	
	 	
	use tailles
	use comon,only:zi,date,t0,t1,c,nt0,nt1,nsujet,nva,ndate,nst, &
	PEN_deri,I_hess,H_hess,Hspl_hess,hess	

	
	implicit none
	
	
	integer::nz1,nz2,i,j,n,np,k,l,effet
	double precision::x1,x2,su,bsup,binf,lam,lbinf, &
	h,lbsup
	double precision,dimension(npmax,npmax)::hes1,hes2
	double precision,dimension(-2:npmax):: the1,the2
	double precision,dimension(npmax):: b
	double precision,dimension(99)::x1Out,x2Out
	double precision,dimension(99,3)::lamOut,suOut,lam2Out,su2Out
	
	n  = nz1+2
	
	if(nst.eq.2)then      
		np  = nz1+2+nz2+2+effet+nva
	else
		np  = nz1+2+effet+nva
	endif
	
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
		
		
	do i=1,99 
		if(i .ne.1)then
			x1 = x1 + h 
		end if
		call cospadd(x1,the1,nz1+2,hes1,zi,binf,su,bsup,lbinf,lam,lbsup)
	
		if(bsup.lt.0.d0)then
			bsup = 0.d0
		endif
		if(binf.gt.1.d0)then
			binf = 1.d0 
		endif
		if(lbinf.lt.0.d0)then
			lbinf = 0.d0
		endif 
!!
!!
!!	
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
			call cospadd(x2,the2,nz2+2,hes2,zi,binf,su,bsup,lbinf,lam,lbsup)
			if(binf.lt.0.d0)then
				binf = 0.d0
			endif
			if(bsup.gt.1.d0)then
				bsup = 1.d0
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
	
	end subroutine distanceadd

!==========================  SUSP  ====================================
	subroutine suspadd(x,the,n,su,lam,zi)
	
	use tailles
	
	implicit none
	integer::j,k,n,i
	double precision::x,ht,ht2,h2,som,lam,su,htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3
	double precision::ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2,h,gl,hh 
	double precision,dimension(-2:npmax)::zi,the
	
	som = 0.d0
	gl=0.d0
	do k = 2,n+1
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
		do i=1,n+1
			som = som+the(i-3)
		end do
		gl = som
	endif
	su  = dexp(-gl)
	
	return
	
	end subroutine suspadd
	
!==========================  COSP  ====================================
	
	subroutine cospadd(x,the,n,y,zi,binf,su,bsup,lbinf,lam,lbsup)
	
	use tailles
	
	implicit none
	
	integer::j,k,n,i
	double precision::x,ht,ht2,h2,som,lam,su,binf,bsup,lbinf,lbsup,pm, &
	htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3,ht3,hht,h4, &
	h3m,hh3,hh2,mm,im3,mm2,h,gl,hh
	double precision,dimension(-2:npmax)::the,zi
	double precision,dimension(npmax,npmax)::y

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
	
	call confadd(x,j,n,y,pm,zi)
	
	binf = dexp(-gl - 1.96d0*pm)
	su  = dexp(-gl)
	bsup = dexp(-gl + 1.96d0*pm)
	
	call conf1add(x,j,n,y,pm,zi)
	lbinf = lam - 1.96d0*pm
	lbsup = lam + 1.96d0*pm
	
	return
	
	end subroutine cospadd
	
!=====================  CONF1  =============================
	subroutine conf1add(x,ni,n,y,pm,zi)
	
	use tailles
	
	implicit none
	
	integer::ni,i,n,j
	double precision::mmspadd,x,pm,res
	double precision,dimension(-2:npmax) :: zi
	double precision,dimension(npmax) :: vecti,aux
	double precision,dimension(npmax,npmax) :: y
	
	do i=1,n
		vecti(i) = mmspadd(x,ni,i,zi)
	end do   
	
	do i=1,n
		aux(i) = 0.d0
		do j=1,n
			aux(i) = aux(i) - y(i,j)*vecti(j)
		end do
	end do   
	
	res = 0.d0
	do i=1,n
		res = res + aux(i)*vecti(i)
	end do
	
	if (res.lt.0)then 
		res=-res
	endif 
	
	pm = dsqrt(res) 
	
	end subroutine conf1add
	

!=====================  CONF  =============================
	subroutine confadd(x,ni,n,y,pm,zi)
	
	use tailles
	
	implicit none
	
	integer::ni,i,n,j
	double precision::ispadd,x,pm,res
	double precision,dimension(-2:npmax) :: zi
	double precision,dimension(52) :: vecti,aux
	double precision,dimension(npmax,npmax) :: y
	
	
	do i=1,n
		vecti(i) = ispadd(x,ni,i,zi)
	end do   
	
	do i=1,n
		aux(i) = 0.d0
		do j=1,n
			aux(i) = aux(i) - y(i,j)*vecti(j)
		end do
	end do   
	
	res = 0.d0
	do i=1,n
		res = res + aux(i)*vecti(i)
	end do
	
	if (res.lt.0)then 
		res=-res
	endif 
	
	pm = dsqrt(res)
	
	end subroutine confadd
	

!==========================   ISP   ==================================
	double precision function ispadd(x,ni,ns,zi)
	
	use tailles
	
	implicit none
	
	integer::ni,ns
	double precision::val,mmspadd,x
	double precision,dimension(-2:npmax)::zi
	
	if(x.eq.zi(ni))then
		if(ni.le.ns-3)then
			val = 0.d0
		else
			if(ni.le.ns-2)then
				val = ((zi(ni)-zi(ni-1))*mmspadd(x,ni,ns,zi))*0.25d0
			else
				if (ni.eq.ns-1)then
					val = ((zi(ni)-zi(ni-2))*mmspadd(x,ni,ns,zi)+ &
					(zi(ni+3)-zi(ni-1))*mmspadd(x,ni,ns+1,zi))*0.25d0
				else
					if(ni.eq.ns)then
						val = ((zi(ni)-zi(ni-3))*mmspadd(x,ni,ns,zi)+ &
						(zi(ni+2)-zi(ni-2))*mmspadd(x,ni,ns+1,zi) &
						+(zi(ni+3)-zi(ni-1))*mmspadd(x,ni,ns+2,zi))*0.25d0
					else
						val = 1.d0
					endif
				endif
			endif   
		endif
	else   
		if(ni.lt.ns-3)then
			val = 0.d0
		else
			if(ni.eq.ns-3)then
				val = (x-zi(ni))*mmspadd(x,ni,ns,zi)*0.25d0
			else  
				if(ni.eq.ns-2)then
					val = ((x-zi(ni-1))*mmspadd(x,ni,ns,zi)+ &
					(zi(ni+4)-zi(ni))*mmspadd(x,ni,ns+1,zi))*0.25d0
				else   
					if (ni.eq.ns-1)then
						val =((x-zi(ni-2))*mmspadd(x,ni,ns,zi)+ &
						(zi(ni+3)-zi(ni-1))*mmspadd(x,ni,ns+1,zi) &
						+(zi(ni+4)-zi(ni))*mmspadd(x,ni,ns+2,zi))*0.25d0
					else
						if(ni.eq.ns)then
							val =((x-zi(ni-3))*mmspadd(x,ni,ns,zi)+ &
							(zi(ni+2)-zi(ni-2))*mmspadd(x,ni,ns+1,zi) &
							+(zi(ni+3)-zi(ni-1))*mmspadd(x,ni,ns+2,zi) &
							+(zi(ni+4)-zi(ni))*mmspadd(x,ni,ns+3,zi))*0.25d0
						else
							val = 1.d0
						endif
					endif
				endif
			endif
		endif 
	endif
	
	ispadd = val
	
	return
	
	end function ispadd
	
!==========================  MMSP   ==================================
	
	double precision function mmspadd(x,ni,ns,zi)
	
	use tailles
	
	implicit none
	
	integer::ni,ns
	double precision::val,x
	double precision,dimension(-2:npmax)::zi

	
	if(ni.lt.ns-3)then
		val = 0.d0
	else
		if(ns-3.eq.ni)then
			if(x.eq.zi(ni))then
				val = 0.d0
			else  
				val = (4.d0*(x-zi(ni))*(x-zi(ni)) &
				*(x-zi(ni)))/((zi(ni+4)-zi(ni))*(zi(ni+3) &
				-zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
			endif
		else 
			if(ns-2.eq.ni)then
				if(x.eq.zi(ni))then
					val = (4.d0*(zi(ni)-zi(ni-1))*(zi(ni)-zi(ni-1))) &
					/((zi(ni+3)-zi(ni-1))*(zi(ni+2)-zi(ni-1)) &
					*(zi(ni+1)-zi(ni-1)))
				else  
					val = (4.d0*(x-zi(ni-1))*(x-zi(ni-1)) &
					*(zi(ni+1)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
					-zi(ni-1))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)-zi(ni))) &
					+   (4.d0*(x-zi(ni-1))*(x-zi(ni)) &
					*(zi(ni+2)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
					-zi(ni))*(zi(ni+1)-zi(ni))*(zi(ni+2)-zi(ni-1)))  &
					+   (4.d0*(x-zi(ni))*(x-zi(ni)) &
					*(zi(ni+3)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+3) &
					-zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
				endif
			else   
				if (ns-1.eq.ni)then
					if(x.eq.zi(ni))then
						val = (4.d0*((zi(ni)-zi(ni-2))*(zi(ni+1) &
						-zi(ni)))/((zi(ni+2)-zi(ni-2))*(zi(ni+1) &
						-zi(ni-1))*(zi(ni+1)-zi(ni-2)))) &
						+((4.d0*((zi(ni)-zi(ni-1))*(zi(ni+2)-zi(ni))) &
						/((zi(ni+2)-zi(ni-2))*(zi(ni+2)-zi(ni-1)) &
						*(zi(ni+1)-zi(ni-1)))))
					else
						val = (4.d0*((x-zi(ni-2))*(zi(ni+1) &
						-x)*(zi(ni+1)-x))/((zi(ni+2) &
						-zi(ni-2))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)- &
						zi(ni))*(zi(ni+1)-zi(ni-2)))) &
						+((4.d0*((x-zi(ni-1))*(zi(ni+2)-x)  &
						*(zi(ni+1)-x))/((zi(ni+2)-zi(ni-2)) &
						*(zi(ni+2)-zi(ni-1))*(zi(ni+1)-zi(ni-1))* &
						(zi(ni+1)-zi(ni))))) &
						+((4.d0*((zi(ni+2)-x)*(zi(ni+2)-x)  &
						*(x-zi(ni)))/((zi(ni+2)-zi(ni-2)) &
						*(zi(ni+2)-zi(ni))*(zi(ni+2)-zi(ni-1))* &
						(zi(ni+1)-zi(ni)))))
					endif 
				else
					if(ni.eq.ns)then
						if(x.eq.zi(ni))then
							val =(4.d0*(x-zi(ni+1))*(x &
							-zi(ni+1))/((zi(ni+1)-zi(ni-1))*(zi(ni+1) &
							-zi(ni-2))*(zi(ni+1)-zi(ni-3))))
						else   
							val =(4.d0*(x-zi(ni+1))*(x &
							-zi(ni+1))*(zi(ni+1)-x)/((zi(ni+1) &
							-zi(ni-1))*(zi(ni+1)-zi(ni-2))*(zi(ni+1) &
							-zi(ni))*(zi(ni+1)-zi(ni-3))))
						endif
					else
						val = 0.d0
					endif
				endif
			endif
		endif
	endif
	
	mmspadd = val
	
	return
	
	end function mmspadd



	
!===================================================================
!========================          FUNCPA       ====================
	double precision function funcpaadd(b,np,id,thi,jd,thj,k0)
	
	use tailles
	use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m, &
	mm3,mm2,mm1,mm,im3,im2,im1,im,date,zi,pe,effet,nz1,nz2,stra, &
	t0,t1,c,nt0,nt1,nsujet,nva,ndate,nst,auxig,ng,ve,g,nig,indictronq,resnonpen
	use additiv,only: &
	dut1,dut2,ut1,ut2,correl,ngexact,ve2,sigma2,tau2,rho,cov,mid,betaaux,invD
	
	implicit none
	
	integer::n,np,id,jd,i,j,k,vj,cptg,ig,ip
	integer,dimension(ngmax)::cpt
	double precision::thi,thj,pe1,pe2,som1,som2, &
	res,vet,h1
	double precision,dimension(-2:npmax)::the1,the2
	double precision,dimension(np)::b,bh
	double precision,dimension(ngmax)::integrale1,funcaux
	double precision,dimension(2)::k0
	
	real,parameter::pi = 3.1415926535
        integer::restar,nf
!****** pour la maximisation avec marq98aux
	integer :: npaux,niaux,ieraux,istopaux
	DOUBLE PRECISION::resaux
	double precision,dimension((npmax*(npmax+3)/2))::vaux
	double precision,dimension(npmax)::baux
!****** derivanal
	double precision , dimension(ngmax)::res1,res2,res3,res4
	double precision , dimension(ngmax)::res5,res6
!      common /derivanal/res1,res2,res3,res4,res5,res6,res8

!****** u_tilde
	double precision  :: u_tilde,v_tilde
!      common /utilde/u_tilde,v_tilde
!******************


	j=0
	res=0.d0
	resnonpen=0.d0
	som2=0.d0
	pe=0.d0
	restar = 0
	nf = 1   
	
	do i=1,np
		bh(i)=b(i)
	end do 
	
	if (id.ne.0) bh(id)=bh(id)+thi 
	if (jd.ne.0) bh(jd)=bh(jd)+thj    
	
	if(nst.eq.2)then
		n = (nz1+2 +nz2+2)/nst
	else
		n = nz1+2
	endif
	
	do i=1,n
		the1(i-3)=(bh(i))*(bh(i))
		j = n+i 
		if (nst.eq.2) then
			the2(i-3)=(bh(j))*(bh(j))
		endif
	end do
	
	if(effet.eq.1) then
!terme de correlation entre 2 frailties/avec contraintes
		sigma2 = (bh(np-nva-1)*bh(np-nva-1)) ! variance intercept
		tau2 = (bh(np-nva)*bh(np-nva))  ! variance traitement * groupe
		cov=dsqrt(sigma2*tau2)*(2.d0 * dexp(bh(np-nva-2))/(1.d0+dexp(bh(np-nva-2)))-1.d0)
	endif
	
	vj = 0
	som1 = 0.d0
	dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))
	ut1(0) = 0.d0
	ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))

	if (nst.eq.2) then
		som2 = 0.d0
		dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))
		ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))
		ut2(0) = 0.d0
	endif
	
	do i=2,ndate-1
		do k = 2,n-2
			if (((date(i)).ge.(zi(k-1))).and.(date(i).lt.zi(k)))then
				j = k-1
				if ((j.gt.1).and.(j.gt.vj))then
					som1 = som1+the1(j-4)
					som2 = som2+the2(j-4)
					vj  = j
				endif   
			endif
		end do 
		
		ut1(i) = som1 +(the1(j-3)*im3(i))+(the1(j-2)*im2(i)) &
		+(the1(j-1)*im1(i))+(the1(j)*im(i))
		dut1(i) = (the1(j-3)*mm3(i))+(the1(j-2)*mm2(i)) &
		+(the1(j-1)*mm1(i))+(the1(j)*mm(i))

		if(nst.eq.2)then
			ut2(i) = som2 +(the2(j-3)*im3(i))+(the2(j-2)*im2(i)) &
			+(the2(j-1)*im1(i))+(the2(j)*im(i))
			dut2(i) = (the2(j-3)*mm3(i))+(the2(j-2)*mm2(i)) &
			+(the2(j-1)*mm1(i))+(the2(j)*mm(i))
		endif         
	end do
	
	
	i = n-2
	h1 = (zi(i)-zi(i-1))
	ut1(ndate)=som1+the1(i-4)+the1(i-3)+the1(i-2)+the1(i-1)
	dut1(ndate) = (4.d0*the1(i-1)/h1)   
!AD:the1(i-4) en th2(i-4)
	if(nst.eq.2)then
		ut2(ndate)=som2+the2(i-4)+the2(i-3)+the2(i-2)+the2(i-1)
		dut2(ndate) = (4.d0*the2(i-1)/h1) 
	endif
!AD:
!--------------------------------------------------------
!----------calcul de la vraisemblance ------------------
!---------------------------------------------------------

!---- avec ou sans variable explicative  ------cc
	
	do ig=1,ngexact!ng!
		res1(ig) = 0.d0
		res2(ig) = 0.d0
	end do
	cpt=0
!*******************************************     
!----- sans effet aleatoire dans le modele
!*******************************************     

	if (effet.eq.0) then
		
		do i=1,nsujet
	
			cpt(g(i))=cpt(g(i))+1
	
			if(nva.gt.0)then
			
				vet = 0.d0  
				 
				do j=1,nva
					vet =vet + bh(np-nva+j)*ve(i,j)
				end do
				
				vet = dexp(vet)
			else
				vet=1.d0
			endif
	
	
			if((c(i).eq.1).and.(stra(i).eq.1))then
				res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet)
			endif  

			if((c(i).eq.1).and.(stra(i).eq.2))then
				res2(g(i)) = res2(g(i))+dlog(dut2(nt1(i))*vet)
			endif

			if(stra(i).eq.1)then
				res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet-ut1(nt0(i))*vet 
			endif
	
			if(stra(i).eq.2)then
				res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet-ut2(nt0(i))*vet 
			endif
	
		end do       

		res = 0.d0

		cptg = 0
	
! k indice les groupes
		do k=1,ngexact 
			if(cpt(k).gt.0)then !nb de sujets dans un gpe=nig()               
				res = res-res1(k)+ res2(k) 
				cptg = cptg + 1 
			endif 
		end do
	
	else 
!*******************************************         
!-----avec deux effets aleatoires dans le modele et correl=0
!*********************************************
		if(effet.eq.1.and.correl.eq.0)then
		
			mid=0
			integrale1=0.d0
			
			do k=1,nsujet
				if(c(k).eq.1)then
					mid(g(k))=mid(g(k))+1
				endif
			end do 
	
			do k=1,nsujet
			
				if(nva.gt.0)then
					vet = 0.d0 
					do ip=1,nva
						vet =vet + bh(np-nva+ip)*ve(k,ip)
					end do
					vet = dexp(vet)
				else
					vet=1.d0
				endif
	
				if((c(k).eq.1).and.(stra(k).eq.1))then
					res2(g(k)) = res2(g(k))+dlog(dut1(nt1(k))*vet)
				endif  
				
				if((c(k).eq.1).and.(stra(k).eq.2))then
					res2(g(k)) = res2(g(k))+dlog(dut2(nt1(k))*vet)
				endif 
	
			end do 

!=========================================================================
!==== calcul des int�grales par transformation de LAPLACE pour chq gpe ===
!=========================================================================
	
			do ig=1,ng!ngexact
			
				res3(ig)=0.d0
				res4(ig)=0.d0
				res5(ig)=0.d0
				res6(ig)=0.d0
				auxig=ig
				baux(1)=0.05d0      !initialisation de u_tilde
				baux(2)=0.05d0      !initialisation de v_tilde

!======================================================================
! maximisation pour deux parametres, les 2 effets aleatoires : u et v 
!======================================================================
				npaux=2
				niaux=0
				ieraux=0
				istopaux=0
				vaux=0.d0!vecteur derivees 2nd et 1ERES
				funcaux=0.d0!vraisemblance
				resaux=0.d0!vraisemblance
				
				do ip=1,nva
					betaaux(ip)= bh(np-nva+ip)
				end do



	call marq98o(baux,npaux,niaux,vaux,resaux,ieraux,istopaux)

				u_tilde = baux(1)!u_tilde
				v_tilde = baux(2)!v_tilde

				if(effet.eq.1) then
				end if
				
				do k=1,nsujet
				
					if(nva.gt.0.and.g(k).eq.ig)then
						vet = 0.d0 
						do ip=1,nva
							vet =vet + bh(np-nva+ip)*ve(k,ip)
						end do
						vet = dexp(vet)
					else
						vet=1.d0
					endif
	
	
					if(g(k).eq.ig)then
					
						if(c(k).eq.1)then
							res3(ig) = res3(ig)+u_tilde+v_tilde*ve2(k,1)
						endif
	
						if(stra(k).eq.1)then
							res4(ig) = res4(ig) &
							+ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))
							res5(ig) = res5(ig) &
							+ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
							res6(ig) = res6(ig) &
							+ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2
						endif
						
						if(stra(k).eq.2)then
							res4(ig) = res4(ig) &
							+ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))
							res5(ig) = res5(ig) &
							+ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
							res6(ig) = res6(ig) &
							+ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2
						endif
	
					endif 
				end do 

!=====fin maximisation aux
	
				integrale1(ig)= -0.5d0*dlog(sigma2*tau2) &
				-0.5d0* &
				dlog(res4(ig)*res6(ig)+res4(ig)/tau2+res6(ig)/sigma2 &
				+1.d0/(sigma2*tau2)-res5(ig)*res5(ig)) & !det(K"(b)) 
				+res2(ig)+res3(ig)-res4(ig) &
				-0.5d0*((u_tilde**2)/sigma2+(v_tilde**2)/tau2)!-k(b)
			end do
	!======= fin calcul de log integrale
!======================================================================
	
	
			res = 0.d0
			do k=1,ng!ngexact  
				if(nig(k).gt.0)then
				
					if(indictronq.eq.0)then
						res = res+integrale1(k)!integrale1 donne le log de I directement
					else !troncature
				!		write(*,*)'***TRAITER TRONCATURE**'
						stop
					endif
				endif
			end do
		endif                !fin boucle effet= 1 and correl = 0

!=======================================================================

!*******************************************         
!-----avec deux effets aleatoires dans le modele et correl=1
!*********************************************
	
		if(effet.eq.1.and.correl.eq.1)then
	
			mid=0
			integrale1=0.d0
	
			do k=1,nsujet
				if(c(k).eq.1)then
					mid(g(k))=mid(g(k))+1
				endif
			end do 
	
			do k=1,nsujet
				if(nva.gt.0)then
					vet = 0.d0 
					do ip=1,nva
						vet =vet + bh(np-nva+ip)*ve(k,ip)
					end do
					vet = dexp(vet) 
				else
					vet=1.d0
				endif
	
				if((c(k).eq.1).and.(stra(k).eq.1))then
					res2(g(k)) = res2(g(k))+dlog(dut1(nt1(k))*vet)
				endif  
				if((c(k).eq.1).and.(stra(k).eq.2))then
					res2(g(k)) = res2(g(k))+dlog(dut2(nt1(k))*vet)
				endif 
	
			end do 

!=========================================================================
!==== calcul des int�grales par transformation de LAPLACE pour chq gpe ===
!=========================================================================
	
			do ig=1,ngexact!ng!
				res3(ig)=0.d0
				res4(ig)=0.d0
				res5(ig)=0.d0
				res6(ig)=0.d0
				auxig=ig
				baux(1)=0.05d0   !initialisation de u_tilde
				baux(2)=0.05d0   !initialisation de v_tilde

!======================================================================
! maximisation pour deux parametres, les 2 effets aleatoires : u et v 
!======================================================================
				npaux=2
				niaux=0
				vaux=0.d0!vecteur derivees 2nd et 1ERES
				funcaux=0.d0!vraisemblance
				resaux=0.d0!vraisemblance
				do ip=1,nva
					betaaux(ip)= bh(np-nva+ip)
				end do
	call marq98o(baux,npaux,niaux,vaux,resaux,ieraux,istopaux)
		
				u_tilde = baux(1)!u_tilde(ig)
				v_tilde = baux(2)!v_tilde
	
				if(effet.eq.1) then
				end if
				do k=1,nsujet
					if(nva.gt.0.and.g(k).eq.ig)then
						vet = 0.d0 
						do ip=1,nva
							vet =vet + bh(np-nva+ip)*ve(k,ip)
						end do
						vet = dexp(vet)
					else
						vet=1.d0
					endif
	
					if(g(k).eq.ig)then
						if(c(k).eq.1)then
							res3(ig) = res3(ig) &
							+u_tilde+v_tilde*ve2(k,1)
						endif
	
						if(stra(k).eq.1)then
							res4(ig) = res4(ig) &
							+ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))
							res5(ig) = res5(ig) &
							+ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
							res6(ig) = res6(ig) &
							+ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2
	
						endif
						if(stra(k).eq.2)then
							res4(ig) = res4(ig) &
							+ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))
							res5(ig) = res5(ig) &
							+ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
							res6(ig) = res6(ig) &
							+ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2
	
						endif
	
					endif 
				end do 
	
!=====fin maximisation aux
	
				cov=dsqrt(sigma2*tau2)* & !avec contrainte
				(2.d0 * dexp(bh(np-nva-2))/(1.d0+dexp(bh(np-nva-2)))-1.d0)
	
				integrale1(ig)= -0.5d0*dlog(sigma2*tau2-cov*cov) &!-0.5logdet(D)
				-0.5d0* &
				dlog((res4(ig)+tau2/(tau2*sigma2-cov**2)) &
				*(res6(ig)+sigma2/(tau2*sigma2-cov**2)) &
				-(res5(ig)-cov/(tau2*sigma2-cov**2))**2) &!-0.5*log(det(ka2)
				+res2(ig)+res3(ig)-res4(ig) &
				-0.5d0*((u_tilde**2)/sigma2+(v_tilde**2)/tau2 &
				-2.d0*u_tilde*v_tilde*cov/(sigma2*tau2)) &
				/(1.d0-(cov**2)/(sigma2*tau2))         !-ka
	
			end do
!======= fin calcul de log integrale
!======================================================================


!======= fin calcul de integrale
!         print*,'** integrale1 **',(integrale1(ig),ig,ig=1,ngexact)
!          stop
!======================================================================
	
	
			res = 0.d0
			do k=1,ngexact  !ng!
				if(nig(k).gt.0)then
					if(indictronq.eq.0)then
						res = res &
						+integrale1(k)!integrale1 donne le log de I directement
					else !troncature
						write(*,*)'***TRAITER TRONCATURE**'
						stop
					endif
				endif
			end do
	
		endif                     !fin boucle effet= 1 and correl = 1
	endif                     !fin boucle globale effet=0 

!----------calcul de la penalisation -------------------
	pe=0.d0
	pe1 = 0.d0
	pe2 = 0.d0
	do i=1,n-3
		pe1 = pe1+(the1(i-3)*the1(i-3)*m3m3(i))+(the1(i-2) &
		*the1(i-2)*m2m2(i))+(the1(i-1)*the1(i-1)*m1m1(i))+( &
		the1(i)*the1(i)*mmm(i))+(2.d0*the1(i-3)*the1(i-2)* &
		m3m2(i))+(2.d0*the1(i-3)*the1(i-1)*m3m1(i))+(2.d0* &
		the1(i-3)*the1(i)*m3m(i))+(2.d0*the1(i-2)*the1(i-1)* &
		m2m1(i))+(2.d0*the1(i-2)*the1(i)*m2m(i))+(2.d0*the1(i-1) &
		*the1(i)*m1m(i))
	
		if (nst.eq.2) then
			pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2) &
			*the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+( &
			the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)* &
			m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0* &
			the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)* &
			m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1) &
			*the2(i)*m1m(i))
		endif
	end do
	
	if (nst.eq.2) then
		pe = k0(1)*pe1 + k0(2)*pe2 
	else
		pe = k0(1)*pe1
	endif 
	resnonpen = res
	res = res - pe
	
	funcpaadd= res 
		
	return
	
	end function funcpaadd
	





!=====================cross validation

!========================          MNBRAK         ===================
	subroutine mnbrakadd(ax,bx,cx,fa,fb,fc,b,n)
	use tailles
	implicit none
	
	double precision::ax,bx,cx,fa,fb,fc,aux,res,dum,fu,q,r,u,ulim
	double precision,dimension(npmax)::b
	double precision,dimension(npmax,npmax)::y
	double precision::estimvadd,gold,glimit,tiny
	parameter (gold=1.618034d0,glimit=100.d0,tiny=1.d-20)
	integer::n,ni
	
	
	fa = estimvadd(ax,n,b,y,aux,ni,res)
	
	fb = estimvadd(bx,n,b,y,aux,ni,res)
	
	
	if(fb.gt.fa)then
		dum = ax
		ax = bx
		bx = dum
		dum = fb
		fb = fa
		fa = dum
	endif
	cx = bx + gold*(bx-ax)
	
	fc = estimvadd(cx,n,b,y,aux,ni,res)
	
1       if(fb.ge.fc)then
		r = (bx-ax)*(fb-fc)
		q = (bx-cx)*(fb-fa)
		u = bx-((bx-cx)*q-(bx-ax)*r)/ &
		(2.d0*sign(max(abs(q-r),tiny),q-r))
		ulim = bx + glimit*(cx-bx)
		if((bx-u)*(u-cx).gt.0.d0)then
			fu = estimvadd(u,n,b,y,aux,ni,res)
			if(fu.lt.fc)then
				ax = bx
				fa = fb
				bx = u
				fb = fu
				return
			else
				if(fu.gt.fb)then
					cx = u
					fc = fu
					return
				endif   
			endif
			u = cx + gold*(cx-bx)
			fu = estimvadd(u,n,b,y,aux,ni,res)
		else
			if((cx-u)*(u-ulim).gt.0.d0)then
				fu = estimvadd(u,n,b,y,aux,ni,res)
				if(fu.lt.fc)then
					bx = cx
					cx = u
					u = cx + gold*(cx-bx)
					fb = fc
					fc = fu
					fu = estimvadd(u,n,b,y,aux,ni,res)
				endif  
			else
				if((u-ulim)*(ulim-cx).ge.0.d0)then
					u = ulim
					fu = estimvadd(u,n,b,y,aux,ni,res)
				else
					u = cx + gold*(cx-bx)
					fu = estimvadd(u,n,b,y,aux,ni,res)
				endif
			endif   
		endif
		ax = bx
		bx = cx
		cx = u
		fa = fb
		fb = fc
		fc = fu
		goto 1
	endif
	
	return 
	
	end subroutine mnbrakadd

!========================      GOLDEN   =========================
	double precision function goldenadd(ax,bx,cx,tol,xmin,n,b,y,aux)
	
	use tailles
	
	implicit none
	
	double precision,dimension(npmax,npmax)::y
	double precision,dimension(npmax)::b
	double precision ax,bx,cx,tol,xmin,r,c,aux,res
	parameter (r=0.61803399d0,c=1.d0-r)
	double precision::f1,f2,x0,x1,x2,x3,estimvadd
	integer::n,ni
	
	x0 = ax
	x3 = cx
	if(abs(cx-bx).gt.abs(bx-ax))then
		x1 = bx
		x2 = bx + c*(cx-bx)
	else
		x2 = bx
		x1 = bx - c*(bx-ax)
	endif
	f1 = estimvadd(x1,n,b,y,aux,ni,res)
	f2 = estimvadd(x2,n,b,y,aux,ni,res)
	
1       if(abs(x3-x0).gt.tol*(abs(x1)+abs(x2)))then
			if(f2.lt.f1)then
				x0 = x1
				x1 = x2
				x2 = r*x1 + c*x3
				f1 = f2
				f2 = estimvadd(x2,n,b,y,aux,ni,res)
	
			else
				x3 = x2
				x2 = x1
				x1 = r*x2+c*x0
				f2 = f1
				f1 = estimvadd(x1,n,b,y,aux,ni,res)
	
			endif
			go to 1
		endif
		if(f1.lt.f2)then
			goldenadd = f1
			xmin = x1
		else
			goldenadd = f2
			xmin = x2
		endif
		return
	
	end function goldenadd
	

!========================          ESTIMV         ===================

	double precision function estimvadd(k00,n,b,y,aux,ni,res)
	
	use tailles
	use comon,only:t0,t1,c,nt0,nt1,nsujet,nva,ndate,nst,&
	date,zi,pe,effet,nz1,nz2,mm3,mm2,mm1,mm,im3,im2,im1,im
	use additiv,only:correl
	use optim
	
	implicit none
	
	double precision,dimension((n*(n+3)/2))::v
	double precision,dimension(n,n)::y
	double precision,dimension(ndatemax)::ut,dut
	double precision,dimension(n)::bh,b
	double precision,dimension(-2:npmax)::the	
	double precision::res,k00,som,h1,aux
	double precision,dimension(2)::k0
	integer n,ij,i,k,j,vj,ier,istop,ni
	double precision::ca,cb,dd,funcpaadd
	external::funcpaadd
	
	ca=0.d0
	cb=0.d0
	dd=0.d0
	j=0
	estimvadd=0.d0
	k0(1) = k00*k00
	k0(2)=0.d0
	
	call marq98J(k0,b,n,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaadd)	
	
	if(k0(1).gt.0.d0)then
		do ij=1,n
			the(ij-3)=(b(ij))*(b(ij))
			bh(ij) = (b(ij))*(b(ij))
		end do
	
		vj = 0
		som = 0.d0
		dut(1) = (the(-2)*4.d0/(zi(2)-zi(1)))
		ut(1) = the(-2)*dut(1)*0.25d0*(zi(1)-zi(-2))
		do i=2,ndate-1
			do k = 2,n-2
				if ((date(i).ge.zi(k-1)).and.(date(i).lt.zi(k)))then
					j = k-1
					if ((j.gt.1).and.(j.gt.vj))then
						som = som+the(j-4)
						vj  = j
					endif   
				endif
			end do 
			ut(i) = som +(the(j-3)*im3(i))+(the(j-2)*im2(i)) &
			+(the(j-1)*im1(i))+(the(j)*im(i))
			dut(i) = (the(j-3)*mm3(i))+(the(j-2)*mm2(i)) &
			+(the(j-1)*mm1(i))+(the(j)*mm(i))
		end do
		i = n-2
		h1 = (zi(i)-zi(i-1))
		ut(ndate) = som+ the(i-4) + the(i-3)+the(i-2)+the(i-1)
		dut(ndate) = (4.d0*the(i-1)/h1)
	
		call testadd(dut,k0,n,aux,y)
		estimvadd = - ((res-pe)) - aux
	
	else
		aux = -n
	endif

	return
	end function estimvadd
      
!=================calcul de la hessienne  et de omega  ==============
	subroutine testadd(dut,k0,n,res,y)
	
	use tailles
	use comon,only:date,zi,t0,t1,c,nt0,nt1,nsujet,nva,ndate,nst
	implicit none
	
	integer::n,i,j,np
	double precision::res,tra,d
	double precision,dimension(npmax,npmax)::hessh,hess,omeg,y
	double precision,dimension(ndatemax)::dut
	integer,dimension(npmax)::indx     
	double precision,dimension(2)::k0
	
	
	
	do i = 1,n
		do j = 1,n
		hess(i,j) = 0.d0 
		end do
	end do
	
	
	do i = 1,n
		do j = i,n
			call matadd(hess(i,j),dut,i,j,n)
		end do
	end do
	do i = 2,n
		do j = 1,i-1
			hess(i,j)=hess(j,i)
		end do
	end do
	
	
	call calcomegadd(n,omeg)
	
	do i = 1,n
		do j = 1,n
			hessh(i,j)=-hess(i,j)
			hess(i,j) = hess(i,j) - (2.d0*k0(1)*omeg(i,j)) 
		end do   
	end do
	
	np = n
	do i=1,n
		do j=1,n
			y(i,j)=0.d0
		end do
		y(i,i)=1.d0
	end do
	
	call ludcmpadd(hess,n,indx,d)
	
	do j=1,n
		call lubksbadd(hess,n,indx,y(1,j))
	end do
	
	tra = 0.d0
	do i=1,n
		do j=1,n
			tra = tra + y(i,j)*hessh(j,i)
		end do
	end do
	
	res = (tra)
	
	end subroutine testadd

!======================  LUBKSB  ======================================
	subroutine lubksbadd(a,n,indx,b)
	
	use tailles
	
	implicit none
	
	integer::n,i,ii,j,ll
	double precision::sum
	integer,dimension(npmax)::indx   
	double precision,dimension(npmax)::b
	double precision,dimension(npmax,npmax)::a
	
	ii = 0
	do i=1,n
		ll = indx(i)
		sum = b(ll)
		b(ll) = b(i)
		if(ii.ne.0)then
			do j=ii,i-1
				sum = sum -a(i,j)*b(j)
			end do
		else
			if(sum.ne.0.d0)then
				ii=i
			endif
		endif
		b(i)=sum
	end do
	do i=n,1,-1
		sum = b(i)
		do j = i+1,n
			sum = sum-a(i,j)*b(j)
		end do
		b(i)=sum/a(i,i)
	end do
	
	return
	
	end subroutine lubksbadd
	
	

!======================  LUDCMP  ======================================
	subroutine ludcmpadd(a,n,indx,d)
	
	use tailles
	
	implicit none
	
	integer::n,nmax,i,imax,j,k
	integer,dimension(npmax)::indx 
	double precision,dimension(npmax,npmax)::a         
	double precision::d,tiny,aamax,dum,sum
	parameter (nmax=500,tiny=1.d-20)
	double precision,dimension(nmax)::vv
	
	d = 1.d0
	imax=0
	do i=1,n
		aamax=0.d0
		do j=1,n
			if (dabs(a(i,j)).gt.aamax)then
				aamax=dabs(a(i,j))
			endif
		end do
		if (aamax.eq.0.d0) then
		end if
		vv(i) = 1.d0/aamax
	end do
	do j = 1,n
		do i=1,j-1
			sum = a(i,j)
			do k=1,i-1
				sum = sum - a(i,k)*a(k,j)
			end do
			a(i,j) = sum
		end do
		aamax = 0.d0
		do i = j,n
			sum = a(i,j)
			do k=1,j-1
				sum = sum -a(i,k)*a(k,j)
			end do
			a(i,j) = sum
			dum = vv(i)*dabs(sum)
			if (dum.ge.aamax) then
				imax = i
				aamax = dum
			endif
		end do
		if(j.ne.imax)then
			do k=1,n
				dum = a(imax,k)
				a(imax,k)=a(j,k)
				a(j,k) = dum
			end do
			d = -d
			vv(imax)=vv(j)
		endif
		indx(j)=imax
		if(a(j,j).eq.0.d0)then
			a(j,j)=tiny
		endif
		if(j.ne.n)then
			dum = 1.d0/a(j,j)
			do i = j+1,n
				a(i,j) = a(i,j)*dum
			end do
		endif
	end do
	
	return
		
	end subroutine ludcmpadd
!=======================  CALOMEG  ===========================
	subroutine calcomegadd(n,omeg)
	
	use tailles
	use comon,only:date,zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m
	
	implicit none
	
	integer::n,i,j
	double precision::calc00add,calc01add,calc02add
	double precision,dimension(npmax,npmax)::omeg
	
	
	do i=1,n
		do j=1,n
		omeg(i,j)=0.d0
		end do
	end do
	
	omeg(1,1)=calc00add(1,n)
	omeg(1,2)=calc01add(1,n)
	omeg(1,3)=calc02add(1,n)
	omeg(1,4)=m3m(1)
	omeg(2,1)=omeg(1,2)
	omeg(2,2)=calc00add(2,n)
	omeg(2,3)=calc01add(2,n)
	omeg(2,4)=calc02add(2,n)
	omeg(2,5)=m3m(2)
	omeg(3,1)=omeg(1,3)
	omeg(3,2)=omeg(2,3)
	omeg(3,3)=calc00add(3,n)
	omeg(3,4)=calc01add(3,n)
	omeg(3,5)=calc02add(3,n)
	omeg(3,6)=m3m(3)
	do i=4,n-3
		omeg(i,i-3)=omeg(i-3,i)
		omeg(i,i-2)=omeg(i-2,i)
		omeg(i,i-1)=omeg(i-1,i)
		omeg(i,i)=calc00add(i,n)
		omeg(i,i+1)=calc01add(i,n)
		omeg(i,i+2)=calc02add(i,n)
		omeg(i,i+3)=m3m(i)
	end do   
	omeg(n-2,n-5)=omeg(n-5,n-2)
	omeg(n-2,n-4)=omeg(n-4,n-2)
	omeg(n-2,n-3)=omeg(n-3,n-2)
	omeg(n-2,n-2)=calc00add(n-2,n)
	omeg(n-2,n-1)=calc01add(n-2,n)
	omeg(n-2,n)=calc02add(n-2,n)
	omeg(n-1,n-4)=omeg(n-4,n-1)
	omeg(n-1,n-3)=omeg(n-3,n-1)
	omeg(n-1,n-2)=omeg(n-2,n-1)
	omeg(n-1,n-1)=calc00add(n-1,n)
	omeg(n-1,n)=calc01add(n-1,n)
	omeg(n,n-3)=omeg(n-3,n)
	omeg(n,n-2)=omeg(n-2,n)
	omeg(n,n-1)=omeg(n-1,n)
	omeg(n,n)=calc00add(n,n)
	
	end subroutine calcomegadd


!====================  MAT  ==================================
	subroutine matadd(res,dut,k,l,n)
	
	use tailles
	use comon,only:date,zi,t0,t1,c,nt0,nt1,nsujet,nva,ndate,nst

	implicit none
	
	integer::k,l,j,ni,n,i
	double precision,dimension(ndatemax)::dut
	double precision::res,res1,mspadd,aux2,u2
	
!---------- calcul de la hessienne ij ------------------
	res = 0.d0
	res1 = 0.d0
	do i=1,nsujet
		if(c(i).eq.1)then  !event
			u2 = dut(nt1(i)) 
			do j = 2,n-2
			if((date(nt1(i)).ge.zi(j-1)).and.(date(nt1(i)).lt.zi(j)))then
				ni = j-1
			endif
		end do 
			if(date(nt1(i)).eq.zi(n-2))then
				ni = n-2
			endif   
!-------attention numero spline 
			aux2 = mspadd(nt1(i),ni,k)*mspadd(nt1(i),ni,l)
			if (u2.le.0.d0)then
				res1 = 0.d0
			else   
				res1 = - aux2/(u2*u2)
			endif  
		else !censure  
			res1 = 0.d0
		endif 
		res = res + res1
	end do   
	
	end subroutine matadd

!==========================  MSP   ==================================
	double precision function mspadd(i,ni,ns)
	
	use tailles
	use comon,only:date,zi
	
	implicit none
	
	integer::ni,ns,i
	double precision::val
	
	
	if(ni.lt.ns-3)then
		val = 0.d0
	else
		if(ns-3.eq.ni)then
			if(date(i).eq.zi(ni))then
				val = 0.d0
			else  
				val = (4.d0*(date(i)-zi(ni))*(date(i)-zi(ni)) &
				*(date(i)-zi(ni)))/((zi(ni+4)-zi(ni))*(zi(ni+3) &
				-zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
			endif
		else 
			if(ns-2.eq.ni)then
				if(date(i).eq.zi(ni))then
					val = (4.d0*(zi(ni)-zi(ni-1))*(zi(ni)-zi(ni-1))) &
					/((zi(ni+3)-zi(ni-1))*(zi(ni+2)-zi(ni-1)) &
					*(zi(ni+1)-zi(ni-1)))
				else  
					val = (4.d0*(date(i)-zi(ni-1))*(date(i)-zi(ni-1)) &
					*(zi(ni+1)-date(i)))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
					-zi(ni-1))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)-zi(ni))) &
					+   (4.d0*(date(i)-zi(ni-1))*(date(i)-zi(ni)) &
					*(zi(ni+2)-date(i)))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
					-zi(ni))*(zi(ni+1)-zi(ni))*(zi(ni+2)-zi(ni-1)))  &
					+   (4.d0*(date(i)-zi(ni))*(date(i)-zi(ni)) &
					*(zi(ni+3)-date(i)))/((zi(ni+3)-zi(ni-1))*(zi(ni+3) &
					-zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
				endif
			else   
				if (ns-1.eq.ni)then
					if(date(i).eq.zi(ni))then
						val = (4.d0*((zi(ni)-zi(ni-2))*(zi(ni+1) &
						-zi(ni)))/((zi(ni+2)-zi(ni-2))*(zi(ni+1) &
						-zi(ni-1))*(zi(ni+1)-zi(ni-2)))) &
						+((4.d0*((zi(ni)-zi(ni-1))*(zi(ni+2)-zi(ni)))  &
						/((zi(ni+2)-zi(ni-2))*(zi(ni+2)-zi(ni-1)) &
						*(zi(ni+1)-zi(ni-1)))))
					else
						val = (4.d0*((date(i)-zi(ni-2))*(zi(ni+1) &
						-date(i))*(zi(ni+1)-date(i)))/((zi(ni+2) &
						-zi(ni-2))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)- &
						zi(ni))*(zi(ni+1)-zi(ni-2)))) &
						+((4.d0*((date(i)-zi(ni-1))*(zi(ni+2)-date(i))  &
						*(zi(ni+1)-date(i)))/((zi(ni+2)-zi(ni-2)) &
						*(zi(ni+2)-zi(ni-1))*(zi(ni+1)-zi(ni-1))* &
						(zi(ni+1)-zi(ni))))) &
						+((4.d0*((zi(ni+2)-date(i))*(zi(ni+2)-date(i))  &
						*(date(i)-zi(ni)))/((zi(ni+2)-zi(ni-2)) &
						*(zi(ni+2)-zi(ni))*(zi(ni+2)-zi(ni-1))* &
						(zi(ni+1)-zi(ni)))))
					endif 
				else
					if(ni.eq.ns)then
						if(date(i).eq.zi(ni))then
							val =(4.d0*(date(i)-zi(ni+1))*(date(i) &
							-zi(ni+1))/((zi(ni+1)-zi(ni-1))*(zi(ni+1) &
							-zi(ni-2))*(zi(ni+1)-zi(ni-3))))
						else   
							val =(4.d0*(date(i)-zi(ni+1))*(date(i) &
							-zi(ni+1))*(zi(ni+1)-date(i))/((zi(ni+1) &
							-zi(ni-1))*(zi(ni+1)-zi(ni-2))*(zi(ni+1) &
							-zi(ni))*(zi(ni+1)-zi(ni-3))))
						endif
					else
						val = 0.d0
					endif
				endif
			endif
		endif
	endif
	
	mspadd = val
		
	return
		
	end function mspadd

!==========================   SP   ==================================
	double precision function spadd(i,ni,ns)
	
	use tailles
	use comon,only:date,zi
	
	implicit none
	
	integer::ni,ns,i
	double precision::val,mspadd
	
	if(date(i).eq.zi(ni))then
		if(ni.le.ns-3)then
			val = 0.d0
		else
			if(ni.le.ns-2)then
				val = ((zi(ni)-zi(ni-1))*mspadd(i,ni,ns))*0.25d0
			else
				if (ni.eq.ns-1)then
					val = ((zi(ni)-zi(ni-2))*mspadd(i,ni,ns)+ &
					(zi(ni+3)-zi(ni-1))*mspadd(i,ni,ns+1))*0.25d0
				else
					if(ni.eq.ns)then
						val = ((zi(ni)-zi(ni-3))*mspadd(i,ni,ns)+ &
						(zi(ni+2)-zi(ni-2))*mspadd(i,ni,ns+1) &
						+(zi(ni+3)-zi(ni-1))*mspadd(i,ni,ns+2))*0.25d0
					else
						val = 1.d0
					endif
				endif
			endif   
		endif
	else   
		if(ni.lt.ns-3)then
			val = 0.d0
		else
			if(ni.eq.ns-3)then
				val = (date(i)-zi(ni))*mspadd(i,ni,ns)*0.25d0
			else  
				if(ni.eq.ns-2)then
					val = ((date(i)-zi(ni-1))*mspadd(i,ni,ns)+ &
					(zi(ni+4)-zi(ni))*mspadd(i,ni,ns+1))*0.25d0
				else   
					if (ni.eq.ns-1)then
						val =((date(i)-zi(ni-2))*mspadd(i,ni,ns)+ &
						(zi(ni+3)-zi(ni-1))*mspadd(i,ni,ns+1) &
						+(zi(ni+4)-zi(ni))*mspadd(i,ni,ns+2))*0.25d0
					else
						if(ni.eq.ns)then
							val =((date(i)-zi(ni-3))*mspadd(i,ni,ns)+ &
							(zi(ni+2)-zi(ni-2))*mspadd(i,ni,ns+1) &
							+(zi(ni+3)-zi(ni-1))*mspadd(i,ni,ns+2) &
							+(zi(ni+4)-zi(ni))*mspadd(i,ni,ns+3))*0.25d0
						else
							val = 1.d0
						endif
					endif
				endif
			endif
		endif 
	endif
	spadd = val
		
	return
		
	end function spadd
!================

!=========================  CALC00  =========================
	double precision function calc00add(j,n) 
	
	use tailles
	use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m
		
	implicit none
	
	double precision::part
	integer::j,n
	
	
	if(j.eq.1)then
		part = m3m3(j)
	else
		if(j.eq.2)then
			part = m3m3(j) + m2m2(j-1)
		else
			if(j.eq.3)then
				part = m3m3(j) + m2m2(j-1) + m1m1(j-2)
			else
				if(j.eq.n-2)then
					part = m2m2(j-1) + m1m1(j-2) + mmm(j-3)
				else   
					if(j.eq.n-1)then
						part = mmm(j-3) + m1m1(j-2)
					else
						if(j.eq.n)then
							part = mmm(j-3)
						else   
							part=mmm(j-3)+m1m1(j-2)+m2m2(j-1)+m3m3(j)
						endif
					endif
				endif   
			endif   
		endif   
	endif 
	
	calc00add = part
	
	return
	
	end function calc00add
	
!=========================  CALC01  =========================
	
	double precision function calc01add(j,n)
	
	use tailles
	use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m
		
	implicit none
	
	double precision::part
	integer::j,n
	
	
	
	if(j.eq.1)then
		part = m3m2(j)
	else   
		if(j.eq.2)then
			part = m3m2(j) + m2m1(j-1) 
		else
			if(j.eq.n-2)then
				part = m1m(j-2) + m2m1(j-1) 
			else
				if(j.ne.n-1)then
					part = m3m2(j) + m2m1(j-1) + m1m(j-2)
				else
					part = m1m(j-2)
				endif
			endif   
		endif
	endif   
	
	calc01add = part
	
	return
	
	end function calc01add

!=========================  CALC02  =========================
	double precision function calc02add(j,n)
	
	use tailles
	use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m
		
	implicit none
	
	
	
	double precision::part
	integer::j,n
	
	
	
	if(j.eq.1)then
		part = m3m1(j)
	else   
		if(j.ne.n-2)then
			part = m3m1(j) + m2m(j-1) 
		else
			part = m2m(j-1)
		endif
	endif   
	
	calc02add = part
	
	return
	
	end function calc02add


!===============================    MARQ98AUX =========================
!

 	
	subroutine marq98o(b,m,ni,v,rl,ier,istop)
	
	use tailles
	use optim
	use parameters,only:maxiter
	use comon,only:t0,t1,c,nt0,nt1,nsujet,nva,ndate,nst,auxig
	
	Implicit none
	
	integer::m,ni,nql,i,ii,nfmax,idpos,ier,istop,igrad,ncount,id,jd,i0
	double precision::da,dm,ga,tr,ca,cb,epsa,epsb,rl &
	,funcpao,step,eps,epsd,vw,fi,z,rl1,th,ep,dd,auxmax     
	double precision,dimension(m*(m+3)/2)::V,fu
	double precision,dimension(m)::b,delta,bh,b1
	double precision,dimension(2)::zero
	double precision::maxtadd
	
	zero(1)=0.d0
	zero(2)=0.d0
	id=0
	jd=0
	z=0.d0
	th=1.d-5
	eps=1.d-7
	epsa=1.d-6
	epsb=1.d-6
	epsd=1.d-6
	
	nfmax=m*(m+1)/2
	ca=epsa+1.d0
	cb=epsb+1.d0
	rl1=-1.d+10
	ni=0
	istop=0
	da=0.01d0
	dm=5.d0
	
	nql=1
	
10   continue

		
	z=0.d0
	i0=0 	
	rl=funcpao(b,m,i0,z,i0,z)
!	write(*,*)'iter',ni,'vrais',rl1
	rl1=rl 
	
	call derivao(b,m,v,rl)
	
	dd = 0.d0
	do i=m*(m+1)/2+1,m*(m+3)/2
		dd = dd + v(i)*v(i)
	end do
	
	dd=dd/dabs(RL)

	if (ca.lt.epsa.and.cb.lt.epsb.and.dd.lt.epsd) goto 100
	
		tr=0.d0
		
		do i=1,m
			ii=i*(i+1)/2
			tr=tr+dabs(v(ii))
		end do 
		
		tr=tr/dble(m)
		ncount=0
		ga=0.01d0
		
400      	do i=1,nfmax+m
			fu(i)=v(i)
		end do
		
		do i=1,m
			ii=i*(i+1)/2
			if (v(ii).ne.0) then
				fu(ii)=v(ii)+da*((1.d0-ga)*dabs(v(ii))+ga*tr)
			else
				fu(ii)=da*ga*tr
			end if
		end do
		
		call dcholej(fu,m,nql,idpos)
	
		if (idpos.ne.0) then  
	
			if(b(1).lt.-1.d0.or.b(1).gt.1.d0.or.b(2).lt.-1.d0.or.b(2).gt.1.d0) then
				b(1)=0.5d0
				b(2)=0.5d0
				goto 110
			endif
				
			ncount=ncount+1
			
			if (ncount.le.3.or.ga.ge.1.d0) then
				da=da*dm
			else
				ga=ga*dm
				if (ga.gt.1.d0) ga=1.d0
			end if
			
			if (ncount > 10) then
				fu=0.d0
				do i=1,m
					ii=i*(i+1)/2
					fu(ii)=1
				end do
				return
			end if
			
			goto 400
			
		else
			do i=1,m
				delta(i)=fu(nfmax+i)
				b1(i)=b(i)+delta(i)
			end do
	
			rl=funcpao(b1,m,id,z,jd,z)
			
			igrad=1
			
			if (rl1.lt.rl) then
				if(da.lt.eps) then
					da=eps
				else
					da=da/(dm+2.d0)
				end if
				goto 800
			end if
		end if
		
		call dmaxt(maxtadd,delta,m)
		auxmax=maxtadd
		
		if(auxmax.eq.0.d0) then
			vw=1.D-5
		else
		!	call dmaxt(maxtadd,delta,m)
			vw=th/auxmax
		end if

		step=dlog(1.5d0)
		
		call searpaso(vw,step,b,bh,m,delta,fi) 
		
		rl=-fi
		
		IF(rl1.gt.rl) then
		end if
		
		do i=1,m
			delta(i)=vw*delta(i)
		end do
		
		da=(dm-3.d0)*da  
	
800       cb=dabs(rl1-rl)

		ca=0.d0
		
		do i=1,m
			ca=ca+delta(i)*delta(i)
		end do

		do i=1,m
			b(i)=b(i)+delta(i)
		end do

		ni=ni+1
		
		if (ni.gt.maxiter) then
			istop=2	
			goto 110
		end if
		goto 10
	
100     continue

	istop=1
	
	ep=10.d-10
	
	call dsinvj(v,m,ep,ier)	
	if (ier.eq.-1) then
	write(6,103)          
103       format(1x,'echec inversion mat information ds marq98aux')	
		istop=3
	end if	


110       continue

	return
	
	end subroutine marq98o
	     
!=================================DERIVAAUX =======================
	
	subroutine derivao(b,m,v,rl)
		
	integer ::i0,iun,m,m1,ll,i,k,j
	double precision :: funcpao,thn,th,z,rl,vl
	double precision :: th2
	double precision ,dimension(m):: fcith
	double precision ,dimension(m):: b
	double precision , dimension(m*(m+3)/2)::V
	
	
	th=1.d-5
	thn=-th
	th2=th*th
	z=0.d0
	i0=0
	iun =1
	
	rl=funcpao(b,m,iun,z,iun,z)
	
	do i=1,m
		fcith(i)=funcpao(b,m,i,th,i0,z)
	end do
	
	k=0
	m1=m*(m+1)/2
	ll=m1
	
	do i=1,m
		ll=ll+1
		vl=(fcith(i)-funcpao(b,m,i,thn,i0,z))/(2.d0*th)
		v(ll)=vl
		do  j=1,i
			k=k+1
			v(k)=-(funcpao(b,m,i,th,j,th)-fcith(j)-fcith(i)+rl)/th2 
		end do
	end do
	
	return
		
	end subroutine derivao




!================================  SEARPAS joly    ==============================

	subroutine searpaso(vw,step,b,bh,m,delta,fim)
	
	implicit none
	
	integer::m,i      
	double precision,dimension(m)::b
	double precision::vw
	double precision,dimension(m)::bh,delta    
	double precision::fim,step   
	double precision::vlw,vlw1,vlw2,vlw3,vm,fi1,fi2,fi3    


	
	vlw1=dlog(vw)
	vlw2=vlw1+step
	call valfpao(vlw1,fi1,b,bh,m,delta)
	call valfpao(vlw2,fi2,b,bh,m,delta)       
	
	if(fi2.ge.fi1) then
		vlw3=vlw2
		vlw2=vlw1
		fi3=fi2
		fi2=fi1
		step=-step
	
		vlw1=vlw2+step
		call valfpao(vlw1,fi1,b,bh,m,delta)   
		if(fi1.gt.fi2) goto 50
	else 
		vlw=vlw1
		vlw1=vlw2
		vlw2=vlw
		fim=fi1
		fi1=fi2
		fi2=fim
	end if
	
	do i=1,40
		vlw3=vlw2
		vlw2=vlw1
		fi3=fi2
		fi2=fi1
	
		vlw1=vlw2+step
		call valfpao(vlw1,fi1,b,bh,m,delta)
		if(fi1.gt.fi2) goto 50
		if(fi1.eq.fi2) then
		fim=fi2
		vm=vlw2 
		goto 100
		end if
	end do
	
50     continue
	
	vm=vlw2-step*(fi1-fi3)/(2.d0*(fi1-2.d0*fi2+fi3))   
	call valfpao(vm,fim,b,bh,m,delta)	
	if(fim.le.fi2) goto 100
	vm=vlw2
	fim=fi2
100   continue
	vw=dexp(vm)
	
	return
	
	end subroutine searpaso

   
!========================   VALFPAAUX   ============================== 
       
	subroutine valfpao(vw,fi,b,bk,m,delta)
	
	implicit none
	
	integer::m  
	double precision,dimension(m)::b,delta  
	double precision,dimension(m)::bk 
	double precision::fi 
	double precision::vw,funcpao,z	
	integer::i0,i
	
	z=0.d0
	i0=1
	do i=1,m
	bk(i)=b(i)+dexp(vw)*delta(i)
	end do
	fi=-funcpao(bk,m,i0,z,i0,z)
	
	return
		
	end subroutine valfpao

!========================    DEBUT FUNCPAAUX       ====================
	double precision function funcpao(b,np,id,thi,jd,thj)
	
	use tailles
	use comon,only:g,nig,t0,t1,c,nt0,nt1,nsujet,nva,ndate,nst, &
	stra,ve,auxig
	use additiv,only:dut1,dut2,ut1,ut2,betaaux,ve2,sigma2,tau2,rho,cov
	implicit none
	
	integer::np,id,jd,i,k,ip
	double precision,dimension(np)::bhaux,b
	double precision::thi,thj,res,vet
!****** u_tilde
	double precision  :: u_tilde,v_tilde
!****** derivanal
	double precision , dimension(ngmax)::res3,res4
	double precision , dimension(ngmax)::res5,res6,res8

!==============================================
!================POUR UN GROUPE AUXIG donn� !
!==============================================
	
	res3=0.d0
	res4=0.d0
	res5=0.d0
	res6=0.d0
	res8=0.d0
	
	do i=1,np
		bhaux(i)=b(i)
	end do 
	
	
	if (id.ne.0) bhaux(id)=bhaux(id)+thi 
	if (jd.ne.0) bhaux(jd)=bhaux(jd)+thj    
		
	u_tilde=bhaux(1)
	v_tilde=bhaux(2)
!--------------------------------------------------------
!----------calcul de la vraisemblance ------------------
!---------------------------------------------------------
	
	do k=1,nsujet
		if(nva.gt.0.and.g(k).eq.auxig)then
			vet = 0.d0 
			do ip=1,nva
				vet =vet + betaaux(ip)*ve(k,ip)
			end do
			vet = dexp(vet)
		else
			vet=1.d0
		endif
		
		if(g(k).eq.auxig)then
		
			if(c(k).eq.1)then
				res3(auxig) = res3(auxig)+u_tilde+v_tilde*ve(k,1)
				res8(auxig) = res8(auxig)+ve(k,1)
			endif
			if(stra(k).eq.1)then
				res4(auxig) =res4(auxig)+ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve(k,1))
				
				res5(auxig) = res5(auxig)+ut1(nt1(k))*vet* &
				dexp(u_tilde+v_tilde*ve(k,1))*ve(k,1)
				
				res6(auxig) = res6(auxig)+ut1(nt1(k))*vet* &
				dexp(u_tilde+v_tilde*ve(k,1))*(ve(k,1))**2
			endif
			if(stra(k).eq.2)then
				res4(auxig) = res4(auxig) &
				+ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve(k,1))
				
				res5(auxig) = res5(auxig)+ut2(nt1(k))*vet* &
				dexp(u_tilde+v_tilde*ve(k,1))*ve(k,1)
				
				res6(auxig) = res6(auxig)+ut2(nt1(k))*vet* &
				dexp(u_tilde+v_tilde*ve(k,1))*(ve(k,1))**2
			endif
		
		endif 
	end do 
		
	res=-res3(auxig) + res4(auxig)+0.5d0*(((u_tilde)**2)/sigma2+((v_tilde)**2)/tau2 &
	-2.d0*u_tilde*v_tilde*cov/(sigma2*tau2))/(1.d0-(cov**2)/(sigma2*tau2))         !-ka
	
	funcpao= -res
	
	return
	
	end function funcpao


	logical function isnann(x)

	implicit none
	
	double precision,intent(in)::x
	
	if (x .ne. x) then
		isnann=.true.
	else
		isnann=.false.
	end if

	end function isnann
