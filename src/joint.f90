


	module donnees
	
	implicit none
	double precision,dimension(6)::cof
	double precision::stp,half,one,fpf
	double precision,dimension(20),save::x,w!The abscissas-weights.
	
	DATA w/0.181080062419,0.422556767879,0.666909546702,0.91535237279, &
	1.1695397071,1.43135498624,1.7029811359,1.98701589585, &
	2.28663576323,2.60583465152,2.94978381794,3.32539569477, &
	3.74225636246,4.21424053477,4.76252016007,5.42172779036, &
	6.25401146407,7.38731523837,9.15132879607,12.8933886244/
	
	DATA x/0.070539889692,0.372126818002,0.916582102483,1.70730653103, &
	2.74919925531,4.04892531384,5.61517497087,7.45901745389, &
	9.59439286749,12.0388025566,14.8142934155,17.9488955686, &
	21.4787881904,25.4517028094,29.9325546634,35.0134341868, &
	40.8330570974,47.6199940299,55.8107957541,66.5244165252/
	
	data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0, &
	-1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
	data half,one,fpf/0.5d0,1.0d0,5.5d0/
	
	end module donnees
      

!=======================================================================================
!=======================================================================================
!                                FIN marq98 version optim
!=======================================================================================
!=======================================================================================


	module splines
	
	double precision,dimension(:),allocatable,save::aux
	double precision,dimension(:),allocatable,save:: v
	double precision,dimension(:,:),allocatable,save::I1_hess,H1_hess
	double precision,dimension(:,:),allocatable,save::I2_hess,H2_hess
	double precision,dimension(:,:),allocatable,save::HI2
	double precision,dimension(:,:),allocatable,save::HIH,IH,HI
	double precision,dimension(:,:),allocatable,save::BIAIS
	double precision,dimension(:),allocatable,save:: vax,vaxdc
	integer,dimension(:),allocatable,save::filtre,filtre2
	character(LEN=20),dimension(:),allocatable,save:: nomvar,nomvar2
	integer,save::ver
	
	end module splines





!--entête pour fortran	
	subroutine joint(nsujet0,ng0,nz0,k0,tt00,tt10,ic0,groupe0      &
	,tt0dc0,tt1dc0,icdc0,nva10,vax0,nva20,vaxdc0,noVar1,noVar2,maxit0   &
	,np,b,H_hessOut,HIHOut,resOut,traceLCV,x1Out,lamOut,suOut,x2Out &
	,lam2Out,su2Out,ni,cpt,cpt_dc,ier)
!AD: add for new marq    
	use parameters	
!AD:end
	use comon
	use additiv,only:aux1,aux2
	use tailles
	use splines
	use optim
	
	IMPLICIT NONE  
	
	integer::maxit0
	integer,intent(in)::nsujet0,ng0,nz0,nva10,nva20
	integer::np
	integer,dimension(nsujet0)::groupe0,ic0
	integer,dimension(ng0),intent(in)::icdc0
	double precision,dimension(ng0)::tt0dc0,tt1dc0
	double precision,dimension(nsujet0)::tt00,tt10
	double precision,dimension(2)::k0
	double precision,dimension(nsujet0,nva10),intent(in):: vax0
	double precision,dimension(ng0,nva20),intent(in):: vaxdc0
	double precision,dimension(np,np)::H_hessOut,HIHOut
	double precision::resOut
	double precision,dimension(99)::x1Out,x2Out
	double precision,dimension(99,3)::lamOut,suOut,lam2Out,su2Out   
	integer::ss,sss
	double precision,dimension(np):: b
	double precision,intent(out)::traceLCV
	
	integer,intent(in)::noVar1,noVar2
	integer,intent(out)::cpt,cpt_dc,ier,ni
	integer::groupe,ij,kk,j,k,nz,n,cptcens,ii,iii,iii2,cptstr1,cptstr2   &
	,i,ic,icdc,istop,cptni,cptni1,cptni2,nb_echec,nb_echecor,id,cptbiais &
	,cptaux,cptauxdc,nb0recu      
	double precision::tt0,tt0dc,moyrecu,tt1,tt1dc,h,res,min,mindc,max, &
	maxdc,maxt,maxtdc,moy_peh0,moy_peh1,lrs,BIAIS_moy
	double precision,dimension(2)::res01
!AD: add for new marq
	double precision::ca,cb,dd,funcpaj
	external::funcpaj		
	indic_alpha=1 
	maxiter=maxit0
!AD:add for new marq
	epsa=1.d-3
	epsb=1.d-3
	epsd=1.d-3
!AD:end 
!----- Type of model joint==1
	model=1
				
	lrs=0.d0
	moy_peh0=0.d0
	moy_peh1=0.d0
	
	nb_echec=0
	nb_echecor=0
	nb0recu =0
	moyrecu =0.d0             
		
	ngmax=ng0
	ng=ng0
	allocate(nig(ngmax),cdc(ngmax),nt0dc(ngmax),nt1dc(ngmax),t0dc(ngmax), &
	t1dc(ngmax),aux1(ngmax),aux2(ngmax),res1(ngmax),res3(ngmax),mi(ngmax))
	
	nsujetmax=nsujet0
	nsujet=nsujet0	    
	allocate(t0(nsujetmax),t1(nsujetmax),c(nsujetmax),nt0(nsujetmax),nt1(nsujetmax) &
	,stra(nsujetmax),g(nsujetmax),aux(2*nsujetmax))
	
	ndatemaxdc=2*ng0     
	allocate(mm3dc(ndatemaxdc),mm2dc(ndatemaxdc),mm1dc(ndatemaxdc),mmdc(ndatemaxdc) &
	,im3dc(ndatemaxdc),im2dc(ndatemaxdc),im1dc(ndatemaxdc),imdc(ndatemaxdc))


	nst=2
	
	ni=0      
!---debut des iterations de simulations       
	id=1
	cptni=0
	cptni1=0
	cptni2=0
	biais_moy=0.d0
	cptbiais=0
	cptaux=0
	cptauxdc=0
	ij=0
	kk=0
!**************************************************
!********************* prog spline****************
	effet=1
	res01(1)=0.d0
	res01(2)=0.d0
!------------  entre non fichier et nombre sujet -----        
	nvarmax=ver
	
	allocate(vax(nva10),vaxdc(nva20))
			
	nva1=nva10
	nva2=nva20
	nva = nva1+nva2
	nvarmax=nva
	
	allocate(ve(nsujetmax,nvarmax),vedc(ngmax,nvarmax))
	allocate(filtre(nva10),filtre2(nva20))
	do i=1,ng
		nig(i) = 0
	end do
      
! AD: recurrent
	if (noVar1.eq.1) then 
!		write(*,*)'filtre 1 desactive'
		filtre=0
		nva1=0
	else
		do i=1,nva10
			filtre(i)=1
		end do	
	end if	
!AD:death
	if (noVar2.eq.1) then 
!		write(*,*)'filtre 2 desactive'
		filtre2=0
		nva2=0
	else
		do i=1,nva20
			filtre2(i)=1
		enddo 
	end if	
	
	if ((noVar1.eq.1).or.(noVar2.eq.1)) then
		nva = nva1+nva2 
	end if
	
		
!AD:end
                  
!------------  lecture fichier -----------------------
	maxt = 0.d0
	maxtdc = 0.d0
	cpt = 0
	cptcens = 0
	cpt_dc = 0
	k = 0
	cptstr1 = 0
	cptstr2 = 0
   
!ccccccccccccccccccccc
! pour le deces
!cccccccccccccccccccc

	do k = 1,ng 
		tt0dc=tt0dc0(k)
		tt1dc=tt1dc0(k)
		icdc=icdc0(k)
		groupe=groupe0(k)
		do j=1,nva20	    
			vaxdc(j)=vaxdc0(k,j)
		enddo	    	    
		if(tt0dc.gt.0.d0)then
			cptauxdc=cptauxdc+1
		endif                  
!------------------   deces c=1 pour donn�es de survie
		if(icdc.eq.1)then
			cpt_dc = cpt_dc + 1
			cdc(k)=1
			t0dc(k) = tt0dc      !/100.d0
			t1dc(k) = tt1dc      !+0.000000001
			iii = 0
			iii2 = 0	       
			do ii = 1,nva20
				if(filtre2(ii).eq.1)then
					iii2 = iii2 + 1
					vedc(k,iii2) = dble(vaxdc(ii))
				endif
			end do   
		else 
!------------------   censure a droite ou event recurr  c=0 
			if(icdc.eq.0)then
				cdc(k) = 0 
				iii = 0
				iii2 = 0
				do ii = 1,nva20
					if(filtre2(ii).eq.1)then
					iii2 = iii2 + 1
					vedc(k,iii2) = dble(vaxdc(ii))
					endif
				end do 
				t0dc(k) =  tt0dc 
				t1dc(k) = tt1dc    
			endif
		endif
		if (maxtdc.lt.t1dc(k))then
			maxtdc = t1dc(k)
		endif
	end do
	
	k = 0
	cptstr1 = 0
	cptstr2 = 0

!cccccccccccccccccccccccccccccccccc
! pour les donn�es recurrentes  
!cccccccccccccccccccccccccccccccccc     	 
	do i = 1,nsujet     !sur les observations
		tt0=tt00(i)
		tt1=tt10(i)
		ic=ic0(i)
		groupe=groupe0(i)	       	       	    
!-------------------	    
		do j=1,nva10	    
			vax(j)=vax0(i,j)
		enddo
!--------------	    	    
		if(tt0.gt.0.d0)then
			cptaux=cptaux+1
		endif  	                      
!-----------------------------------------------------
            
!     essai sans troncature
!     tt0=0.
!------------------   observation c=1 pour donn�es recurrentes
		if(ic.eq.1)then
			cpt = cpt + 1
			c(i)=1
			t0(i) = tt0 
			t1(i) = tt1  
			t1(i) = t1(i)
			g(i) = groupe
			nig(groupe) = nig(groupe)+1 ! nb d event recurr dans un groupe
			iii = 0
			iii2 = 0
!                  do ii = 1,ver		  
			do ii = 1,nva10
				if(filtre(ii).eq.1)then
					iii = iii + 1
					ve(i,iii) = dble(vax(ii)) !ici sur les observations
				endif
			end do   
		else 
!------------------   censure a droite  c=0 pour donn�es recurrentes
			if(ic.eq.0)then
				cptcens=cptcens+1
				c(i) = 0 
				iii = 0
				iii2 = 0
		!                     do ii = 1,ver
				do ii = 1,nva10
					if(filtre(ii).eq.1)then
					iii = iii + 1
					ve(i,iii) = dble(vax(ii))
					endif
				end do 
				t0(i) =  tt0 
				t1(i) = tt1 
				t1(i) = t1(i) 
				g(i) = groupe
			endif
		endif
		if (maxt.lt.t1(i))then
			maxt = t1(i)
		endif
	end do 
         	
	nsujet=i-1
	
	nz=nz0
	nz1=nz
	nz2=nz
	    
	if(nz.gt.20)then
		nz = 20
	endif
	if(nz.lt.4)then
		nz = 4
	endif
	 
	ndatemax=2*nsujet

	allocate(date(ndatemax),datedc(ndatemax),mm3(ndatemax),mm2(ndatemax) &
	,mm1(ndatemax),mm(ndatemax),im3(ndatemax),im2(ndatemax),im1(ndatemax),im(ndatemax))
        
!!!  DONNEES DECES

	mindc = 0.d0
	maxdc = maxtdc
	do i = 1,2*ng	 
		do k = 1,ng
			if((t0dc(k).ge.mindc))then
				if(t0dc(k).lt.maxdc)then
					maxdc = t0dc(k)
				endif
			endif
			if((t1dc(k).ge.mindc))then
				if(t1dc(k).lt.maxdc)then
					maxdc = t1dc(k)
				endif
			endif
		end do   
		aux(i) = maxdc
		mindc = maxdc + 1.d-12
		maxdc = maxtdc
	end do
	
	datedc(1) = aux(1)
	k = 1
	do i=2,2*ng
		if(aux(i).gt.aux(i-1))then
			k = k+1
			datedc(k) = aux(i)
		endif 
	end do 
	ndatedc = k      
!!         write(*,*)'** ndatemax,ndatemaxdc',ndatemax,ndatemaxdc	        
!--------------- zi- ----------------------------------

!      construire vecteur zi (des noeuds)

!!! DONNEES RECURRENTES

	min = 0.d0
	aux =0.d0
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
	ndate = k
	nzmax=nz+3
		
	allocate(zi(-2:nzmax))
	zi(-2) =0
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

!---------- affectation nt0dc,nt1dc DECES ----------------------------
	indictronqdc=0
	do k=1,ng 
		if(nig(k).eq.0.d0)then
			nb0recu = nb0recu + 1 !donne nb sujet sans event recu
		endif
		moyrecu =  moyrecu + dble(nig(k))
	
		if(t0dc(k).eq.0.d0)then
			nt0dc(k) = 0
		endif
		if(t0dc(k).ne.0.d0)then
			indictronqdc=1
		endif
		do j=1,ndatedc
			if(datedc(j).eq.t0dc(k))then
				nt0dc(k)=j
			endif
			if(datedc(j).eq.t1dc(k))then
				nt1dc(k)=j
			endif
		end do
	end do 

!---------- affectation nt0,nt1 RECURRENTS----------------------------

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
!AD:add argument:ndatedc
	call vecspliJ(n,ndate,ndatedc) 
!AD:end	
	allocate(m3m3(nzmax),m2m2(nzmax),m1m1(nzmax),mmm(nzmax),m3m2(nzmax) &
	,m3m1(nzmax),m3m(nzmax),m2m1(nzmax),m2m(nzmax),m1m(nzmax)) 
			
	call vecpenJ(n)
	
	npmax=np
	
	allocate(I_hess(npmax,npmax),H_hess(npmax,npmax),Hspl_hess(npmax,npmax) &
	,PEN_deri(npmax,1),hess(npmax,npmax),v((npmax*(npmax+3)/2)),I1_hess(npmax,npmax) &
	,H1_hess(npmax,npmax),I2_hess(npmax,npmax),H2_hess(npmax,npmax),HI2(npmax,npmax) & 
	,HIH(npmax,npmax),IH(npmax,npmax),HI(npmax,npmax),BIAIS(npmax,1))


               
!------- initialisation des parametres
                   
	do i=1,npmax
	b(i)=5.d-1
	end do

	b(np-nva-indic_alpha)=1.d0 ! pour theta
	b(np-nva)=1.d0

	call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaj)
	resOut=res
	
	if ((istop .eq. 4).or.(res .eq. -1.d9)) then
		goto 1400
	end if
	
	call multiJ(I_hess,H_hess,np,np,np,IH)
	call multiJ(H_hess,IH,np,np,np,HIH)   
		
	if(effet.eq.1.and.ier.eq.-1)then
	v((np-nva-indic_alpha)*(np-nva-indic_alpha+1)/2)=10.d10
	endif          
	res01(effet+1)=res
	 
! --------------  Lambda and survival estimates JRG January 05

	call distanceJ(nz1,nz2,b,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)      

	
	
	do ss=1,npmax
		do sss=1,npmax
			HIHOut(ss,sss) = HIH(ss,sss)
			H_hessOut(ss,sss)= H_hess(ss,sss)
		end do  
	end do
!AD:add LCV
!     calcul de la trace, pour le LCV (likelihood cross validation)
	traceLCV=0.d0
	call multiJ(H_hess,I_hess,np,np,np,HI)
	
	do i =1,np
		traceLCV = traceLCV + HI(i,i)
	end do 
	
	traceLCV = (traceLCV - resnonpen) / nsujet
	
	
1400 continue	
!AD:end

	deallocate(nig,cdc,nt0dc,nt1dc,t0dc,t1dc,aux1,aux2,res1,res3,mi) 
	deallocate(t0,t1,c,nt0,nt1,stra,g,aux)  
!AD: add deallocate  
	deallocate(mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc)
!AD:end
	deallocate(ve,vedc,filtre,filtre2,vax,vaxdc,date,datedc)   
	deallocate(mm3,mm2,mm1,mm,im3,im2,im1,im,zi)
	deallocate(m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m) 
	deallocate(I_hess,H_hess,Hspl_hess,PEN_deri,hess,v,I1_hess,H1_hess &
	,I2_hess,H2_hess,HI2,HIH,IH,HI,BIAIS) 
		
	return
				
	end subroutine joint
      

!========================== VECSPLI ==============================
!AD:add argument:ndatedc 
	subroutine vecspliJ(n,ndate,ndatedc) 
!AD:end	
	use tailles
!AD:	
	use comon,only:date,datedc,zi,mm3,mm2,mm1,mm,im3,im2,im1,im &
	,mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc
!AD:end	
	IMPLICIT NONE
	
	integer,intent(in)::n,ndate,ndatedc
	integer::i,j,k
	double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2
      
!----------  calcul de u(ti) :  STRATE1 ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)
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
		im2(i) = (0.25d0*hht*mm2(i))+(h3m*mm1(i)*0.25d0) &
			+(h4*mm(i)*0.25d0)
		im1(i) = (htm*mm1(i)*0.25d0)+(h4*mm(i)*0.25d0)
		im(i)  = ht*mm(i)*0.25d0

	end do
!AD: add for death 
!----------  calcul de u(ti) :  STRATE2 ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

	do i=1,ndatedc-1
		do k = 2,n-2
			if ((datedc(i).ge.zi(k-1)).and.(datedc(i).lt.zi(k)))then
				j = k-1
			endif
		end do 
		ht = datedc(i)-zi(j)
		htm= datedc(i)-zi(j-1)
		h2t= datedc(i)-zi(j+2)
		ht2 = zi(j+1)-datedc(i)
		ht3 = zi(j+3)-datedc(i)
		hht = datedc(i)-zi(j-2)
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
		mm3dc(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
		mm2dc(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
		*ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
		mm1dc(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
		h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
		mmdc(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
		im3dc(i) = (0.25d0*(datedc(i)-zi(j-3))*mm3dc(i))+(0.25d0*hh2 &
		*mm2dc(i))+(0.25d0*h3m*mm1dc(i))+(0.25d0*h4*mmdc(i))
		im2dc(i) = (0.25d0*hht*mm2dc(i))+(h3m*mm1dc(i)*0.25d0) &
			+(h4*mmdc(i)*0.25d0)
		im1dc(i) = (htm*mm1dc(i)*0.25d0)+(h4*mmdc(i)*0.25d0)
		imdc(i)  = ht*mmdc(i)*0.25d0

	end do
!AD:end	    
	end subroutine vecspliJ  

!========================== VECPEN ==============================
	subroutine vecpenJ(n) 
	
	use tailles
	
	use comon,only:date,datedc,zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m
	
	IMPLICIT NONE
	
	integer,intent(in)::n
	integer::i
	double precision::h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2,a3,a2,b2 &
	,c2,a1,b1,c1,a0,x3,x2,x


!*********************************************************************
         
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
		m2m2(i) = m2m2(i) + 64.d0*(((3.d0*x3-(3.d0*x2*(zi(i+2)  &
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
		m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i-1)+zi(i) &     
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
		*zi(i)+zi(i-1)+zi(i+2)))+x*(zi(i-1)*zi(i+3)+2.d0*zi(i-1) &    
		*zi(i)+zi(i+3)*zi(i)+2.d0*zi(i)*zi(i)+zi(i+2)*zi(i+3) &
		+2.d0*zi(i+2)*zi(i)))/(b1*c1))
		mmm(i) = (192.d0*h/(h4*h3*h2*h4*h3*h2))
		m3m2(i) = 192.d0*(((-x3+(0.5d0*x2*(5.d0*zi(i+1)+zi(i-2) &
		))-x*(2.d0*zi(i+1)*zi(i+1)+zi(i+1)*zi(i-2)))/(a3*a2)) &
		+((-x3+(0.5d0*x2*(4.d0*zi(i+1)+zi(i-1)+zi(i+2)))-x* &
		(zi(i+1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*b2)) &
		+((-x3+(0.5d0*x2*(3.d0*zi(i+1)+2.d0*zi(i+2)+zi(i)))-x* &
		(2.d0*zi(i+1)*zi(i+2)+zi(i+1)*zi(i)))/(a3*c2)) )
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

	end subroutine vecpenJ

!========================          FUNCPA NEW         ====================
	double precision function funcpaJ(b,np,id,thi,jd,thj,k0)
	
	use tailles
	use comon
	use additiv,only:aux1,aux2
	
	IMPLICIT NONE

! *** NOUVELLLE DECLARATION F90 :
	
	integer,intent(in)::id,jd,np
	double precision,dimension(np),intent(in)::b
	double precision,dimension(2),intent(in)::k0
	double precision,intent(in)::thi,thj
	
	integer::n,i,j,k,vj,ig,choix
	integer,dimension(ngmax)::cpt
	double precision::pe1,pe2,sum,inv,som1,som2,res,vet,vet2,h1
	
	double precision,dimension(-2:npmax):: the1,the2
	double precision,dimension(np)::bh
	double precision,dimension(ngmax)::res2,res1dc,res2dc &
	,res3dc,integrale1,integrale2,integrale3
!AD: for death,change dimension 
	double precision,dimension(ndatemax)::dut1
	double precision,dimension(ndatemaxdc)::dut2
	logical::isnan
!AD:end
	double precision,dimension(0:ndatemax)::ut1
	double precision,dimension(0:ndatemaxdc)::ut2
	double precision::int,gammaJ

	choix=0
	ig=0
	k=0
	vj=0
	n=0
	j=0
	do i=1,np
	bh(i)=b(i)
	end do 

	if (id.ne.0) bh(id)=bh(id)+thi
	if (jd.ne.0) bh(jd)=bh(jd)+thj    



	n = (np-nva-effet-indic_ALPHA)/nst

	do i=1,n
		the1(i-3)=(bh(i))*(bh(i))
		j = n+i 
		if (nst.eq.2) then
			the2(i-3)=(bh(j))*(bh(j))
		endif
	end do

	
	if(effet.eq.1) then
		theta = bh(np-nva-indic_ALPHA)*bh(np-nva-indic_ALPHA)
		alpha = bh(np-nva)
	endif

!----------  calcul de ut1(ti) et ut2(ti) ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

!AD:modify
	dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))
	dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))
	
	ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))
	ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))
	
	ut1(0) = 0.d0
	ut2(0) = 0.d0
	 
!//// NEW AMADOU vvv :
!--- strate1
	som1 = 0.d0
	vj = 0
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

		if(nst.eq.2)then
			ut2(i) = som2 +(the2(j-3)*im3dc(i))+(the2(j-2)*im2dc(i)) &
			+(the2(j-1)*im1dc(i))+(the2(j)*imdc(i))
			dut2(i) = (the2(j-3)*mm3dc(i))+(the2(j-2)*mm2dc(i)) &
			+(the2(j-1)*mm1dc(i))+(the2(j)*mmdc(i))
		endif
            
	end do

!-------------fin strate2  
	i = n-2
	h1 = (zi(i)-zi(i-1))
	
	ut1(ndate)=som1+the1(i-4)+the1(i-3)+the1(i-2)+the1(i-1)
	dut1(ndate) = (4.d0*the1(i-1)/h1)

		
	ut2(ndatedc)=som2+the2(i-4)+the2(i-3)+the2(i-2)+the2(i-1)!am the1(i-4)
	dut2(ndatedc) = (4.d0*the2(i-1)/h1)
         
!//// fin NEW AMADOU-vvv
!AD:end
!-------------------------------------------------------
!---------- calcul de la vraisemblance ------------------
!---------------------------------------------------------

!---- avec ou sans variable explicative  ------

	do k=1,ng
		res1(k) = 0.d0
		res2(k) = 0.d0
		res3(k) = 0.d0
		res1dc(k) = 0.d0
		res2dc(k) = 0.d0
		res3dc(k) = 0.d0
		cpt(k) = 0
		integrale1(k) = 0.d0
		integrale2(k) = 0.d0
		integrale3(k) = 0.d0
		aux1(k)=0.d0
		aux2(k)=0.d0
	end do

!*******************************************         
!-----avec un effet aleatoire dans le modele
!*********************************************

	inv = 1.d0/theta

!ccccccccccccccccccccccccccccccccccccccccc
!     pour les donnees recurrentes
!ccccccccccccccccccccccccccccccccccccccccc

	
	do i=1,nsujet 
		cpt(g(i))=cpt(g(i))+1  
		if(nva1.gt.0)then
			vet = 0.d0   
			do j=1,nva1
				vet =vet + bh(np-nva+j)*dble(ve(i,j))
!      write(*,*)'*** funcpa vet',vet,ve(i,j),i,j
			end do
			vet = dexp(vet)
		else
			vet=1.d0
		endif
            
		if((c(i).eq.1))then
			res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet) 
!AD:
			if (isnan(res2(g(i)))) then
				funcpaj=-1.d9
				goto 600	
			end if
!AD: 
		endif  
!     nouvelle version
		res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet            
!     modification pour nouvelle vraisemblance / troncature:
		res3(g(i)) = res3(g(i)) + ut1(nt0(i))*vet 
	end do
	 
!           stop
!ccccccccccccccccccccccccccccccccccccccccc
! pour le deces 
!ccccccccccccccccccccccccccccccccccccccccc 

	do k=1,ng  
		if(nva2.gt.0)then
			vet2 = 0.d0   
			do j=1,nva2
				vet2 =vet2 + bh(np-nva2+j)*dble(vedc(k,j))
			end do
			vet2 = dexp(vet2)
		else
			vet2=1.d0
		endif
		if(cdc(k).eq.1)then
			res2dc(k) = dlog(dut2(nt1dc(k))*vet2)  
!AD:
			if (isnan(res2dc(k))) then
				funcpaj=-1.d9
				goto 600	
			end if
!AD: 
		endif 
             
! pour le calcul des integrales / pour la survie, pas les donn�es recurrentes:
		aux1(k)=ut2(nt1dc(k))*vet2
		aux2(k)=aux2(k)+ut2(nt0(k))*vet2 !vraie troncature
	end do

!**************INTEGRALES ****************************
	do ig=1,ng 
		auxig=ig 
		choix = 3  
		call gaulagJ(int,choix)
		integrale3(ig) = int !moins bon
	end do
!************* FIN INTEGRALES **************************
                      
	res = 0.d0 
	do k=1,ng  
		sum=0.d0
		if(cpt(k).gt.0)then
			if(theta.gt.(1.d-8)) then
!cccc ancienne vraisemblance : pour calendar sans vrai troncature cccccccc
                   
				res= res + res2(k) &
!--      pour le deces:
				+ res2dc(k)  &
				- gammaJ(1./theta)-dlog(theta)/theta  &
				+ dlog(integrale3(k))
!AD:
				if (isnan(res)) then
					funcpaj=-1.d9
					goto 600	
				end if
!AD: 
			else
!*************************************************************************
!     developpement de taylor d ordre 3
!*************************************************************************
!                   write(*,*)'************** TAYLOR *************'                   
				res= res + res2(k) &
				+ res2dc(k)  &
				- gammaJ(1./theta)-dlog(theta)/theta  &
				+ dlog(integrale3(k)) 
!AD:
				if (isnan(res)) then
					funcpaj=-1.d9
					goto 600	
				end if
!AD: 
			endif
		endif 
	end do
               
!---------- calcul de la penalisation -------------------

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
	
	funcpaJ = res 
!Ad:
600     continue

	if (isnan(funcpaj)) then
		funcpaj=-1.d9
	end if
	
!AD:
	return
	
	
	
	end function funcpaJ





!==========================  DISTANCE   =================================
         
	subroutine distanceJ(nz1,nz2,b,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)
	
	use tailles
	use comon,only:date,datedc,zi,t0,t1,t0dc,t1dc,c,cdc &
	,nt0,nt1,nt0dc,nt1dc,nsujet,nva,nva1,nva2,ndate,ndatedc,nst &
	,PEN_deri,I_hess,H_hess,Hspl_hess,hess
	IMPLICIT NONE  
	
	
	integer,intent(in):: nz1,nz2
	double precision ,dimension(npmax),intent(in):: b
	integer::i,j,n,k,l      
	double precision::x1,x2,h,su,bsup,binf,lam,lbinf,lbsup
	double precision ,dimension(npmax,npmax)::hes1,hes2
	double precision ,dimension(-2:npmax):: the1,the2     
	double precision,dimension(99)::x1Out,x2Out
	double precision,dimension(99,3):: lamOut,suOut,lam2Out,su2Out 
		
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
            
	do i=1,99
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
		suOut(i,2)=binf
		suOut(i,3)=bsup  
	              
		if(nst.eq.2)then
			if(i.ne.1)then
				x2 = x2 + h 
			endif
 
			call cospJ(x2,the2,nz2+2,hes2,zi,binf,su,bsup,lbinf,lam,lbsup)
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
		
	end subroutine distanceJ





!==========================  SUSP  ====================================
	subroutine suspJ(x,the,n,su,lam,zi)
	
	use tailles
	
	IMPLICIT NONE 
	
	integer,intent(in)::n
	double precision,intent(out)::lam,su
	double precision,dimension(-2:npmax),intent(in)::zi,the
	double precision,intent(in)::x
	integer::j,k,i
	double precision::ht,ht2,h2,som,htm,h2t,h3,h2n,hn, &
	im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2 &
	,h,gl,hh

	gl=0.d0
	som = 0.d0
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
	 
	end subroutine suspJ

!==========================  COSP  ====================================
! calcul les points pour les fonctions 
! et leur bandes de confiance

	subroutine cospJ(x,the,n,y,zi,binf,su,bsup,lbinf,lam,lbsup)
	
	use tailles
	
	IMPLICIT NONE
	
	integer,intent(in)::n
	double precision,intent(in)::x
	double precision,intent(out)::lam,su
	double precision,intent(out)::binf,bsup,lbinf,lbsup
	double precision,dimension(npmax,npmax),intent(in)::y
	double precision,dimension(-2:npmax),intent(in)::the,zi
	integer::j,k,i
	double precision::ht,ht2,h2,som,pm,htm,h2t,h3,h2n,hn, &
	im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2, &
	h,gl,hh
      
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

	call confJ(x,j,n,y,pm,zi)

	binf = dexp(-gl - 1.96d0*pm)
	su  = dexp(-gl)
	bsup = dexp(-gl + 1.96d0*pm)

	call conf1J(x,j,n,y,pm,zi)
	lbinf = lam - 1.96d0*pm
	lbsup = lam + 1.96d0*pm
!         write(*,*)'lbinf apres conf1',lbinf,lam,pm

	return

	end subroutine cospJ
	 
	 
!=====================  CONF1  =============================


	subroutine  conf1J(x,ni,n,y,pm,zi)
	
	use tailles
	
	IMPLICIT NONE  
	
	integer,intent(in)::ni,n
	double precision,intent(in)::x
	double precision,dimension(-2:npmax),intent(in)::zi
	double precision,dimension(npmax,npmax),intent(in)::y
	double precision,intent(out)::pm
	integer::i,j
	double precision::res,mmspJ
	double precision,dimension(npmax)::vecti,aux

      
           
	do i=1,n
		vecti(i) = mmspJ(x,ni,i,zi)
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
	
	res=-res 
	pm = dsqrt(res)
	
	end subroutine  conf1J
	 
!=====================  CONF  =============================

	subroutine  confJ(x,ni,n,y,pm,zi)
	
	use tailles
	
	IMPLICIT NONE  
	
	integer,intent(in)::ni,n
	double precision,intent(in)::x
	double precision,dimension(-2:npmax),intent(in)::zi
	double precision,dimension(npmax,npmax),intent(in)::y
	double precision,intent(out)::pm
	integer::i,j
	double precision::res,ispJ
	double precision,dimension(52)::vecti,aux
	
	do i=1,n
	vecti(i) = ispJ(x,ni,i,zi)
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
	res=-res
	pm = dsqrt(res)
               
	end subroutine  confJ


!==========================   ISP   ==================================

	double precision function ispJ(x,ni,ns,zi)
	
	use tailles
	
	IMPLICIT NONE  
	
	integer,intent(in)::ni,ns
	double precision,intent(in)::x
	double precision,dimension(-2:npmax),intent(in)::zi
	double precision::val,mmspJ



	if(x.eq.zi(ni))then
		if(ni.le.ns-3)then
			val = 0.d0
			else
				if(ni.le.ns-2)then
					val = ((zi(ni)-zi(ni-1))*mmspJ(x,ni,ns,zi))*0.25d0
				else
					if (ni.eq.ns-1)then
						val = ((zi(ni)-zi(ni-2))*mmspJ(x,ni,ns,zi)+ &
						(zi(ni+3)-zi(ni-1))*mmspJ(x,ni,ns+1,zi))*0.25d0
					else
						if(ni.eq.ns)then
							val = ((zi(ni)-zi(ni-3))*mmspJ(x,ni,ns,zi)+ &
							(zi(ni+2)-zi(ni-2))*mmspJ(x,ni,ns+1,zi) &
							+(zi(ni+3)-zi(ni-1))*mmspJ(x,ni,ns+2,zi))*0.25d0
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
				val = (x-zi(ni))*mmspJ(x,ni,ns,zi)*0.25d0
			else  
				if(ni.eq.ns-2)then
					val = ((x-zi(ni-1))*mmspJ(x,ni,ns,zi)+ &
					(zi(ni+4)-zi(ni))*mmspJ(x,ni,ns+1,zi))*0.25d0
				else   
					if (ni.eq.ns-1)then
						val =((x-zi(ni-2))*mmspJ(x,ni,ns,zi)+ &
						(zi(ni+3)-zi(ni-1))*mmspJ(x,ni,ns+1,zi) &
						+(zi(ni+4)-zi(ni))*mmspJ(x,ni,ns+2,zi))*0.25d0
					else
						if(ni.eq.ns)then
							val =((x-zi(ni-3))*mmspJ(x,ni,ns,zi)+ &
							(zi(ni+2)-zi(ni-2))*mmspJ(x,ni,ns+1,zi) &
							+(zi(ni+3)-zi(ni-1))*mmspJ(x,ni,ns+2,zi) &
							+(zi(ni+4)-zi(ni))*mmspJ(x,ni,ns+3,zi))*0.25d0
						else
							val = 1.d0
						endif
					endif
				endif
			endif
		endif 
	endif
	  
	ispJ = val
	
	return
	
	end function ispJ
	
!==========================  MMSP   ==================================

	double precision function mmspJ(x,ni,ns,zi)
	
	use tailles
	
	IMPLICIT NONE 
	
	integer,intent(in)::ni,ns
	double precision,intent(in)::x
	double precision,dimension(-2:npmax),intent(in)::zi
	double precision::val

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
					-zi(ni))*(zi(ni+1)-zi(ni))*(zi(ni+2)-zi(ni-1))) &
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
						+((4.d0*((zi(ni+2)-x)*(zi(ni+2)-x) &
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

	mmspJ = val
	
	return
	
	end function mmspJ

!================== multiplication de matrice  ==================

! multiplie A par B avec le resultat dans C

	subroutine multiJ(A,B,IrowA,JcolA,JcolB,C)
!     remarque :  jcolA=IrowB
	use tailles
	
	IMPLICIT NONE
	
	integer,intent(in)::IrowA,JcolA,JcolB
	double precision,dimension(npmax,npmax),intent(in):: A,B
	double precision,dimension(npmax,npmax),intent(out)::C       
	integer::i,j,k
	double precision::sum

	do I=1,IrowA
		do J=1,JcolB
			sum=0
			do K=1,JcolA
				sum=sum+A(I,K)*B(K,J)
			end do
			C(I,J)=sum
		end do
	end do
	
	return
	
	end subroutine multiJ
!====================================================================
!============================    GAMMA      ==============================

!       function qui calcule le log de  Gamma
	double precision function gammaJ(xx)
	
	use donnees,only:cof,stp,half,one,fpf 
	
	implicit none
	
	integer::j
	double precision,intent(in)::xx
            
	
	double precision::x,tmp,ser
	
	x = xx - one
	tmp = x + fpf
	tmp = (x+half)*dlog(tmp) - tmp
	ser = one
	do j = 1,6
		x = x + one
		ser = ser + cof(j)/x
	end do
	gammaJ = tmp + dlog(stp*ser)
	
	return
	
	end function gammaJ
   
!======================================================

!==================================================================
!==================================================================
! 20 points / sur (0,+infty)
	SUBROUTINE gaulagJ(ss,choix) 
	
	use tailles
	use comon,only:auxig
	use donnees,only:w,x
	
	IMPLICIT NONE
	
	integer,intent(in)::choix
	double precision,intent(out):: ss
	double precision ::auxfunca,func1J,func2J,func3J
	external :: func1J,func2J,func3J
! gauss laguerre
! func1 est l int�grant, ss le r�sultat de l integrale sur 0 ,  +infty
	integer :: j,nan

      
	ss=0.d0 
! Will be twice the average value of the function,since the ten
! wei hts (five numbers above each used twice) sum to 2.
	do j=1,20
		if (choix.eq.1) then !integrale 1
			auxfunca=func1J(x(j))
			ss = ss+w(j)*(auxfunca)
			if(ss.eq.nan)then 
!               write(*,*)'----int1',ss ,w(j),auxfunca,j
!               stop
			endif

		else                   !choix=2, survie marginale, vraie troncature
			if (choix.eq.2) then 
				auxfunca=func2J(x(j))
				ss = ss+w(j)*(auxfunca)
			else                   !choix=3, AG model
				if (choix.eq.3) then
					auxfunca=func3J(x(j))
					ss = ss+w(j)*(auxfunca)
				endif
			endif
		endif
	end do
	
	return
	
	END SUBROUTINE gaulagJ
  
!================================================

	double precision function func1J(frail) 
! calcul de l integrant, pour un effet aleatoire donn� frail et un groupe donne auxig (cf funcpa)      
	use tailles
	
	use comon,only:auxig,alpha,theta,res1,res3,mi
	use additiv,only:aux1,aux2
	
	IMPLICIT NONE

	double precision,intent(in)::frail

!============================================================
	
	func1J = (frail**(alpha*mi(auxig)+1./theta-1.))* &
		dexp(-(frail**alpha) *aux1(auxig))*dexp(-frail/theta)
	
	return
	
	end function func1J
!==================================================================

!================================================

	double precision function func2J(frail) 
! calcul de l integrant, pour un effet aleatoire donn� frail et un groupe donne auxig (cf funcpa)      
	use tailles
	
	use comon,only:auxig,ALPHA,theta,mi,res1,res3
	use additiv,only:aux1,aux2
	IMPLICIT NONE
	
	double precision,intent(in)::frail
	double precision::gammaJ
            
!============================================================
      
	func2J = dexp(-(frail**alpha)*aux2(auxig))*dexp(-frail/theta)*(frail) &
	/(exp(gammaJ(1.d0/theta))*(theta**(1./theta)))
		
	return
	
	end function func2J
	
!==================================================================

	double precision function func3J(frail) 
! calcul de l integrant, pour un effet aleatoire donn� frail et un groupe donne auxig (cf funcpa)      
      	use tailles
	use comon,only:g,nig,auxig,alpha,theta,res1,res3,t0,t1,t0dc,t1dc,c, &
        cdc,nt0,nt1,nt0dc,nt1dc,nsujet,nva,nva1,nva2,ndate,ndatedc,nst,AG 
	use additiv,only:aux1,aux2   
	
	IMPLICIT NONE 

	double precision,intent(in)::frail
!      double precision func3



	
	func3J = (nig(auxig)+ alpha*cdc(auxig)+ 1./theta-1.)*dlog(frail) &
		- frail*(res1(auxig)-res3(auxig)) &!res3=0 si AG=0
		- (frail**alpha)*(aux1(auxig))- frail/theta
	
	func3J = exp(func3J)
	
	return
	
	end function func3J
!==================================================================
!AD: IS NAN

	logical function isnan(x)

	implicit none
	
	double precision,intent(in)::x
	
	if (x .ne. x) then
		isnan=.true.
	else
		isnan=.false.
	end if

	end function isnan


!AD:end
!====================================================================
