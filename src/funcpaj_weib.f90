
!!!!_____________________________________________________
!========================          FUNCPA NEW         ====================
	double precision function funcpaj_weib(b,np,id,thi,jd,thj,k0)
	
	use tailles
	use comon
        use residusM
		
	implicit none

! *** NOUVELLLE DECLARATION F90 :
	
	integer,intent(in)::id,jd,np
	double precision,dimension(np),intent(in)::b
	double precision,dimension(2)::k0
	double precision,intent(in)::thi,thj
	
	integer::n,i,j,k,vj,ig,choix
	integer,dimension(ngmax)::cpt
	double precision::sum,res,vet,vet2
	
	double precision,dimension(np)::bh
	double precision,dimension(ngmax)::res2,res1dc,res2dc &
	,res3dc,integrale1,integrale2,integrale3
	double precision::int,gammaJ
	
	kkapa=k0
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

	betaR= bh(1)**2
	etaR= bh(2)**2
	betaD= bh(3)**2
	etaD= bh(4)**2
	
	if(effet.eq.1) then
		theta = bh(np-nva-indic_ALPHA)*bh(np-nva-indic_ALPHA)
		alpha = bh(np-nva)
	endif

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

!	inv = 1.d0/theta

!ccccccccccccccccccccccccccccccccccccccccc
!     pour les donnees recurrentes
!ccccccccccccccccccccccccccccccccccccccccc

	
	do i=1,nsujet 
		cpt(g(i))=cpt(g(i))+1  
		if(nva1.gt.0)then
			vet = 0.d0   
			do j=1,nva1
				vet =vet + bh(np-nva+j)*dble(ve(i,j))
			end do
			vet = dexp(vet)
		else
			vet=1.d0
		endif
            
		if((c(i).eq.1))then
			
			res2(g(i)) = res2(g(i))+(betaR-1.d0)*dlog(t1(i))+dlog(betaR)-betaR*dlog(etaR)+dlog(vet) 
			if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
				funcpaj_weib=-1.d9
				goto 123
			end if	
		endif  
!     nouvelle version
		res1(g(i)) = res1(g(i))+((t1(i)/etaR)**betaR)*vet 
 		if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
			funcpaj_weib=-1.d9
			goto 123
		end if	         
!     modification pour nouvelle vraisemblance / troncature:
		res3(g(i)) = res3(g(i))+((t0(i)/etaR)**betaR)*vet
		if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
			funcpaj_weib=-1.d9
			goto 123
		end if	
	!	write(*,*)'***res2',res3(g(i)),vet,i,g(i) 
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
			res2dc(k) = (betaD-1.d0)*dlog(t1dc(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet2) 
			if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
				funcpaj_weib=-1.d9
				goto 123
			end if	
		endif 
             
! pour le calcul des integrales / pour la survie, pas les donn�es recurrentes:
		aux1(k)=((t1dc(k)/etaD)**betaD)*vet2
		if ((aux1(k).ne.aux1(k)).or.(abs(aux1(k)).ge. 1.d30)) then
			funcpaj_weib=-1.d9
			goto 123
		end if		
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
				+ res2dc(k)- gammaJ(1./theta)-dlog(theta)/theta+dlog(integrale3(k))
			else
!*************************************************************************
!     developpement de taylor d ordre 3
!*************************************************************************
!                   write(*,*)'************** TAYLOR *************'                   
				res= res + res2(k)+res2dc(k)-gammaJ(1./theta)-dlog(theta)/theta  &
				+ dlog(integrale3(k)) 
			endif
			if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
				funcpaj_weib=-1.d9
				goto 123
			end if	
		endif 
	end do
               
!---------- calcul de la penalisation -------------------
	
	if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
		funcpaj_weib =-1.d9
		Rrec = 0.d0
		Nrec = 0.d0
		Rdc = 0.d0
		Ndc = 0.d0
		goto 123

	else
		funcpaj_weib = res 
		
		do k=1,ng
			Rrec(k)=res1(k)
			Nrec(k)=nig(k)
			Rdc(k)=aux1(k)
			Ndc(k)=cdc(k)
		end do
	end if

!Ad:
123     continue

	return
		
	end function funcpaj_weib
	
!=================================================================================================
