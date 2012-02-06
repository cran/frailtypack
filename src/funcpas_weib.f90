

!========================          FUNCPA_WEIB          ====================
	double precision function funcpas_weib(b,np,id,thi,jd,thj,k0)
	use tailles
	use comon,only:t0,t1,c,nsujet,nva, &
	nst,stra,ve,effet,ng,g,nig,AG,etaR,etaD,betaR,betaD,kkapa,theta
        use residusM

	implicit none
	
! *** NOUVELLLE DECLARATION F90 :
	
	integer::nb,np,id,jd,i,j,k,cptg,l
	integer,dimension(ngmax)::cpt
	double precision::thi,thj,dnb,sum,inv,res,vet
	double precision,dimension(np)::b,bh
	double precision,dimension(ngmax)::res1,res2,res3
	double precision,dimension(2)::k0 

	kkapa=k0
	j=0
	theta=0.d0
	do i=1,np
		bh(i)=b(i)
	end do 

	if (id.ne.0) bh(id)=bh(id)+thi
	if (jd.ne.0) bh(jd)=bh(jd)+thj    
	
	if (nst == 1) then
		betaR= bh(1)**2
		etaR= bh(2)**2
		etaD= 0.d0
		betaD= 0.d0
	else
		betaR= bh(1)**2
		etaR= bh(2)**2
		betaD= bh(3)**2
		etaD= bh(4)**2		
	end if

	if(effet.eq.1) then
		theta = bh(np-nva)*bh(np-nva)
	endif

!-------------------------------------------------------
!--------- calcul de la vraisemblance ------------------
!--------------------------------------------------------

!--- avec ou sans variable explicative  ------cc

	do k=1,ng
		res1(k) = 0.d0
		res2(k) = 0.d0
		res3(k) = 0.d0
		cpt(k) = 0
	end do

!*******************************************     
!---- sans effet aleatoire dans le modele
!*******************************************     

	if (effet.eq.0) then
		do i=1,nsujet
			cpt(g(i))=cpt(g(i))+1
			
			if(nva.gt.0)then
				vet = 0.d0   
				do j=1,nva
					vet =vet + bh(np-nva+j)*dble(ve(i,j))
				end do
				vet = dexp(vet)
			else
				vet=1.d0
			endif
			
			if((c(i).eq.1).and.(stra(i).eq.1))then
				res2(g(i)) = res2(g(i))+(betaR-1.d0)*dlog(t1(i))+dlog(betaR)-betaR*dlog(etaR)+dlog(vet)
			endif  
	
			if((c(i).eq.1).and.(stra(i).eq.2))then
				res2(g(i)) = res2(g(i))+(betaD-1.d0)*dlog(t1(i))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet)
			endif
	               if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                          funcpas_weib=-1.d9
                          goto 123
                       end if	
			if(stra(i).eq.1)then
				res1(g(i)) = res1(g(i)) + ((t1(i)/etaR)**betaR)*vet - ((t0(i)/etaR)**betaR)*vet  
				RisqCumul(i) = ((t1(i)/etaR)**betaR)*vet
			endif
	
			if(stra(i).eq.2)then
				res1(g(i)) = res1(g(i)) + ((t1(i)/etaD)**betaD)*vet - ((t0(i)/etaD)**betaD)*vet 
				RisqCumul(i) = ((t1(i)/etaD)**betaD)*vet
			endif
	               if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                          funcpas_weib=-1.d9
                          goto 123
                       end if			
		end do       
		res = 0.d0         
		cptg = 0
		
! k indice les groupes
		do k=1,ng   
			if(cpt(k).gt.0)then
				nb = nig(k)
				dnb = dble(nig(k))               
				res = res-res1(k)+res2(k) 
				cptg = cptg + 1 
			       if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
		                  funcpas_weib=-1.d9
		                  goto 123
		               end if	
			endif 
		end do

!*******************************************         
!----avec un effet aleatoire dans le modele
!*********************************************

	else
!      write(*,*)'AVEC EFFET ALEATOIRE'
		inv = 1.d0/theta
!     i indice les sujets
		do i=1,nsujet 
			
			cpt(g(i))=cpt(g(i))+1 
		
			if(nva.gt.0)then
				vet = 0.d0   
				do j=1,nva
					vet =vet + bh(np-nva+j)*dble(ve(i,j))
				end do
				vet = dexp(vet)
			else
				vet=1.d0
			endif
			if((c(i).eq.1).and.(stra(i).eq.1))then
				res2(g(i)) = res2(g(i))+(betaR-1.d0)*dlog(t1(i))+dlog(betaR)-betaR*dlog(etaR)+dlog(vet)
			endif  
			if((c(i).eq.1).and.(stra(i).eq.2))then
				res2(g(i)) = res2(g(i))+(betaD-1.d0)*dlog(t1(i))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet)
			endif  
	               if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                          funcpas_weib=-1.d9
                          goto 123
                       end if	
			if(stra(i).eq.1)then
				res1(g(i)) = res1(g(i)) + ((t1(i)/etaR)**betaR)*vet
			endif
			
			if(stra(i).eq.2)then
				res1(g(i)) = res1(g(i)) + ((t1(i)/etaD)**betaD)*vet
			endif
	               if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                          funcpas_weib=-1.d9
                          goto 123
                       end if	
! modification pour nouvelle vraisemblance / troncature:
			if(stra(i).eq.1)then
				res3(g(i)) = res3(g(i)) + ((t0(i)/etaR)**betaR)*vet
			endif
			
			if(stra(i).eq.2)then
				res3(g(i)) = res3(g(i)) + ((t0(i)/etaD)**betaD)*vet
			endif
	               if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                          funcpas_weib=-1.d9
                          goto 123
                       end if			
		end do 

		res = 0.d0
		cptg = 0
!     gam2 = gamma(inv)
! k indice les groupes
		do k=1,ng  
			sum=0.d0
			if(cpt(k).gt.0)then
				nb = nig(k)
				dnb = dble(nig(k))
				
				if (dnb.gt.1.d0) then
					do l=1,nb
						sum=sum+dlog(1.d0+theta*dble(nb-l))
					end do
				endif
				if(theta.gt.(1.d-5)) then
!ccccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
					if(AG.EQ.1)then
						res= res-(inv+dnb)*dlog(theta*(res1(k)-res3(k))+1.d0) &
						+ res2(k) + sum  
!ccccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
					else
						res = res-(inv+dnb)*dlog(theta*res1(k)+1.d0)  &
						+(inv)*dlog(theta*res3(k)+1.d0)+ res2(k) + sum  
					endif
				else              
!     developpement de taylor d ordre 3
!                   write(*,*)'************** TAYLOR *************'
!cccc ancienne vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
					if(AG.EQ.1)then
						res = res-dnb*dlog(theta*(res1(k)-res3(k))+1.d0) &
						-(res1(k)-res3(k))*(1.d0-theta*(res1(k)-res3(k))/2.d0 &
						+theta*theta*(res1(k)-res3(k))*(res1(k)-res3(k))/3.d0)+res2(k)+sum

!cccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
					else
						res = res-dnb*dlog(theta*res1(k)+1.d0)-res1(k)*(1.d0-theta*res1(k) &
						/2.d0+theta*theta*res1(k)*res1(k)/3.d0) &
						+res2(k)+sum &
						+res3(k)*(1.d0-theta*res3(k)/2.d0 &
							+theta*theta*res3(k)*res3(k)/3.d0)
					endif
				endif 
			       if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
		                  funcpas_weib=-1.d9
		                  goto 123
		               end if	
			endif 
		end do


       	endif !fin boucle effet=0

!--------- calcul de la penalisation -------------------

!    Changed JRG 25 May 05
	if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
		funcpas_weib=-1.d9
		goto 123
	end if	

	funcpas_weib = res 

	do k=1,ng
		cumulhaz(k)=res1(k)
	end do

123     continue

	return

	end function funcpas_weib


!=================================================================================================
