

!========================          FUNCPA_SPLINES          ====================
	double precision function funcpassplines(b,np,id,thi,jd,thj,k0)
	use tailles
	use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m, &
	mm3,mm2,mm1,mm,im3,im2,im1,im,date,zi,t0,t1,c,nt0,nt1,nsujet,nva,ndate, &
	nst,stra,ve,pe,effet,nz1,nz2,ng,g,nig,AG,resnonpen,theta
        use residusM
	
	
	implicit none
	
! *** NOUVELLLE DECLARATION F90 :
	
	integer::nb,n,np,id,jd,i,j,k,vj,cptg,l
	integer,dimension(ngmax)::cpt
	double precision::thi,thj,pe1,pe2,dnb,sum,inv,som1,som2,res,vet,h1
	double precision,dimension(-2:npmax)::the1,the2
	double precision,dimension(np)::b,bh
	double precision,dimension(ngmax)::res1,res2,res3
	double precision,dimension(2)::k0 
	double precision,dimension(ndatemax)::dut1,dut2
	double precision,dimension(0:ndatemax)::ut1,ut2
	
	j=0
	theta=0.d0
	do i=1,np
		bh(i)=b(i)
	end do 
	dut1=0.d0
	dut2=0.d0
	ut1=0.d0
	ut2=0.d0
	if (id.ne.0) bh(id)=bh(id)+thi
	if (jd.ne.0) bh(jd)=bh(jd)+thj    
	
	n = (np-nva-effet)/nst
          
	do i=1,n
		the1(i-3)=(bh(i))*(bh(i))
		j = n+i 
		if (nst.eq.2) then
			the2(i-3)=(bh(j))*(bh(j))
		endif
	end do


	if(effet.eq.1) then
		theta = bh(np-nva)*bh(np-nva)
	endif
!---------  calcul de ut1(ti) et ut2(ti) ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

	vj = 0
	som1 = 0.d0
	som2 = 0.d0
	dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))

	dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))
	ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))
	ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))
	ut1(0) = 0.d0
	ut2(0) = 0.d0
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
	ut2(ndate)=som2+the2(i-4)+the2(i-3)+the2(i-2)+the2(i-1)
	dut1(ndate) = (4.d0*the1(i-1)/h1)
	dut2(ndate) = (4.d0*the2(i-1)/h1)

	
!-------------------------------------------------------
!--------- calcul de la vraisemblance ------------------
!--------------------------------------------------------

!--- avec ou sans variable explicative  ------cc
	res1= 0.d0
	res2= 0.d0
	res3= 0.d0
	cpt = 0

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
				res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet)
			endif  
	
			if((c(i).eq.1).and.(stra(i).eq.2))then
				res2(g(i)) = res2(g(i))+dlog(dut2(nt1(i))*vet)
			endif
	               if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                          funcpassplines=-1.d9
                          goto 123
                       end if	
			if(stra(i).eq.1)then
				res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet-ut1(nt0(i))*vet 
				RisqCumul(i) = ut1(nt1(i))*vet
			endif
	
			if(stra(i).eq.2)then
				res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet-ut2(nt0(i))*vet 
				RisqCumul(i) = ut2(nt1(i))*vet
			endif
	               if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                          funcpassplines=-1.d9
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
		                  funcpassplines=-1.d9
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
				res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet)
			endif  
			if((c(i).eq.1).and.(stra(i).eq.2))then
				res2(g(i)) = res2(g(i))+dlog(dut2(nt1(i))*vet)
			endif  
	               if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                          funcpassplines=-1.d9
                          goto 123
                       end if
			if(stra(i).eq.1)then
				res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet
			endif
			
			if(stra(i).eq.2)then
				res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet
			endif
	               if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                          funcpassplines=-1.d9
                          goto 123
                       end if
! modification pour nouvelle vraisemblance / troncature:
			if(stra(i).eq.1)then
				res3(g(i)) = res3(g(i)) + ut1(nt0(i))*vet
			endif
			
			if(stra(i).eq.2)then
				res3(g(i)) = res3(g(i)) + ut2(nt0(i))*vet
			endif
	               if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                          funcpassplines=-1.d9
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
				       if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
				          funcpassplines=-1.d9
				          goto 123
				       end if
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
				       if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
				          funcpassplines=-1.d9
				          goto 123
				       end if
				endif 
			endif 
		end do
	
       	endif !fin boucle effet=0

!--------- calcul de la penalisation -------------------

	pe1 = 0.d0
	pe2 = 0.d0
	pe=0.d0
	do i=1,n-3
	
		pe1 = pe1+(the1(i-3)*the1(i-3)*m3m3(i))+(the1(i-2) &
		*the1(i-2)*m2m2(i))+(the1(i-1)*the1(i-1)*m1m1(i))+( &
		the1(i)*the1(i)*mmm(i))+(2.d0*the1(i-3)*the1(i-2)* &
		m3m2(i))+(2.d0*the1(i-3)*the1(i-1)*m3m1(i))+(2.d0* &
		the1(i-3)*the1(i)*m3m(i))+(2.d0*the1(i-2)*the1(i-1)* &
		m2m1(i))+(2.d0*the1(i-2)*the1(i)*m2m(i))+(2.d0*the1(i-1) &
		*the1(i)*m1m(i))
		pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2) &
		*the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+( &
		the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)* &
		m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0* &
		the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)* &
		m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1) &
		*the2(i)*m1m(i))
	
	end do
    

!    Changed JRG 25 May 05
	if (nst.eq.1) then
		pe2=0.d0
	end if

	pe = k0(1)*pe1 + k0(2)*pe2 
	
	resnonpen = res

	res = res - pe

	if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
		funcpassplines=-1.d9
		goto 123
	end if

	funcpassplines = res 

	do k=1,ng
		cumulhaz(k)=res1(k)
	end do

123     continue

	return

	end function funcpassplines


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

	h = (zi(n)-zi(1))*0.01d0
	x1 = zi(1)
	x2 = zi(1)     
	
	do i=1,mt 
		
		if(i .ne.1)then
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

