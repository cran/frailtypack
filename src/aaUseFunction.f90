	
! ******************** BGOS ********************************

	subroutine bgos(sx,id,x1,x2,ro)

! ID=1:U(0,SX); ID DIFF DE 1 :N(0,SX)
	implicit none
	
	integer::id
	double precision::ro,sx,f,v1,v2,s,dls,ro2,x1,x2,uniran

5     continue
      
	x1=uniran()
	x2=uniran()
    
	if (id.ne.1) goto 10
	f=2.*sqrt(3.)
	x1=(x1-0.5)*f
	x2=(x2-0.5)*f
	goto 20
10    continue
	v1=2.*x1-1
	v2=2.*x2-1
	s=v1*v1+v2*v2
	if (s.ge.1.) goto 5
	dls=sqrt(-2.*log(s)/s)
	x1=v1*dls
	x2=v2*dls
20    continue
	ro2=ro*ro
	if (abs(ro).gt.1.E-10) x2=(x1+x2*sqrt(1./ro2-1.))*ro
	x1=x1*sx
	x2=x2*sx
	
	return
	
	end subroutine bgos


!------------------- FIN SUBROUTINE BGOS -----------------

! ------------------------------------------------------

	double precision function uniran()
!
!     Random number generator(RCARRY), adapted from F. James
!     "A Review of Random Number Generators"
!      Comp. Phys. Comm. 60(1990), pp. 329-344.
!
	double precision,save::carry
	double precision,dimension(24),save::seeds 
	double precision,parameter::one=1 
	double precision,parameter::twom24 = ONE/16777216
	integer,save::i,j
	data i, j, carry / 24, 10, 0.0 /
	data seeds / &
	0.8804418, 0.2694365, 0.0367681, 0.4068699, 0.4554052, 0.2880635, &
	0.1463408, 0.2390333, 0.6407298, 0.1755283, 0.7132940, 0.4913043, &
	0.2979918, 0.1396858, 0.3589528, 0.5254809, 0.9857749, 0.4612127, &
	0.2196441, 0.7848351, 0.4096100, 0.9807353, 0.2689915, 0.5140357/
      
	uniran = seeds(i) - seeds(j) - carry
	
	if (uniran .lt. 0) then
		uniran = uniran + 1
		carry = twom24
	else
		carry = 0
	end if
	
	seeds(I) = uniran
	I = 24 - MOD( 25-I, 24 )
	J = 24 - MOD( 25-J, 24 )
	
	end function uniran	

	
	
      
	subroutine percentile(t,t25,t975) ! pour le MCMC
	    
	integer::indd,i    
	double precision::t25,t975,temp
	double precision,dimension(1000)::t
      
	indd=1
	do while (indd.eq.1)
		indd=0
		do i=1,999
			if (t(i).gt.t(i+1)) then
				temp=t(i)
				t(i)=t(i+1)
				t(i+1)=temp
				indd=1        
			end if
		end do  
	end do 
		t25=t(25)
		t975=t(975) 
	
	end subroutine percentile
	
