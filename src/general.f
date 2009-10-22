
c=================================  DMFSD  ===================================

      SUBROUTINE DMFSD(A,N,EPS,IER)
C
C   FACTORISATION DE CHOLESKY D'UNE MATRICE SDP
C   MATRICE = TRANSPOSEE(T)*T
C   ENTREE : TABLEAU A CONTENANT LA PARTIE SUPERIEURE STOCKEE COLONNE
C            PAR COLONNE DE LA METRICE A FACTORISER
C   SORTIE : A CONTIENT LA PARTIE SUPPERIEURE DE LA MATRICE SYMETRIQUE T
C 
C   SUBROUTINE APPELE PAR DSINV
C  
C   N : DIM. MATRICE 
C   EPS : SEUIL DE TOLERANCE
C   IER = 0 PAS D'ERREUR
C   IER = -1 ERREUR
C   IER = K COMPRIS ENTRE 1 ET N, WARNING, LE CALCUL CONTINUE
C
      DOUBLE PRECISION A(N*(N+1)/2)
      DOUBLE PRECISION DPIV,DSUM
      DOUBLE PRECISION EPS,TOL
      INTEGER I,K,L,N,IER,KPIV,IND,LEND,LANF,LIND
C
C   TEST ON WRONG INPUT PARAMETER N
C
      IF(N-1) 12,1,1
1     IER=0
C
C   INITIALIZE DIAGONAL-LOOP
C
      KPIV=0
      DO 11 K=1,N
      KPIV=KPIV+K
      IND=KPIV
      LEND=K-1
C
C   CALCULATE TOLERANCE
C
      TOL=DABS(EPS*SNGL(A(KPIV)))
C
C   START FACTORIZATION-LOOP OVER K-TH ROW
C
      DO 11 I=K,N
      DSUM=0.D0
      IF(LEND) 2,4,2
C
C   START INNER LOOP
C
2     DO 3 L=1,LEND
      LANF=KPIV-L
      LIND=IND-L
3     DSUM=DSUM+A(LANF)*A(LIND)
C
C   END OF INNEF LOOP
C
C   TRANSFORM ELEMENT A(IND)
C
4     DSUM=A(IND)-DSUM
      IF(I-K)10,5,10
C
C   TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE
C
5     IF(SNGL(DSUM)-TOL)6,6,9
6     IF(DSUM)12,12,7
7     IF(IER)8,8,9
8     IER=K-1
C
C   COMPUTE PIVOT ELEMENT
C
9     DPIV=DSQRT(DSUM)
      A(KPIV)=DPIV
      DPIV=1.D0/DPIV
      GO TO 11
C
C   CALCULATE TERMS IN ROW
C
10    A(IND)=DSUM*DPIV
11    IND=IND+I
C
C   END OF DIAGONAL-LOOP
C
      RETURN
12    IER=-1
      RETURN
 
      END

c==============================   DSINV  ======================================

      SUBROUTINE DSINV(A,N,EPS,IER)
C
C     INVERSION D'UNE MATRICE SYMETRIQUE DEFINIE POSITIVE :
C
C     MATRICE = TRANSPOSEE(T)*T
C     INERSE(MATRICE) = INVERSE(T)*INVERSE(TRANSPOSEE(T))
C
C     A : TABLEAU CONTENANT LA PARTIE SUPERIEURE DE LA MATRICE A INVERSER
C         STOCKEE COLONNE PAR COLONNE
C     DIM. MATRICE A INVERSER = N 
C     DIM. TABLEAU A = N*(N+1)/2
C
C     EPS : SEUIL DE TOLERANCE AU-DESSOUS DUQUEL UN PIVOT EST CONSIDERE
C           COMME NUL
C
C     IER : CODE D'ERREUR
C         IER=0 PAS D'ERREUR
C         IER=-1 ERREUR SUR LA DIM.N OU MATRICE PAS DEFINIE POSITIVE
C         IER=1 PERTE DE SIGNIFICANCE, LE CALCUL CONTINUE
C
      DOUBLE PRECISION A(N*(N+1)/2)
      DOUBLE PRECISION DIN,WORK
      DOUBLE PRECISION EPS
      INTEGER N,IER,IND,IPIV,I,J,K,L,MIN,KEND
      INTEGER LHOR,LVER,LANF
C
C     FACTORIZE GIVEN MATRIX BY MEANS OF SUBROUTINE DMFSD
C     A=TRANSPOSE(T) * T
C
      CALL DMFSD(A,N,EPS,IER)
      IF(IER) 9,1,1
C
C     INVERT UPPER TRIANGULAR MATRIX T
C     PREPARE INVERSION-LOOP
C
1     IPIV=N*(N+1)/2
      IND=IPIV
C
C     INITIALIZE INVERSION-LOOP
C
      DO 6 I=1,N
      DIN=1.D0/A(IPIV)
      A(IPIV)=DIN
      MIN=N
      KEND=I-1
      LANF=N-KEND
      IF(KEND) 5,5,2
2     J=IND
C
C     INITIALIZE ROW-LOOP
C
      DO 4 K=1,KEND
      WORK=0.D0
      MIN=MIN-1
      LHOR=IPIV
      LVER=J
C
C     START INNER LOOP
C
      DO 3 L=LANF,MIN
      LVER=LVER+1
      LHOR=LHOR+L
3     WORK=WORK+A(LVER)*A(LHOR)
C
C     END OF INNER LOOP
C
      A(J)=-WORK*DIN
4     J=J-MIN
C
C     END OF ROW-LOOP
C
5     IPIV=IPIV-MIN
6     IND=IND-1
C
C     END OF INVERSION-LOOP
C
C     CALCULATE INVERSE(A) BY MEANS OF INVERSE(T)
C     INVERSE(A) = INVERSE(T) * TRANSPOSE(INVERSE(T))
C     INITIALIZE MULTIPLICATION-LOOP
C
      DO 8 I=1,N
      IPIV=IPIV+I
      J=IPIV
C
C     INITIALIZE ROW-LOOP
C
      DO 8 K=I,N
      WORK=0.D0
      LHOR=J
C
C     START INNER LOOP
C
      DO 7 L=K,N
      LVER=LHOR+K-I
      WORK=WORK+A(LHOR)*A(LVER)
7     LHOR=LHOR+L
C
C     END OF INNER LOOP
C
      A(J)=WORK
8     J=J+K
C
C     END OF ROW-AND MULTIPLICATION-LOOP
C
9     RETURN
 
      END

c===============================    MAXT    =============================

      DOUBLE PRECISION FUNCTION MAXT(DELTA,M)
      DOUBLE PRECISION DELTA(m)
      INTEGER  I,M
c
      MAXT=DELTA(1)
      DO 2 I=2,M
      IF(DELTA(I).GT.MAXT) MAXT=DELTA(I)
2     CONTINUE
      RETURN
      END

c================================  DCHOLE  ===========================
      subroutine dchole(a,k,nq,idpos)

      parameter(npmax=50)      
      integer  k,nq,i,ii,i1,i2,i3,m,is,j,k2,jmk
      integer  ijm,irm,jji,jjj,l,jj,iil,jjl,il,idpos
      double precision a((npmax*(npmax+3)/2)) 
      double precision  term,xn,diag,p
      dimension is(500)
      equivalence (term,xn)
c
c      ss programme de resolution d'un systeme lineaire symetrique
c
c       k ordre du systeme /
c       nq nombre de seconds membres
c
c       en sortie les seconds membres sont remplaces par les solutions
c       correspondantes
c
      idpos=0
      k2=k+nq
c     calcul des elements de la matrice
      do 13 i=1,k
      ii=i*(i+1)/2
c     elements diagonaux
      diag=a(ii)
      i1=ii-i
      if(i-1) 1,4,1
1     i2=i-1
      do 3 l=1,i2
      m=i1+l
      p=a(m)
      p=p*p
      if(is(l)) 2,3,3
2     p=-p
3     diag=diag-p
4     if(diag) 5,50,6
5     is(i)=-1
      idpos=idpos+1
      diag=-dsqrt(-diag)
      a(ii)=-diag
      go to 7
6     is(i)=1
      diag=dsqrt(diag)
      a(ii)=diag
c       elements non diagonaux
7     i3=i+1
      do 13 j=i3,k2
      jj=j*(j-1)/2+i
      jmk=j-k-1
      if(jmk) 9,9,8
8     jj=jj-jmk*(jmk+1)/2
9     term=a(jj)
      if(i-1) 10,13,10
10    do 12 l=1,i2
      iil=ii-l
      jjl=jj-l
      p=a(iil)*a(jjl)
      il=i-l
      if(is(il)) 11,12,12
11    p=-p
12    term=term-p
13    a(jj)=term/diag
c       calcul des solutions
      jj=ii-k+1
      do 45 l=1,nq
      jj=jj+k
      i=k-1
14    jji=jj+i
      xn=a(jji)
      if(i-k+1) 20,22,22
20    j=k-1
21    jjj=jj+j
      ijm=i+1+j*(j+1)/2
      xn=xn-a(jjj)*a(ijm)
      if(j-i-1) 22,22,30
30    j=j-1
      go to 21
22    irm=(i+1)*(i+2)/2
      a(jji)=xn/a(irm)
      if(i) 45,45,40
40    i=i-1
      go to 14
45    continue
50    continue
      return 
      end


c================== multiplication de matrice  ==================

c multiplie A par B avec le resultat dans C

      subroutine multi(A,B,IrowA,JcolA,JcolB,C)
c     remarque :  jcolA=IrowB

      integer , parameter :: npmax=50
      
      integer IrowA,JcolA,JcolB,i,j,k
      double precision sum
      double precision A(npmax,npmax) ,B(npmax,npmax) ,C(npmax,npmax) 
      
      do 1 I=1,IrowA
         do 2 J=1,JcolB
            sum=0
            do 3 K=1,JcolA
               sum=sum+A(I,K)*B(K,J)
 3          continue
            C(I,J)=sum
 2       continue
 1    continue
      return
      end
c====================================================================



c======================  LUBKSB  ======================================
      subroutine lubksb(a,n,np,indx,b)
      parameter(npmax=50)
         integer n,np,indx(npmax)
         double precision a(npmax,npmax),b(npmax)
         integer i,ii,j,ll
         double precision sum

         ii = 0
         do 12 i=1,n
            ll = indx(i)
            sum = b(ll)
            b(ll) = b(i)
            if(ii.ne.0)then
               do 11 j=ii,i-1
                  sum = sum -a(i,j)*b(j)
 11            continue
            else
               if(sum.ne.0.d0)then
                  ii=i
               endif
            endif
            b(i)=sum
 12      continue
         do 14 i=n,1,-1
            sum = b(i)
            do 13 j = i+1,n
               sum = sum-a(i,j)*b(j)
 13         continue
            b(i)=sum/a(i,i)
 14      continue
         return

         end
c==================
c======================  LUDCMP  ======================================
       subroutine ludcmp(a,n,np,indx,d)
       parameter(npmax=50)

         integer n,np,indx(n),nmax
         double precision d,a(npmax,npmax),tiny
         parameter (nmax=500,tiny=1.d-20)
         integer i,imax,j,k
         double precision aamax,dum,sum,vv(nmax)

         d = 1.d0
         do 12 i=1,n
            aamax=0.d0
            do 11 j=1,n
               if (dabs(a(i,j)).gt.aamax)then
                  aamax=dabs(a(i,j))
               endif
 11         continue
c          if (aamax.eq.0.d0) 
c               pause 'matrice singuliere'
            vv(i) = 1.d0/aamax
 12      continue
         do 19 j = 1,n
            do 14 i=1,j-1
               sum = a(i,j)
               do 13 k=1,i-1
                  sum = sum - a(i,k)*a(k,j)
 13            continue
               a(i,j) = sum
 14         continue
            aamax = 0.d0
            do 16 i = j,n
               sum = a(i,j)
               do 15 k=1,j-1
                  sum = sum -a(i,k)*a(k,j)
 15            continue
               a(i,j) = sum
               dum = vv(i)*dabs(sum)
               if (dum.ge.aamax) then
                  imax = i
                  aamax = dum
               endif
 16         continue
            if(j.ne.imax)then
               do 17 k=1,n
                  dum = a(imax,k)
                  a(imax,k)=a(j,k)
                  a(j,k) = dum
 17            continue
               d = -d
               vv(imax)=vv(j)
            endif
            indx(j)=imax
            if(a(j,j).eq.0.d0)then
               a(j,j)=tiny
            endif
            if(j.ne.n)then
               dum = 1.d0/a(j,j)
               do 18 i = j+1,n
                  a(i,j) = a(i,j)*dum
 18            continue
            endif
 19      continue
         return
         end


