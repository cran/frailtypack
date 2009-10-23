c avec crossvalidation


c programme modifié le 21 mars 2005 
c on a la possibilité de choisir deux formaulations pour la vraisemblance :
c soit la vraisemblance avec troncature a gauche : cf LIDA2003
c soit l'ancienne version de la troncature, adaptée aux données récurrentes de type 
c andersen gill (calendar time????? --> counting process formulation), qui permet aussi 
c de traiter des variables explicatievs dépendantes du temps .


c********************************************************************
c********  version F77 pour le passer sous R ************************
c********************************************************************
c avec ou sans stratification
c programme repris en janvier 2003 pour le comparer 
c a frailtycarre.f qui utilise l algo numerique de DC
c  sur la vraisemblance partielle marginale
c on supprime dans marquard le chgt de variable ,
c juste le considerer au mment de l affichage de se(theta)=2*b(np-nva)*sqrt(v(1))
c Rem: en modifiant le frailty_nonstrat.inf on retombe sur le programme frailty.f
c  avec stratification
c********************************************************************
c ancien nouveau programme ...
c********************************************************************
c bis le 7/11/2001, modification de funcpa, erreur sur ddl d ordre 3
c********************************************************************
c programme modifie le 01juin2001, car modification de la vraisemblance
c en tenant compte de la loi de la frailty conditionnelle a la survie

c********************************************************************
c  puis, modification de marquard (erreur) le 8juin2001
c     car dans dchol le vecteur v() est modifie, boucle 4001 rajoutee 
c  lorsque la matrice des derivees 2 n'est pas inversible,
c  il faut reinitialiser completement la matrice FU et
c  pas seulement modifier la diagonale car FU a ete modifie dans Dchole.
c********************************************************************

c BANDES DE CONFIANCES 
c bayesiennes (H-1) ou non baysienne H-1.I.H-1

c TRACE = nb de degres de liberte pour les fonctions de risque
c  a partir du kappa fixe on peut en deduire le mdf

c BIAIS
c on essaye d obtenir un estimateur sans biais pour theta,
c en corrigeant par -H-1(dPn/dtheta) = biais
c mais , le terme de penalisation ne depend pas de theta !!!
c faire correction uniquement sur les parametres des splines

c VARIANCE des parametres :
c en utisant la matrice sandwich H-1 I H-1
c ici dans marquard, on modifie la matrice hessienne 
c (en tenant compte du changement de variable) avant de l'inverser 
c donc apres marquard on a directement la variance par H_hess.I_hess.H_hess
c TEST
c de Wald en utilisant la variance corrigee ou non 
c avec test du LRS penalise ?

c dans chaque simulation, on n'initialise pas les parametres avec valeurs precedentes

CCCC vraisemblance calculee avec ou sans effet aleatoire

CCCC programme qui permet des donnees par strates (2 uniquement)
CCCC la fonction de risque de base est donc différente
CCCC  pour chaque strate (2 fonctions differentes)
CCCC  estimees par 2 bases de splines differentes

c I_hess : -hessienne non inversee sur vraisemblance non penalisee
c H_hess : inverse de -hessienne  sur vraisemblance penalisee

c  np:     total number of parameters

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c    estimation des coefficient de regression  ... spline d'ordre 4
c    avec censure a droite, troncature a gauche, noeuds aux choix


c                        avec fraitly avec loi Gamma

      subroutine frailpenal(nsujetAux,ngAux,icenAux,nstAux,effetAux,
     &             nzAux,ax1,ax2,tt0Aux,tt1Aux,icAux,groupeAux,
     &             nvaAux,strAux,vaxAux,AGAux,noVar,maxitAux,irep1,
     &             np,b,H_hessOut,HIHOut,resOut,x1Out,lamOut,suOut,
     &             x2Out,lam2Out,su2Out,ni,cpt,ier,k0,ddl)

c
c  Obs: noVar=0 indicates no variables but I need to pass all 0's in vaxAux

       parameter(npmax=50,nsujetmax=15000,nvarmax=50,ngmax=1000)
       parameter(nboumax=1000,ndatemax=30000)


       integer  groupe,ij,kk,j,k,nz,n,np,cpt,ii,iii,ver
       integer  cptstr1,cptstr2,trace,trace1,trace2
       integer  i,ic,ni,ier,istop,ef
       integer  cptni,cptni1,cptni2,nb_echec,nb_echecor
       integer  ibou,nbou2,id,cptbiais,l
       integer  m,nbou,idum,icen
       integer  filtre(nvarmax) 
       integer  cpt1(nboumax) 
       integer  cpt2(nboumax) 
       integer  cpt3(nboumax) 
       integer  ind(nboumax)
       
       integer  stracross(nsujetmax)
       integer irep1,nvacross,nstcross,effetcross

       real vax(nvarmax)
       double precision tt0
       double precision tt1
      
       double precision h
       double precision ro,wres,csi,csi1
       double precision ax1,ax2,res,min,max,maxt
       double precision bi,bs,wald
       double precision moyvar,moyse,moyse_cor,str
       double precision varxij,eca,varsmarg,smoy,smoyxij
       double precision pe1,pe2,moy_peh0,moy_peh1,lrs


       double precision auxkappa(2)
       double precision aux(2*nsujetmax)
       double precision res_tab(ngmax) 
       double precision v((npmax*(npmax+3)/2))
       double precision k0(2),res01(2)
       double precision b(npmax)
       double precision se(nboumax),se_cor(nboumax),tvars(nboumax)
       double precision biais_theta(nboumax)
       double precision I1_hess(npmax,npmax),H1_hess(npmax,npmax)
       double precision I2_hess(npmax,npmax),H2_hess(npmax,npmax)
       double precision HI1,HI2(npmax,npmax)
       double precision HIH(npmax,npmax),IH(npmax,npmax),HI(npmax,npmax)
 





c******************************************   Add JRG January 05

         double precision  x1,x2,lam,lbinf,lbsup,margi
         double precision  lam2,lbinf2,lbsup2,margi2
         integer nsujetAux,ngAux,icenAux,nstAux,effetAux,nzAux,nvaAux 
         double precision tt0Aux(nsujetAux),tt1Aux(nsujetAux)
         integer icAux(nsujetAux),groupeAux(nsujetAux),ss,sss
         double precision strAux(nsujetAux)
         double precision vaxAux(nsujetAux,nvaAux)
         double precision resOut,H_hessOut(npmax,npmax)
         double precision HIHOut(npmax,npmax)


         double precision  x1Out(99),lamOut(99,3),suOut(99,3),
     &                     x2Out(99),lam2Out(99,3),su2Out(99,3)

         integer noVar,AGAux,maxitAux

c*******************************************  Add JRG May 05 (Cross-validation)

       double precision auxi,ax,bx,cx,tol,ddl
       double precision fa,fb,fc,golden,estimv
       double precision y(npmax,npmax)  
       double precision xmin1,xmin2




c*****************************************************************
c*****dace1 
      double precision date(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,zi
c*****dace2
      double precision t0(nsujetmax),t1(nsujetmax)
      double precision  ncsrl
      integer c(nsujetmax)
      integer nt0(nsujetmax),nt1(nsujetmax)
      integer  nsujet,nsim,nva,ndate,nst,maxit
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,
     &              nst,maxit
c*****dace4
      integer  stra(nsujetmax)
      common /dace4/stra
c*****ve1
      double precision ve(nsujetmax,nvarmax)
      common /ve1/ve
c*****dace3
      double precision  pe
      integer  effet,nz1,nz2
      common /dace3/pe,effet,nz1,nz2
c*****dace7
      double precision I_hess(npmax,npmax),H_hess(npmax,npmax)
      double precision Hspl_hess(npmax,npmax)
      double precision PEN_deri(npmax,1)
      double precision hess(npmax,npmax)
      common /dace7/PEN_deri,I_hess,H_hess,Hspl_hess,hess
c*****contrib
      integer ng          !nb de gpes
      common /contrib/ng       
c*****groupe
      integer g(nsujetmax)
      integer nig(ngmax) 
      common /gpe/g,nig

c*****mem1
      double precision mm3(ndatemax),mm2(ndatemax)
      double precision mm1(ndatemax),mm(ndatemax)
      common /mem1/mm3,mm2,mm1,mm
c %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
      integer AG
      common /andersengill/AG



c************ FIN COMMON ***********************************
     
c nst: nombre de strates
c ist: appartenance aux strates
c icen=1: censure
c ib: matrice de variance de beta chapeau
c I_hess : -hessienne non inversee sur vraisemblance non penalisee
c H_hess : inverse de -hessienne  sur vraisemblance penalisee
        
    
      id=1

      cptni=0
      cptni1=0
      cptni2=0
      biais_moy=0.d0
      cptbiais=0

       

            

c
c  --------- Added January 05 JRG  
c
 
      ncsrl=0.d0
      pe=0.d0

      do i=1,nsujetmax
        c(i)=0
        nt0(i)=0
        nt1(i)=0
        g(i)=0
        do j=1,nvarmax 
           ve(i,j)=0.d0
        end do 
      enddo

      
      nsujet=0
      nsim=0
      ndate=0
      nst=0    

      do i=1,ngmax
        nig(i)=0
      enddo 
       
      do i=-2,npmax
        zi(i)=0.d0
      enddo  

      do i=1,ndatemax
        date(i)=0.d0
        mm(i)=0.d0
        mm1(i)=0.d0
        mm2(i)=0.d0
        mm3(i)=0.d0
      end do   

      ibou=1  
      ind(ibou)=0
      ij=0
      kk=0

      np=0 
      ni=0
      cpt=0  
      do i=1,npmax
        b(i)=0.d0
      end do  
      resOut=0.d0     

      nsujet=nsujetAux
      ng=ngAux
      icen=icenAux
      nst=nstAux
      effet=effetAux 
      nz=nzAux

      if (noVar.eq.1) then 
        do i=1,nvaAux
         filtre(i)=0
        enddo  
        nva=0  
      else
        do i=1,nvaAux
         filtre(i)=1
        enddo  
        nva=nvaAux  
      end if  

      ver=nvaAux


      res=0.d0
      
      do i=1,npmax*(npmax+3)/2
        v(i)=0.d0
      end do 

      do i=1,npmax
        Pen_deri(i,1)=0.d0
        do j=1,npmax
          H_hess(i,j)=0.d0
          I_hess(i,j)=0.d0
          Hspl_hess(i,j)=0.d0
          hess(i,j)=0.d0
        end do
      end do
      
      AG=AGAux
      maxit=maxitAux 

      istop=0
c      ier=0     

       

        


C**************************************************
C**************************************************
C********************* prog spline****************


        res01(1)=0.d0
        res01(2)=0.d0
      
           cpt1(ibou)=0
           cpt2(ibou)=0 
           cpt3(ibou)=0 


       
       
         
c------------  lecture fichier -----------------------

         maxt = 0.d0
         
         cpt = 0
         k = 0
         cptstr1 = 0
         cptstr2 = 0

         do 10 i = 1,nsujet 
            str=1
            if(nst.eq.2)then                     
               tt0=tt0Aux(i)
               tt1=tt1Aux(i)
               ic=icAux(i)
               groupe=groupeAux(i)
               str=strAux(i)
               do j=1,nva
                 vax(j)=vaxAux(i,j)  
               enddo
            else 
               tt0=tt0Aux(i)
               tt1=tt1Aux(i)
               ic=icAux(i)
               groupe=groupeAux(i)
               do j=1,nva
                 vax(j)=vaxAux(i,j)  
               enddo
            endif
            k = k +1


c     essai sans troncature
c            tt0=0.
c------------------   observation c=1
               if(ic.eq.1)then
                  cpt = cpt + 1
                  c(k)=1
                  if(str.eq.1)then
                  stra(k) = 1
                  cptstr1 = cptstr1 + 1
                  else
                     if(str.eq.2)then
                        stra(k) = 2
                        cptstr2 = cptstr2 + 1
                     endif
                  endif
                  t0(k) = tt0
                  t1(k) = tt1
                  g(k) = groupe

c ! nb de dc dans un groupe
                  nig(groupe) = nig(groupe)+1 


                  iii = 0
                  do 6 ii = 1,ver
                     if(filtre(ii).eq.1)then
                        iii = iii + 1
                        ve(i,iii) = dble(vax(ii))
c                        write(7,*)'** ve **', ve(i,iii),i,iii
                     endif
 6                continue   


                

               else 
c------------------   censure a droite  c=0
                  if(ic.eq.0)then
                     c(k) = 0 
                     if(str.eq.1)then
                        stra(k) = 1
                        cptstr1 = cptstr1 + 1
                     else
                        if(str.eq.2)then
                           stra(k) = 2
                           cptstr2 = cptstr2 + 1
                        endif
                     endif
                     iii = 0
                         

                     do 8 ii = 1,ver
                        if(filtre(ii).eq.1)then
                           iii = iii + 1
                           ve(i,iii) = dble(vax(ii))
                        endif
 8                    continue 
                      
              
                      t0(k) =  tt0
                      t1(k) = tt1
                      g(k) = groupe
                   endif
                endif
                if (maxt.lt.t1(k))then
                   maxt = t1(k)
                endif



             
 10         continue 

             

c %%%%%%%%%%%%% SANS EFFET ALEATOIRE %%%%%%%%%%%%%%%%%%%%%%%%% 


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

         nsujet = k
                
            nz1=nz
            nz2=nz
         if(nz.gt.20)then
            nz = 20
         endif
         if(nz.lt.4)then
            nz = 4
         endif



c***************************************************
c--------------- zi- ----------------------------------

c      construire vecteur zi (des noeuds)

         min = 1.d-10
c          min = 0.d0
          
 
         max = maxt
         

         do 15 i = 1,2*nsujet
            do 16 k = 1,nsujet
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
 16            continue   
            aux(i) = max
            min = max + 1.d-12
            max = maxt
 15      continue

         date(1) = aux(1)
         k = 1
         do 17 i=2,2*nsujet
               if(aux(i).gt.aux(i-1))then
                  k = k+1
                  date(k) = aux(i)
               endif 
 17         continue 
         ndate = k

         
         zi(-2) = date(1)
         zi(-1) = date(1)
         zi(0) = date(1)
         zi(1) = date(1)
         h = (date(ndate)-date(1))/dble(nz-1)

         do 18 i=2,nz-1
            zi(i) =zi(i-1) + h
 18      continue
         
         zi(nz) = date(ndate)
         zi(nz+1)=zi(nz)
         zi(nz+2)=zi(nz)
         zi(nz+3)=zi(nz)

 
c---------- affectation nt0,nt1----------------------------

            do 50 i=1,nsujet 
               if(t0(i).eq.0.d0)then
                  nt0(i) = 0
               endif
               do 45 j=1,ndate
                  if(date(j).eq.t0(i))then
                     nt0(i)=j
                  endif
                  if(date(j).eq.t1(i))then
                     nt1(i)=j
                  endif
 45            continue
 50          continue 

c         write(*,*)'nt0',(nt0(i),i=1,nsujet)
c         write(*,*)'nt1',(nt1(i),i=1,nsujet)  
c         stop

c---------- affectation des vecteurs de splines -----------------
           
             n  = nz+2
           
            
             call vecspli(n,ndate) 
             call vecpen(n)
             
             np = nst*n + nva + effet


c          write(*,*)'nombre total de paramètres',np,nst,n,nva,effet,nz
               
               
c------- initialisation des parametres
                   
                   do 75 i=1,np
                      b(i)=1.d-1
 75                continue


          
C deux parametres de lissage :

c         k0(1) = ax1
c         if(nst.eq.2)then
c            k0(2) = ax2
c         else
c            k0(2) = 0.d0 
c         endif
          


c    Esto se cambia ya que ahora xmin1 es kappa1
 
      
         xmin1=ax1 
         if(nst.eq.2)then
            xmin2 = ax2
         else
            xmin2 = 0.d0 
         endif


c***********************************************************
c************** NEW : cross validation  ***********************
ccccc sur une seule strate, sans var expli , sans frailties ****
c***********************************************************
         
         nvacross=nva !pour la recherche du parametre de lissage sans var expli
         nva=0
         effetcross=effet
         effet=0
         nstcross=nst
         nst=1

         do 751 l=1,nsujet  
            stracross(l)=stra(l)
 751     continue
         do 752 l=1,nsujet  
            stra(l)=1
 752     continue
        

         if(irep1.eq.1)then   !pas recherche du parametre de lissage

            xmin1 = dsqrt(xmin1)
c            write(*,*) 'auxi',n,nst,effet,effetcross,np,K0,xmin1
c            write(*,*)'==avant premier estimv',n,k0,xmin1,effet,
c     & effetcross,nst,nvacross
c            stop

            auxi = estimv(xmin1,n,b,y,ddl,ni,res)

            if (ni.ge.maxit) then
c              write(*,*) ' '
c              write(*,*) 'no convergence 
c     &             with the chosen smoothing parameter'

              do 175 i=1,nz+2
                 b(i)=1.d-1
 175           continue     
              xmin1 = sqrt(10.d0)*xmin1
              auxi = estimv(xmin1,n,b,y,ddl,ni,res)
              if (ni.lt.maxit) then
c                 write(*,*)' '
c            write(*,*)'Value of the smoothing parameter :'
c     &       ,real(xmin1*xmin1), '  DoF :',-ddl
c                 write(*,*)' '
c                 write(*,*)'Log-vraisemblance :',res
c                 stop
              else
                 do 176 i=1,nz+2
                    b(i)=1.d-1
 176             continue     
                 xmin1 = sqrt(10.d0)*xmin1
                 auxi = estimv(xmin1,n,b,y,ddl,ni,res)
                 if (ni.lt.maxit) then
c                    write(*,*)' '
c            write(*,*)'Value of the smoothing parameter :'
c     &       ,real(xmin1*xmin1), '  DoF :',-ddl
c                    write(*,*)' '
c                    write(*,*)'Log-vraisemblance :',res
                 endif   
               endif
            else
c               write(*,*)' '
c            write(*,*)'Value of the smoothing parameter :'
c     &       ,real(xmin1*xmin1), '  DoF :',-ddl
c               write(4,*)' '
c            write(4,*)'Value of the smoothing parameter :'
c     &       ,real(xmin1*xmin1), '  DoF :',-ddl
c               write(*,*)' '
c               write(*,*)'Log-vraisemblance :',res
            endif

            
c----------------------------------------------------
         else                   !recherche du parametre de lissage
            
c            write(*,*)' b ',b
c            stop
            if(xmin1.le.0.d0)then
               xmin1 = 1.d0
            endif  
c            write(*,*)' '
c            write(*,*)'                Searching smoothing parameter'
c            write(*,*)' '
            xmin1 = dsqrt(xmin1)
            auxi = estimv(xmin1,n,b,y,ddl,ni,res)
c            write(*,*)'estimv1',ddl
c               stop
            if(ddl.gt.-2.5d0)then
c            write(*,*)'OK non'
c               stop
               xmin1 = dsqrt(xmin1)
c               stop
               auxi = estimv(xmin1,n,b,y,ddl,ni,res)
c               write(*,*)'estimv2'
c               stop !pb
               if(ddl.gt.-2.5d0)then
                  xmin1 = dsqrt(xmin1)
                  auxi = estimv(xmin1,n,b,y,ddl,ni,res)
c                  write(*,*)'estimv3'
                  if(ddl.gt.-2.5d0)then
                     xmin1 = dsqrt(xmin1)
                     auxi = estimv(xmin1,n,b,y,ddl,ni,res)
c                     write(*,*)'estimv4'
                     if(ddl.gt.-2.5d0)then
                        xmin1 = dsqrt(xmin1)
                        auxi = estimv(xmin1,n,b,y,ddl,ni,res)
c                        write(*,*)'estimv5'
                        if(ddl.gt.-2.5d0)then
                           xmin1 = dsqrt(xmin1)
                        endif   
                     endif   
                  endif   
               endif
            endif 
c            write(*,*)'estimv1 bis',ddl
            if (ni.ge.maxit) then
              do 275 i=1,nz+2
                 b(i)=1.d-1
 275          continue     
              xmin1 = sqrt(10.d0)*xmin1
              auxi = estimv(xmin1,n,b,y,ddl,ni,res)
c              write(*,*)'estimv5'
              if (ni.ge.maxit) then
                 do 276 i=1,nz+2
                    b(i)=1.d-1
 276             continue     
                 xmin1 = sqrt(10.d0)*xmin1
              endif
            endif 
            ax = xmin1
            bx = xmin1*dsqrt(1.5d0)  
c            write(*,*)'avant mnbrak' ,ax ,bx,cx,fa,fb,fc,n
            
            call mnbrak(ax,bx,cx,fa,fb,fc,b,n)
c            write(*,*)'OK 2' 
c            stop   
            
            tol = 0.001d0
c            write(*,*)'avant golden',ax,bx,cx,tol,xmin1,n,b,ddl
c            stop
            res = golden(ax,bx,cx,tol,xmin1,n,b,y,ddl)
c            write(*,*)'apres golden'
c            stop

c            write(4,*)'******************************************* '
c            write(4,*)'Best smoothing parameter',real(xmin1*xmin1)
c     &                , '  DoF :',-ddl
c            write(4,*)'******************************************* '
c            write(*,*)'Best smoothing parameter',real(xmin1*xmin1)
c     &                , '  DoF :',-ddl
           
c            stop

            auxkappa(1)=xmin1*xmin1
            auxkappa(2)=0.d0

            call marq98(auxkappa,b,n,ni,v,res,ier,istop)

            if (ni.ge.maxit) then
c               write(*,*)' '
c              write(*,*) 'non convergence'
            else
c               write(*,*)' '
c               write(*,*)'Log-vraisemblance :',res
            endif
         endif   
        

ccccccccc
         nva=nvacross ! pour la recherche des parametres de regression
         nst=nstcross ! avec stratification si nécessaire
         effet=effetcross ! avec effet initial
         do 753 l=1,nsujet  
            stra(l)=stracross(l) !rétablissement stratification
 753     continue
c        stop

ccccc********************************************************************


         k0(1) = xmin1*xmin1
         if(nst.eq.2)then
            k0(2) = xmin2
         endif

         


C-------------  indicateur d'effet aleatoire ou non dans le modele
      
c          write(*,*) k0          
          call marq98(k0,b,np,ni,v,res,ier,istop)
  


           j=(np-nva)*(np-nva+1)/2       

      
             call multi(I_hess,H_hess,np,np,np,IH)
             call multi(H_hess,IH,np,np,np,HIH)

             if(effet.eq.1)then
               j=(np-nva)*(np-nva+1)/2
   	     endif
          
            if(effet.eq.1.and.ier.eq.-1)then
               v((np-nva)*(np-nva+1)/2)=10.d10
            endif

c vraisemblance non penalisee:            
c         res01(effet+1)=res+pe
c vraisemblance penalisee:            
         res01(effet+1)=res
         
c         write(4,*)'valeur de ni',ni 

          j=(np-nva)*(np-nva+1)/2





c --------------  Lambda and survival estimates JRG January 05

       

c             if(effet.eq.0)then
c               call distance(nz1,nz2,b,effet,
c     &	       x1,lam,lbinf,lbsup,margi,x2,lam2,lbinf2,lbsup2,margi2)
c	     else
c               call distance(nz1,nz2,b,effet,
c     &	       x1,lam,lbinf,lbsup,margi,x2,lam2,lbinf2,lbsup2,margi2)
c	     endif

             call distance(nz1,nz2,b,effet,
     &	       x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)      
             
 
      cpt3(ibou)=0
      if(ni.ge.999.and.effet.eq.0)then
         cpt1(ibou)=1 
      endif
      if(ni.ge.999.and.effet.eq.1)then
         cpt2(ibou)=1
      endif

c 765  continue

c-----fin de la boucle sur  effet aleatoire ou non 

      if(cpt1(ibou).eq.1.or.cpt2(ibou).eq.1)then
         cpt3(ibou)=1
c------cpt3=1 lorsque le nombre d''iterations est max pour un des 2 modeles
      endif
      
      
CCCCC*******************************
      if (cpt3(ibou).eq.1) then
         cptni=cptni+1
      endif  
      if (cpt1(ibou).eq.1) then
         cptni1=cptni1+1
      endif
 
      if (cpt2(ibou).eq.1) then
         cptni2=cptni2+1
         ind(ibou)=1
      else
         nbou2=nbou2+1  
      endif
c         cp = sp/sqrt(varsp)
      if (effet.eq.1)then
         tvars(ibou)= b(np-nva)*b(np-nva)
         j=(np-nva)*(np-nva+1)/2
c         write(*,*)'vj',v(j)
c         write(*,*)'bnp-nva 2',b(np-nva)*b(np-nva)
      if(H_hess(np-nva,np-nva).lt.2.and.H_hess(np-nva,np-nva).gt.-2)then
            se(ibou)=dsqrt(H_hess(np-nva,np-nva))
         endif
         if(HIH(np-nva,np-nva).lt.2.and.HIH(np-nva,np-nva).gt.-2)then
            se_cor(ibou)=dsqrt(HIH(np-nva,np-nva))
         endif
      endif
       
c     on retourne a l'iteration suivante
 1000 continue


      resOut=res

      do ss=1,npmax
       do sss=1,npmax
         HIHOut(ss,sss) = HIH(ss,sss)
         H_hessOut(ss,sss)= H_hess(ss,sss)
       end do  
      end do


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
C    Bias and Var eliminated  
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


                    
      return     
      
      end subroutine frailpenal
      
  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC





CCCCCCCCCCCCC**********SUBROUTINES******  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



c========================== VECSPLI ==============================
      subroutine vecspli(n,ndate) 
      
      integer   n,ndate,i,j,k
      double precision  ht,htm,h2t,ht2,ht3,hht,h,hh,h2
      double precision  h3,h4,h3m,h2n,hn,hh3,hh2
         
      parameter(ndatemax=30000,npmax=50,nsujetmax=15000)
      parameter(nvarmax=50)

c*****ve1
      double precision ve(nsujetmax,nvarmax)
      common /ve1/ve
c*****dace1 
      double precision date(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,zi
c*****mem1
      double precision mm3(ndatemax),mm2(ndatemax)
      double precision mm1(ndatemax),mm(ndatemax)
      common /mem1/mm3,mm2,mm1,mm
c*****mem2
      double precision im3(ndatemax),im2(ndatemax)
      double precision im1(ndatemax),im(ndatemax)
      common /mem2/im3,im2,im1,im

c----------  calcul de u(ti) ---------------------------
c    attention the(1)  sont en nz=1
c        donc en ti on a the(i)

         do 8 i=1,ndate-1
            do 6 k = 2,n-2
               if ((date(i).ge.zi(k-1)).and.(date(i).lt.zi(k)))then
                  j = k-1
               endif
 6          continue 
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
            mm2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm
     &       *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht*
     &       h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mm(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3(i) = (0.25d0*(date(i)-zi(j-3))*mm3(i))+(0.25d0*hh2
     &       *mm2(i))+(0.25d0*h3m*mm1(i))+(0.25d0*h4*mm(i))
            im2(i) = (0.25d0*hht*mm2(i))+(h3m*mm1(i)*0.25d0)
     &               +(h4*mm(i)*0.25d0)
            im1(i) = (htm*mm1(i)*0.25d0)+(h4*mm(i)*0.25d0)
            im(i)  = ht*mm(i)*0.25d0


 8       continue
         end
c========================== VECPEN ==============================
      subroutine vecpen(n) 

      parameter(ndatemax=30000,npmax=50)

      integer   n,i
      double precision  h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2
      double precision  a3,a2,b2,c2,a1,b1,c1,a0,x3,x2,x

c*****dace1 
      double precision date(ndatemax) 
      double precision zi(-2:npmax) 
      common /dace1/date,zi
c*****pen1
      double precision  m3m3(npmax),m2m2(npmax),m1m1(npmax)
      double precision  mmm(npmax),m3m2(npmax)
      common /pen1/m3m3,m2m2,m1m1,mmm,m3m2
c*****pen2
      double precision m3m1(npmax),m3m(npmax),m2m1(npmax)
      double precision m2m(npmax),m1m(npmax)
      common /pen2/m3m1,m3m,m2m1,m2m,m1m
c*********************************************************************
         
         do 20 i=1,n-3
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
            m2m2(i) = 64.d0*(((3.d0*x3-(3.d0*x2*(2.d0*zi(i+1)+zi(i-2)
     &      ))+x*(4.d0*zi(i+1)*zi(i+1)+zi(i-2)*zi(i-2)+4.d0*zi(i+1)
     &      *zi(i-2)))/(a2*a2)))
            m2m2(i) = m2m2(i) + 64.d0*(((3.d0*x3-(3.d0*x2*(zi(i+2)
     &      +zi(i-1)+zi(i+1)))+x*(zi(i+2)*zi(i+2)+zi(i-1)*zi(i-1)
     &      +zi(i+1)*zi(i+1)+2.d0*zi(i+2)*zi(i-1)+2.d0*zi(i+2)
     &      *zi(i+1)+2.d0*zi(i-1)*zi(i+1)))/(b2*b2)))
            m2m2(i) = m2m2(i) +64.d0*((3.d0*x3-(3.d0*x2*(2.d0*zi(i+2)
     &     +zi(i)))+x*(4.d0*zi(i+2)*zi(i+2)+zi(i)*zi(i)+4.d0*zi(i+2)
     &      *zi(i)))/(c2*c2))
            m2m2(i) = m2m2(i) +128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+2)
     &      +zi(i-1)+3.d0*zi(i+1)+zi(i-2)))+x*(2.d0*zi(i+1)*zi(i+2)
     &      +2.d0*zi(i+1)*zi(i-1)+2.d0*zi(i+1)*zi(i+1)+zi(i-2)*zi(i+2)
     &      +zi(i-2)*zi(i-1)+zi(i-2)*zi(i+1)))/(a2*b2))
            m2m2(i) = m2m2(i) + 128.d0*((3.d0*x3-(1.5d0*
     &      x2*(2.d0*zi(i+2)+zi(i)+2.d0*zi(i+1)+zi(i-2)))+x*
     &      (4.d0*zi(i+1)*zi(i+2)+2.d0*zi(i+1)*zi(i)+2.d0*zi(i-2)
     &      *zi(i+2)+zi(i-2)*zi(i)))/(a2*c2))
            m2m2(i) = m2m2(i) + 128.d0*((3.d0*x3-(1.5d0*x2
     &     *(3.d0*zi(i+2)+zi(i)+zi(i-1)+zi(i+1)))+x*(zi(i+2)*zi(i)+
     &      2.d0*zi(i-1)*zi(i+2)+zi(i)*zi(i-1)+2.d0*zi(i+1)*zi(i+2)
     &      +zi(i+1)*zi(i)+2.d0*zi(i+2)*zi(i+2)))/(b2*c2))
            m1m1(i) = 64.d0*((3.d0*x3-(3.d0*x2*(2.d0*zi(i-1)+zi(i+1)))
     &      +x*(4.d0*zi(i-1)*zi(i-1)+zi(i+1)*zi(i+1)+4.d0*zi(i-1)
     &      *zi(i+1)))/(a1*a1))
            m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i-1)+zi(i)     
     &      +zi(i+2)))+x*(zi(i-1)*zi(i-1)+zi(i)*zi(i)+zi(i+2)*
     &      zi(i+2)+2.d0*zi(i-1)*zi(i)+2.d0*zi(i-1)*zi(i+2)+2.d0*
     &      zi(i)*zi(i+2)))/(b1*b1))
            m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i+3)
     &      +2.d0*zi(i)))+x*(zi(i+3)*zi(i+3)+4.d0*zi(i)*zi(i)
     &      +4.d0*zi(i+3)*zi(i)))/(c1*c1)) 
            m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(3.d0
     &      *zi(i-1)+zi(i)+zi(i+2)+zi(i+1)))+x*(2.d0*zi(i-1)*zi(i-1)
     &      +2.d0*zi(i-1)*zi(i)+2.d0*zi(i-1)*zi(i+2)+zi(i+1)*zi(i-1)
     &      +zi(i+1)*zi(i)+zi(i+1)*zi(i+2)))/(a1*b1))
            m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+3)+
     &      2.d0*zi(i)+2.d0*zi(i-1)+zi(i+1)))+x*(2.d0*zi(i-1)*zi(i+3)
     &      +4.d0*zi(i-1)*zi(i)+zi(i+1)*zi(i+3)+2.d0*zi(i+1)*zi(i)))
     &       /(a1*c1))    
            m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+3)+3.d0
     &      *zi(i)+zi(i-1)+zi(i+2)))+x*(zi(i-1)*zi(i+3)+2.d0*zi(i-1)    
     &      *zi(i)+zi(i+3)*zi(i)+2.d0*zi(i)*zi(i)+zi(i+2)*zi(i+3)
     &      +2.d0*zi(i+2)*zi(i)))/(b1*c1))
            mmm(i) = (192.d0*h/(h4*h3*h2*h4*h3*h2))
            m3m2(i) = 192.d0*(((-x3+(0.5d0*x2*(5.d0*zi(i+1)+zi(i-2)
     &      ))-x*(2.d0*zi(i+1)*zi(i+1)+zi(i+1)*zi(i-2)))/(a3*a2))
     &      +((-x3+(0.5d0*x2*(4.d0*zi(i+1)+zi(i-1)+zi(i+2)))-x*
     &       (zi(i+1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*b2))
     &      +((-x3+(0.5d0*x2*(3.d0*zi(i+1)+2.d0*zi(i+2)+zi(i)))-x*
     &      (2.d0*zi(i+1)*zi(i+2)+zi(i+1)*zi(i)))/(a3*c2)) )
            m3m1(i) = 192.d0*(((x3-(0.5d0*x2*(4.d0*zi(i+1)+2.d0*zi(i-1)
     &      ))+x*(2.d0*zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*a1))
     &      +((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+2)+zi(i-1)+zi(i)))
     &      +x*(zi(i+1)*zi(i-1)+zi(i+1)*zi(i)+zi(i+1)*zi(i+2)))/(b1*a3))
     &     +((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+3)+2.d0*zi(i)))+x*(zi(i+1)
     &      *zi(i+3)+2.d0*zi(i+1)*zi(i)))/(c1*a3)) )
            m3m(i) = 576.d0*((-(x3/3.d0)+(0.5d0*x2*(zi(i+1)+zi(i)))
     &       -x*zi(i+1)*zi(i))/(a3*a0))
            m2m1(i) = 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0*zi(i-1)+3.d0*
     &      zi(i+1)+zi(i-2)))-x*(4.d0*zi(i+1)*zi(i-1)+2.d0*zi(i+1)
     &      *zi(i+1)+2.d0*zi(i-2)*zi(i-1)+zi(i-2)*zi(i+1)))/(a2*a1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i-1)+
     &      zi(i)+zi(i+2)+2.d0*zi(i+1)+zi(i-2)))-x*(2.d0*zi(i+1)*zi(i-1)
     &      +2.d0*zi(i+1)*zi(i)+2.d0*zi(i+1)*zi(i+2)+zi(i-2)*zi(i-1)+
     &      zi(i-2)*zi(i)+zi(i-2)*zi(i+2)))/(a2*b1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3)+2.d0
     &      *zi(i)+2.d0*zi(i+1)+zi(i-2)))-x*(2.d0*zi(i+1)*zi(i+3)+4.d0
     &      *zi(i+1)*zi(i)+zi(i-2)*zi(i+3)+2.d0*zi(i-2)*zi(i)))/(a2*c1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*
     &      (3.d0*zi(i-1)+2.d0*zi(i+1)+zi(i+2)))-x*(2.d0*zi(i+2)*zi(i-1)
     &      +zi(i+2)*zi(i+1)+2.d0*zi(i-1)*zi(i-1)+3.d0
     &       *zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(b2*a1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0
     &      *zi(i-1)+zi(i)+2.d0*zi(i+2)+zi(i+1)))-x*(zi(i+2)*zi(i-1)
     &      +zi(i+2)*zi(i)+zi(i+2)*zi(i+2)+zi(i-1)*zi(i-1)+zi(i-1)
     &      *zi(i)+zi(i-1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i)
     &      +zi(i+1)*zi(i+2)))/(b2*b1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3)
     &      +2.d0*zi(i)+zi(i+2)+zi(i-1)+zi(i+1)))-x*(zi(i+2)*zi(i+3)
     &      +2.d0*zi(i+2)*zi(i)+zi(i-1)*zi(i+3)+2.d0*zi(i-1)*zi(i)
     &      +zi(i+1)*zi(i+3)+2.d0*zi(i+1)*zi(i)))/(b2*c1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0*zi(i-1)
     &      +zi(i+1)+2.d0*zi(i+2)+zi(i)))-x*(4.d0*zi(i+2)*zi(i-1)+2.d0*
     &      zi(i+2)*zi(i+1)+2.d0*zi(i)*zi(i-1)+zi(i)*zi(i+1)))/(c2*a1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i-1)
     &      +2.d0*zi(i)+3.d0*zi(i+2)))-x*(2.d0*zi(i+2)*zi(i-1)+2.d0
     &      *zi(i+2)*zi(i)+2.d0*zi(i+2)*zi(i+2)+zi(i)*zi(i-1)+zi(i)
     &      *zi(i)+zi(i)*zi(i+2)))/(c2*b1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3)
     &      +3.d0*zi(i)+2.d0*zi(i+2)))-x*(2.d0*zi(i+2)*zi(i+3)+4.d0
     &       *zi(i+2)*zi(i)+zi(i)*zi(i+3)+2.d0*zi(i)*zi(i)))/(c2*c1))
            m2m(i) = 192.d0*(((x3-(0.5d0*x2*(3.d0*zi(i)+2.d0*zi(i+1)
     &      +zi(i-2)))+x*(2.d0*zi(i+1)*zi(i)+zi(i-2)*zi(i)))/(a2*a0))
     &      +((x3-(0.5d0*x2*(3.d0*zi(i)+zi(i+2)+zi(i-1)+zi(i+1)))
     &       +x*(zi(i+2)*zi(i)+zi(i-1)*zi(i)+zi(i+1)*zi(i)))/(b2*a0))
     &      +((x3-(0.5d0*x2*(4.d0*zi(i)+2.d0*zi(i+2)))+x*(2.d0*zi(i+2)
     &      *zi(i)+zi(i)*zi(i)))/(c2*a0)) )
            m1m(i) = 192.d0*(((-x3+(0.5d0*x2*(3.d0*zi(i)+2.d0*zi(i-1)
     &      +zi(i+1)))-x*(2.d0*zi(i-1)*zi(i)+zi(i+1)*zi(i)))/(a1*a0))
     &      +((-x3+(0.5d0*x2*(4.d0*zi(i)+zi(i-1)+zi(i+2)))
     &       -x*(zi(i-1)*zi(i)+zi(i)*zi(i)+zi(i+2)*zi(i)))/(b1*a0))
     &      +((-x3+(0.5d0*x2*(5.d0*zi(i)+zi(i+3)))-x*(zi(i+3)*zi(i)
     &      +2.d0*zi(i)*zi(i)))/(c1*a0)) )
 20       continue

          end




c========================          FUNCPA          ====================
      double precision function funcpa(b,np,id,thi,jd,thj,k0)

c *** NOUVELLLE DECLARATION F90 :

      parameter(npmax=50,nsujetmax=15000,nvarmax=50)
      parameter(ngmax=1000,nssgmax=1000)
      parameter(ndatemax=30000)

      integer  nb,n,np,id,jd,i,j,k,vj,cptg,l
      integer cpt(ngmax)
c      integer nan

      double precision  thi,thj,pe1,pe2,dnb,sum
      double precision  theta,inv
      double precision  som1,som2
      double precision  res,vet,h1
      double precision the1(-2:npmax),the2(-2:npmax)
      double precision b(np),bh(np)
      double precision res1(ngmax),res2(ngmax),res3(ngmax)
      double precision k0(2) 
      double precision dut1(ndatemax),dut2(ndatemax)
      double precision ut1(0:ndatemax),ut2(0:ndatemax)

c ************************* les commons :
c*****pen1
      double precision m3m3(npmax),m2m2(npmax)
      double precision m1m1(npmax),mmm(npmax),m3m2(npmax)
      common /pen1/m3m3,m2m2,m1m1,mmm,m3m2
c*****pen2
      double precision m3m1(npmax),m3m(npmax)
      double precision m2m1(npmax),m2m(npmax),m1m(npmax)
      common /pen2/m3m1,m3m,m2m1,m2m,m1m
c*****mem1
      double precision mm3(ndatemax),mm2(ndatemax)
      double precision mm1(ndatemax),mm(ndatemax)
      common /mem1/mm3,mm2,mm1,mm
c*****mem2
      double precision im3(ndatemax),im2(ndatemax),im1(ndatemax)
      double precision im(ndatemax)
      common /mem2/im3,im2,im1,im
c*****dace1 
      double precision date(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,zi
c*****dace2
      double precision t0(nsujetmax),t1(nsujetmax)
      double precision  ncsrl
      integer c(nsujetmax)
      integer nt0(nsujetmax),nt1(nsujetmax)
      integer  nsujet,nsim,nva,ndate,nst,maxit
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,
     &              nst,maxit
c*****dace4
      integer stra(nsujetmax)
      common /dace4/stra
c*****ve1
      double precision ve(nsujetmax,nvarmax)
      common /ve1/ve
c*****dace3
      double precision  pe
      integer  effet,nz1,nz2
      common /dace3/pe,effet,nz1,nz2
c*****contrib
      integer ng          !nb de gpes
      common /contrib/ng    
c*****groupe
      integer g(nsujetmax) 
      integer nig(ngmax) 
      common /gpe/g,nig
c %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
      integer AG
      common /andersengill/AG
c************************************************************        
c     write(*,*)'%%%%%%%%%% ENTREE DANS FUNCPA %%%%%%%%%%'

c      write(*,*)'* funcm3*',m3m(1),m2m1(1),m2m(1),m1m(1)
c      write(*,*)'* m3m3* ',m3m3(1),m2m2(1),m1m1(1),mmm(1),m3m2(1)
        
c      write(*,*)'* mm3* ',mm3(1),mm2(1),mm1(1),mm(1)
c      write(*,*)'* im3* ',im3(4)
c      stop
c      write(*,*)'nig',nig
         do 3 i=1,np
            bh(i)=b(i)
c             write(*,*)'funcpa bh(i)',bh(i),i,np
 3       continue 
c         stop
         if (id.ne.0) bh(id)=bh(id)+thi
         if (jd.ne.0) bh(jd)=bh(jd)+thj    
         
         n = (np-nva-effet)/nst
          
c         write(*,*)'funcpa nombre total de paramètres',np
c         write(*,*)'funcpa nb noeuds+2',nz1,nz2
         do 4 i=1,n
            the1(i-3)=(bh(i))*(bh(i))
            j = n+i 
c           write(*,*)'the1',the1(i-3)
             if (nst.eq.2) then
                the2(i-3)=(bh(j))*(bh(j))
             endif
 4        continue


        if(effet.eq.1) then
         theta = bh(np-nva)*bh(np-nva)
c         write(*,*)'theta',theta
         endif
c----------  calcul de ut1(ti) et ut2(ti) ---------------------------
c    attention the(1)  sont en nz=1
c        donc en ti on a the(i)

         vj = 0
         som1 = 0.d0
         som2 = 0.d0
         dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))

         dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))
         ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))
c      print*,'-- AUX --',the1(-2),dut1(1),zi(1),zi(-2)
c      stop
         ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))
         ut1(0) = 0.d0
         ut2(0) = 0.d0
         do 8 i=2,ndate-1
            do 6 k = 2,n-2
               if (((date(i)).ge.(zi(k-1))).and.(date(i).lt.zi(k)))then
                  j = k-1
                  if ((j.gt.1).and.(j.gt.vj))then
                     som1 = som1+the1(j-4)
                     som2 = som2+the2(j-4)
                     vj  = j
                  endif   
               endif
 6          continue 
            ut1(i) = som1 +(the1(j-3)*im3(i))+(the1(j-2)*im2(i))
     &       +(the1(j-1)*im1(i))+(the1(j)*im(i))
            dut1(i) = (the1(j-3)*mm3(i))+(the1(j-2)*mm2(i))
     &       +(the1(j-1)*mm1(i))+(the1(j)*mm(i))
c         write(*,*)'** ut1',ut1(i),som1,i,the1(j-3),im3(i),im2(i)
c             write(*,*)'** dut1',mm3(i),i
c         stop

c            write(4,*)date(i),dut1(i),ut1(i)

            if(nst.eq.2)then
            ut2(i) = som2 +(the2(j-3)*im3(i))+(the2(j-2)*im2(i))
     &       +(the2(j-1)*im1(i))+(the2(j)*im(i))
            dut2(i) = (the2(j-3)*mm3(i))+(the2(j-2)*mm2(i))
     &       +(the2(j-1)*mm1(i))+(the2(j)*mm(i)) 
c            if(i.eq.157.or.i.eq.158)then
c               write(*,*)'** dut2',dut2(i),mm2(i),mm1(i),mm(i)
c            endif
            endif
            
 8       continue
         i = n-2
         h1 = (zi(i)-zi(i-1))
      ut1(ndate)=som1+the1(i-4)+the1(i-3)+the1(i-2)+the1(i-1)
      ut2(ndate)=som2+the1(i-4)+the2(i-3)+the2(i-2)+the2(i-1)
         dut1(ndate) = (4.d0*the1(i-1)/h1)
         dut2(ndate) = (4.d0*the2(i-1)/h1)

c         write(*,*)'** ut2',(ut2(i),i=0,ndate)
c         write(*,*)'** dut2',(dut2(i),i,i=0,ndate)
c         stop
c         write(*,*)'dut1(ndate)',dut1(ndate),dut1(1)
c         write(*,*)'dut2(ndate)',dut2(ndate),dut2(1)
c_-------------------------------------------------------
c---------- calcul de la vraisemblance ------------------
c---------------------------------------------------------

cc---- avec ou sans variable explicative  ------cc

         do 89 k=1,ng
            res1(k) = 0.d0
            res2(k) = 0.d0
            res3(k) = 0.d0
            cpt(k) = 0
 89      continue

c*******************************************     
C----- sans effet aleatoire dans le modele
c*******************************************     

        if (effet.eq.0) then
           do 110 i=1,nsujet
                   cpt(g(i))=cpt(g(i))+1
                   
                   if(nva.gt.0)then
                      vet = 0.d0   
                      do 91 j=1,nva
                         vet =vet + bh(np-nva+j)*dble(ve(i,j))
c       write(*,*)'funcpa ve ******',ve(i,j),i,j
 91                   continue
                      vet = dexp(vet)
                   else
                      vet=1.d0
                   endif
                  
                   if((c(i).eq.1).and.(stra(i).eq.1))then
                     res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet)
c               write(*,*)'FUNC res2',res2(g(i)),g(i),stra(i),nt1(i) 
                  endif  

                   if((c(i).eq.1).and.(stra(i).eq.2))then
                      res2(g(i)) = res2(g(i))+dlog(dut2(nt1(i))*vet)
                  endif

                   if(stra(i).eq.1)then
                      res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet
     &                     -ut1(nt0(i))*vet 
c               write(*,*)'FUNC res1',res1(g(2)),g(2),stra(2),nt1(2)
                   endif

                   if(stra(i).eq.2)then
                      res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet
     &                     -ut2(nt0(i))*vet 
c                   write(*,*)'FUNC res1',res1(g(i)),g(i),stra(i),nt1(i)
                   endif
           
 110           continue       
c           write(*,*)'FUNC fin res1',res1(g(2)),g(2),stra(2),nt1(2)
c           write(*,*)'FUNC',ut1(nt1(i)) ,ut1(nt0(i)) 
c           write(*,*)'FUNC res2',res2(g(1)),vet
c           write(*,*)'FUNC vet',vet,bh
c           stop
           res = 0.d0
           
          cptg = 0
          
C k indice les groupes
          do 115 k=1,ng   
             if(cpt(k).gt.0)then
                nb = nig(k)
                dnb = dble(nig(k))
                
                  res = res-res1(k) 
     &               + res2(k) 
                cptg = cptg + 1 
             endif 
 115       continue
  
c*******************************************         
C-----avec un effet aleatoire dans le modele
c*********************************************

       else
c      write(*,*)'AVEC EFFET ALEATOIRE'
          inv = 1.d0/theta
c     write(*,*)'INV 1============',inv
C     i indice les sujets
          do 10 i=1,nsujet 
             
             cpt(g(i))=cpt(g(i))+1 
c             write(4,*)'funcpa cpt ******',cpt(g(i)),g(i)
 
             if(nva.gt.0)then
                vet = 0.d0   
                do 19 j=1,nva
                   vet =vet + bh(np-nva+j)*dble(ve(i,j))
c                     write(*,*)'funcpa ve ******',ve(i,j),i,j
 19             continue
                vet = dexp(vet)
             else
                vet=1.d0
             endif
             if((c(i).eq.1).and.(stra(i).eq.1))then
                res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet)
c           write(*,*)'***res2',res2(g(i)),dut1(nt1(i)),nt1(i),vet,i,g(i)
             endif  
             if((c(i).eq.1).and.(stra(i).eq.2))then
                res2(g(i)) = res2(g(i))+dlog(dut2(nt1(i))*vet)
             endif  
             if(stra(i).eq.1)then
                res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet
c                write(*,*)'***res1',res1(g(i)),ut1(nt1(i)),nt1(i),i
             endif
             
             if(stra(i).eq.2)then
                res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet
             endif
c modification pour nouvelle vraisemblance / troncature:
             if(stra(i).eq.1)then
                res3(g(i)) = res3(g(i)) + ut1(nt0(i))*vet
             endif
             
             if(stra(i).eq.2)then
                res3(g(i)) = res3(g(i)) + ut2(nt0(i))*vet
             endif
             
 10       continue 
    
c           write(*,*)'***FUNC res2',(res2(i),i=1,ng)
                  
          res = 0.d0
          cptg = 0
c     gam2 = gamma(inv)
C k indice les groupes
          do 15 k=1,ng  
             sum=0.d0
             if(cpt(k).gt.0)then
                nb = nig(k)
                dnb = dble(nig(k))
c                write(*,*)'nb,dnb',nb,dnb,k
c     gam1 = gamma(dnb + inv) 
                
                if (dnb.gt.1.d0) then
                   do 16 l=1,nb
                      sum=sum+dlog(1.d0+theta*dble(nb-l))
 16                continue
                endif
                if(theta.gt.(1.d-5)) then
ccccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                   if(AG.EQ.1)then
                   res= res-(inv+dnb)*dlog(theta*(res1(k)-res3(k))+1.d0) 
     &                  + res2(k) + sum  
ccccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                else
                   res = res-(inv+dnb)*dlog(theta*res1(k)+1.d0) 
     &                  +(inv)*dlog(theta*res3(k)+1.d0) 
     &                  + res2(k) + sum  
                   endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc              
c                   write(*,*)'res',res
c                   write(*,*)'vet',vet
c                   write(*,*)'inv+dnb',inv+dnb
c                   write(*,*)'theta',theta
c                   write(*,*)'res1',res1(k)
c                   write(*,*)'dlog',dlog(theta*res1(k)+1.d0) 
c                   write(*,*)'inv',inv
c                   write(*,*)'dlog2',dlog(theta*res3(k)+1.d0) 
c                   write(*,*)'res2',res2(k)
c                   write(*,*)'sum',sum
c                   stop
                else              
c     developpement de taylor d ordre 3
c                   write(*,*)'************** TAYLOR *************'
ccccc ancienne vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                     if(AG.EQ.1)then
                   res = res-dnb*dlog(theta*(res1(k)-res3(k))+1.d0)
     &             -(res1(k)-res3(k))*(1.d0-theta*(res1(k)-res3(k))/2.d0
     &            +theta*theta*(res1(k)-res3(k))*(res1(k)-res3(k))/3.d0)
     &                  +res2(k)+sum

ccccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                     else
                   res = res-dnb*dlog(theta*res1(k)+1.d0)
     &                  -res1(k)*(1.d0-theta*res1(k)/2.d0
     &                            +theta*theta*res1(k)*res1(k)/3.d0)
     &                  +res2(k)+sum
     &                  +res3(k)*(1.d0-theta*res3(k)/2.d0
     &                            +theta*theta*res3(k)*res3(k)/3.d0)
                endif
             endif 
          endif 
 15    continue

       endif !fin boucle effet=0

c---------- calcul de la penalisation -------------------

         pe1 = 0.d0
         pe2 = 0.d0
         do 20 i=1,n-3
            pe1 = pe1+(the1(i-3)*the1(i-3)*m3m3(i))+(the1(i-2)
     &    *the1(i-2)*m2m2(i))+(the1(i-1)*the1(i-1)*m1m1(i))+(
     &    the1(i)*the1(i)*mmm(i))+(2.d0*the1(i-3)*the1(i-2)*
     &    m3m2(i))+(2.d0*the1(i-3)*the1(i-1)*m3m1(i))+(2.d0*
     &    the1(i-3)*the1(i)*m3m(i))+(2.d0*the1(i-2)*the1(i-1)*
     &    m2m1(i))+(2.d0*the1(i-2)*the1(i)*m2m(i))+(2.d0*the1(i-1)
     &    *the1(i)*m1m(i))
            pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2)
     &    *the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+(
     &    the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)*
     &    m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0*
     &    the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)*
     &    m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1)
     &    *the2(i)*m1m(i))
 20       continue
    

c    Changed JRG 25 May 05
          if (nst.eq.1) then
            pe2=0.d0
          end if
          
          pe = k0(1)*pe1 + k0(2)*pe2 
              

c          write(*,*)'avant  penalisation ',res,pe,pe1,pe2
          res = res - pe
c           write(*,*)'vraisemblance penalisee',res,pe,pe1,pe2
c           write(*,*)'vraisemblance penalisee',k0(1),k0(2)
          funcpa = res 

c          if(funcpa.eq.nan)then 
c             write(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù'
c             write(*,*)'%%%%%%%%%%%% PB FUNCPA :NAN %%%%%%%%%ù'
c             write(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù'
c             stop
c          endif
        
c      write(*,*)'%%%%%%%%%% SORTIE DANS FUNCPA %%%%%%%%%%' ,res
c      stop
          return
          end



c================================  SEARPAS joly    ==============================

      SUBROUTINE SEARPAS(VW,STEP,B,BH,M,DELTA,FIM,EPSV,k0)
C
C  MINIMISATION UNIDIMENSIONNELLE
C
       INTEGER  I,M
       DOUBLE PRECISION  VLW,VLW1,VLW2,VLW3,VW,VM
       DOUBLE PRECISION  FI1,FI2,FI3,FIM,EPSV
       DOUBLE PRECISION  STEP
       DOUBLE PRECISION k0(2)
       DOUBLE PRECISION B(M),BH(M),DELTA(M)
C
       VLW1=DLOG(VW)
       VLW2=VLW1+STEP
       CALL VALFPA(VLW1,FI1,B,BH,M,DELTA,k0)
       CALL VALFPA(VLW2,FI2,B,BH,M,DELTA,k0)
C
       IF(FI2.GE.FI1) THEN
          VLW3=VLW2
          VLW2=VLW1
          FI3=FI2
          FI2=FI1
C
          STEP=-STEP
C
          VLW1=VLW2+STEP
          CALL VALFPA(VLW1,FI1,B,BH,M,DELTA,k0)   
          IF (FI1.GT.FI2) GOTO 50
       ELSE
          VLW=VLW1
          VLW1=VLW2
          VLW2=VLW
          FIM = FI1
          FI1 = FI2
          FI2 = FIM
       ENDIF
C
       DO 20 I=1,40
          VLW3=VLW2
          VLW2=VLW1
          FI3=FI2
          FI2=FI1
C
          VLW1=VLW2+STEP
          CALL VALFPA(VLW1,FI1,B,BH,M,DELTA,k0)
          IF(FI1.GT.FI2) GO TO 50
c          IF (dabs(FI1-FI2).LT.EPSV) THEN
          IF (FI1.eq.FI2) THEN
             FIM=FI2
             VM=VLW2
             GO TO 100
          ENDIF
 20    CONTINUE
C
C  PHASE 2 APPROXIMATION PAR QUADRIQUE
C
50     CONTINUE
C
C  CALCUL MINIMUM QUADRIQUE
C
         VM=VLW2-STEP*(FI1-FI3)/(2.d0*(FI1-2.d0*FI2+FI3))
         CALL VALFPA(VM,FIM,B,BH,M, DELTA,k0)
         IF (FIM.LE.FI2) GO TO 100
         VM=VLW2
         FIM=FI2
100   CONTINUE
      VW=DEXP(VM)
      RETURN

      END




c===================================   VALFPA joly   ==============================

        subroutine valfpa(vw,fi,b,bk,m,delta,k0)
        integer  m,i
        double precision  vw,fi
        double precision  funcpa,z        
        double precision k0(2)
        double precision b(M),bk(M),delta(M)


         z=0.d0
         do 1 i=1,m
            bk(i)=b(i)+dexp(vw)*delta(i)
1        continue
c
         fi=-funcpa(bk,m,1,z,1,z,k0)
c
         return
         end   




c===============================    MARQ98  HJ =========================

   	 subroutine marq98(k0,b,m,ni,v,rl,ier,istop)
c
c
c  fu = matrice des derivees secondes et premieres
c
c  istop: raison de l'arret
c  1: critere d'arret satisfait (prm=ca, vraisblce=cb, derivee=dd)
c  2: nb max d'iterations atteints
c  3: 1 mais echec inversion matrice d'info (ier=1)
C  4: non amelioration vraisblce apres searpas
c      
         parameter(npmax=50,nsujetmax=15000,nvarmax=50)

         integer  m,ni,nq,i,ii,nfmax,idpos,ier,istop,igrad,j
         integer  ncount,id,jd,i0,kkk


         double precision  da,dm,tr,g0
         double precision  ca,cb,epsa,epsb,rl
         double precision  funcpa,det, step,eps,epsd
         double precision  vw,fi,maxt,z
         double precision  rl1,th,ep,dd

         double precision v((npmax*(npmax+3)/2))
         double precision v1((npmax*(npmax+3)/2))
         double precision vnonpen((npmax*(npmax+3)/2))
         double precision b(npmax),delta(npmax)
         double precision fu((npmax*(npmax+3)/2))
         double precision k0(2)
         double precision bh(npmax),b1(npmax)
         double precision zero(2)

         double precision ga
c*****dace2
      double precision t0(nsujetmax),t1(nsujetmax)
      double precision  ncsrl
      integer c(nsujetmax)
      integer nt0(nsujetmax),nt1(nsujetmax)
      integer  nsujet,nsim,nva,ndate,nst,maxit
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,
     &              nst,maxit
c*****dace3
      double precision  pe
      integer  effet,nz1,nz2
      common /dace3/pe,effet,nz1,nz2  
c*****dace7
      double precision I_hess(npmax,npmax),H_hess(npmax,npmax)
      double precision Hspl_hess(npmax,npmax)
      double precision PEN_deri(npmax,1) 
      double precision hess(npmax,npmax)
      common /dace7/PEN_deri,I_hess,H_hess,Hspl_hess,hess
c***********************************************
     
c      ga=0.d0
c
c      write(*,*)'Marq'
c      write(*,*)'ga 1=',ga
c      stop
      zero(1)=0.d0
      zero(2)=0.d0
      id=0
      jd=0
      z=0.d0
      th=1.d-5
      eps=1.d-7
      epsa=1.d-3
      epsb=1.d-3
      epsd=1.d-3
      nfmax=m*(m+1)/2
      ca=epsa+1.d0
      cb=epsb+1.d0
      rl1=-1.d+10
      ni=0
      istop=0
      da=0.01d0
      dm=5.d0
c     nq=nbre de seconds membres pour chole
      nq=1
 10   continue
c     
c      write(6,*)
c      write(*,*)'** Iteration :', ni 
c      if (ni.eq.6) stop
c      write(*,*)'** Parametres:'
c      write(*,*)(b(i),i=1,m)
c      write(*,*)'da',da,' ga=',ga
      
      z=0.d0
      i0=0 

      
      rl=funcpa(b,m,i0,z,i0,z,k0)
      rl1=rl 

c      write(*,*)'log Vrais rl1,nsujet',rl1,nsujet
c      stop

      if(ni.eq.0) write(8,*)'vrais init=',rl1
     
      if (RL1.gt.0)stop
     
      call deriva(b,m,v,rl,k0)
       
c          if (ni.eq.1) stop
          dd = 0.d0
          do 13 i=m*(m+1)/2+1,m*(m+3)/2
c             write(*,*)'d',v(i)
             dd = dd + v(i)*v(i)
 13       continue 
          dd=dd/dabs(RL)

c          pause
c          write(6,*) (b(ii),ii=1,m)
          if (ca.lt.epsa.and.cb.lt.epsb.and.dd.lt.epsd) goto 100
c

          tr=0.d0
          do 300 i=1,m
             ii=i*(i+1)/2
             tr=tr+dabs(v(ii))
 300      continue 
         tr=tr/dble(m)
c
          ncount=0
          ga=0.01d0
 400      do 1 i=1,nfmax+m
             fu(i)=v(i)
c             write(*,*)'fu',fu(i),i
1         continue
          do 500 i=1,m
             ii=i*(i+1)/2
             if (v(ii).ne.0) then
             fu(ii)=v(ii)+da*((1.d0-ga)*dabs(v(ii))+ga*tr) 
             
             else
             fu(ii)=da*ga*tr
             endif
 500      continue
          call dchole(fu,m,nq,idpos)
c
          if (idpos.ne.0) then  
c             write(6,*) 'echec inversion choleski'
c modif ligne suivante le 25/10/99 (da*5 et non *25 a chaque iter)
c             da=dm*da
             ncount=ncount+1
             if (ncount.le.3.or.ga.ge.1.d0) then
                da=da*dm
             else
                ga=ga*dm
                if (ga.gt.1.d0) ga=1.d0
             endif
             goto 400
          else
c              write(6,*) 'matrice reguliere, idpos=',idpos
c              write(*,*)'ncount=',ncount, 'da=',da, 'ga=',ga
              do 2 i=1,m
                 delta(i)=fu(nfmax+i)
                 b1(i)=b(i)+delta(i)
2             continue
c              write(6,*) '** ** avant func.....',b1,delta(1)
              rl=funcpa(b1,m,id,z,jd,z,k0)
c              write(6,*) 'rl1 =',rl1,' rl =',rl
              igrad=1
              if (rl1.lt.rl) then
                 if(da.lt.eps) then
                    da=eps
                 else
                    da=da/(dm+2.d0)
                 endif
                 goto 800
              endif
           endif
c           write(6,*) ' non amelioration de la vraisemblance '
           
           if(dabs(maxt(delta,m)).eq.0.d0) then
              vw=1.D-5
           else
              vw=th/maxt(delta,m)
           endif
c           
           step=dlog(1.5d0)
c           write(*,*)'vw ',vw,'th maxt ',th,maxt(delta,m)
c           write(*,*)'delta',(delta(i),i=1,m)
c           pause
c           write(*,*) 'searpas'
           call searpas(vw,step,b,bh,m,delta,fi,eps,k0) 
c           write(*,*) 'apres searpas',-fi,rl1
           rl=-fi
           IF(rl1.gt.rl) then
c              write(6,*) 'non amelioration vraisblce apres searpas'
c              write(6,*) 'ancienne rl1 =',rl1,' nouvelle rl =',rl
c              write(6,*) 'anciens prms =', (B(i),i=1,m)
c              write(6,*) 'delta =', (delta(i),i=1,m)
c              write(6,*)'ncount',ncount
c              pause
c              rl=rl1
c              istop=4
c              goto 110
            endif
           do 4 i=1,m
              delta(i)=vw*delta(i)
 4         continue
           da=(dm-3.d0)*da  

 800       cb =dabs(rl1-rl)
           ca=0.d0
           do 5 i=1,m
              ca=ca+delta(i)*delta(i)
 5         continue
c           write(6,*)'rl1=',rl1,'rl=',rl,'(rl1-rl)',(rl1-rl)
c           write(6,*) 'ca =',ca,' cb =',cb,' dd =',dd
    
           do 6 i=1,m
              b(i)=b(i)+delta(i)
c              write(*,*)'delta',delta(i),i
 6         continue

           ni=ni+1
           if (ni.gt.maxit) then
              istop=2
c              write(6,*) 'nombre iteration max atteinte'
              goto 110
           end if
           goto 10
c********************
c     inversion matrice d'information
c********************
100      continue
           istop=1
          
c================ pour les bandes de confiance
           call deriva(b,m,v,rl,k0)
           do 303 i=1,(m*(m+3)/2)
              v1(i)=0.d0
 303       continue
           
           m1=m-nva-effet

           kkk=m1*(m1+1)/2
           do 3005 i=1,m1
              kkk=kkk+1
              do 3004 j=i,m1
                 k = (((j-1)*j)/2) +i
c                 write(*,*)'kkk',kkk,k,i,j,b(i),b(j)
c                 hess(i,j)=v(k)/(4.d0*b(i)*b(j))
                 v1(k)=v(k)/(4.d0*b(i)*b(j))
 3004         continue
              v1(kkk)=v1(kkk)+(v(kkk)/(4.d0*b(i)*b(i)*b(i)))
c              hess(i,i)=hess(i,i)+(v(kkk)/(4.d0*b(i)*b(i)*b(i)))
c              write(*,*)'hess(i,i)',v1(kkk),v(kkk),kkk
 3005      continue 
c           stop
           do 304 i=1,(m1*(m1+3)/2)
c              v1(i)=-v1(i)
 304       continue
c           write(*,*)'hess conf:v1 avant dsinv',v1 !(1),v1(2),v1(3),v1(4),v1(5) ! qqs NAN 
c           stop
           ep=10.d-10
           call dsinv(v1,m1,ep,ier)
           if (ier.eq.-1) then
c              write(*,*)   'echec inversion matrice information'
              istop=3
           endif

           do 3014 i=1,m1
              do 3015 j=i,m1
                 hess(i,j)=v1((j-1)*j/2+i)
 3015         continue
 3014      continue

           do 3016 i=2,m1
              do 3017 j=1,i-1
                 hess(i,j)=hess(j,i)
 3017         continue
 3016      continue

c           do 3018 i=1,m
c              do 3019 j=1,m
c                 hess(i,j)=-hess(i,j)
c                 write(*,*)'hess',hess(i,j),i,j
c 3019         continue
c 3018      continue
c           write(*,*)'hess dans marq',hess(1,1)
c           stop
c           stop
c======== fin 


c           DO k=1,npl*(npl+1)/2
c             v1(k)=v(k)
c           END DO
c           write(*,*)'ca,cb,dd final',ca,cb,dd
c           write(8,*)'ca,cb,dd final',ca,cb,dd
           ep=10.d-10
           call dsinv(v,m,ep,ier)
           if (ier.eq.-1) then
c             write(6,103)          
103          format(1x,'echec inversion matrice information')
             istop=3
c             call dsinv(v1,npl,ep,ier)
c               if (ier.eq.-1) then
c             write(*,*)'echec inversion matrice information
c     & prms fixes'
c                istop=31
c               else
c                 DO k=1,npl*(npl+1)/2
c                  v(k)=v1(k)
c                 END DO
c               end if
          endif
c******************************************
c            do 104 i=1,m
c               ii=i*(i+1)/2
c               write(6,*) 'variance parametre avant deriva ',i,'=',v(ii)
c 104        continue

c            write(6,*)
c            ii=0
c            do 105 i=m*(m+1)/2+1,m*(m+3)/2
c               ii=ii+1
c               write(6,*) 'derivee premiere % parametre',ii,' = ',v(i)
c 105        continue
c========================================================
c pour pouvoir calculer la matrice I_hess  = 
c (-) la hessienne sur la vraisemblances non penalisee (information).
          
         ep=10.d-10
         call deriva(b,m,vnonpen,rl,zero)
         
         do 101 i=1,m
            do 102 j=i,m
               I_hess(i,j)=vnonpen((j-1)*j/2+i)
 102        continue
 101     continue
         
         do 111 i=2,m
            do 122 j=1,i-1
               I_hess(i,j)=I_hess(j,i)
 122        continue
 111     continue
           
c========================================================

c   H_hess est moins la hessienne inverse sur la vraisemblance penalisee
       do 106 i=1,m
              do 107 j=i,m
                 H_hess(i,j)=v((j-1)*j/2+i)
 107          continue
 106       continue
c       write(*,*) 'H_hess(16,16) fin marq',H_hess(16,16),m,((j-1)*j/2+i)
           do 108 i=2,m
              do 109 j=1,i-1
                 H_hess(i,j)=H_hess(j,i)
 109          continue
 108       continue

 110       continue

           return
       end

       
c=================================    DERIVA pjoly =======================

      subroutine deriva(b,m,v,rl,k0)

      parameter(npmax=50)
       
      integer i0,iun,m,m1,ll,i,k,j
      double precision  funcpa
      double precision thn,th
      double precision z
      double precision rl
      double precision vl
      double precision vaux
      double precision  th2
      double precision k0(2)
      double precision fcith(m)
      double precision b(m)
      double precision v((npmax*(npmax+3)/2))
c      double precision v(1000)

c
c     v:matrice d'information+score
c     calcul de la derivee premiere
c
c      print*,'entree deriva'
      th=1.d-5
      thn=-th
      th2=th*th
      z=0.d0
      i0=0
      iun =1
      rl=funcpa(b,m,iun,z,iun,z,k0)
      do 2 i=1,m
	    fcith(i)=funcpa(b,m,i,th,i0,z,k0)
c            print*,'fcith',fcith(i),i
2     continue

      k=0
      m1=m*(m+1)/2
      ll=m1
      do 1 i=1,m
      ll=ll+1
      vaux=funcpa(b,m,i,thn,i0,z,k0)
c      write(*,*)'vaux ok',vaux
      vl=(fcith(i)-vaux)/(2.d0*th)
c      write(*,*)'vl / deriva',vl,ll,i
c      write(*,*)'detail',fcith(i),vaux,(fcith(i)-vaux)
c      write(*,*)'vl',(fcith(i)-vaux)/(2.d0*th)
c      vl=(fcith(i)-funcpa(b,m,i,thn,i0,z,k0))/(2.d0*th)
      v(ll)=vl
      do 1 j=1,i
      k=k+1
      v(k)=-(funcpa(b,m,i,th,j,th,k0)-fcith(j)-fcith(i)+rl)/th2
1     continue
      
c      stop
      return
      end


c==========================  DISTANCE   =================================
    
   
         subroutine distance(nz1,nz2,b,effet,
     &	       x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)



         parameter(ndatemax=30000,npmax=50,nsujetmax=15000)

         integer  nz1,nz2,nz3,i,j,n,np,k,l,effet
         integer indx(71) 

         double precision  su1,ri1,su2,ri2,x1,x2,h1,h2
         double precision  d,h
         double precision  su,bsup,binf,lam,lbinf,lbsup,margi
         double precision HIH(npmax,npmax),IH(npmax,npmax)
          double precision hes1(npmax,npmax),hes2(npmax,npmax)
         double precision the1(-2:npmax),the2(-2:npmax)
         double precision b(npmax)
         character*14 fic1,fic2,fic1b,fic2b


         double precision  x1Out(99),lamOut(99,3),suOut(99,3),
     &                     x2Out(99),lam2Out(99,3),su2Out(99,3)

c*****dace1 
         double precision date(ndatemax)
         double precision zi(-2:npmax)
         common /dace1/date,zi
c*****dace2
      double precision t0(nsujetmax),t1(nsujetmax)
      double precision  ncsrl
      integer c(nsujetmax)
      integer nt0(nsujetmax),nt1(nsujetmax)
      integer  nsujet,nsim,nva,ndate,nst,maxit
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,
     &              nst,maxit
c*****dace7
         double precision I_hess(npmax,npmax),H_hess(npmax,npmax)
         double precision Hspl_hess(npmax,npmax)
         double precision PEN_deri(npmax,1)
         double precision hess(npmax,npmax) 
         common /dace7/PEN_deri,I_hess,H_hess,Hspl_hess,hess
c***********************************************
  
         n  = nz1+2
         if(nst.eq.2)then      
            np  = nz1+2+nz2+2+effet+nva
         else
            np  = nz1+2+effet+nva
         endif
            
        
             do 116 i=1,nz1+2
                do 115 j=1,nz1+2
                   hes1(i,j)=hess(i,j)
 115            continue
 116         continue 

         if(nst.eq.2)then  
            k = 0
            do 118 i=nz1+3,nz1+2+nz2+2
               k = k + 1 
               l = 0
               do 117 j=nz1+3,nz1+2+nz2+2
                  l = l + 1
                  hes2(k,l)=hess(i,j)
 117           continue
 118        continue   
         endif

         do 4 i=1,nz1+2
            the1(i-3)=(b(i))*(b(i))
 4       continue
         if(nst.eq.2)then  
            do 6 i=1,nz2+2
               j = nz1+2+i
               the2(i-3)=(b(j))*(b(j))
 6          continue
         endif

         h = (zi(n)-zi(1))*0.01d0
         x1 = zi(1)
         x2 = zi(1)     
         
         do 10 i=1,99 
             x1 = x1 + h 
             call cosp(x1,the1,nz1+2,hes1,zi,date,binf,su,bsup
     &              ,lbinf,lam,lbsup)
                         
               if(binf.lt.0.d0)then
                 binf = 0.d0
               endif
               if(bsup.gt.1.d0)then
                 bsup = 1.d0
               endif
              if(lbinf.lt.0.d0)then
                 lbinf = 0.d0
              endif

ccc   Replaced by next sentences and add new ones JRG January 05
c
c
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
               call cosp(x2,the2,nz2+2,hes2,zi,date,binf,su,bsup
     &              ,lbinf,lam,lbsup)
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
 10      continue
              
      return
         
      end subroutine



c==========================  SUSP  ====================================
      subroutine susp(x,the,n,su,lam,zi)
      parameter(ndatemax=30000,npmax=50)

      integer  j,k,n,i
      double precision  x,ht,ht2,h2,som,lam,su
      double precision  htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3
      double precision  ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2
      double precision  h,gl,hh
      double precision zi(-2:npmax),the(-2:npmax)


         som = 0.d0
         gl = 0.d0
          


         do 58 k = 2,n+1
            if ((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
               j = k-1
               if (j.gt.1)then
                  do 55 i=2,j
                     som = som+the(i-4)
 55               continue  
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
            mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm
     &       *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht*
     &       h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3 = (0.25d0*(x-zi(j-3))*mm3)+(0.25d0*hh2*mm2)
     &       +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
            im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
            im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
            im  = ht*mm*0.25d0
            gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)
     &       +(the(j)*im)
            lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)
     &       +(the(j)*mm)
            endif
 58      continue
   
         if(x.ge.zi(n))then
            som = 0.d0
            do 59 i=1,n+1
               som = som+the(i-3)
 59         continue
            gl = som
         endif
         su  = dexp(-gl)
         return
         end

c==========================  COSP  ====================================
c calcul les points pour les fonctions 
c et leur bandes de confiance

      subroutine cosp(x,the,n,y,zi,date,binf,su,bsup,lbinf,lam,lbsup)
      parameter(ndatemax=30000,npmax=50)
      integer  j,k,n,i
      double precision  x,ht,ht2,h2,som,lam,su
      double precision  binf,bsup,lbinf,lbsup,pm
      double precision  htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3
      double precision  ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2
      double precision  h,gl,hh
      double precision the(-2:npmax),zi(-2:npmax)
      double precision y(npmax,npmax)
      double precision date(ndatemax)

c      do 4 i=1,n
c         the(i-3)=(b(i))*(b(i))
c 4    continue

         som = 0.d0
         gl = 0.d0 
 
         do 58 k = 2,n-1
            if ((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
               j = k-1
               if (j.gt.1)then
                  do 55 i=2,j
                     som = som+the(i-4)
 55               continue  
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
            mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm
     &       *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht*
     &       h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3 = (0.25d0*(x-zi(j-3))*mm3)+(0.25d0*hh2*mm2)
     &       +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
            im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
            im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
            im  = ht*mm*0.25d0
            gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)
     &       +(the(j)*im)
            lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)
     &       +(the(j)*mm)
            endif
 58      continue
   
         if(x.ge.zi(n))then
            som = 0.d0
            do 59 i=1,n
               som = som+the(i-3)
 59         continue
            gl = som
         endif

         call conf(x,j,n,y,pm,zi,date)

         binf = dexp(-gl + 1.96d0*pm)
         su  = dexp(-gl)
         bsup = dexp(-gl - 1.96d0*pm)

         call conf1(x,j,n,y,pm,zi,date)
         lbinf = lam - 1.96d0*pm
         lbsup = lam + 1.96d0*pm
c         write(*,*)'lbinf apres conf1',lbinf,lam,pm

         return

         end
c=====================  CONF1  =============================
      subroutine  conf1(x,ni,n,y,pm,zi,date) 
       parameter(npmax=50)
       parameter(ndatemax=30000)

      integer  ni,i,n,j

      double precision  mmsp,x,pm
      double precision  res
      double precision date(ndatemax) 
      double precision zi(-2:npmax)
      double precision  vecti(npmax),aux(npmax)
      double precision y(npmax,npmax) 
                 
            do 10 i=1,n
               vecti(i) = mmsp(x,ni,i,zi,date)
c               write(*,*)' vecti(i)' ,vecti(i),i / ok idem
 10         continue
          
            do 20 i=1,n
               aux(i) = 0.d0
               do 15 j=1,n
                  aux(i) = aux(i) - y(i,j)*vecti(j)
c                  write(*,*)'aux(i)', y(i,j),i,j !/non y
 15            continue
 20         continue 
c            stop

            res = 0.d0
            do 30 i=1,n
               res = res + aux(i)*vecti(i)
c            write(*,*)'conf1 : res ',aux(i),vecti(i),i
 30         continue
c            write(*,*)'conf1 : res ',res
c            stop
c            if (res.lt.0)then 
               res=-res
c            endif 
               pm = dsqrt(res)
             
c            write(*,*)'conf1 pm',pm 
         end
c=====================  CONF  =============================
      subroutine  conf(x,ni,n,y,pm,zi,date)
       parameter(npmax=50)
       parameter(ndatemax=30000)
       
       integer  ni,i,n,j
       double precision isp,x,pm
       double precision res
       double precision zi(-2:npmax)
       double precision date(ndatemax)
       double precision vecti(npmax),aux(npmax)
       double precision y(npmax,npmax)

            do 10 i=1,n
               vecti(i) = isp(x,ni,i,zi,date)
 10         continue   

            do 20 i=1,n
               aux(i) = 0.d0
               do 15 j=1,n
                  aux(i) = aux(i) - y(i,j)*vecti(j)
c                 write(*,*)'aux(i) conf1', y(i,j),i,j
 15            continue
 20         continue   

            res = 0.d0
            do 30 i=1,n
               res = res + aux(i)*vecti(i)
 30         continue

c            write(*,*)'conf : res ',res
c            if (res.lt.0)then 
               res=-res
c            endif 
c               pm = dsqrt(2.d0*res)
               pm = dsqrt(res)
              
c            write(4,*)'conf : pm ',pm
               
         end


c==========================   ISP   ==================================
          double precision function isp(x,ni,ns,zi,date)
          parameter(ndatemax=30000,npmax=50)
          integer  ni,ns
          double precision  val,mmsp,x
          double precision zi(-2:npmax)
          double precision date(ndatemax)


          if(x.eq.zi(ni))then
             if(ni.le.ns-3)then
                val = 0.d0
             else
                if(ni.le.ns-2)then
           val = ((zi(ni)-zi(ni-1))*mmsp(x,ni,ns,zi,date))*0.25d0
                else
                   if (ni.eq.ns-1)then
                  val = ((zi(ni)-zi(ni-2))*mmsp(x,ni,ns,zi,date)+
     &        (zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+1,zi,date))*0.25d0
                   else
                      if(ni.eq.ns)then
                  val = ((zi(ni)-zi(ni-3))*mmsp(x,ni,ns,zi,date)+
     &                (zi(ni+2)-zi(ni-2))*mmsp(x,ni,ns+1,zi,date)
     &       +(zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+2,zi,date))*0.25d0
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
                 val = (x-zi(ni))*mmsp(x,ni,ns,zi,date)*0.25d0
             else  
             if(ni.eq.ns-2)then
                   val = ((x-zi(ni-1))*mmsp(x,ni,ns,zi,date)+
     &       (zi(ni+4)-zi(ni))*mmsp(x,ni,ns+1,zi,date))*0.25d0
             else   
                if (ni.eq.ns-1)then
                   val =((x-zi(ni-2))*mmsp(x,ni,ns,zi,date)+
     &             (zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+1,zi,date)
     &       +(zi(ni+4)-zi(ni))*mmsp(x,ni,ns+2,zi,date))*0.25d0
                else
                   if(ni.eq.ns)then
                      val =((x-zi(ni-3))*mmsp(x,ni,ns,zi,date)+
     &             (zi(ni+2)-zi(ni-2))*mmsp(x,ni,ns+1,zi,date)
     &             +(zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+2,zi,date)
     &        +(zi(ni+4)-zi(ni))*mmsp(x,ni,ns+3,zi,date))*0.25d0
                   else
                      val = 1.d0
                   endif
                endif
             endif
             endif
          endif 
          endif
             isp = val
             return
             end
c==========================  MMSP   ==================================
      double precision function mmsp(x,ni,ns,zi,date)
      parameter(ndatemax=30000,npmax=50)
      integer  ni,ns
      double precision  val,x
      double precision zi(-2:npmax)
      double precision date(ndatemax)

          if(ni.lt.ns-3)then
             val = 0.d0
          else
             if(ns-3.eq.ni)then
                if(x.eq.zi(ni))then
                   val = 0.d0
                else  
                   val = (4.d0*(x-zi(ni))*(x-zi(ni))
     &             *(x-zi(ni)))/((zi(ni+4)-zi(ni))*(zi(ni+3)
     &             -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
                endif
             else 
                if(ns-2.eq.ni)then
                   if(x.eq.zi(ni))then
                      val = (4.d0*(zi(ni)-zi(ni-1))*(zi(ni)-zi(ni-1)))
     &              /((zi(ni+3)-zi(ni-1))*(zi(ni+2)-zi(ni-1))
     &              *(zi(ni+1)-zi(ni-1)))
                   else  
                      val = (4.d0*(x-zi(ni-1))*(x-zi(ni-1))
     &              *(zi(ni+1)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+2)
     &              -zi(ni-1))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)-zi(ni)))
     &              +   (4.d0*(x-zi(ni-1))*(x-zi(ni))
     &              *(zi(ni+2)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+2)
     &              -zi(ni))*(zi(ni+1)-zi(ni))*(zi(ni+2)-zi(ni-1))) 
     &              +   (4.d0*(x-zi(ni))*(x-zi(ni))
     &              *(zi(ni+3)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+3)
     &              -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
                   endif
                else   
                   if (ns-1.eq.ni)then
                      if(x.eq.zi(ni))then
                         val = (4.d0*((zi(ni)-zi(ni-2))*(zi(ni+1)
     &                  -zi(ni)))/((zi(ni+2)-zi(ni-2))*(zi(ni+1)
     &                  -zi(ni-1))*(zi(ni+1)-zi(ni-2))))
     &                 +((4.d0*((zi(ni)-zi(ni-1))*(zi(ni+2)-zi(ni))) 
     &                 /((zi(ni+2)-zi(ni-2))*(zi(ni+2)-zi(ni-1))
     &                 *(zi(ni+1)-zi(ni-1)))))
                      else
                        val = (4.d0*((x-zi(ni-2))*(zi(ni+1)
     &                  -x)*(zi(ni+1)-x))/((zi(ni+2)
     &                  -zi(ni-2))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)-
     &                   zi(ni))*(zi(ni+1)-zi(ni-2))))
     &                 +((4.d0*((x-zi(ni-1))*(zi(ni+2)-x) 
     &                 *(zi(ni+1)-x))/((zi(ni+2)-zi(ni-2))
     &                 *(zi(ni+2)-zi(ni-1))*(zi(ni+1)-zi(ni-1))*
     &                 (zi(ni+1)-zi(ni)))))
     &                 +((4.d0*((zi(ni+2)-x)*(zi(ni+2)-x) 
     &                 *(x-zi(ni)))/((zi(ni+2)-zi(ni-2))
     &                 *(zi(ni+2)-zi(ni))*(zi(ni+2)-zi(ni-1))*
     &                 (zi(ni+1)-zi(ni)))))
                      endif 
                   else
                      if(ni.eq.ns)then
                         if(x.eq.zi(ni))then
                            val =(4.d0*(x-zi(ni+1))*(x
     &                    -zi(ni+1))/((zi(ni+1)-zi(ni-1))*(zi(ni+1)
     &                    -zi(ni-2))*(zi(ni+1)-zi(ni-3))))
                         else   
                           val =(4.d0*(x-zi(ni+1))*(x
     &                      -zi(ni+1))*(zi(ni+1)-x)/((zi(ni+1)
     &                      -zi(ni-1))*(zi(ni+1)-zi(ni-2))*(zi(ni+1)
     &                      -zi(ni))*(zi(ni+1)-zi(ni-3))))
                         endif
                      else
                         val = 0.d0
                      endif
                   endif
                endif
             endif
          endif

             mmsp = val
             return
             end



c========================          MNBRAK         ===================
      subroutine mnbrak(ax,bx,cx,fa,fb,fc,b,n)
      parameter(npmax=50)
         double precision ax,bx,cx,fa,fb,fc,aux,res
         double precision b(npmax),y(npmax,npmax)
         double precision estimv,gold,glimit,tiny
         parameter (gold=1.618034d0,glimit=100.d0,tiny=1.d-20)
         double precision dum,fu,q,r,u,ulim
         integer n,ni

c     write(*,*)' DEBUT mnbrak',fa,fb
         fa = estimv(ax,n,b,y,aux,ni,res)
ccc   write(*,*)'fa mnbrak',fa
c     stop
         fb = estimv(bx,n,b,y,aux,ni,res)
c         write(*,*)'mnbrak',fa,fb
c         stop

         if(fb.gt.fa)then
            dum = ax
            ax = bx
            bx = dum
            dum = fb
            fb = fa
            fa = dum
         endif
         cx = bx + gold*(bx-ax)
c         write(*,*)'ok yes mnbrak'
c         stop
         fc = estimv(cx,n,b,y,aux,ni,res)
c         write(*,*)'ok mnbrak'
c         stop
 1       if(fb.ge.fc)then
            r = (bx-ax)*(fb-fc)
            q = (bx-cx)*(fb-fa)
            u = bx-((bx-cx)*q-(bx-ax)*r)/
     &       (2.d0*sign(max(abs(q-r),tiny),q-r))
            ulim = bx + glimit*(cx-bx)
            if((bx-u)*(u-cx).gt.0.d0)then
               fu = estimv(u,n,b,y,aux,ni,res)
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
               fu = estimv(u,n,b,y,aux,ni,res)
            else
               if((cx-u)*(u-ulim).gt.0.d0)then
                  fu = estimv(u,n,b,y,aux,ni,res)
                  if(fu.lt.fc)then
                     bx = cx
                     cx = u
                     u = cx + gold*(cx-bx)
                     fb = fc
                     fc = fu
                     fu = estimv(u,n,b,y,aux,ni,res)
                  endif  
               else
                  if((u-ulim)*(ulim-cx).ge.0.d0)then
                     u = ulim
                     fu = estimv(u,n,b,y,aux,ni,res)
                  else
                     u = cx + gold*(cx-bx)
                     fu = estimv(u,n,b,y,aux,ni,res)
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
c         write(*,*)'fin'
c         stop
         return 
         end

c========================      GOLDEN   =========================
      double precision function golden(ax,bx,cx,tol,xmin,n,b,y,aux)

      parameter(ndatemax=30000,npmax=50)
      
      double precision y(npmax,npmax)
      double precision ax,bx,cx,tol,xmin,b(npmax)
      double precision r,c,aux,res
      parameter (r=0.61803399d0,c=1.d0-r)
      double precision f1,f2,x0,x1,x2,x3,estimv
      integer n,ni
      
         x0 = ax
         x3 = cx
         if(abs(cx-bx).gt.abs(bx-ax))then
            x1 = bx
            x2 = bx + c*(cx-bx)
         else
            x2 = bx
            x1 = bx - c*(bx-ax)
         endif
c         write(*,*)'DANS golden',x2,x1,n,aux,ni,res
         f1 = estimv(x1,n,b,y,aux,ni,res)
         f2 = estimv(x2,n,b,y,aux,ni,res)
         
c         write(*,*)'f2 f1',f2,f1
c         write(*,*)'abs',(x3-x0),tol*(abs(x1)+abs(x2))
c         stop
 1       if(abs(x3-x0).gt.tol*(abs(x1)+abs(x2)))then
            if(f2.lt.f1)then
               x0 = x1
               x1 = x2
               x2 = r*x1 + c*x3
               f1 = f2
c         write(*,*)'f2',f1
c         stop
               f2 = estimv(x2,n,b,y,aux,ni,res)
               
c         write(*,*)'f2 DANS golden',f2,f1
c         stop
            else
               x3 = x2
               x2 = x1
               x1 = r*x2+c*x0
               f2 = f1
c         write(*,*)'f2 else',f1
               f1 = estimv(x1,n,b,y,aux,ni,res)
c         write(*,*)'f1 else',f1
c         stop
            endif
            go to 1
          endif
          if(f1.lt.f2)then
             golden = f1
             xmin = x1
          else
             golden = f2
             xmin = x2
          endif
          return
          end


c========================          ESTIMV         ===================

      double precision function estimv(k00,n,b,y,aux,ni,res)

      parameter(ndatemax=30000,npmax=50,nsujetmax=15000)

      double precision v((npmax*(npmax+3)/2)),y(npmax,npmax)
      double precision res,k0(2),k00,ut(ndatemax),bh(npmax),som,h1
      double precision b(npmax),the(-2:npmax),aux,dut(ndatemax)
      integer n,ij,i,k,j,vj
      integer ier,istop,ni


c*****dace2
      double precision t0(nsujetmax),t1(nsujetmax)
      double precision ncsrl
      integer c(nsujetmax)
      integer nt0(nsujetmax),nt1(nsujetmax)
      integer  nsujet,nsim,nva,ndate,nst,maxit
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,
     &              nst,maxit
c*****dace1 
      double precision date(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,zi
      
c*****dace3
      double precision  pe
      integer  effet,nz1,nz2
      common /dace3/pe,effet,nz1,nz2
c*****mem1
      double precision mm3(ndatemax),mm2(ndatemax)
      double precision mm1(ndatemax),mm(ndatemax)
      common /mem1/mm3,mm2,mm1,mm
      
c*****mem2
      double precision im3(ndatemax),im2(ndatemax)
      double precision im1(ndatemax),im(ndatemax)
      common /mem2/im3,im2,im1,im
      
c************
      
      k0(1) = k00*k00
      k0(2) = 0.d0
 
c      write(*,*)'DEBUT ESTIMV',effetcross ,k0
c      stop
      call marq98(k0,b,n,ni,v,res,ier,istop)
      
      
      if(k0(1).gt.0.d0)then
         do 4 ij=1,n
            the(ij-3)=(b(ij))*(b(ij))
            bh(ij) = (b(ij))*(b(ij))
 4       continue
         
         vj = 0
         som = 0.d0
         dut(1) = (the(-2)*4.d0/(zi(2)-zi(1)))
         ut(1) = the(-2)*dut(1)*0.25d0*(zi(1)-zi(-2))
         do 8 i=2,ndate-1
            do 6 k = 2,n-2
               if ((date(i).ge.zi(k-1)).and.(date(i).lt.zi(k)))then
                  j = k-1
                  if ((j.gt.1).and.(j.gt.vj))then
                     som = som+the(j-4)
                     vj  = j
                  endif   
               endif
 6          continue 
            ut(i) = som +(the(j-3)*im3(i))+(the(j-2)*im2(i))
     &           +(the(j-1)*im1(i))+(the(j)*im(i))
            dut(i) = (the(j-3)*mm3(i))+(the(j-2)*mm2(i))
     &           +(the(j-1)*mm1(i))+(the(j)*mm(i))
 8       continue
         i = n-2
         h1 = (zi(i)-zi(i-1))
         ut(ndate) = som+ the(i-4) + the(i-3)+the(i-2)+the(i-1)
         dut(ndate) = (4.d0*the(i-1)/h1)
         
         call test(bh,ut,dut,k0,n,aux,v,y)
         estimv = - ((res-pe)) - aux
c         write(*,*)'estimv',k0,-aux,ni,estimv
      else
         aux = -n
c     write(*,*)'estimv2',k0,-aux,ni
      endif
      
      return
      end
      
c=================calcul de la hessienne  et de omega  ==============
      subroutine test(b,ut,dut,k0,n,res,v,y)
      parameter(ndatemax=30000,npmax=50,nsujetmax=15000)

      double precision hessh(npmax,npmax),hess(npmax,npmax)
      double precision omeg(npmax,npmax),y(npmax,npmax)
      integer n,i,j,np,indx(npmax)
      double precision b(npmax),k0(2),d,dut(ndatemax),ut(ndatemax)
      double precision res,tra ,v((npmax*(npmax+3)/2))

c*****dace1 
      double precision date(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,zi
c*****dace2
      double precision t0(nsujetmax),t1(nsujetmax)
      double precision ncsrl
      integer c(nsujetmax)
      integer nt0(nsujetmax),nt1(nsujetmax)
      integer  nsujet,nsim,nva,ndate,nst,maxit
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,
     &              nst,maxit
c*****
         do 10 i = 1,n
            do 5 j = 1,n
               hess(i,j) = 0.d0 
 5         continue
 10      continue
   
 
         do 20 i = 1,n
            do 15 j = i,n
               call mat(hess(i,j),ut,dut,i,j,n)
c               write(*,*)'hess test',hess(i,j),i,j
 15         continue
 20      continue
        do 40 i = 2,n
            do 35 j = 1,i-1
               hess(i,j)=hess(j,i)
 35         continue
 40      continue


         call calcomeg(n,omeg)

         do 100 i = 1,n
            do 90 j = 1,n
               hessh(i,j)=-hess(i,j)
               hess(i,j) = hess(i,j) - (2.d0*k0(1)*omeg(i,j)) 
 90         continue   
 100     continue

         np = n
         do 112 i=1,n
            do 111 j=1,n
               y(i,j)=0.d0
 111        continue
            y(i,i)=1.d0
 112     continue
         call ludcmp(hess,n,np,indx,d)

         do 113 j=1,n
            call lubksb(hess,n,np,indx,y(1,j))
 113     continue

c         write(*,*)'hess test',y(1,1)

         tra = 0.d0
         do 150 i=1,n
            do 140 j=1,n
               tra = tra + y(i,j)*hessh(j,i)
c         write(*,*)'test tra',tra,y(i,j),hessh(j,i),i,j
 140        continue
 150     continue
   
c         write(*,*)'test',k0,-tra
         res = (tra)

         end


c=======================  CALOMEG  ===========================
      subroutine calcomeg(n,omeg)
c        remplissage de la matrice omega n*n
c          elle a 7 diagonales
      parameter(ndatemax=30000,npmax=50)

      double precision omeg(npmax,npmax)
      integer n
      double precision calc00,calc01,calc02
      integer i,j
c*****dace1 
      double precision date(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,zi
c*****pen1
      double precision  m3m3(npmax),m2m2(npmax),m1m1(npmax)
      double precision  mmm(npmax),m3m2(npmax)
      common /pen1/m3m3,m2m2,m1m1,mmm,m3m2
c*****pen2
      double precision m3m1(npmax),m3m(npmax),m2m1(npmax)
      double precision m2m(npmax),m1m(npmax)
      common /pen2/m3m1,m3m,m2m1,m2m,m1m
c**************************
      do 5 i=1,n
         do 3 j=1,n
            omeg(i,j)=0.d0
 3       continue
 5    continue
   
      omeg(1,1)=calc00(1,n)
      omeg(1,2)=calc01(1,n)
      omeg(1,3)=calc02(1,n)
      omeg(1,4)=m3m(1)
      omeg(2,1)=omeg(1,2)
      omeg(2,2)=calc00(2,n)
      omeg(2,3)=calc01(2,n)
      omeg(2,4)=calc02(2,n)
      omeg(2,5)=m3m(2)
      omeg(3,1)=omeg(1,3)
      omeg(3,2)=omeg(2,3)
      omeg(3,3)=calc00(3,n)
      omeg(3,4)=calc01(3,n)
      omeg(3,5)=calc02(3,n)
      omeg(3,6)=m3m(3)
      do 10 i=4,n-3
         omeg(i,i-3)=omeg(i-3,i)
         omeg(i,i-2)=omeg(i-2,i)
         omeg(i,i-1)=omeg(i-1,i)
         omeg(i,i)=calc00(i,n)
         omeg(i,i+1)=calc01(i,n)
         omeg(i,i+2)=calc02(i,n)
         omeg(i,i+3)=m3m(i)
 10   continue   
      omeg(n-2,n-5)=omeg(n-5,n-2)
      omeg(n-2,n-4)=omeg(n-4,n-2)
      omeg(n-2,n-3)=omeg(n-3,n-2)
      omeg(n-2,n-2)=calc00(n-2,n)
      omeg(n-2,n-1)=calc01(n-2,n)
      omeg(n-2,n)=calc02(n-2,n)
      omeg(n-1,n-4)=omeg(n-4,n-1)
      omeg(n-1,n-3)=omeg(n-3,n-1)
      omeg(n-1,n-2)=omeg(n-2,n-1)
      omeg(n-1,n-1)=calc00(n-1,n)
      omeg(n-1,n)=calc01(n-1,n)
      omeg(n,n-3)=omeg(n-3,n)
      omeg(n,n-2)=omeg(n-2,n)
      omeg(n,n-1)=omeg(n-1,n)
      omeg(n,n)=calc00(n,n)

      end


c====================  MAT  ==================================
      subroutine mat(res,ut,dut,k,l,n)
      parameter(ndatemax=30000,npmax=50,nsujetmax=15000)

      double precision res,dut(ndatemax),ut(ndatemax)
      integer k,l,j,ni,n,ni1
      double precision res1,msp,aux2,sp
      double precision u2
      integer i

c*****dace1 
      double precision date(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,zi
c*****dace2
      double precision t0(nsujetmax),t1(nsujetmax)
      double precision  ncsrl
      integer c(nsujetmax)
      integer nt0(nsujetmax),nt1(nsujetmax)
      integer  nsujet,nsim,nva,ndate,nst,maxit
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,
     &              nst,maxit
c**

             
c---------- calcul de la hessienne ij ------------------
          res = 0.d0
          res1 = 0.d0
          do 10 i=1,nsujet
             if(c(i).eq.1)then  !event
                u2 = dut(nt1(i)) 
c                write(*,*)'mat u2',u2,nt1(i),i
                do 6 j = 2,n-2
                   if((date(nt1(i)).ge.zi(j-1)).and.
     &                  (date(nt1(i)).lt.zi(j)))then
                      ni = j-1
c                     write(*,*)'mat ni',ni
                   endif
 6              continue 
                if(date(nt1(i)).eq.zi(n-2))then
                   ni = n-2
                endif   
c-------attention numero spline 
                aux2 = msp(nt1(i),ni,k)*msp(nt1(i),ni,l)
c                write(*,*)'mat aux2',aux2,nt1(i),ni,k
                if (u2.le.0.d0)then
                   res1 = 0.d0
                else   
                   res1 = - aux2/(u2*u2)
                endif  
c                write(*,*)'mat',res1
c                stop
             else !censure  
                res1 = 0.d0
             endif 
          res = res + res1
 10    continue   
       
       end

c==========================  MSP   ==================================
          double precision function msp(i,ni,ns)
          parameter(ndatemax=30000,npmax=50)
             integer ni,ns,i
             double precision val
c*****dace1 
      double precision date(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,zi
c************************************

          if(ni.lt.ns-3)then
             val = 0.d0
          else
             if(ns-3.eq.ni)then
                if(date(i).eq.zi(ni))then
                   val = 0.d0
                else  
                   val = (4.d0*(date(i)-zi(ni))*(date(i)-zi(ni))
     &             *(date(i)-zi(ni)))/((zi(ni+4)-zi(ni))*(zi(ni+3)
     &             -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
                endif
             else 
                if(ns-2.eq.ni)then
                   if(date(i).eq.zi(ni))then
                      val = (4.d0*(zi(ni)-zi(ni-1))*(zi(ni)-zi(ni-1)))
     &              /((zi(ni+3)-zi(ni-1))*(zi(ni+2)-zi(ni-1))
     &              *(zi(ni+1)-zi(ni-1)))
                   else  
                      val = (4.d0*(date(i)-zi(ni-1))*(date(i)-zi(ni-1))
     &              *(zi(ni+1)-date(i)))/((zi(ni+3)-zi(ni-1))*(zi(ni+2)
     &              -zi(ni-1))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)-zi(ni)))
     &              +   (4.d0*(date(i)-zi(ni-1))*(date(i)-zi(ni))
     &              *(zi(ni+2)-date(i)))/((zi(ni+3)-zi(ni-1))*(zi(ni+2)
     &              -zi(ni))*(zi(ni+1)-zi(ni))*(zi(ni+2)-zi(ni-1))) 
     &              +   (4.d0*(date(i)-zi(ni))*(date(i)-zi(ni))
     &              *(zi(ni+3)-date(i)))/((zi(ni+3)-zi(ni-1))*(zi(ni+3)
     &              -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
                   endif
                else   
                   if (ns-1.eq.ni)then
                      if(date(i).eq.zi(ni))then
                         val = (4.d0*((zi(ni)-zi(ni-2))*(zi(ni+1)
     &                  -zi(ni)))/((zi(ni+2)-zi(ni-2))*(zi(ni+1)
     &                  -zi(ni-1))*(zi(ni+1)-zi(ni-2))))
     &                 +((4.d0*((zi(ni)-zi(ni-1))*(zi(ni+2)-zi(ni))) 
     &                 /((zi(ni+2)-zi(ni-2))*(zi(ni+2)-zi(ni-1))
     &                 *(zi(ni+1)-zi(ni-1)))))
                      else
                        val = (4.d0*((date(i)-zi(ni-2))*(zi(ni+1)
     &                  -date(i))*(zi(ni+1)-date(i)))/((zi(ni+2)
     &                  -zi(ni-2))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)-
     &                   zi(ni))*(zi(ni+1)-zi(ni-2))))
     &                 +((4.d0*((date(i)-zi(ni-1))*(zi(ni+2)-date(i)) 
     &                 *(zi(ni+1)-date(i)))/((zi(ni+2)-zi(ni-2))
     &                 *(zi(ni+2)-zi(ni-1))*(zi(ni+1)-zi(ni-1))*
     &                 (zi(ni+1)-zi(ni)))))
     &                 +((4.d0*((zi(ni+2)-date(i))*(zi(ni+2)-date(i)) 
     &                 *(date(i)-zi(ni)))/((zi(ni+2)-zi(ni-2))
     &                 *(zi(ni+2)-zi(ni))*(zi(ni+2)-zi(ni-1))*
     &                 (zi(ni+1)-zi(ni)))))
                      endif 
                   else
                      if(ni.eq.ns)then
                         if(date(i).eq.zi(ni))then
                            val =(4.d0*(date(i)-zi(ni+1))*(date(i)
     &                    -zi(ni+1))/((zi(ni+1)-zi(ni-1))*(zi(ni+1)
     &                    -zi(ni-2))*(zi(ni+1)-zi(ni-3))))
                         else   
                           val =(4.d0*(date(i)-zi(ni+1))*(date(i)
     &                      -zi(ni+1))*(zi(ni+1)-date(i))/((zi(ni+1)
     &                      -zi(ni-1))*(zi(ni+1)-zi(ni-2))*(zi(ni+1)
     &                      -zi(ni))*(zi(ni+1)-zi(ni-3))))
                         endif
                      else
                         val = 0.d0
                      endif
                   endif
                endif
             endif
          endif

             msp = val
             return
             end
c==========================   SP   ==================================
          double precision function sp(i,ni,ns)
          parameter(ndatemax=30000,npmax=50)
          
             integer ni,ns,i
             double precision val,msp
c*****dace1 
      double precision date(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,zi
c-----------

          if(date(i).eq.zi(ni))then
             if(ni.le.ns-3)then
                val = 0.d0
             else
                if(ni.le.ns-2)then
                   val = ((zi(ni)-zi(ni-1))*msp(i,ni,ns))*0.25d0
                else
                   if (ni.eq.ns-1)then
                      val = ((zi(ni)-zi(ni-2))*msp(i,ni,ns)+
     &                (zi(ni+3)-zi(ni-1))*msp(i,ni,ns+1))*0.25d0
                   else
                      if(ni.eq.ns)then
                         val = ((zi(ni)-zi(ni-3))*msp(i,ni,ns)+
     &                       (zi(ni+2)-zi(ni-2))*msp(i,ni,ns+1)
     &              +(zi(ni+3)-zi(ni-1))*msp(i,ni,ns+2))*0.25d0
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
                   val = (date(i)-zi(ni))*msp(i,ni,ns)*0.25d0
             else  
             if(ni.eq.ns-2)then
                   val = ((date(i)-zi(ni-1))*msp(i,ni,ns)+
     &             (zi(ni+4)-zi(ni))*msp(i,ni,ns+1))*0.25d0
             else   
                if (ni.eq.ns-1)then
                   val =((date(i)-zi(ni-2))*msp(i,ni,ns)+
     &             (zi(ni+3)-zi(ni-1))*msp(i,ni,ns+1)
     &             +(zi(ni+4)-zi(ni))*msp(i,ni,ns+2))*0.25d0
                else
                   if(ni.eq.ns)then
                      val =((date(i)-zi(ni-3))*msp(i,ni,ns)+
     &             (zi(ni+2)-zi(ni-2))*msp(i,ni,ns+1)
     &             +(zi(ni+3)-zi(ni-1))*msp(i,ni,ns+2)
     &             +(zi(ni+4)-zi(ni))*msp(i,ni,ns+3))*0.25d0
                   else
                      val = 1.d0
                   endif
                endif
             endif
             endif
          endif 
          endif
             sp = val
             return
             end
c================

c=========================  CALC00  =========================
          double precision function calc00(j,n) 
          parameter(npmax=50)
          double precision part
          integer j,n
c*****pen1
          double precision  m3m3(npmax),m2m2(npmax),m1m1(npmax)
          double precision  mmm(npmax),m3m2(npmax)
          common /pen1/m3m3,m2m2,m1m1,mmm,m3m2
c*****pen2
          double precision m3m1(npmax),m3m(npmax),m2m1(npmax)
          double precision m2m(npmax),m1m(npmax)
          common /pen2/m3m1,m3m,m2m1,m2m,m1m
c----------------------------
c        entre i et i+1 ---> m*m et       i = j-3

c       entre i+1 et i+2 ---> m1*m1 et      i = j-2

c       entre i+2 et i+3 ---> m2*m2 et    i = j-1

c       entre i+3 et i+4 ---> m3*m3  et    i = j

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

             calc00 = part
             return
             end
c=========================  CALC01  =========================
      double precision function calc01(j,n)
      parameter(npmax=50)
      double precision part
      integer j,n

c*****pen1
          double precision  m3m3(npmax),m2m2(npmax),m1m1(npmax)
          double precision  mmm(npmax),m3m2(npmax)
          common /pen1/m3m3,m2m2,m1m1,mmm,m3m2
c*****pen2
          double precision m3m1(npmax),m3m(npmax),m2m1(npmax)
          double precision m2m(npmax),m1m(npmax)
          common /pen2/m3m1,m3m,m2m1,m2m,m1m
c----------------------------
c        entre j+1 et j+2 ---> m1*m et       i = j-2

c        entre j+2 et j+3 ---> m2*m1 et    i = j -1

c        entre j+3 et j+4 ---> m3*m2 et     i = j 


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

             calc01 = part
             return
             end
c=========================  CALC02  =========================
          double precision function calc02(j,n)
          parameter(npmax=50)

          double precision part
          integer j,n

c*****pen1
          double precision  m3m3(npmax),m2m2(npmax),m1m1(npmax)
          double precision  mmm(npmax),m3m2(npmax)
          common /pen1/m3m3,m2m2,m1m1,mmm,m3m2
c*****pen2
          double precision m3m1(npmax),m3m(npmax),m2m1(npmax)
          double precision m2m(npmax),m1m(npmax)
          common /pen2/m3m1,m3m,m2m1,m2m,m1m
c====================
c        entre j+3 et j+4 ---> m2*m et     i = j-1
 
c        entre j+2 et j+3 ---> m3*m1 et       i = j


             if(j.eq.1)then
                part = m3m1(j)
             else   
                if(j.ne.n-2)then
                   part = m3m1(j) + m2m(j-1) 
                else
                   part = m2m(j-1)
                endif
             endif   

             calc02 = part
             return
             end
