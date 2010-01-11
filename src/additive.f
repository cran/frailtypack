c==ici : cov , le terme de covariance est estimé directement
c=== sans aucune contraintes
c== au lieu de rho 
c====================================================
c 13 fevrier 2007
c===================================================
c idem à progcorrel_LAPLACE.f mais sans contraintes sur rho
c car dans progcorrel_LAPALCE.f on a souvent un probleme 
c d inversion de la mat3(rice des derivées secondes
c
c====================================================
c 26 avril 2006
c===================================================
c avec en plus un terme de correlation non nul = rho
c entre intercept et pente aleatoire
c====================================================
c 16-mars 2006
c===================================================
c PROGRAMME POUR intercept et pente aléatoire
c dans un modèle de survie
c application : méta analyse =
c effet "essai" (géographique) en effet simple puis en interaction avec le traitement
c
c G (nb d'essais) intégration doubles à calculer => GENZ = quadrature gaussienne adaptative
c vraisemblance pénalisée
c parametres à estimer =
c cov= covariance entre(intercept,pente aleatoire) =rho*sigma*tau
c rho= correlation(intercept,pente aleatoire) (1)
c sigma2 = variance pour interecept aleatoire (2)
c tau2 = variance pour interaction avec TTT, pente aleatoire (3)
c beta = variables explicatives (4....)
c
c====================================================
c====================================================
c modifie en avril 2004 pour etudier des fichiers de donnees
c et non des donnees simulees
c attention a trier le fichier de donnees par groupe et sous groupe !!!
c attention : integration sur [0;5.]
c=====================================================
c SIMULATIONS sur plusieurs fichiers generes
c a partir du programme : nested_essai3.f (1 seule simulation)
c     avril /2004
c=========================================================
c pb sur VET dans funcpa corrigé le 28/01/2004 
c===========================================
c ATTENTION dans le fichier de données, 
c les données doivent être classées par sous groupe et groupe.
c et les données doivent etre nested :
c  important au moment de l attribution d un numero de groupe
c**************************************************************
c prog test ou on reprend le funcpa précédent des shared frailty.

c integration par gauss-laguerre
c====================================================================
C 17/04/2003
c NESTED FRAILTY MODEL
c DEUX EFFETS ALÉATOIRES EMBOITÉS
c EXEMPLE : UN SUR LES VILLES ET UN PAR ZONE
c INTEGRATION PAR MONTE CARLO
C (2 INTEGRALES DANS LA LOG VRAISEMBLANCE MARGINALE avec troncature)
c ON NE CONSIDERE QUE LE CAS NESTED AVEC 2 FRAILTY ET VAR EXPLICATIVES
c********************************************************************
c ancien nouveau programme ...
c********************************************************************
c bis le 12/05/2003, modification de funcpa, erreur sur ddl d ordre 3
c********************************************************************
c bis le 7/11/2001, modification de funcpa, erreur sur ddl d ordre 3
c********************************************************************
c programme modifie le 01juin2001, car modification de la vraisemblance
c en tenant compte de la loi de la frailty conditionnelle a la survie

c********************************************************************
c  puis, modification de marquard (erreur) le 8juin2001
c     car dans dchol le vecteur v() est modifie, boucle 4001 rajoutee 
c  lorsque la mat3(rice des derivees 2 n'est pas inversible,
c  il faut reinitialiser completement la mat3(rice FU et
c  pas seulement modifier la diagonale car FU a ete modifie dans Dchole.
c********************************************************************

c BANDES DE CONFIANCES 
c bayesiennes (H-1) ou non baysienne H-1.I.H-1

c TRACE = nb de degres de liberte pour les fonctions de risque
c  a partir du kappa fixe on peut en deduire le mdf

c delta method poour la variance de theta= b(1)**2
c TEST
c de Wald en utilisant la variance corrigee ou non 
c avec test du LRS penalise ?

c dans chaque simulation, on ninitialise pas les parametres avec valeurs precedentes
CCC-----stratification------CCCCCCCCCCCCCCCCCCCCCCCC
CCc sans troncature

CCCC vraisemblance calculee avec ou sans effet aleatoire

CCCC programme qui permet des donnees par strates (2 uniquement)
CCCC la fonction de risque de base est donc différente
CCCC  pour chaque strate (2 fonctions differentes)
CCCC  estimees par 2 bases de splines differentes

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c    estimat3(ion des coefficient de regression  ... spline d'ordre 4
c    avec censure a droite, troncature a gauche, noeuds aux choix


c     avec fraitly avec loi Gamma

       subroutine additive(nsujetAux,ngexactAux,icenAux,nstAux,
     &     nzAux,ax1,ax2,tt0Aux,tt1Aux,icAux,groupeAux,nvaAux,
     &     strAux,vaxAux,interaction,AGAux,noVar,maxitAux,irep1,
     &     correlAux,
     &     b,coef,varcoef,varcoef2,rhoEnd,covEnd,varcovEnd,
     &     varSigma2,varTau2,ni,res,k0,x1Out,lamOut,suOut,
     &     x2Out,lam2Out,su2Out)

     
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************


c  Added JRG August '06
      integer nsujetAux,ngexactAux,icenAux,nstAux, interaction
      integer nzAux,nvaAux,AGaux,noVar,maxitAux,correlAux 
      double precision ax1,ax2
      double precision tt0Aux(nsujetAux),tt1Aux(nsujetAux)
      integer icAux(nsujetAux),groupeAux(nsujetAux)
      integer strAux(nsujetAux)
      double precision vaxAux(nsujetAux,nvaAux)
      double precision varSigma2(2),varTau2(2)
      double precision rhoEnd,covEnd,varcovEnd
      double precision coef(nva),varcoef(nva),varcoef2(nva)

      double precision  x1Out(99),lamOut(99,3),suOut(99,3),
     &                  x2Out(99),lam2Out(99,3),su2Out(99,3)


      integer :: isujet,AG
      integer :: groupe
      integer :: j,k,nz,n,np,cpt,ii,iii,iiii,ver
      integer :: cptstr1,cptstr2
      integer :: i,ic,ni,ier,istop,ef,l,str
      integer :: cptni,cptni1,cptni2
      integer :: nb_echec,nb_echecor
      integer :: nis,m,idum,id
      integer :: auxng  ,indic_codage 
      integer :: ig

c .............. yas déb        
	   integer , dimension(nvarmax) :: filtre,filtre2 
c filtre = dans le modele ou pas ? filtre2 = en  interaction ou pas
  
      integer :: stracross(nsujetmax) !pour crossvalidation
      integer :: irep1,nvat,nvacross,nstcross,effetcross !pour crossvalidation
c .............. yas fin     	  

c yas     integer , dimension(nvarmax) :: filtre,filtre2 
c filtre = dans le modele ou pas ? filtre2 = en  interaction ou pas 
c yas      integer irep1,nvat
 
      real , dimension(nvarmax) :: vax
      double precision tt0 
      double precision tt1

      double precision :: h
      double precision :: eca,varsmarg
      double precision :: smoy
      double precision :: wres
      double precision :: xmin1,xmin2,res,min,max,maxt
      double precision :: bi,bs,wald
      double precision :: lrs,trace,trace1,trace2
      double precision :: uniran
      double precision :: moyvar1,moyvar2
      double precision :: moyse1,moyse2, moyse_cor1,moyse_cor2
      double precision moy_peh0,moy_peh1
      double precision BIAIS_moy
      double precision auxi,ax,bx,cx,tol,ddl
      double precision fa,fb,fc,golden3,estimv3
      double precision f1,f2,f3,varcov
      
      
      double precision , dimension(2*nsujetmax)::aux
      double precision , dimension(npmax,npmax):: y
      double precision , dimension((npmax*(npmax+3)/2))::v
      double precision , dimension(2)::k0
      double precision , dimension(npmax)::b
      double precision , dimension(npmax,npmax):: I1_hess,H1_hess
      double precision , dimension(npmax,npmax):: I2_hess,H2_hess
      double precision , dimension(npmax,npmax):: HI1,HI2
      double precision , dimension(npmax,npmax):: HIH,IH,HI
      

c*****************************************************************
      
c*****nmax
      integer :: nmax
      common /nmax/nmax
c*****dace 1 
      double precision , dimension(ndatemax)::date
      double precision , dimension(-2:npmax)::zi
      common /dace1/date,zi
      
c*****dace2
      double precision , dimension(nsujetmax) ::t0,t1
      double precision :: ncsrl
      integer , dimension(nsujetmax) ::c
      integer, dimension(nsujetmax) :: nt0,nt1
      integer :: nsujet,nsim,nva,ndate,nst
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,nst
c*****dace4
      integer , dimension(nsujetmax) :: stra
      common /dace4/stra
c*****ve1
      double precision , dimension(nsujetmax,nvarmax):: ve,ve2
      common /ve1/ve,ve2
c*****dace3
      double precision :: pe
      integer :: effet,nz1,nz2,correl
      common /dace3/pe,effet,nz1,nz2,correl
c***  
      integer correlini
c*****dace7
      double precision I_hess(npmax,npmax),H_hess(npmax,npmax)
      double precision Hspl_hess(npmax,npmax)
      double precision PEN_deri(npmax,1) 
      double precision hess(npmax,npmax)
      common /dace7/PEN_deri,I_hess,H_hess,Hspl_hess,hess
c*****contrib
      integer :: ngexact,ng        !nb EXACT de gpes 
      common /contrib/ngexact,ng
      
c*****groupe
      integer , dimension(nsujetmax) :: g
      integer , dimension(ngmax) :: nig ! nb de sujet par groupe
      common /gpe/g,nig
      
c******indicateur de troncature
      integer :: indictronq     ! =0 si donnees non tronquées
      common /troncature/indictronq
      
c*****mij
      integer , dimension(ngmax) :: mi ! nb de dc dans gpe i
      common /mij/mi
      
c******indicateur du nb de parametres
      integer :: nbpara
      common /nbpara/nbpara

c****** integration
      integer :: NDIM,MINPTS,MAXPTS,npts,type_int
      DOUBLE PRECISION :: EPSABS,EPSREL
      common /integration/NDIM,MINPTS,MAXPTS,npts,EPSABS,EPSREL,type_int
      
c*****inversion
      integer :: indic_sousv,sousm
      common /invers/indic_sousv,sousm
c*******sigma2tau2rho
      double precision :: sigma2,tau2,rho,cov
      common /sigma2tau2/sigma2,tau2,rho,cov
      
c******u_tildefinal
      double precision , dimension(ngmax) :: u_tildefinal,v_tildefinal
c      common /utildefinal/u_tildefinal,v_tildefinal
c******************
c********************************************************************
       
      character (len=20):: nomvarl,nominter
      character (len=20),dimension(25):: nomvar
      character (len=20):: donnees
      character (len=20)::fich1
      character (len=20)::fich2
      character (len=20)::fich3
      character (len=20)::fich4
      character (len=20)::fich1b
      character (len=20)::fich2b
      character (len=20)::fich3b
      character (len=20)::fich4b
      character (len=20)::fich1c
      character (len=20)::fich2c
      character (len=20)::fich3c
      character (len=20)::fich4c
      character*20 dateamj
      character*20 zone
      character*20 heure1
      character*20 heure2
      integer , dimension(8) :: values
      
      
c     nst: nombre de strates
c     ist: appartenance aux strates
c     ib: mat3(rice de variance de beta chapeau
c     I_hess : -hessienne non inversee sur vraisemblance non penalisee
c     H_hess : inverse de -hessienne  sur vraisemblance penalisee
      
c     nig(ngmax) : nb sujet dans groupe (different dans frailty.f !!!)
c     mi(ngmax) : nb deces dans groupe  ! double precision
c     g(nsujetmax) : numero du groupe pour un  sujet donné 
      
c     ngexact : nb exact de groupes
      


c      write(*,*)'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
c      write(*,*)'     PROG.f additive(correl) gaussien frailty model '
c      write(*,*)'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'

c      open(4,file='outcorrel3')

c      write(4,*)'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
c      write(4,*)'     PROG.f additive(correl) gaussien frailty model '
c      write(4,*)'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
c      write(4,*)''

c      call date_and_time(dateamj,heure1,zone,values)
      

      nmax = maxitAux           !nb iterations max dans marquard



c=== initialisations 
        lrs=0.d0    
        nb_echec=0 
        nb_echecor=0
c=== fin initialisations 


    
c	open(2,file='mach2.inf')
c	open(2,file='simul.inf')
c 1900	continue (pour faire plusieurs jeux de simulations)
c        read(2,*)nsujet         !nb de sujets total
c        read(2,*)ng!ngexact
c        write(*,*)'  nsujet',nsujet,ng!ngexact
c        read(2,*)nst
c        read(2,*)correl  
c        correlini=correl
c        read(2,*)ver
c         nva = 0
c         if(ver.gt.0)then 
c            do 44 j=1,ver
c               read(2,*)nomvarl,filtre(j),filtre2(j) 
c               PRINT*,'filtre ',nomvarl,filtre(j)  ,filtre2(j)  
c               nva = nva + filtre(j)
c               if(filtre(j).eq.1)then
c                  nomvar(nva) = nomvarl
c               endif 
c               if(filtre2(j).eq.1)then
c                  nominter = nomvarl
c               endif     
c 44         continue
c         endif
c        read(2,*)donnees
c        open(9,file=donnees)
c        read(2,*)correl !=1 si terme de correlation à estimer

c initialisations pour integration
c        write(4,*)'======================================'
              
c        read(2,*)indic_codage
c        if(indic_codage.eq.1)then 
c           write(4,*)'** codage traitement en + ou - 1/2 **'
c        else
c           write(4,*)'** codage traitement en + ou - 0/1 **'
c        endif

c        read(2,*)AG
c        read(2,*)nz


c Juan  Added August'07

        nsujet=nsujetAux
        ng=ngexactAux
        nssgbyg=nssgbygAux
        icen=icenAux
        nst=nstAux

        xmin1=ax1
        xmin2=ax2   

        if (noVar.eq.1) then 
          do i=1,nvaAux
           filtre(i)=0
          enddo  
          nva=0  
        else
          do i=1,nvaAux
           filtre(i)=1
           filtre2(i)=0
          enddo  
          filtre2(interaction)=1
          nva=nvaAux  
        end if  

        ver=nvaAux
 

        AG=AGAux
        nz=nzAux  

        correl=correlAux
        correlini=correl

c        write(*,*) strAux 


c******************************************
c******************************************
c---debut des iterations des simulations
       	
        effet = 1
        nig = 0
        g=0
        
c        write(4,*)'nb de variable explicative dans le modele:',nva
c        write(*,*)'nb de variable explicative dans le modele:',nva
    
c------------  lecture fichier -----------------------

         maxt = 0.d0
         cpt = 0
         k = 0
         cptstr1 = 0
         cptstr2 = 0

c         write(*,*)'** Fichier de données = ',donnees
c         write(4,*)'** Fichier de données = ',donnees
c         open(9,file=donnees)
         
         do 10 i = 1,nsujet    
c           str=1             
            if(nst.eq.2)then
c Juan         read(9,*)tt0,tt1,groupe,ic,str,(vax(j),j=1,ver)
               tt0=tt0Aux(i)
               tt1=tt1Aux(i)
               ic=icAux(i)
               groupe=groupeAux(i)
               str=strAux(i)
               do j=1,ver
                 vax(j)=vaxAux(i,j)  
               enddo


            else
c Juan        read(9,*)tt0,tt1,groupe,ic,(vax(j),j=1,ver) 
               tt0=tt0Aux(i)
               tt1=tt1Aux(i)
               ic=icAux(i)
               groupe=groupeAux(i)
               str=strAux(i)
               do j=1,ver
                 vax(j)=vaxAux(i,j)  
               enddo
            endif

c           write(*,*)'**',i,tt0,tt1,groupe,ic,(vax(j),j=1,ver) 

            k = k +1
            if(k.eq.1)then
               auxng=groupe
               ngexact=1
               g(k)=1
            else
               g(k)=groupe
            endif
c------------------   observation c=1
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
c     t0(k) = 0.d0pro
                  t1(k) = tt1
                
                  if(auxng.ne.groupe)then !chgt de groupes
                     ngexact=ngexact+1
                     auxng=groupe

                     if(k.ne.1)then
                        g(k)=g(k-1)+1
c             write(*,*)'**  obser groupe **',g(k),g(k-1),ngexact,auxng,k
                     endif
                 goto 100             
                  endif
                  
 100              continue

                  nig(g(k)) = nig(g(k))+1
                  iii = 0
                  iiii = 0
                  do 6 ii = 1,ver
                     if(filtre(ii).eq.1)then
                        iii = iii + 1
                        ve(i,iii) = dble(vax(ii))
                        ve2(i,iii) = ve(i,iii) 
c==================recodage en +-1/2 ===============================
            !si indic_codage=1 ==> recodage en +-1/2, sinon, 0/1
c                        if(indic_codage.eq.1)then
c                           ve(k,iii) = ve(k,iii) -0.5d0
c                           ve2(k,iii) = ve(k,iii)
c                        endif
c====================================================================
                     endif
                     if(filtre2(ii).eq.1)then
                        iiii = iiii + 1
                        ve(k,iiii) = dble(vax(ii))
                        ve2(k,iiii) = ve(k,iiii)  
c==================recodage en +-1/2 =============================== 
c                        if(indic_codage.eq.1)then
c                           ve(k,iiii) =ve(k,iiii) -0.5d0
c                           ve2(k,iiii) = ve(k,iiii) 
c                        endif
                     endif
c====================================================================
 6                continue   
               else 
c------------------   censure a droite  c=0
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
                     do 8 ii = 1,ver
                        if(filtre(ii).eq.1)then
                           iii = iii + 1
                           ve(k,iii) = dble(vax(ii))
                           ve2(k,iii) =  ve(k,iii) 
c==================recodage en +-1/2 ===============================
c                        if(indic_codage.eq.1)then
c                           ve(k,iii) =ve(k,iii) -0.5d0
c                           ve2(k,iii) =  ve(k,iii) 
c                        endif
c====================================================================
                        endif
                        if(filtre2(ii).eq.1)then
                           iiii = iiii + 1
                           ve(k,iiii) = dble(vax(ii))
                           ve2(k,iiii) =  ve(k,iiii) 
c==================recodage en +-1/2 =============================== 
c                           if(indic_codage.eq.1)then
c                              ve(k,iiii) =ve(k,iiii) -0.5d0
c                              ve2(k,iiii) =  ve(k,iiii) 
c                           endif
c====================================================================
                        endif
 8                   continue 

                      t0(k) =  tt0
                      t1(k) = tt1
                      if(auxng.ne.groupe)then !chgt de groupes
                         ngexact=ngexact+1
                         auxng=groupe
                         if(k.ne.1)then
                            g(k)=g(k-1)+1
c                write(*,*)'******  cen groupe **',g(k),ngexact,auxng,k

                         endif
                         goto 101
                      endif
                      
 101                  continue
                      nig(g(k)) = nig(g(k))+1
                   endif
                endif
                if (maxt.lt.t1(k))then
                   maxt = t1(k)
                endif
 10          continue 


             

c--------------------------- fin lecture du fichier
c             write(*,*)'** nig **',(nig(i),i,i=1,ngexact)
c          write(*,*)'**  groupe **',(g(i),i,i=1,nsujet)
c          write(*,*)'**  *********************************'

c         write(*,*)'max',maxt
c         write(*,*)'** Fichier de données = ',donnees
c         write(4,*)'** Fichier de données = ',donnees
c         write(*,*)'** nb de groupes annonce =',ng
c         write(*,*)'** nb exact de groupes =',ngexact
c     write(*,*)'** num de groupe =',(g(i),i=1,nsujet)
         
c         write(4,*)'** nb de groupes annonce =',ng
c         write(4,*)'** nb de groupes exact=',ngexact 
c     write(4,*)'nb individu par groupe et par strate = ',nis
c         write(4,*)'** nb de strates = ',nst
c         write(4,*)'** nb d observations',nsujet
c         write(4,*)'** nb de deces ',cpt
c         write(4,*)'** nb de censures ',(nsujet-cpt)
c         write(*,*)'** nb d observations',nsujet
c         write(*,*)'** nb de deces ',cpt
c         write(*,*)'** nb de censures ',(nsujet-cpt)
        

            nz1=nz
            nz2=nz
            if(nz.gt.20)then
               nz = 20
            endif
            if(nz.lt.4)then
               nz = 4
            endif
c     write(*,*)'nombre de noeuds :',nz
c     write(*,*)'nombre de strates:',nst
    


c***************************************************
c--------------- zi- ----------------------------------

c      construire vecteur zi (des noeuds)

         min = 1.d-10
c         min = 0.d0
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
c            print*,'** aux **',aux(i)
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
c         write(*,*)'*-------------------------------------'
c         write(*,*)'** ndate,maxtemps',ndate,maxt
c         write(*,*)'** date **',date(1),date(2),date(3)
c         write(*,*)'** date **',(date(i),i,i=1,ndate)
c         stop
   
         zi(-2) = date(1)
         zi(-1) = date(1)
         zi(0) = date(1)
         zi(1) = date(1)
         h = (date(ndate)-date(1))/dble(nz-1)
         do 18 i=2,nz-1
            zi(i) =zi(i-1) + h  
c         write(*,*)'',zi(i),zi(i-1),i  
 18      continue
         
         zi(nz) = date(ndate)
         zi(nz+1)=zi(nz)
         zi(nz+2)=zi(nz)
         zi(nz+3)=zi(nz)


c         do 188 i=1,ndate     
c         write(*,*)'date(i)',date(i) 
c         write(*,*)'i',i
c 188     continue
          
          
c         write(*,*)'** h **',h
c         write(*,*)'** zi **',zi
c         stop

c---------- affectation nt0,nt1----------------------------

         indictronq=0
            do 50 i=1,nsujet 
               if(t0(i).eq.0.d0)then
                  nt0(i) = 0
               endif
               if(t0(i).ne.0.d0)then
                  indictronq=1
               endif
               do 45 j=1,ndate
                  if(date(j).eq.t0(i))then
                     nt0(i)=j
                  endif
                  if(date(j).eq.t1(i))then
                     nt1(i)=j
                  endif
 45            continue
c                write(*,*)'*** nt1',nt1(1),1
 50       continue   
c          stop

c test sans troncature
c             indictronq=0

c---------- affectation des vecteurs de splines -----------------
             n  = nz+2
             
c      write(*,*)'*** vet avant prog ppal 2** ',(ve(l,1),l,l=1,nsujet) ! ok
             call vecspli3(n,ndate)
c      write(*,*)'*** vet prog ppal 2** ',(ve(l,1),l,l=1,nsujet) ! ok
c      stop
             call vecpen3(n)
             
c             np= nst*n + nva + effet         
c             np = nst*n + nva + 2*effet ! 2 termes de var de la meme frailty        
             np = nst*n + nva + correl + 2*effet 
! term de correlation entre frailty + 2 termes de var de la meme frailty   
             nbpara =np

c             write(*,*)'nombre total de paramètres',np,nst,n,nva,effet
c             write(*,*)'KO****',k0(1),k0(2),ax1,ax2

c------- initialisation des parametres
                   
             do 75 i=1,np
                b(i)= 1.d-1!5.d-1!
 75          continue
c traitement, intercept et pente initialisés plus loin


c***********************************************************
c************** NEW : cross validation  ***********************
c***********************************************************

c yas             effet=0
c yas             correl=0
c         read(2,*)irep1  !=0 pour une recherche du parametre de lissage

c         if(irep1.eq.0.and.nst.ge.2)then
c            write(*,*)'Error : you need only one stratum'
c            stop
c         endif


c         read(2,*)xmin1


c         write(4,*)'STARTING VALUE FOR KAPPA',xmin1
c         write(*,*)'STARTING VALUE FOR KAPPA',xmin1

c yas déb
             nvacross=nva !pour la recherche du parametre de lissage sans var expli
             nva=0
             effetcross=effet
             effet=0
             nstcross=nst
             nst=1

             correl=0

         do 751 l=1,nsujet  
            stracross(l)=stra(l)
 751     continue
         do 752 l=1,nsujet  
            stra(l)=1
 752     continue

c Juan Sep'09         read(2,*)irep1  !=0 pour une recherche du parametre de lissage
         if(irep1.eq.0.and.nst.ge.2)then
            write(*,*)'Error : you need only one stratum'
c ne se produit jamais maintenant (1 seule strate pour CV): changement de juin 2009
            stop
         endif
c Juan Sep'09        read(2,*)xmin1
c Juan Sep'09        if(nst.eq.2)then
c Juan Sep'09         read(2,*)xmin2
c Juan Sep'09         endif

c Juan Sep'09         write(4,*)'STARTING VALUE FOR KAPPA',xmin1
c Juan Sep'09         write(*,*)'STARTING VALUE FOR KAPPA',xmin1
c yas fin
         if(xmin1.le.0.d0)then
            xmin1 = 0.d0
         endif  

c*************************************************
         nvat=nva
         nva=0
  
c yas         if(irep1.eq.1)then   !recherche du parametre de lissage
c yas            xmin1 = dsqrt(xmin1)
c            write(*,*) 'auxi',n,y,ddl,ni,res
c yas            auxi = estimv3(xmin1,n,b,y,ddl,ni,res)

c yas            if (ni.ge.250) then
c              write(*,*) ' '
c              write(*,*) 'no convergence 
c     &             with the chosen smoothing parameter'

c yas             do 175 i=1,nz+2
c yas                b(i)=1.d-1
c yas 175           continue     
c yas             xmin1 = sqrt(10.d0)*xmin1
c yas              auxi = estimv3(xmin1,n,b,y,ddl,ni,res)
c yas              if (ni.lt.250) then
c                 write(*,*)' '
c            write(*,*)'Value of the smoothing parameter :'
c     &       ,real(xmin1*xmin1), '  DoF :',-ddl
c                 write(*,*)' '
c                 write(*,*)'Log-vraisemblance :',res
c                 stop
c yas              else
c yas                 do 176 i=1,nz+2
c yas                    b(i)=1.d-1
c yas 176             continue     
c yas                 xmin1 = sqrt(10.d0)*xmin1
c yas                 auxi = estimv3(xmin1,n,b,y,ddl,ni,res)
c yas                 if (ni.lt.250) then
c                    write(*,*)' '
c            write(*,*)'Value of the smoothing parameter :'
c     &       ,real(xmin1*xmin1), '  DoF :',-ddl
c                    write(*,*)' '
c                    write(*,*)'Log-vraisemblance :',res
c yas                 endif   
c yas               endif
c yas            else
c               write(*,*)' '
c            write(*,*)'Value of the smoothing parameter :'
c     &       ,real(xmin1*xmin1), '  DoF :',-ddl
c               write(4,*)' '
c            write(4,*)'Value of the smoothing parameter :'
c     &       ,real(xmin1*xmin1), '  DoF :',-ddl
c               write(*,*)' '
c               write(*,*)'Log-vraisemblance :',res
c yas            endif

c yas déb

!on  travaille d'abord sur une seule strate
c         nstbis=nst
c         nst=1

         if(irep1.eq.1)then   !pas recherche du parametre de lissage

            xmin1 = dsqrt(xmin1)
c            write(*,*) 'auxi',n,y,ddl,ni,res,xmin1
c        write(*,*)'==avant premier estimv',np,n,k0,xmin1,xmin2,effet,nst
                        
            auxi = estimv3(xmin1,n,b,y,ddl,ni,res)

            if (ni.ge.250) then
c              write(*,*) ' '
c              write(*,*) 'no convergence 
c     &             with the chosen smoothing parameter'

              do 175 i=1,nz+2
                 b(i)=1.d-1
 175           continue     
              xmin1 = sqrt(10.d0)*xmin1
              auxi = estimv3(xmin1,n,b,y,ddl,ni,res)
              if (ni.lt.250) then
c                 write(*,*)' '
c            write(*,*)'Value of the smoothing parameter :'
c     &       ,real(xmin1*xmin1), '  DoF :',-ddl
c                 write(*,*)'Log-vraisemblance :',res
c                 stop
              else
                 do 176 i=1,nz+2
                    b(i)=1.d-1
 176             continue     
                 xmin1 = sqrt(10.d0)*xmin1
                 auxi = estimv3(xmin1,n,b,y,ddl,ni,res)
                 if (ni.lt.250) then
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
c              write(*,*)'Log-vraisemblance :',res
            endif

c yas fin
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
            auxi = estimv3(xmin1,n,b,y,ddl,ni,res)
c            write(*,*)'estimv1',ddl
c               stop
            if(ddl.gt.-2.5d0)then
c            write(*,*)'OK non'
c               stop
               xmin1 = dsqrt(xmin1)
c               stop
               auxi = estimv3(xmin1,n,b,y,ddl,ni,res)
c               write(*,*)'estimv2'
c               stop !pb
               if(ddl.gt.-2.5d0)then
                  xmin1 = dsqrt(xmin1)
                  auxi = estimv3(xmin1,n,b,y,ddl,ni,res)
c                  write(*,*)'estimv3'
                  if(ddl.gt.-2.5d0)then
                     xmin1 = dsqrt(xmin1)
                     auxi = estimv3(xmin1,n,b,y,ddl,ni,res)
c                     write(*,*)'estimv4'
                     if(ddl.gt.-2.5d0)then
                        xmin1 = dsqrt(xmin1)
                        auxi = estimv3(xmin1,n,b,y,ddl,ni,res)
c                        write(*,*)'estimv5'
                        if(ddl.gt.-2.5d0)then
                           xmin1 = dsqrt(xmin1)
                        endif   
                     endif   
                  endif   
               endif
            endif 
c            write(*,*)'estimv1 bis',ddl
            if (ni.ge.250) then
              do 275 i=1,nz+2
                 b(i)=1.d-1
 275          continue     
              xmin1 = sqrt(10.d0)*xmin1
              auxi = estimv3(xmin1,n,b,y,ddl,ni,res)
c             write(*,*)'estimv5'
              if (ni.ge.250) then
                 do 276 i=1,nz+2
                    b(i)=1.d-1
 276             continue     
                 xmin1 = sqrt(10.d0)*xmin1
              endif
            endif 
            ax = xmin1
            bx = xmin1*dsqrt(1.5d0)  
c            write(*,*)'avant mnbrak' ,ax ,bx,cx,fa,fb,fc,n
            
            call mnbrak3(ax,bx,cx,fa,fb,fc,b,n)
c            write(*,*)'OK 2' 
c            stop   
            
            tol = 0.001d0
c            write(*,*)'avant golden',ax,bx,cx,tol,xmin1,n,b,ddl
c            stop
            res = golden3(ax,bx,cx,tol,xmin1,n,b,y,ddl)
c            write(*,*)'apres golden'
c            stop
c            write(4,*)'******************************************* '
c            write(4,*)'Best smoothing parameter',real(xmin1*xmin1)
c     &                , '  DoF :',-ddl
c           write(4,*)'******************************************* '
c            write(*,*)'******************************************* '
c            write(*,*)'Best smoothing parameter',real(xmin1*xmin1)
c     &                , '  DoF :',-ddl
c            write(*,*)'******************************************* '
            effet=0
            correl=0
c            stop


            call marq983((xmin1*xmin1),b,n,ni,v,res,ier,istop)

            if (ni.ge.250) then
c               write(*,*)' '
c               write(*,*) 'non convergence'
            else
c               write(*,*)' '
c               write(*,*)'Log-vraisemblance :',res
            endif
         endif 
c yas déb
         nva=nvacross ! pour la recherche des parametres de regression
         nst=nstcross ! avec stratification si nécessaire
         effet=effetcross ! avec effet initial
         do 753 l=1,nsujet  
            stra(l)=stracross(l) !rétablissement stratification
 753     continue
c yas fin
		 
c yas         nva=nvat    ! pour la recherche des parametres de regression
c        write(*,*)'---nva',nva
c        stop
ccccc********************************************************************


c         if(nst.eq.2)then
c            read(2,*)xmin2
c         endif


c         write(4,*)'nombre de noeuds:',nz
c         read(2,*)fich1 !strate1 sans frailty
c         if(nst.eq.2)then
c            read(2,*)fich2 !strate2 sans frailty
c         endif
c         read(2,*)fich3 !ef_strate1 AVEC frailty
c         if(nst.eq.2)then
c            read(2,*)fich4!ef_strate2 AVEC frailty
c         endif

c         read(2,*)fich1b!surv1 sans frailty
c         if(nst.eq.2)then
c            read(2,*)fich2b!surv2 sans frailty
c         endif
c         read(2,*)fich3b !ef_surv1 AVEC frailty
c         if(nst.eq.2)then
c            read(2,*)fich4b!ef_surv2 AVEC frailty
c         endif

c         read(2,*)fich1c!cumulhazard1 sans frailty
c         if(nst.eq.2)then
c            read(2,*)fich2c!cumulhazard2 sans frailty
c         endif
c         read(2,*)fich3c !ef_cumulhazard1 AVEC frailty
c         if(nst.eq.2)then
c            read(2,*)fich4c!ef_cumulhazard2 AVEC frailty
c         endif

         k0(1) = xmin1*xmin1
         k0(2) = 0.d0
         if(nst.eq.2)then
            k0(2) = xmin2
         endif
c         write(*,*),'lissage',k0(1)
c         write(*,*),'lissage2',k0(2),xmin2
c         stop

C=============================- fin cross validation
         
c===== initialisation des parametres de regression/pas effets aleatopires
c         write(*,*),'====================================='
c         write(*,*),'== avec var explicatives============='
c        write(*,*),'====================================='
         effet=0
         np = nst*n + nva
         b(nst*n+1)=-0.15d0           !initialisation traitementc

c         write(*,*),'===avant marq98==========',n,np,nva,nst,effet

         call marq983(k0,b,np,ni,v,res,ier,istop,effet)
        
c===== recherche de l'ensemble des parametres
         
c         write(*,*),'====================================='
c         write(*,*),'== ensemble des parametres ======='
c         write(*,*),'====================================='
         effet=1
         correl=correlini
         do 255 i=1,nva
            b(np-i+2+correl+1)=b(np-i+1)
 255     continue
         np=nbpara 
c         b(np)= -0.11d0          !initialisation traitement
         if(correl.eq.1)then
c            b(np-nva-2)=-0.25d0!initialisation cov sans contrainte
            b(np-nva-2)=0.1d0!-0.69d0!-1.61d0!-0.64d0!initialisation cov avec contrainte 
         endif
         b(np-nva-1)=0.5d0!0.5477d0      !      initialisation intercept
         b(np-nva)=0.15d0!0.5477d0  !0.13d0!  !   0.15d0!   initialisation pente

c       write(4,*)'======================================'
c       write(4,*)'=== Initialisation parametres ========'
c       if(correl.eq.1)then
c          write(4,*)' pour correl => b(np-nva-2)',b(np-nva-2)
c       endif
c       write(4,*)'sigma2  = ',b(np-nva-1)*b(np-nva-1)
c       write(4,*)'tau2 = ',b(np-nva)*b(np-nva)
c       write(4,*)'======================================'

        
c         write(*,*),'== ensemble des parametres = ',np,nbpara,
c     &                      correl, effet,k0

        
        call marq983(k0,b,np,ni,v,res,ier,istop)
         
c        write(4,*)'*********************************'
c        write(4,*)'nombre d iteration', ni
c        write(4,*)'*********************************'
c        write(4,*)'** log vraisemblance ',res

        j=(np-nva)*(np-nva+1)/2
         
        trace=0
        trace1=0
        trace2=0

c%%%% trace(H-1 * I)=mdf :

         if(indic_sousv.eq.0)then 
cc! que lorsque la mat3(rice totale est inversible, ie sans echec inversion
c strate1 : 
             do 772 i=1,nz1+2
                do 773 j=1,nz1+2
                   H1_hess(i,j)=H_hess(i,j)
                   I1_hess(i,j)=I_hess(i,j)
 773            continue
 772         continue 
             call multi3(H1_hess,I1_hess,nz1+2,nz1+2,nz1+2,HI1)
             do 776 i =1,nz1+2
                trace1=trace1+HI1(i,i)
 776         continue
c             write(4,*)'trace1',trace1
c             write(*,*)'trace1',trace1
c strate2 :
             if(nst.eq.2)then
                do 774 i=1,nz2+2
                   k=nz1+2+i
                   do 775 j=1,nz2+2
                      l=nz1+2+j                   
                      H2_hess(i,j)=H_hess(k,l)
                      I2_hess(i,j)=I_hess(k,l)
 775               continue
 774            continue
                call multi3(H2_hess,I2_hess,nz2+2,nz2+2,nz2+2,HI2)
                do 777 i =1,nz2+2
                   trace2=trace2+HI2(i,i)
 777            continue
c                write(4,*)'trace2',trace2
             endif ! pour la deuxieme strate

c             call multi3(H_hess,I_hess,np,np,np,HI)
c             do 771 i =1,nz1+2
c                trace=trace+HI(i,i)
c 771         continue
c             write(4,*)''
c             write(4,*)'trace1',trace
c             trace=0.d0
c             do 7711 i =1,nz2+2
c                trace=trace+HI(nz1+2+i,nz1+2+i)
c 7711        continue
c             write(4,*)'trace2',trace
          endif
c fin calcul trace

          if(indic_sousv.eq.0)then
             call multi3(I_hess,H_hess,np,np,np,IH)
             call multi3(H_hess,IH,np,np,np,HIH)
          endif
          if(indic_sousv.eq.1)then
             call multi3(I_hess,H_hess,sousm,sousm,sousm,IH)
             call multi3(H_hess,IH,sousm,sousm,sousm,HIH)
          endif

c          if (effet.eq.1.and.correl.eq.1)then
c             write(*,*)'=============================================='
c            write(4,*)'=============================================='
c            write(*,*)'=> cov =', cov
c            write(*,*)'=> b(np-nva-2) =', b(np-nva-2)
c     &         (2.d0*(dexp(b(np-nva-2))/(1.d0+dexp(b(np-nva-2))))-1)
c            write(4,*)'=> cov =', cov
c            write(4,*)'=> b(np-nva-2) =', b(np-nva-2)
c     &         (2.d0*(dexp(b(np-nva-2))/(1.d0+dexp(b(np-nva-2))))-1)
c            if(indic_sousv.eq.0)then
c        write(*,*)'var(b(np-nva-2))=H-1:',H_hess(np-nva-2,np-nva-2)            
c        write(*,*)'var(b(np-nva-2))=H-1 I H -1 :',HIH(np-nva-2,np-nva-2)
c        write(4,*)'var(b(np-nva-2))=H-1:',H_hess(np-nva-2,np-nva-2)            
c        write(4,*)'var(b(np-nva-2))=H-1 I H -1 :',HIH(np-nva-2,np-nva-2)

c        write(4,*)'---------'

c         write(4,*)'var(cov)=:',varcov
c        write(4,*)'et SD(COV)=',dsqrt(varcov)
c        write(4,*)'WALD(COV)=',cov/dsqrt(varcov)
c        write(*,*)'var(cov)=:',varcov
c        write(*,*)'et SD(COV)=',dsqrt(varcov)
c        write(*,*)'WALD(COV)=',cov/dsqrt(varcov)

c        write(4,*)'---------'
c            endif 
c            if(indic_sousv.eq.1)then
c       write(*,*)'var(b(np-nva-2))=H-1:',H_hess(sousm-nva-2,sousm-nva-2)            
c       write(*,*)'var(b(np-nva-2))=H-1IH -1:'
c     &      ,HIH(sousm-nva-2,sousm-nva-2)
c       write(4,*)'var(b(np-nva-2))=H-1:',H_hess(sousm-nva-2,sousm-nva-2)            
c       write(4,*)'var(b(np-nva-2))=H-1IH -1:'
c     &      ,HIH(sousm-nva-2,sousm-nva-2)
c            endif
            


c Salida Juan Aug'07

       f1 = (b(np-nva))*(2.d0* dexp(b(np-nva-2))/
     &       (1.d0+dexp(b(np-nva-2)))-1.d0)
        f2 = (b(np-nva-1))*(2.d0* dexp(b(np-nva-2))/
     &       (1.d0+dexp(b(np-nva-2)))-1.d0)
        f3= 2.d0*(b(np-nva-1))*(b(np-nva))*
     &       (dexp(b(np-nva-2))/
     &       (1.d0+dexp(b(np-nva-2)))-
     &       (dexp(b(np-nva-2))/(1.d0+dexp(b(np-nva-2))))**2)

        varcov=
     &       f1*f1*H_hess(np-nva-1,np-nva-1)+
     &       f2*f2*H_hess(np-nva,np-nva)+
     &       f3*f3*H_hess(np-nva-2,np-nva-2)+
     &       2.d0*f1*f3*H_hess(np-nva-1,np-nva-2)+
     &       2.d0*f2*f3*H_hess(np-nva,np-nva-2)+
     &       2.d0*f1*f2*H_hess(np-nva-1,np-nva)

      if(correl.eq.1)then
            rho=cov/(dsqrt(b(np-nva-1)*b(np-nva-1)*b(np-nva)*b(np-nva)))
      else
           rho=-1
      end if
  
      rhoEnd=rho
      covEnd=cov
      varcovEnd=varcov


      if(indic_sousv.eq.0)then
        varSigma2(1)=((2.d0*b(np-nva-1))**2)*
     &               H_hess(np-nva-1,np-nva-1)
        varSigma2(2)=((2.d0*b(np-nva-1))**2)*HIH(np-nva-1,np-nva-1)

        varTau2(1)=((2.d0*b(np-nva))**2)*H_hess(np-nva,np-nva)
        varTau2(2)=((2.d0*b(np-nva))**2)*HIH(np-nva,np-nva)
      end if

      if(indic_sousv.eq.1)then
        varSigma2(1)=((2.d0*b(np-nva-1))**2)*
     &                H_hess(sousm-nva-1,sousm-nva-1)
        varSigma2(2)=((2.d0*b(np-nva-1))**2)*
     &                HIH(sousm-nva-1,sousm-nva-1)

        varTau2(1)=((2.d0*b(np-nva))**2)*H_hess(sousm-nva,sousm-nva)
        varTau2(2)=((2.d0*b(np-nva))**2)*HIH(sousm-nva,sousm-nva)
      end if  


c Covariates

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


c            if(correl.eq.1)then
c               write(*,*)''
c               write(*,*)'=> CORRELATION=rho',rho
c
c               write(4,*)''
c               write(4,*)'=> CORRELATION=rho',rho     
c            endif
c
c            write(*,*)''
c            write(*,*)'=> sigma2 (var pour intercept aleatoire)=',
c     &           b(np-nva-1)*b(np-nva-1) 
c            write(*,*)'var(sigma2)=H-1:'
c            write(4,*)''
c            write(4,*)'=> sigma2 (var pour intercept aleatoire)=',
c     &           b(np-nva-1)*b(np-nva-1) ,sigma2
c            write(4,*)'var(sigma2)=H-1:'
c            if(indic_sousv.eq.0)then
c            write(*,*)((2.d0*b(np-nva-1))**2)*H_hess(np-nva-1,np-nva-1) !delta method
c           !   =((2.d0*b(np-nva-1))**2)*v(j) avec j=(np-nva)(np-nva+1)/2
c            write(4,*)((2.d0*b(np-nva-1))**2)*H_hess(np-nva-1,np-nva-1) !delta method
c            endif
c            if(indic_sousv.eq.1)then
c            write(*,*)
c     &      ((2.d0*b(np-nva-1))**2)*H_hess(sousm-nva-1,sousm-nva-1) !delta method
c            write(4,*)
c     &      ((2.d0*b(np-nva-1))**2)*H_hess(sousm-nva-1,sousm-nva-1) !delta method
c            endif
c           
c            write(*,*)'var(sigma2)=H-1 I H -1 :'
c            write(4,*)'var(sigma2)=H-1 I H -1 :'
c            if(indic_sousv.eq.0)then
c            write(*,*)((2.d0*b(np-nva-1))**2)*HIH(np-nva-1,np-nva-1)
c            write(4,*)((2.d0*b(np-nva-1))**2)*HIH(np-nva-1,np-nva-1)
c            endif
c            if(indic_sousv.eq.1)then
c            write(*,*)
c     &        ((2.d0*b(np-nva-1))**2)*HIH(sousm-nva-1,sousm-nva-1)
c            write(4,*)
c     &        ((2.d0*b(np-nva-1))**2)*HIH(sousm-nva-1,sousm-nva-1)
c            endif
            
c            write(*,*)''
c            write(*,*)'Interaction avec  : ',nominter
c            write(*,*)'=> tau2( var pour pente aleatoire)=',
c     &           b(np-nva)*b(np-nva),tau2 
c            write(*,*)'var(tau2)=H-1:'
c            write(4,*)''
c            write(4,*)'Interaction avec  : ',nominter
c            write(4,*)'=> tau2( var pour pente aleatoire)=',
c     &           b(np-nva)*b(np-nva) 
c            write(4,*)'var(tau2)=H-1:'

c            if(indic_sousv.eq.0)then
c               write(*,*)((2.d0*b(np-nva))**2)*H_hess(np-nva,np-nva)
c               write(4,*)((2.d0*b(np-nva))**2)*H_hess(np-nva,np-nva)
c            endif
c            if(indic_sousv.eq.1)then
c             write(*,*)((2.d0*b(np-nva))**2)*H_hess(sousm-nva,sousm-nva)
c             write(4,*)((2.d0*b(np-nva))**2)*H_hess(sousm-nva,sousm-nva)
c            endif
c            write(*,*)'var(tau2)=H-1 I H -1 :'
c            write(4,*)'var(tau2)=H-1 I H -1 :'
c            if(indic_sousv.eq.0)then
c            write(*,*)((2.d0*b(np-nva))**2)*HIH(np-nva,np-nva)
c            write(4,*)((2.d0*b(np-nva))**2)*HIH(np-nva,np-nva)
c            endif
c            if(indic_sousv.eq.1)then
c            write(*,*)((2.d0*b(np-nva))**2)*HIH(sousm-nva,sousm-nva)
c            write(4,*)((2.d0*b(np-nva))**2)*HIH(sousm-nva,sousm-nva)
c            endif
            
            
            
c            if(indic_codage.eq.0)then
c            write(4,*)'==========='
c            write(4,*)'===> mat3(RICE DE VARIANCE COVARIANCE
c     &           (u_i+v_i x_ij1)'
c            IF(correl.eq.1)then
c               write(4,*)'sigma2=',real(b(np-nva-1)*b(np-nva-1)),
c     &real(sigma2) ,
c     &              'sigma2+cov=',
c     &           real(b(np-nva-1)*b(np-nva-1)+  
c     &           cov)
c               write(4,*)'sigma2+cov=',
c     &           real(b(np-nva-1)*b(np-nva-1)+  
c     &           cov),
c     &              'sigma2+2cov +tau2=',
c     &           real(b(np-nva-1)*b(np-nva-1)+  
c     &        2.d0*cov+
c     &              b(np-nva)*b(np-nva))
c            else
c               write(4,*)'sigma2=',real(b(np-nva-1)*b(np-nva-1)) ,
c     &              'sigma2=',real(b(np-nva-1)*b(np-nva-1) )
c               write(4,*)'sigma2=',real(b(np-nva-1)*b(np-nva-1)) ,
c    & 'sigma2+tau2=',real(b(np-nva-1)*b(np-nva-1) +b(np-nva)*b(np-nva))
c            endif
c            write(4,*)''
c            endif
c            endif

c            if(indic_codage.eq.1)then
c            write(4,*)'==========='
c            write(4,*)'===> mat3(RICE DE VARIANCE COVARIANCE
c     &           (u_i+v_i x_ij1)'
c            IF(correl.eq.1)then
c               write(4,*)'sigma2-cov+tau2/4='
c     &              ,real(b(np-nva-1)*b(np-nva-1)-cov
c     &+b(np-nva)*b(np-nva)/4.0),
c     &              'sigma2-tau2/4=',
c     &           real(b(np-nva-1)*b(np-nva-1)-  
c     &            b(np-nva)*b(np-nva)/4.0)
c               write(4,*)'sigma2-tau2/4=',
c     &           real(b(np-nva-1)*b(np-nva-1)-  
c     &            b(np-nva)*b(np-nva)/4.0),
c     &              'sigma2+cov +tau2/4=',
c     &              real(b(np-nva-1)*b(np-nva-1)+cov
c     &+b(np-nva)*b(np-nva)/4.0)
c            else
c               write(4,*)'sigma2+tau2/4=',real(b(np-nva-1)*b(np-nva-1)
c     &+b(np-nva)*b(np-nva)/4.0),
c     &              'sigma2-tau2/4=',real(b(np-nva-1)*b(np-nva-1)
c     &-b(np-nva)*b(np-nva)/4.0)
c               write(4,*)'sigma2-tau2/4=',real(b(np-nva-1)*b(np-nva-1)
c     &-b(np-nva)*b(np-nva)/4.0),
c     &              'sigma2+tau2/4=',real(b(np-nva-1)*b(np-nva-1)
c c    &+b(np-nva)*b(np-nva)/4.0)
c c           endif
c c           write(4,*)''
c            endif


         if(effet.eq.1.and.ier.eq.-1)then
            v((np-nva)*(np-nva+1)/2)=10.d10
          endif
                  
c          write(*,*)'valeur de ni',ni 
       
          j=(np-nva)*(np-nva+1)/2


c ------------------------------------ Covariates

c         if(nva.gt.0)then
c            write(4,*)'=============================================='
c            write(4,*)'=============================================='
c         do 124 i=1,nva
c            j=(np-nva+i)*(np-nva+i+1)/2
c            if(indic_sousv.eq.0)then
c            bi = b(np-nva+i) - 1.96*dsqrt(H_hess(np-nva+i,np-nva+i))
c            bs = b(np-nva+i) + 1.96*dsqrt(H_hess(np-nva+i,np-nva+i))
c            endif
c            if(indic_sousv.eq.1)then
c          bi = b(np-nva+i) - 1.96*dsqrt(H_hess(sousm-nva+i,sousm-nva+i))
c          bs = b(np-nva+i) + 1.96*dsqrt(H_hess(sousm-nva+i,sousm-nva+i))
c            endif
c            if(indic_sousv.eq.1)then
c       write(4,*)i,')',' SE (=H)',dsqrt(H_hess(sousm-nva+i,sousm-nva+i))
c            wres=(b(np-nva+i))/dsqrt(H_hess(sousm-nva+i,sousm-nva+i))
c            endif
c            write(4,*)'Variable : ',nomvar(i)
c            write(4,*)i,')','beta=',b(np-nva+i)
c            write(4,*)' '
c            if(indic_sousv.eq.0)then
c            write(4,*)i,')',' SE (=H)',dsqrt(H_hess(np-nva+i,np-nva+i))
c            wres=(b(np-nva+i))/dsqrt(H_hess(np-nva+i,np-nva+i))
c            endif
c            write(4,*)'---> WALD',wres
c            write(4,*)' '
c            if(indic_sousv.eq.0)then
c            write(4,*)i,')',' SE (=HIH)',dsqrt(HIH(np-nva+i,np-nva+i))
c       write(4,*)'---> WALD',(b(np-nva+i))/dsqrt(HIH(np-nva+i,np-nva+i))
c       endif
c            if(indic_sousv.eq.1)then
c        write(4,*)i,')',' SE (=HIH)',dsqrt(HIH(sousm-nva+i,sousm-nva+i))
c       write(4,*)'---> WALD',
c     &       (b(np-nva+i))/dsqrt(HIH(sousm-nva+i,sousm-nva+i))
c      endif
c            write(4,*)' '
c            write(4,*)'RR : ',dexp(b(np-nva+i)),'  IC',dexp(bi),dexp(bs) 
c            write(4,*)'=============================================='
c            write(4,*)'=============================================='
c 124     continue 
c         write(4,*)'---> log vraisemb marginale complete pénalisée',res
c         write(4,*)'=============================================='
c         write(4,*)'=============================================='
c      endif
      
c --------------  Lambda and survival and cumulative hazard estimat3(es 
      
c      if(effet.eq.2)then
c         call distance3(nz1,nz2,b,fich3,fich4,fich3b,fich4b,
c     & fich3c,fich4c,effet)
c      endif

       call distance3(nz1,nz2,b,effet,
     &	       x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)      

      

c      call date_and_time(dateamj,heure2,zone,values)     
c      write(4,*) '***************************************************'   
c      write(4,*) '*** starting time: ***', dateamj,heure1
c      write(4,*) '*** Ending time:(hhmmss.sss) ***', dateamj,heure2

      end subroutine 


    


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC                             CCCCCCCCCCCCCCC
CCCCCCCCCCCCC          SUBROUTINES        CCCCCCCCCCCCCCC
CCCCCCCCCCCCC                             CCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c========================== VECSPLI =====================
      subroutine vecspli3(n,ndate) 
      
      integer ::  n,ndate,i,j,k
      double precision :: ht,htm,h2t,ht2,ht3,hht,h,hh,h2
      double precision :: h3,h4,h3m,h2n,hn,hh3,hh2
         
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************

c*****dace1 
      double precision , dimension(ndatemax) ::date
      double precision , dimension(-2:npmax) :: zi
      common /dace1/date,zi
c*****mem1
      double precision ,dimension(ndatemax) :: mm3,mm2,mm1,mm
      common /mem1/mm3,mm2,mm1,mm
c*****mem2
      double precision ,dimension(ndatemax) :: im3,im2,im1,im
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
      subroutine vecpen3(n) 
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************

      integer ::  n,i
      double precision :: h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2
      double precision :: a3,a2,b2,c2,a1,b1,c1,a0,x3,x2,x

c*****dace1 
      double precision , dimension(ndatemax) ::date
      double precision , dimension(-2:npmax) :: zi
      common /dace1/date,zi

c*****pen1
      double precision ,dimension(npmax) ::m3m3,m2m2,m1m1,mmm,m3m2
      common /pen1/m3m3,m2m2,m1m1,mmm,m3m2
c*****pen2
      double precision ,dimension(npmax) ::m3m1,m3m,m2m1,m2m,m1m
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
 20      continue
         end

c================================  DCHOLE  ===========================
      subroutine dchole3(a,k,nq,idpos)

c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
      
      integer :: k,nq,i,ii,i1,i2,i3,m,is,j,k2,jmk
      integer :: ijm,irm,jji,jjj,l,jj,iil,jjl,il,idpos
      double precision ,dimension(npmax*(npmax+3)/2) :: a
      double precision :: term,xn,diag,p
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
c     calcul des elements de la mat3(rice
      do 13 i=1,k
      ii=i*(i+1)/2
c       elements diagonaux
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
c=================================  DMFSD  ===================================

      SUBROUTINE DMFSD3(A,N,EPS,IER)
C
C   FACTORISATION DE CHOLESKY D'UNE mat3(RICE SDP
C   mat3(RICE = TRANSPOSEE(T)*T
C   ENTREE : TABLEAU A CONTENANT LA PARTIE SUPERIEURE STOCKEE COLONNE
C            PAR COLONNE DE LA METRICE A FACTORISER
C   SORTIE : A CONTIENT LA PARTIE SUPPERIEURE DE LA mat3(RICE SYMETRIQUE T
C 
C   SUBROUTINE APPELE PAR DSINV
C  
C   N : DIM. mat3(RICE 
C   EPS : SEUIL DE TOLERANCE
C   IER = 0 PAS D'ERREUR
C   IER = -1 ERREUR
C   IER = K COMPRIS ENTRE 1 ET N, WARNING, LE CALCUL CONTINUE
C
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
      integer , parameter ::npmax2=(Npmax*(Npmax+1)/2) 

      DOUBLE PRECISION  A(npmax2)!((N*(N+1)/2)) :: A
      DOUBLE PRECISION :: DPIV,DSUM,EPS,TOL
      INTEGER :: I,K,L,N,IER,KPIV,IND,LEND,LANF,LIND
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

      SUBROUTINE DSINV3(A,N,EPS,IER)
C
C     INVERSION D'UNE mat3(RICE SYMETRIQUE DEFINIE POSITIVE :
C
C     mat3(RICE = TRANSPOSEE(T)*T
C     INERSE(mat3(RICE) = INVERSE(T)*INVERSE(TRANSPOSEE(T))
C
C     A : TABLEAU CONTENANT LA PARTIE SUPERIEURE DE LA mat3(RICE A INVERSER
C         STOCKEE COLONNE PAR COLONNE
C     DIM. mat3(RICE A INVERSER = N 
C     DIM. TABLEAU A = N*(N+1)/2
C
C     EPS : SEUIL DE TOLERANCE AU-DESSOUS DUQUEL UN PIVOT EST CONSIDERE
C           COMME NUL
C
C     IER : CODE D'ERREUR
C         IER=0 PAS D'ERREUR
C         IER=-1 ERREUR SUR LA DIM.N OU mat3(RICE PAS DEFINIE POSITIVE
C         IER=1 PERTE DE SIGNIFICANCE, LE CALCUL CONTINUE
C
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
      integer , parameter ::npmax2=(Npmax*(Npmax+1)/2) 

      DOUBLE PRECISION,dimension(npmax2) :: A
      DOUBLE PRECISION :: DIN,WORK
      DOUBLE PRECISION :: EPS
      INTEGER ::  N,IER,IND,IPIV,I,J,K,L,MIN,KEND
      INTEGER :: LHOR,LVER,LANF
C
C     FACTORIZE GIVEN mat3(RIX BY MEANS OF SUBROUTINE DMFSD
C     A=TRANSPOSE(T) * T
C
      CALL DMFSD3(A,N,EPS,IER)
      IF(IER) 9,1,1
C
C     INVERT UPPER TRIANGULAR mat3(RIX T
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

      DOUBLE PRECISION FUNCTION MAXT3(DELTA,M)
c************** definition commune des parameter ***********************

c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
      DOUBLE PRECISION ,dimension(npmax)::DELTA
      INTEGER :: I,M
c
      MAXT3=DELTA(1)
      DO 2 I=2,M
      IF(DELTA(I).GT.MAXT3) MAXT3=DELTA(I)
2     CONTINUE
      RETURN
      END



c================================  SEARPAS joly    ==============================

      SUBROUTINE SEARPAS3(VW,STEP,B,BH,M,DELTA,FIM,EPSV,k0)
C
C  MINIMISATION UNIDIMENSIONNELLE
C
       INTEGER :: I,M
       DOUBLE PRECISION :: VLW,VLW1,VLW2,VLW3,VW,VM
       DOUBLE PRECISION :: FI1,FI2,FI3,FIM,EPSV
       DOUBLE PRECISION :: STEP
       DOUBLE PRECISION , dimension(2):: k0
       DOUBLE PRECISION , dimension(M):: B,BH,DELTA 
C
       VLW1=DLOG(VW)
       VLW2=VLW1+STEP
       CALL valfpa3(VLW1,FI1,B,BH,M,DELTA,k0)
       CALL valfpa3(VLW2,FI2,B,BH,M,DELTA,k0)
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
          CALL valfpa3(VLW1,FI1,B,BH,M,DELTA,k0)   
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
          CALL valfpa3(VLW1,FI1,B,BH,M,DELTA,k0)
          IF(FI1.GT.FI2) GO TO 50
c          IF (dabs(FI1-FI2).LT.EPSV) THEN
          IF (FI1.eq.FI2) THEN
             FIM=FI2
             VM=VLW2
             GO TO 100
          ENDIF
 20    CONTINUE
C
C  PHASE 2 APPROXImat3(ION PAR QUADRIQUE
C
50     CONTINUE
C
C  CALCUL MINIMUM QUADRIQUE
C
         VM=VLW2-STEP*(FI1-FI3)/(2.d0*(FI1-2.d0*FI2+FI3))
         CALL valfpa3(VM,FIM,B,BH,M, DELTA,k0)
         IF (FIM.LE.FI2) GO TO 100
         VM=VLW2
         FIM=FI2
100   CONTINUE
      VW=DEXP(VM)
      RETURN

      END


   
c===================================   VALFPA joly   ==============================

        subroutine valfpa3(vw,fi,b,bk,m,delta,k0)
        integer :: m,i
        double precision :: vw,fi
        double precision :: funcpa3,z        
        double precision ,dimension(2) :: k0
        double precision ,dimension(M) :: b,bk,delta

        z=0.d0
        do 1 i=1,m
           bk(i)=b(i)+dexp(vw)*delta(i)
 1      continue
c     
        fi=-funcpa3(bk,m,1,z,1,z,k0)
c
         return
         end   
       
c=================================    DERIVA pjoly =======================

      subroutine deriva3(b,m,v,rl,k0)
            
      integer ::i0,iun,m,m1,ll,i,k,j
      double precision :: funcpa3,thn,th,z,rl,vl
      double precision :: th2
      double precision ,dimension(2):: k0
      double precision ,dimension(m):: fcith
      double precision ,dimension(m):: b
      double precision , dimension(m*(m+3)/2)::V
      
c     v:mat3(rice d'informat3(ion+score
c     calcul de la derivee premiere
c     
c      print*,'entree deriva'
      th=5.d-3!1.d-4!1.d-5!
      thn=-th
      th2=th*th
      z=0.d0
      i0=0
      iun =1
      rl=funcpa3(b,m,iun,z,iun,z,k0)
      do 2 i=1,m
         fcith(i)=funcpa3(b,m,i,th,i0,z,k0)
c         print*,'*** fcith :',fcith(i),b(i),i,m
 2    continue
      
      k=0
      m1=m*(m+1)/2
      ll=m1
      do 1 i=1,m
         ll=ll+1
         vl=(fcith(i)-funcpa3(b,m,i,thn,i0,z,k0))/(2.d0*th)
c         print*,'*** - funcpa :', funcpa3(b,m,i,thn,i0,z,k0),b(i),i,m
         v(ll)=vl
         do 1 j=1,i
            k=k+1
            v(k)=-(funcpa3(b,m,i,th,j,th,k0)-fcith(j)-fcith(i)+rl)/th2 
c            print*,'*** v(k)',v(k),k,j,i
 1       continue
         
         return
         end


c==========================  DISTANCE   =================================
    
   
         subroutine distance3(nz1,nz2,b,effet,
     &	       x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)



c*************** SAME **********************************************
      integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
      integer , parameter ::ngmax=1000,nboumax=1000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1000
c******************************************************************

         integer :: nz1,nz2,i,j,n,np,k,l,effet

         double precision :: su1,ri1,su2,ri2,x1,x2,h1,h2
         double precision :: d,h
         double precision :: su,bsup,binf,lam,lbinf,lbsup,margi
         double precision ,dimension(npmax,npmax)::HIH,IH,hes1,hes2
         double precision ,dimension(-2:npmax):: the1,the2
         double precision ,dimension(npmax):: b
         character*20 fic1,fic2,fic1b,fic2b,fic1c,fic2c

         double precision  x1Out(99),lamOut(99,3),suOut(99,3),
     &                     x2Out(99),lam2Out(99,3),su2Out(99,3)

c*****dace1 
         double precision , dimension(ndatemax) ::date
         double precision , dimension(-2:npmax) :: zi
         common /dace1/date,zi

c*****dace2
      double precision , dimension(nsujetmax) ::t0,t1
      double precision :: ncsrl
      integer , dimension(nsujetmax) ::c
      integer, dimension(nsujetmax) :: nt0,nt1
      integer :: nsujet,nsim,nva,ndate,nst
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,nst
c*****dace7
         double precision , dimension(npmax,npmax) :: I_hess,H_hess
         double precision , dimension(npmax,npmax) :: Hspl_hess
         double precision , dimension(npmax,1) ::PEN_deri
         common /dace7/PEN_deri,I_hess,H_hess,Hspl_hess
c***********************************************

  
         n  = nz1+2
         if(nst.eq.2)then      
            np  = nz1+2+nz2+2+effet+nva
         else
            np  = nz1+2+effet+nva
         endif
            
        
             do 116 i=1,nz1+2
                do 115 j=1,nz1+2
                   hes1(i,j)=h_Hess(i,j)
 115            continue
 116         continue 

         if(nst.eq.2)then  
            k = 0
            do 118 i=nz1+3,nz1+2+nz2+2
               k = k + 1 
               l = 0
               do 117 j=nz1+3,nz1+2+nz2+2
                  l = l + 1
                  hes2(k,l)=h_Hess(i,j)
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
             call cosp3(x1,the1,nz1+2,hes1,zi,date,binf,su,bsup
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
               call cosp3(x2,the2,nz2+2,hes2,zi,date,binf,su,bsup
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
      subroutine SUPS3(x,the,n,su,lam,zi)

c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
      integer :: j,k,n,i
      double precision :: x,ht,ht2,h2,som,lam,su
      double precision :: htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3
      double precision :: ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2
      double precision :: h,gl,hh
      double precision ,dimension(-2:npmax) :: zi,the

      som = 0.d0
      do 58 k = 2,n+1
         if ((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
            j = k-1
            if (j.gt.1)then
               do 55 i=2,j
                  som = som+the(i-4)
 55            continue  
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

      subroutine cosp3(x,the,n,y,zi,date,binf,su,bsup,lbinf,lam,lbsup)
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************

      integer :: j,k,n,i
      double precision :: x,ht,ht2,h2,som,lam,su
      double precision :: binf,bsup,lbinf,lbsup,pm
      double precision :: htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3
      double precision :: ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2
      double precision :: h,gl,hh
      double precision ,dimension(-2:npmax):: the,zi
      double precision ,dimension(npmax,npmax):: y
      double precision ,dimension(ndatemax):: date

c         do 4 i=1,n
c            the(i-3)=(b(i))*(b(i))
c 4       continue

         som = 0.d0
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

         call conf3(x,j,n,y,pm,zi,date)

         binf = dexp(-gl + 1.96d0*pm)
         su  = dexp(-gl)
         bsup = dexp(-gl - 1.96d0*pm)

         call conf13(x,j,n,y,pm,zi,date)
         lbinf = lam - 1.96d0*pm
         lbsup = lam + 1.96d0*pm

         return

         end
c=====================  CONF1  =============================
      subroutine  conf13(x,ni,n,y,pm,zi,date)

c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************



      integer :: ni,i,n,j

      double precision :: mmsp3,x,pm
      double precision :: res
      double precision,dimension(ndatemax) :: date
      double precision,dimension(-2:npmax) :: zi
      double precision,dimension(npmax) :: vecti,aux
      double precision,dimension(npmax,npmax) :: y
      
            do 10 i=1,n
               vecti(i) = mmsp3(x,ni,i,zi,date)
 10         continue   

            do 20 i=1,n
               aux(i) = 0.d0
               do 15 j=1,n
                  aux(i) = aux(i) - y(i,j)*vecti(j)
 15            continue
 20         continue   

            res = 0.d0
            do 30 i=1,n
               res = res + aux(i)*vecti(i)
 30         continue
c            write(4,*)'conf1 : res ',res
            if (res.lt.0)then 
               res=-res
            endif 
c               pm = dsqrt(2.d0*res)
               pm = dsqrt(res)
             
c            write(4,*)'conf1 pm',pm 
         end
c=====================  CONF  =============================
      subroutine  conf3(x,ni,n,y,pm,zi,date)
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
       
       integer :: ni,i,n,j
       double precision :: isp3,x,pm
       double precision :: res
       double precision,dimension(-2:npmax) :: zi
       double precision,dimension(ndatemax) :: date
       double precision,dimension(52) :: vecti,aux
       double precision,dimension(npmax,npmax) :: y


            do 10 i=1,n
               vecti(i) = isp3(x,ni,i,zi,date)
 10         continue   

            do 20 i=1,n
               aux(i) = 0.d0
               do 15 j=1,n
                  aux(i) = aux(i) - y(i,j)*vecti(j)
 15            continue
 20         continue   

            res = 0.d0
            do 30 i=1,n
               res = res + aux(i)*vecti(i)
 30         continue

c            write(4,*)'conf : res ',res
            if (res.lt.0)then 
               res=-res
            endif 
c               pm = dsqrt(2.d0*res)
               pm = dsqrt(res)
              
c            write(4,*)'conf : pm ',pm
               
         end


c==========================   ISP   ==================================
          double precision function isp3(x,ni,ns,zi,date)
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
          integer :: ni,ns
          double precision :: val,mmsp3,x
          double precision ,dimension(-2:npmax)::zi
          double precision ,dimension(ndatemax)::date

          if(x.eq.zi(ni))then
             if(ni.le.ns-3)then
                val = 0.d0
             else
                if(ni.le.ns-2)then
           val = ((zi(ni)-zi(ni-1))*mmsp3(x,ni,ns,zi,date))*0.25d0
                else
                   if (ni.eq.ns-1)then
                  val = ((zi(ni)-zi(ni-2))*mmsp3(x,ni,ns,zi,date)+
     &        (zi(ni+3)-zi(ni-1))*mmsp3(x,ni,ns+1,zi,date))*0.25d0
                   else
                      if(ni.eq.ns)then
                  val = ((zi(ni)-zi(ni-3))*mmsp3(x,ni,ns,zi,date)+
     &                (zi(ni+2)-zi(ni-2))*mmsp3(x,ni,ns+1,zi,date)
     &       +(zi(ni+3)-zi(ni-1))*mmsp3(x,ni,ns+2,zi,date))*0.25d0
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
                 val = (x-zi(ni))*mmsp3(x,ni,ns,zi,date)*0.25d0
             else  
             if(ni.eq.ns-2)then
                   val = ((x-zi(ni-1))*mmsp3(x,ni,ns,zi,date)+
     &       (zi(ni+4)-zi(ni))*mmsp3(x,ni,ns+1,zi,date))*0.25d0
             else   
                if (ni.eq.ns-1)then
                   val =((x-zi(ni-2))*mmsp3(x,ni,ns,zi,date)+
     &             (zi(ni+3)-zi(ni-1))*mmsp3(x,ni,ns+1,zi,date)
     &       +(zi(ni+4)-zi(ni))*mmsp3(x,ni,ns+2,zi,date))*0.25d0
                else
                   if(ni.eq.ns)then
                      val =((x-zi(ni-3))*mmsp3(x,ni,ns,zi,date)+
     &             (zi(ni+2)-zi(ni-2))*mmsp3(x,ni,ns+1,zi,date)
     &             +(zi(ni+3)-zi(ni-1))*mmsp3(x,ni,ns+2,zi,date)
     &        +(zi(ni+4)-zi(ni))*mmsp3(x,ni,ns+3,zi,date))*0.25d0
                   else
                      val = 1.d0
                   endif
                endif
             endif
             endif
          endif 
          endif
             isp3 = val
             return
             end
c==========================  MMSP   ==================================
          double precision function mmsp3(x,ni,ns,zi,date)
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
      integer :: ni,ns
      double precision :: val,x
      double precision ,dimension(-2:npmax)::zi
      double precision ,dimension(ndatemax)::date

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

             mmsp3 = val
             return
             end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



C****************************************************************
      
	SUBROUTINE BGOS(SX,ID,X1,X2,RO)
c ID=1:U(0,SX); ID DIFF DE 1 :N(0,SX)
c bgos peut renvoyer deux nb aleat X1 et X2 
c rO est leur correlation, si r0=0, var independantes
c sx = ecart type de X1 ou X2        
      double precision :: uniran,sx,x1,x2,RO,RO2
      double precision :: f,v1,v2,s,dls
      integer :: id
      real,dimension(2) :: rx

 5 	CONTINUE
        call RANDOM_number(rx)
        X1=rx(1)
        X2=rx(2)
c 	X1=rand(idum)
c	X2=rand(idum) 
c	X1=uniran()
c	X2=uniran()
	 IF(ID.NE.1) GO TO 10
	F=2.d0*DSQRT(3.d0)
	X1=(X1-0.5)*F
	X2=(X2-0.5)*F
	GO TO 20
10 	CONTINUE
	V1=2.*X1-1
	V2=2.*X2-1
	S=V1*V1+V2*V2
	IF(S.GE.1.d0) GO TO 5
	DLS=SQRT(-2.*LOG(S)/S)
	X1=V1*DLS
	X2=V2*DLS
20	CONTINUE
	RO2=RO*RO
	IF(ABS(RO).GT.1.E-10) X2=(X1+X2*SQRT(1./RO2-1.))*RO
	X1=X1*SX
	X2=X2*SX
	RETURN
      END

c***************************************************************


	SUBROUTINE VAREXP2(Y)
        DOUBLE PRECISION :: Y,uniran
        
        y=uniran()

c	Y=-DLOG(1.d0-Y)/1.414213d0
	Y=-DLOG(1.d0-Y)
	RETURN
      	END

c================== multiplication de mat3(rice  ==================

c multiplie A par B avec le resultat dans C

	subroutine multi3(A,B,IrowA,JcolA,JcolB,C)
c remarque :  jcolA=IrowB
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************

        integer ::IrowA,JcolA,JcolB,i,j,k
        double precision ::sum
        double precision ,dimension(npmax,npmax) ::A,B,C
                
        do 1 I=1,IrowA
           do 2 J=1,JcolB
              sum=0
              do 3 K=1,JcolA
                 sum=sum+A(I,K)*B(K,J)
 3            continue
              C(I,J)=sum
 2            continue
 1      continue
        return
        end


c===============================    MARQ98  HJ =========================

   	 subroutine marq983(k0,b,m,ni,v,rl,ier,istop)

c  fu = matrice des derivees secondes et premieres
c
c  istop: raison de l'arret
c  1: critere d'arret satisfait (prm=ca, vraisblce=cb, derivee=dd)
c  2: nb max d'iterations atteints
c  3: 1 mais echec inversion matrice d'info (ier=1)
C  4: non amelioration vraisblce apres searpas


c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
 
         integer :: m,ni,nql,i,ii,nfmax,idpos,ier,istop,igrad,j
         integer :: ncount,id,jd,i0

         double precision :: da,dm,ga,tr,g0
         double precision :: ca,cb,epsa,epsb,rl
         double precision :: funcpa3,det, step,eps,epsd
         double precision :: vw,fi,maxt,z
         double precision :: rl1,th,ep,dd

         double precision , dimension(npmax*(npmax+3)/2)::V
         double precision , dimension(npmax*(npmax+3)/2)::V1
         double precision , dimension(npmax*(npmax+3)/2)::vnonpen
         double precision , dimension(m)::b,delta
         double precision , dimension(npmax*(npmax+3)/2) :: fu
         double precision , dimension(2) :: k0
         double precision , dimension(m) :: bh,b1
         double precision , dimension(2) :: zero
         character*20 heure

c***** nmax
         integer :: nmax
         common /nmax/nmax

c*****dace2
      double precision , dimension(nsujetmax) ::t0,t1
      double precision :: ncsrl
      integer , dimension(nsujetmax) ::c
      integer, dimension(nsujetmax) :: nt0,nt1
      integer :: nsujet,nsim,nva,ndate,nst
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,nst

c*****dace3
      double precision :: pe
      integer :: effet,nz1,nz2,correl
      common /dace3/pe,effet,nz1,nz2,correl
c*****dace7
      double precision , dimension(npmax,npmax) :: I_hess,H_hess
      double precision , dimension(npmax,npmax) :: Hspl_hess
      double precision , dimension(npmax,1) ::PEN_deri
      common /dace7/PEN_deri,I_hess,H_hess,Hspl_hess
c***********************************************

c*****dace1 
c      double precision , dimension(ndatemax) ::date
c      double precision , dimension(-2:npmax) :: zi
c     common /dace1/date,zi


      zero(1)=0.d0
      zero(2)=0.d0
      id=0
      jd=0
      z=0.d0
      th=1.d-5
      eps=1.d-7
      epsa=1.d-4!1.d-3
      epsb=1.d-4!1.d-3
      epsd=1.d-3!1.d-2
      nfmax=m*(m+1)/2
      ca=epsa+1.d0
      cb=epsb+1.d0
      rl1=-1.d+10
      ni=0
      istop=0
      da=0.01d0
      dm=5.d0
c     nql=nbre de seconds membres pour chole
      nql=1

 10   continue
c     
c      write(6,*)
c      write(*,*)'** ITERATION :**', ni 
c         write(*,*)'** Parametres:'
c         write(*,*)(b(i),i=1,m)
c         write(6,*)'da',da,' ga=',ga
 
      z=0.d0
      i0=0 

c      write(*,*) "hola VR......"
c      write(*,*) m, i0, z, k0
c      stop

      rl=funcpa3(b,m,i0,z,i0,z,k0)
      rl1=rl 

c      write(*,*)'log Vrais rl1,nsujet',rl1,nsujet
c      stop
      if (RL1.gt.0) then
        write(*,*) "problems ..." 
        stop
      end if
      
      call deriva3(b,m,v,rl,k0)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      do 144 i=1,m*(m+1)/2
c          ii=i*(i+1)/2
c          write(*,*) 'toutes derivees 2nds marq ',i,' = ',v(i)
c  144  continue 
c       ii=0
c       do 155 i=m*(m+1)/2+1,m*(m+3)/2
c          ii=ii+1
c       write(*,*) 'derivee premiere % parametre, marq',ii,' = ',v(i)
c  155  continue
cccccccccccccccccccccccc
c      stop
          dd = 0.d0
          do 13 i=m*(m+1)/2+1,m*(m+3)/2
c             write(*,*)'d',v(i)
             dd = dd + v(i)*v(i)
 13       continue 
          dd=dd/dabs(RL)

c          pause
c          write(6,*) (b(ii),ii=1,m)
          if (ca.lt.epsa.and.cb.lt.epsb.and.dd.lt.epsd) goto 100
c          if (ca.lt.epsa.and.cb.lt.epsb) goto 100
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
c          write(*,*)'avant dchole',fu(nfmax+1),m,nql,idpos
          call dchole3(fu,m,nql,idpos)
c          write(*,*)'apres dchole',fu(nfmax+1),m,nql,idpos
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
c              write(*,*)'ncount=',ncount,'da=',da,'ga=',ga
              do 2 i=1,m
                 delta(i)=fu(nfmax+i)
                 b1(i)=b(i)+delta(i)
2             continue
c              write(6,*) '** ** avant func.....B1',b1
c              write(6,*) '**'
c              write(6,*) '** ** avant func.....B:',b
c              write(6,*) '** ** avant func..',b1,m,id,z,jd,z,k0
              rl=funcpa3(b1,m,id,z,jd,z,k0)
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
           call searpas3(vw,step,b,bh,m,delta,fi,eps,k0) 
c           write(*,*) 'apres searpas'
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

 800       cb=dabs(rl1-rl)
           ca=0.d0
           do 5 i=1,m
              ca=ca+delta(i)*delta(i)
 5         continue

c           write(6,*) 'ca =',ca,' cb =',cb,' dd =',dd
    
           do 6 i=1,m
              b(i)=b(i)+delta(i)
c              write(*,*)'delta',delta(i),i
 6         continue

           ni=ni+1
           if (ni.gt.nmax) then
              istop=2
c              write(6,*) 'nombre iteration max atteinte'
              goto 110
           end if
           goto 10
c********************
c     inversion matrice d'information
c********************
 100       continue
           istop=1

c           write(*,*)'ca,cb,dd final',ca,cb,dd


c================ pour les bandes de confiance/delta method

           call deriva3(b,m,v,rl,k0)
           do 303 i=1,(m*(m+3)/2)
              v1(i)=0.d0
 303       continue
           
           m1=m-nva-2*effet

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
           call dsinv3(v1,m1,ep,ier)
           if (ier.eq.-1) then
c              write(*,*)   'echec inversion matrice information'
              istop=3
           endif

           do 3014 i=1,m1
              do 3015 j=i,m1
                 h_hess(i,j)=v1((j-1)*j/2+i)
 3015         continue
 3014      continue

           do 3016 i=2,m1
              do 3017 j=1,i-1
                 h_hess(i,j)=h_hess(j,i)
 3017         continue
 3016      continue

c======== fin 



           ep=10.d-10      
           call dsinv3(v,m,ep,ier)
           if (ier.eq.-1) then
c             write(6,103)          
c 103         format(1x,'echec inversion matrice information')

             istop=3
          endif

c            do 104 i=1,m
c               ii=i*(i+1)/2
c               write(6,*) 'variance parametre ',i,' = ',v(ii)
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
         call deriva3(b,m,vnonpen,rl,zero)
         
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

           do 108 i=2,m
              do 109 j=1,i-1
                 H_hess(i,j)=H_hess(j,i)
 109          continue
 108       continue

 110       continue
           return
       end
c======== fin marquard HJ


c===============================    MARQ98  HJ =========================

   	 subroutine marq983old(k0,b,m,ni,v,rl,ier,istop)

c  fu = mat3(rice des derivees secondes et premieres
c
c  istop: raison de l'arret
c  1: critere d'arret satisfait (prm=ca, vraisblce=cb, derivee=dd)
c  2: nb max d'iterations atteints
c  3: 1 mais echec inversion mat3(rice d'info (ier=1)
C  4: non amelioration vraisblce apres searpas
 
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************

         integer :: m,ni,nql,i,ii,nfmax,idpos,ier,istop,igrad,j
         integer :: ncount,id,jd,i0

         double precision :: da,dm,ga,tr,g0
         double precision :: ca,cb,epsa,epsb,rl
         double precision :: funcpa3,det, step,eps,epsd
         double precision :: vw,fi,maxt3,z
         double precision :: rl1,th,ep,dd

         double precision , dimension(npmax*(npmax+3)/2)::V
         double precision , dimension(npmax*(npmax+3)/2)::V1
         double precision , dimension(npmax*(npmax+3)/2)::vnonpen
         double precision , dimension(m)::b,delta
         double precision , dimension(npmax*(npmax+3)/2) :: fu
         double precision , dimension(2) :: k0
         double precision , dimension(m) :: bh,b1
         double precision , dimension(2) :: zero
         character*20 heure

c***** nmax
         integer :: nmax
         common /nmax/nmax

c*****dace2
      double precision , dimension(nsujetmax) ::t0,t1
      double precision :: ncsrl
      integer , dimension(nsujetmax) ::c
      integer, dimension(nsujetmax) :: nt0,nt1
      integer :: nsujet,nsim,nva,ndate,nst
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,nst
c*****dace3
      double precision :: pe
      integer :: effet,nz1,nz2,correl
      common /dace3/pe,effet,nz1,nz2,correl
c*****dace7
      double precision , dimension(npmax,npmax) :: I_hess,H_hess
      double precision , dimension(npmax,npmax) :: Hspl_hess
      double precision , dimension(npmax,1) ::PEN_deri
      common /dace7/PEN_deri,I_hess,H_hess,Hspl_hess
c***********************************************

c*****dace1 
c      double precision , dimension(ndatemax) ::date
c      double precision , dimension(-2:npmax) :: zi
c     common /dace1/date,zi


      zero(1)=0.d0
      zero(2)=0.d0
      id=0
      jd=0
      z=0.d0
      th=1.d-5
      eps=1.d-7
      epsa=1.d-4!1.d-3
      epsb=1.d-4!1.d-3
      epsd=1.d-3!1.d-2
      nfmax=m*(m+1)/2
      ca=epsa+1.d0
      cb=epsb+1.d0
      rl1=-1.d+10
      ni=0
      istop=0
      da=0.01d0
      dm=5.d0
c     nql=nbre de seconds membres pour chole
      nql=1

 10   continue
c     
c      write(6,*)
c      write(*,*)'** ITERATION :**', ni 
c         write(*,*)'** Parametres:'
c         write(*,*)(b(i),i=1,m)
c         write(6,*)'da',da,' ga=',ga
 
      z=0.d0
      i0=0 
      rl=funcpa3(b,m,i0,z,i0,z,k0)
      rl1=rl 

c      write(*,*)'log Vrais rl1,nsujet',rl1,nsujet
c      stop
      if (RL1.gt.0)stop
      
      call deriva3(b,m,v,rl,k0)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      do 144 i=1,m*(m+1)/2
c          ii=i*(i+1)/2
c          write(*,*) 'toutes derivees 2nds marq ',i,' = ',v(i)
c  144  continue 
c       ii=0
c       do 155 i=m*(m+1)/2+1,m*(m+3)/2
c          ii=ii+1
c       write(*,*) 'derivee premiere % parametre, marq',ii,' = ',v(i)
c  155  continue
cccccccccccccccccccccccc
c      stop
          dd = 0.d0
          do 13 i=m*(m+1)/2+1,m*(m+3)/2
c             write(*,*)'d',v(i)
             dd = dd + v(i)*v(i)
 13       continue 
          dd=dd/dabs(RL)

c          pause
c          write(6,*) (b(ii),ii=1,m)
          if (ca.lt.epsa.and.cb.lt.epsb.and.dd.lt.epsd) goto 100
c          if (ca.lt.epsa.and.cb.lt.epsb) goto 100
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
c          write(*,*)'avant dchole',fu(nfmax+1),m,nql,idpos
          call dchole3(fu,m,nql,idpos)
c          write(*,*)'apres dchole',fu(nfmax+1),m,nql,idpos
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
c              write(6,*) 'mat3(rice reguliere, idpos=',idpos
c              write(*,*)'ncount=',ncount,'da=',da,'ga=',ga
              do 2 i=1,m
                 delta(i)=fu(nfmax+i)
                 b1(i)=b(i)+delta(i)
2             continue
c              write(6,*) '** ** avant func.....B1',b1
c              write(6,*) '**'
c              write(6,*) '** ** avant func.....B:',b
c              write(6,*) '** ** avant func..',b1,m,id,z,jd,z,k0
              rl=funcpa3(b1,m,id,z,jd,z,k0)
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
           if(dabs(maxt3(delta,m)).eq.0.d0) then
              vw=1.D-5
           else
              vw=th/maxt3(delta,m)
           endif
c           
           step=dlog(1.5d0)
c           write(*,*)'vw ',vw,'th maxt ',th,maxt3(delta,m)
c           write(*,*)'delta',(delta(i),i=1,m)
c           pause
c           write(*,*) 'searpas'
           call SEARPAS3(vw,step,b,bh,m,delta,fi,eps,k0) 
c           write(*,*) 'apres searpas'
           rl=-fi
           IF(rl1.gt.rl) then
c              write(6,*) 'non amelioration vraisblce apres searpas'
c             write(6,*) 'ancienne rl1 =',rl1,' nouvelle rl =',rl
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

 800       cb=dabs(rl1-rl)
           ca=0.d0
           do 5 i=1,m
              ca=ca+delta(i)*delta(i)
 5         continue

c           write(6,*) 'ca =',ca,' cb =',cb,' dd =',dd
    
           do 6 i=1,m
              b(i)=b(i)+delta(i)
c              write(*,*)'delta',delta(i),i
 6         continue

           ni=ni+1
           if (ni.gt.nmax) then
              istop=2
c              write(6,*) 'nombre iteration max atteinte'
              goto 110
           end if
           goto 10
c********************
c     inversion mat3(rice d'informat3(ion
c********************
 100       continue
           istop=1

c           write(*,*)'ca,cb,dd final',ca,cb,dd


c================ pour les bandes de confiance/delta method

           call deriva3(b,m,v,rl,k0)
           do 303 i=1,(m*(m+3)/2)
              v1(i)=0.d0
 303       continue
           
           m1=m-nva-2*effet

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
           call DSINV3(v1,m1,ep,ier)
           if (ier.eq.-1) then
c              write(*,*)   'echec inversion mat3(rice informat3(ion'
              istop=3
           endif

           do 3014 i=1,m1
              do 3015 j=i,m1
                 h_hess(i,j)=v1((j-1)*j/2+i)
 3015         continue
 3014      continue

           do 3016 i=2,m1
              do 3017 j=1,i-1
                 h_hess(i,j)=h_hess(j,i)
 3017         continue
 3016      continue

c======== fin 



           ep=10.d-10      
           call DSINV3(v,m,ep,ier)
           if (ier.eq.-1) then
c             write(6,103)          
c 103         format(1x,'echec inversion mat3(rice informat3(ion')

             istop=3
          endif

c            do 104 i=1,m
c               ii=i*(i+1)/2
c               write(6,*) 'variance parametre ',i,' = ',v(ii)
c 104        continue
c            write(6,*)
c            ii=0
c            do 105 i=m*(m+1)/2+1,m*(m+3)/2
c               ii=ii+1
c               write(6,*) 'derivee premiere % parametre',ii,' = ',v(i)
c 105        continue

c========================================================
c pour pouvoir calculer la mat3(rice I_hess  = 
c (-) la hessienne sur la vraisemblances non penalisee (informat3(ion).
          
         ep=10.d-10
         call deriva3(b,m,vnonpen,rl,zero)
         
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

           do 108 i=2,m
              do 109 j=1,i-1
                 H_hess(i,j)=H_hess(j,i)
 109          continue
 108       continue

 110       continue
           return
       end
c======== fin marquard HJ





c===================================================================
c========================          FUNCPA       ====================
      double precision function funcpa3(b,np,id,thi,jd,thj,k0)

c      USE nr, ONLY : gammln

c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************

      integer :: nb,n,np,id,jd,i,j,k,vj,cptg,l
      integer :: ig,ip,choix
      integer , dimension(ngmax) :: cpt

      real , dimension(nptsmax) :: points,poids
      REAL :: gammln,essai

      double precision :: thi,thj,pe1,pe2,dnb,sum,hermite
      double precision :: som1,som2
      double precision :: res,vet,h1
      double precision :: int1,int2,int

      double precision , dimension(-2:npmax) ::the1,the2
      double precision , dimension(npmax) :: b,bh
      double precision , dimension(ngmax) :: integrale1,funcaux
      double precision , dimension(2) :: k0
      double precision , dimension(ngmax) :: sum1

      real, parameter :: pi = 3.1415926535
      double precision , dimension(3) :: vecD

c****** pour la maximisation avec marq98aux
      integer :: npaux,niaux,ieraux,istopaux
      DOUBLE PRECISION :: resaux
      double precision , dimension((npmax*(npmax+3)/2))::vaux
      double precision , dimension(npmax)::baux


c ************************* les commons :
c****** integration
      integer :: NDIM,MINPTS,MAXPTS,npts,type_int
      DOUBLE PRECISION :: EPSABS,EPSREL
      common /integration/NDIM,MINPTS,MAXPTS,npts,EPSABS,EPSREL,type_int

 
c*****ut1ut2
      double precision , dimension(ndatemax) :: dut1,dut2
      double precision , dimension(0:ndatemax) :: ut1,ut2
      common /ut1ut2/dut1,dut2,ut1,ut2
c*****pen1
      double precision ,dimension(npmax) ::m3m3,m2m2,m1m1,mmm,m3m2
      common /pen1/m3m3,m2m2,m1m1,mmm,m3m2
c*****pen2
      double precision ,dimension(npmax) ::m3m1,m3m,m2m1,m2m,m1m
      common /pen2/m3m1,m3m,m2m1,m2m,m1m

c*****mem1
      double precision ,dimension(ndatemax) :: mm3,mm2,mm1,mm
      common /mem1/mm3,mm2,mm1,mm
c*****mem2
      double precision ,dimension(ndatemax) :: im3,im2,im1,im
      common /mem2/im3,im2,im1,im
c*****dace1 
      double precision , dimension(ndatemax) ::date
      double precision , dimension(-2:npmax) :: zi
      common /dace1/date,zi

c*****dace2
      double precision , dimension(nsujetmax) ::t0,t1
      double precision :: ncsrl
      integer , dimension(nsujetmax) ::c
      integer, dimension(nsujetmax) :: nt0,nt1
      integer :: nsujet,nsim,nva,ndate,nst
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,nst
c*****dace4
      integer , dimension(nsujetmax) :: stra
      common /dace4/stra
c*****dace3
      double precision :: pe
      integer :: effet,nz1,nz2,correl
      common /dace3/pe,effet,nz1,nz2,correl

c*****ve1
      double precision , dimension(nsujetmax,nvarmax):: ve,ve2
      common /ve1/ve,ve2
c*****contrib
      integer :: ngexact,ng    !nb EXACT de gpes 
      common /contrib/ngexact,ng
c*****groupe
      integer , dimension(nsujetmax) :: g
      integer , dimension(ngmax) :: nig
      common /gpe/g,nig
c****** indicateur de troncature
      integer :: indictronq     ! =0 si donnees non tronquées
      common /troncature/indictronq
c*****auxig
      integer :: auxig
      common /auxig/auxig
c*******sigma2tau2rho
      double precision :: sigma2,tau2,rho,cov
      common /sigma2tau2/sigma2,tau2,rho,cov
c*****mij
      integer , dimension(ngmax) :: mi ! nb de dc dans gpe i
      common /mij/mi
c**** betaaux
      double precision , dimension(nvarmax) :: betaaux
      common /betaaux/betaaux
c*****invD
      double precision , dimension(2,2) :: invD
      common /invD/invD
c****** derivanal
      double precision , dimension(ngmax)::res1,res2,res3,res4
      double precision , dimension(ngmax)::res5,res6,res8
c      common /derivanal/res1,res2,res3,res4,res5,res6,res8

c****** u_tilde
      double precision  :: u_tilde,v_tilde
c      common /utilde/u_tilde,v_tilde
c******************

      restar = 0
      nf = 1   

c      write(*,*)'%%%%%%%% ENTREE DANS FUNCPA %%%%%%%%%%',effet,np,nst,K0
c      stop

c     on force sigma2 a 0.17 = 0.41 * 0.41
c      bh(np-nva-1)=0.41

      do 3 i=1,np
         bh(i)=b(i)
c         write(*,*)'%%%%%%%%%% ENTREE: bh ',bh(i),i
 3    continue 

      if (id.ne.0) bh(id)=bh(id)+thi 
      if (jd.ne.0) bh(jd)=bh(jd)+thj    
         
      if(nst.eq.2)then
c yas déb         n = nz1+2 +nz2+2
		 n = (nz1+2 +nz2+2)/nst
c yas fin		 
      else
         n = nz1+2
      endif
    
      do 4 i=1,n
         the1(i-3)=(bh(i))*(bh(i))
         j = n+i 
         if (nst.eq.2) then
            the2(i-3)=(bh(j))*(bh(j))
         endif
 4    continue

      if(effet.eq.1) then
 !terme de correlation entre 2 frailties/avec contraintes
         sigma2 = (bh(np-nva-1)*bh(np-nva-1)) ! variance intercept
         tau2 = (bh(np-nva)*bh(np-nva))  ! variance traitement * groupe
c         cov=bh(np-nva-2)!terme de covariance
         cov=dsqrt(sigma2*tau2)*
     &        (2.d0 * dexp(bh(np-nva-2))/(1.d0+dexp(bh(np-nva-2)))-1.d0)
c         write(*,*)'cov',cov
c         rho=bh(np-nva-2) !terme de correlation entre 2 frailties/sans contraintes
c         rho=2.d0 * dexp(bh(np-nva-2))/(1.d0+dexp(bh(np-nva-2)))-1.d0
      endif

c----------calcul de ut1(ti) et ut2(ti) ---------------------------
c     attention the(1)  sont en nz=1
c     donc en ti on a the(i)
      
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
 6       continue 
         ut1(i) = som1 +(the1(j-3)*im3(i))+(the1(j-2)*im2(i))
     &        +(the1(j-1)*im1(i))+(the1(j)*im(i))
         dut1(i) = (the1(j-3)*mm3(i))+(the1(j-2)*mm2(i))
     &        +(the1(j-1)*mm1(i))+(the1(j)*mm(i))

         if(nst.eq.2)then
            ut2(i) = som2 +(the2(j-3)*im3(i))+(the2(j-2)*im2(i))
     &           +(the2(j-1)*im1(i))+(the2(j)*im(i))
            dut2(i) = (the2(j-3)*mm3(i))+(the2(j-2)*mm2(i))
     &           +(the2(j-1)*mm1(i))+(the2(j)*mm(i))
         endif         
 8    continue

      i = n-2
      h1 = (zi(i)-zi(i-1))
      ut1(ndate)=som1+the1(i-4)+the1(i-3)+the1(i-2)+the1(i-1)
      dut1(ndate) = (4.d0*the1(i-1)/h1)    
c      write(*,*)'** dut1(ndate)',dut1(ndate),zi(i),zi(i-1),i,h1
c     & ,nva,np,n,effet
c
      if(nst.eq.2)then
         ut2(ndate)=som2+the1(i-4)+the2(i-3)+the2(i-2)+the2(i-1)
         dut2(ndate) = (4.d0*the2(i-1)/h1)
      endif
c         write(*,*)'** ut1',(ut1(i),i,i=0,ndate)!non
c         write(*,*)'** dut1',(dut1(i),i,i=0,ndate)
c         write(*,*)'** date **',(date(i),i,i=1,ndate)!ok idem 
c      if (nva.gt.0)then
c         write(*,*)'** dut2',dut2(ndate),the2(i-1),h1,ndate
c         stop
c         endif
c--------------------------------------------------------
c----------calcul de la vraisemblance ------------------
c---------------------------------------------------------
            
cc---- avec ou sans variable explicative  ------cc

      do 89 ig=1,ngexact!ng!
         res1(ig) = 0.d0
         res2(ig) = 0.d0
         cpt(ig)=0
 89   continue
c      stop
c*******************************************     
C----- sans effet aleatoire dans le modele
c*******************************************     

        if (effet.eq.0) then
           do 110 i=1,nsujet

c              if (nva.gt.0)then
c                 write(*,*)'***OK',nva,stra(i),c(i),i
c              endif
                   cpt(g(i))=cpt(g(i))+1
                   
                   if(nva.gt.0)then
                      vet = 0.d0   
                      do 91 j=1,nva
                         vet =vet + bh(np-nva+j)*ve(i,j)
c       write(*,*)'funcpa ve ******',ve(i,j),i,j
 91                   continue
                      vet = dexp(vet)
                   else
                      vet=1.d0
                   endif
                  
                   if((c(i).eq.1).and.(stra(i).eq.1))then
c                      if (nva.gt.0)  write(*,*)'res2 avant*',res2(g(i))
                     res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet)

c                      if (nva.gt.0) then
c       write(*,*)'res2 *',res2(g(i)),dlog(dut1(nt1(i))*vet),vet
c     &                    ,g(i),stra(i),dut1(nt1(i)),nt1(i),i
c       stop
c       endif

                  endif  

                   if((c(i).eq.1).and.(stra(i).eq.2))then
                      res2(g(i)) = res2(g(i))+dlog(dut2(nt1(i))*vet)
                      
c                      if (nva.gt.0)then
c                         write(*,*)'** dut2',dut2(1),dut2(ndate)
c(dut2(l),l,l=1,ndate)  
c                write(*,*)'res2 *',res2(g(i)),dlog(dut2(nt1(i))*vet),vet
c     &                        ,g(i),stra(i),dut2(nt1(i)),nt1(i),i
c                         stop
c                      endif
                  endif

                   if(stra(i).eq.1)then
                      res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet
     &                     -ut1(nt0(i))*vet 
                   endif

                   if(stra(i).eq.2)then
                      res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet
     &                     -ut2(nt0(i))*vet 
                   endif
           
 110           continue       
               res = 0.d0
           
               cptg = 0
          
C k indice les groupes
          do 115 k=1,ngexact 
             if(cpt(k).gt.0)then !nb de sujets dans un gpe=nig()
                             
                res = res-res1(k) 
     &               + res2(k) 
                cptg = cptg + 1  
c                if (nva.gt.0)then
c                  write(*,*)' **rescpt',res,res1(k),res2(k),cpt(k),k,nst
c                   stop
c                endif

             endif 
 115       continue
c           stop

        else !effet not =0
           

c*******************************************         
C-----avec deux effets aleatoires dans le modele et correl=0
c*********************************************


           if(effet.eq.1.and.correl.eq.0)then
              mi=0
              integrale1=0.d0

              do 116 k=1,nsujet
                 if(c(k).eq.1)then
                    mi(g(k))=mi(g(k))+1
c     write(*,*)'nb de deces mi', mi(g(k)),k,'(groupe=',g(k),')'
                 endif
 116          continue 

          do 106 k=1,nsujet
             if(nva.gt.0)then
                vet = 0.d0 
                do 196 ip=1,nva
                   vet =vet + bh(np-nva+ip)*ve(k,ip)
 196            continue
                vet = dexp(vet)
             else
                vet=1.d0
             endif

             if((c(k).eq.1).and.(stra(k).eq.1))then
                res2(g(k)) = res2(g(k))+dlog(dut1(nt1(k))*vet)
c     &                 ,nva
             endif  
             if((c(k).eq.1).and.(stra(k).eq.2))then
                res2(g(k)) = res2(g(k))+dlog(dut2(nt1(k))*vet)
             endif 
 
 106      continue 

c=========================================================================
c==== calcul des intégrales par transformat3(ion de LAPLACE pour chq gpe ===
c=========================================================================

          do 1706 ig=1,ng!ngexact
             res3(ig)=0.d0
             res4(ig)=0.d0
             res5(ig)=0.d0
             res6(ig)=0.d0
             auxig=ig
             baux(1)=0.05d0      !initialisation de u_tilde
             baux(2)=0.05d0      !initialisation de v_tilde
c             write(*,*)'***res2',res2(ig),ig
c======================================================================
c maximisation pour deux parametres, les 2 effets aleatoires : u et v 
c======================================================================
             npaux=2
             niaux=0
             ieraux=0
             istopaux=0
             vaux=0.d0!vecteur derivees 2nd et 1ERES
             funcaux=0.d0!vraisemblance
             resaux=0.d0!vraisemblance
             do 1916 ip=1,nva
                betaaux(ip)= bh(np-nva+ip)
 1916        continue


             
         call marq98AUX3(k0,baux,npaux,niaux,vaux,resaux,ieraux,
     &                  istopaux)

             u_tilde = baux(1)!u_tilde
             v_tilde = baux(2)!v_tilde

             if(effet.eq.1) then
c                write(*,*),'aux',baux(1),baux(2)
c                pause
                end if
                 do 1016 k=1,nsujet
                if(nva.gt.0.and.g(k).eq.ig)then
                   vet = 0.d0 
                   do 19126 ip=1,nva
                      vet =vet + bh(np-nva+ip)*ve(k,ip)
19126              continue
                   vet = dexp(vet)
                else
                   vet=1.d0
                endif

c              write(*,*),'res4 avant',res4(ig),ig
                if(g(k).eq.ig)then
                   if(c(k).eq.1)then
                     res3(ig) = res3(ig)
     &                +u_tilde+v_tilde*ve2(k,1)
                  endif

                  if(stra(k).eq.1)then
                      res4(ig) = res4(ig)
     &                 +ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))
                      res5(ig) = res5(ig)
     &        +ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                      res6(ig) = res6(ig)
     &     +ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2

                   endif
                  if(stra(k).eq.2)then
                      res4(ig) = res4(ig)
     &                 +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))
                      res5(ig) = res5(ig)
     &        +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                      res6(ig) = res6(ig)
     &     +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2

                   endif

                 endif 
 1016         continue 

c=====fin maximisation aux

             integrale1(ig)= -0.5d0*dlog(sigma2*tau2)
     &      -0.5d0*
     & dlog(res4(ig)*res6(ig)+res4(ig)/tau2+res6(ig)/sigma2
     &      +1.d0/(sigma2*tau2)-res5(ig)*res5(ig))!det(K"(b))
     &  +res2(ig)+res3(ig)-res4(ig)
     &-0.5d0*((u_tilde**2)/sigma2+(v_tilde**2)/tau2)!-k(b)

c          write(*,*)'*** integ1**',integrale1(ig),sigma2,tau2,
c     &res2(ig),res3(ig),res4(ig),res5(ig),res6(ig),u_tilde,v_tilde
c          stop
 1706     continue
c======= fin calcul de log integrale
c======================================================================


          res = 0.d0
          do 156 k=1,ng!ngexact  
             if(nig(k).gt.0)then
                if(indictronq.eq.0)then
                   res = res
c     &                  -4.d0*mi(ig)!correction faite dans integrant
c     &                  +res2(k) !+dlog(1.d0/(2.d0*dble(pi)))
     &                  +integrale1(k)!integrale1 donne le log de I directement
c          write(*,*)'*** res**',res,res2(k),integrale1(k)
c          stop
                else !troncature
                   write(*,*)'***TRAITER TRONCATURE**'
                   stop
                endif
             endif
 156      continue
c          write(*,*)'*** res**',res
c          stop

           endif                !fin boucle effet= 1 and correl = 0

c=======================================================================

c*******************************************         
C-----avec deux effets aleatoires dans le modele et correl=1
c*********************************************

           if(effet.eq.1.and.correl.eq.1)then

              mi=0
              integrale1=0.d0

              do 11 k=1,nsujet
                 if(c(k).eq.1)then
                    mi(g(k))=mi(g(k))+1
c             write(*,*)'nb de deces mi', mi(g(k)),k,'(groupe=',g(k),')'
                 endif
 11           continue 

          do 10 k=1,nsujet
             if(nva.gt.0)then
                vet = 0.d0 
                do 19 ip=1,nva
                   vet =vet + bh(np-nva+ip)*ve(k,ip)
c                   write(*,*)'ici vet**',vet,bh(np-nva+ip),ve(k,ip),k,ip
 19             continue
                vet = dexp(vet)
c                write(*,*)'**vet**',vet 
             else
                vet=1.d0
             endif

             if((c(k).eq.1).and.(stra(k).eq.1))then
                res2(g(k)) = res2(g(k))+dlog(dut1(nt1(k))*vet)
             endif  
             if((c(k).eq.1).and.(stra(k).eq.2))then
                res2(g(k)) = res2(g(k))+dlog(dut2(nt1(k))*vet)
             endif 
 
 10       continue 

c=========================================================================
c==== calcul des intégrales par transformat3(ion de LAPLACE pour chq gpe ===
c=========================================================================

          do 170 ig=1,ngexact!ng!
             res3(ig)=0.d0
             res4(ig)=0.d0
             res5(ig)=0.d0
             res6(ig)=0.d0
             auxig=ig
c     IF(IG.EQ.1)THEN
                baux(1)=0.05d0   !initialisation de u_tilde
                baux(2)=0.05d0   !initialisation de v_tilde
c     ENDIF
c======================================================================
c maximisation pour deux parametres, les 2 effets aleatoires : u et v 
c======================================================================
             npaux=2
             niaux=0
             vaux=0.d0!vecteur derivees 2nd et 1ERES
             funcaux=0.d0!vraisemblance
             resaux=0.d0!vraisemblance
             do 191 ip=1,nva
                betaaux(ip)= bh(np-nva+ip)
 191         continue
             
         call marq98AUX3(k0,baux,npaux,niaux,vaux,resaux,ieraux,
     &              istopaux)
             
             u_tilde = baux(1)!u_tilde(ig)
             v_tilde = baux(2)!v_tilde

             if(effet.eq.1) then
c                write(*,*),'aux',baux(1),baux(2)
c                pause
                end if
                 do 101 k=1,nsujet
                if(nva.gt.0.and.g(k).eq.ig)then
                   vet = 0.d0 
                   do 1912 ip=1,nva
                      vet =vet + bh(np-nva+ip)*ve(k,ip)
 1912              continue
                   vet = dexp(vet)
                else
                   vet=1.d0
                endif

c              write(*,*),'res4 avant',res4(ig),ig
                if(g(k).eq.ig)then
                   if(c(k).eq.1)then
                     res3(ig) = res3(ig)
     &                +u_tilde+v_tilde*ve2(k,1)
                  endif

                  if(stra(k).eq.1)then
                      res4(ig) = res4(ig)
     &                 +ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))
                      res5(ig) = res5(ig)
     &        +ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                      res6(ig) = res6(ig)
     &     +ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2

                   endif
                  if(stra(k).eq.2)then
                      res4(ig) = res4(ig)
     &                 +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))
                      res5(ig) = res5(ig)
     &        +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                      res6(ig) = res6(ig)
     &     +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2

                   endif

                 endif 
 101          continue 

c=====fin maximisation aux
             
c             cov=bh(np-nva-2)!sans contrainte
         cov=dsqrt(sigma2*tau2)*!avec contrainte
     &        (2.d0 * dexp(bh(np-nva-2))/(1.d0+dexp(bh(np-nva-2)))-1.d0)

             integrale1(ig)= -0.5d0*dlog(sigma2*tau2-cov*cov)!-0.5logdet(D)
     &      -0.5d0*
     & dlog((res4(ig)+tau2/(tau2*sigma2-cov**2))
     &      *(res6(ig)+sigma2/(tau2*sigma2-cov**2))
     &      -(res5(ig)-cov/(tau2*sigma2-cov**2))**2)!-0.5*log(det(ka2)
     &  +res2(ig)
     &+res3(ig)-res4(ig)
     &  -0.5d0*((u_tilde**2)/sigma2+(v_tilde**2)/tau2
     &            -2.d0*u_tilde*v_tilde*cov/(sigma2*tau2))
     &            /(1.d0-(cov**2)/(sigma2*tau2))         !-ka

c          write(*,*)'*** integ1**',integrale1(k),sigma2,tau2,
c     &res2(ig),res3(ig),res4(ig),res5(ig),res6(ig),u_tilde,v_tilde
c          stop
 170      continue
c======= fin calcul de log integrale
c          stop
c======================================================================


c======= fin calcul de integrale
c         print*,'** integrale1 **',(integrale1(ig),ig,ig=1,ngexact)
c          stop
c======================================================================


          res = 0.d0
          do 15 k=1,ngexact  !ng!
             if(nig(k).gt.0)then
                if(indictronq.eq.0)then
                   res = res
c     &                  +res2(k) !+dlog(1.d0/(2.d0*dble(pi)))
     &                  +integrale1(k)!integrale1 donne le log de I directement
c          write(*,*)'*** res**',res,res2(k),integrale1(k)

                else !troncature
                   write(*,*)'***TRAITER TRONCATURE**'
                   stop
                endif
             endif
 15   continue
c          write(*,*)'*** res**',res
c          stop
c             stop
      endif                     !fin boucle effet= 1 and correl = 1
      endif                     !fin boucle globale effet=0 

c----------calcul de la penalisation -------------------
             
             pe1 = 0.d0
             pe2 = 0.d0
             do 20 i=1,n-3
                pe1 = pe1+(the1(i-3)*the1(i-3)*m3m3(i))+(the1(i-2)
     &               *the1(i-2)*m2m2(i))+(the1(i-1)*the1(i-1)*m1m1(i))+(
     &               the1(i)*the1(i)*mmm(i))+(2.d0*the1(i-3)*the1(i-2)*
     &               m3m2(i))+(2.d0*the1(i-3)*the1(i-1)*m3m1(i))+(2.d0*
     &              the1(i-3)*the1(i)*m3m(i))+(2.d0*the1(i-2)*the1(i-1)*
     &          m2m1(i))+(2.d0*the1(i-2)*the1(i)*m2m(i))+(2.d0*the1(i-1)
     &    *the1(i)*m1m(i))

c                 write(*,*)'===PE1===',i,m3m3(i),m2m2(i),m1m1(i)

            if (nst.eq.2) then
               pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2)
     &              *the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+(
     &              the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)*
     &              m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0*
     &              the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)*
     &          m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1)
     &              *the2(i)*m1m(i))
            endif
 20      continue
         
c         print*,'ok6'
         if (nst.eq.2) then
            pe = k0(1)*pe1 + k0(2)*pe2 
c           write(*,*)'avant  penalisation res et pe1 et pe2',res,pe1,pe2
c           write(*,*)'vraisemblance penalisee',k0(1),k0(2)
          
         else
            pe = k0(1)*pe1
c           write(*,*)'avant  penalisation res et pe1 ',res,pe1
         endif 
         
         res = res - pe

c     write(*,*)'vraisemblance penalisee',pe,pe1,pe2
c     write(*,*)'vraisemblance penalisee',k0(1),k0(2)
        
         funcpa3= res 
         
c      print*,'*** res :',res,bh(id),id,np         
c      write(*,*)'%%%%%%%%%% SORTIE DANS FUNCPA %%%%%',res,pe,pe1,k0(1)
c      stop
        return
         end



c===================================   WEIGUI2    ============================

	subroutine weigui2(a,b,betau2,x)

        double precision ::a,b,x,u,v,betau2
        double precision ::uniran

        call RANDOM_number(u)
        v = (1.d0-u)
        x = (1.d0/b)*((-dexp(-betau2)*dlog(v))**(1.d0/a))

	return
	end

c***************************************************************

      DOUBLE PRECISION FUNCTION UNIRAN()
*
*     Random number generator(RCARRY), adapted from F. James
*     "A Review of Random Number Generators"
*      Comp. Phys. Comm. 60(1990), pp. 329-344.
*
      DOUBLE PRECISION SEEDS(24), TWOM24, CARRY, ONE
      PARAMETER ( ONE = 1, TWOM24 = ONE/16777216 )
      INTEGER I, J
      SAVE I, J, CARRY, SEEDS
      DATA I, J, CARRY / 24, 10, 0.0 /
      DATA SEEDS /
     & 0.8804418, 0.2694365, 0.0367681, 0.4068699, 0.4554052, 0.2880635,
     & 0.1463408, 0.2390333, 0.6407298, 0.1755283, 0.7132940, 0.4913043,
     & 0.2979918, 0.1396858, 0.3589528, 0.5254809, 0.9857749, 0.4612127,
     & 0.2196441, 0.7848351, 0.4096100, 0.9807353, 0.2689915, 0.5140357/
      UNIRAN = SEEDS(I) - SEEDS(J) - CARRY
      IF ( UNIRAN .LT. 0 ) THEN
         UNIRAN = UNIRAN + 1
         CARRY = TWOM24
      ELSE
         CARRY = 0
      ENDIF
      SEEDS(I) = UNIRAN
      I = 24 - MOD( 25-I, 24 )
      J = 24 - MOD( 25-J, 24 )
      END

c===============================================================
c======== INTEGRATION DOUBLE / genz ============================
c===============================================================

      SUBROUTINE HRMSYM( NDIM, NF, MINPTS, MAXPTS, FUNSUB, EPSABS, 
     &     EPSREL, RESTAR, RESULT, ABSERR, NEVAL, IFAIL, WORK)

c CALCUL D'UNE INTÉGRALE DE LA FORME :
c \int 1/\sqrt(2 pi) * exp (-x^2/2) *f(x) dx
c avec une quadrature gaussienne non adaptative (ie sans faire de transformat3(ions sur f(x))

****BEGIN PROLOGUE HRMSYM
*
****AUTHOR
*            Alan Genz, 
*            Department of mat3(hemat3(ics 
*            Washington State University
*            Pullman, WA 99164-3113, USA
*            Email: alangenz@wsu.edu
*
*      Reference: Genz, A., and Keister, B. (1996), Fully Symmetric 
*      Interpolatory Rules for Multiple Integrals over Infinite 
*      Regions with Gaussian Weight, J. Comp. Appl. mat3(h., 71, 299-309.
*
****KEYWORDS automat3(ic multidimensional integrator,
*            n-dimensional region (-infin, infin)^n 
*            Gaussian weight
****PURPOSE  The routine calculates an approximat3(ion to a given
*            vector of definite integrals
*
*      infin     infin 
*     I    ...  I     w(X)(F ,F ,...,F  ) DX(NDIM)...DX(2)DX(1),
*     -infin    -infin      1  2      NF
*
*       where F = F (X ,X ,...,X    ), I = 1,2,...,NF,
*              I   I  1  2      NDIM
*
*       w(X) = EXP(-( X(1)**2 + ... + X(NDIM)**2 )/2)/SQRT(2*PI)**NDIM.
*
****DESCRIPTION Computation of integrals over infinite regions with
*               Gaussian weight.
*
*   ON ENTRY
*
*     NDIM   Integer, number of variables, with 1 < NDIM <= 1000.
*     NF     Integer, number of components of the integral.
*     MINPTS Integer, minimum number of function evaluations.
*     MAXPTS Integer, maximum number of function evaluations.
*     FUNSUB Externally declared subroutine for computing
*            all components of the integrand at the given
*            evaluation point.
*            It must have parameters (NDIM,X,NF,FUNVLS)
*            Input parameters:
*              NDIM   Integer that defines the dimension of the
*                     integral.
*              X      Real array of dimension NDIM
*                     that defines the evaluation point.
*              NF Integer that defines the number of components of I.
*            Output parameter:
*              FUNVLS Real array of dimension NF
*                     that defines NF components of the integrand.
*
*     EPSABS Real, requested absolute accuracy.
*     EPSREL Real, requested relative accuracy.
*     RESTAR Integer.
*            If RESTAR = 0, this is the first attempt to compute
*            the integral.
*            If RESTAR = 1, then we restart a previous attempt.
*     WORK   Real array of working storage, must have dimension at
*            least 2*NF + (NF+1)*NUMSMS. NUMSMS is the number 
*            of NDIM partitions of the integers 0,...,RULE, after
*            RULE+1 approximat3(ions to the integrals have been computed.
*            The required value for NUMSMS will never exceed 10000.
*
*   ON RETURN
*
*     RESULT Real array of dimension NF.
*            Approximat3(ions to all components of the integral.
*     ABSERR Real array of dimension NF.
*            Estimat3(es of absolute accuracies.
*     NEVAL  Integer, number of function evaluations used.
*     IFAIL  Integer.
*            IFAIL = 0 for normal exit, when 
*              ABSERR(K) <= MAX( EPSABS, ABS(RESULT(K))*EPSREL ) for
*              all K, 0 < K <= NF, with <= MAXPTS function values. 
*            IFAIL = 1 if MAXPTS was too small to obtain the required 
*              accuracy. In this case values of RESULT are returned
*              with estimat3(ed absolute accuracies ABSERR.
*     WORK   Real array of working storage, that contains informat3(ion
*            that might be needed for subsequent calls of HRMSYM. This
*            array should not be modified. 
*
****ROUTINES CALLED 
****END PROLOGUE HRMSYM
*
*   Global variables.
*

      EXTERNAL FUNSUB
      INTEGER NDIM, NF, MINPTS, MAXPTS, RESTAR
      INTEGER NEVAL, IFAIL
      DOUBLE PRECISION EPSABS, EPSREL
      DOUBLE PRECISION RESULT(*), ABSERR(*), WORK(*)
*
*   Local variables.
*
      INTEGER I, MAXRUL, RULE, MNRULE, INTCLS, NUMSMS
      PARAMETER ( MAXRUL = 25 )
      DOUBLE PRECISION DIFFER, CLTOTL, RNTCLS, EPOWER
      PARAMETER ( EPOWER = 1.5 )
      SAVE RULE, MNRULE, CLTOTL
      IFAIL = 1
      IF ( NDIM .LE. 4 ) THEN 
*
*     Call product Gauss-Hermite Rule
*
         CALL HERMIT( NDIM, NF, MINPTS, MAXPTS, FUNSUB, EPSABS, EPSREL, 
     &                RESTAR, RESULT, ABSERR, NEVAL, IFAIL, WORK )
      ELSE
         IF ( RESTAR .EQ. 0 ) THEN
            MNRULE = -1
            RULE = 0
            CLTOTL = 0
            DO I = 1, NF
               WORK(I) = 0
               WORK(I+NF) = 0
            END DO
         END IF
         NEVAL = 0
         DO WHILE ( NEVAL .LE. MAXPTS .AND. RULE .LE. MAXRUL 
     &        .AND. ( IFAIL .GT. 0 .OR. NEVAL .LT. MINPTS ) ) 
            CALL HRMTRL( NDIM, NF, FUNSUB, MNRULE, RULE, RESULT, 
     &                   INTCLS, WORK(2*NF+1), WORK(3*NF+1), NUMSMS ) 
            RNTCLS = REAL( INTCLS )**EPOWER 
            CLTOTL = CLTOTL + RNTCLS
            DO I = 1, NF
               DIFFER = ( RESULT(I) - WORK(I) )/CLTOTL
               WORK(I) = WORK(I) + RNTCLS*DIFFER 
               WORK(NF+I) = ( CLTOTL - RNTCLS )
     &                     *( WORK(NF+I)/CLTOTL + RNTCLS*DIFFER**2 )
            END DO
            IFAIL = 0
            DO I = 1,NF
               IF ( RULE .GT. 0 ) THEN
                  ABSERR(I) = SQRT( WORK(NF+I) )
               ELSE
                  ABSERR(I) = ABS( RESULT(I) )
               END IF
               IF ( ABSERR(I) .GT. MAX(EPSABS,EPSREL*ABS(RESULT(I))) )
     &              IFAIL = 1
            END DO
            NEVAL = NEVAL + INTCLS
            RULE = RULE + 1
         END DO
      END IF
*
****END HRMSYM
*
      END
      SUBROUTINE HRMTRL( S, N, F, MINORD, MAXORD, INTVAL, INTCLS, 
     &                   WORK, FULSMS, NUMSMS )
*  MULTIDIMENSIONAL FULLY SYMMETRIC RULE INTEGRATION SUBROUTINE
*
*   THIS SUBROUTINE COMPUTES A SEQUENCE OF FULLY SYMMETRIC RULE
*   APPROXImat3(IONS TO A FULLY SYMMETRIC MULTIPLE INTEGRAL. WRITTEN BY 
*     Alan Genz
*     Department of mat3(hemat3(ics
*     Washington State University
*     Pullman, Washington 99164-3113  USA
*     Telephone: 509-335-2131 
*     Electronic Mail: alangenz@wsu.edu
*     Fax: 509-335-1188 
*
***************  PARAMETERS FOR HRMTRL  ********************************
******INPUT PARAMETERS
*  S       INTEGER number of variables, with 0 < S <= 1000.
*  N       Integer, number of components of the integral.
*  F       EXTERNALly declared user defined integrand subroutine.
*          it must have parameters (S,X,N,F), where X is a REAL S-array
*          and F is a real N-array.
*  MINORD  INTEGER minimum order parameter.  On entry MINORD specifies
*          the current highest order approximat3(ion to the integral,
*          available in the array INTVAL.  For the first call
*          MINORD should be set to -1. Otherwise a previous call is
*          assumed that computed INTVAL. On exit MINORD is set to MAXORD.
*  MAXORD  INTEGER maximum order parameter, must be greater than MINORD
*          and not exceed 25. The subroutine computes approximat3(ions of 
*          polynomial degree 2*MAXORD+1.
*******OUTPUT PARAMETERS
*  INTVAL  REAL array of length N.  Upon successful exit
*          INTVAL(1),..., INTVAL(N) are approximat3(ions to the components
*          of the integral.  These are all approximat3(ions of polynomial 
*          degree 2*MAXORD+1.
*  INTCLS  INTEGER number of F values needed for INTVAL
*  WORK    REAL working storage array. 
*  FULSMS  REAL working storage array with dimension (N+1,*). On exit
*          FULSMS(I,J) contains the fully symmetric basic rule sum
*          indexed by the jth s-partition of the integers 0,...,MAXORD, 
*          for the Ith component of the integrand.
*          FULSMS(N+1,J) contains number of points for the fully 
*          symmetric basic rule sum indexed by the Jth S-partition of 
*          the integers 0,...,MAXORD.
*  NUMSMS  INTEGER number of S-partitions of the integers 0,...,MAXORD.
************************************************************************
      EXTERNAL F
      INTEGER S, N, MAXDIM, MINORD, MAXORD, MAXRDM, NUMSMS
      PARAMETER ( MAXDIM = 1000, MAXRDM = 25 ) 
      DOUBLE PRECISION X(MAXDIM), FULWGT, WEIGHT, 
     &                 INTVAL(N), FULSMS(N+1,*), WORK(*) 
      INTEGER D, I, L, MODOFM, M(MAXDIM), K(MAXDIM), INTCLS, PRTCNT
      D = MINORD + 1
      INTCLS = 0
      IF ( D .EQ. 0 ) THEN
         DO I = 1,N
            INTVAL(I) = 0
         END DO
      END IF
*
****  Begin loop for each D
*      for each D find all distinct partitions M with |M| <= D
*
      DO WHILE( D .LE. MIN ( MAXORD, MAXRDM ) )
         PRTCNT = 0
         CALL NXPART( PRTCNT, S, M, MODOFM )
         DO WHILE( MODOFM .LE. D )
*     
****  Calculate updated weight for partition M alnd 
****     fully symmetric sums ( when necessary )
*
            FULWGT = WEIGHT( S, X, M, K, MODOFM, D )
            IF ( D .EQ. MODOFM ) THEN
               DO I = 1,N
                  FULSMS(I,PRTCNT) = 0
               END DO
               FULSMS(N+1,PRTCNT) = 0
            END IF
            IF ( FULSMS(N+1,PRTCNT) .EQ. 0 .AND. FULWGT .NE. 0 ) THEN
               CALL FULSMH( S, M, N, F, FULSMS(1,PRTCNT), X, WORK ) 
     &                      
               INTCLS = INTCLS + FULSMS(N+1,PRTCNT)
            END IF
            DO I = 1,N 
               INTVAL(I) = INTVAL(I) + FULWGT*FULSMS(I,PRTCNT)
            END DO
            CALL NXPART( PRTCNT, S, M, MODOFM )
         END DO
*     
****  End loop for each D
*
         D = D + 1
      END DO
      MINORD = MAXORD
      NUMSMS = PRTCNT - 1
      END
*
      SUBROUTINE FULSMH( S, M, N, F, FULSMS, X, FUNVAL )
*
****  To compute fully symmetric basic rule sums
*
      INTEGER S, M(*), SUMCLS, IX, LX, I,L, MI, ML, IL, N, MX
      DOUBLE PRECISION X(*), INTWGT, FULSMS(*), FUNVAL(*)
*
*        Generators for 1 + 2 + 6 + 10 + 16 = 35 point 
*        degree 51 rule, with degree 1, 5, 15 and 29 imbedded rules.     
*
      PARAMETER ( MX = 17 )
      DOUBLE PRECISION G(0:MX)
      SAVE G
      DATA G( 0),G( 1)/0D0                  , 0.17320508075688773D1/
      DATA G( 2),G( 3)/0.41849560176727319D1, 0.74109534999454084D0/
      DATA G( 4),G( 5)/0.28612795760570581D1, 0.63633944943363700D1/
      DATA G( 6),G( 7)/0.12304236340273060D1, 0.51870160399136561D1/
      DATA G( 8),G( 9)/0.25960831150492022D1, 0.32053337944991945D1/
      DATA G(10),G(11)/0.90169397898903025D1, 0.24899229757996061D0/
      DATA G(12),G(13)/0.79807717985905609D1, 0.22336260616769417D1/
      DATA G(14),G(15)/0.71221067008046167D1, 0.36353185190372782D1/
      DATA G(16),G(17)/0.56981777684881096D1, 0.47364330859522971D1/
      INTWGT = 1
      DO I = 1, S
        IF ( M(I) .NE. 0 ) INTWGT = INTWGT/2
      END DO
      SUMCLS = 0
      DO I = 1,N
         FULSMS(I) = 0
      END DO
*
********  Compute centrally symmetric sum for permutation M
*
 10   DO I = 1, S
         X(I) = -G(M(I))
      END DO
 20   SUMCLS = SUMCLS + 1
      CALL F( S, X, N, FUNVAL )
      DO I = 1,N
         FULSMS(I) = FULSMS(I) + INTWGT*FUNVAL(I)
      END DO
      DO I = 1, S
         X(I) = -X(I)
         IF ( X(I) .GT. 0 ) GO TO 20
      END DO
*
********  END Integration loop for M
*
*
********  Find next distinct permutation of M and loop back
*          to compute next centrally symmetric sum
*
      DO I = 2, S
         IF ( M(I-1) .GT. M(I) ) THEN
            MI = M(I)
            IX = I - 1
            IF ( I .GT. 2 ) THEN
               DO L = 1, IX/2
                  ML = M(L)
                  IL = I - L
                  M(L) = M(IL)
                  M(IL) = ML
                  IF ( ML .LE. MI ) IX = IX - 1
                  IF ( M(L) .GT. MI ) LX = L
               END DO
               IF ( M(IX) .LE. MI ) IX = LX
            END IF
            M(I) = M(IX)
            M(IX) = MI
            GO TO 10
         END IF
      END DO
*
****  END Loop for permutations of M and associated sums
*
*
**** Restore original order to M.
*
      DO I = 1, S/2
         MI = M(I)
         M(I) = M(S-I+1)
         M(S-I+1) = MI
      END DO
      FULSMS(N+1) = SUMCLS
      END
*
      DOUBLE PRECISION FUNCTION WEIGHT( S, INTRPS, M,K, MODOFM, D )
*
****  Subroutine to update weight for partition m
*
      INTEGER S, M(S), K(S), I, L, D, MAXRDM, NZRMAX, MODOFM
      PARAMETER ( MAXRDM = 25 )
      DOUBLE PRECISION INTRPS(S), MOMPRD(0:MAXRDM,0:MAXRDM), MOMNKN,
     &     INTMPA, INTMPB
      SAVE MOMPRD
*
*     Modified moments and generators for 1 + 2 + 6 + 10 + 16 = 35 point 
*        degree 51 rule, with degree 1, 5 and 19 imbedded rules.     
*
      PARAMETER ( NZRMAX = 17 )
      DOUBLE PRECISION A(0:NZRMAX), G(0:NZRMAX)
      DATA A /  2*1D0, 0D0, 6D0, 
     &     -0.48378475125832451D2, 3*0D0, 34020D0, 
     &     -0.98606453173677489D6, 5*0D0, 0.12912054173706603D13, 
     &     -0.11268664521456168D15, 0.29248520348796280D16 /  
      DATA G( 0),G( 1)/0D0                  , 0.17320508075688773D1/
      DATA G( 2),G( 3)/0.41849560176727319D1, 0.74109534999454084D0/
      DATA G( 4),G( 5)/0.28612795760570581D1, 0.63633944943363700D1/
      DATA G( 6),G( 7)/0.12304236340273060D1, 0.51870160399136561D1/
      DATA G( 8),G( 9)/0.25960831150492022D1, 0.32053337944991945D1/
      DATA G(10),G(11)/0.90169397898903025D1, 0.24899229757996061D0/
      DATA G(12),G(13)/0.79807717985905609D1, 0.22336260616769417D1/
      DATA G(14),G(15)/0.71221067008046167D1, 0.36353185190372782D1/
      DATA G(16),G(17)/0.56981777684881096D1, 0.47364330859522971D1/
      DATA MOMPRD(0,0) / 0D0 /
      IF ( MOMPRD(0,0) .EQ. 0 ) THEN
*
****  Calculate moments 
*
         DO L = 0, MAXRDM 
            DO I = 0, MAXRDM 
               MOMPRD(L,I) = 0
            END DO
         END DO
         MOMPRD(0,0) = A(0)
         DO L = 0, NZRMAX
            MOMNKN = 1
            DO I = 1, NZRMAX
               IF ( I .LE. L ) THEN
                  MOMNKN = MOMNKN*( G(L)**2 - G(I-1)**2 )
               ELSE
                  MOMNKN = MOMNKN*( G(L)**2 - G(I)**2 )
               END IF
               IF ( I .GE. L ) MOMPRD(L,I) = A(I)/MOMNKN
            END DO
         END DO
      END IF
*
*     Determine Updated Weight Contribution
*
      DO I = 2,S
        INTRPS(I) = 0
        K(I) = M(I)
      END DO
      K(1) = D - MODOFM + M(1)
 10   INTRPS(1) = MOMPRD( M(1), K(1) )
      DO I = 2, S
        INTRPS(I) = INTRPS(I) + MOMPRD( M(I), K(I) )*INTRPS(I-1)
        INTRPS(I-1) = 0
        K(1) = K(1) - 1
        K(I) = K(I) + 1
        IF ( K(1) .GE. M(1) ) GO TO 10
        K(1) = K(1) + K(I) - M(I)
        K(I) = M(I)
      END DO
      WEIGHT = INTRPS(S)
      END
      SUBROUTINE NXPART( PRTCNT, S, M, MODOFM )
*
****  SUBROUTINE TO COMPUTE THE NEXT S PARTITION
*
      INTEGER S, M(S), PRTCNT, MODOFM, I, MSUM, L
      IF ( PRTCNT .EQ. 0 ) THEN
         DO I = 1, S
            M(I) = 0
         END DO
         PRTCNT = 1
         MODOFM = 0
      ELSE
         PRTCNT = PRTCNT + 1
         MSUM = M(1)
         DO I = 2, S
            MSUM = MSUM + M(I)
            IF ( M(1) .LE. M(I) + 1 ) THEN
               M(I) = 0
            ELSE
               M(1) = MSUM - (I-1)*(M(I)+1)
               DO L = 2, I
                  M(L) = M(I) + 1
               END DO
               RETURN
            END IF
         END DO
         M(1) = MSUM + 1
         MODOFM = M(1)
      END IF
      END
*
      SUBROUTINE MLTRUL( NDIM, NUMFUN, FUNSUB, NP, POINT, WEIGHT, 
     &                   INTVAL, FUNS, X, IC )
*
*     Computes product integration rule
*
      EXTERNAL FUNSUB
      INTEGER NDIM, NP, NUMFUN, ICI, I
      DOUBLE PRECISION POINT(*), WEIGHT(*), INTVAL(*)
      DOUBLE PRECISION WTPROD, X(*), FUNS(*), IC(*)
      DO I = 1, NDIM
         IC(I) = 1
      END DO
      DO I = 1, NUMFUN
         INTVAL(I) = 0
      END DO
 10   WTPROD = 1
      DO I = 1, NDIM
         ICI = IC(I)
         X(I) = POINT(ICI)
         WTPROD = WTPROD*WEIGHT(ICI)
      END DO
      CALL FUNSUB( NDIM, X, NUMFUN, FUNS ) 
      DO I = 1, NUMFUN
         INTVAL(I) = INTVAL(I) + WTPROD*FUNS(I)
      END DO
      DO I = 1, NDIM
         IC(I) = IC(I) + 1
         IF ( IC(I) .LE. NP ) GO TO 10
         IC(I) = 1
      END DO
      END
*
      SUBROUTINE HERMIT( NDIM, NUMFUN, MINPTS, MAXPTS, FUNSUB,
     &     EPSABS, EPSREL, RESTAR, RESULT, ABSERR, NEVAL, IFAIL, WORK ) 
****BEGIN PROLOGUE HERMIT
****CATEGORY NO. H2B1A1
****AUTHOR
*
*     Alan Genz
*     Department of mat3(hemat3(ics
*     Washington State University
*     Pullman, Washington 99164-3113  USA
*     Telephone: 509-335-2131 
*     Electronic Mail: alangenz@wsu.edu
*     Fax: 509-335-1188 
*
****KEYWORDS automat3(ic multidimensional integrator, Gaussian weight
*            n-dimensional region (-infin, infin)^n 
****PURPOSE  The routine calculates an approximat3(ion to a vector of
*            definite integrals using Gauss-Hermite product rules.
*
*      infin     infin
*     I    ...  I       w(X) (F ,F ,...,F      ) DX(NDIM)...DX(2)DX(1),
*     -infin    -infin         1  2      NUMFUN
*
*       where F = F (X ,X ,...,X    ), I = 1,2,...,NUMFUN,
*              I   I  1  2      NDIM
*
*       w(X) = EXP(-( X(1)**2 + ... + X(NDIM)**2 )/2)/SQRT(2*PI)**NDIM.
*
****DESCRIPTION Computation of integrals over infinite regions with
*               Gaussian weight using Gauss-Hermite product rules.
*
*   ON ENTRY
*
*     NDIM   Integer, number of variables.
*     NUMFUN Integer, number of components of the integral.
*     MINPTS Integer, minimum number of function evaluations.
*     MAXPTS Integer, maximum number of function evaluations.
*     FUNSUB Externally declared subroutine for computing components of 
*            the integrand at the given evaluation point.
*            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
*            Input parameters:
*              NDIM   Integer, number of dimensions for the integral.
*              X      Real NDIM-array for the evaluation point.
*              NUMFUN Integer, number of components of I.
*            Output parameter:
*              FUNVLS Real NUMFUN-array for integrand components.
*
*     EPSABS Real, requested absolute accuracy.
*     EPSREL Real, requested relative accuracy.
*     RESTAR Integer.
*            If RESTAR = 0, this is the first attempt to compute
*            the integral.
*            If RESTAR = 1, then we restart a previous attempt.
*     WORK   Real array of working storage, must have dimensions at
*            least 2*NUMFUN+2*NDIM
*
*   ON RETURN
*
*     RESULT Real NUMFUN-array of approximat3(ions to all components 
*            of the integral.
*     ABSERR Real NUMFUN-array of estimat3(es of absolute accuracies.
*     NEVAL  Integer, number of function evaluations used.
*     IFAIL  Integer.
*            IFAIL = 0 for normal exit, when 
*              ABSERR(K) <= MAX( EPSABS, ABS(RESULT(K))*EPSREL ) for
*              all K, 0 < K <= NUMFUN, with <= MAXPTS function values. 
*            IFAIL = 1 if MAXPTS was too small to obtain the required 
*              accuracy. In this case values of RESULT are returned
*              with estimat3(ed absolute accuracies ABSERR.
*
****ROUTINES CALLED MLTRUL
****END PROLOGUE HERMIT
*
*   Global variables.
*
      EXTERNAL FUNSUB
      INTEGER NDIM, NUMFUN, MINPTS, MAXPTS, RESTAR
      INTEGER NEVAL, IFAIL
      DOUBLE PRECISION EPSABS, EPSREL
      DOUBLE PRECISION RESULT(*), ABSERR(*), WORK(*)
*
*   Local variables.
*
      INTEGER I, K, MAXRUL, RULE
      PARAMETER ( MAXRUL = 50 )
      DOUBLE PRECISION POINT(MAXRUL), WEIGHT(MAXRUL)
      SAVE RULE
      DOUBLE PRECISION W(25,50), T(25,50)
      SAVE W, T
      IF ( RESTAR .EQ. 0 ) RULE = 1
      NEVAL = 0
 10   IF ( NEVAL + RULE**NDIM .LE. MAXPTS .AND. RULE .LT. MAXRUL ) THEN
         DO I = 1, RULE/2
            POINT( I) = -T(I,RULE)
            WEIGHT(I) =  W(I,RULE)
            POINT( RULE-I+1) = T(I,RULE)
            WEIGHT(RULE-I+1) = W(I,RULE)
         END DO
         IF ( MOD( RULE, 2 ) .EQ. 1 ) THEN
            POINT(  RULE/2 + 1 ) = 0
            WEIGHT( RULE/2 + 1 ) = W( RULE/2 + 1, RULE ) 
         END IF
         CALL MLTRUL( NDIM, NUMFUN, FUNSUB, RULE, POINT, WEIGHT, 
     &        RESULT, WORK, WORK(NUMFUN+1), WORK(NUMFUN+NDIM+1) )
         NEVAL = NEVAL + RULE**NDIM
         IFAIL = 0
         DO I = 1,NUMFUN
            IF ( RULE .GT. 1 ) THEN
               ABSERR(I) = ABS( RESULT(I) - WORK(2*NDIM+NUMFUN+I)  )
            ELSE
               ABSERR(I) = ABS( RESULT(I) )
            ENDIF
            WORK(2*NDIM+NUMFUN+I) = RESULT(I)
            IF ( ABSERR(I) .GT. MAX( EPSABS, EPSREL*ABS(RESULT(I)) ) )
     &           IFAIL = 1
         END DO
         RULE = RULE + 1
         IF ( IFAIL .GT. 0 .OR. NEVAL .LT. MINPTS ) GO TO 10 
      END IF
*
*     Gauss Hermite Weights and Points, N = 1,50
*
      DATA   W(1, 1), T(1, 1) / 1D0, 0D0 /
      DATA ( W(I, 2), T(I, 2), I = 1, 1) /
     &  0.5000000000000001D+00, 0.1000000000000000D+01/
      DATA ( W(I, 3), T(I, 3), I = 1, 2) /
     &  0.1666666666666667D+00, 0.1732050807568877D+01,
     &  0.6666666666666664D+00, 0.1107367643833737D-15/
      DATA ( W(I, 4), T(I, 4), I = 1, 2) /
     &  0.4587585476806855D-01, 0.2334414218338977D+01,
     &  0.4541241452319317D+00, 0.7419637843027258D+00/
      DATA ( W(I, 5), T(I, 5), I = 1, 3) /
     &  0.1125741132772071D-01, 0.2856970013872805D+01,
     &  0.2220759220056126D+00, 0.1355626179974265D+01,
     &  0.5333333333333342D+00, 0.9386691848789097D-16/
      DATA ( W(I, 6), T(I, 6), I = 1, 3) /
     &  0.2555784402056243D-02, 0.3324257433552119D+01,
     &  0.8861574604191447D-01, 0.1889175877753710D+01,
     &  0.4088284695560291D+00, 0.6167065901925933D+00/
      DATA ( W(I, 7), T(I, 7), I = 1, 4) /
     &  0.5482688559722184D-03, 0.3750439717725742D+01,
     &  0.3075712396758645D-01, 0.2366759410734542D+01,
     &  0.2401231786050126D+00, 0.1154405394739968D+01,
     &  0.4571428571428575D+00, 0.2669848554723344D-16/
      DATA ( W(I, 8), T(I, 8), I = 1, 4) /
     &  0.1126145383753679D-03, 0.4144547186125893D+01,
     &  0.9635220120788268D-02, 0.2802485861287542D+01,
     &  0.1172399076617590D+00, 0.1636519042435109D+01,
     &  0.3730122576790775D+00, 0.5390798113513754D+00/
      DATA ( W(I, 9), T(I, 9), I = 1, 5) /
     &  0.2234584400774664D-04, 0.4512745863399781D+01,
     &  0.2789141321231769D-02, 0.3205429002856470D+01,
     &  0.4991640676521780D-01, 0.2076847978677829D+01,
     &  0.2440975028949394D+00, 0.1023255663789133D+01,
     &  0.4063492063492066D+00, 0.0000000000000000D+00/
      DATA ( W(I,10), T(I,10), I = 1, 5) /
     &  0.4310652630718282D-05, 0.4859462828332311D+01,
     &  0.7580709343122187D-03, 0.3581823483551927D+01,
     &  0.1911158050077027D-01, 0.2484325841638954D+01,
     &  0.1354837029802680D+00, 0.1465989094391158D+01,
     &  0.3446423349320194D+00, 0.4849357075154977D+00/
      DATA ( W(I,11), T(I,11), I = 1, 6) /
     &  0.8121849790214922D-06, 0.5188001224374871D+01,
     &  0.1956719302712241D-03, 0.3936166607129977D+01,
     &  0.6720285235537304D-02, 0.2865123160643646D+01,
     &  0.6613874607105794D-01, 0.1876035020154847D+01,
     &  0.2422402998739701D+00, 0.9288689973810635D+00,
     &  0.3694083694083690D+00, 0.0000000000000000D+00/
      DATA ( W(I,12), T(I,12), I = 1, 6) /
     &  0.1499927167637166D-06, 0.5500901704467746D+01,
     &  0.4837184922590630D-04, 0.4271825847932281D+01,
     &  0.2203380687533207D-02, 0.3223709828770096D+01,
     &  0.2911668791236414D-01, 0.2259464451000800D+01,
     &  0.1469670480453302D+00, 0.1340375197151617D+01,
     &  0.3216643615128298D+00, 0.4444030019441390D+00/
      DATA ( W(I,13), T(I,13), I = 1, 7) /
     &  0.2722627642805887D-07, 0.5800167252386502D+01,
     &  0.1152659652733391D-04, 0.4591398448936520D+01,
     &  0.6812363504429268D-03, 0.3563444380281636D+01,
     &  0.1177056050599653D-01, 0.2620689973432215D+01,
     &  0.7916895586044999D-01, 0.1725418379588239D+01,
     &  0.2378715229641365D+00, 0.8566794935194499D+00,
     &  0.3409923409923412D+00, 0.2011511664336819D-15/
      DATA ( W(I,14), T(I,14), I = 1, 7) /
     &  0.4868161257748367D-08, 0.6087409546901291D+01,
     &  0.2660991344067620D-05, 0.4896936397345567D+01,
     &  0.2003395537607445D-03, 0.3886924575059772D+01,
     &  0.4428919106947401D-02, 0.2963036579838668D+01,
     &  0.3865010882425336D-01, 0.2088344745701943D+01,
     &  0.1540833398425136D+00, 0.1242688955485464D+01,
     &  0.3026346268130198D+00, 0.4125904579546022D+00/
      DATA ( W(I,15), T(I,15), I = 1, 8) /
     &  0.8589649899633300D-09, 0.6363947888829836D+01,
     &  0.5975419597920602D-06, 0.5190093591304780D+01,
     &  0.5642146405189029D-04, 0.4196207711269018D+01,
     &  0.1567357503549958D-02, 0.3289082424398766D+01,
     &  0.1736577449213763D-01, 0.2432436827009758D+01,
     &  0.8941779539984458D-01, 0.1606710069028730D+01,
     &  0.2324622936097322D+00, 0.7991290683245483D+00,
     &  0.3182595182595181D+00, 0.0000000000000000D+00/
      DATA ( W(I,16), T(I,16), I = 1, 8) /
     &  0.1497814723161838D-09, 0.6630878198393126D+01,
     &  0.1309473216286842D-06, 0.5472225705949343D+01,
     &  0.1530003216248727D-04, 0.4492955302520013D+01,
     &  0.5259849265739089D-03, 0.3600873624171548D+01,
     &  0.7266937601184742D-02, 0.2760245047630703D+01,
     &  0.4728475235401395D-01, 0.1951980345716333D+01,
     &  0.1583383727509496D+00, 0.1163829100554964D+01,
     &  0.2865685212380120D+00, 0.3867606045005573D+00/
      DATA ( W(I,17), T(I,17), I = 1, 9) /
     &  0.2584314919374932D-10, 0.6889122439895331D+01,
     &  0.2808016117930569D-07, 0.5744460078659410D+01,
     &  0.4012679447979839D-05, 0.4778531589629983D+01,
     &  0.1684914315513387D-03, 0.3900065717198010D+01,
     &  0.2858946062284621D-02, 0.3073797175328194D+01,
     &  0.2308665702571097D-01, 0.2281019440252989D+01,
     &  0.9740637116272111D-01, 0.1509883307796740D+01,
     &  0.2267063084689769D+00, 0.7518426007038956D+00,
     &  0.2995383701266057D+00, 0.0000000000000000D+00/
      DATA ( W(I,18), T(I,18), I = 1, 9) /
     &  0.4416588769358736D-11, 0.7139464849146476D+01,
     &  0.5905488478836554D-08, 0.6007745911359599D+01,
     &  0.1021552397636983D-05, 0.5054072685442739D+01,
     &  0.5179896144116204D-04, 0.4188020231629400D+01,
     &  0.1065484796291652D-02, 0.3374736535778089D+01,
     &  0.1051651775194131D-01, 0.2595833688911239D+01,
     &  0.5489663248022256D-01, 0.1839779921508646D+01,
     &  0.1606853038935128D+00, 0.1098395518091501D+01,
     &  0.2727832346542882D+00, 0.3652457555076979D+00/
      DATA ( W(I,19), T(I,19), I = 1,10) /
     &  0.7482830054057162D-12, 0.7382579024030434D+01,
     &  0.1220370848447449D-08, 0.6262891156513252D+01,
     &  0.2532220032092866D-06, 0.5320536377336039D+01,
     &  0.1535114595466674D-04, 0.4465872626831029D+01,
     &  0.3785021094142701D-03, 0.3664416547450636D+01,
     &  0.4507235420342067D-02, 0.2898051276515753D+01,
     &  0.2866669103011841D-01, 0.2155502761316934D+01,
     &  0.1036036572761442D+00, 0.1428876676078373D+01,
     &  0.2209417121991433D+00, 0.7120850440423796D+00,
     &  0.2837731927515210D+00, 0.4118522463420039D-15/
      DATA ( W(I,20), T(I,20), I = 1,10) /
     &  0.1257800672437914D-12, 0.7619048541679760D+01,
     &  0.2482062362315163D-09, 0.6510590157013660D+01,
     &  0.6127490259983006D-07, 0.5578738805893195D+01,
     &  0.4402121090230841D-05, 0.4734581334046057D+01,
     &  0.1288262799619300D-03, 0.3943967350657311D+01,
     &  0.1830103131080496D-02, 0.3189014816553389D+01,
     &  0.1399783744710099D-01, 0.2458663611172367D+01,
     &  0.6150637206397690D-01, 0.1745247320814126D+01,
     &  0.1617393339840001D+00, 0.1042945348802752D+01,
     &  0.2607930634495551D+00, 0.3469641570813557D+00/
      DATA ( W(I,21), T(I,21), I = 1,11) /
     &  0.2098991219565665D-13, 0.7849382895113822D+01,
     &  0.4975368604121770D-10, 0.6751444718717456D+01,
     &  0.1450661284493093D-07, 0.5829382007304472D+01,
     &  0.1225354836148259D-05, 0.4994963944782024D+01,
     &  0.4219234742551696D-04, 0.4214343981688420D+01,
     &  0.7080477954815349D-03, 0.3469846690475375D+01,
     &  0.6439697051408779D-02, 0.2750592981052372D+01,
     &  0.3395272978654278D-01, 0.2049102468257161D+01,
     &  0.1083922856264195D+00, 0.1359765823211230D+01,
     &  0.2153337156950595D+00, 0.6780456924406435D+00,
     &  0.2702601835728773D+00, 0.0000000000000000D+00/
      DATA ( W(I,22), T(I,22), I = 1,11) /
     &  0.3479460647877136D-14, 0.8074029984021710D+01,
     &  0.9841378982346056D-11, 0.6985980424018808D+01,
     &  0.3366514159458310D-08, 0.6073074951122888D+01,
     &  0.3319853749814059D-06, 0.5247724433714421D+01,
     &  0.1334597712680954D-04, 0.4476361977310866D+01,
     &  0.2622833032559635D-03, 0.3741496350266517D+01,
     &  0.2808761047577212D-02, 0.3032404227831676D+01,
     &  0.1756907288080571D-01, 0.2341759996287707D+01,
     &  0.6719631142889003D-01, 0.1664124839117906D+01,
     &  0.1619062934136754D+00, 0.9951624222712152D+00,
     &  0.2502435965869353D+00, 0.3311793157152742D+00/
      DATA ( W(I,23), T(I,23), I = 1,12) /
     &  0.5732383167802038D-15, 0.8293386027417354D+01,
     &  0.1922935311567786D-11, 0.7214659435051866D+01,
     &  0.7670888862399855D-09, 0.6310349854448401D+01,
     &  0.8775062483861979D-07, 0.5493473986471793D+01,
     &  0.4089977244992140D-05, 0.4730724197451473D+01,
     &  0.9340818609031275D-04, 0.4004775321733304D+01,
     &  0.1167628637497855D-02, 0.3305040021752963D+01,
     &  0.8579678391465647D-02, 0.2624323634059181D+01,
     &  0.3886718370348111D-01, 0.1957327552933424D+01,
     &  0.1120733826026210D+00, 0.1299876468303978D+01,
     &  0.2099596695775429D+00, 0.6484711535344957D+00,
     &  0.2585097408088385D+00, 0.0000000000000000D+00/
      DATA ( W(I,24), T(I,24), I = 1,12) /
     &  0.9390193689041782D-16, 0.8507803519195264D+01,
     &  0.3714974152762395D-12, 0.7437890666021664D+01,
     &  0.1718664927964866D-09, 0.6541675005098631D+01,
     &  0.2267461673480609D-07, 0.5732747175251204D+01,
     &  0.1217659745442582D-05, 0.4978041374639117D+01,
     &  0.3209500565274598D-04, 0.4260383605019904D+01,
     &  0.4647187187793975D-03, 0.3569306764073560D+01,
     &  0.3976608929181313D-02, 0.2897728643223314D+01,
     &  0.2112634440896754D-01, 0.2240467851691752D+01,
     &  0.7206936401717838D-01, 0.1593480429816420D+01,
     &  0.1614595128670001D+00, 0.9534219229321088D+00,
     &  0.2408701155466405D+00, 0.3173700966294525D+00/
      DATA ( W(I,25), T(I,25), I = 1,13) /
     &  0.1530038997998690D-16, 0.8717597678399592D+01,
     &  0.7102103037003980D-13, 0.7656037955393078D+01,
     &  0.3791150000477161D-10, 0.6767464963809719D+01,
     &  0.5738023868899356D-08, 0.5966014690606704D+01,
     &  0.3530152560245470D-06, 0.5218848093644280D+01,
     &  0.1067219490520254D-04, 0.4508929922967284D+01,
     &  0.1777669069265268D-03, 0.3825900569972490D+01,
     &  0.1757850405263803D-02, 0.3162775679388193D+01,
     &  0.1085675599146230D-01, 0.2514473303952205D+01,
     &  0.4337997016764489D-01, 0.1877058369947839D+01,
     &  0.1148809243039517D+00, 0.1247311975616789D+01,
     &  0.2048510256503405D+00, 0.6224622791860757D+00,
     &  0.2481693511764858D+00, 0.0000000000000000D+00/
      DATA ( W(I,26), T(I,26), I = 1,13) /
     &  0.2480694260393664D-17, 0.8923051727828243D+01,
     &  0.1344547649663596D-13, 0.7869426697637738D+01,
     &  0.8242809443163844D-11, 0.6988088770623415D+01,
     &  0.1424293237988014D-08, 0.6193693483796317D+01,
     &  0.9986755573314568D-07, 0.5453615383857833D+01,
     &  0.3443413612308114D-05, 0.4750947483085378D+01,
     &  0.6557558694333818D-04, 0.4075427214412228D+01,
     &  0.7442025763604303D-03, 0.3420156373999979D+01,
     &  0.5302198015682246D-02, 0.2780138499509748D+01,
     &  0.2459766565712125D-01, 0.2151530090121648D+01,
     &  0.7622953220630281D-01, 0.1531215708695402D+01,
     &  0.1605865456137948D+00, 0.9165450413386282D+00,
     &  0.2324707356300776D+00, 0.3051559707592978D+00/
      DATA ( W(I,27), T(I,27), I = 1,14) /
     &  0.4003364766550257D-18, 0.9124421250672931D+01,
     &  0.2522363250873417D-14, 0.8078349274534165D+01,
     &  0.1768236511219616D-11, 0.7203876611910644D+01,
     &  0.3472604702845840D-09, 0.6416154934562174D+01,
     &  0.2761933918447925D-07, 0.5682760761629052D+01,
     &  0.1080581533683832D-05, 0.4986906410679802D+01,
     &  0.2339557671566820D-04, 0.4318417671936682D+01,
     &  0.3028398259361645D-03, 0.3670472986492407D+01,
     &  0.2471872445961970D-02, 0.3038150251871036D+01,
     &  0.1321102584046355D-01, 0.2417683983162542D+01,
     &  0.4748957556274264D-01, 0.1806045213138672D+01,
     &  0.1169962651750102D+00, 0.1200683354549981D+01,
     &  0.2000149701605136D+00, 0.5993548807899852D+00,
     &  0.2389778937255043D+00, 0.0000000000000000D+00/
      DATA ( W(I,28), T(I,28), I = 1,14) /
     &  0.6432547438801930D-19, 0.9321937814408766D+01,
     &  0.4691765569500354D-15, 0.8283069540861421D+01,
     &  0.3745901035176660D-12, 0.7415125286176065D+01,
     &  0.8326609843882241D-10, 0.6633731493950435D+01,
     &  0.7479362584613589D-08, 0.5906656325824994D+01,
     &  0.3304864449926482D-06, 0.5217223673447450D+01,
     &  0.8093584057145153D-05, 0.4555340384596974D+01,
     &  0.1188285381401779D-03, 0.3914253725963635D+01,
     &  0.1104305927857598D-02, 0.3289106970171833D+01,
     &  0.6752459709030160D-02, 0.2676201879526944D+01,
     &  0.2793578476788097D-01, 0.2072582674144621D+01,
     &  0.7977336601159966D-01, 0.1475781736957922D+01,
     &  0.1594181936613094D+00, 0.8836525629929802D+00,
     &  0.2248886297506769D+00, 0.2942517144887133D+00/
      DATA ( W(I,29), T(I,29), I = 1,15) /
     &  0.1029341808721942D-19, 0.9515812006947357D+01,
     &  0.8657491667957282D-16, 0.8483826557846555D+01,
     &  0.7842840425658472D-13, 0.7622102722480985D+01,
     &  0.1965709944734762D-10, 0.6846722135707994D+01,
     &  0.1986123546067053D-08, 0.6125635348243716D+01,
     &  0.9868968543560684D-07, 0.5442271089178505D+01,
     &  0.2721127828058099D-05, 0.4786611062352805D+01,
     &  0.4508394026980976D-04, 0.4151964855100983D+01,
     &  0.4743663738893483D-03, 0.3533533770990993D+01,
     &  0.3297972210833669D-02, 0.2927678154322886D+01,
     &  0.1559400577786720D-01, 0.2331504884065565D+01,
     &  0.5121083528871912D-01, 0.1742616232610662D+01,
     &  0.1185603192669036D+00, 0.1158946149400189D+01,
     &  0.1954459569679010D+00, 0.5786461780331649D+00,
     &  0.2307372767004873D+00, 0.0000000000000000D+00/
      DATA ( W(I,30), T(I,30), I = 1,15) /
     &  0.1640807008117853D-20, 0.9706235997359524D+01,
     &  0.1585560944966296D-16, 0.8680837722732207D+01,
     &  0.1624080129972436D-13, 0.7825051744352813D+01,
     &  0.4573425871326147D-11, 0.7055396866960296D+01,
     &  0.5178459467189710D-09, 0.6339997686869597D+01,
     &  0.2882175154047618D-07, 0.5662381850082873D+01,
     &  0.8909088868621158D-06, 0.5012600596486518D+01,
     &  0.1657998163067346D-04, 0.4384020365898051D+01,
     &  0.1965129439848249D-03, 0.3771894423159236D+01,
     &  0.1544707339866097D-02, 0.3172634639420402D+01,
     &  0.8295747557723240D-02, 0.2583402100229274D+01,
     &  0.3111177018350134D-01, 0.2001858612956431D+01,
     &  0.8278683671562172D-01, 0.1426005658374115D+01,
     &  0.1580469532090208D+00, 0.8540733517109733D+00,
     &  0.2179999718155776D+00, 0.2844387607362094D+00/
      DATA ( W(I,31), T(I,31), I = 1,16) /
     &  0.2605973854893011D-21, 0.9893385708986649D+01,
     &  0.2883352367857899D-17, 0.8874301409488794D+01,
     &  0.3328468324148409D-14, 0.8024193227361653D+01,
     &  0.1049603362311349D-11, 0.7260000488890867D+01,
     &  0.1327251483589731D-09, 0.6550014268765684D+01,
     &  0.8243931619119761D-08, 0.5877855885986261D+01,
     &  0.2845610088162858D-06, 0.5233641511712708D+01,
     &  0.5923202317686233D-05, 0.4610789797323995D+01,
     &  0.7871624069602249D-04, 0.4004600901491224D+01,
     &  0.6960312713792868D-03, 0.3411532415843158D+01,
     &  0.4221717767270697D-02, 0.2828792768157509D+01,
     &  0.1796787584344161D-01, 0.2254095000754410D+01,
     &  0.5456725889447496D-01, 0.1685497905069052D+01,
     &  0.1196831096958545D+00, 0.1121297374047009D+01,
     &  0.1911320047746435D+00, 0.5599475878410030D+00,
     &  0.2232941387424060D+00, 0.9191380905810332D-16/
      DATA ( W(I,32), T(I,32), I = 1,16) /
     &  0.4124607489018384D-22, 0.1007742267422945D+02,
     &  0.5208449591960853D-18, 0.9064399210702408D+01,
     &  0.6755290223670036D-15, 0.8219728765382246D+01,
     &  0.2378064855777808D-12, 0.7460755754121516D+01,
     &  0.3347501239801238D-10, 0.6755930830540704D+01,
     &  0.2312518412074224D-08, 0.6088964309076983D+01,
     &  0.8881290713105934D-07, 0.5450033273623426D+01,
     &  0.2059622103953437D-05, 0.4832604613244488D+01,
     &  0.3055980306089618D-04, 0.4232021109995410D+01,
     &  0.3025570258170642D-03, 0.3644781249880835D+01,
     &  0.2062051051307883D-02, 0.3068135169013122D+01,
     &  0.9903461702320572D-02, 0.2499840415187396D+01,
     &  0.3410984772609194D-01, 0.1938004905925718D+01,
     &  0.8534480827208071D-01, 0.1380980199272144D+01,
     &  0.1565389937575984D+00, 0.8272849037797656D+00,
     &  0.2117055698804795D+00, 0.2755464192302757D+00/
      DATA ( W(I,33), T(I,33), I = 1,17) /
     &  0.6506889970402893D-23, 0.1025849562613868D+02,
     &  0.9349155921250728D-19, 0.9251297851734609D+01,
     &  0.1358447599822066D-15, 0.8411842935668213D+01,
     &  0.5323023871225154D-13, 0.7657866034784416D+01,
     &  0.8316046910555386D-11, 0.6957971061087896D+01,
     &  0.6369261724538402D-09, 0.6295953125159368D+01,
     &  0.2712480030928315D-07, 0.5662046690089215D+01,
     &  0.6982937054010693D-06, 0.5049763451908822D+01,
     &  0.1152282979948381D-04, 0.4454485185983318D+01,
     &  0.1271924628565943D-03, 0.3872747224621657D+01,
     &  0.9695342064240463D-03, 0.3301836743259241D+01,
     &  0.5227605765843533D-02, 0.2739550282026145D+01,
     &  0.2030404475704196D-01, 0.2184038256077331D+01,
     &  0.5758631173578715D-01, 0.1633699795932273D+01,
     &  0.1204510521056565D+00, 0.1087107916669903D+01,
     &  0.1870581852279859D+00, 0.5429533766656418D+00,
     &  0.2165276496896071D+00, 0.0000000000000000D+00/
      DATA ( W(I,34), T(I,34), I = 1,17) /
     &  0.1023327129805438D-23, 0.1043674187100505D+02,
     &  0.1668144375578546D-19, 0.9435150833760799D+01,
     &  0.2708054709262291D-16, 0.8600705233983431D+01,
     &  0.1177937790622538D-13, 0.7851517590699072D+01,
     &  0.2036657679770825D-11, 0.7156339259134842D+01,
     &  0.1724305566745258D-09, 0.6499046353927969D+01,
     &  0.8117409040122166D-08, 0.5869927588596779D+01,
     &  0.2312034264322868D-06, 0.5262536481734633D+01,
     &  0.4227725748387717D-05, 0.4672290691400979D+01,
     &  0.5182712643366873D-04, 0.4095758971954162D+01,
     &  0.4399649667746255D-03, 0.3530261634074169D+01,
     &  0.2650824238310194D-02, 0.2973629650303989D+01,
     &  0.1155073894677711D-01, 0.2424050509756231D+01,
     &  0.3692346295804431D-01, 0.1879964366420186D+01,
     &  0.8751168135862399D-01, 0.1339990625522619D+01,
     &  0.1549420965148600D+00, 0.8028734607837125D+00,
     &  0.2059249366691133D+00, 0.2674391839515571D+00/
      DATA ( W(I,35), T(I,35), I = 1,18) /
     &  0.1604619191790137D-24, 0.1061228847764259D+02,
     &  0.2959542628709020D-20, 0.9616099851106188D+01,
     &  0.5354094198066113D-17, 0.8786471736571588D+01,
     &  0.2578597004420350D-14, 0.8041881508402966D+01,
     &  0.4921158318497817D-12, 0.7351222593955577D+01,
     &  0.4592895709255784D-10, 0.6698448669188526D+01,
     &  0.2383123333705684D-08, 0.6073899909835108D+01,
     &  0.7486460178518943D-07, 0.5471169042499284D+01,
     &  0.1511924078161164D-05, 0.4885706922446520D+01,
     &  0.2050985370117561D-04, 0.4314112806712937D+01,
     &  0.1931447146139865D-03, 0.3753736849678044D+01,
     &  0.1294831077498348D-02, 0.3202440651527510D+01,
     &  0.6300195959720374D-02, 0.2658444295466304D+01,
     &  0.2258156121393648D-01, 0.2120223604836576D+01,
     &  0.6029658086774051D-01, 0.1586437891088980D+01,
     &  0.1209324521970301D+00, 0.1055876792225099D+01,
     &  0.1832085621911518D+00, 0.5274192342262778D+00,
     &  0.2103411454127607D+00, 0.2220535009031490D-16/
      DATA ( W(I,36), T(I,36), I = 1,18) /
     &  0.2509037634634927D-25, 0.1078525331238753D+02,
     &  0.5222366200862934D-21, 0.9794276019583023D+01,
     &  0.1050294193474738D-17, 0.8969286534562617D+01,
     &  0.5587113978649920D-15, 0.8229115367471579D+01,
     &  0.1174029255105150D-12, 0.7542793039211395D+01,
     &  0.1204744586548017D-10, 0.6894347646173911D+01,
     &  0.6871083387212228D-09, 0.6274168326809511D+01,
     &  0.2373822777365743D-07, 0.5675884710106672D+01,
     &  0.5278315192800947D-06, 0.5094978513857614D+01,
     &  0.7896978086723001D-05, 0.4528076990600175D+01,
     &  0.8220025562410442D-04, 0.3972557341929988D+01,
     &  0.6107548335511020D-03, 0.3426308595129129D+01,
     &  0.3304134538435289D-02, 0.2887579695004719D+01,
     &  0.1321657401560215D-01, 0.2354877715992540D+01,
     &  0.3955236976655980D-01, 0.1826896577986744D+01,
     &  0.8934249750438648D-01, 0.1302464954480165D+01,
     &  0.1532910133997116D+00, 0.7805064920524665D+00,
     &  0.2005920064390222D+00, 0.2600079252490002D+00/
      DATA ( W(I,37), T(I,37), I = 1,19) /
     &  0.3912701900272739D-26, 0.1095574594356165D+02,
     &  0.9167997950194993D-22, 0.9969800945691102D+01,
     &  0.2045033189096993D-18, 0.9149282977696753D+01,
     &  0.1198838092763842D-15, 0.8413364679488879D+01,
     &  0.2767202859289950D-13, 0.7731209035776908D+01,
     &  0.3114550633754001D-11, 0.7086915685017345D+01,
     &  0.1947529536776704D-09, 0.6470920475390299D+01,
     &  0.7379429831202629D-08, 0.5876887892594201D+01,
     &  0.1801392729842884D-06, 0.5300328473991422D+01,
     &  0.2963204697508371D-05, 0.4737895299703953D+01,
     &  0.3397941877127311D-04, 0.4186990225627911D+01,
     &  0.2788064134714879D-03, 0.3645526993515297D+01,
     &  0.1670452623119233D-02, 0.3111780274991816D+01,
     &  0.7424836460756516D-02, 0.2584284718210775D+01,
     &  0.2478561700046269D-01, 0.2061764482625703D+01,
     &  0.6272612686101582D-01, 0.1543082026656558D+01,
     &  0.1211816464812620D+00, 0.1027199336691835D+01,
     &  0.1795672590244481D+00, 0.5131472568284307D+00,
     &  0.2046562495907946D+00, 0.0000000000000000D+00/
      DATA ( W(I,38), T(I,38), I = 1,19) /
     &  0.6086019568424894D-27, 0.1112386843494987D+02,
     &  0.1601586834974089D-22, 0.1014278766112967D+02,
     &  0.3953752210235847D-19, 0.9326584757395571D+01,
     &  0.2548653789376282D-16, 0.8594764136386107D+01,
     &  0.6447906524854284D-14, 0.7916616928480361D+01,
     &  0.7941720700519787D-12, 0.7276311665555140D+01,
     &  0.5431469647397207D-10, 0.6664328864414711D+01,
     &  0.2251484789660455D-08, 0.6074366042059564D+01,
     &  0.6017570329946914D-07, 0.5501960756811018D+01,
     &  0.1085199304765951D-05, 0.4943790029802265D+01,
     &  0.1366659848918500D-04, 0.4397278309932714D+01,
     &  0.1234207566124386D-03, 0.3860361738163258D+01,
     &  0.8159914125966156D-03, 0.3331338060544006D+01,
     &  0.4014342501780207D-02, 0.2808766413901917D+01,
     &  0.1488362222069118D-01, 0.2291397563624214D+01,
     &  0.4200051997416258D-01, 0.1778123425136563D+01,
     &  0.9088415555258827D-01, 0.1267939114780749D+01,
     &  0.1516111462601078D+00, 0.7599132481367392D+00,
     &  0.1956519870413640D+00, 0.2531636353359348D+00/
      DATA ( W(I,39), T(I,39), I = 1,20) /
     &  0.9443344575063092D-28, 0.1128971604447632D+02,
     &  0.2784787505225604D-23, 0.1031334144275699D+02,
     &  0.7592450542206611D-20, 0.9501306853782864D+01,
     &  0.5370701458462819D-17, 0.8773438698067656D+01,
     &  0.1486129587733307D-14, 0.8099152213152410D+01,
     &  0.1998726680623669D-12, 0.7462682377866057D+01,
     &  0.1491720010448932D-10, 0.6854552519790964D+01,
     &  0.6748646364787737D-09, 0.6268491549685291D+01,
     &  0.1969872920599359D-07, 0.5700062454271859D+01,
     &  0.3884118283713278D-06, 0.5145964544094186D+01,
     &  0.5356584901373566D-05, 0.4603643074288992D+01,
     &  0.5307742112417270D-04, 0.4071054595947238D+01,
     &  0.3859381697690359D-03, 0.3546517669768644D+01,
     &  0.2093837438880663D-02, 0.3028613314092798D+01,
     &  0.8588083029362219D-02, 0.2516115854339318D+01,
     &  0.2690624148395815D-01, 0.2007943065265498D+01,
     &  0.6490157458316884D-01, 0.1503118900228415D+01,
     &  0.1212421280631244D+00, 0.1000744572356018D+01,
     &  0.1761190277014499D+00, 0.4999751896072681D+00,
     &  0.1994086534474408D+00, 0.1536204102353874D-15/
      DATA ( W(I,40), T(I,40), I = 1,20) /
     &  0.1461839873869467D-28, 0.1145337784154873D+02,
     &  0.4820467940200524D-24, 0.1048156053467427D+02,
     &  0.1448609431551587D-20, 0.9673556366934033D+01,
     &  0.1122275206827074D-17, 0.8949504543855559D+01,
     &  0.3389853443248306D-15, 0.8278940623659475D+01,
     &  0.4968088529197761D-13, 0.7646163764541459D+01,
     &  0.4037638581695192D-11, 0.7041738406453829D+01,
     &  0.1989118526027766D-09, 0.6459423377583766D+01,
     &  0.6325897188548972D-08, 0.5894805675372016D+01,
     &  0.1360342421574886D-06, 0.5344605445720084D+01,
     &  0.2048897436081474D-05, 0.4806287192093873D+01,
     &  0.2221177143247582D-04, 0.4277826156362752D+01,
     &  0.1770729287992397D-03, 0.3757559776168985D+01,
     &  0.1055879016901825D-02, 0.3244088732999869D+01,
     &  0.4773544881823334D-02, 0.2736208340465433D+01,
     &  0.1653784414256937D-01, 0.2232859218634873D+01,
     &  0.4427455520227679D-01, 0.1733090590631720D+01,
     &  0.9217657917006089D-01, 0.1236032004799159D+01,
     &  0.1499211117635710D+00, 0.7408707252859313D+00,
     &  0.1910590096619904D+00, 0.2468328960227240D+00/
      DATA ( W(I,41), T(I,41), I = 1,21) /
     &  0.2257863956583089D-29, 0.1161493725433746D+02,
     &  0.8308558938782992D-25, 0.1064753678631932D+02,
     &  0.2746891228522292D-21, 0.9843433249157988D+01,
     &  0.2326384145587187D-18, 0.9123069907984480D+01,
     &  0.7655982291966812D-16, 0.8456099083269388D+01,
     &  0.1220334874202772D-13, 0.7826882004053867D+01,
     &  0.1077818394935909D-11, 0.7226022663732790D+01,
     &  0.5769853428092003D-10, 0.6647308470747191D+01,
     &  0.1994794756757339D-08, 0.6086349164878472D+01,
     &  0.4667347708107243D-07, 0.5539884440458126D+01,
     &  0.7658186077982435D-06, 0.5005396683404125D+01,
     &  0.9058608622433030D-05, 0.4480878331594004D+01,
     &  0.7894719319504627D-04, 0.3964684028033266D+01,
     &  0.5158014443431912D-03, 0.3455432217780992D+01,
     &  0.2561642428649777D-02, 0.2951937016381193D+01,
     &  0.9777902738208298D-02, 0.2453159345907049D+01,
     &  0.2893721174793441D-01, 0.1958170711977291D+01,
     &  0.6684765935446599D-01, 0.1466125457295967D+01,
     &  0.1211489170115104D+00, 0.9762387671800500D+00,
     &  0.1728495310506020D+00, 0.4877685693194347D+00,
     &  0.1945450277536008D+00, 0.2585532684499631D-15/
      DATA ( W(I,42), T(I,42), I = 1,21) /
     &  0.3479841758734498D-30, 0.1177447255645880D+02,
     &  0.1426197845863333D-25, 0.1081135621818894D+02,
     &  0.5178070329449428D-22, 0.1001103095231325D+02,
     &  0.4785541849652557D-19, 0.9294235815925036D+01,
     &  0.1712817711028008D-16, 0.8630736540442662D+01,
     &  0.2963878417982294D-14, 0.8004954459331870D+01,
     &  0.2839400066530831D-12, 0.7407531683161432D+01,
     &  0.1648408975918176D-10, 0.6832282984221202D+01,
     &  0.6182418956905496D-09, 0.6274839704458262D+01,
     &  0.1570405641380752D-07, 0.5731959941924930D+01,
     &  0.2800444926136691D-06, 0.5201142761234950D+01,
     &  0.3605310211304410D-05, 0.4680396489703305D+01,
     &  0.3425737772910867D-04, 0.4168091525525806D+01,
     &  0.2445266890869028D-03, 0.3662862441045986D+01,
     &  0.1329885905572134D-02, 0.3163540283549010D+01,
     &  0.5573924561218487D-02, 0.2669104116603810D+01,
     &  0.1816809011555160D-01, 0.2178645208762502D+01,
     &  0.4638274179115778D-01, 0.1691339732834912D+01,
     &  0.9325376860459540D-01, 0.1206427277827926D+01,
     &  0.1482345543992958D+00, 0.7231933449391704D+00,
     &  0.1867743488620199D+00, 0.2409545348665914D+00/
      DATA ( W(I,43), T(I,43), I = 1,22) /
     &  0.5352075224079761D-31, 0.1193205730105631D+02,
     &  0.2438509785586889D-26, 0.1097309952495776D+02,
     &  0.9705961174387400D-23, 0.1017643700187913D+02,
     &  0.9772188868027905D-20, 0.9463096735522324D+01,
     &  0.3797453851791020D-17, 0.8802954705722227D+01,
     &  0.7121219095950012D-15, 0.8180490511448660D+01,
     &  0.7386428797509359D-13, 0.7586383052570858D+01,
     &  0.4641647576106650D-11, 0.7014473354182730D+01,
     &  0.1884798249393598D-09, 0.6460413330895460D+01,
     &  0.5186705460818634D-08, 0.5920978461471519D+01,
     &  0.1003006675541783D-06, 0.5393683423125381D+01,
     &  0.1402088827949991D-05, 0.4876551284706595D+01,
     &  0.1448871182961238D-04, 0.4367966934806672D+01,
     &  0.1126824016383807D-03, 0.3866579655692208D+01,
     &  0.6691683916969021D-03, 0.3371235817148052D+01,
     &  0.3070012395787847D-02, 0.2880930777450376D+01,
     &  0.1098373316809114D-01, 0.2394773428842528D+01,
     &  0.3087516724563331D-01, 0.1911959274951712D+01,
     &  0.6858704943109335D-01, 0.1431749364451692D+01,
     &  0.1209303683096624D+00, 0.9534532756297385D+00,
     &  0.1697454597838771D+00, 0.4764148861084949D+00,
     &  0.1900207247825863D+00, 0.0000000000000000D+00/
      DATA ( W(I,44), T(I,44), I = 1,22) /
     &  0.8215250893174055D-32, 0.1208776070905845D+02,
     &  0.4153634860809548D-27, 0.1113284252424393D+02,
     &  0.1809481070701256D-23, 0.1033973350764328D+02,
     &  0.1981515599992759D-20, 0.9629741154666469D+01,
     &  0.8346676981526894D-18, 0.8972848703628870D+01,
     &  0.1693427110357361D-15, 0.8353592294951005D+01,
     &  0.1898518041004563D-13, 0.7762686386163553D+01,
     &  0.1289064815963923D-11, 0.7193997236488432D+01,
     &  0.5656565204102218D-10, 0.6643196399727595D+01,
     &  0.1683045805404500D-08, 0.6107075817029749D+01,
     &  0.3522080139402352D-07, 0.5583164829645622D+01,
     &  0.5334152162808288D-06, 0.5069500234589782D+01,
     &  0.5980506497871435D-05, 0.4564480302199131D+01,
     &  0.5055024343502252D-04, 0.4066767790495557D+01,
     &  0.3269029726329101D-03, 0.3575222999425195D+01,
     &  0.1636879726834301D-02, 0.3088855994142003D+01,
     &  0.6407990218344719D-02, 0.2606791462987426D+01,
     &  0.1976568083382755D-01, 0.2128242119185467D+01,
     &  0.4833422729538417D-01, 0.1652487989479395D+01,
     &  0.9414471679496565D-01, 0.1178859803328855D+01,
     &  0.1465614441714095D+00, 0.7067252381063451D+00,
     &  0.1827650568597313D+00, 0.2354771181719222D+00/
      DATA ( W(I,45), T(I,45), I = 1,23) /
     &  0.1258601266761721D-32, 0.1224164801738294D+02,
     &  0.7049452438028993D-28, 0.1129065655801674D+02,
     &  0.3355898215846036D-24, 0.1050099761933589D+02,
     &  0.3990921764567262D-21, 0.9794252095357150D+01,
     &  0.1819415421781025D-18, 0.9140507651218687D+01,
     &  0.3987388394988911D-16, 0.8524355348620958D+01,
     &  0.4823870286948576D-14, 0.7936544056928108D+01,
     &  0.3532976789571211D-12, 0.7370964332170609D+01,
     &  0.1672373578472200D-10, 0.6823306517591415D+01,
     &  0.5370107061670978D-09, 0.6290378188903940D+01,
     &  0.1213731971819947D-07, 0.5769722503692224D+01,
     &  0.1987368663904533D-06, 0.5259389088450044D+01,
     &  0.2412167676023244D-05, 0.4757788618871960D+01,
     &  0.2210680714129846D-04, 0.4263596248773437D+01,
     &  0.1552889103199588D-03, 0.3775684997351202D+01,
     &  0.8463565255385183D-03, 0.3293078264001083D+01,
     &  0.3614802525284821D-02, 0.2814914962782011D+01,
     &  0.1219644981036976D-01, 0.2340423205074869D+01,
     &  0.3271890357088906D-01, 0.1868899890086888D+01,
     &  0.7014033220973423D-01, 0.1399694432579061D+01,
     &  0.1206095596750258D+00, 0.9321954002025563D+00,
     &  0.1667945553649214D+00, 0.4658191781783683D+00,
     &  0.1857980420096402D+00, 0.0000000000000000D+00/
      DATA ( W(I,46), T(I,46), I = 1,23) /
     &  0.1924664627046275D-33, 0.1239378079201923D+02,
     &  0.1192246272132870D-28, 0.1144660885260523D+02,
     &  0.6192859698073460D-25, 0.1066030193428052D+02,
     &  0.7986127260533585D-22, 0.9956707572498630D+01,
     &  0.3934545090204791D-19, 0.9306015173119082D+01,
     &  0.9300353507700858D-17, 0.8692869193236842D+01,
     &  0.1212244415430154D-14, 0.8108051845054231D+01,
     &  0.9561363845435019D-13, 0.7545477116064648D+01,
     &  0.4874112963243041D-11, 0.7000853362533224D+01,
     &  0.1686106290702903D-09, 0.6471003045526279D+01,
     &  0.4108222480852658D-08, 0.5953482378238384D+01,
     &  0.7258481863429432D-07, 0.5446353016216582D+01,
     &  0.9517578349775918D-06, 0.4948037176754664D+01,
     &  0.9436598821520223D-05, 0.4457221460031489D+01,
     &  0.7183252612530643D-04, 0.3972790546235975D+01,
     &  0.4250585899081316D-03, 0.3493779976317873D+01,
     &  0.1975260370663049D-02, 0.3019341533379349D+01,
     &  0.7268729480993403D-02, 0.2548717171165183D+01,
     &  0.2132401414735712D-01, 0.2081218863808104D+01,
     &  0.5013852830109192D-01, 0.1616212619554741D+01,
     &  0.9487419222577377D-01, 0.1153105444663487D+01,
     &  0.1449090160374568D+00, 0.6913343908693306D+00,
     &  0.1790029030973513D+00, 0.2303570440689780D+00/
      DATA ( W(I,47), T(I,47), I = 1,24) /
     &  0.2937991918931770D-34, 0.1254421721021957D+02,
     &  0.2009633017503842D-29, 0.1160076284239976D+02,
     &  0.1137324853690828D-25, 0.1081771486308368D+02,
     &  0.1588168143815981D-22, 0.1011718100450289D+02,
     &  0.8443830093861547D-20, 0.9469449861345275D+01,
     &  0.2149644881290608D-17, 0.8859217846077781D+01,
     &  0.3014376026444397D-15, 0.8277299513807538D+01,
     &  0.2556494006562142D-13, 0.7717631482085797D+01,
     &  0.1401228890589135D-11, 0.7175939408289193D+01,
     &  0.5213279817493255D-10, 0.6649059958212412D+01,
     &  0.1366926136038253D-08, 0.6134561715294918D+01,
     &  0.2601160240819878D-07, 0.5630517648304252D+01,
     &  0.3677447240564392D-06, 0.5135360748767347D+01,
     &  0.3936411635330397D-05, 0.4647788224870149D+01,
     &  0.3239911617951789D-04, 0.4166695488690111D+01,
     &  0.2076565749155176D-03, 0.3691129181403769D+01,
     &  0.1047278482271871D-02, 0.3220252777152719D+01,
     &  0.4191754574650512D-02, 0.2753320728634848D+01,
     &  0.1340827900837288D-01, 0.2289658541588672D+01,
     &  0.3446881631958411D-01, 0.1828647032406020D+01,
     &  0.7152609366400937D-01, 0.1369709566513547D+01,
     &  0.1202053639683183D+00, 0.9123014234783601D+00,
     &  0.1639855806134269D+00, 0.4559006613049348D+00,
     &  0.1818448921796476D+00, 0.0000000000000000D+00/
      DATA ( W(I,48), T(I,48), I = 1,24) /
     &  0.4477155473876067D-35, 0.1269301231543958D+02,
     &  0.3376456143135106D-30, 0.1175317846161662D+02,
     &  0.2079058971418782D-26, 0.1097330095851353D+02,
     &  0.3139476644799137D-23, 0.1027574158173351D+02,
     &  0.1798854916237414D-20, 0.9630885686949055D+01,
     &  0.4925465109698832D-18, 0.9023480280417658D+01,
     &  0.7419993598065164D-16, 0.8444371322523445D+01,
     &  0.6756677274657091D-14, 0.7887517316500372D+01,
     &  0.3975805958471478D-12, 0.7348660565908911D+01,
     &  0.1588360981246298D-10, 0.6824651320745980D+01,
     &  0.4474287153436318D-09, 0.6313069914905501D+01,
     &  0.9154055746313672D-08, 0.5811999987702533D+01,
     &  0.1392791689559962D-06, 0.5319884620419087D+01,
     &  0.1606393695542611D-05, 0.4835430885878957D+01,
     &  0.1426609232424675D-04, 0.4357544108790293D+01,
     &  0.9881804917611217D-04, 0.3885281117248394D+01,
     &  0.5395865846271863D-03, 0.3417816048684776D+01,
     &  0.2343082111760394D-02, 0.2954414688880120D+01,
     &  0.8149696855351926D-02, 0.2494414743040025D+01,
     &  0.2283821210009132D-01, 0.2037210303259033D+01,
     &  0.5180518354956036D-01, 0.1582239319674964D+01,
     &  0.9546340056142952D-01, 0.1128973231509649D+01,
     &  0.1432824569931762D+00, 0.6769081422063483D+00,
     &  0.1754635418118663D+00, 0.2255570729801638D+00/
      DATA ( W(I,49), T(I,49), I = 1,25) /
     &  0.6811389123116583D-36, 0.1284021824817448D+02,
     &  0.5655218285268972D-31, 0.1190391240788757D+02,
     &  0.3783656387961836D-27, 0.1112712121198839D+02,
     &  0.6170433961110509D-24, 0.1043245459795100D+02,
     &  0.3805275010242669D-21, 0.9790392369513130D+01,
     &  0.1119149962220974D-18, 0.9185730837014976D+01,
     &  0.1808768646225968D-16, 0.8609346484896321D+01,
     &  0.1765991900170034D-14, 0.8055219008757504D+01,
     &  0.1113983833798787D-12, 0.7519106753949869D+01,
     &  0.4771647040124784D-11, 0.6997872987021261D+01,
     &  0.1441769692730330D-09, 0.6489109229832087D+01,
     &  0.3166143293671931D-08, 0.5990909213118121D+01,
     &  0.5175289227690812D-07, 0.5501725495034997D+01,
     &  0.6419581525894787D-06, 0.5020274351287227D+01,
     &  0.6139401305210235D-05, 0.4545470293840825D+01,
     &  0.4586356164456607D-04, 0.4076379533242183D+01,
     &  0.2705406889774859D-03, 0.3612185969516418D+01,
     &  0.1271499268760875D-02, 0.3152165704938866D+01,
     &  0.4796631268346652D-02, 0.2695667488564915D+01,
     &  0.1461268197482865D-01, 0.2242097365824569D+01,
     &  0.3612646737149015D-01, 0.1790906348772901D+01,
     &  0.7276104720165333D-01, 0.1341580271678923D+01,
     &  0.1197332834080540D+00, 0.8936312256129486D+00,
     &  0.1613082628631473D+00, 0.4465901175404579D+00,
     &  0.1781337719310835D+00, 0.0000000000000000D+00/
      DATA ( W(I,50), T(I,50), I = 1,25) /
     &  0.1034607500576990D-36, 0.1298588445541555D+02,
     &  0.9443414659584510D-32, 0.1205301838092448D+02,
     &  0.6856280758924735D-28, 0.1127923332148262D+02,
     &  0.1206044550761014D-24, 0.1058738174919177D+02,
     &  0.7995094477915292D-22, 0.9948035709637500D+01,
     &  0.2522482807168144D-19, 0.9346039593575728D+01,
     &  0.4368171816201588D-17, 0.8772299579514598D+01,
     &  0.4566698246800344D-15, 0.8220815907982127D+01,
     &  0.3083828687005300D-13, 0.7687362406712500D+01,
     &  0.1414228936126661D-11, 0.7168814837853899D+01,
     &  0.4576636712310442D-10, 0.6662775399018720D+01,
     &  0.1077060789389039D-08, 0.6167347388659921D+01,
     &  0.1888225976835208D-07, 0.5680992291033284D+01,
     &  0.2514609880838772D-06, 0.5202434993399912D+01,
     &  0.2584937658949391D-05, 0.4730598550228594D+01,
     &  0.2078485175734569D-04, 0.4264557843038109D+01,
     &  0.1321726328668984D-03, 0.3803505741742012D+01,
     &  0.6708280619787080D-03, 0.3346727774732429D+01,
     &  0.2738160896935348D-02, 0.2893582727707738D+01,
     &  0.9045054154849623D-02, 0.2443487452654017D+01,
     &  0.2430481286424306D-01, 0.1995904709795124D+01,
     &  0.5334352453170102D-01, 0.1550333214338771D+01,
     &  0.9593054035810168D-01, 0.1106299289397183D+01,
     &  0.1416854132499443D+00, 0.6633496795082918D+00,
     &  0.1721258519924433D+00, 0.2210451816445435D+00/
*
*
****END HERMIT
*
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c==============programme Jeremie

       FUNCTION hermite(npts,funsub,ier) RESULT(fn_val)

!************************************************************************

      IMPLICIT NONE
       INTEGER, INTENT(IN)  :: npts
       INTEGER, INTENT(IN)  :: ier
      double precision            :: fn_val,funvls1,funvls2
      double precision,dimension(2)::X
      double precision,dimension(25,50)::Wg,Tg
      double precision,dimension(2,50)::gauss
      integer::i,ndim,j,nf
      external::funsub

      DATA ( Wg(I, 5), Tg(I, 5), I = 1, 3) /
     &  0.1125741132772071D-01, 0.2856970013872805D+01,
     &  0.2220759220056126D+00, 0.1355626179974265D+01,
     &  0.5333333333333342D+00, 0.9386691848789097D-16/
      DATA ( Wg(I, 7), Tg(I, 7), I = 1, 4) /
     &  0.5482688559722184D-03, 0.3750439717725742D+01,
     &  0.3075712396758645D-01, 0.2366759410734542D+01,
     &  0.2401231786050126D+00, 0.1154405394739968D+01,
     &  0.4571428571428575D+00, 0.2669848554723344D-16/
      DATA ( Wg(I, 9), Tg(I, 9), I = 1, 5) /
     &  0.2234584400774664D-04, 0.4512745863399781D+01,
     &  0.2789141321231769D-02, 0.3205429002856470D+01,
     &  0.4991640676521780D-01, 0.2076847978677829D+01,
     &  0.2440975028949394D+00, 0.1023255663789133D+01,
     &  0.4063492063492066D+00, 0.0000000000000000D+00/
      DATA ( Wg(I,15), Tg(I,15), I = 1, 8) /
     &  0.8589649899633300D-09, 0.6363947888829836D+01,
     &  0.5975419597920602D-06, 0.5190093591304780D+01,
     &  0.5642146405189029D-04, 0.4196207711269018D+01,
     &  0.1567357503549958D-02, 0.3289082424398766D+01,
     &  0.1736577449213763D-01, 0.2432436827009758D+01,
     &  0.8941779539984458D-01, 0.1606710069028730D+01,
     &  0.2324622936097322D+00, 0.7991290683245483D+00,
     &  0.3182595182595181D+00, 0.0000000000000000D+00/
      DATA ( Wg(I,20), Tg(I,20), I = 1,10) /
     &  0.1257800672437914D-12, 0.7619048541679760D+01,
     &  0.2482062362315163D-09, 0.6510590157013660D+01,
     &  0.6127490259983006D-07, 0.5578738805893195D+01,
     &  0.4402121090230841D-05, 0.4734581334046057D+01,
     &  0.1288262799619300D-03, 0.3943967350657311D+01,
     &  0.1830103131080496D-02, 0.3189014816553389D+01,
     &  0.1399783744710099D-01, 0.2458663611172367D+01,
     &  0.6150637206397690D-01, 0.1745247320814126D+01,
     &  0.1617393339840001D+00, 0.1042945348802752D+01,
     &  0.2607930634495551D+00, 0.3469641570813557D+00/
       DATA ( Wg(I,30), Tg(I,30), I = 1,15) /
     &  0.1640807008117853D-20, 0.9706235997359524D+01,
     &  0.1585560944966296D-16, 0.8680837722732207D+01,
     &  0.1624080129972436D-13, 0.7825051744352813D+01,
     &  0.4573425871326147D-11, 0.7055396866960296D+01,
     &  0.5178459467189710D-09, 0.6339997686869597D+01,
     &  0.2882175154047618D-07, 0.5662381850082873D+01,
     &  0.8909088868621158D-06, 0.5012600596486518D+01,
     &  0.1657998163067346D-04, 0.4384020365898051D+01,
     &  0.1965129439848249D-03, 0.3771894423159236D+01,
     &  0.1544707339866097D-02, 0.3172634639420402D+01,
     &  0.8295747557723240D-02, 0.2583402100229274D+01,
     &  0.3111177018350134D-01, 0.2001858612956431D+01,
     &  0.8278683671562172D-01, 0.1426005658374115D+01,
     &  0.1580469532090208D+00, 0.8540733517109733D+00,
     &  0.2179999718155776D+00, 0.2844387607362094D+00/
      DATA ( Wg(I,40), Tg(I,40), I = 1,20) /
     &  0.1461839873869467D-28, 0.1145337784154873D+02,
     &  0.4820467940200524D-24, 0.1048156053467427D+02,
     &  0.1448609431551587D-20, 0.9673556366934033D+01,
     &  0.1122275206827074D-17, 0.8949504543855559D+01,
     &  0.3389853443248306D-15, 0.8278940623659475D+01,
     &  0.4968088529197761D-13, 0.7646163764541459D+01,
     &  0.4037638581695192D-11, 0.7041738406453829D+01,
     &  0.1989118526027766D-09, 0.6459423377583766D+01,
     &  0.6325897188548972D-08, 0.5894805675372016D+01,
     &  0.1360342421574886D-06, 0.5344605445720084D+01,
     &  0.2048897436081474D-05, 0.4806287192093873D+01,
     &  0.2221177143247582D-04, 0.4277826156362752D+01,
     &  0.1770729287992397D-03, 0.3757559776168985D+01,
     &  0.1055879016901825D-02, 0.3244088732999869D+01,
     &  0.4773544881823334D-02, 0.2736208340465433D+01,
     &  0.1653784414256937D-01, 0.2232859218634873D+01,
     &  0.4427455520227679D-01, 0.1733090590631720D+01,
     &  0.9217657917006089D-01, 0.1236032004799159D+01,
     &  0.1499211117635710D+00, 0.7408707252859313D+00,
     &  0.1910590096619904D+00, 0.2468328960227240D+00/
      DATA ( Wg(I,50), Tg(I,50), I = 1,25) /
     &  0.1034607500576990D-36, 0.1298588445541555D+02,
     &  0.9443414659584510D-32, 0.1205301838092448D+02,
     &  0.6856280758924735D-28, 0.1127923332148262D+02,
     &  0.1206044550761014D-24, 0.1058738174919177D+02,
     &  0.7995094477915292D-22, 0.9948035709637500D+01,
     &  0.2522482807168144D-19, 0.9346039593575728D+01,
     &  0.4368171816201588D-17, 0.8772299579514598D+01,
     &  0.4566698246800344D-15, 0.8220815907982127D+01,
     &  0.3083828687005300D-13, 0.7687362406712500D+01,
     &  0.1414228936126661D-11, 0.7168814837853899D+01,
     &  0.4576636712310442D-10, 0.6662775399018720D+01,
     &  0.1077060789389039D-08, 0.6167347388659921D+01,
     &  0.1888225976835208D-07, 0.5680992291033284D+01,
     &  0.2514609880838772D-06, 0.5202434993399912D+01,
     &  0.2584937658949391D-05, 0.4730598550228594D+01,
     &  0.2078485175734569D-04, 0.4264557843038109D+01,
     &  0.1321726328668984D-03, 0.3803505741742012D+01,
     &  0.6708280619787080D-03, 0.3346727774732429D+01,
     &  0.2738160896935348D-02, 0.2893582727707738D+01,
     &  0.9045054154849623D-02, 0.2443487452654017D+01,
     &  0.2430481286424306D-01, 0.1995904709795124D+01,
     &  0.5334352453170102D-01, 0.1550333214338771D+01,
     &  0.9593054035810168D-01, 0.1106299289397183D+01,
     &  0.1416854132499443D+00, 0.6633496795082918D+00,
     &  0.1721258519924433D+00, 0.2210451816445435D+00/
c      INTERFACE
c         subroutine funsub(NDIM,X,NF,FUNVLS)
c         IMPLICIT NONE 
c         integer :: ndim,nf
c	double precision, dimension(2)::X
c    	double precision :: funvls


c        end subroutine funsub
c      END INTERFACE


      if(npts.ne.5.and.npts.ne.7.and.npts.ne.9.and.npts.ne.15.and.
     & npts.ne.20.and.npts.ne.30.and.npts.ne.40.and.npts.ne.50) then
c       write(*,*)'nb pts GH = 5,7,9,15,20,30,40, ou 50'
       stop
       end if

ccccccccccccccccccccccccccccccccccccccccccccccc
c        print*,"npts",npts
    
       fn_val=0.d0

c       open(2,file="integrant.dat")

        DO J = 1, NPTS- MOD( NPTS, 2 )
          
                      if (j .le. npts/2) then
            GAUSS(1,j) = -Tg(j,NPTS)
            GAUSS(2,j) =  Wg(j,NPTS)
            X(2)=gauss(1,j)
c             print*,"x(2) ds hermite",X(2)          
            end if
         
           if (j .ge. npts/2 +1 ) then
            GAUSS(1,j) = Tg(NPTS/2-NPTS+j+MOD( NPTS, 2 ),NPTS)
            GAUSS(2,j) = Wg(NPTS/2-NPTS+j+MOD( NPTS, 2 ),NPTS)   
           X(2)=gauss(1,j)
c            print*,"x(2) ds hermite",X(2)  
               end if
              
        DO I = 1, NPTS- MOD( NPTS, 2 ) 
c           print*,"i",i,NPTS- MOD( NPTS, 2 )-1 
           if (i .le. npts/2) then
            GAUSS(1,I) = -Tg(I,NPTS)
            GAUSS(2,I) =  Wg(I,NPTS)
        
            X(1)=gauss(1,i)
            end if
            if (i .ge. npts/2 +1 ) then
            GAUSS(1,i) = Tg(NPTS/2-NPTS+i+MOD( NPTS, 2 ) ,NPTS)
            GAUSS(2,i) = Wg(NPTS/2-NPTS+i+MOD( NPTS, 2 ),NPTS)   
           X(1)=gauss(1,i)
               end if
c              print*,"X(1)",X(1),GAUSS(2,I),"X(2)",X(2),GAUSS(2,j)
               call funsub(NDIM,X,NF,FUNVLS1)
            
           fn_val=fn_val+ gauss(2,j)*GAUSS(2,I)*funvls1 
c           write(2,'(f15.7,1x,f23.7,f15.7,1x,f23.7,*)'), 

c ===========pour la représentation graphique de l'intégrant :
c           write(2,*),
c     &          X(1),X(2),funvls1*exp(-(5.d-1)*dot_product(X,X))
         
c            print*,"gauss(2,i)",gauss(2,i)
         END DO
       
         IF ( MOD( NPTS, 2 ) .EQ. 1 ) THEN
            GAUSS(1, NPTS/2 + 1 ) = 0.D0
            GAUSS(2, NPTS/2 + 1 ) = Wg( NPTS/2 + 1, NPTS ) 
            X(1)=0
            call funsub(NDIM,X,NF,FUNVLS1)      
             fn_val=fn_val+ gauss(2,j)*GAUSS(2,NPTS/2+1)*funvls1
cc            write(2,'(f15.7,1x,f23.7,f15.7,1x,f23.7)'),
cc     &X(1),X(2),funvls1*exp(-(5.d-1)*dot_product(X,X))
         END IF
      end do

       IF ( MOD( NPTS, 2 ) .EQ. 1 ) THEN
          gauss(2,j)=GAUSS(2,NPTS/2+1)
            X(2)=0

            DO I = 1, NPTS/2
            GAUSS(1,I) = -Tg(I,NPTS)
            GAUSS(2,I) =  Wg(I,NPTS)
            GAUSS(1,NPTS-I+1) = Tg(I,NPTS)
            GAUSS(2,NPTS-I+1) = Wg(I,NPTS)
            X(1)=gauss(1,i)

            call funsub(NDIM,X,NF,FUNVLS1)
cc        write(2,'(f15.7,1x,f23.7,f15.7,1x,f23.7)'),
cc     &X(1),X(2),funvls1*exp(-(5.d-1)*dot_product(X,X))    
            X(1)=-gauss(1,i)
            call funsub(NDIM,X,NF,FUNVLS2)
cc         write(2,'(f15.7,1x,f23.7,f15.7,1x,f23.7)'),
cc     &X(1),X(2),funvls1*exp(-(5.d-1)*dot_product(X,X))
 
              fn_val=fn_val+ gauss(2,i)*GAUSS(2,NPTS/2+1)*(funvls1+
     &funvls2)
              end do

              X=0.d0
              call funsub(NDIM,X,NF,FUNVLS1)
            fn_val=fn_val+ 4*(GAUSS(2,NPTS/2+1)**2)*funvls1  
cc         write(2,'(f15.7,1x,f23.7,f15.7,1x,f23.7)'),
cc     &X(1),X(2),funvls1*exp(-(5.d-1)*dot_product(X,X))    
           end if
            fn_val=fn_val
cc            close(2)
                RETURN

                END FUNCTION hermite

c=====================cross validation

c========================          MNBRAK         ===================
      subroutine mnbrak3(ax,bx,cx,fa,fb,fc,b,n)
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************

         double precision ax,bx,cx,fa,fb,fc,aux,res
         double precision b(npmax),y(npmax,npmax)
         double precision estimv3,gold,glimit,tiny
         parameter (gold=1.618034d0,glimit=100.d0,tiny=1.d-20)
         double precision dum,fu,q,r,u,ulim
         integer n,ni

c     write(*,*)' DEBUT mnbrak',fa,fb
         fa = estimv3(ax,n,b,y,aux,ni,res)
ccc   write(*,*)'fa mnbrak',fa
c     stop
         fb = estimv3(bx,n,b,y,aux,ni,res)
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
         fc = estimv3(cx,n,b,y,aux,ni,res)
c         write(*,*)'ok mnbrak'
c         stop
 1       if(fb.ge.fc)then
            r = (bx-ax)*(fb-fc)
            q = (bx-cx)*(fb-fa)
            u = bx-((bx-cx)*q-(bx-ax)*r)/
     &       (2.d0*sign(max(abs(q-r),tiny),q-r))
            ulim = bx + glimit*(cx-bx)
            if((bx-u)*(u-cx).gt.0.d0)then
               fu = estimv3(u,n,b,y,aux,ni,res)
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
               fu = estimv3(u,n,b,y,aux,ni,res)
            else
               if((cx-u)*(u-ulim).gt.0.d0)then
                  fu = estimv3(u,n,b,y,aux,ni,res)
                  if(fu.lt.fc)then
                     bx = cx
                     cx = u
                     u = cx + gold*(cx-bx)
                     fb = fc
                     fc = fu
                     fu = estimv3(u,n,b,y,aux,ni,res)
                  endif  
               else
                  if((u-ulim)*(ulim-cx).ge.0.d0)then
                     u = ulim
                     fu = estimv3(u,n,b,y,aux,ni,res)
                  else
                     u = cx + gold*(cx-bx)
                     fu = estimv3(u,n,b,y,aux,ni,res)
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
      double precision function golden3(ax,bx,cx,tol,xmin,n,b,y,aux)

c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
      
      double precision y(npmax,npmax)
      double precision ax,bx,cx,tol,xmin,b(npmax)
      double precision r,c,aux,res
      parameter (r=0.61803399d0,c=1.d0-r)
      double precision f1,f2,x0,x1,x2,x3,estimv3
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
         f1 = estimv3(x1,n,b,y,aux,ni,res)
         f2 = estimv3(x2,n,b,y,aux,ni,res)
         
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
               f2 = estimv3(x2,n,b,y,aux,ni,res)
               
c         write(*,*)'f2 DANS golden',f2,f1
c         stop
            else
               x3 = x2
               x2 = x1
               x1 = r*x2+c*x0
               f2 = f1
c         write(*,*)'f2 else',f1
               f1 = estimv3(x1,n,b,y,aux,ni,res)
c         write(*,*)'f1 else',f1
c         stop
            endif
            go to 1
          endif
          if(f1.lt.f2)then
             golden3 = f1
             xmin = x1
          else
             golden3 = f2
             xmin = x2
          endif
          return
          end


c========================          ESTIMV         ===================

      double precision function estimv3(k00,n,b,y,aux,ni,res)

c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
 
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
      integer  nsujet,nsim,nva,ndate,nst
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,nst
c*****dace1 
      double precision date(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,zi
      
c*****dace3
      double precision :: pe
      integer :: effet,nz1,nz2,correl
      common /dace3/pe,effet,nz1,nz2,correl
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
      k0(2)=0.d0
c      write(*,*)'DEBUT ESTIMV',k0,effet,n 

c      stop
      call marq983(k0,b,n,ni,v,res,ier,istop)
      
      
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
         
         call test3(bh,ut,dut,k0,n,aux,v,y)
         estimv3 = - ((res-pe)) - aux
c         write(*,*)'estimv',k0,-aux,ni,estimv3
      else
         aux = -n
c     write(*,*)'estimv2',k0,-aux,ni
      endif
c      stop
      return
      end
      
c=================calcul de la hessienne  et de omega  ==============
      subroutine test3(b,ut,dut,k0,n,res,v,y)
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
 
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
      integer  nsujet,nsim,nva,ndate,nst
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,nst
c*****
         do 10 i = 1,n
            do 5 j = 1,n
               hess(i,j) = 0.d0 
 5         continue
 10      continue
   
 
         do 20 i = 1,n
            do 15 j = i,n
               call mat3(hess(i,j),ut,dut,i,j,n)
c               write(*,*)'hess test',hess(i,j),i,j
 15         continue
 20      continue
        do 40 i = 2,n
            do 35 j = 1,i-1
               hess(i,j)=hess(j,i)
 35         continue
 40      continue


         call calcomeg3(n,omeg)

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
         call ludcmp3(hess,n,np,indx,d)

         do 113 j=1,n
            call lubksb3(hess,n,np,indx,y(1,j))
 113     continue

c         write(*,*)'hess test',y(1,1)

         tra = 0.d0
         do 150 i=1,n
            do 140 j=1,n
               tra = tra + y(i,j)*hessh(j,i)
c         write(*,*)'test tra',tra,y(i,j),hessh(j,i),i,j
 140        continue
 150     continue
   
c         write(*,*)'test',k0(1),-tra
         res = (tra)

         end

c======================  LUBKSB  ======================================
      subroutine lubksb3(a,n,np,indx,b)


c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
 
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
       subroutine ludcmp3(a,n,np,indx,d)

c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
 

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


            if (aamax.eq.0.d0) then
c              pause 'mat3(rice singuliere'
c
c  do not remove!
             write(*,*) 'singular matrix'
            end if 
 
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
c=======================  CALOMEG  ===========================
      subroutine calcomeg3(n,omeg)
c        remplissage de la mat3(rice omega n*n
c          elle a 7 diagonales
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
 

      double precision omeg(npmax,npmax)
      integer n
      double precision calc003,calc013,calc023
      integer i,j
c*****dace1 
      double precision date(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,zi

c*****pen1
      double precision ,dimension(npmax) ::m3m3,m2m2,m1m1,mmm,m3m2
      common /pen1/m3m3,m2m2,m1m1,mmm,m3m2
c*****pen2
      double precision ,dimension(npmax) ::m3m1,m3m,m2m1,m2m,m1m
      common /pen2/m3m1,m3m,m2m1,m2m,m1m

c**************************
      do 5 i=1,n
         do 3 j=1,n
            omeg(i,j)=0.d0
 3       continue
 5    continue
   
      omeg(1,1)=calc003(1,n)
      omeg(1,2)=calc013(1,n)
      omeg(1,3)=calc023(1,n)
      omeg(1,4)=m3m(1)
      omeg(2,1)=omeg(1,2)
      omeg(2,2)=calc003(2,n)
      omeg(2,3)=calc013(2,n)
      omeg(2,4)=calc023(2,n)
      omeg(2,5)=m3m(2)
      omeg(3,1)=omeg(1,3)
      omeg(3,2)=omeg(2,3)
      omeg(3,3)=calc003(3,n)
      omeg(3,4)=calc013(3,n)
      omeg(3,5)=calc023(3,n)
      omeg(3,6)=m3m(3)
      do 10 i=4,n-3
         omeg(i,i-3)=omeg(i-3,i)
         omeg(i,i-2)=omeg(i-2,i)
         omeg(i,i-1)=omeg(i-1,i)
         omeg(i,i)=calc003(i,n)
         omeg(i,i+1)=calc013(i,n)
         omeg(i,i+2)=calc023(i,n)
         omeg(i,i+3)=m3m(i)
 10   continue   
      omeg(n-2,n-5)=omeg(n-5,n-2)
      omeg(n-2,n-4)=omeg(n-4,n-2)
      omeg(n-2,n-3)=omeg(n-3,n-2)
      omeg(n-2,n-2)=calc003(n-2,n)
      omeg(n-2,n-1)=calc013(n-2,n)
      omeg(n-2,n)=calc023(n-2,n)
      omeg(n-1,n-4)=omeg(n-4,n-1)
      omeg(n-1,n-3)=omeg(n-3,n-1)
      omeg(n-1,n-2)=omeg(n-2,n-1)
      omeg(n-1,n-1)=calc003(n-1,n)
      omeg(n-1,n)=calc013(n-1,n)
      omeg(n,n-3)=omeg(n-3,n)
      omeg(n,n-2)=omeg(n-2,n)
      omeg(n,n-1)=omeg(n-1,n)
      omeg(n,n)=calc003(n,n)

      end


c====================  mat3(  ==================================
      subroutine mat3(res,ut,dut,k,l,n)

c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
 

      double precision res,dut(ndatemax),ut(ndatemax)
      integer k,l,j,ni,n
      double precision res1,msp3,aux2
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
      integer  nsujet,nsim,nva,ndate,nst
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,nst
c**

             
c---------- calcul de la hessienne ij ------------------
          res = 0.d0
          res1 = 0.d0
          do 10 i=1,nsujet
             if(c(i).eq.1)then  !event
                u2 = dut(nt1(i)) 
c                write(*,*)'mat3( u2',u2,nt1(i),i
                do 6 j = 2,n-2
                   if((date(nt1(i)).ge.zi(j-1)).and.
     &                  (date(nt1(i)).lt.zi(j)))then
                      ni = j-1
c                     write(*,*)'mat3( ni',ni
                   endif
 6              continue 
                if(date(nt1(i)).eq.zi(n-2))then
                   ni = n-2
                endif   
c-------attention numero spline 
                aux2 = msp3(nt1(i),ni,k)*msp3(nt1(i),ni,l)
c                write(*,*)'mat3( aux2',aux2,nt1(i),ni,k
                if (u2.le.0.d0)then
                   res1 = 0.d0
                else   
                   res1 = - aux2/(u2*u2)
                endif  
c                write(*,*)'mat3(',res1
c                stop
             else !censure  
                res1 = 0.d0
             endif 
          res = res + res1
 10    continue   
       
       end

c==========================  MSP   ==================================
          double precision function msp3(i,ni,ns)
 
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************

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

             msp3 = val
             return
             end
c==========================   SP   ==================================
          double precision function sp3(i,ni,ns)
 
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
  
             integer ni,ns,i
             double precision val,msp3
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
                   val = ((zi(ni)-zi(ni-1))*msp3(i,ni,ns))*0.25d0
                else
                   if (ni.eq.ns-1)then
                      val = ((zi(ni)-zi(ni-2))*msp3(i,ni,ns)+
     &                (zi(ni+3)-zi(ni-1))*msp3(i,ni,ns+1))*0.25d0
                   else
                      if(ni.eq.ns)then
                         val = ((zi(ni)-zi(ni-3))*msp3(i,ni,ns)+
     &                       (zi(ni+2)-zi(ni-2))*msp3(i,ni,ns+1)
     &              +(zi(ni+3)-zi(ni-1))*msp3(i,ni,ns+2))*0.25d0
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
                   val = (date(i)-zi(ni))*msp3(i,ni,ns)*0.25d0
             else  
             if(ni.eq.ns-2)then
                   val = ((date(i)-zi(ni-1))*msp3(i,ni,ns)+
     &             (zi(ni+4)-zi(ni))*msp3(i,ni,ns+1))*0.25d0
             else   
                if (ni.eq.ns-1)then
                   val =((date(i)-zi(ni-2))*msp3(i,ni,ns)+
     &             (zi(ni+3)-zi(ni-1))*msp3(i,ni,ns+1)
     &             +(zi(ni+4)-zi(ni))*msp3(i,ni,ns+2))*0.25d0
                else
                   if(ni.eq.ns)then
                      val =((date(i)-zi(ni-3))*msp3(i,ni,ns)+
     &             (zi(ni+2)-zi(ni-2))*msp3(i,ni,ns+1)
     &             +(zi(ni+3)-zi(ni-1))*msp3(i,ni,ns+2)
     &             +(zi(ni+4)-zi(ni))*msp3(i,ni,ns+3))*0.25d0
                   else
                      val = 1.d0
                   endif
                endif
             endif
             endif
          endif 
          endif
             sp3 = val
             return
             end
c================

c=========================  CALC00  =========================
          double precision function calc003(j,n) 
  
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************

          double precision part
          integer j,n

c*****pen1
      double precision ,dimension(npmax) ::m3m3,m2m2,m1m1,mmm,m3m2
      common /pen1/m3m3,m2m2,m1m1,mmm,m3m2
c*****pen2
      double precision ,dimension(npmax) ::m3m1,m3m,m2m1,m2m,m1m
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

             calc003 = part
             return
             end
c=========================  CALC01  =========================
      double precision function calc013(j,n)

c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************

      double precision part
      integer j,n

c*****pen1
      double precision ,dimension(npmax) ::m3m3,m2m2,m1m1,mmm,m3m2
      common /pen1/m3m3,m2m2,m1m1,mmm,m3m2
c*****pen2
      double precision ,dimension(npmax) ::m3m1,m3m,m2m1,m2m,m1m
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

             calc013 = part
             return
             end
c=========================  CALC02  =========================
          double precision function calc023(j,n)
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************


          double precision part
          integer j,n

c*****pen1
      double precision ,dimension(npmax) ::m3m3,m2m2,m1m1,mmm,m3m2
      common /pen1/m3m3,m2m2,m1m1,mmm,m3m2
c*****pen2
      double precision ,dimension(npmax) ::m3m1,m3m,m2m1,m2m,m1m
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

             calc023 = part
             return
             end


c===============================    MARQ98AUX =========================
c pour la maximisation auxiliaire
   	 subroutine marq98AUX3(k0,b,m,ni,v,rl,ier,istop)
c
c
c  fu = mat3(rice des derivees secondes et premieres
c
c  istop: raison de l'arret
c  1: critere d'arret satisfait (prm=ca, vraisblce=cb, derivee=dd)
c  2: nb max d'iterations atteints
c  3: 1 mais echec inversion mat3(rice d'info (ier=1)
C  4: non amelioration vraisblce apres searpas
c   
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
   
         integer :: m,ni,nql,i,ii,nfmax,idpos,ier,istop,igrad,j
         integer :: ncount,id,jd,i0
           
         double precision :: da,dm,ga,tr
         double precision :: ca,cb,epsa,epsb,rl
         double precision :: funcpaAUX,step,eps,epsd
         double precision :: vw,fi,maxt3,z
         double precision :: rl1,th,ep,dd,auxmax

         double precision , dimension(npmax*(npmax+3)/2)::V
         double precision , dimension(npmax*(npmax+3)/2)::vnonpen
         double precision , dimension(npmax)::b,delta
         double precision , dimension(npmax*(npmax+3)/2) :: fu
         double precision , dimension(2) :: k0
         double precision , dimension(npmax) :: bh,b1
         double precision , dimension(2) :: zero

c***** nmax
         integer :: nmax
         common /nmax/nmax
c*****dace2
      double precision , dimension(nsujetmax) ::t0,t1
      double precision :: ncsrl
      integer , dimension(nsujetmax) ::c
      integer, dimension(nsujetmax) :: nt0,nt1
      integer :: nsujet,nsim,nva,ndate,nst
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,nst
c*****dace3
      double precision :: pe
      integer :: effet,nz1,nz2,correl
      common /dace3/pe,effet,nz1,nz2,correl
c*****auxig
      integer :: auxig
      common /auxig/auxig

c***********************************************
c***** dace1 ************************************

c      double precision , dimension(ndatemax) ::date
c      double precision , dimension(-2:npmax) :: zi
c     common /dace1/date,zi

      zero(1)=0.d0
      zero(2)=0.d0
      id=0
      jd=0
      z=0.d0
      th=1.d-5
      eps=1.d-7
      epsa=1.d-4!1.d-3
      epsb=1.d-4!1.d-3
      epsd=1.d-3!1.d-2
      nfmax=m*(m+1)/2
      ca=epsa+1.d0
      cb=epsb+1.d0
      rl1=-1.d+10
      ni=0
      istop=0
      da=0.01d0
      dm=5.d0
c     nql=nbre de seconds membres pour chole
      nql=1

 10   continue
c     
c      write(6,*)
c      write(*,*)'** ITERATION :**', ni 
c         write(*,*)'** Parametres:'
c         write(*,*)'/////dans MARQ98AUX',(b(i),i=1,m)
c         write(6,*)'da',da,' ga=',ga
 
      z=0.d0
      i0=0 
      rl=funcpaAUX(b,m,i0,z,i0,z,k0)
      rl1=rl 

c      write(*,*)'log Vrais rl1 AUX,nsujet',rl1,auxig,nsujet
c      stop
c      if (RL1.gt.0)stop
      
      call derivaAUX3(b,m,v,rl,k0)

c      write(*,*)'FIN derivaaux',(v(i),i=1,5)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      do 144 i=1,m*(m+1)/2
c          ii=i*(i+1)/2
c           write(*,*) 'toutes derivees 2nds marq ',i,' = ',v(i)
c  144  continue 
c       ii=0
c       do 155 i=m*(m+1)/2+1,m*(m+3)/2
c          ii=ii+1
c        write(*,*) 'derivee premiere % parametre, marq',ii,' = ',v(i)
c  155  continue
cccccccccccccccccccccccc
c      stop
          dd = 0.d0
          do 13 i=m*(m+1)/2+1,m*(m+3)/2
c             write(*,*)'--------------dd v(i)',v(i),i
             dd = dd + v(i)*v(i)
 13       continue 
          dd=dd/dabs(RL)

c          pause
c          write(6,*) (b(ii),ii=1,m)
          if (ca.lt.epsa.and.cb.lt.epsb.and.dd.lt.epsd) goto 100
c          if (ca.lt.epsa.and.cb.lt.epsb) goto 100
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
c          write(*,*)'avant dcholeaux',fu(nfmax+1),m,nql,idpos
          call dcholeaux(fu,m,nql,idpos) !ok idem
c         write(*,*)'apres dcholeaux',(fu(i),i=1,nfmax+m),m,nql,idpos
c          stop
          if (idpos.ne.0) then  
             write(6,*) 'echec  choleskiAUX',b(1),b(2),ncount 

             if(b(1).lt.-1.d0.or.b(1).gt.1.d0.or.
     &            b(2).lt.-1.d0.or.b(2).gt.1.d0) then
             b(1)=0.5d0
             b(2)=0.5d0
             goto 110
             endif
             
             
c             write(*,*)'ncount',ncount
c             stop
c modif ligne suivante le 25/10/99 (da*5 et non *25 a chaque iter)
c             da=dm*da
             ncount=ncount+1
             if (ncount.le.3.or.ga.ge.1.d0) then
                da=da*dm
             else
                ga=ga*dm
                if (ga.gt.1.d0) ga=1.d0
             endif
CCC   ajout 240305 JEREMIE
               if (ncount > 10) then
                  fu=0.d0
                  do i=1,m
                     ii=i*(i+1)/2
                     fu(ii)=1
                  end do
                  print*,"pas de choleski, on a pris I"
                  return
               end if
ccccc FIN AJOUT         
             goto 400
          else
c              write(6,*) 'mat3(rice reguliere, idpos=',idpos
c              write(*,*)'ncount=',ncount,'da=',da,'ga=',ga
              do 2 i=1,m
                 delta(i)=fu(nfmax+i)
                 b1(i)=b(i)+delta(i)
2             continue
c              write(6,*) '** ** avant func.....delta  ',  delta
c              write(6,*) '**'
c              write(6,*) '** ** avant func.....B:',b
c              write(6,*) '** ** avant func..',b1,m,id,z,jd,z,k0
              rl=funcpaAUX(b1,m,id,z,jd,z,k0)
c              write(6,*) 'rl1 =',rl1,' rl =',rl
c              stop
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
c           
c----------rajouté juin 2006
c sinon on a des maxt3(delta,m) < 0 donc le passage au log dans searpas pose pb)
c attention ne pas mettre directement maxt3() dans le calcul de vw
           auxmax=dabs(maxt3(delta,m))
           if(auxmax.eq.0.d0) then 
              vw=1.D-5
           else
              vw=th/auxmax
           endif
c----------
c           
           step=dlog(1.5d0)
c            write(*,*)'vw ',vw,'th maxt ',th,maxt3(delta,m)
c           write(*,*)'delta',(delta(i),i=1,m)
c           pause
c           write(*,*) 'searpas'
           call SEARPASAUX3(vw,step,b,bh,m,delta,fi,eps,k0) 
c           write(*,*) 'apres searpas'
           rl=-fi
           IF(rl1.gt.rl) then
c              write(6,*) 'non amelioration vraisblce apres searpasAUX'
c              write(6,*) 'ancienne rl1AUX =',rl1,' nouvelle rlAUX =',rl
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

 800       cb=dabs(rl1-rl)
           ca=0.d0
           do 5 i=1,m
              ca=ca+delta(i)*delta(i)
 5         continue

c           write(6,*) 'ca =',ca,' cb =',cb,' dd =',dd
    
           do 6 i=1,m
              b(i)=b(i)+delta(i)
c              write(*,*)'delta',delta(i),i
 6         continue

           ni=ni+1
           if (ni.gt.nmax) then
              istop=2
              write(6,*) 'nombre iteration max atteinte'
              goto 110
           end if
           goto 10
c********************
c     inversion mat3(rice d'informat3(ion
c********************
100      continue
           istop=1

c           write(*,*)'ca,cb,dd final',ca,cb,dd
           ep=10.d-10
      
           call DSINV3(v,m,ep,ier)
           if (ier.eq.-1) then
             write(6,103)          
103       format(1x,'echec inversion mat3( informat3(ion ds marq98aux')

          istop=3
       endif

 110   continue
c     stop
       return
       end
c========fin marq98AUX

      
c=================================DERIVAAUX =======================

      subroutine derivaAUX3(b,m,v,rl,k0)
            
      integer ::i0,iun,m,m1,ll,i,k,j
      double precision :: funcpaAUX,thn,th,z,rl,vl
      double precision :: th2
      double precision ,dimension(2):: k0
      double precision ,dimension(m):: fcith
      double precision ,dimension(m):: b
      double precision , dimension(m*(m+3)/2)::V
      
c     v:mat3(rice d'informat3(ion+score
c     calcul de la derivee premiere
c     
c      print*,''
c      print*,'entree deriva'
c      print*,''
      th=1.d-5!1.d-4!
      thn=-th
      th2=th*th
      z=0.d0
      i0=0
      iun =1
      rl=funcpaAUX(b,m,iun,z,iun,z,k0)
      do 2 i=1,m
         fcith(i)=funcpaAUX(b,m,i,th,i0,z,k0)
c         print*,'*** fcith :',fcith(i),b(i),i,m
 2    continue

      k=0
      m1=m*(m+1)/2
      ll=m1
      do 1 i=1,m
         ll=ll+1
         vl=(fcith(i)-funcpaAUX(b,m,i,thn,i0,z,k0))/(2.d0*th)
c       print*,'*funcpaAUX :',vl,b(i),i,m
         v(ll)=vl
         do 1 j=1,i
            k=k+1
            v(k)=-(funcpaAUX(b,m,i,th,j,th,k0)-fcith(j)-fcith(i)+rl)/th2 
c            print*,'*** v(k)',v(k),k,j,i
 1       continue
c         print*,'FIN DERIVA$$$$$$'
c         stop
         return
         end
c fin derivaAUX



c================================  SEARPAS joly    ==============================

      SUBROUTINE SEARPASAUX3(VW,STEP,B,BH,M,DELTA,FIM,EPSV,k0)
C
C  MINIMISATION UNIDIMENSIONNELLE
C
       INTEGER :: I,M
       DOUBLE PRECISION :: VLW,VLW1,VLW2,VLW3,VW,VM
       DOUBLE PRECISION :: FI1,FI2,FI3,FIM,EPSV
       DOUBLE PRECISION :: STEP
       DOUBLE PRECISION , dimension(2):: k0
       DOUBLE PRECISION , dimension(M):: B,BH,DELTA 
C
       VLW1=DLOG(VW)
       VLW2=VLW1+STEP
       CALL VALFPAAUX(VLW1,FI1,B,BH,M,DELTA,k0)
c       write(*,*)'**dans searpasaux',VLW1,vw
       CALL VALFPAAUX(VLW2,FI2,B,BH,M,DELTA,k0)
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
          CALL VALFPAAUX(VLW1,FI1,B,BH,M,DELTA,k0)   
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
          CALL VALFPAAUX(VLW1,FI1,B,BH,M,DELTA,k0)
          IF(FI1.GT.FI2) GO TO 50
c          IF (dabs(FI1-FI2).LT.EPSV) THEN
          IF (FI1.eq.FI2) THEN
             FIM=FI2
             VM=VLW2
             GO TO 100
          ENDIF
 20    CONTINUE
C
C  PHASE 2 APPROXImat3(ION PAR QUADRIQUE
C
50     CONTINUE
C
C  CALCUL MINIMUM QUADRIQUE
C
         VM=VLW2-STEP*(FI1-FI3)/(2.d0*(FI1-2.d0*FI2+FI3))
         CALL VALFPAAUX(VM,FIM,B,BH,M,DELTA,k0)
         IF (FIM.LE.FI2) GO TO 100
         VM=VLW2
         FIM=FI2
100   CONTINUE
      VW=DEXP(VM)
      RETURN

      END


c fin searpasaux
   
c========================   VALFPAAUX   ==============================

        subroutine valfpaAUX(vw,fi,b,bk,m,delta,k0)
        integer :: m,i
        double precision :: vw,fi
        double precision :: funcpaAUX,z        
        double precision ,dimension(2) :: k0
        double precision ,dimension(m) :: b,bk,delta

        z=0.d0
        do 1 i=1,m
           bk(i)=b(i)+dexp(vw)*delta(i)
c        print*,'***  VALFPAAUX ***',b(i),dexp(vw),vw,delta(i),i
       
 1      continue
c     
        fi=-funcpaAUX(bk,m,1,z,1,z,k0)
c
         return
         end   
       
c fin valfpaAUX

c================================  DCHOLE  ===========================
      subroutine dcholeaux(a,k,nq,idpos)

c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************
      
      integer :: k,nq,i,ii,i1,i2,i3,m,is,j,k2,jmk
      integer :: ijm,irm,jji,jjj,l,jj,iil,jjl,il,idpos
      double precision ,dimension(npmax*(npmax+3)/2) :: a
      double precision :: term,xn,diag,p
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
c     calcul des elements de la mat3(rice
      do 13 i=1,k
      ii=i*(i+1)/2
c       elements diagonaux
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
c=====================

c========================    DEBUT FUNCPAAUX       ====================
      double precision function funcpaAUX(b,np,id,thi,jd,thj,k0)
 
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************


      integer :: n,np,id,jd,i,j,k,vj,cptg
      integer :: ig,ip
      integer , dimension(ngmax) :: cpt
      double precision , dimension(npmax) :: bhaux,b

      double precision :: thi,thj

      double precision ::res,vet
      double precision , dimension(2) :: k0

c ************************* les commons :
c*****groupe
      integer , dimension(nsujetmax) :: g
      integer , dimension(ngmax) :: nig! nb de sujet par groupe
      common /gpe/g,nig
c*****dace2
      double precision , dimension(nsujetmax) ::t0,t1
      double precision :: ncsrl
      integer , dimension(nsujetmax) ::c
      integer, dimension(nsujetmax) :: nt0,nt1
      integer :: nsujet,nsim,nva,ndate,nst
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,nst
c*****dace4
      integer , dimension(nsujetmax) :: stra
      common /dace4/stra
c*****ut1ut2
      double precision , dimension(ndatemax) :: dut1,dut2
      double precision , dimension(0:ndatemax) :: ut1,ut2
      common /ut1ut2/dut1,dut2,ut1,ut2
c**** betaaux
      double precision , dimension(nvarmax) :: betaaux
      common /betaaux/betaaux
c*****ve1
      double precision , dimension(nsujetmax,nvarmax):: ve,ve2
      common /ve1/ve,ve2
c*****auxig
      integer :: auxig
      common /auxig/auxig
c*******sigma2tau2rho
      double precision :: sigma2,tau2,rho,cov
      common /sigma2tau2/sigma2,tau2,rho,cov
c****** u_tilde
      double precision  :: u_tilde,v_tilde
c      common /utilde/u_tilde,v_tilde
c****** derivanal
      double precision , dimension(ngmax)::res1,res2,res3,res4
      double precision , dimension(ngmax)::res5,res6,res8
c      common /derivanal/res1,res2,res3,res4,res5,res6,res8
c******************
c==============================================
c================POUR UN GROUPE AUXIG donné !
c==============================================
c ici que 2 parametres  à estimer : u et v
c pour des valeurs fixées des splines, des beta et sigma2 et tau2 
c sans pénalisation
   
      res3=0.d0
      res4=0.d0
      res5=0.d0
      res6=0.d0
      res8=0.d0

      do 3 i=1,np
         bhaux(i)=b(i)
c         write(*,*)'%% ENTREE funcapaaux: bhaux ',bhaux(i),thi,thj,i,np
 3    continue 
c      stop

      if (id.ne.0) bhaux(id)=bhaux(id)+thi 
      if (jd.ne.0) bhaux(jd)=bhaux(jd)+thj    
         
      u_tilde=bhaux(1)
      v_tilde=bhaux(2)
c--------------------------------------------------------
c----------calcul de la vraisemblance ------------------
c---------------------------------------------------------

      do 101 k=1,nsujet
         if(nva.gt.0.and.g(k).eq.auxig)then
            vet = 0.d0 
            do 1912 ip=1,nva
               vet =vet + betaaux(ip)*ve(k,ip)
 1912       continue
            vet = dexp(vet)
         else
            vet=1.d0
         endif
         
         if(g(k).eq.auxig)then
            if(c(k).eq.1)then
               res3(auxig) = res3(auxig)
     &              +u_tilde+v_tilde*ve(k,1)
c                      write(*,*)'))))))))res3',res3(auxig),ve(k,1),
c     &                     u_tilde,v_tilde,auxig
               res8(auxig) = res8(auxig)
     &              +ve(k,1)
            endif
            if(stra(k).eq.1)then
               res4(auxig) = res4(auxig)
     &              +ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve(k,1))
               res5(auxig) = res5(auxig)
     &              +ut1(nt1(k))*vet*
     &              dexp(u_tilde+v_tilde*ve(k,1))*ve(k,1)
               res6(auxig) = res6(auxig)
     &              +ut1(nt1(k))*vet*
     &              dexp(u_tilde+v_tilde*ve(k,1))*(ve(k,1))**2
            endif
            if(stra(k).eq.2)then
               res4(auxig) = res4(auxig)
     &              +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve(k,1))
               res5(auxig) = res5(auxig)
     &              +ut2(nt1(k))*vet*
     &              dexp(u_tilde+v_tilde*ve(k,1))*ve(k,1)
               res6(auxig) = res6(auxig)
     &              +ut2(nt1(k))*vet*
     &              dexp(u_tilde+v_tilde*ve(k,1))*(ve(k,1))**2
            endif
           
         endif 
 101  continue 
c              write(*,*),'fin funcpaux res4,res5',
c     &           res4(auxig),res5(auxig),res6(auxig),auxig
         
      res= !-res2(auxig)
     & -res3(auxig) + res4(auxig) 
     &  +0.5d0*(((u_tilde)**2)/sigma2+((v_tilde)**2)/tau2
     &      -2.d0*u_tilde*v_tilde*cov/(sigma2*tau2))
     &        /(1.d0-(cov**2)/(sigma2*tau2))         !-ka


      funcpaAUX= -res
c      write(*,*)'%%% SORTIE funcapaaux',-res,! =K(.)
c     &            res2(auxig),res3(auxig),res4(auxig),sigma2,tau2,
c     &            invD(1,1),invD(1,2),invD(2,2),
c     &            betaaux(1),u_tilde,v_tilde,auxig
c      stop
      return
      end


c=================================    DERIVATH  =======================

      subroutine derivaTH(b,m,v,rl,k0)
          
c*************** SAME **********************************************
         integer , parameter ::npmax=50,NSUJETMAX=20000,nvarmax=50
         integer , parameter ::ngmax=1000,nboumax=1000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1000
c******************************************************************

      integer ::i0,iun,m,m1,ll,i,k,j
      double precision :: funcpa3,thn,th,z,rl,vl
      double precision :: th2
      double precision ,dimension(2):: k0
      double precision ,dimension(npmax):: fcith
      double precision ,dimension(npmax):: b
      double precision , dimension(npmax*(npmax+3)/2)::V
      
c     v:mat3(rice d'informat3(ion+score
c     calcul de la derivee premiere
c     
c      print*,'entree deriva'
      th=1.d-4!1.d-5
      thn=-th
      th2=th*th
      z=0.d0
      i0=0
      iun =1
      rl=funcpa3(b,m,iun,z,iun,z,k0)
      do 2 i=1,m
         fcith(i)=funcpa3(b,m,i,th,i0,z,k0)
c         print*,'*** fcith :',fcith(i),b(i),i,m
 2    continue
      
      k=0
      m1=m*(m+1)/2
      ll=m1
      do 1 i=1,m
         ll=ll+1
         vl=(fcith(i)-funcpa3(b,m,i,thn,i0,z,k0))/(2.d0*th)
c         print*,'*** - funcpa :', funcpa3(b,m,i,thn,i0,z,k0),b(i),i,m
         v(ll)=vl
         do 1 j=1,i
            k=k+1
            v(k)=-(funcpa3(b,m,i,th,j,th,k0)-fcith(j)-fcith(i)+rl)/th2 
c            print*,'*** v(k)',v(k),k,j,i
 1       continue
         
         return
         end


c====================
