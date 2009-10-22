c rem : pour le faire avec du gap time : utiliser recurr_gap.txt
c NEW : 03/11/2005
c Nous avons un seul temps de censure et de deces par sujet
c pour garder la meme base de temps pour les events recurr et le deces 
c on decide de travailler sur les "calendars time" de type ANdersen -Gill 
c (// left truncation)
c l'ecriture de la vraisemblance a été reprise, 
c il ne fait intervenir qu une seule integrale
c on ne limite pas le nb de données recurrentes, 
c ces temps sont générés tant que l on ne dépasse pas les temps de cens ou deces
c========================================================


c 17 / 10 / 2005
c legere modification du calcul de func1 pour eviter les pbs numeriques
c------
c JOINT FRAILTY MODEL : ex: recurrent events and death
c dans l esprit de Liu, Biometrics,2004
c mais on propose plusieurs echelles de temps (gap, calendar time),
c et donc comme precedemment plusieurs types de troncature,
c et plusieurs vraisemblances pour estimer les parametres
c les modeles sont:
c** recuurent : \lambda(t_ij|Z_i)=\lambda_0(t_ij) Z_i exp(beta X_ij)
c** event     : \lambda(t_ij*|Z_i)=\lambda_0(t_ij*) Z_i^alpha exp(beta X_ij)
c** en fait t_ij et t_ij* sont identiques, ce qui change c les indicateurs de censure

c deux bases de splines identiques mais coefficientsdes splines differents
c  pour les deux risques de base : ici NST = 2
c on ne traite plus ici la stratification
c attention : les effets des var expli sont differents pour recurr event et death
c*****************************************************************************
c programme modifié le 21 mars 2005 
c on a la possibilité de choisir deux formaulations pour la vraisemblance :
c soit la vraisemblance avec troncature a gauche : cf LIDA2003
c soit l'ancienne version de la troncature, adaptée aux données récurrentes de type andersen gill (calendar time), qui permet aussi de traiter des variables explicatievs dépendantes du temps .

c********************************************************************
c********  version F77 pour le passer sous R ************************
c********************************************************************
c****le 25/01/2005 modification
c BANDES DE CONFIANCES avec H-1
c dans marq / on tient compte du changement de variable 
c effectué sur les parametres des splines .
c cf v1 et hess pour le resultat : 
c inversion de la matrice comprenant uniquement les para des splines
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
c bis le 7/11/2001, modification de 
c erreur sur ddl d ordre 3
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


       subroutine frailpenalJoint(nsujetAux,ngAux,icenAux,
     &      nzAux,k0,tt0Aux,tt1Aux,icAux,groupeAux,
     &      tt0dcAux,tt1dcAux,icdcAux,nvaAux,vaxAux,
     &      nvadcAux,vaxdcAux,AGAux,noVar,maxitAux,irep1,
     &      np,b,H_hessOut,HIHOut,resOut,x1Out,lamOut,suOut,
     &      x2Out,lam2Out,su2Out,ni,cpt,cpt_dc,ier)


      
c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************


       integer  groupe,ij,kk,j,k,nz,n,np,cpt,cptcens,cpt_dc,ii,iii,iii2
       integer  cptstr1,cptstr2,trace,trace1,trace2,ver
       integer  i,ic,icdc,ni,ier,istop,ef
       integer  cptni,cptni1,cptni2,nb_echec,nb_echecor
       integer  id,cptbiais,l
       integer  m,idum
       integer filtre(nvarmax), filtre2(nvarmax) 
       integer  cptaux , cptauxdc 
       integer  nb0recu
       

       real vax(nvarmax),vaxdc(nvarmax)
       double precision tt0,tt0dc,moyrecu
       double precision tt1,tt1dc
      
       double precision h
       double precision ro,wres,csi,csi1
       double precision res,min,mindc,max,maxdc,maxt,maxtdc
       double precision bi,bs,wald
       double precision moyvar,moyse,moyse_cor,str
       double precision varxij,eca,varsmarg,smoy,smoyxij
       double precision pe1,pe2,moy_peh0,moy_peh1,lrs
       double precision BIAIS_moy

       double precision aux(2*nsujetmax)
       double precision v((npmax*(npmax+3)/2))
       double precision k0(2),res01(2)
       double precision b(npmax)
       double precision I1_hess(npmax,npmax),H1_hess(npmax,npmax)
       double precision I2_hess(npmax,npmax),H2_hess(npmax,npmax)
       double precision HI1,HI2(npmax,npmax)
       double precision HIH(npmax,npmax),IH(npmax,npmax),HI(npmax,npmax)
       double precision BIAIS(npmax,1)


       character*20 nomvarl
       character*20 nomvar(nvarmax),nomvar2(nvarmax)
       character*20 donnees
       character*20 death
       character*24 ficpar
       character*20 fich1
       character*20 fich2
       character*20 fich3
       character*20 fich4
       character*20 fich1b
       character*20 fich2b
       character*20 fich3b
       character*20 fich4b
      character (len=20)::fich1c
      character (len=20)::fich2c
      character (len=20)::fich3c
      character (len=20)::fich4c
       character*20 dateamj
       character*20 zone
       character*20 heure1
       character*20 heure2
       integer values(8)





c******************************************   Add JRG January 05

         double precision  x1,x2,lam,lbinf,lbsup,margi
         double precision  lam2,lbinf2,lbsup2,margi2
         integer nsujetAux,ngAux,icenAux,nstAux,effetAux,nzAux,nvaAux
         integer nvadcAux 
         double precision tt0Aux(nsujetAux),tt1Aux(nsujetAux),
     &                    tt0dcAux(ngAux),tt1dcAux(ngAux)

         integer icAux(nsujetAux),groupeAux(nsujetAux),ss,sss,
     &           icdcAux(ngAux)   

         double precision strAux(nsujetAux)
         double precision vaxAux(nsujetAux,nvaAux),
     &                    vaxdcAux(ngAux,nvadcAux)
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
      double precision date(ndatemax),datedc(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,datedc,zi

c*****dace2
      double precision t0dc(ngmax),t1dc(ngmax)
      double precision t0(nsujetmax),t1(nsujetmax)
      integer c(nsujetmax), cdc(ngmax)
      integer nt0(nsujetmax),nt1(nsujetmax)
      integer nt0dc(ngmax),nt1dc(ngmax)
      integer  nsujet,nva,nva1,nva2,ndate,ndatedc,nst
      common /dace2/t0,t1,t0dc,t1dc,c,cdc,nt0,nt1,nt0dc,nt1dc,nsujet
     &     ,nva,nva1,nva2,ndate,ndatedc,nst
c*****dace4
      integer  stra(nsujetmax)
      common /dace4/stra
c*****ve1
      double precision ve(nsujetmax,nvarmax)
      double precision vedc(ngmax,nvarmax)
      common /ve1/ve,vedc
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
      integer nig(ngmax)  ! nb d events recurrents par sujet
      common /gpe/g,nig

c*****mem1
      double precision mm3(ndatemax),mm2(ndatemax)
      double precision mm1(ndatemax),mm(ndatemax)
      common /mem1/mm3,mm2,mm1,mm
c %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
      integer AG
      common /andersengill/AG
c %%%%%%%%%%%%% indic ALPHA %%%%%%%%%%%%%%%%%%%%%%%%% 
      integer indic_ALPHA
      common /alpha/indic_ALPHA ! pour preciser un para en plus 
c****  theta/alpha
      double precision  theta,alpha !en exposant pour la frailty deces 
      common /thetaalpha/ALPHA,theta
c****** indicateur de troncature
      integer :: indictronq,indictronqdc ! =0 si donnees non tronquées reellement
      common /troncature/indictronq,indictronqdc

c************ FIN COMMON ***********************************
     
c nst: deux finctions de risque a estimer (meme bases de splines)
c ist: appartenance aux strates
c ib: matrice de variance de beta chapeau
c I_hess : -hessienne non inversee sur vraisemblance non penalisee
c H_hess : inverse de -hessienne  sur vraisemblance penalisee
        
CCCCCCCCCCCCCCCCC hosur9.f CCCCCCCCCCCCCCCCCCCCCCCC
c      write(*,*)'    ******************************************'
c      write(*,*)'  ****** DEBUT PROGRAMME FRAILTY.F**********'
c      write(*,*)'******************************************'

      indic_alpha=1 ! on precise que l on a un parametre en plus estimer
      
c      call date_and_time(dateamj,heure1,zone,values)                  
c      write(*,*)'Starting time: ', dateamj,heure1,zone,values
      
      lrs=0.d0
      moy_peh0=0.d0
      moy_peh1=0.d0
      
      nb_echec=0
      nb_echecor=0
      nb0recu =0
      moyrecu =0.d0             ! nb de cens
c      write(*,*)'name of the parameter file ??'
c      read(*,*)ficpar
c      write(*,*)'name',ficpar
c      ficpar='joint.inf'
c      open(2,file=ficpar)
c      open(4,file='outjoint')


c     1900	continue (pour faire plusieurs jeux de simulations)
      
c      read(2,*)nsujet
c      read(2,*)ng
c      read(2,*)nst


      nst=2

c JRG aug'07 
      nsujet=nsujetAux
      ng=ngAux
      if (noVar.eq.1) then 
        do i=1,nvaAux
         filtre(i)=0
         filtre2(i)=0
        enddo  
        nva=0  
      else
        do i=1,nvaAux
         filtre(i)=1
        end do
        do i=1,nvadcAux
         filtre2(i)=1
        enddo  
      end if  

      ni=0
      

c      write(4,*)'**************************************************'
c      write(4,*)'************ JOINT MODEL *************************'
c      write(4,*)'*** RECURRENT EVENTS and TERMINATING EVENT *******'
c      write(4,*)'**************************************************'
c      write(4,*)'** nb de groupes+ =',ng 
c      write(*,*)'** nb de groupes+ =',ng 
c      write(*,*)'** deux fonctions de risque de base  = ',nst
      

c******************************************
c******************************************
c---debut des iterations de simulations
       
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

        AG=AGAux
        maxit=maxitAux 

        istop=0




C**************************************************
C**************************************************
C********************* prog spline****************

c        read(2,*)effet
        effet=1
        res01(1)=0.d0
        res01(2)=0.d0
        
        
c     do 765 ef=0,1
c           if (effet.eq.1)then
c              write(*,*)'-----modele AVEC effet aleatoire'
c              write(4,*)'-----modele AVEC effet aleatoire'
c           else
c              write(*,*)'-----modele SANS effet aleatoire'
c              write(4,*)'-----modele SANS effet aleatoire'
c           endif


c------------  entre non fichier et nombre sujet -----
        
c         read(2,*)ver

c         nva1 = 0 ! nb de var expli pour donnees recurrentes
c         nva2 = 0 ! nb de var expli pour deces

c         if(ver.gt.0)then
c            do 44 j=1,ver
c               read(2,*)nomvarl,filtre(j),filtre2(j)
c               nva1 = nva1 + filtre(j) ! adjustment for recurrent events
c               nva2 = nva2 + filtre2(j) ! adjustment for survival
c               if(filtre(j).eq.1)then
c                  nomvar(nva1) = nomvarl
c               endif  
c               if(filtre2(j).eq.1)then
c                  nomvar2(nva2) = nomvarl
c               endif   
c 44         continue
c         endif

c         nva = nva1+nva2


c JRG aug'07 
c  Suponemos que las variables son para recurrencias y supervivencia

c  ver ahora quitado para controlar que sean ambas distintas
c          ver=nvaAux 
          nva1=nvaAux
          nva2=nvadcAux
          nva=nva1+nva2

         do 4 i=1,ng
            nig(i) = 0
 4       continue
         
c     write(4,*)'** explanatory variables for recurrent events:',nva1
c      write(*,*)'** explanatory variables for recurrent events:',nva1
c     write(4,*)'** explanatory variables for deaths:',nva2
c      write(*,*)'** explanatory variables for deaths:',nva2
         
c------------  lecture fichier -----------------------

         maxt = 0.d0
         maxtdc = 0.d0
         
         cpt = 0
         cptcens = 0
         cpt_dc = 0
         k = 0
         cptstr1 = 0
         cptstr2 = 0

c         read(2,*)donnees
c         read(2,*)death
c         write(*,*)'** Fichiers de données = ',donnees,'et ' , death
c         write(4,*)'** Fichiers de données = ',donnees,'et ' , death

c         open(9,file=donnees)
c         open(10,file=death)
         
cccccccccccccccccccccc
c pour le deces
ccccccccccccccccccccc

          do 101 k=1,ng

c         do 101 k = 1,ng !sur les groupes uniquement (= sujet)
c            read(10,*)tt0dc,tt1dc,icdc,groupe,(vaxdc(j),j=1,ver) 


           tt0dc=tt0dcAux(k)
           tt1dc=tt1dcAux(k)
           icdc=icdcAux(k)


c Not necessary JRG aug'07
c           groupe=groupeAux(k)
 
          do j=1,nva2
             vaxdc(j)=vaxdcAux(k,j)
          enddo
           
            if(tt0dc.gt.0.d0)then
               cptauxdc=cptauxdc+1
            endif              

c            if(tt0dc.gt.tt1dc)then
c               write(*,*)'FAUX : TTOdc greater than TT1dc at line ', i
c               stop
c            endif    



c------------------   deces c=1 pour données de survie
            if(icdc.eq.1)then
               cpt_dc = cpt_dc + 1
               cdc(k)=1
               t0dc(k) = tt0dc      !/100.d0
               t1dc(k) = tt1dc      ! +0.001
               iii = 0
               iii2 = 0
c                  do 66 ii = 1,ver
                  do 66 ii = 1,nva2
                     if(filtre2(ii).eq.1)then
                        iii2 = iii2 + 1
                        vedc(k,iii2) = dble(vaxdc(ii))
                     endif
 66               continue   
               else 
c------------------   censure a droite ou event recurr  c=0 
                  if(icdc.eq.0)then
                     cdc(k) = 0 
                     iii = 0
                     iii2 = 0
c                     do 88 ii = 1,ver
                     do 88 ii = 1,nva2
                        if(filtre2(ii).eq.1)then
                           iii2 = iii2 + 1
                           vedc(k,iii2) = dble(vaxdc(ii))
                        endif
 88                   continue 
                     t0dc(k) =  tt0dc 
                     t1dc(k) = tt1dc 
                  endif
               endif
               if (maxtdc.lt.t1dc(k))then
                  maxtdc = t1dc(k)
               endif


c               write(*,*)'temps to et t1',t0(k),t1(k),k
 101     continue

         k = 0
         cptstr1 = 0
         cptstr2 = 0

ccccccccccccccccccccccccccccccccccc
c pour les données recurrentes  
ccccccccccccccccccccccccccccccccccc    
   

         do 10 i = 1,nsujet     !sur les observations
            if(nst.eq.2)then
c               read(9,*)tt0,tt1,ic,groupe,(vax(j),j=1,ver)
 
               tt0=tt0Aux(i)
               tt1=tt1Aux(i)
               ic=icAux(i)
               groupe=groupeAux(i)
c               do j=1,nva
               do j=1,nva1
                 vax(j)=vaxAux(i,j)  
               enddo
             endif 
c            vax(13)=vax(13)/100.d0 !pour reduire l age

            if(tt0.gt.0.d0)then
               cptaux=cptaux+1
            endif               

c            if(tt0.gt.tt1)then
c               write(*,*)'FAUX : TTO greater than TT1 at line ', i
c               stop
c            endif           
             
c             write(*,*)'--2--',i,tt0,tt1,ic,icdc,groupe,(vax(j),j=1,ver)
              
          
c-----------------------------------------------------
            
c     essai sans troncature
c     tt0=0.
c------------------   observation c=1 pour données recurrentes
               if(ic.eq.1)then
                  cpt = cpt + 1
                  c(i)=1
                  t0(i) = tt0 !/100.d0
                  t1(i) = tt1  ! +0.001
                  t1(i) = t1(i)!/100.d0
                  g(i) = groupe
                  nig(groupe) = nig(groupe)+1 ! nb d event recurr dans un groupe
                  iii = 0
                  iii2 = 0
c                  do 6 ii = 1,ver
                  do 6 ii = 1,nva1
                     if(filtre(ii).eq.1)then
                        iii = iii + 1
                        ve(i,iii) = dble(vax(ii)) !ici sur les observations
c                        write(*,*)'** ve **', ve(i,iii),i,iii
                     endif
 6                continue   
               else 
c------------------   censure a droite  c=0 pour données recurrentes
                  if(ic.eq.0)then
                     cptcens=cptcens+1
                     c(i) = 0 
                     iii = 0
                     iii2 = 0
c                     do 8 ii = 1,ver
                     do 8 ii = 1,nva1
                        if(filtre(ii).eq.1)then
                           iii = iii + 1
                           ve(i,iii) = dble(vax(ii))
                        endif
 8                   continue 
                     t0(i) =  tt0 !/100.d0
                     t1(i) = tt1 ! +0.001
                     t1(i) = t1(i) !/100.d0
                     g(i) = groupe
                  endif
               endif
               if (maxt.lt.t1(i))then
                  maxt = t1(i)
               endif
c                write(*,*)'temps to et t1',t0(i),t1(i),i
 10          continue 
c             write(*,*)'** max',maxt,maxtdc,cptaux,cptauxdc

c             write(4,*)'** max',maxt,maxtdc,cptaux,cptauxdc


c            stop
c %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
c         read(2,*)AG
c         write(*,*)'*AG**',AG

          

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
         nsujet = i-1
         if (effet.eq.1) then
c            write(*,*)'** nombre d observations',nsujet
c            write(*,*)'** nombre de données recurrentes ',cpt
c            write(*,*)'** nombre de données censureees ',cptcens
c            write(*,*)'** nombre de deces ',cpt_dc
c            write(*,*)'** nombre de données tronquées gauche',cptaux
c            write(4,*)'** nombre d observations',nsujet
c            write(4,*)'** nombre de données recurrentes ',cpt
c            write(4,*)'** nombre de données censureees ',cptcens
c            write(4,*)'** nombre de deces ',cpt_dc
c            write(4,*)'** nombre de données tronquées gauche',cptaux
         endif




c         read(2,*)nz

         nz=nzAux
         nz1=nz
         nz2=nz
         if(nz.gt.20)then
            nz = 20
         endif
         if(nz.lt.4)then
            nz = 4
         endif
c      write(4,*)'nombre de noeuds :',nz

c***************************************************
ccc DONNEES DECES

         mindc = 0.d0
         maxdc = maxtdc

         do 151 i = 1,2*ng
            do 161 k = 1,ng
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
 161           continue   
            aux(i) = maxdc
            mindc = maxdc + 1.d-12
            maxdc = maxtdc
c            print*,'** aux **',aux(i),t1dc(i)
 151     continue

         datedc(1) = aux(1)
         k = 1
         do 171 i=2,2*ng
               if(aux(i).gt.aux(i-1))then
                  k = k+1
                  datedc(k) = aux(i)
c         print*,'** ndatedc **',datedc(k),k,aux(i),i
               endif 
 171        continue 
         ndatedc = k      
c         write(*,*)'** ndatedc,maxtdc',ndatedc,maxdc
c         stop
         

c--------------- zi- ----------------------------------

c      construire vecteur zi (des noeuds)

ccc DONNEES RECURRENTES

         min = 1.d-10
c         min = 0.d0
         aux =0.d0
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
c         write(*,*)'** ndate,maxt',ndate,maxt
c         write(*,*)'** date **',(date(i),i,i=1,ndate)

         
         zi(-2) = date(1)
         zi(-1) = date(1)
         zi(0) = date(1)
         zi(1) = date(1)
         h = (date(ndate)-date(1))/dble(nz-1)
c         write(*,*)'** h **',h,date(1),nz,dble(nz-1)
c         stop
         do 18 i=2,nz-1
            zi(i) =zi(i-1) + h
 18      continue
         
         zi(nz) = date(ndate)
         zi(nz+1)=zi(nz)
         zi(nz+2)=zi(nz)
         zi(nz+3)=zi(nz)

c         do 188 i=1,nsujet    
c         write(*,*)'date(i)',date(i),i 
c         write(*,*)'t1(i)',t1(i),i 
c         write(*,*)'t0(i)',t0(i),i 
c 188     continue
c          stop
c         write(*,*)'** zi **',(zi(i),i,i=-2,nz+3)
 
c---------- affectation nt0dc,nt1dc DECES ----------------------------

         indictronqdc=0
            do 501 k=1,ng 
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
               do 451 j=1,ndatedc
                  if(datedc(j).eq.t0dc(k))then
                     nt0dc(k)=j
                  endif
                  if(datedc(j).eq.t1dc(k))then
                     nt1dc(k)=j
                  endif
 451           continue
 501        continue 
c            write(4,*)'** nombre de sujet SANS données recurrentes '
c     &           ,nb0recu
c            write(4,*)'** nombre de données recurrentes ',moyrecu
c            moyrecu=moyrecu/ng
c            write(4,*)'** nombre moyen de données recurrentes '
c     &           ,moyrecu
c            write(*,*)'** nombre de sujet SANS données recurrentes '
c    &           ,nb0recu
c            write(*,*)'** nombre de données recurrentes ',moyrecu
c            moyrecu=moyrecu/ng
c            write(*,*)'** nombre moyen de données recurrentes '
c     &           ,moyrecu
c         write(*,*)'datedc',(datedc(j),j=1,ndatedc)
c          write(*,*)'nt0dc',(nt0dc(k),k=1,ng)
c          write(*,*)'nt1dc',(nt1dc(k),k=1,ng)  
c          stop
c---------- affectation nt0,nt1 RECURRENTS----------------------------

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
c               write(*,*)'*** var expli',(ve(i,j),j=1,nva)
 50          continue 
c         write(*,*)'date',(date(j),j=1,ndate)
c         write(*,*)'nt0',(nt0(i),i=1,nsujet)
c         write(*,*)'nt1',(nt1(i),i=1,nsujet)  
c         stop

c---------- affectation des vecteurs de splines -----------------
           
             n  = nz+2
             call vecspli1(n,ndate) 
c      write(*,*)'** fin vecspli mm3 ** ',mm3(5)
c      stop
             call vecpen1(n)
        
             np = nst*n + nva + effet + indic_alpha
c          write(*,*)'nombre total de paramètres',np
c          write(4,*)'nombre total de paramètres',np
               
c------- initialisation des parametres
                   
                   do 75 i=1,np
                      b(i)=5.d-1
 75                continue
c                   write(4,*)'Initialisation des paramètres à ',b(1)
                   b(np-nva-indic_alpha)=1.d0 ! pour theta
                   b(np-nva)=1.d0!2.2d0 ! pour alpha

C deux parametres de lissage :


c JRG aug'07 they are passed through the subroutine
c                   read(2,*)ax1
c                   read(2,*)ax2




c                   write(4,*)'kappa1',ax1
c                   write(4,*)'kappa2',ax2
c                   write(4,*)'nombre de noeuds:',nz
                   

c                   if(effet.eq.1)then
c                      read(2,*)fich3 !ef_hazard recurrent AVEC frailty
c                      read(2,*)fich4 !ef_hazard death  AVEC frailty
                   
c                      read(2,*)fich3b !ef_surv1 AVEC frailty
c                      read(2,*)fich4b !ef_surv2 AVEC frailty
                      
c                      read(2,*)fich3c !ef_cumulhazard1 AVEC frailty
c                      read(2,*)fich4c !ef_cumulhazard2 AVEC frailty
c                   endif
                   
c     close(2)


c                   k0(1) = ax1
c                   k0(2) = ax2

c       write(*,*)'lissage',k0(1)
c       write(*,*)'lissage2',k0(2),nva1,nva2,nva   


c         write(*,*),'==================================='
c         write(*,*),'== ensemble des parametres ========',k0(1),k0(2)
c         write(*,*),'==================================='

            
         call marq981(k0,b,np,ni,v,res,ier,istop,effet)
         

         call multi3(I_hess,H_hess,np,np,np,IH)
         call multi3(H_hess,IH,np,np,np,HIH)


c         write(*,*) "....", I_hess(1,1)
c         write(*,*) "hola....", H_hess(1,1)
c         write(*,*) "hola....", IH(1,1)
c         write(*,*) "hola....", HIH(1,1)


c         write(*,*) ier
c         write(*,*) istop

                   
         if(effet.eq.1.and.ier.eq.-1)then
           v((np-nva-indic_alpha)*(np-nva-indic_alpha+1)/2)=10.d10
         endif

c     vraisemblance non penalisee:            
c     res01(effet+1)=res+pe
c     vraisemblance penalisee:            
       res01(effet+1)=res
       
c       write(4,*)'valeur de ni',ni 
c       write(*,*)'valeur de ni',ni 
       
       
c     do 77 i=1,1000
c     write(*,*) 'derivee',v(i)
c     77      continue
       
c       if (effet.eq.1)then
c          if(AG.eq.1)then
c             write(4,*)'*************************************** ' 
c             write(4,*)'**** ANDERSEN-GILL APPROACH *********** ' 
c             write(*,*)'**** ANDERSEN-GILL APPROACH *********** ' 
c             write(4,*)'*************************************** ' 
c          endif
          
c          write(4,*)'**************** ' 
c          write(4,*)'THETA = variance de Z dans Z.exp(bX)'
c          write(4,*)'(Z suit une GAMMA)'
c          write(4,*)b(np-nva-indic_alpha)*b(np-nva-indic_alpha)
c          write(6,*)'**************** ' 
c       write(6,*)'THETA = ',b(np-nva-indic_alpha)*b(np-nva-indic_alpha)
c     pour tenir compte du changt de var : delta methode
c          write(6,*)'SE theta (=H)'
c     &         ,dsqrt(((2.d0*b(np-nva-indic_alpha))**2)*
c     &         H_hess(np-nva-indic_alpha,np-nva-indic_alpha))
c          write(4,*)'SE theta(=H)'
c     &         ,dsqrt(((2.d0*b(np-nva-indic_alpha))**2)*
c     &         H_hess(np-nva-indic_alpha,np-nva-indic_alpha))
c          write(6,*)'SE theta (=HIH)'
c     &         ,dsqrt(((2.d0*b(np-nva-indic_alpha))**2)*
c     &         HIH(np-nva-indic_alpha,np-nva-indic_alpha))
c          write(4,*)'SE theta(=HIH)'
c     &         ,dsqrt(((2.d0*b(np-nva-indic_alpha))**2)*
c     &         HIH(np-nva-indic_alpha,np-nva-indic_alpha))
c          write(4,*)'**************** ' 
c          write(6,*)'**************** ' 
          
c          write(4,*)'**************** ' 
c          write(4,*)'ALPHA dans z_i ^ alpha, pour le deces'
c          write(4,*)b(np-nva)
c          write(6,*)'**************** ' 
c          write(6,*)'ALPHA = ',b(np-nva)
c     pour tenir compte du changt de var : delta methode
c          write(6,*)'SE ALPHA (=H)',
c     &         dsqrt(H_hess(np-nva,np-nva))
c          write(4,*)'SE ALPHA (=H)',
c     &         dsqrt(H_hess(np-nva,np-nva))
c          write(6,*)'SE alpha (=HIH)',
c     &         dsqrt(HIH(np-nva,np-nva))
c          write(4,*)'SE alpha (=HIH)',
c     &         dsqrt(HIH(np-nva,np-nva))
          
c          write(4,*)'**************** ' 
c          write(6,*)'**************** ' 
          
c       endif
c       if(nva.gt.0)then
c         do 124 i=1,nva
c            j=(np-nva+i)*(np-nva+i+1)/2
c            bi = b(np-nva+i) - 1.96*dsqrt(H_hess(np-nva+i,np-nva+i))
c            bs = b(np-nva+i) + 1.96*dsqrt(H_hess(np-nva+i,np-nva+i))
c            if(i.eq.1)then
c               write(4,*)'*** FOR RECURRENT EVENTS ***'
c            endif
c            if(i.eq.nva1+1)then
c               write(4,*)'*** FOR DEATH ***',nva,i
c            endif
c            write(4,*)'**************** '
c            if(i.gt.nva1)then
c               write(4,*)'Variable : ',nomvar2(i-nva1)
c            else
c               write(4,*)'Variable : ',nomvar(i)
c            endif
c            write(4,*)i,')','beta=',b(np-nva+i)
c            write(4,*)' '
c            write(4,*)i,')',' SE (=H)',dsqrt(H_hess(np-nva+i,np-nva+i))
c            wres=(b(np-nva+i))/dsqrt(H_hess(np-nva+i,np-nva+i))
c            write(4,*)'---> WALD',wres
c            write(4,*)' '
c            write(4,*)i,')',' SE (=HIH)',dsqrt(HIH(np-nva+i,np-nva+i))
c       write(4,*)'---> WALD',(b(np-nva+i))/dsqrt(HIH(np-nva+i,np-nva+i))
c            write(4,*)' '
c            write(4,*)'RR : ',dexp(b(np-nva+i)),'  IC',dexp(bi),dexp(bs) 
c            write(4,*)'**************** ' 
           
c 124     continue 
c         write(4,*)'**************** ' 
c         write(4,*)'---> log vraisemb marginale complete pénalisée',res
c         write(4,*)'**************** ' 
c      endif
      
c             if(effet.eq.0)then
c               call distance1(nz1,nz2,b,fich1,fich2,fich1b,fich2b,effet)
c             else


c         call distance1(nz1,nz2,b,fich3,fich4,fich3b,fich4b,
c     & fich3c,fich4c,effet)


c             endif
          
c      call date_and_time(dateamj,heure2,zone,values)     
c      write(4,*) '***************************************************'   
c      write(4,*) '*** starting time: ***', dateamj,heure1
c      write(4,*) '*** Ending time:(hhmmss.sss) ***', dateamj,heure2
      



c --------------  Lambda and survival estimates JRG January 05

             call distance1(nz1,nz2,b,effet,
     &	       x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)      


      resOut=res


      do ss=1,npmax
       do sss=1,npmax
         HIHOut(ss,sss) = HIH(ss,sss)
         H_hessOut(ss,sss)= H_hess(ss,sss)
       end do  
      end do


      return

      end subroutine frailpenalJoint
      
  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCC**********SUBROUTINES******  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



c========================== VECSPLI ==============================
      subroutine vecspli1(n,ndate) 
      
      integer   n,ndate,i,j,k
      double precision  ht,htm,h2t,ht2,ht3,hht,h,hh,h2
      double precision  h3,h4,h3m,h2n,hn,hh3,hh2
         
c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************

c*****dace1 
      double precision date(ndatemax),datedc(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,datedc,zi
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
      subroutine vecpen1(n) 
      
c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************

      integer   n,i
      double precision  h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2
      double precision  a3,a2,b2,c2,a1,b1,c1,a0,x3,x2,x

c*****dace1 
      double precision date(ndatemax),datedc(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,datedc,zi
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


c================================  SEARPAS joly    ==============================

      SUBROUTINE searpas1(VW,STEP,B,BH,M,DELTA,FIM,EPSV,k0)
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
       CALL valfpa1(VLW1,FI1,B,BH,M,DELTA,k0)
       CALL valfpa1(VLW2,FI2,B,BH,M,DELTA,k0)
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
          CALL valfpa1(VLW1,FI1,B,BH,M,DELTA,k0)   
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
          CALL valfpa1(VLW1,FI1,B,BH,M,DELTA,k0)
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
         CALL valfpa1(VM,FIM,B,BH,M, DELTA,k0)
         IF (FIM.LE.FI2) GO TO 100
         VM=VLW2
         FIM=FI2
100   CONTINUE
      VW=DEXP(VM)
      RETURN

      END


   
c===================================   VALFPA joly   ==============================

        subroutine valfpa1(vw,fi,b,bk,m,delta,k0)
        integer  m,i
        double precision  vw,fi
        double precision  funcpa1,z        
        double precision k0(2)
        double precision b(M),bk(M),delta(M)


         z=0.d0
         do 1 i=1,m
            bk(i)=b(i)+dexp(vw)*delta(i)
1        continue
c
         fi=-funcpa1(bk,m,1,z,1,z,k0)
c
         return
         end   



c========================          FUNCPA NEW         ====================
      double precision function funcpa1(b,np,id,thi,jd,thj,k0)

c *** NOUVELLLE DECLARATION F90 :

c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************

      integer  nb,n,np,id,jd,i,j,k,vj,cptg,l,ig
      integer cpt(ngmax),choix

      double precision  thi,thj,pe1,pe2,dnb,sum
      double precision  inv
      double precision  som1,som2
      double precision  res,vet,vet2,h1
      double precision the1(-2:npmax),the2(-2:npmax)
      double precision b(np),bh(np)
      double precision res2(ngmax)
      double precision res1dc(ngmax),res2dc(ngmax),res3dc(ngmax)
      double precision k0(2) 
      double precision dut1(ndatemax),dut2(ndatemax)
      double precision ut1(0:ndatemax),ut2(0:ndatemax)
      double precision integrale1(ngmax),integrale2(ngmax)
      double precision integrale3(ngmax),int
      double precision gamma,nan

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
      double precision date(ndatemax),datedc(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,datedc,zi

c*****dace2
      double precision t0dc(ngmax),t1dc(ngmax)
      double precision t0(nsujetmax),t1(nsujetmax)
      integer c(nsujetmax), cdc(ngmax)
      integer nt0(nsujetmax),nt1(nsujetmax)
      integer nt0dc(ngmax),nt1dc(ngmax)
      integer  nsujet,nva,nva1,nva2,ndate,ndatedc,nst
      common /dace2/t0,t1,t0dc,t1dc,c,cdc,nt0,nt1,nt0dc,nt1dc,nsujet
     &     ,nva,nva1,nva2,ndate,ndatedc,nst
c*****dace4
      integer stra(nsujetmax)
      common /dace4/stra
c*****ve1
      double precision ve(nsujetmax,nvarmax)
      double precision vedc(ngmax,nvarmax)
      common /ve1/ve,vedc
c*****dace3
      double precision  pe
      integer  effet,nz1,nz2
      common /dace3/pe,effet,nz1,nz2
c*****contrib
      integer ng          !nb de gpes
      common /contrib/ng    
c*****groupe
      integer g(nsujetmax) 
      integer nig(ngmax) ! nb d events recurrents 
      common /gpe/g,nig
c %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
      integer AG
      common /andersengill/AG
c %%%%%%%%%%%%% indic ALPHA %%%%%%%%%%%%%%%%%%%%%%%%% 
      integer indic_ALPHA
      common /alpha/indic_ALPHA
c****  theta/alpha
      double precision  theta,alpha !en exposant pour la frailty deces 
      common /thetaalpha/ALPHA,theta
c*****auxig
      integer :: auxig
      common /auxig/auxig
c******  aux1 aux2
      double precision , dimension(ngmax) ::  aux1,aux2
      double precision , dimension(ngmax) :: res1,res3
      common /risqcumul/aux1,aux2,res1,res3
c****** indicateur de troncature
      integer :: indictronq,indictronqdc ! =0 si donnees non tronquées reellement
      common /troncature/indictronq,indictronqdc
c************************************************************  

     
c     write(*,*)'%% DANS FUNCPA%%',indictronq,AG,nva,nva1,nva2,g(1)

         do 3 i=1,np
            bh(i)=b(i)
 3       continue 
         if (id.ne.0) bh(id)=bh(id)+thi
         if (jd.ne.0) bh(jd)=bh(jd)+thj    
         
         n = (np-nva-effet-indic_ALPHA)/nst

          
         do 4 i=1,n
            the1(i-3)=(bh(i))*(bh(i))
            j = n+i 
c           write(*,*)'---the1',the1(i-3),i
             if (nst.eq.2) then
                the2(i-3)=(bh(j))*(bh(j))
c                write(*,*)'the2',the2(i-3),i
             endif
 4        continue


        if(effet.eq.1) then
         theta = bh(np-nva-indic_ALPHA)*bh(np-nva-indic_ALPHA)
c         write(*,*)'theta',theta
         alpha = bh(np-nva)
c         write(*,*)'alpha',alpha !pas de chgt de variable sur alpha
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
                     som1 = som1 + the1(j-4)
                     som2 = som2 + the2(j-4)
                     vj  = j
                  endif   
               endif
 6          continue 
            ut1(i) = som1 +(the1(j-3)*im3(i))+(the1(j-2)*im2(i))
     &       +(the1(j-1)*im1(i))+(the1(j)*im(i))
            dut1(i) = (the1(j-3)*mm3(i))+(the1(j-2)*mm2(i))
     &       +(the1(j-1)*mm1(i))+(the1(j)*mm(i))

            if(nst.eq.2)then
            ut2(i) = som2 +(the2(j-3)*im3(i))+(the2(j-2)*im2(i))
     &       +(the2(j-1)*im1(i))+(the2(j)*im(i))
            dut2(i) = (the2(j-3)*mm3(i))+(the2(j-2)*mm2(i))
     &       +(the2(j-1)*mm1(i))+(the2(j)*mm(i))
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
            res1dc(k) = 0.d0
            res2dc(k) = 0.d0
            res3dc(k) = 0.d0
            cpt(k) = 0
            integrale1(k) = 0.d0
            integrale2(k) = 0.d0
            integrale3(k) = 0.d0
            aux1(k)=0.d0
            aux2(k)=0.d0
 89      continue

c*******************************************         
C-----avec un effet aleatoire dans le modele
c*********************************************

         inv = 1.d0/theta

cccccccccccccccccccccccccccccccccccccccccc
c     pour les donnees recurrentes
cccccccccccccccccccccccccccccccccccccccccc
         do 10 i=1,nsujet 
            cpt(g(i))=cpt(g(i))+1  
            if(nva1.gt.0)then
               vet = 0.d0   
               do 19 j=1,nva1
                  vet =vet + bh(np-nva+j)*dble(ve(i,j))
c      write(*,*)'*** funcpa vet',vet,ve(i,j),i,j
 19            continue
               vet = dexp(vet)
            else
               vet=1.d0
            endif
            
            if((c(i).eq.1))then
               res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet) 
c     write(*,*)'***res2',res2(g(i)),dut1(nt1(i)),nt1(i),vet,i,g(i)
            endif  
c     nouvelle version
            res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet            
c     modification pour nouvelle vraisemblance / troncature:
              res3(g(i)) = res3(g(i)) + ut1(nt0(i))*vet 
c     write(*,*)'**res123',res1(g(i)),res2(g(i)),res3(g(i)),g(i),i
 10        continue

cccccccccccccccccccccccccccccccccccccccccc
c pour le deces 
cccccccccccccccccccccccccccccccccccccccccc 
         do 1000 k=1,ng 
             if(nva2.gt.0)then
                vet2 = 0.d0   
                do 191 j=1,nva2
                   vet2 =vet2 + bh(np-nva2+j)*dble(vedc(k,j))
c      write(*,*)'*** vet2',vet2,vedc(k,j),k,j
 191            continue
                vet2 = dexp(vet2)
             else
                vet2=1.d0
             endif
             if(cdc(k).eq.1)then
                res2dc(k) = dlog(dut2(nt1dc(k))*vet2)
c     write(*,*)'*** res2dc',res2dc(k),dut2(nt1dc(k)),vet2,k,ng
             endif 
             
c pour le calcul des integrales / pour la survie, pas les données recurrentes:
             aux1(k)=ut2(nt1dc(k))*vet2
             aux2(g(i))=aux2(g(i))+ut2(nt0(i))*vet2 !vraie troncature
 1000     continue
c          write(*,*)'*** aux1',aux1(ng)
c          stop
  
c**************INTEGRALES ****************************
          do 101   ig=1,ng 
             auxig=ig 
c              choix=1
c              call gaulag1(int,choix)
c              integrale1(ig) = int
c        write(*,*)'integrale1',integrale1(ig),frail,theta,
c     &             alpha,aux1(auxig),ig
              
c              if(indictronq.eq.1.and.AG.eq.0)then
c                 choix = 2
c                 call gaulag1(int,choix)
c                 integrale2(ig) = int
c                 write(*,*)'integrale2',integrale2(ig),ig
c              else 
c                 integrale2(ig) = 1.d0
c              endif

c              if(AG.eq.1)then
                 choix = 3  

                 call gaulag1(int,choix)

c                 integrale3(ig) = dexp(gamma(1.d0/theta +nig(ig)))*
c     &                (1.d0/theta + res1(ig)- res3(ig))**
c     &                (-nig(ig) - 1.d0 / theta)

                 integrale3(ig) = int !moins bon
c                 write(*,*)'integrale3',integrale3(ig),ig
c             endif
 101       continue
c                 write(*,*)'integrale3',integrale3(ng)
c      stop
c************* FIN INTEGRALES **************************
                      
          res = 0.d0 
          do 15 k=1,ng  
             sum=0.d0
             if(cpt(k).gt.0)then
                if(theta.gt.(1.d-8)) then
ccccc ancienne vraisemblance : pour calendar sans vrai troncature cccccccc
                   
                   res= res + res2(k)
c--      pour le deces:
     &                  + res2dc(k) 
     &                  - gamma(1./theta)-dlog(theta)/theta 
     &                  + dlog(integrale3(k))
c                   if(res.eq.nan)then 
c      write(*,*)'--func',res2(k),k!,integrale3(k),res,k
c                      stop
c                   endif 
                else
c*************************************************************************
c     developpement de taylor d ordre 3
c*************************************************************************
c                   write(*,*)'************** TAYLOR *************'                   
                   res= res + res2(k)
     &                  + res2dc(k) 
     &                  - gamma(1./theta)-dlog(theta)/theta 
     &                  + dlog(integrale3(k)) 
                endif
             endif 
 15       continue
               
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
            if(nst.eq.1)then
               pe2=0.d0
            else
            pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2)
     &    *the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+(
     &    the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)*
     &    m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0*
     &    the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)*
     &    m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1)
     &    *the2(i)*m1m(i))
            endif
 20       continue
          pe = k0(1)*pe1 + k0(2)*pe2 
c          write(*,*)'avant  penalisation ',res,pe,pe1,pe2
          res = res - pe
c           write(*,*)'vraisemblance penalisee',res,pe,pe1,pe2
c           write(*,*)'vraisemblance penalisee',k0(1),k0(2)
          funcpa1 = res 
c          write(*,*)''
c      write(*,*)'%%%%%%%%%% SORTIE DANS FUNCPA %%%%%%%%%%' 
c     &         res,integrale3(ng)
c 

c               stop
          return
          end


c================================  DCHOLE  ===========================
      subroutine dchole1(a,k,nq,idpos)

c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************  
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


c===============================    MARQ98  HJ =========================

   	 subroutine marq981(k0,b,m,ni,v,rl,ier,istop,effet)
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
c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************  
      
         integer  m,ni,nql,i,ii,nfmax,idpos,ier,istop,igrad,j
         integer  ncount,id,jd,i0,kkk
         integer effet,m1,k

         double precision  da,dm,ga,tr,g0
         double precision  ca,cb,epsa,epsb,rl
         double precision  funcpa1,det, step,eps,epsd
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

c*****dace2
      double precision t0dc(ngmax),t1dc(ngmax)
      double precision t0(nsujetmax),t1(nsujetmax)
      integer c(nsujetmax), cdc(ngmax)
      integer nt0(nsujetmax),nt1(nsujetmax)
      integer nt0dc(ngmax),nt1dc(ngmax)
      integer  nsujet,nva,nva1,nva2,ndate,ndatedc,nst
      common /dace2/t0,t1,t0dc,t1dc,c,cdc,nt0,nt1,nt0dc,nt1dc,nsujet
     &     ,nva,nva1,nva2,ndate,ndatedc,nst
      
c*****dace7
      double precision I_hess(npmax,npmax),H_hess(npmax,npmax)
      double precision Hspl_hess(npmax,npmax)
      double precision PEN_deri(npmax,1) 
      double precision hess(npmax,npmax)
      common /dace7/PEN_deri,I_hess,H_hess,Hspl_hess,hess
c %%%%%%%%%%%%% indic ALPHA %%%%%%%%%%%%%%%%%%%%%%%%% 
      integer indic_ALPHA
      common /alpha/indic_ALPHA ! pour preciser un para en plus 


c***********************************************
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
c     nql=nbre de seconds membres pour chole
      nql=1
 10   continue
c     
c     write(6,*)
c      write(*,*)'** Iteration :', ni 
c      if (ni.eq.5) stop
c      write(*,*)'** Parametres:'
c      write(*,*)(b(i),i=1,m)
c      write(6,*)'da',da,' ga=',ga
      
      z=0.d0
      i0=0 
      rl=funcpa1(b,m,i0,z,i0,z,k0)
      rl1=rl 

c     write(*,*)'log Vrais rl1,nsujet',rl1,nsujet,m
c      stop
c      if(ni.eq.0) write(8,*)'vrais init=',rl1
     
c      if (RL1.gt.0)stop
     
      call deriva1(b,m,v,rl,k0)
       
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
          call dchole1(fu,m,nql,idpos)
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
c              write(6,*) '** ** avant func.....',b1(i),delta(i),i
2             continue
              rl=funcpa1(b1,m,id,z,jd,z,k0)
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
           call searpas1(vw,step,b,bh,m,delta,fi,eps,k0) 
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
           if (ni.gt.350) then
              istop=2
c              write(6,*) 'nombre iteration max atteinte'
c              write(*,*) 'non convergence'
              goto 110
           end if
           goto 10
c********************
c     inversion matrice d'information
c********************
100      continue
           istop=1
          
c================ pour les bandes de confiance
c==== on ne retient que les para des splines
           call deriva1(b,m,v,rl,k0)
           do 303 i=1,(m*(m+3)/2)
              v1(i)=0.d0
 303       continue
           
           m1=m-nva-effet-indic_alpha

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
c             write(*,*)   'echec inversion matrice information'
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
c 103          format(1x,'echec inversion matrice information')
             istop=3
             call dsinv(v1,npl,ep,ier)
               if (ier.eq.-1) then
c             write(*,*)'echec inversion matrice information
c     & prms fixes'
                istop=31
               else
                 DO k=1,npl*(npl+1)/2
                  v(k)=v1(k)
                 END DO
               end if
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
         call deriva1(b,m,vnonpen,rl,zero)
         
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

      subroutine deriva1(b,m,v,rl,k0)

c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************  
       
      integer i0,iun,m,m1,ll,i,k,j
      double precision  funcpa1
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
c      print*,'***** entree deriva',m
c      stop
      th=1.d-5
      thn=-th
      th2=th*th
      z=0.d0
      i0=0
      iun =1
      rl=funcpa1(b,m,iun,z,iun,z,k0)
      do 2 i=1,m
	    fcith(i)=funcpa1(b,m,i,th,i0,z,k0)
c            print*,'fcith',fcith(i),i
2     continue
c      stop
      k=0
      m1=m*(m+1)/2
      ll=m1
      do 1 i=1,m
      ll=ll+1
      vaux=funcpa1(b,m,i,thn,i0,z,k0)
c      write(*,*)'vaux ok',vaux
      vl=(fcith(i)-vaux)/(2.d0*th)
c      write(*,*)'vl / deriva',vl,ll,i
c      write(*,*)'detail',fcith(i),vaux,(fcith(i)-vaux)
c      write(*,*)'vl',(fcith(i)-vaux)/(2.d0*th)
c      vl=(fcith(i)-funcpa1(b,m,i,thn,i0,z,k0))/(2.d0*th)
      v(ll)=vl
      do 1 j=1,i
      k=k+1
      v(k)=-(funcpa1(b,m,i,th,j,th,k0)-fcith(j)-fcith(i)+rl)/th2
1     continue
      
c      stop
      return
      end


c==========================  DISTANCE   =================================
    
   
         subroutine distance1(nz1,nz2,b,effet,
     &	       x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)

c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************  


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
      double precision date(ndatemax),datedc(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,datedc,zi
       
c*****dace2
      double precision t0dc(ngmax),t1dc(ngmax)
      double precision t0(nsujetmax),t1(nsujetmax)
      integer c(nsujetmax), cdc(ngmax)
      integer nt0(nsujetmax),nt1(nsujetmax)
      integer nt0dc(ngmax),nt1dc(ngmax)
      integer  nsujet,nva,nva1,nva2,ndate,ndatedc,nst
      common /dace2/t0,t1,t0dc,t1dc,c,cdc,nt0,nt1,nt0dc,nt1dc,nsujet
     &     ,nva,nva1,nva2,ndate,ndatedc,nst
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
                   hes1(i,j)=H_hess(i,j)
 115            continue
 116         continue 

         if(nst.eq.2)then  
            k = 0
            do 118 i=nz1+3,nz1+2+nz2+2
               k = k + 1 
               l = 0
               do 117 j=nz1+3,nz1+2+nz2+2
                  l = l + 1
                  hes2(k,l)=H_hess(i,j)
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
             call cosp1(x1,the1,nz1+2,hes1,zi,date,binf,su,bsup
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
               call cosp1(x2,the2,nz2+2,hes2,zi,date,binf,su,bsup
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
      subroutine susp1(x,the,n,su,lam,zi)
c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************  

      integer  j,k,n,i
      double precision  x,ht,ht2,h2,som,lam,su
      double precision  htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3
      double precision  ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2
      double precision  h,gl,hh
      double precision zi(-2:npmax),the(-2:npmax)


         som = 0.d0
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

      subroutine cosp1(x,the,n,y,zi,date,binf,su,bsup,lbinf,lam,lbsup)
c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************  
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

         call conf01(x,j,n,y,pm,zi,date)

         binf = dexp(-gl + 1.96d0*pm)
         su  = dexp(-gl)
         bsup = dexp(-gl - 1.96d0*pm)

         call conf11(x,j,n,y,pm,zi,date)
         lbinf = lam - 1.96d0*pm
         lbsup = lam + 1.96d0*pm
c         write(*,*)'lbinf apres conf1',lbinf,lam,pm

         return

         end
c=====================  CONF1  =============================
      subroutine  conf11(x,ni,n,y,pm,zi,date)
c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************  

      integer  ni,i,n,j

      double precision  mmsp1,x,pm
      double precision  res
      double precision date(ndatemax) 
      double precision zi(-2:npmax)
      double precision  vecti(npmax),aux(npmax)
      double precision y(npmax,npmax) 
      
           
            do 10 i=1,n
               vecti(i) = mmsp1(x,ni,i,zi,date)
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
      subroutine  conf01(x,ni,n,y,pm,zi,date)
c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************  
       
       integer  ni,i,n,j
       double precision isp1,x,pm
       double precision res
       double precision zi(-2:npmax)
       double precision date(ndatemax)
       double precision vecti(52),aux(52)
       double precision y(npmax,npmax)

            do 10 i=1,n
               vecti(i) = isp1(x,ni,i,zi,date)
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
          double precision function isp1(x,ni,ns,zi,date)
c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************  
          integer  ni,ns
          double precision  val,mmsp1,x
          double precision zi(-2:npmax)
          double precision date(ndatemax)


          if(x.eq.zi(ni))then
             if(ni.le.ns-3)then
                val = 0.d0
             else
                if(ni.le.ns-2)then
           val = ((zi(ni)-zi(ni-1))*mmsp1(x,ni,ns,zi,date))*0.25d0
                else
                   if (ni.eq.ns-1)then
                  val = ((zi(ni)-zi(ni-2))*mmsp1(x,ni,ns,zi,date)+
     &        (zi(ni+3)-zi(ni-1))*mmsp1(x,ni,ns+1,zi,date))*0.25d0
                   else
                      if(ni.eq.ns)then
                  val = ((zi(ni)-zi(ni-3))*mmsp1(x,ni,ns,zi,date)+
     &                (zi(ni+2)-zi(ni-2))*mmsp1(x,ni,ns+1,zi,date)
     &       +(zi(ni+3)-zi(ni-1))*mmsp1(x,ni,ns+2,zi,date))*0.25d0
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
                 val = (x-zi(ni))*mmsp1(x,ni,ns,zi,date)*0.25d0
             else  
             if(ni.eq.ns-2)then
                   val = ((x-zi(ni-1))*mmsp1(x,ni,ns,zi,date)+
     &       (zi(ni+4)-zi(ni))*mmsp1(x,ni,ns+1,zi,date))*0.25d0
             else   
                if (ni.eq.ns-1)then
                   val =((x-zi(ni-2))*mmsp1(x,ni,ns,zi,date)+
     &             (zi(ni+3)-zi(ni-1))*mmsp1(x,ni,ns+1,zi,date)
     &       +(zi(ni+4)-zi(ni))*mmsp1(x,ni,ns+2,zi,date))*0.25d0
                else
                   if(ni.eq.ns)then
                      val =((x-zi(ni-3))*mmsp1(x,ni,ns,zi,date)+
     &             (zi(ni+2)-zi(ni-2))*mmsp1(x,ni,ns+1,zi,date)
     &             +(zi(ni+3)-zi(ni-1))*mmsp1(x,ni,ns+2,zi,date)
     &        +(zi(ni+4)-zi(ni))*mmsp1(x,ni,ns+3,zi,date))*0.25d0
                   else
                      val = 1.d0
                   endif
                endif
             endif
             endif
          endif 
          endif
             isp1 = val
             return
             end
c==========================  MMSP   ==================================
      double precision function mmsp1(x,ni,ns,zi,date)
c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************ 
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

             mmsp1 = val
             return
             end

c============================    GAMMA      ==============================

c       function qui calcule le log de  Gamma
         double precision function gamma(xx)
            double precision xx
            integer j
            double precision cof(6),stp,half,one,fpf,x,tmp,ser
            data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     &       -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
            data half,one,fpf/0.5d0,1.0d0,5.5d0/

            x = xx - one
            tmp = x + fpf
            tmp = (x+half)*dlog(tmp) - tmp
            ser = one
            do 11 j = 1,6
               x = x + one
               ser = ser + cof(j)/x
 11         continue
            gamma = tmp + dlog(stp*ser)
            return
            end
   
c======================================================

c==================================================================
c==================================================================
c 20 points / sur (0,+infty)
      SUBROUTINE gaulag1(ss,choix) 

c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************

      double precision :: ss
      double precision ::auxfunca
      double precision :: func1,func2,func3
      external :: func1,func2,func3
c gauss laguerre
c func1 est l intégrant, ss le résultat de l integrale sur 0 ,  +infty
      
      INTEGER :: j,choix,nan

c*****auxig
      integer :: auxig
      common /auxig/auxig 
c*****************************
      DOUBLE PRECISION :: x(20),w(20)!The abscissas-weights.
      SAVE  w,x
c      DATA w/0.354009738607,0.831902301044,1.33028856175,1.86306390311,
c     & 2.45025555808,3.12276415514,3.9341526956,4.99241487226,
c     & 6.57220248513,9.78469584034/
c      DATA x/0.13779347054,0.729454549503,1.80834290174,
c     & 3.40143369785,5.55249614006,8.33015274676,11.8437858379,
c     & 16.2792578314,21.996585812,29.9206970123/ 
c/*** 20 points ***/
      DATA w/0.181080062419,0.422556767879,0.666909546702,0.91535237279,
     &  1.1695397071,1.43135498624,1.7029811359,1.98701589585,
     &  2.28663576323,2.60583465152,2.94978381794,3.32539569477,
     &  3.74225636246,4.21424053477,4.76252016007,5.42172779036,
     &  6.25401146407,7.38731523837,9.15132879607,12.8933886244/
    
      DATA x/0.070539889692,0.372126818002,0.916582102483,1.70730653103,
     &  2.74919925531,4.04892531384,5.61517497087,7.45901745389,
     &  9.59439286749,12.0388025566,14.8142934155,17.9488955686,
     &  21.4787881904,25.4517028094,29.9325546634,35.0134341868,
     &  40.8330570974,47.6199940299,55.8107957541,66.5244165252/
     
c****************************************
      ss=0.d0 
c Will be twice the average value of the function,since the ten
c wei hts (five numbers above each used twice) sum to 2.
      do 11 j=1,20
         if (choix.eq.1) then !integrale 1
            auxfunca=func1(x(j))
            ss = ss+w(j)*(auxfunca)
            if(ss.eq.nan)then 
c               write(*,*)'----int1',ss ,w(j),auxfunca,j
c               stop
            endif

         else                   !choix=2, survie marginale, vraie troncature
            if (choix.eq.2) then 
               auxfunca=func2(x(j))
               ss = ss+w(j)*(auxfunca)
               else                   !choix=3, AG model
                  if (choix.eq.3) then
                     auxfunca=func3(x(j))
                     ss = ss+w(j)*(auxfunca)
                  endif
               endif
            endif
 11         continue
c      print*,'** ss',ss
c      stop
      return
      END
  
c================================================

      double precision  function func1(frail) 
! calcul de l integrant, pour un effet aleatoire donné frail et un groupe donne auxig (cf funcpa)      
      
c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************

      double precision  frail
c      double precision func1
      integer :: ig,issg,k,i,nan
c*****auxig
      integer :: auxig
      common /auxig/auxig
c****  theta/alpha
      double precision  theta,alpha !en exposant pour la frailty deces 
      common /thetaalpha/alpha,theta
c*****mi
      integer , dimension(ngmax) :: mi ! nb de dc dans gpe i/pas donnees recurrentes
      common /mi/mi
c******  aux1 aux2
      double precision , dimension(ngmax) ::  aux1,aux2
      double precision , dimension(ngmax) :: res1,res3
      common /risqcumul/aux1,aux2,res1,res3
            
c============================================================
      
      func1 = (frail**(alpha*mi(auxig)+1./theta-1.))*
     &     dexp(-(frail**alpha) *aux1(auxig))*
     &     dexp(-frail/theta)
               if(func1.le.0.d0)then
c      print*,"** func1 **",func1,frail,theta,alpha,mi(auxig),aux1(auxig)
c      print*,"** **",(frail**(alpha*mi(auxig)+1./theta-1.))
c      print*,"** **",dexp(-(frail**alpha) *aux1(auxig))
c      print*,"** **",dexp(-frail/theta)
c** pb : quand frail grand (35.01) et theta petit 0.0422 --> exp(-836)=0
c** et func1 = 0
c      print*,"** **",dexp(-(frail**alpha) *aux1(auxig) - frail/theta)
c       stop
      endif
c            if(func1.eq.nan)then 
c              write(*,*)'----func1',func1,frail,theta,dexp(-frail/theta)
c               stop
c            endif
      return
      end
c==================================================================

c================================================

      double precision  function func2(frail) 
! calcul de l integrant, pour un effet aleatoire donné frail et un groupe donne auxig (cf funcpa)      
      
c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************

      double precision  frail,gamma
c      double precision func2
      integer :: ig,issg,k,i
c*****auxig
      integer :: auxig
      common /auxig/auxig
c****  theta/alpha
      double precision  theta,alpha !en exposant pour la frailty deces 
      common /thetaalpha/ALPHA,theta
c*****mi
      integer , dimension(ngmax) :: mi ! nb de dc dans gpe i/pas donnees recurrentes
      common /mi/mi
c******  aux1 aux2
      double precision , dimension(ngmax) ::  aux1,aux2
      double precision , dimension(ngmax) :: res1,res3
      common /risqcumul/aux1,aux2,res1,res3
            
c============================================================
      
      func2 = dexp(-(frail**alpha)*aux2(auxig))*dexp(-frail/theta)
     & *(frail)
     & /(exp(gamma(1.d0/theta))*(theta**(1./theta)))
               
c     print*,"** func2 **",func2,frail,eta,alpha,mi(auxig)
c     stop
      return
      end
c==================================================================

      double precision  function func3(frail) 
! calcul de l integrant, pour un effet aleatoire donné frail et un groupe donne auxig (cf funcpa)      
      
c************** definition commune des parameter ***********************
      integer , parameter :: npmax=50,nsujetmax=15000,nvarmax=50
      integer , parameter :: ngmax=1000
      integer , parameter :: ndatemax=30000
c************************************************************************ 

      double precision  frail
c      double precision func3
      double precision gamma
      integer :: ig,issg,k,i ,nan    
c*****groupe
      integer g(nsujetmax)
      integer nig(ngmax)  ! nb d events recurrents par sujet
      common /gpe/g,nig
c*****auxig
      integer :: auxig
      common /auxig/auxig
c****  theta/alpha
      double precision  theta,alpha !en exposant pour la frailty deces 
      common /thetaalpha/alpha,theta
c******  aux1 aux2
      double precision , dimension(ngmax) ::  aux1,aux2
      double precision , dimension(ngmax) :: res1,res3
      common /risqcumul/aux1,aux2,res1,res3
c*****dace2
      double precision t0dc(ngmax),t1dc(ngmax)
      double precision t0(nsujetmax),t1(nsujetmax)
      integer c(nsujetmax), cdc(ngmax)
      integer nt0(nsujetmax),nt1(nsujetmax)
      integer nt0dc(ngmax),nt1dc(ngmax)
      integer  nsujet,nva,nva1,nva2,ndate,ndatedc,nst
      common /dace2/t0,t1,t0dc,t1dc,c,cdc,nt0,nt1,nt0dc,nt1dc,nsujet
     &     ,nva,nva1,nva2,ndate,ndatedc,nst    
            
c %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
      integer AG
      common /andersengill/AG
c============================================================
      


      func3 = (nig(auxig)+ alpha*cdc(auxig)+ 1./theta-1.)*dlog(frail)
     &     - frail*(res1(auxig)-res3(auxig)) !res3=0 si AG=0
     &     - (frail**alpha)*(aux1(auxig))
     &     - frail/theta
     
      func3 = exp(func3)

      return
      end
c==================================================================
