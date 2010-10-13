c gaulag NEW avec 20 points / modifie 25 fevrier
c
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
c=====================================
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
c  lorsque la matrice des derivees 2 n'est pas inversible,
c  il faut reinitialiser completement la matrice FU et
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

c    estimation des coefficient de regression  ... spline d'ordre 4
c    avec censure a droite, troncature a gauche, noeuds aux choix


c     avec fraitly avec loi Gamma

       subroutine nested(nsujetAux,ngexactAux,nssgbygAux,icenAux,
     &     nstAux,nzAux,ax1,ax2,tt0Aux,tt1Aux,icAux,groupeAux,
     &     ssgroupeAux,nvaAux,strAux,vaxAux,AGAux,noVar,
     &     maxitAux,irep1,np,b,H_hessOut,HIHOut,resOut,x1Out,
     &     lamOut,suOut,x2Out,lam2Out,su2Out,ni,cpt,ier,k0,ddl)


      

c PARAMETER
c***************SAME **********************************************
      integer npmax,nsujetmax,nvarmax,ngmax,nssgmax,
     &   nboumax,ndatemax,nptsmax 
      parameter(npmax=50,nsujetmax=20000,nvarmax=50)
      parameter(ngmax=2500,nssgmax=5000,nboumax=2000)
      parameter(ndatemax=30000,nptsmax=1500)
c******************************************************************

     
c  Added JRG August '06
      integer nsujetAux,ngexactAux,nssgbygAux,icenAux,nstAux
      integer nzAux,nvaAux,AGaux,noVar,maxitAux 
      double precision ax1,ax2
      double precision tt0Aux(nsujetAux),tt1Aux(nsujetAux)
      integer icAux(nsujetAux),groupeAux(nsujetAux)
      integer ssgroupeAux(nsujetAux)
      double precision strAux(nsujetAux)
      double precision vaxAux(nsujetAux,nvaAux)
      double precision resOut,H_hessOut(npmax,npmax)
      double precision HIHOut(npmax,npmax)

      double precision  x1Out(99),lamOut(99,3),suOut(99,3),
     &                  x2Out(99),lam2Out(99,3),su2Out(99,3)

      

      integer :: isujet
      integer :: groupe,ssgroupe
      integer :: j,k,nz,n,np,cpt,ii,iii,ver
      integer :: cptstr1,cptstr2
      integer :: i,ic,ni,ier,istop,ef,l
      integer :: cptni,cptni1,cptni2
      integer :: nb_echec,nb_echecor
      integer :: nis,m,idum,icen,id,ibou,nbou,nbou2
      integer :: auxng,auxssng  
      integer :: ig,issg
c yas déb
      integer  stracross(nsujetmax) !pour crossvalidation
      integer irep1,nvacross,nstcross,effetcross !pour crossvalidation
      integer nvat,ss,sss,ngaux	  
c yas fin      
      integer  filtre(nvarmax)
      integer  cpt1(nboumax)
      integer  cpt2(nboumax)
      integer  cpt3(nboumax)
      integer ind(nboumax)
c yas      integer irep1,nvat,ss,sss
      
      real  vax(nvarmax),tt0,tt1

      double precision :: h
      double precision :: varxij,eca,varsmarg
      double precision  :: smoy
      double precision :: wres
      double precision :: xmin1,xmin2,res,min,max,maxt
      double precision :: bi,bs,wald
      double precision ::str
      double precision :: pe1,pe2,lrs,trace,trace1,trace2
      double precision :: uniran
      double precision :: moyvar1,moyvar2
      double precision :: moyse1,moyse2, moyse_cor1,moyse_cor2
      double precision auxi,ax,bx,cx,tol,ddl
      double precision fa,fb,fc,golden2,estimv2
      
      double precision aux(2*nsujetmax)
      double precision y(npmax,npmax)
      double precision v((npmax*(npmax+3)/2))
      double precision k0(2),res01(2)
      double precision b(npmax)
      double precision se1(nboumax),se_cor1(nboumax)
      double precision se2(nboumax),se_cor2(nboumax)
      double precision tvars1(nboumax),tvars2(nboumax)
      double precision I1_hess(npmax,npmax),H1_hess(npmax,npmax)
      double precision I2_hess(npmax,npmax),H2_hess(npmax,npmax)
      double precision HI1(npmax,npmax),HI2(npmax,npmax),
     &    HIH(npmax,npmax),IH(npmax,npmax),HI(npmax,npmax)


c *** Added Virginie Aug'07
      integer , dimension(nsujetmax) :: gaux,gnew
      double precision :: bgpe,bssgpe




c************ declarations pour donnees generees **********
         integer :: no,nbn,nb1,nb0,delta,grp,sgrp
         integer :: sg,nb
         real :: v1,v2
         real :: piece,rien,demi
         double precision :: ui,uij,temps1
         double precision :: x,tronc,cens,beta1,beta2,beta
         double precision :: bg1(2),bg2(2),bw1(2)
c*****************************************************************

c***** nmax
         integer :: nmax
         common /nmax/nmax
c*****dace1 
      double precision date(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,zi

c*****dace2
      double precision t0(nsujetmax),t1(nsujetmax)
      double precision ncsrl
      integer c(nsujetmax)
      integer nt0(nsujetmax),nt1(nsujetmax)
      integer nsujet,nsim,nva,ndate,nst
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,nst
c*****dace4
      integer stra(nsujetmax)
      common /dace4/stra
c*****ve1
      double precision ve(nsujetmax,npmax)
      common /ve1/ve
c*****dace3
      double precision pe
c yas      integer :: effet,nz1,nz2,effetcross
   	  integer effet,nz1,nz2
      common /dace3/pe,effet,nz1,nz2
c*****dace7
      double precision I_hess(npmax,npmax),H_hess(npmax,npmax)
      double precision Hspl_hess(npmax,npmax)
      double precision PEN_deri(npmax,1)
      double precision hess(npmax,npmax)
      common /dace7/PEN_deri,I_hess,H_hess,Hspl_hess,hess
c*****contrib
      integer ngexact,nssgexact    !nb EXACT de gpes et nb exact de ss gpes au total
      integer nssgbyg  ! nb de ss gpes par groupe MAXimum 
      common /contrib/ngexact,nssgexact,nssgbyg 
c*****groupe
      integer g(nsujetmax) 
      integer ssg(nsujetmax,ngmax)
      integer nig(ngmax)
      common /gpe/g,nig,ssg
c****** indicateur de troncature
      integer :: indictronq ! =0 si donnees non tronquées
      common /troncature/indictronq
c*****mij
c ! nb de dc dans gpe i
      integer mi(ngmax) 
c ! nb de dc dans ssgpe ij
      integer  mij(ngmax,nssgmax)
      common /mij/mij,mi
 
c %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
      integer AG
      common /andersengill/AG

c****** indicateur du nb de parametres
      integer  nbpara
      common /nbpara/nbpara

c************ FIN COMMON ***********************************
       
      character (len=10):: nomvarl
      character (len=10):: nomvar(25)
      character (len=14):: donnees
      character (len=14)::fich1
      character (len=14)::fich2
      character (len=14)::fich3
      character (len=14)::fich4
      character (len=14)::fich1b
      character (len=14)::fich2b
      character (len=14)::fich3b
      character (len=14)::fich4b
      character*20 dateamj
      character*20 zone
      character*20 heure1
      character*20 heure2
      integer  values(8)

        
c nst: nombre de strates
c ist: appartenance aux strates
c icen=1: censure
c ib: matrice de variance de beta chapeau
c I_hess : -hessienne non inversee sur vraisemblance non penalisee
c H_hess : inverse de -hessienne  sur vraisemblance penalisee

c nig(ngmax) : nb sujet dans groupe (different dans frailty.f !!!)
c mi(ngmax) : nb deces dans groupe  ! double precision
c mij(ngmax,nssgmax) : nb deces dans ss groupe
c g(nsujetmax) : numero du groupe pour un  sujet donné 
c ssg(nsujetmax,ngmax): numero du sous-groupe pour 1 sujet donnée d'1 groupe donné

c ngexact : nb exact de groupes
c nssgexact : nb exact de sous groupes au total
c nssgbyg : nb de sous-groupes par groupe au max



c Changed by Juan August '07
!nb iterations max dans marquard

       nmax = maxitAux 


 
c Juan     write(*,*)'  *************************************************'
c Juan      write(*,*)' ****** DEBUT PROGRAMME NESTED.F **********'
c Juan      write(*,*)'***********************************************'

c Juan      open(4,file='output')

c Juan      write(4,*)'  *************************************************'
c Juan      write(4,*)' ****** DEBUT PROGRAMME NESTED.F **********'
c Juan      write(4,*)'***********************************************'

c Juan      call date_and_time(dateamj,heure1,zone,values)      

c=== initialisations 
        lrs=0.d0    
        nb_echec=0 
        nb_echecor=0
c=== fin initialisations 
    


c Juan  	open(2,file='nested.inf')
c Juan  c 1900	continue (pour faire plusieurs jeux de simulations)
c Juan          read(2,*)nsujet
c Juan          read(2,*)ngexact
c Juan          read(2,*)nssgbyg
c Juan          read(2,*)icen
c Juan          read(2,*)nst
c Juan          read(2,*)ver
c Juan           nva = 0
c Juan           if(ver.gt.0)then 
c Juan              do 44 j=1,ver
c Juan                 read(2,*)nomvarl,filtre(j) 
c Juan  c               PRINT*,'filtre ',nomvarl,filtre(j)   
c Juan                 nva = nva + filtre(j)
c Juan                 if(filtre(j).eq.1)then
c Juan                    nomvar(nva) = nomvarl
c Juan                 endif   
c Juan   44         continue
c Juan           endif
         




c Juan  Added August'07

        nsujet=nsujetAux
        ngexact=ngexactAux
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
          enddo  
          nva=nvaAux  
        end if  

        ver=nvaAux
       
        
c Juan        nbou=1
c Juan        read(2,*)donnees
c Juan        open(9,file=donnees)


c %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
c Juan        read(2,*)AG
        AG=AGAux

c Juan        read(2,*)nz
        nz=nzAux  

c        nz1=nz
c        if(nz.gt.20)then
c            nz = 20
c        endif
c        if(nz.lt.4)then
c            nz = 4
c        endif

        cptni=0
        cptni1=0
        cptni2=0

        res01(1)=0.d0
        res01(2)=0.d0
       
               
        effet=1 
        do i=1,ngmax   
         nig(i) = 0
         do j=1,nsujetmax
           ssg(j,i)=0
           g(j)=0
         enddo
        enddo

       
c Juan        write(4,*)'nb de variable explicative dans le modele:',nva
c Juan        write(*,*)'nb de variable explicative dans le modele:',nva
    
c------------  lecture fichier -----------------------

         maxt = 0.d0
         cpt = 0
         k = 0
         cptstr1 = 0
         cptstr2 = 0

c Juan         write(*,*)'** Fichier de données = ',donnees
c Juan         write(4,*)'** Fichier de données = ',donnees
c Juan         open(9,file=donnees)
         
         do 10 i = 1,nsujet    
c            str=1
            if(nst.eq.2)then
c Juan               read(9,*)tt0,tt1,ic,ssgroupe,groupe,str,(vax(j),j=1,ver)
               tt0=tt0Aux(i)
               tt1=tt1Aux(i)
               ic=icAux(i)
               ssgroupe=ssgroupeAux(i)
               groupe=groupeAux(i)
               str=strAux(i)
               do j=1,ver
                 vax(j)=vaxAux(i,j)  
               enddo
            else
c Juan               read(9,*)tt0,tt1,ic,ssgroupe,groupe,(vax(j),j=1,ver) 
c              write(*,*)'**',tt0,tt1,ic,groupe,ssgroupe,(vax(j),j=1,ver) 

               tt0=tt0Aux(i)
               tt1=tt1Aux(i)
               ic=icAux(i)
               ssgroupe=ssgroupeAux(i)
               groupe=groupeAux(i)
               str=strAux(i)
               do j=1,ver
                 vax(j)=vaxAux(i,j)  
               enddo
            endif
            k = k +1

            if(k.eq.1)then
               auxng=groupe
               auxssng=ssgroupe
               ngexact=1
               nssgexact=1
               g(k)=1
               ssg(k,g(k))=1
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
                  t0(k) = dble(tt0)   
c                  t0(k) = 0.d0
                  t1(k) = dble(tt1)   
                
                  if(auxng.ne.groupe)then !chgt de groupes
                     ngexact=ngexact+1
                     auxng=groupe
                     nssgexact=nssgexact+1
                     auxssng=ssgroupe
                     if(k.ne.1)then
                        g(k)=g(k-1)+1
                        ssg(k,g(k))=1  !ssg(k-1,g(k-1))+1 
                     endif
c     si on suppose vraiment une structure NESTED
c        write(*,*)'** chgt gp',k,ssg(k,g(k)),nssgexact,g(k),ngexact,ic
                 goto 100             
                  endif
                  
                  if(auxssng.ne.ssgroupe.and.auxng.eq.groupe)then 
c     chgt de ssgroupe mais pas de groupe
                     nssgexact=nssgexact+1
                     auxssng=ssgroupe
                     if(k.ne.1)then
                        g(k)=g(k-1)
                        ssg(k,g(k))=ssg(k-1,g(k-1))+1
c        write(*,*)'** auxA ',k,g(k),ssg(k,g(k)),ssg(k-1,g(k-1))
                     endif
c      write(*,*)'*chgA ',k,ssg(k-1,g(k-1)),auxssng,ssgroupe,auxng,groupe
c        write(*,*)'** chgA ssgp',k,ssg(k,g(k)),nssgexact,g(k),ngexact,ic
                 goto 100             
                  endif
                  
                  if(k.ne.1)then
                     if(auxssng.eq.ssgroupe.and.auxng.eq.groupe)then 
c     aucun chgt de ssgroupe ni de gpe
                     if(k.ne.1)then
                        g(k)=g(k-1)
                        ssg(k,g(k))=ssg(k-1,g(k-1))
                     endif
c      write(*,*)'** pas chgt ',k,ssg(k,g(k)),nssgexact,g(k),ngexact,ic
                        goto 100             
                     endif
                  endif
 100              continue

                  nig(g(k)) = nig(g(k))+1
                  iii = 0
                  do 6 ii = 1,ver
                     if(filtre(ii).eq.1)then
                        iii = iii + 1
                        ve(i,iii) = dble(vax(ii))
                     endif
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
                     do 8 ii = 1,ver
                        if(filtre(ii).eq.1)then
                           iii = iii + 1
                           ve(i,iii) = dble(vax(ii))
                        endif
 8                    continue 
                      t0(k) =  dble(tt0)
c                      t0(k) = 0.d0
                      t1(k) = dble(tt1)
                      if(auxng.ne.groupe)then !chgt de groupes
                         ngexact=ngexact+1
                         auxng=groupe
                         nssgexact=nssgexact+1
                         auxssng=ssgroupe
                         if(k.ne.1)then
                            g(k)=g(k-1)+1
                            ssg(k,g(k))=1   !ssg(k-1,g(k-1))+1 
                         endif
c     si on suppose vraiment une structure NESTED
c      write(*,*)'** chgt gp',k,ssg(k,g(k)),nssgexact,g(k),ngexact,ic
                 goto 101
                      endif
                      
                      if(auxssng.ne.ssgroupe.and.auxng.eq.groupe)then 
c     chgt de ssgroupe mais pas de groupe
                         nssgexact=nssgexact+1
                         auxssng=ssgroupe
                         if(k.ne.1)then
                            g(k)=g(k-1)
                            ssg(k,g(k))=ssg(k-1,g(k-1))+1
                         endif
c      write(*,*)'** chgt ssgp',k,ssg(k,g(k)),nssgexact,g(k),ngexact,ic
               goto 101             
                      endif
                      
                      if(k.ne.1)then
                         if(auxssng.eq.ssgroupe.and.auxng.eq.groupe)then 
c     aucun chgt de ssgroupe ni de gpe
                         if(k.ne.1)then
                            g(k)=g(k-1)
                            ssg(k,g(k))=ssg(k-1,g(k-1))
                         endif
c      write(*,*)'** pas chgt',k,ssg(k,g(k)),nssgexact,g(k),ngexact,ic
                 goto 101             
                         endif
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
c             write(*,*)'** nig **',(nig(i),i=1,ngexact)
c          write(*,*)'**  groupe **',(g(i),i,i=1,nsujet)
c          write(*,*)'**  *********************************'
c             do 103 i=1,nsujet
c          write(*,*)'** sous groupe **',ssg(i,g(i)),i,g(i),ngexact
c 103   continue
c          stop

c         write(*,*)'max',maxtemps
c         write(*,*)'** exact nb of groups =',ngexact
c     write(*,*)'** num de groupe =',(g(i),i=1,nsujet)
c         write(*,*)'** exact nb of sub-groups =',nssgexact
c     write(*,*)'** num de ssgroupe =',(ssg(i,g(i)),i=1,nsujet)
         
c         write(4,*)'** exact nb of groups =',ngexact 
c         write(4,*)'** exact nb of sub-groups =',nssgexact 
c     write(4,*)'nb individu par groupe et par strate = ',nis
c         write(4,*)'** nb of strata = ',nst
c         write(4,*)'** nb of observations',nsujet
c         write(4,*)'** nb of failures ',cpt
c         write(4,*)'** nb of censoring ',(nsujet-cpt)
c         write(*,*)'** nb of observations',nsujet
c         write(*,*)'** nb of failures ',cpt
c         write(*,*)'** nb of censuring ',(nsujet-cpt)



            nz1=nz
            if(nst.eq.2)then
               nz2=nz
            endif
            if(nz.gt.20)then
               nz = 20
            endif
            if(nz.lt.4)then
               nz = 4
            endif
c      write(4,*)'nombre de noeuds :',nz
c     write(*,*)'nombre de strates:',nst
    

c***************************************************
c--------------- zi- ----------------------------------


c      construire vecteur zi (des noeuds)



         min = 1.d-10
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
c            write(*,*) '** aux **',aux(i)
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
c         stop
   
         zi(-2) = date(1)
         zi(-1) = date(1)
         zi(0) = date(1)
         zi(1) = date(1)
         h = (date(ndate)-date(1))/(nz-1)
         do 18 i=2,nz-1
            zi(i) =zi(i-1) + h
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
 50       continue   

c test sans troncature
c             indictronq=0

c---------- affectation des vecteurs de splines -----------------
             n  = nz+2
            
             call vecspli2(n,ndate)
             call vecpen2(n)
                    
             np = nst*n + nva + 2*effet
             nbpara =np

c             write(*,*)'nombre total de paramètres',np,nst,n,nva,effet
               
c------- initialisation des parametres
                   
             do 75 i=1,np
                b(i)=3.d-1
 75          continue
c-------------------------------------------------


         

        
c***********************************************************
c************** NEW : cross validation  ***********************
c***********************************************************

c yas             effet=0

c   Juan             read(2,*)irep1     !=0 pour une recherche du parametre de lissage


c  Juan: This must be considered in the R program
c             if(irep1.eq.0.and.nst.ge.2)then
c                write(*,*)'Error : you need only one stratum'
c                stop
c             endif


c   Juan             read(2,*)xmin1

c   Juan    irep1 e xmin1 se pasan en la subroutina

c   Juan             write(4,*)'STARTING VALUE FOR KAPPA',xmin1
c              write(*,*)'STARTING VALUE FOR KAPPA',xmin1

c yas déb
             nvacross=nva   !pour la recherche du parametre de lissage sans var expli
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

c         read(2,*)irep1  !=0 pour une recherche du parametre de lissage
         if(irep1.eq.0.and.nst.ge.2)then
c            write(*,*)'Error : you need only one stratum'
c            stop
         endif
c         read(2,*)xmin1
c         write(4,*)'STARTING VALUE FOR KAPPA',xmin1
c         write(*,*)'STARTING VALUE FOR KAPPA',xmin1
c yas fin

             if(xmin1.le.0.d0)then
                xmin1 = 0.d0
             endif  




c*************************************************

         nvat=nva
         nva=0
  
         if(irep1.eq.1)then   !recherche du parametre de lissage
            xmin1 = dsqrt(xmin1)
c            write(*,*) 'auxi',n,y,ddl,ni,res
            auxi = estimv2(xmin1,n,b,y,ddl,ni,res)

            if (ni.ge.250) then
c              write(*,*) ' '
c              write(*,*) 'no convergence 
c     &             with the chosen smoothing parameter'

              do 175 i=1,nz+2
                 b(i)=1.d-1
 175           continue     
              xmin1 = sqrt(10.d0)*xmin1
              auxi = estimv2(xmin1,n,b,y,ddl,ni,res)
              if (ni.lt.250) then
c                 write(*,*)' '
c            write(*,*)'Value of the smoothing parameter :'
c     &       ,real(xmin1*xmin1), '  DoF :',-ddl
c                write(*,*)' '
c                 write(*,*)'Log-vraisemblance :',res
c                 stop
              else
                 do 176 i=1,nz+2
                    b(i)=1.d-1
 176             continue     
                 xmin1 = sqrt(10.d0)*xmin1
                 auxi = estimv2(xmin1,n,b,y,ddl,ni,res)
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
c    &       ,real(xmin1*xmin1), '  DoF :',-ddl
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
            auxi = estimv2(xmin1,n,b,y,ddl,ni,res)
c            write(*,*)'estimv1',ddl
c               stop
            if(ddl.gt.-2.5d0)then
c            write(*,*)'OK non'
c               stop
               xmin1 = dsqrt(xmin1)
c               stop
               auxi = estimv2(xmin1,n,b,y,ddl,ni,res)
c               write(*,*)'estimv2'
c               stop !pb
               if(ddl.gt.-2.5d0)then
                  xmin1 = dsqrt(xmin1)
                  auxi = estimv2(xmin1,n,b,y,ddl,ni,res)
c                  write(*,*)'estimv3'
                  if(ddl.gt.-2.5d0)then
                     xmin1 = dsqrt(xmin1)
                     auxi = estimv2(xmin1,n,b,y,ddl,ni,res)
c                     write(*,*)'estimv4'
                     if(ddl.gt.-2.5d0)then
                        xmin1 = dsqrt(xmin1)
                        auxi = estimv2(xmin1,n,b,y,ddl,ni,res)
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
              auxi = estimv2(xmin1,n,b,y,ddl,ni,res)
c              write(*,*)'estimv5'
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
            
            call mnbrak2(ax,bx,cx,fa,fb,fc,b,n)
c            write(*,*)'OK 2' 
c            stop   
            
            tol = 0.001d0
c            write(*,*)'avant golden',ax,bx,cx,tol,xmin1,n,b,ddl
c            stop
            res = golden2(ax,bx,cx,tol,xmin1,n,b,y,ddl)
c            write(*,*)'apres golden'
c            stop
c            write(4,*)'******************************************* '
c            write(4,*)'Best smoothing parameter',real(xmin1*xmin1)
c     &                , '  DoF :',-ddl
c            write(4,*)'******************************************* ' 
c            write(*,*)'******************************************* '
c            write(*,*)'Best smoothing parameter',real(xmin1*xmin1)
c     &                , '  DoF :',-ddl
c            write(*,*)'******************************************* '
            effet=0
c            stop
c yas            call marq982((xmin1*xmin1),b,n,ni,v,res,ier,istop,
c yas     &           effet)
c yas déb
            call marq982((xmin1*xmin1),b,n,ni,v,res,ier,istop)	 
c yas fin	 
c            if (ni.ge.250) then
c               write(*,*)' '
c              write(*,*) 'non convergence'
c            else
c               write(*,*)' '
c               write(*,*)'Log-vraisemblance :',res
c            endif
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
c        write(*,*)'nva',nva
c        stop
ccccc********************************************************************



c   Juan         if(nst.eq.2)then
c   Juan            read(2,*)xmin2
c   Juan         endif


c   Juan         read(2,*)fich1
c   Juan         if(nst.eq.2)then
c   Juan            read(2,*)fich2
c   Juan         endif
c   Juan         read(2,*)fich3 
c   Juan         if(nst.eq.2)then
c   Juan            read(2,*)fich4
c   Juan         endif
c   Juan         read(2,*)fich1b
c   Juan         if(nst.eq.2)then
c   Juan            read(2,*)fich2b
c   Juan         endif
c   Juan         read(2,*)fich3b 
c   Juan         if(nst.eq.2)then
c   Juan            read(2,*)fich4b
c   Juan         endif
c         close(2)

                 k0(1) = xmin1*xmin1
                 if(nst.eq.2)then
                   k0(2) = xmin2
                 endif

c         write(*,*),'lissage',k0(1)
c         write(*,*),'lissage2',k0(2)




C=============================- fin cross validation

c===== initialisation des parametres de regression/pas effets aleatopires


c         write(*,*),'====================================='
c         write(*,*),'== avec var explicatives !t ============='
c         write(*,*),'====================================='

c         write(*,*),'a ver si ok',k0

         effet=0
         np = nst*n + nva

         call marq982(k0,b,np,ni,v,res,ier,istop)
         
c         write(*,*),'================================================'
c         write(*,*),'== avec var explicatives + effet groupe ========='
c         write(*,*),'================================================'


c         write(*,*),'a ver si ok',k0

         effet=1
         do 256 i=1,nva
            b(np-i+2)=b(np-i+1)
 256     continue
         b(np-nva+1)=0.1d0
         np = nst*n + nva +effet
        call marq982(k0,b,np,ni,v,res,ier,istop)

c         write(*,*),'================================================'
c         write(*,*),'== avec var explicatives + effet sub group ====='
c         write(*,*),'================================================'

c         write(*,*),'a ver si ok',k0

         bgpe=b(np-nva)
         gaux=g !numéro de groupe stokés
         gnew=0
         do 257 i=1,ngexact
            do 258 j=1,nsujet
               if(g(j).eq.i)then
                  gnew(j)=ssg(j,i)
               endif
 258        continue
 257     continue   
         g=gnew

         effet=1
         ngaux=ngexact
         ngexact=nssgexact
         b(np-nva)=0.5d0
         np = nst*n + nva +effet      
        call marq982(k0,b,np,ni,v,res,ier,istop)
         bssgpe=b(np-nva)
         
c         write(*,*),'=========================================='
c         write(*,*),'== ensemble des parametres =============',k0(1)
c         write(*,*),'====================================='
         effet=2
         ngexact = ngaux
         g=gaux
         do 255 i=1,nva
            b(np-i+2)=b(np-i+1)
 255     continue

         np=nbpara 
         b(np-nva-1)=bgpe!0.1d0!0.15d0  !initialisation alpha(groupe)
         b(np-nva)=bssgpe!1.0d0!0.15d0   !initialisation eta(sous groupe)

c        write(4,*)'======================================'
c        write(4,*)'=== Initialisation parametres ========'
c        write(4,*)'alpha  = ',b(np-nva-1)*b(np-nva-1)
c        write(4,*)'eta = ',b(np-nva)*b(np-nva)
c        write(4,*)'======================================'

         call marq982(k0,b,np,ni,v,res,ier,istop)
         
c         write(4,*) '*********************************'
c        write(4,*) 'nombre d iteration', ni
c         write(4,*) '*********************************'

         j=(np-nva)*(np-nva+1)/2
         
         trace=0
         trace1=0
         trace2=0

c%%%% trace(H-1 * I)=mdf :
c strate1 : 
             do 772 i=1,nz1+2
                do 773 j=1,nz1+2
                   H1_hess(i,j)=H_hess(i,j)
                   I1_hess(i,j)=I_hess(i,j)
 773            continue
 772         continue 
             call multi(H1_hess,I1_hess,nz1+2,nz1+2,nz1+2,HI1)
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
                call multi(H2_hess,I2_hess,nz2+2,nz2+2,nz2+2,HI2)
                do 777 i =1,nz2+2
                   trace2=trace2+HI2(i,i)
 777            continue

c                write(4,*)'trace2',trace2
             endif ! pour la deuxieme strate

c             call multi(H_hess,I_hess,np,np,np,HI)
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
c fin calcul trace

         call multi3(I_hess,H_hess,np,np,np,IH)
         call multi3(H_hess,IH,np,np,np,HIH)


c         if (effet.eq.2)then

c            if(AG.EQ.1)then
c            write(4,*)'*************************************** ' 
c            write(4,*)'**** ANDERSEN-GILL APPROACH *********** ' 
c            write(4,*)'*************************************** ' 
c            endif

c            write(*,*)'=============================================='
c            write(4,*)'=============================================='
c            write(*,*)'=> alpha(effet groupe)=',b(np-nva-1)*b(np-nva-1) 
c            write(*,*)'SD(alpha)=H-1:'
c            write(*,*)
c     &       dsqrt(((2.d0*b(np-nva-1))**2)*H_hess(np-nva-1,np-nva-1))
c     &,b(np-nva-1),H_hess(np-nva-1,np-nva-1)
            
c           write(*,*)'SD(alpha)=H-1 I H -1 :'
c           write(*,*)dsqrt(((2.d0*b(np-nva-1))**2)*HIH(np-nva-1,np-nva-1))
c            write(*,*)''
c            write(*,*)'=> eta(effet sous groupe)=',b(np-nva)*b(np-nva) 
c            write(*,*)'SD(eta)=H-1:'
c            write(*,*)dsqrt(((2.d0*b(np-nva))**2)*H_hess(np-nva,np-nva))
c     &           ,b(np-nva),H_hess(np-nva,np-nva)
c            write(*,*)'SD(eta)=H-1 I H -1 :'
c            write(*,*)dsqrt(((2.d0*b(np-nva))**2)*HIH(np-nva,np-nva))
c            write(4,*)'=============================================='
            
c            write(4,*)'=> alpha(effet groupe)=',b(np-nva-1)*b(np-nva-1) 
c            write(4,*)'SD(alpha)=H-1:'
c         write(4,*)dsqrt
c     &(((2.d0*b(np-nva-1))**2)*H_hess(np-nva-1,np-nva-1))
c               !  =((2.d0*b(np-nva-1))**2)*v(j) avec j=(np-nva)(np-nva+1)/2
c            write(4,*)'SD(alpha)=H-1 I H -1 :'
c         write(4,*)dsqrt(((2.d0*b(np-nva-1))**2)*HIH(np-nva-1,np-nva-1))
c            write(4,*)''
c            write(4,*)'=> eta(effet sous groupe)=',b(np-nva)*b(np-nva) 
c            write(4,*)'SD(eta)=H-1:'
c            write(4,*)dsqrt(((2.d0*b(np-nva))**2)*H_hess(np-nva,np-nva))
c            write(4,*)'SD(eta)=H-1 I H -1 :'
c            write(4,*)dsqrt(((2.d0*b(np-nva))**2)*HIH(np-nva,np-nva))
           
c            do 800 ig=1,ngexact
c               do 801 issg=1,nssgbyg
c            write(*,*)'nb deces par gpe et ss gpe',mij(ig,issg),ig,issg
c            write(4,*)'nb deces par gpe et ss gpe',mij(ig,issg),ig,issg
c 801     continue 
c 800  continue 

c         endif
         

         if(effet.eq.2.and.ier.eq.-1)then
            v((np-nva)*(np-nva+1)/2)=10.d10
          endif
         
c          write(*,*)'valeur de ni',ni 
       

c          j=(np-nva)*(np-nva+1)/2

c         wald=(b(np-nva)*b(np-nva))
c         wald=wald/dsqrt(H_hess(np-nva,np-nva))
c         write(8,*)'WALD',wald
c         write(8,*)' '
c         write(8,*)' '
c         endif
c         if(nva.gt.0)then
c         do 123 i=1,nva
c            j=(np-nva+i)*(np-nva+i+1)/2
c            bi = b(np-nva+i) - 1.96*dsqrt(v(j))
c            bs = b(np-nva+i) + 1.96*dsqrt(v(j))
c            write(8,*)' '
c            write(8dsqrt(v(j))
c            bs = b(np-nva+i) + 1.96*dsqrt(v(j))
c            write(8,*)' '
c            write(8,*)' '
c            write(8,*)i, '  beta = ',b(np-nva+i),'  SE',dsqrt(v(j))
c            write(8,*)'WALD',(b(np-nva+i))/dsqrt(H_hess(np-nva,np-nva))
c            write(8,*)' '
c       write(8,*)'RR : ',dexp(b(,*)' '
c            write(8,*)i, '  beta = ',b(np-nva+i),'  SE',dsqrt(v(j))
c            write(8,*)'WALD',(b(np-nva+i))/dsqrt(H_hess(np-nva,np-nva))
c            write(8,*)' '
c       write(8,*)'RR : ',dexp(b(np-nva+i)),'  IC',dexp(bi),dexp(bs)
c 123  continue  
c      endif
        
c         write(*,*)' '
c         j=(np-nva)*(np-nva+1)/2
c         if (effet.eq.1)then
c         write(*,*)'theta', b(np-nva)*b(np-nva)
c         write(*,*)'SE',dsqrt(v(j))
         
c         endif
c         if(nva.gt.0)then
c         do 124 i=1,nva
c            j=(np-nva+i)*(np-nva+i+1)/2
c            bi = b(np-nva+i) - 1.96*dsqrt(H_hess(np-nva+i,np-nva+i))
c            bs = b(np-nva+i) + 1.96*dsqrt(H_hess(np-nva+i,np-nva+i))
c            write(4,*)'**************** '
c            write(4,*)'Variable : ',nomvar(i)
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



c  Output JRG Aug'07

           resOut=res

           do ss=1,npmax
            do sss=1,npmax
             HIHOut(ss,sss) = HIH(ss,sss)
             H_hessOut(ss,sss)= H_hess(ss,sss)
            end do  
           end do



c --------------  Lambda and survival estimates 
      
c      if(effet.eq.2)then
c         call distance2(nz1,nz2,b,fich3,fich4,fich3b,fich4b,
c     & fich3c,fich4c,effet)
c      endif

       call distance2(nz1,nz2,b,effet,
     &	       x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)      


  
      end subroutine nested









  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC                             CCCCCCCCCCCCCCC
CCCCCCCCCCCCC          SUBROUTINES        CCCCCCCCCCCCCCC
CCCCCCCCCCCCC                             CCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c========================== VECSPLI =====================
      subroutine vecspli2(n,ndate) 
      
c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************

      integer ::  n,ndate,i,j,k
      double precision :: ht,htm,h2t,ht2,ht3,hht,h,hh,h2
      double precision :: h3,h4,h3m,h2n,hn,hh3,hh2
         

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
      subroutine vecpen2(n) 

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
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


c================================  SEARPAS joly    ==============================

      SUBROUTINE searpas2(VW,STEP,B,BH,M,DELTA,FIM,EPSV,k0)
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
       CALL valfpa2(VLW1,FI1,B,BH,M,DELTA,k0)
       CALL valfpa2(VLW2,FI2,B,BH,M,DELTA,k0)
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
          CALL valfpa2(VLW1,FI1,B,BH,M,DELTA,k0)   
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
          CALL valfpa2(VLW1,FI1,B,BH,M,DELTA,k0)
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
         CALL valfpa2(VM,FIM,B,BH,M, DELTA,k0)
         IF (FIM.LE.FI2) GO TO 100
         VM=VLW2
         FIM=FI2
100   CONTINUE
      VW=DEXP(VM)
      RETURN

      END


   
c===================================   VALFPA joly   ==============================

        subroutine valfpa2(vw,fi,b,bk,m,delta,k0)
        integer :: m,i
        double precision :: vw,fi
        double precision :: funcpa2,z        
        double precision ,dimension(2) :: k0
        double precision ,dimension(M) :: b,bk,delta

        z=0.d0
        do 1 i=1,m
           bk(i)=b(i)+dexp(vw)*delta(i)
 1      continue
c     
        fi=-funcpa2(bk,m,1,z,1,z,k0)
c
         return
         end   
 
c=================================    DERIVA pjoly =======================

      subroutine deriva2(b,m,v,rl,k0)

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************

      integer ::i0,iun,m,m1,ll,i,k,j
      double precision :: funcpa2,thn,th,z,rl,vl
      double precision :: th2
      double precision ,dimension(2):: k0
      double precision ,dimension(m):: fcith
      double precision ,dimension(m):: b
      double precision , dimension(m*(m+3)/2)::V
c      double precision ,dimension(npmax):: fcith
c      double precision ,dimension(npmax):: b
c      double precision , dimension(npmax*(npmax+3)/2)::V
      
c     v:matrice d'information+score
c     calcul de la derivee premiere
c     
c      print*,'entree deriva'
      th=1.d-5 !5.d-4 !5.d-5!
      thn=-th
      th2=th*th
      z=0.d0
      i0=0
      iun =1
      rl=funcpa2(b,m,iun,z,iun,z,k0)
      do 2 i=1,m
         fcith(i)=funcpa2(b,m,i,th,i0,z,k0)
c         print*,'*** fcith :',fcith(i),b(i),i,m,k0(1)
 2    continue
      
      k=0
      m1=m*(m+1)/2
      ll=m1
      do 1 i=1,m
         ll=ll+1
         vl=(fcith(i)-funcpa2(b,m,i,thn,i0,z,k0))/(2.d0*th)
c         print*,'*** - funcpa :', funcpa2(b,m,i,thn,i0,z,k0),b(i),i,m
         v(ll)=vl
         do 1 j=1,i
            k=k+1
            v(k)=-(funcpa2(b,m,i,th,j,th,k0)-fcith(j)-fcith(i)+rl)/th2 
c            print*,'*** v(k)',v(k),k,j,i,k0
c     &,funcpa2(b,m,i,th,j,th,k0),fcith(j),fcith(i),rl,th2
 1       continue
         
         return
         end


c==========================  DISTANCE   =================================
    
   
         subroutine distance2(nz1,nz2,b,effet,
     &	       x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)

c***************SAME **********************************************
      integer npmax,nsujetmax,nvarmax,ngmax,nssgmax,
     &   nboumax,ndatemax,nptsmax 
      parameter(npmax=50,nsujetmax=20000,nvarmax=50)
      parameter(ngmax=2500,nssgmax=5000,nboumax=2000)
      parameter(ndatemax=30000,nptsmax=1500)
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
             call cosp2(x1,the1,nz1+2,hes1,zi,date,binf,su,bsup
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
               call cosp2(x2,the2,nz2+2,hes2,zi,date,binf,su,bsup
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
      subroutine susp2(x,the,n,su,lam,zi)

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
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

      subroutine cosp2(x,the,n,y,zi,date,binf,su,bsup,lbinf,lam,lbsup)

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
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

         call conf02(x,j,n,y,pm,zi,date)

         binf = dexp(-gl + 1.96d0*pm)
         su  = dexp(-gl)
         bsup = dexp(-gl - 1.96d0*pm)

         call conf12(x,j,n,y,pm,zi,date)
         lbinf = lam - 1.96d0*pm
         lbsup = lam + 1.96d0*pm

         return

         end
c=====================  CONF1  =============================
      subroutine  conf12(x,ni,n,y,pm,zi,date)

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************

      integer :: ni,i,n,j

      double precision :: mmsp2,x,pm
      double precision :: res
      double precision,dimension(ndatemax) :: date
      double precision,dimension(-2:npmax) :: zi
      double precision,dimension(npmax) :: vecti,aux
      double precision,dimension(npmax,npmax) :: y
      
            do 10 i=1,n
               vecti(i) = mmsp2(x,ni,i,zi,date)
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
      subroutine  conf02(x,ni,n,y,pm,zi,date)

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************
       
       integer :: ni,i,n,j
       double precision :: isp2,x,pm
       double precision :: res
       double precision,dimension(-2:npmax) :: zi
       double precision,dimension(ndatemax) :: date
       double precision,dimension(52) :: vecti,aux
       double precision,dimension(npmax,npmax) :: y


            do 10 i=1,n
               vecti(i) = isp2(x,ni,i,zi,date)
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
          double precision function isp2(x,ni,ns,zi,date)

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************
  
          integer :: ni,ns
          double precision :: val,mmsp2,x
          double precision ,dimension(-2:npmax)::zi
          double precision ,dimension(ndatemax)::date

          if(x.eq.zi(ni))then
             if(ni.le.ns-3)then
                val = 0.d0
             else
                if(ni.le.ns-2)then
           val = ((zi(ni)-zi(ni-1))*mmsp2(x,ni,ns,zi,date))*0.25d0
                else
                   if (ni.eq.ns-1)then
                  val = ((zi(ni)-zi(ni-2))*mmsp2(x,ni,ns,zi,date)+
     &        (zi(ni+3)-zi(ni-1))*mmsp2(x,ni,ns+1,zi,date))*0.25d0
                   else
                      if(ni.eq.ns)then
                  val = ((zi(ni)-zi(ni-3))*mmsp2(x,ni,ns,zi,date)+
     &                (zi(ni+2)-zi(ni-2))*mmsp2(x,ni,ns+1,zi,date)
     &       +(zi(ni+3)-zi(ni-1))*mmsp2(x,ni,ns+2,zi,date))*0.25d0
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
                 val = (x-zi(ni))*mmsp2(x,ni,ns,zi,date)*0.25d0
             else  
             if(ni.eq.ns-2)then
                   val = ((x-zi(ni-1))*mmsp2(x,ni,ns,zi,date)+
     &       (zi(ni+4)-zi(ni))*mmsp2(x,ni,ns+1,zi,date))*0.25d0
             else   
                if (ni.eq.ns-1)then
                   val =((x-zi(ni-2))*mmsp2(x,ni,ns,zi,date)+
     &             (zi(ni+3)-zi(ni-1))*mmsp2(x,ni,ns+1,zi,date)
     &       +(zi(ni+4)-zi(ni))*mmsp2(x,ni,ns+2,zi,date))*0.25d0
                else
                   if(ni.eq.ns)then
                      val =((x-zi(ni-3))*mmsp2(x,ni,ns,zi,date)+
     &             (zi(ni+2)-zi(ni-2))*mmsp2(x,ni,ns+1,zi,date)
     &             +(zi(ni+3)-zi(ni-1))*mmsp2(x,ni,ns+2,zi,date)
     &        +(zi(ni+4)-zi(ni))*mmsp2(x,ni,ns+3,zi,date))*0.25d0
                   else
                      val = 1.d0
                   endif
                endif
             endif
             endif
          endif 
          endif
             isp2 = val
             return
             end
c==========================  MMSP   ==================================
      double precision function mmsp2(x,ni,ns,zi,date)

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
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

             mmsp2 = val
             return
             end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



c=======cross validation

c========================          MNBRAK         ===================
      subroutine mnbrak2(ax,bx,cx,fa,fb,fc,b,n)

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************

         double precision ax,bx,cx,fa,fb,fc,aux,res
         double precision b(npmax),y(npmax,npmax)
         double precision estimv2,gold,glimit,tiny
         parameter (gold=1.618034d0,glimit=100.d0,tiny=1.d-20)
         double precision dum,fu,q,r,u,ulim
         integer n,ni

c     write(*,*)' DEBUT mnbrak',fa,fb
         fa = estimv2(ax,n,b,y,aux,ni,res)
ccc   write(*,*)'fa mnbrak',fa
c     stop
         fb = estimv2(bx,n,b,y,aux,ni,res)
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
         fc = estimv2(cx,n,b,y,aux,ni,res)
c         write(*,*)'ok mnbrak'
c         stop
 1       if(fb.ge.fc)then
            r = (bx-ax)*(fb-fc)
            q = (bx-cx)*(fb-fa)
            u = bx-((bx-cx)*q-(bx-ax)*r)/
     &       (2.d0*sign(max(abs(q-r),tiny),q-r))
            ulim = bx + glimit*(cx-bx)
            if((bx-u)*(u-cx).gt.0.d0)then
               fu = estimv2(u,n,b,y,aux,ni,res)
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
               fu = estimv2(u,n,b,y,aux,ni,res)
            else
               if((cx-u)*(u-ulim).gt.0.d0)then
                  fu = estimv2(u,n,b,y,aux,ni,res)
                  if(fu.lt.fc)then
                     bx = cx
                     cx = u
                     u = cx + gold*(cx-bx)
                     fb = fc
                     fc = fu
                     fu = estimv2(u,n,b,y,aux,ni,res)
                  endif  
               else
                  if((u-ulim)*(ulim-cx).ge.0.d0)then
                     u = ulim
                     fu = estimv2(u,n,b,y,aux,ni,res)
                  else
                     u = cx + gold*(cx-bx)
                     fu = estimv2(u,n,b,y,aux,ni,res)
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
      double precision function golden2(ax,bx,cx,tol,xmin,n,b,y,aux)

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************
     
      double precision y(npmax,npmax)
      double precision ax,bx,cx,tol,xmin,b(npmax)
      double precision r,c,aux,res
      parameter (r=0.61803399d0,c=1.d0-r)
      double precision f1,f2,x0,x1,x2,x3,estimv2
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
         f1 = estimv2(x1,n,b,y,aux,ni,res)
         f2 = estimv2(x2,n,b,y,aux,ni,res)
         
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
               f2 = estimv2(x2,n,b,y,aux,ni,res)
               
c         write(*,*)'f2 DANS golden',f2,f1
c         stop
            else
               x3 = x2
               x2 = x1
               x1 = r*x2+c*x0
               f2 = f1
c         write(*,*)'f2 else',f1
               f1 = estimv2(x1,n,b,y,aux,ni,res)
c         write(*,*)'f1 else',f1
c         stop
            endif
            go to 1
          endif
          if(f1.lt.f2)then
             golden2 = f1
             xmin = x1
          else
             golden2 = f2
             xmin = x2
          endif
          return
          end


c========================          ESTIMV         ===================

      double precision function estimv2(k00,n,b,y,aux,ni,res)

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
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
      double precision pe
c yas      integer :: effet,nz1,nz2,effetcross
   	  integer effet,nz1,nz2
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
c      write(*,*)'DEBUT ESTIMV',effet ,k0
c      stop
      call marq982(k0,b,n,ni,v,res,ier,istop)
      
      
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
         
         call test2(bh,ut,dut,k0,n,aux,v,y)
         estimv2 = - ((res-pe)) - aux
c         write(*,*)'estimv',k0,-aux,ni,estimv
      else
         aux = -n
c     write(*,*)'estimv2',k0,-aux,ni
      endif
      
      return
      end
      
c=================calcul de la hessienne  et de omega  ==============
      subroutine test2(b,ut,dut,k0,n,res,v,y)

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
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
               call mat2(hess(i,j),ut,dut,i,j,n)
c               write(*,*)'hess test',hess(i,j),i,j
 15         continue
 20      continue
        do 40 i = 2,n
            do 35 j = 1,i-1
               hess(i,j)=hess(j,i)
 35         continue
 40      continue


         call calcomeg2(n,omeg)

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
   
c         write(*,*)'test',k0(1),-tra
         res = (tra)

         end




c=======================  CALOMEG  ===========================
      subroutine calcomeg2(n,omeg)
c        remplissage de la matrice omega n*n
c          elle a 7 diagonales


c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************


      double precision omeg(npmax,npmax)
      integer n
      double precision calc002,calc012,calc022
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
   
      omeg(1,1)=calc002(1,n)
      omeg(1,2)=calc012(1,n)
      omeg(1,3)=calc022(1,n)
      omeg(1,4)=m3m(1)
      omeg(2,1)=omeg(1,2)
      omeg(2,2)=calc002(2,n)
      omeg(2,3)=calc012(2,n)
      omeg(2,4)=calc022(2,n)
      omeg(2,5)=m3m(2)
      omeg(3,1)=omeg(1,3)
      omeg(3,2)=omeg(2,3)
      omeg(3,3)=calc002(3,n)
      omeg(3,4)=calc012(3,n)
      omeg(3,5)=calc022(3,n)
      omeg(3,6)=m3m(3)
      do 10 i=4,n-3
         omeg(i,i-3)=omeg(i-3,i)
         omeg(i,i-2)=omeg(i-2,i)
         omeg(i,i-1)=omeg(i-1,i)
         omeg(i,i)=calc002(i,n)
         omeg(i,i+1)=calc012(i,n)
         omeg(i,i+2)=calc022(i,n)
         omeg(i,i+3)=m3m(i)
 10   continue   
      omeg(n-2,n-5)=omeg(n-5,n-2)
      omeg(n-2,n-4)=omeg(n-4,n-2)
      omeg(n-2,n-3)=omeg(n-3,n-2)
      omeg(n-2,n-2)=calc002(n-2,n)
      omeg(n-2,n-1)=calc012(n-2,n)
      omeg(n-2,n)=calc022(n-2,n)
      omeg(n-1,n-4)=omeg(n-4,n-1)
      omeg(n-1,n-3)=omeg(n-3,n-1)
      omeg(n-1,n-2)=omeg(n-2,n-1)
      omeg(n-1,n-1)=calc002(n-1,n)
      omeg(n-1,n)=calc012(n-1,n)
      omeg(n,n-3)=omeg(n-3,n)
      omeg(n,n-2)=omeg(n-2,n)
      omeg(n,n-1)=omeg(n-1,n)
      omeg(n,n)=calc002(n,n)

      end


c====================  MAT  ==================================
      subroutine mat2(res,ut,dut,k,l,n)

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************

      double precision res,dut(ndatemax),ut(ndatemax)
      integer k,l,j,ni,n
      double precision res1,msp2,aux2
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
                aux2 = msp2(nt1(i),ni,k)*msp2(nt1(i),ni,l)
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
          double precision function msp2(i,ni,ns)
 
c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
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

             msp2 = val
             return
             end
c==========================   SP   ==================================
          double precision function sp2(i,ni,ns)

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************
         
             integer ni,ns,i
             double precision val,msp2
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
                   val = ((zi(ni)-zi(ni-1))*msp2(i,ni,ns))*0.25d0
                else
                   if (ni.eq.ns-1)then
                      val = ((zi(ni)-zi(ni-2))*msp2(i,ni,ns)+
     &                (zi(ni+3)-zi(ni-1))*msp2(i,ni,ns+1))*0.25d0
                   else
                      if(ni.eq.ns)then
                         val = ((zi(ni)-zi(ni-3))*msp2(i,ni,ns)+
     &                       (zi(ni+2)-zi(ni-2))*msp2(i,ni,ns+1)
     &              +(zi(ni+3)-zi(ni-1))*msp2(i,ni,ns+2))*0.25d0
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
                   val = (date(i)-zi(ni))*msp2(i,ni,ns)*0.25d0
             else  
             if(ni.eq.ns-2)then
                   val = ((date(i)-zi(ni-1))*msp2(i,ni,ns)+
     &             (zi(ni+4)-zi(ni))*msp2(i,ni,ns+1))*0.25d0
             else   
                if (ni.eq.ns-1)then
                   val =((date(i)-zi(ni-2))*msp2(i,ni,ns)+
     &             (zi(ni+3)-zi(ni-1))*msp2(i,ni,ns+1)
     &             +(zi(ni+4)-zi(ni))*msp2(i,ni,ns+2))*0.25d0
                else
                   if(ni.eq.ns)then
                      val =((date(i)-zi(ni-3))*msp2(i,ni,ns)+
     &             (zi(ni+2)-zi(ni-2))*msp2(i,ni,ns+1)
     &             +(zi(ni+3)-zi(ni-1))*msp2(i,ni,ns+2)
     &             +(zi(ni+4)-zi(ni))*msp2(i,ni,ns+3))*0.25d0
                   else
                      val = 1.d0
                   endif
                endif
             endif
             endif
          endif 
          endif
             sp2 = val
             return
             end
c================

c=========================  calc00  =========================
          double precision function calc002(j,n) 

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************
      
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

             calc002 = part
             return
             end
c=========================  CALC01  =========================
      double precision function calc012(j,n)
  
c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************
 
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

             calc012 = part
             return
             end
c=========================  CALC02  =========================
          double precision function calc022(j,n)

 
c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************


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

             calc022 = part
             return
             end

c===============================    MARQ98  HJ =========================

   	 subroutine marq982(k0,b,m,ni,v,rl,ier,istop)

c  fu = matrice des derivees secondes et premieres
c
c  istop: raison de l'arret
c  1: critere d'arret satisfait (prm=ca, vraisblce=cb, derivee=dd)
c  2: nb max d'iterations atteints
c  3: 1 mais echec inversion matrice d'info (ier=1)
C  4: non amelioration vraisblce apres searpas
 
c*************** SAME **********************************************
         integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
         integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
         integer , parameter ::ndatemax=30000
         integer , parameter ::nptsmax=1500
c******************************************************************

         integer :: m,ni,nql,i,ii,nfmax,idpos,ier,istop,igrad,j
         integer :: ncount,id,jd,i0

         double precision :: da,dm,ga,tr,g0
         double precision :: ca,cb,epsa,epsb,rl
         double precision :: funcpa2,det, step,eps,epsd
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
c*****dace7
      double precision , dimension(npmax,npmax) :: I_hess,H_hess
      double precision , dimension(npmax,npmax) :: Hspl_hess
      double precision , dimension(npmax,1) ::PEN_deri
      common /dace7/PEN_deri,I_hess,H_hess,Hspl_hess

c add JRG after Virginie's e-mail
c*****dace3
      double precision pe
c yas      integer :: effet,nz1,nz2,effetcross
   	  integer effet,nz1,nz2
      common /dace3/pe,effet,nz1,nz2



c***********************************************

c*****dace1 
c      double precision , dimension(ndatemax) ::date
c      double precision , dimension(-2:npmax) :: zi
c     common /dace1/date,zi


c      write(*,*) "bamsmsms: ", date
c      write(*,*) "bamsmsms: ", zi


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
      rl=funcpa2(b,m,i0,z,i0,z,k0)
      rl1=rl 

c      write(*,*) 'hola:',b,m,i0,z


c      write(*,*)'log Vrais rl1,nsujet',rl1,nsujet
c      stop
      if (RL1.gt.0)stop
      
      call deriva2(b,m,v,rl,k0)

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
          call dchole(fu,m,nql,idpos)
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
              rl=funcpa2(b1,m,id,z,jd,z,k0)
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
           call searpas2(vw,step,b,bh,m,delta,fi,eps,k0) 
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

           call deriva2(b,m,v,rl,k0)
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
           call dsinv(v,m,ep,ier)
           if (ier.eq.-1) then
c             write(6,103)          
 103         format(1x,'echec inversion matrice information')

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
         call deriva2(b,m,vnonpen,rl,zero)
         
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


c==============================================================

c========================          FUNCPA  NESTED       ====================
      double precision function funcpa2(b,np,id,thi,jd,thj,k0)

c      USE nr, ONLY : gammln

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************


      integer :: nb,n,np,id,jd,i,j,k,vj,cptg,l
      integer :: ig,ip,issg,cptaux,choix
      integer , dimension(ngmax) :: cpt

      real , dimension(nptsmax) :: points,poids
      REAL :: gammln,essai

c      double precision :: func12,func22
      double precision :: thi,thj,pe1,pe2,dnb,sum
      double precision :: theta,inv
      double precision :: som1,som2
      double precision :: res,vet,h1
      double precision :: int1,int2,int

      double precision , dimension(-2:npmax) ::the1,the2
      double precision , dimension(np) :: b,bh
      double precision , dimension(ngmax) :: res1,res2,res3
      double precision , dimension(ngmax) :: integrale1,integrale2
      double precision , dimension(ngmax) :: integrale3 !andersen gill
      double precision , dimension(2) :: k0
      double precision , dimension(ndatemax) :: dut1,dut2
      double precision , dimension(0:ndatemax) :: ut1,ut2
      double precision , dimension(ngmax) :: sum1

c ************************* les commons :
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
      double precision pe
c yas      integer :: effet,nz1,nz2,effetcross
   	  integer effet,nz1,nz2
      common /dace3/pe,effet,nz1,nz2
c*****ve1
      double precision , dimension(nsujetmax,npmax):: ve
      common /ve1/ve
c*****contrib
      integer :: ngexact,nssgexact    !nb EXACT de gpes et ss gpes au total
      integer :: nssgbyg  !nb de ss gpes par groupe au MAX
      common /contrib/ngexact,nssgexact,nssgbyg 
c*****groupe
      integer , dimension(nsujetmax) :: g
      integer , dimension(nsujetmax,ngmax) :: ssg
! numero du ssgpe pour un sujet et un gpe donné
      integer , dimension(ngmax) :: nig
      common /gpe/g,nig,ssg
c****** indicateur de troncature
      integer :: indictronq ! =0 si donnees non tronquées
      common /troncature/indictronq
c******  aux1 aux2
      double precision , dimension(ngmax,nssgmax) ::  aux1,aux2
      common /risqcumul/aux1,aux2
c*****auxig
      integer :: auxig
      common /auxig/auxig
c*******alphaeta
      double precision :: alpha,eta
      common /alphaeta/alpha,eta
c*****mij
      integer , dimension(ngmax) :: mi ! nb de dc dans gpe i
      integer  , dimension(ngmax,nssgmax) ::  mij ! nb de dc dans ssgpe ij
      common /mij/mij,mi
c %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
      integer AG
      common /andersengill/AG
c***********************************************************
      
    
c      write(*,*)'%% ENTREE DANS FUNCPA %'
c     on force alpha a 0.17 = 0.41 * 0.41
c      bh(np-nva-1)=0.41

      do 3 i=1,np
         bh(i)=b(i)
 3    continue 

      if (id.ne.0) bh(id)=bh(id)+thi 
      if (jd.ne.0) bh(jd)=bh(jd)+thj    
      
      n = nz1+2
      
      do 4 i=1,n
         the1(i-3)=(bh(i))*(bh(i))
         j = n+i 
         if (nst.eq.2) then
            the2(i-3)=(bh(j))*(bh(j))
         endif
 4    continue

      if(effet.eq.1) then
         theta = (bh(np-nva)*bh(np-nva)) ! variance effet groupe
      endif

      if(effet.eq.2) then
         alpha = (bh(np-nva-1)*bh(np-nva-1)) ! variance effet groupe
         eta = (bh(np-nva)*bh(np-nva))  ! variance effet sous groupe
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
      if(nst.eq.2)then
         ut2(ndate)=som2+the1(i-4)+the2(i-3)+the2(i-2)+the2(i-1)
         dut2(ndate) = (4.d0*the2(i-1)/h1)
      endif

c         write(*,*)'** ut2',(ut2(i),i,i=0,ndate)
c         write(*,*)'** dut1',(dut1(i),i,i=0,ndate)
c         write(*,*)'** date **',(date(i),i,i=1,ndate)

c--------------------------------------------------------
c----------calcul de la vraisemblance ------------------
c---------------------------------------------------------
            
      do 89 ig=1,ngexact
         res1(ig) = 0.d0
         res2(ig) = 0.d0
         cpt(ig) = 0
 89   continue
      
c*******************************************     
C----- sans effet aleatoire dans le modele
c*******************************************     

        if (effet.eq.0) then


           do 110 i=1,nsujet
                   cpt(g(i))=cpt(g(i))+1
c                   write(*,*)' CPT',cpt(g(i)),(g(i)),i
                   
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
                     res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet)
c               write(*,*)'FUNC res2',res2(g(i)),g(i),stra(i),nt1(i) 
                  endif  

                   if((c(i).eq.1).and.(stra(i).eq.2))then
                      res2(g(i)) = res2(g(i))+dlog(dut2(nt1(i))*vet)
                  endif

                   if(stra(i).eq.1)then
                      res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet
     &                     -ut1(nt0(i))*vet 
c             write(*,*)'FUNC res1',res1(g(i)),g(i),stra(i),nt1(i),nt0(i)
                   endif

                   if(stra(i).eq.2)then
                      res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet
     &                     -ut2(nt0(i))*vet 
c                   write(*,*)'FUNC res1',res1(g(i)),g(i),stra(i),nt1(i)
                   endif
           
 110           continue     
c               stop
c           write(*,*)'FUNC fin res1',res1(g(2)),g(2),stra(2),nt1(2)
c           write(*,*)'FUNC',ut1(nt1(i)) ,ut1(nt0(i)) 
c           write(*,*)'FUNC res2',res2(g(1)),vet
c           write(*,*)'FUNC vet',vet,bh,cpt
c           stop
           res = 0.d0
           
          cptg = 0
          
C k indice les groupes
          do 115 k=1,ngexact
             if(cpt(k).gt.0)then !nb de sujets dans un gpe=nig()
                             
                res = res-res1(k) 
     &               + res2(k) 
                cptg = cptg + 1 
c             write(*,*)' **rescpt',res,res1(k),res2(k),cpt(k),k
             endif 
 115       continue

        endif !fin boucle effet=0

c*******************************************
C-----avec un seul  effet aleatoire dans le modele
c*********************************************

        if (effet.eq.1) then
c      write(*,*)'AVEC 1 EFFET ALEATOIRE'
          inv = 1.d0/theta
          cpt=0.d0
          res1=0.d0
          res2=0.d0
          res3=0.d0
c      write(*,*)'INV 1============',inv
C     i indice les sujets
          do 109 i=1,nsujet 
             
             cpt(g(i))=cpt(g(i))+1 
c             write(*,*)'funcpa cpt ******',cpt(g(i)),g(i)
 
             if(nva.gt.0)then
                vet = 0.d0   
                do 199 j=1,nva
                   vet =vet + bh(np-nva+j)*dble(ve(i,j))
c                     write(*,*)'funcpa ve ******',ve(i,j),i,j
 199            continue
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
c              res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet -ut1(nt0(i))*vet
c nouvelle version
              res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet 
             endif
             
             if(stra(i).eq.2)then
c              res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet-ut2(nt0(i))*vet
c nouvelle version
              res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet
             endif
c modification pour nouvelle vraisemblance / troncature:
             if(stra(i).eq.1)then
                res3(g(i)) = res3(g(i)) + ut1(nt0(i))*vet 
             endif
             
             if(stra(i).eq.2)then
                res3(g(i)) = res3(g(i)) + ut2(nt0(i))*vet 
             endif
             
 109      continue 
c          stop
c           write(*,*)'***FUNC res2',(res2(i),i=1,ngexact)
c          stop
          res = 0.d0
          cptg = 0
          mi =0.d0
c     gam2 = gamma(inv)
C k indice les groupes


          do 119 k=1,nsujet
             if(c(k).eq.1)then
                mi(g(k))=mi(g(k))+1
             endif
 119      continue 

          do 159 k=1,ngexact  
             sum=0.d0
             if(cpt(k).gt.0)then
                nb = mi(k)!nb de deces par groupe
                dnb = dble(nb)
c                write(*,*)'nb,dnb',nb,dnb,k
c     gam1 = gamma(dnb + inv) 
                
                if (dnb.gt.1.d0) then
                   do 169 l=1,nb
                      sum=sum+dlog(1.d0+theta*dble(nb-l))
 169               continue
                endif 
                if(theta.gt.(1.d-5)) then
ccccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                   if(AG.EQ.1)then
                   res= res-(inv+dnb)*dlog(theta*(res1(k)-res3(k))+1.d0) 
     &                  + res2(k) + sum  
ccccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                   else
                   res= res-(inv+dnb)*dlog(theta*(res1(k))+1.d0) 
     &                  +(inv)*dlog(theta*res3(k)+1.d0) 
     &                  + res2(k) + sum  
c                   write(*,*)'***',res2(k),sum,inv,theta,k
                   endif

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
ccccccccccccccccccccccccccccccccccccccccccccccc

                   endif
             endif 
 159      continue

c       stop
       endif !fin boucle effet=1

c*******************************************
C-----avec deux effets aleatoires dans le modele
c*********************************************

        if (effet.eq.2) then

         mi=0
         mij=0
         res1=0.d0
         res2=0.d0
         aux1=0.d0
         aux2=0.d0
         integrale1=0.d0
         integrale2=0.d0

c     === MODIFICATION DE LA VRAISEMBLANCE POUR LE NESTED FRAILTY MODEL
          
          do 11 k=1,nsujet
             if(c(k).eq.1)then
                mi(g(k))=mi(g(k))+1
                mij(g(k),ssg(k,g(k)))=mij(g(k),ssg(k,g(k)))+1 
!nb de dc ds ss gpe ssg(k)                
c         write(*,*)'**par gpe/ssgpe',
c     &               mij(g(k),ssg(k,g(k))),g(k),ssg(k,g(k)),k
             endif
 11       continue 

c          do 800 i=1,nsujet
c           write(*,*)'num de groupe/sujet',g(i),i          
c           write(*,*)'num de ss groupe/sujet et gpe',ssg(i,g(i)),i,g(i)
c 800      continue 
c          stop

          do 10 k=1,nsujet
             if(nva.gt.0)then
                vet = 0.d0 
                do 19 ip=1,nva
                   vet =vet + bh(np-nva+ip)*dble(ve(k,ip))
c                   write(*,*)'ici vet**',vet,bh(np-nva+ip),ve(k,ip),k,ip
 19             continue
                vet = dexp(vet)
c                write(*,*)'**vet**',vet 
             else
                vet=1.d0
             endif

             if((c(k).eq.1).and.(stra(k).eq.1))then
                res2(g(k)) = res2(g(k))+dlog(dut1(nt1(k))*vet)
c           write(*,*)'***res2',res2(g(k)),dut1(nt1(k)),nt1(k),vet,k,g(k)
             endif  
             if((c(k).eq.1).and.(stra(k).eq.2))then
                res2(g(k)) = res2(g(k))+dlog(dut2(nt1(k))*vet)
c           write(*,*)'***res2',res2(g(k)),dut1(nt1(k)),nt1(k),vet,k,g(k)
             endif 
 
             if(stra(k).eq.1)then
           aux1(g(k),ssg(k,g(k)))=aux1(g(k),ssg(k,g(k)))+ut1(nt1(k))*vet
           aux2(g(k),ssg(k,g(k)))=aux2(g(k),ssg(k,g(k)))+ut1(nt0(k))*vet  
c       write(*,*)'*** aux1***',aux1(g(k),ssg(k,g(k))),g(k),ssg(k,g(k)),k
             endif
             if(stra(k).eq.2)then
           aux1(g(k),ssg(k,g(k)))=aux1(g(k),ssg(k,g(k)))+ut2(nt1(k))*vet
           aux2(g(k),ssg(k,g(k)))=aux2(g(k),ssg(k,g(k)))+ut2(nt0(k))*vet
             endif 
 10       continue  
c          write(*,*)'*** FIN aux1***',aux1(10,10),aux1(10,9)
c          stop
c================== calcul des intégrales par Gauss LAGUERRE
c     memes points et poids dans chq groupe

          do 170 ig=1,ngexact
             auxig=ig
             choix=1
             call gaulag2(int,choix)
             integrale1(auxig)=int
c             call qgauss1(0.d0,5.d0,int1)
c             integrale1(auxig)=int1
c             write(*,*)'**** CHOIX',choix

c     integrale sur la troncature: 
             if(indictronq.eq.1)then
                if(AG.eq.1)then !andersen gill
                   choix=3
                   call gaulag2(int,choix)
                   integrale3(auxig)=int   
c                   write(*,*)'**** CHOIX',choix
                else !troncature classique
                   choix=2
c                   write(*,*)'**** CHOIX',choix
                   call gaulag2(int,choix)
                   integrale2(auxig)=int    
c                call qgauss2(0.d0,5.d0,int2)
c                integrale2(auxig)=int2 
c      write(*,*)'**** qgauss2 ****',integrale2(auxig),auxig
                endif
             endif
 170      continue

c=======fin calcul de integrale
c      print*,'** integrale1 **',(integrale1(ig),ig,ig=1,ngexact)
c     &,choix
c      print*,'**** integrale2 **',(integrale2(ig),ig,ig=1,ngexact)
c      stop
c======================================================================

          do 150 ig=1,ngexact  
             sum1(ig)=0.d0
             do 151 issg=1,nssgbyg !!! NON ICI NSSGBYG
c             print*,'**Asum1',nssgbyg  ,mij(ig,issg),ig,issg 
                if(mij(ig,issg).gt.1) then
                   do 16 l=1,mij(ig,issg)
             sum1(ig)=sum1(ig)+dlog(1.d0+eta*dble(mij(ig,issg)-l)) 
c             print*,'**sum1',sum1(ig),mij(ig,issg),issg,ig
 16                continue
                endif
 151         continue
 150       continue

          res = 0.d0
          do 15 k=1,ngexact  
             if(nig(k).gt.0)then
                if(indictronq.eq.0)then
                   res = res
     &                  +res2(k)
     &                  +sum1(k)
     &                -dlog(alpha)/(alpha)-dble(gammln(real(1./alpha)))
     &                  +dlog(integrale1(k))
c        print*,'*** res',res ,integrale1(k),k
                endif
                if(indictronq.eq.1)then
                   if(AG.eq.1)then
ccccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                      res = res
     &                     +res2(k)
     &                     +sum1(k)
     &            -dlog(alpha)/(alpha)-dble(gammln(real(1./alpha)))
     &                     +dlog(integrale3(k))
                   else
c vraisemblance pr donnees censurées dte et tronquées a gauche
                      res = res
     &                     +res2(k)
     &                     +sum1(k)
     &                     +dlog(integrale1(k))
     &                     -dlog(integrale2(k)) 
c
c       print*,'***'
c        print*,'*** 1',    res,+res2(k),k,ngexact
c        print*,'*** 2',   sum1(k)
c        print*,'*** int1',  dlog(integrale1(k)),integrale1(k),k
c        print*,'*** int2',  dlog(integrale2(k)),integrale2(k),k
                   endif
                endif
             endif
c       stop
 15   continue
c      stop

       endif !fin boucle effet=2
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
c           write(*,*)'avant  penalisation res et pe1 ',res,pe1,k0(1)
         endif 
         
         res = res - pe

c     write(*,*)'vraisemblance penalisee',pe,pe1,pe2
c     write(*,*)'vraisemblance penalisee',k0(1),k0(2)
        
         funcpa2= res 
         
c     print*,'*** res :',res,bh(id),id,np
         
c         write(*,*)'%%% SORTIE DANS FUNCPA %%%%',res,pe
c     stop
         return
         end


c===================================================================

      FUNCTION gammln(xx)
      REAL gammln,xx
!     Returns the value ln[gamma(xx)] for xx > 0.
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
!     Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure accuracy is good enough.
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     &    24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     &     -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11      j=1,6
         y=y+1.d0
         ser=ser+cof(j)/y
 11   continue
      gammln=tmp+log(stp*ser/x)
      return
      END

c==================================================================
c==================================================================

      SUBROUTINE gaulag2(ss,choix) 

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************


      double precision :: ss
      double precision ::auxfunca
      double precision ::func0,func12,func22,func3,func4,func5,func6
      external :: func0,func12,func22,func3,func4,func5,func6
c gauss laguerre
c func1 est l intégrant, ss le résultat de l integrale sur 0 ,  +infty
      
      INTEGER :: j,choix

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
c weights (five numbers above each used twice) sum to 2.
      do 11 j=1,20
         if (choix.eq.1) then 
            auxfunca=func12(x(j))
            ss = ss+w(j)*(auxfunca)
         else                   !choix=2, troncature
            if (choix.eq.2) then 
               auxfunca=func22(x(j))
               ss = ss+w(j)*(auxfunca)
               else                   !choix=3,essai, res = -1
                  if (choix.eq.3) then
                     auxfunca=dexp(-x(j))
                     ss = ss+w(j)*(auxfunca)
                  endif
               endif
            endif
 11         continue
c      print*,'** ss',ss,xr
c      stop
      return
      END
  
c================================================



c==================================================================

      SUBROUTINE qgauss1(a,b,ss) ! sans troncature

c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************

      double precision :: a,b,ss,i
      double precision ::auxfunc1a,auxfunc1b
c      external :: func1
c      Returns as ss the integral of the function func1 between a and b, by ten-point Gauss-Legendre integration: the function is evaluated exactly ten times at interior points in the range of integration.
c func1 est l intégrant, ss le résultat de l integrale
      
      INTEGER :: j
      double precision :: func12

c*****auxig
      integer :: auxig
      common /auxig/auxig
c*****mij
      integer , dimension(ngmax) :: mi ! nb de dc dans gpe i
      integer  , dimension(ngmax,nssgmax) ::  mij ! nb de dc dans ssgpe ij
      common /mij/mij,mi
c*****contrib
      integer :: ngexact,nssgexact    !nb EXACT de gpes et ss gpes au total
      integer :: nssgbyg  !nb de ss gpes par groupe
      common /contrib/ngexact,nssgexact,nssgbyg 
 
c*****************************
      DOUBLE PRECISION :: dx,xm,xr,w(5),x(5)  !The abscissas and weights.
      SAVE  w,x
      DATA w/.2955242247,.2692667193,.2190863625,.1494513491,
     &     .0666713443/
      DATA x/.1488743389,.4333953941,.6794095682,.8650633666,
     &     .9739065285/ 

c****************************************

c      write(*,*)'**** qgauss1 ****',(mi(i),i=1,ngexact)
c      write(*,*)'******'

      xm=5.d-1*(b+a)
      xr=5.d-1*(b-a)
      ss=0.d0 ! Will be twice the average value of the function,since the ten
      do 11 j=1,5  !weights (five numbers above each used twice) sum to 2.
      dx=xr*x(j)
c      print*,'** début qgauss1'

c      auxfunc1a=(xm+dx)*(xm+dx) ! essai sur x2
c      auxfunc1b=(xm-dx)*(xm-dx)

      auxfunc1a=func12(xm+dx)
      auxfunc1b=func12(xm-dx)

c      print*,'***** début qgauss1 ',auxfunc1a,auxfunc1b,xm+dx,xm-dx
c      print*,'**prod1 ** ',auxfunc1a,xm+dx,auxig
      ss = ss+w(j)*( auxfunc1a+ auxfunc1b)
c      print*,'** qgauss1 **', ss, auxfunc1a, auxfunc1b,w(j),j
 11   continue
      ss=xr*ss            !  Scale the answer to the range of integration.
c      print*,'** qgauss1 **', ss
c      stop
      return
      END
  
c================================================
c================================================
      double precision  function func12(frail) 
! calcul de l integrant, pour un effet aleatoire donné frail et un groupe donne auxig (cf funcpa)      
 
c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************


      double precision  frail
c      double precision func1
      integer :: ig,issg,k,i
      double precision , dimension(ngmax) :: prod1
c*****auxig
      integer :: auxig
      common /auxig/auxig
c*****ut1ut2
      double precision , dimension(ndatemax) :: dut1,dut2
      double precision , dimension(0:ndatemax) :: ut1,ut2
      common /ut1ut2/dut1,dut2,ut1,ut2
c*****groupe
      integer , dimension(nsujetmax) :: g
      integer , dimension(nsujetmax,ngmax) :: ssg
! numero du ssgpe pour un sujet et un gpe donné
      integer , dimension(ngmax) :: nig
      common /gpe/g,nig,ssg
c*****dace4
      integer , dimension(nsujetmax) :: stra
      common /dace4/stra
c*******alphaeta
      double precision :: alpha,eta
      common /alphaeta/alpha,eta
c*****mij
      integer , dimension(ngmax) :: mi ! nb de dc dans gpe i
      integer  , dimension(ngmax,nssgmax) ::  mij ! nb de dc dans ssgpe ij
      common /mij/mij,mi
c****** indicateur de troncature
      integer :: indictronq ! =0 si donnees non tronquées
      common /troncature/indictronq
c*****contrib
      integer :: ngexact,nssgexact    !nb EXACT de gpes et ss gpes au total
      integer :: nssgbyg  !nb de ss gpes par groupe
      common /contrib/ngexact,nssgexact,nssgbyg 
c*****dace2
      double precision , dimension(nsujetmax) ::t0,t1
      double precision :: ncsrl
      integer , dimension(nsujetmax) ::c
      integer, dimension(nsujetmax) :: nt0,nt1
      integer :: nsujet,nsim,nva,ndate,nst
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,nst
c******  aux1 aux2
      double precision , dimension(ngmax,nssgmax) ::  aux1,aux2
      common /risqcumul/aux1,aux2
            
c============================================================
     
c      write(*,*)'**** func1 ****',(mi(i),i=1,ngexact)
c      write(*,*)'******'

      ig=auxig
c     initialisation de prod1 et prod2 par le numerateur de l integrant
c     prod1(ig)=1.d0
      
             prod1(ig)=
     &            (dexp(dble(-frail/alpha)))
     &            *(dble(frail)**(1.d0/alpha-1.d0+mi(ig)))

c            print*,"** init1 **", prod1(ig),frail,ig
c             stop
      do 1412 issg=1,nssgbyg !!! attention sous gpe pour un gpe donné
         do 1413 k=1,nsujet
            if((g(k).eq.ig).and.(ssg(k,g(k)).eq.issg))then
c      if(eta.gt.1.d-2)then   
               prod1(ig)=prod1(ig)
     &              *(1.d0+eta*dble(frail)*aux1(g(k),ssg(k,g(k))))
     &              **(-(1.d0/eta)-dble(mij(g(k),ssg(k,g(k))))) 
c               print*,"** prod1 **",prod1(ig),ig,issg,k
               goto 1412
            endif
 1413    continue
 1412 continue
c      stop
             func12=prod1(ig)    
c             print*,"** prod1 **",prod1(ig),frail,eta,alpha,ig
c             stop
             return
             end
c==================================================================

      SUBROUTINE qgauss2(a,b,ss) ! avec troncature
      DOUBLE PRECISION a,b,ss
      double precision func22
      double precision ::auxfunc2a,auxfunc2b

c*****auxig
      integer :: auxig
      common /auxig/auxig

c      Returns as ss the integral of the function func between a and b, by ten-point Gauss-Legendre integration: the function is evaluated exactly ten times at interior points in the range of integration.
c func2 est l intégrant
      INTEGER j
      
      DOUBLE PRECISION dx,xm,xr,w(5),x(5)   !     The abscissas and weights.
      SAVE w,x
      DATA w/.2955242247,.2692667193,.2190863625,.1494513491,
     &     .0666713443/
      DATA x/.1488743389,.4333953941,.6794095682,.8650633666,
     &     .9739065285/
      xm=5.d-1*(b+a)
      xr=5.d-1*(b-a)
      ss=0.d0 ! Will be twice the average value of the function,since the ten
      do 11 j=1,5  !weights (five numbers above each used twice) sum to 2.
      dx=xr*x(j)
      auxfunc2a=func22(xm+dx)
      auxfunc2b=func22(xm-dx)
      ss = ss+w(j)*( auxfunc2a+ auxfunc2b) 
 11   continue
      ss=xr*ss            !  Scale the answer to the range of integration.
      return
      END
c================================================

c================================================
c================================================
      double precision  function func22(frail) ! calcul de l integrant, pour le calcul d intregrale avec troncature
! calcul de l integrant, pour un effet aleatoire donné frail et un groupe donne auxig (cf funcpa)
 
c*************** SAME **********************************************
      integer , parameter ::npmax=50,nsujetmax=20000,nvarmax=50
      integer , parameter ::ngmax=2500,nssgmax=5000,nboumax=2000
      integer , parameter ::ndatemax=30000
      integer , parameter ::nptsmax=1500
c******************************************************************


      integer :: ig,issg,k

      double precision :: frail
c      double precision :: func2
      double precision , dimension(ngmax) :: prod2

c*****auxig
      integer :: auxig
      common /auxig/auxig
c*****ut1ut2
      double precision , dimension(ndatemax) :: dut1,dut2
      double precision , dimension(0:ndatemax) :: ut1,ut2
      common /ut1ut2/dut1,dut2,ut1,ut2
c*****groupe
      integer , dimension(nsujetmax) :: g
      integer , dimension(nsujetmax,ngmax) :: ssg
! numero du ssgpe pour un sujet et un gpe donné
      integer , dimension(ngmax) :: nig
      common /gpe/g,nig,ssg
c*****dace4
      integer , dimension(nsujetmax) :: stra
      common /dace4/stra
c*******alphaeta
      double precision :: alpha,eta
      common /alphaeta/alpha,eta
c*****mij
      integer , dimension(ngmax) :: mi ! nb de dc dans gpe i
      integer  , dimension(ngmax,nssgmax) ::  mij ! nb de dc dans ssgpe ij
      common /mij/mij,mi
c****** indicateur de troncature
      integer :: indictronq ! =0 si donnees non tronquées
      common /troncature/indictronq
c*****contrib
      integer :: ngexact,nssgexact    !nb EXACT de gpes et ss gpes au total
      integer :: nssgbyg  !nb de ss gpes par groupe
      common /contrib/ngexact,nssgexact,nssgbyg 
c*****dace2
      double precision , dimension(nsujetmax) ::t0,t1
      double precision :: ncsrl
      integer , dimension(nsujetmax) ::c
      integer, dimension(nsujetmax) :: nt0,nt1
      integer :: nsujet,nsim,nva,ndate,nst
      common /dace2/t0,t1,ncsrl,c,nt0,nt1,nsujet,nsim,nva,ndate,nst
c******  aux1 aux2
      double precision , dimension(ngmax,nssgmax) ::  aux1,aux2
      common /risqcumul/aux1,aux2
            
c============================================================
      ig=auxig
c      prod2(ig)=1.d0

c     initialisation de prod1 et prod2 par le numerateur de l integrant
      PROD2(IG)=   
     &              (DEXP(DBLE(-frail/ALPHA)))
     &             *(DBLE(frail)**(1.d0/alpha-1.d0))

                  do 1412 issg=1,nssgbyg
                     do 1413 k=1,nsujet
                        if((g(k).eq.ig).and.(ssg(k,g(k)).eq.issg))then
                           if(indictronq.eq.1)then
c                              if(eta.gt.1.d-2)then 
                                 prod2(ig)=prod2(ig)*((1.d0
     &       +(eta*dble(frail)*aux2(g(k),ssg(k,g(k)))))**(-1.d0/eta)) 
c                              endif
                           endif
                           goto 1412
                        endif
 1413                continue
 1412             continue


      func22= prod2(ig)
      return
      end

c=====================
c====== END ==========
c=====================
