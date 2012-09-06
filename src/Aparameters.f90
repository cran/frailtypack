	module tailles
	integer,save:: npmax,NSUJETMAX,nvarmax
	integer,save:: ngmax                    !AD:,maxiter
	integer,save:: ndatemax,ndatemaxdc,nzmax
	integer,save::nssgbyg,nssgmax
	end module tailles

	module parameters
		double precision,save::epsa,epsb,epsd
		integer,save::maxiter
	end module parameters
	

		
	module comon
	implicit none
!*****************************************************************
!*****dace1 
	double precision,dimension(:),allocatable,save::date,datedc
	double precision,dimension(:),allocatable,save::zi

!*****dace2
	double precision,dimension(:),allocatable,save::t0dc,t1dc
	double precision,dimension(:),allocatable,save::t0,t1,t2,t3,tU
	integer,dimension(:),allocatable,save:: c, cdc
	integer,dimension(:),allocatable,save:: nt0,nt1,ntU
	integer,dimension(:),allocatable,save:: nt0dc,nt1dc
	integer,save::nsujet,nva,nva1,nva2,ndate,ndatedc,nst
!*****dace4
	integer,dimension(:),allocatable,save::stra
!*****ve1
	double precision,dimension(:,:),allocatable,save::ve
	double precision,dimension(:,:),allocatable,save::vedc
!*****dace3
	double precision,save::pe
	integer,save::effet,nz1,nz2
!*****dace7
	double precision,dimension(:,:),allocatable,save::I_hess,H_hess
	double precision,dimension(:,:),allocatable,save::Hspl_hess!(npmax,npmax)
	double precision,dimension(:,:),allocatable,save::PEN_deri!(npmax,1)
	double precision,dimension(:,:),allocatable,save::hess!(npmax,npmax)
!*****contrib
	integer,save::ng          !nb de gpes     
!*****groupe
	integer,dimension(:),allocatable,save::g!(nsujetmax)
	integer,dimension(:),allocatable,save::nig!(ngmax)  ! nb d events recurrents par sujet
!*****mem1
	double precision,dimension(:),allocatable,save::mm3,mm2!(ndatemax)
	double precision,dimension(:),allocatable,save::mm1,mm!(ndatemax)
!*****mem2
	double precision,dimension(:),allocatable,save::im3,im2,im1,im
!*****pen1
	double precision,dimension(:),allocatable,save::m3m3,m2m2,m1m1,mmm,m3m2
!*****pen2
	double precision,dimension(:),allocatable,save:: m3m1,m3m,m2m1,m2m,m1m,mi
     
!************************************************************ 
!AD: add for death
	double precision, dimension(:),allocatable,save::mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc  
!AD:end
!************************************************************
! %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
	integer,save::AG
! %%%%%%%%%%%%% indic ALPHA %%%%%%%%%%%%%%%%%%%%%%%%% 
	integer,save::indic_ALPHA
!****  theta/alpha
	double precision,save::theta,alpha,eta !en exposant pour la frailty deces 
!****** indicateur de troncature
	integer,save:: indictronq,indictronqdc ! =0 si donnees non tronquées reellement
!*****auxig
	integer,save :: auxig
!******  aux1 aux2
	double precision,dimension(:),allocatable,save::res1,res3,aux1,aux2
	double precision,save::resnonpen
	double precision,dimension(:),allocatable::kkapa
!****** Type du modele
	integer,save::model  
	double precision,dimension(:),allocatable::vvv 
!cpm
	double precision ::cens
	integer,save:: nbrecu,nbdeces,nbintervR,nbintervDC
	integer,save::indic_eta
!	double precision,save::eta !en exposant pour la frailty deces 
	double precision,dimension(:),allocatable,save::res4
	integer,save::typeof,typeof2
        double precision,dimension(:),allocatable,save::ttt,tttdc
        double precision,dimension(:),allocatable,save::betacoef
!Weib
	double precision,save::etaR,etaD,betaR,betaD
	integer,save::indic_tronc
!censure par intervalle
	integer,dimension(:),allocatable,save::d
	integer,save::dmax
	integer::intcens
	
	end module comon
	


	module additiv	
	implicit none	
		integer,save::correlini,correl		
!*****contrib
		integer,save::ngexact       !nb EXACT de gpes   	
!*****mij
		integer,dimension(:),allocatable,save::mid ! nb de dc dans gpe i  
!******indicateur du nb de parametres
		integer,save::nbpara	
		double precision,dimension(:,:),allocatable,save::ve2	
!*****inversion
		integer,save::indic_sousv,sousm	
!*******sigma2tau2rho
		double precision,save::sigma2,tau2,rho,cov  
!*****ut1ut2
		double precision,dimension(:),allocatable,save::dut1,dut2
		double precision,dimension(:),allocatable,save::ut1,ut2	
!**** betaaux
		double precision,dimension(:),allocatable,save::betaaux
!*****invD
		double precision,dimension(:,:),allocatable,save::invD	
		double precision,dimension(:),allocatable,save::aux1,aux2
		double precision,dimension(:,:),allocatable,save::Xbeta
		
	end module additiv	


	module residusM
		double precision,dimension(:),allocatable,save::Residus &
		,varResidus,cumulhaz,vecuiRes,post_esp,post_SD,som_Xbeta
		double precision,dimension(:),allocatable,save::ResidusRec,&
		Residusdc,Rrec,Nrec,Rdc,Ndc,vecviRes,RisqCumul
		double precision,save::cares,cbres,ddres
		double precision,dimension(:),allocatable,save:: vres
		integer , save :: ierres,nires,istopres,effetres,indg
		double precision,save::rlres,varuiR,moyuiR,varviR,moyviR,corruiviR
		double precision,dimension(:),allocatable::vuu,b_temp
		integer,save::indic_cumul
		integer,dimension(:),allocatable::n_ssgbygrp
		double precision,dimension(:,:),allocatable,save::cumulhaz1,invsigma
		double precision,save::detSigma
		
		
	end module residusM 
		