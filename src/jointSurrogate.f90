!================================================================
! Remarque tres importante: Ne pas nommer les subroutines fortran
! a appeler dans R en usant des majuscules
!================================================================
subroutine jointsurrogate(nsujet1,ng,ntrials1,maxiter,nst,nparamfrail,indice_a_estime,param_risque_base,nbrevar,&
                          filtre0,donnees,death,p,prop_i,n_sim1,EPS2,kappa0,vect_kappa,logNormal,nsim_node,Param_kendall_boot,&
                          vrai_val_init,param_init,revision_echelle,random_generator0,sujet_equi,prop_trait,paramSimul,&
                          autreParamSim,fichier_kendall,fichier_R2, param_estimes,sizeVect,b, H_hessOut,HIHOut,&
                          resOut,LCV,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out,ni,ier,istop,ziOut, affiche_itter,Varcov,&
                          dataHessian,dataHessianIH,datab,vbetast,vbetastinit,param_res_rt,res_rt,knotsurro)

    ! programme principal permettant le traitement des donnees et l'appel du joint_surogate pour l'estimation des parametres

    use sortie
    use Autres_fonctions
    use double_precision
    use var_surrogate, only: graine,aleatoire,nbre_sim,nbre_itter_PGH,nb_procs,random_generator,&
                             affiche_itteration,copula_function,mediation
    use Autres_fonctions, only:pos_proc_domaine
    use var_mediation,only:nmc,nmcboot,nboot,ntimes,nknots,splines_ord,&
                              nsplines!,zmed,zdcmed
    use natural_effects,only:compute_rt

    
    !use mpi ! module pour l'environnement MPI
    !$ use OMP_LIB

    implicit none

    ! =======debut declaration des variables================

    ! =====Parametres prises en entree de la subroutine=====
    integer, dimension(3),intent(in)::nbrevar
    integer,dimension(13), intent(inout)::nsim_node
    integer,intent(in)::nsujet1,ng,ntrials1,nst,maxiter,nparamfrail,n_sim1,logNormal,vrai_val_init,random_generator0,sujet_equi,&
                        affiche_itter
    integer::sujet_equi2
    integer,dimension(5),intent(in)::indice_a_estime
    integer,dimension(5),intent(in):: param_risque_base
    integer,dimension(3),intent(in):: Param_kendall_boot
    integer,dimension(16),intent(in):: autreParamSim ! indique si l'on estime (1=oui, 0=non) estime ou pas zeta(1), covST(2), alpha(3), gammaST(4). indique en (5) si on prend en compte l'heterogeneite sur le risque de base
    integer,dimension(nbrevar(3),2), intent(in)::filtre0
    double precision,dimension(nsujet1,5+nbrevar(1)), intent(in):: donnees
    double precision,dimension(ng,5+nbrevar(2)), intent(in):: death
    double precision,dimension(2), intent(in):: kappa0
    double precision,dimension(9), intent(in):: param_init
    double precision,dimension(23+nbrevar(1) + nbrevar(2)-2), intent(in):: paramSimul
    !double precision,dimension(:), intent(in):: paramSimul
    double precision,dimension(3), intent(inout)::EPS2
    character(len=30),dimension(5)::NomFichier
    double precision, intent(in)::prop_trait,revision_echelle
    double precision::prop_trait2
    integer, dimension(5), intent(in)::sizeVect
    double precision, dimension(ntrials1), intent(in)::p,prop_i
    double precision,dimension(n_sim1,2), intent(in):: vect_kappa
    double precision, dimension(nbrevar(3),2), intent(in)::vbetast,vbetastinit


    ! ! =====Parametres fournies en sortie par la subroutine=====
    integer, intent(out):: ni, istop, ier
    double precision,dimension(n_sim1,3), intent(out):: fichier_kendall,fichier_R2
    double precision,dimension(n_sim1,nsim_node(13)),intent(out):: param_estimes
    !double precision,dimension(:,:),intent(inout):: param_estimes
    double precision,dimension(sizeVect(1)), intent(out)::b
    double precision,dimension(2), intent(out)::LCV
    double precision,dimension(sizeVect(2)), intent(out)::x1Out
    double precision,dimension(sizeVect(3)), intent(out)::x2Out
    double precision,dimension(sizeVect(1),sizeVect(1)), intent(out)::H_hessOut,HIHOut ! H_hessOut = matrice hesienne (des variance-covariance), HIHOut= matrice hessienne corrigee
    double precision,dimension(sizeVect(1)*n_sim1,sizeVect(1)), intent(out)::dataHessian, dataHessianIH ! sauvegarde des matrices hessiennes des differentes simulations
    double precision,dimension(n_sim1,sizeVect(1)), intent(out)::datab ! sauvegarde des vecteurs de parametres de toutes les simulations
    double precision,dimension(sizeVect(2),3), intent(out)::lamOut
    double precision,dimension(sizeVect(3),3), intent(out)::lam2Out
    double precision,dimension(sizeVect(4),3), intent(out)::suOut
    double precision,dimension(sizeVect(5),3), intent(out)::su2Out
    double precision,dimension(param_risque_base(5)+6), intent(out)::ziOut
    double precision,dimension(sizeVect(4)), intent(out)::xSu1
    double precision,dimension(sizeVect(5)), intent(out)::xSu2
    double precision, intent(out)::resOut
    double precision,dimension(3,3), intent(out):: Varcov ! pour la matrice de variance covariance de (sigma_S,sigma_ST_,sigma_T) par la delta-metode

    ! =====Autres variables utilisees dans la subroutine
    !character(len=30)::donnees
    character(len=10), dimension(nbrevar(3))::nomvarl
    character(len=30)::dateamj,zone,heure1,heure2,param_estime, param_empirique,param_empirique_NC,tableau_rejet
    integer::i,j,effet,ver,nva1,nva2,nva,ag,nz, cpt,cpt_dc,noVar1,noVar2,k,typeJoint,np !ii,jj,ncur
    double precision::ax1,ax2,tp1,tp2 !tempon
    integer, dimension(:),allocatable::vdeces,vsurrogate !contient les dates devenement: deces et progression
    character(len=20),dimension(:),allocatable::nomvart,nomvar2t,nomvar,nomvar2
    double precision,dimension(:),allocatable::tt0dc,tt1dc
    integer,dimension(:),allocatable::icdc,groupe,pourtrial,ic,trials,nigs,cdcs,nigts,cdcts
    integer,dimension(:,:),allocatable::nig_Ts,cdc_Ts
    double precision,dimension(:,:),allocatable::vaxdc,vaxdct
    double precision,dimension(:),allocatable::tt0,tt1,ttU
    double precision,dimension(:,:),allocatable::vax,vaxt
    double precision,dimension(2)::k0,k0_save,k01_save,ckappa
    double precision,dimension(3)::EPS
    integer,dimension(8)::values
    integer, dimension(0:1)::randomisation,deces,surrogate
    double precision::bi,bs,wres !wald
    character(len=30)::kapa !aaa !les fichiers de sortie
    integer,dimension(:),allocatable::filtre,filtre2
    !cpm
    integer::mt11,mt12,mt1,mt2,n_sim,ntrials,nsujet
    double precision,dimension(2)::shape_weib,scale_weib
    integer::typeof,nbintervR,nbintervDC,equidistant !nbrecu,nbdeces
    !double precision::Xgamma
    !predictor
    double precision,dimension(:,:),allocatable::MartinGales,v_chap_kendall,v_chap_R2,theta_chap_kendall,theta_chap_R2, &
                                                theta_chap_copula, v_chap_copula
    double precision,dimension(:),allocatable::linearpred,vect_kendall_tau,t_chap_kendall,t_chap_R2,vect_R2
    double precision,dimension(:),allocatable::linearpreddc
    double precision,dimension(:),allocatable::time
    double precision,dimension(:),allocatable::timedc,vbetas,vbetat, vbetas_intit, vbetat_intit
    integer,dimension(4)::mtaille
    integer,dimension(3)::paratps
    double precision,dimension(4)::paraweib
    double precision,dimension(3)::descripSurr,descripDeces
    double precision,dimension(:,:),allocatable:: paGH,matrice_generation ! parametre pour l'adaptative: en ligne les individus, en colone on a respectivement: les ui_cham,
    !racine carree du determinant de l'inverse de la cholesky,variance des ui_chap,les covariances estimees des fragilites pour chaque individu, sachant que la matrice de variances covariance est bien la cholesky
    !parametres de simulation
    integer::n_col,mode_cens,n_essai,n_obs,weib,frailty_cor,affiche_stat,s_i,indice_eta,indice_theta&
                ,rangparam,rangparam2,nbre_rejet,ind_temp,seed_,une_donnee,gener_only,kapa_use,&
                i_min,i_max,i_min_t,i_max_t,ind_rech,Rech_kappa,incre_kappa,statut_kappa,statut_kappa1&
                ,indice_varS,indice_varT,indice_covST,rangparam_sigs,rangparam_sigt,rangparam_sigst,&
                np_save,control_kappa,ind_premier_kappa,control,control2,frailt_base,&
                indice_gamma,indice_alpha,rangparam_gamma,nbre_rejet_0,indice_gamma_st,indice_theta_t,&
                indice_theta_st,indice_gamma_t,rangparam_thetat,rangparam_thetast,rangparam_gammat,&
                rangparam_gammast,rangparam_alpha,decoup_simul,incre_decoup,method_int_kendal,N_MC_kendall,&
                param_weibull,donne_reel,indice_seed,npoint1,npoint2,rangparam_eta,nboot_kendal,nparam_kendall,&
                rangparam_theta,erreur_fichier,indicCP,controlgoto,remplnsim,indice_kapa, type_joint_estim,&
                rangparam_copula,pfs


    double precision::theta,eta,betas,alpha,betat,lambdas,nus,lambdat,nut,temps_cens,&
                        cens0,rsqrt,sigma_s,sigma_t,moy_theta,theta_sim,moy_dec,moy_cens,&
                        moy_pros,moy_theta_est,moy_se_theta,moy_eta,moy_se_eta,&
                        bi2,bs2,n_sim_exact,&
                        moy_trt,taux_couverture_theta,taux_couverture_eta,bi_theta,bs_theta,bi_eta,bs_eta,&
                        gamma1,gamma2,theta2,sigmas_sim,sigmat_sim,rho_sim,varS1,varT1,covST1,varS_es,varT_es,covST_es,&
                        moy_sigmas,moy_sigmat,moy_rho,moy_sigmas_est,moy_sigmat_est,moy_sigmast_est,moy_se_sigmas,moy_se_sigmat,&
                        moy_se_sigmast,bi_sigmas,bi_sigmat,bi_sigmast,bs_sigmas,bs_sigmat,bs_sigmast,taux_couverture_sigmas,&
                        taux_couverture_sigmat,taux_couverture_sigmast,sigmast_vrai,sigma_ss_init,sigma_tt_init,sigma_st_init,&
                        theta_init,betas_init,betat_init,moy_ni, gamma_ui,alpha_ui,gamma_sim,moy_gamma,gamma_init,alpha_init,&
                        moy_gamma_est,moy_se_gamma,bi_gamma,bs_gamma,taux_couverture_gamma,moy_ni_0,moy_trt_0,moy_theta_0,&
                        moy_sigmas_0,moy_sigmat_0,moy_rho_0,moy_gamma_0,tab_var_sigma_0,moy_se_theta_0,taux_couverture_theta_0,&
                        moy_gamma_est_0,moy_se_gamma_0,taux_couverture_gamma_0,moy_sigmas_est_0,moy_se_sigmas_0,&
                        taux_couverture_sigmas_0,moy_sigmat_est_0,moy_se_sigmat_0,taux_couverture_sigmat_0,&
                        moy_sigmast_est_0,moy_se_sigmast_0,taux_couverture_sigmast_0,moy_eta_0,moy_se_eta_0,&
                        taux_couverture_eta_0,moy_betaS_0,moy_betaS_se_0,taux_couvertureS_0,moy_betaT_0,&
                        moy_betaT_se_0,taux_couvertureT_0,se_theta_sim,se_sigmas_sim,se_sigmat_sim,se_rho_sim,&
                        se_gamma_sim,& !se_theta_sim_0,se_sigmas_sim_0,se_sigmat_sim_0,se_rho_sim_0,se_gamma_sim_0
                        n_sim_exact_0,moy_theta_est_0,moy_pros_0,moy_dec_0,theta2_t,rsqrt_theta,gamma_uit,rsqrt_gamma_ui,&
                        thetat_init,thetast_init,gammat_init,gammast_init,theta_simt,rho_sim_wij,gamma_simt,rho_sim_ui,&
                        moy_thetat,moy_rho_wij,moy_gammat,moy_rho_ui,thetast_vrai,gammast_vrai,moy_thetat_est,moy_se_thetat,&
                        taux_couverture_thetat,moy_thetast_est,moy_se_thetast,taux_couverture_thetast,moy_gammat_est,moy_se_gammat,&
                        taux_couverture_gammat,moy_gammast_est,moy_se_gammast,taux_couverture_gammast,thetaS1,thetaT1,thetaST1,&
                        thetaT_es,thetaST_es,gammaS1,gammaT1,gammaST1,gammaS_es,gammaT_es,gammaST_es,moy_alpha,moy_se_alpha,&
                        taux_couverture_alpha,se_theta_simt,se_rho_sim_wij,se_gamma_simt,se_rho_sim_ui,se_theta_est,se_eta_est,&
                        se_beta_s,se_beta_t,se_sigmas_est,se_sigmat_est,se_cov_est,se_gamma_est,se_alpha_est,se_thetat_est,&
                        se_cov_est_wij,se_gamma_estt,se_cov_est_ui,se_theta_est_0,se_eta_est_0,se_beta_s_0,se_beta_t_0,&
                        se_sigmat_est_0,se_cov_est_0,se_gamma_est_0,se_alpha_est_0,se_thetat_est_0,se_gamma_estt_0,thetaS_es,&
                        se_cov_est_ui_0,moy_thetat_0,se_theta_simt_0,moy_thetat_est_0,moy_se_thetat_0,taux_couverture_thetat_0,&
                        moy_rho_wij_0,se_rho_sim_wij_0,moy_thetast_est_0,se_cov_est_wij_0,moy_se_thetast_0,&
                        moy_gammat_0,se_gamma_simt_0,moy_gammat_est_0,moy_se_gammat_0,taux_couverture_gammat_0,&
                        moy_rho_ui_0,se_rho_sim_ui_0,moy_gammast_est_0,moy_se_gammast_0,taux_couverture_gammast_0,&
                        moy_alpha_0,moy_se_alpha_0,taux_couverture_alpha_0,R2_trial,se_R2_trial,moy_R2_trial,&
                        taux_couverture_R2_trial,moy_se_R2_trial,moy_bi_R2_trial,moy_bs_R2_trial,moy_kendal_11,tau_kendal_11,&
                        moy_kendal_10,tau_kendal_10,moy_kendal_01,tau_kendal_01,moy_kendal_00,tau_kendal_00,se_kendal_11,&
                        se_kendal_01,se_kendal_00,moy_tau_boots,IC_Inf,IC_sup,zeta_init,moy_R2_boots,IC_Inf_R2,IC_sup_R2,&
                        CP_R2_boot,CP_ktau_boot,se_sigmas_est_0,taux_couverture_thetast_0,se_kendal_10,&
                        bi_R2_trial,bs_R2_trial,thetacopule, thetacopula_init, moy_param_cop, moy_se_param_cop,&
                        bi_param_cop, bs_param_cop, pour_ic, taux_couverture_param_cop,taux_couverture_tauk, vrai_tau_copula

    double precision, dimension(:,:),allocatable::don_simul,don_simulS, don_simultamp,don_simulStamp,don_simulS1,&
                        parametre_empirique, parametre_estimes,parametre_empirique_NC,parametre_estimes_MPI,&
                        parametre_estimes_MPI_T,result_bootstrap
    double precision, dimension(:),allocatable::tab_var_theta,tampon,tampon_all, moy_betaS, moy_betaT,moy_betaS_se, moy_betaT_se,&
                        taux_couvertureS, taux_couvertureT
    double precision, dimension(:,:),allocatable::tab_var_sigma
    integer,parameter ::trt1=1,v_s1=2,v_t1=3,trialref1=4,w_ij1=5,timeS1=6,timeT1=7,&
                      timeC1=8,statusS1=9,statusT1=10,initTime1=11,Patienref1=12,u_i1=13,&
                          w_ijt=14,u_it=15 ! definissent les indices du tableau de donnee simulees
    integer, dimension(:),allocatable::tab_rejet ! pour les rangs des jeux de donnees rejetees
    double precision, dimension(:,:),allocatable::kappa,d_S,d_T !jeu de donnees reelle pour le test
    integer,dimension(:),allocatable::tableEssai,tableNsim ! tableNsim: indique le nombre de simulation a effectuer par chaque processus
    double precision,dimension(:,:),allocatable::donnee_essai,theta_st_2,gamma_st_2,theta_st0_2,gamma_st0_2 !sigma_st_2,sigma_st0_2
    double precision,dimension(2,2)::chol,sigma_st,theta_st,gamma_st,sigma_st0,theta_st0,gamma_st0,Chol_R2,mat_A
    integer, dimension(4)::indice_esti
    integer::nb_processus,rang,n_sim_total,suplement,init_i,max_i,debut_exe,indice_sim_proc, &
            rang_proc, code_print ! je redefini ces indices car les precedentes sont utilisees autrement: cas OpenMP
    !double precision,dimension(10)::t
    double precision,dimension(3,3):: sigmac ! pour la mtrice de variance covariance de Sigma par la delta-metode
    double precision,dimension(3,3):: hb

    integer,dimension(9),intent(in)::param_res_rt
    !param_res_rt(1) : mediation (t/f)
    !            (2) : nknots (innerknots, gamma(s))
    !            (3) : splines_ord for gamma(s)
    !            (4) : ntimes, number of times at which R(t) should be computed
    !            (5) : nmc, monte carlo simulations to be used in R(t) computation
    !            (6) : boot : shall be boostraped sample of R(t) computed?
    !            (7) : nboot : if so, how many bootstrap?
    !            (8) : nmcboot : nmc but for the bootstrap samples
    !            (9) : integ_type : integration method over S :
    !                     1 : simpson's rule with 200 evaluation
    !                     2 : gauss-laguerre with 30 nodes
    !                     3 : adaptive gauss kronrod
    double precision,dimension(2+param_res_rt(2)),intent(out)::knotsurro

    double precision,dimension(param_res_rt(4),4,1+param_res_rt(6)*&
                               param_res_rt(7)),intent(out)::res_rt
    !double precision,dimension(1+param_res_rt(6)*param_res_rt(7),&
    !                           param_res_rt(4),4),intent(inout)::res_rt


    integer::rt_boot,integ_type

    !character(len=1025)::filename
    !double precision::dummy_xx
    !integer::file_id
    !integer::posrecl !rec to write in direct access to binary data file
    covST_es = 0.0d0 
    !dummy if statement to remove warning regarding unused arguments prop_trait & sujet_equi
    if(.false.) then 
        prop_trait2=prop_trait
        sujet_equi2=sujet_equi
    end if 
    if(param_res_rt(1).eq.1) then
      mediation = .TRUE.
      nknots=param_res_rt(2)
      splines_ord=param_res_rt(3)
      ntimes=param_res_rt(4)
      nmc=param_res_rt(5)
      integ_type=param_res_rt(9)
      nsplines=splines_ord+nknots
      if(param_res_rt(6).eq.1) then
        rt_boot = 1
        nboot = param_res_rt(7)
        nmcboot = param_res_rt(8)
      else
        rt_boot = 0
        nboot = 0
        nmcboot = 0
      end if
    else
      mediation = .FALSE.
      nknots = 0
      splines_ord = 0
      ntimes = 0
      nmc = 0
      nboot = 0
      nmcboot = 0
      integ_type=0
      nsplines=0
    end if
    !=====================================================================================
    !*********fin declaration des variables et debut du programme principale**************
    !=====================================================================================

    copula_function = nsim_node(12) ! the copula function, can be 1 for clayton or 2 for Gumbel-Hougaard
    type_joint_estim = nsim_node(8) ! type of estimated model
    ! affectation de certains parametres
    nomvarl(1) = "trt"
    NomFichier(1) = "kappa_valid_crois.txt"
    NomFichier(2) = "Parametre_estime.txt"
    NomFichier(3) = "Parametre_empirique.txt"
    NomFichier(4) = "Parametre_empirique_NC.txt"
    NomFichier(5) = "tab_rejet.txt"
    rang_proc = 0
    rangparam_eta = 0
    rangparam_gammast = 0
    rangparam_thetat = 0
    rangparam_thetast = 0
    rangparam_gammat = 0
    statut_kappa = 0

    pour_ic = 0.d0
    R2_trial = 0.d0
    se_R2_trial = 0.d0
    tau_kendal_10 = 0.d0
    tau_kendal_01 = 0.d0
    tau_kendal_00 = 0.d0
    tau_kendal_11 = 0.d0
    thetaS1 = 0.d0
    thetaT1 = 0.d0
    thetaST1 = 0.d0
    thetaT_es = 0.d0
    thetaST_es = 0.d0
    thetat_init = 0.d0
    thetast_init = 0.d0
    varS_es = 0.d0
    varT_es = 0.d0
    gammaS1 = 0.d0
    gammaT1 = 0.d0
    gammaST1 = 0.d0
    gammaS_es = 0.d0
    gammaT_es = 0.d0
    gammaST_es = 0.d0
    gammat_init = 0.d0
    gammast_init = 0.d0
    moy_tau_boots = 0.d0
    moy_R2_boots = 0.d0



    affiche_itteration = affiche_itter
    n_sim = n_sim1 ! nombre de simulations
    ntrials = ntrials1
    nsujet = nsujet1
    typeof = param_risque_base(1)    !type de function de risque  0:Splines,  1:Cpm  2:weib
    nbintervR = param_risque_base(2) !Nombre intervalle surrogate
    nbintervDC = param_risque_base(3) !Nombre intervalle deces
    equidistant = param_risque_base(4) !cpm (1:equidistant, 0:percentile)')
    nz = param_risque_base(5) ! nombre de noeud pour les spline

    ! indisue si l'on estime (1=oui, 0=non) estime ou pas zeta(1), covST(2), alpha(3), gammaST(4). indique en (5) si on prend en compte l'heterogeneite sur le risque de base
    indice_eta=indice_a_estime(1)    !dit si l'on estime eta (1) ou pas (0)
    indice_covST=indice_a_estime(2) !dit si l'on estime la covariance des frailties essais ou pas
    indice_alpha=indice_a_estime(3)    !dit si l'on estime alpha (1) ou non(0)
    indice_gamma_st=indice_a_estime(4) !dit si l'on estime sigma_us_ut (1) ou non(0)
    frailt_base=indice_a_estime(5) !dit si l'on prend en compte l'heterogeneite sur le risque de base aussi bien dans la generation des donnes que dans l'estimation(1) ou non (0)
    !call intpr("I'm there scl 13:", -1, code_print, 1)
    ver=nbrevar(3) ! nombre total de variables explicatives
    allocate(vbetas(ver),vbetat(ver),vbetas_intit(ver), vbetat_intit(ver))

    AG=0 ! andersen-gill approach(1=oui)

    ! Gestion des kappas pour le jeux de donnees reelles
    ax1 = kappa0(1)    ! pour surrogate
    ax2 = kappa0(2)    ! pour true
    kapa= NomFichier(1)

    ! Parametres associes au taux de kendall et au bootstrap
    method_int_kendal = Param_kendall_boot(1)
    N_MC_kendall = Param_kendall_boot(2)
    nboot_kendal = Param_kendall_boot(3)

    ! Autres noms de fichiers:
    param_estime = NomFichier(2)
    param_empirique = NomFichier(3)
    param_empirique_NC = NomFichier(4)
    tableau_rejet = NomFichier(5)

    ! Parametres initiaux
    theta_init = param_init(1) ! if we are estimating the joint surrogate model
    thetacopula_init = param_init(1) ! if we are estimating copula modele
    sigma_ss_init = param_init(2)
    sigma_tt_init = param_init(3)
    sigma_st_init = param_init(4)
    gamma_init = param_init(5)
    alpha_init = param_init(6)
    zeta_init = param_init(7)
    betas_init = param_init(8)
    betat_init = param_init(9)
    random_generator = random_generator0

    ! parametres de simulation
    gamma1 = paramSimul(1)
    gamma2 = paramSimul(2)
    theta2 = paramSimul(3)
    eta = paramSimul(4)
    gamma_ui = paramSimul(5)
    alpha_ui = paramSimul(6)
    theta2_t = paramSimul(7)
    rsqrt_theta = paramSimul(8)
    gamma_uit = paramSimul(9)
    rsqrt_gamma_ui = paramSimul(10)
    betas = paramSimul(11)
    betat = paramSimul(12)
    lambdas = paramSimul(13)
    nus = paramSimul(14)
    lambdat = paramSimul(15)
    nut = paramSimul(16)
    mode_cens = int(paramSimul(17))
    temps_cens = paramSimul(18)
    cens0 = paramSimul(19)
    rsqrt = paramSimul(20)
    sigma_s = paramSimul(21)
    sigma_t = paramSimul(22)
    thetacopule = paramSimul(23)
    ! les les parametres de simulation pour les autres covariables se trouvent a la fin du tableau paramSimul

    if(nsim_node(11) == 3) then ! joint frailty copula, remplissage des vecteurs des variables explicatives
        do i = 1, size(vbetast,1)
            vbetas(i) = vbetast(i,1) ! beta_s
            vbetat(i) = vbetast(i,2) ! beta_t
            vbetas_intit(i) = vbetastinit(i,1) ! beta_s
            vbetat_intit(i) = vbetastinit(i,2) ! beta_t
        enddo
    endif

    if(nsim_node(11) == 1) then ! joint surrogate, remplissage des vecteurs des variables explicatives
        vbetas(1) = betas ! beta_s
        vbetat(1) = betat ! beta_t
    endif

    ! autres parametres de simulation
    weib = autreParamSim(1)
    param_weibull = autreParamSim(2)
    frailty_cor = autreParamSim(3)
    affiche_stat = autreParamSim(4)
    seed_ = autreParamSim(5)
    une_donnee = autreParamSim(6)
    donne_reel = autreParamSim(7)
    gener_only = autreParamSim(8)
    kapa_use = autreParamSim(9)
    decoup_simul = autreParamSim(10)
    aleatoire = autreParamSim(11)
    nbre_sim = autreParamSim(12)
    graine = autreParamSim(13)
    ckappa(1) =dble(autreParamSim(14))
    ckappa(2) =dble(autreParamSim(15))
    pfs = autreParamSim(16)
    ! call intpr("avant appel joint:pfs", -1,pfs, 1)

    np=sizeVect(1)
    call date_and_time(dateamj,heure1,zone,values) ! pour la date de debut du programme

    rang=0
    nb_processus=1
    nb_procs=1
    controlgoto=0

    typeJoint=8 ! 8 pour les models conjoint surrogate (sans reccurence) et true endpoint


    allocate(groupe(nsujet),ic(nsujet),tt0(nsujet),tt1(nsujet),ttU(nsujet),linearpred(nsujet),pourtrial(nsujet))

    allocate(tt0dc(ng),tt1dc(ng),icdc(ng),trials(ntrials),nigs(ntrials),cdcs(ntrials),nigts(ntrials),cdcts(ntrials),&
    nig_Ts(ntrials,2),cdc_Ts(ntrials,2))
    allocate(MartinGales(ng,4),linearpreddc(ng))

    if(indice_eta==0 .and. (nparamfrail==2 .or. nparamfrail==5)) then
        !!print*,"Attention si indice_eta=0, nparamfrail==1 ou >2"
        !stop
    endif

    indice_esti(1)=indice_eta
    indice_esti(2)=frailt_base
    indice_esti(3)=indice_alpha
    indice_esti(4)=indice_gamma_st

    !indicateur de presence deffet aleatoire (0=Non, 1=Oui)
    if(nparamfrail.eq.0) then
        effet=0
    else
        effet=1
        indice_theta=1
        indice_varS=1
        indice_varT=1
        !indice_covST=1
        if(frailt_base==1) then
            indice_gamma=1
            !indice_alpha=1
        endif

        if(nsim_node(8)==2) then
            indice_gamma=1
            indice_theta_t=1
            indice_theta_st=1
            indice_gamma_t=1
        endif
    end if

    if(type_joint_estim == 3) then ! joint frailty-copula model
        indice_theta = 0
        indice_theta_t = 0
        indice_theta_st = 0
    endif
    mt11 = 100 !scl  A quoi sert ce parametre?
    mt12 = 100

    if(typeof == 1) then
        mt1 = 3*nbintervR    !scl: pourquoi *3
        mt2 = 3*nbintervDC
    else
        if(typeof == 2) then
            ! read(2,*)nbintervR,nbintervDC
        end if
        mt1 = 100
        mt2 = 100
    end if

    allocate(time(nbintervR+1),timedc(nbintervDC+1))


    nva1 = 0         ! nb de var expli pour donnees recurrentes
    nva2 = 0         ! nb de var expli pour deces

    allocate(filtre(ver),filtre2(ver),nomvart(ver),nomvar2t(ver))!scl filtre indique si la variable est prise en compte pour les reccures et filtre2 dit si elle est prise en compte pour les deces
    allocate(vaxt(nsujet,ver),vaxdct(ng,ver)) ! matrice de toutes les variables explicatives

    ! indicatrice de prise en compte de variable
    do i=1,ver
        filtre(i)=filtre0(1,i)
        filtre2(i)=filtre0(2,i)
    enddo


    if(ver.gt.0)then ! aumoins une variable explicative comme c'est le cas
        do j=1,ver
           ! read(2,*)nomvarl,filtre(j),filtre2(j) ! nom de la variable + indicateur d'appartenance aux deux fichiers surrogate et true
            !if(rang_proc==0) !write(*,*)"       ",nomvarl(j)," ",filtre(j)," ",filtre2(j)
            nva1 = nva1 + filtre(j) ! adjustment for recurrent events
            nva2 = nva2 + filtre2(j) ! adjustment for survival
            if(filtre(j).eq.1)then
                nomvart(nva1) = nomvarl(j)
            endif
            if(filtre2(j).eq.1)then
                nomvar2t(nva2) = nomvarl(j)
            endif
        end do
    endif

    !On active ou pas les filtres
    !------> filtre 1
    if (nva1 .ne. 0) then
        noVar1 = 0 !
      !  if(rang_proc==0) !write(*,*)'****** Présence de variables explicatives pour surrogate   *******'
    else
        noVar1 = 1
    end if

    !------> filtre 2
    if (nva2 .ne. 0) then
        noVar2 = 0
       ! if(rang_proc==0) !write(*,*)'********* Présence de variables explicatives pour décès **********'
    else
        noVar2 = 1
    end if


    allocate(vax(nsujet,nva1)) ! matrice des variables explicatives presentes dans le jeux de donnees surrogate
    allocate(nomvar(nva1),nomvar2(nva2)) !Nnom des variables explicatives presentes dans le jeux de donnees surrogate et deces

    if(ver.gt.0)then
        nomvar=nomvart(1:nva1)
        nomvar2=nomvar2t(1:nva2)
    end if

    nva = nva1+nva2

    allocate(vaxdc(ng,nva2)) ! matrice des variables explicatives presentes dans le jeux de donnees deces
    allocate(moy_betaS(nva1), moy_betaT(nva2),moy_betaS_se(nva1), moy_betaT_se(nva2),&
            taux_couvertureS(nva1), taux_couvertureT(nva2))


    npoint1=nsim_node(2) ! nombre de point de quadrature a privilegier initialement
    npoint2 = nsim_node(9) ! nombre de point de quadrature a utiliser en cas de non convergence de prefenrence 7 ou 9 pour la pseudo adaptative et 32 pour la non adaptative
    ! read(2,*)nsim_node(3) ! doit-on faire de l'adaptative(1) ou de la non-adaptative(0)
    nbre_itter_PGH = nsim_node(10) !nombre d'itteration aubout desquelles reestimer les effects aleatoires a posteriori pour la pseude adaptative. si 0 pas de resestimation

    np_save=nsim_node(2)

    nparam_kendall=4 ! on a 4 parametres qui rentrent dans le calcul du tau de kendall: theta, alpha, gamma, zeta

    if(method_int_kendal==4) then
        if(indice_alpha==0) nparam_kendall=nparam_kendall-1
        if(indice_eta==0) nparam_kendall=nparam_kendall-1
    endif
    if(method_int_kendal==5) then
        nparam_kendall=2 ! dans ce cas on a fixe u_i =0, ce qui annule gamma et alpha
        if(indice_eta==0) nparam_kendall=nparam_kendall-1
    endif

    allocate(vect_kendall_tau(nboot_kendal),v_chap_kendall(nparam_kendall,nparam_kendall),vect_R2(nboot_kendal),&
            theta_chap_kendall(1,nparam_kendall),t_chap_kendall(nparam_kendall),v_chap_R2(3,3),t_chap_R2(3),&
            theta_chap_R2(1,3),result_bootstrap(n_sim,6), theta_chap_copula(1,1), v_chap_copula(1,1))

    indice_kapa = 1
    n_essai=ntrials
    n_obs=nsujet
    ! le jeu de donnee doit contenir 13 ou 17 colonnes
    !read(5,*)n_col    ! je ne lis plus cette variable
    if(nsim_node(8)==2)then
        n_col=15 + nbrevar(3)-1 ! j'ajoute le surplus des covariables, -1 pour le traitement qui est deja pris en compte
    else
        n_col=13 + nbrevar(3)-1 ! j'ajoute le surplus des covariables
    endif
    alpha = eta    ! alpha associe a u_i chez les deces

    !generation des donnees par joint failty-copula
    if(nsim_node(11)==3) then 
        allocate(don_simultamp(n_obs,n_col),don_simulStamp(n_obs,n_col))
    else 
        allocate(don_simultamp(1,1),don_simulStamp(1,1))
    end if 
    allocate(don_simul(n_obs,n_col),don_simulS1(n_obs,n_col))

    if(nsim_node(8)==2)then
        allocate(donnee_essai(n_essai,5))
    else
        allocate(donnee_essai(n_essai,4))
    endif

    don_simul=0.d0
    sigmast_vrai=rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
    thetast_vrai=rsqrt_theta*dsqrt(theta2)*dsqrt(theta2_t)
    gammast_vrai=rsqrt_gamma_ui*dsqrt(gamma_ui)*dsqrt(gamma_uit)

    if(nsim_node(11) == 3) then ! si joint frailty copula, alors on ajout les covariable aux jeux de donnees
      allocate(d_S(nsujet*n_sim,6 +size(vbetast,1)-1),d_T(ng*n_sim,6+size(vbetast,1)-1))
    else
      allocate(d_S(nsujet*n_sim,6),d_T(ng*n_sim,6))
    endif

    if(une_donnee==1) then
        ! on recupere le jeu de donnees reelles
        d_S=donnees
        d_T=death
    else
        ! on recupere le jeu de donnees reelles
    endif

    moy_theta=0.d0
    moy_theta=0.d0
    moy_sigmas=0.d0
    moy_sigmat=0.d0
    moy_rho=0.d0
    moy_gamma=0.d0
    moy_trt=0.d0
    moy_dec=0.d0
    moy_pros=0.d0
    moy_cens=0.d0
    moy_theta_est=0.d0
    moy_se_theta=0.d0
    moy_sigmas_est=0.d0
    moy_sigmat_est=0.d0
    moy_sigmast_est=0.d0
    moy_eta=0.d0
    moy_se_eta=0.d0
    moy_se_sigmas=0.d0
    moy_se_sigmat=0.d0
    moy_se_sigmast=0.d0
    moy_betaS=0.d0
    moy_betaS_se=0.d0
    moy_betaT=0.d0
    moy_betaT_se=0.d0
    moy_thetat=0.d0
    moy_rho_wij=0.d0
    moy_gammat=0.d0
    moy_rho_ui=0.d0
    taux_couvertureS=0.d0
    taux_couvertureT=0.d0
    taux_couverture_theta=0.d0
    taux_couverture_eta=0.d0
    taux_couverture_sigmas=0.d0
    taux_couverture_sigmat=0.d0
    taux_couverture_sigmast=0.d0
    nbre_rejet=0.d0
    moy_ni=0
    moy_gamma_est=0.d0
    moy_se_gamma=0.d0
    taux_couverture_gamma=0.d0
    moy_thetat_est=0.d0
    moy_se_thetat=0.d0
    taux_couverture_thetat=0.d0
    moy_thetast_est=0.d0
    moy_se_thetast=0.d0
    taux_couverture_thetast=0.d0
    moy_gammat_est=0.d0
    moy_se_gammat=0.d0
    taux_couverture_gammat=0.d0
    moy_gammast_est=0.d0
    moy_se_gammast=0.d0
    taux_couverture_gammast=0.d0
    moy_alpha=0.d0
    moy_se_alpha=0.d0
    taux_couverture_alpha=0.d0
    moy_R2_trial=0.d0
    moy_se_R2_trial=0.d0
    moy_bi_R2_trial=0.d0
    moy_bs_R2_trial=0.d0
    taux_couverture_R2_trial=0.d0
    moy_kendal_11=0.d0
    moy_kendal_10=0.d0
    moy_kendal_01=0.d0
    moy_kendal_00=0.d0
    debut_exe=1  ! pour dire on peut commencer l'execution pour ce processus
    indice_sim_proc=1 ! pour indicer le table des parametres estimes MPI
    indicCP=1 ! pour indicer le vecteur des IC par bootstrap
    moy_param_cop = 0.d0
    moy_se_param_cop = 0.d0
    taux_couverture_param_cop = 0.d0
    taux_couverture_tauk = 0.d0

    allocate(parametre_empirique(n_sim,12),parametre_empirique_NC(n_sim,12)) !parametre_empirique: contient les parametres empirique pour chaque simul (trt,theta,surrrogate,deces,vars,vart,rhost,gamma)
    if(nsim_node(8)==0)then !model conjoint surrogate classique
        allocate(parametre_estimes(n_sim,8)) !parametres estimes: contient les parametres estimes(theta_chap+sd,zeta+sd,beta_s+sd,beta8t+sd)
        !allocate(param_estimes(n_sim,24))
        allocate(parametre_estimes_MPI(n_sim,1),parametre_estimes_MPI_T(n_sim,1))
    else
        if(nsim_node(8)==1)then ! modele avec les effets aleatoires partages
            allocate(parametre_estimes(n_sim,24)) !parametres estimes: contient les parametres estimes(theta_chap+sd,zeta+sd,beta_s+sd,beta8t+sd,sigma_s+sd,sigma_t+sd,sigmast+sd,gamma_ui+sd,alpha_ui+sd, R2 reduit et sd, taux de kendall)
            allocate(parametre_estimes_MPI(n_sim,24))! contient les parametres estimes par chaque processus dans MPI
            allocate(parametre_estimes_MPI_T(n_sim,24)) ! contient tous les parametres, de tous les processus
        else if(nsim_node(8)==2)then ! modele complet avec les effets correles. on a 11 parametre a estimer dans le pire des cas avec leur SE
            allocate(parametre_estimes(n_sim,32)) !parametres estimes: contient les parametres estimes(theta_chap+sd,zeta+sd,beta_s+sd,beta8t+sd,sigma_s+sd,sigma_t+sd,sigmast+sd,gamma_ui+sd,alpha_ui+sd, R2 reduit et sd, taux de kendall)
            allocate(parametre_estimes_MPI(n_sim,32))! contient les parametres estimes par chaque processus dans MPI
            allocate(parametre_estimes_MPI_T(n_sim,32)) ! contient tous les parametres, de tous les processus
        else! joint frailty-copula model
            allocate(parametre_estimes(n_sim,25 + nva -2)) !parametres estimes: contient les parametres estimes(theta_copula+sd,zeta+sd,beta_s+sd,beta8t+sd,sigma_s+sd,sigma_t+sd,sigmast+sd,gamma_ui+sd,alpha_ui+sd, R2 reduit et sd, taux de kendall + sd) + variables explicative supplementaires
            allocate(parametre_estimes_MPI(n_sim,25 + nva -2))! contient les parametres estimes par chaque processus dans MPI
            allocate(parametre_estimes_MPI_T(n_sim,25 + nva -2)) ! contient tous les parametres, de tous les processus
        endif
        !allocate(param_estimes(n_sim,size(parametre_estimes,2)))
    endif



    allocate(tab_rejet(n_sim),tab_var_theta(n_sim),tab_var_sigma(n_sim,8))
    !allocate(Vect_sim_MC(nsim_node(1),1))
    parametre_empirique=0.d0
    parametre_estimes=0.d0
    parametre_estimes_MPI=0.d0
    parametre_estimes_MPI_T=0.d0
    tab_rejet=0
    allocate(kappa(n_sim,2))

    incre_kappa=1 ! pour contenir le nombre de kappa ayant permis la convergence
    statut_kappa1=0
    s_i=1

    n_sim_total=n_sim ! nombre de simulations a effectuer au total
    n_sim=INT(n_sim/nb_processus)
    suplement=n_sim_total-n_sim*nb_processus ! donne le nombre de simulation a partager entre les premiers processus seulement

    ! remplissage du table du nombre de simulation a effectuer par processus
    !!print*,nb_processus,n_sim,suplement
    allocate(tableNsim(nb_processus))
    tableNsim(1:nb_processus)=n_sim
    !!print*,tableNsim
    tableNsim(1:suplement)=n_sim+1 ! tous les essais jusqu'au rang supplement-1 recoivent une tâche supplementaire a realiser
    n_sim=tableNsim(rang+1) ! nombre de simulations a effectuer par le processus courant

    if (rang==0) then
       ! if(rang_proc==0) !print*, "tableNsim=",tableNsim
        init_i=1 ! ce processus commence a la premiere simulation
    else
        init_i=sum(tableNsim(1:rang))+1 ! ce processus commence a la simulation qui respecte son ordre et doit s'arreter au nombre de simultation dont il a le droit d'executer
    endif

    max_i=init_i+tableNsim(rang+1)-1!rang maximale de la simulation a executer (-1 car on a deja incrementer init_i de 1)
    !!call MPI_ABORT(MPI_COMM_WORLD,erreur,code)! on stop tous les programmes appartenant au communicateur code, equivalent de l'instruction stop en sequantiel
    incre_decoup=0
    do while (s_i<=n_sim_total)
        if(une_donnee.eq.1) then ! dans ce cas on utilise un seul jeu de donnees pour les simul, de preference les donnees reel

            ! indice des donnees a utiliser pour la simulation encours
            i_min=s_i*nsujet-nsujet+1
            i_max=s_i*nsujet
            i_min_t=s_i*ng-ng+1
            i_max_t=s_i*ng
            if(donne_reel==1) then
                don_simulS1(:,initTime1)=d_S(i_min:i_max,4)
                don_simulS1(:,timeS1)=d_S(i_min:i_max,5)
                don_simulS1(:,statusS1)=d_S(i_min:i_max,6)
                don_simulS1(:,trialref1)=d_S(i_min:i_max,1)
                don_simulS1(:,Patienref1)=d_S(i_min:i_max,2)
                don_simulS1(:,trt1)=d_S(i_min:i_max,3)
                don_simul(:,initTime1)=d_T(i_min_t:i_max_t,4)
                don_simul(:,timeT1)=d_T(i_min_t:i_max_t,5)
                don_simul(:,statusT1)=d_T(i_min_t:i_max_t,6)
                don_simul(:,trialref1)=d_T(i_min_t:i_max_t,1)
                don_simul(:,Patienref1)=d_T(i_min_t:i_max_t,2)
                don_simul(:,trt1)=d_T(i_min_t:i_max_t,3)
            else ! alors la position des variables n'est plus la meme
                don_simulS1(:,initTime1)=d_S(i_min:i_max,1)
                don_simulS1(:,timeS1)=d_S(i_min:i_max,2)
                don_simulS1(:,statusS1)=d_S(i_min:i_max,3)
                don_simulS1(:,trialref1)=d_S(i_min:i_max,4)
                don_simulS1(:,Patienref1)=d_S(i_min:i_max,5)
                don_simulS1(:,trt1)=d_S(i_min:i_max,6)
                don_simul(:,initTime1)=d_T(i_min_t:i_max_t,1)
                don_simul(:,timeT1)=d_T(i_min_t:i_max_t,2)
                don_simul(:,statusT1)=d_T(i_min_t:i_max_t,3)
                don_simul(:,trialref1)=d_T(i_min_t:i_max_t,4)
                don_simul(:,Patienref1)=d_T(i_min_t:i_max_t,5)
                don_simul(:,trt1)=d_T(i_min_t:i_max_t,6)
            endif

            ! j'ajoute les autres variables a la fin
            do i = 2,nbrevar(3)
                if(filtre(i).eq.1)then
                    don_simulS1(:,size(don_simulS1,2) - nbrevar(3) + i - 1)=d_S(i_min:i_max,size(don_simulS1,2) &
                    - nbrevar(3) + i - 1)
                endif
                if(filtre2(i).eq.1)then
                    don_simul(:,size(don_simul,2)- nbrevar(3) + i - 1)=d_T(i_min_t:i_max_t,size(don_simul,2)- nbrevar(3) + i - 1)
                endif
            enddo
            ! on met à jour le nombre d'essais
            20041 continue

            ! pour la gestion des paquets de simulation, avance dans le fichier des kappas pour se placer au bon endroit
            if(incre_decoup<decoup_simul) then !incre_decoup<decoup_simul c'est pour gerer le cas des simpulations par paquet
                ax1 = vect_kappa(indice_kapa,1)
                ax2 = vect_kappa(indice_kapa,2)
                indice_kapa = indice_kapa +1
                incre_decoup=incre_decoup+1
                goto 20041 ! pour etre sur qu'on n'utilise pas les kappas des jeux de donnee a ne pas considerer dans ce paquet de simulation
            endif


            if((s_i<init_i).or.s_i>max_i) then
                debut_exe=0 ! pour dire le processus ne considere pas ce jeu de donnee
            else
                debut_exe=1 !jeux de donnees a considerer pour le processus courant
            endif
            allocate(tableEssai(ng))

            tableEssai=table_essai(nint(don_simul(:,trialref1)))

            n_essai=0
            do i=1,ng
                if(tableEssai(i).ne.0)then
                    n_essai=n_essai+1
                endif
            enddo
            ntrials=n_essai
            deallocate(tableEssai)
            ind_temp=ng
            theta=theta2
            !stop
        else
            theta=0.d0
            !call intpr("nsim_node(8) =", -1, nsim_node(8), 1)
            if(nsim_node(8)==0) then
                indice_seed=0
                call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation
                10 continue
                call generation_Gamma(don_simul,don_simulS1,ng,n_col,logNormal,affiche_stat,theta,&
                    ng,ver,alpha,cens0,temps_cens,gamma1,gamma2,theta2,lambdas,nus,lambdat,nut,betas,betat)
                indice_seed=indice_seed+1
                if(indice_seed<s_i) goto 10
            else
                indice_seed=0
                call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation
                11 continue
                if(nsim_node(11)==1) then !modele avec effets aleatoires partages
                    call Generation_surrogate(don_simul,don_simulS1,ng,n_col,logNormal,affiche_stat,theta,&
                        ng,ver,alpha,cens0,temps_cens,gamma1,gamma2,theta2,lambdas,nus,lambdat,nut,vbetas,vbetat,&
                        n_essai,rsqrt,sigma_s,sigma_t,p,prop_i,gamma_ui,alpha_ui,frailt_base,pfs)
                endif
                if(nsim_node(11)==3) then ! joint frailty copula model
                    call Generation_surrogate_copula(don_simultamp,don_simulStamp,ng,n_col,logNormal,affiche_stat,theta,&
                        ng,ver,alpha,cens0,temps_cens,gamma1,gamma2,theta2,lambdas,nus,lambdat,nut,vbetas,vbetat,&
                        n_essai,rsqrt,sigma_s,sigma_t,p,prop_i,gamma_ui,alpha_ui,frailt_base,thetacopule, filtre,&
                        filtre2,pfs)
                    don_simul(:,1:4) = don_simultamp(:,1:4)
                    don_simul(:,5) = 0.d0 ! on le met a 0 car je ne prends pas en compte les w_ij au moment de generation avec les copule. ducoup matricce avec -1 colone, par rapport a la generation a partir du modele joint surrogate
                    don_simul(:,6:size(don_simul,2)) = don_simultamp(:,5:(size(don_simultamp,2)-1)) ! -1 car la derniere colonne n'est pas rempli
                    don_simulS1(:,1:4) = don_simulStamp(:,1:4)
                    don_simulS1(:,5) = 0.d0
                    don_simulS1(:,6:size(don_simulS1,2)) = don_simulStamp(:,5:(size(don_simulStamp,2)-1)) ! -1 car la derniere colonne n'est pas rempli

                endif

                if(nsim_node(11)==2) then
                endif

                !==================22/05/2019 Fin commentaires==== =================================
                !==================================================================

                ! si je suis avec des paquets de simulation, alors quand je suis la, je genere les premier jeu de donnees a ne pas considere, apres reinitialisation de l'environnement de generation
                indice_seed=indice_seed+1
                if(indice_seed<(s_i+decoup_simul)) goto 11

                if((s_i<init_i).or.s_i>max_i) then
                    debut_exe=0 ! pour dire le processus ne considere pas ce jeu de donnee
                else
                    debut_exe=1 !jeux de donnees a considerer pour le processus courant
                endif

                ! jeux de donnees a considerer pour le processus courant
            endif
            ! ind_temp=ng
            !sauvegarde de tous le jeu de donnees simulees
            do k=1,ng

            enddo
        endif


        ind_temp=ng
        allocate(don_simulS(ind_temp,n_col))
        allocate(tableEssai(n_essai))
        tableEssai=0
        don_simulS=don_simulS1(1:ind_temp,:) ! on recupere les significatives du tableaux don_simulS1
        nsujet=ind_temp ! on met à jour le nombre d'observations sur les surrogates

        !---------------------->  Deces
        trials=0 ! vecteur contenant le nombre de sujet par etude
        vaxdc=0.d0
        Deces=0
        randomisation=0
        surrogate=0
        nigs=0
        cdcs=0
        nigts=0
        cdcts=0

        do i= 1,ng !sur les groupes uniquement (= sujet)
            !read(10,*)
            tt0dc(i)=don_simul(i,initTime1)/revision_echelle ! temps initial
            tt1dc(i)=don_simul(i,timeT1)/revision_echelle! temps de deces
            icdc(i)=Int(don_simul(i,statusT1)) ! delta_star
            pourtrial(i)=INT(don_simul(i,trialref1)) !indice de l'essai
            groupe(i)=INT(don_simul(i,Patienref1)) ! numero de l'individu
            vaxdct(i,1)=don_simul(i,trt1) ! vecteur des variables explicatives
            tableEssai(pourtrial(i))=tableEssai(pourtrial(i))+1
            ! j'ajoute les autres variables a la fin

            if(type_joint_estim == 3) then! joint frailty copula
                if(ver > 1) then ! I add the rest of covariates
                    vaxdct(i,2:ver) = don_simul(i,(size(don_simul,2) - ver +2):size(don_simul,2))
                endif
            endif

            if(une_donnee==0 .and. gener_only==1)then
                if(seed_.eq.0) then
                endif
                        !stop
                if(s_i.eq.seed_) then
                endif
            endif
            k=0
            do j=1,ver
                if (filtre2(j).eq.1)then
                    k=k+1
                    vaxdc(i,k)=vaxdct(i,j) ! on maintient parmi les variables explicatives de deces celles indiquees dans le filtre
                end if
            end do
            ! elements de statistique
            trials(int(pourtrial(i)))=trials(int(pourtrial(i)))+1 ! on incremente le nombre de personnes dans letude
            Deces(int(icdc(i)))=Deces(int(icdc(i)))+1 ! compte le nombre de decedes et de vivant
            randomisation(int(vaxdct(i,1)))=randomisation(int(vaxdct(i,1)))+1 ! compte le nombre de personne sous traitements et non traites
            if(icdc(i).eq.1) then !on incremente le nombre de personnes decedees par essai
                cdcs(int(pourtrial(i)))=cdcs(int(pourtrial(i)))+1 ! nombre de cedes par essai
                cdcts(int(pourtrial(i)))=cdcts(int(pourtrial(i)))+int(vaxdct(i,1))! nombre de deces traites par essai
            endif
        end do
        deallocate(groupe)
        allocate(groupe(nsujet))
        !close(10)

        !---------------------->  Surrogate
        vax=0.d0
        do i = 1,nsujet     !sur les observations
            k=0
            !read(9,*)tt0(i),tt1(i),ic(i),pourtrial(i),groupe(i),(vaxt(i,j),j=1,ver)
            tt0(i)=don_simulS(i,initTime1)/revision_echelle ! temps initial
            tt1(i)=don_simulS(i,timeS1)/revision_echelle! temps de deces
            ic(i)=Int(don_simulS(i,statusS1)) ! delta_star
            pourtrial(i)=INT(don_simulS(i,trialref1)) !indice de l'essai
            groupe(i)=INT(don_simulS(i,Patienref1)) ! numero de l'individu
            vaxt(i,1)=don_simulS(i,trt1) ! vecteur des variables explicatives
            ! j'ajoute les autres variables a la fin

            if(type_joint_estim == 3) then! joint frailty copula
                if(ver > 1) then ! I add the rest of covariates
                    vaxt(i,2:ver) = don_simulS(i,(size(don_simulS,2) - ver +2):size(don_simulS,2))
                endif
            endif

            !ecriture des donnees dans le fichier (juste pour le premier jeux de donnee)
            ! on sauvegarde seulement si on est dans la simple generation des donnees

            if(une_donnee==0  .and. gener_only==1)then
                if(seed_==0) then
                 !   if(rang_proc==0) !write(9,*)tt0(i),tt1(i),ic(i),pourtrial(i),groupe(i),(vaxt(i,j),j=1,ver)
                endif

                if(s_i==seed_) then
                 !   if(rang_proc==0) !write(9,*)tt0(i),tt1(i),ic(i),pourtrial(i),groupe(i),(vaxt(i,j),j=1,ver)
                endif
            endif

            do j=1,ver
                if (filtre(j).eq.1)then
                    k=k+1
                    vax(i,k)=vaxt(i,j) ! on maintient parmi les variables explicatives de surrogacy celles indiquees dans le filtre
                end if
            end do
            surrogate(int(ic(i)))=surrogate(int(ic(i)))+1 ! compte le nombre de personne avec et sans progression
            if(ic(i).eq.1) then !on incremente le nombre de personnes avec une progression par essai
                nigs(int(pourtrial(i)))=nigs(int(pourtrial(i)))+1 ! nbre personnes avec une progression
                nigts(int(pourtrial(i)))=nigts(int(pourtrial(i)))+int(vaxt(i,1)) ! nbre personnes avec une progression traites
            endif
        end do
        !close(9)

        ! construction de la matrice du nombre d'observations
        nig_Ts(:,1)=nigs
        nig_Ts(:,2)=nigts
        cdc_Ts(:,1)=cdcs
        cdc_Ts(:,2)=cdcts

        allocate(vdeces(deces(1)),vsurrogate(surrogate(1)))
        ! on recupere les dates de deces et de progression: cas des progression sans recidive

        k=1
        j=1
        do i=1,nsujet
            if(ic(i).eq.1) then
                vsurrogate(k)=int(tt1(i))
                k=k+1
            end if
            if(icdc(i).eq.1) then
                vdeces(j)=int(tt1dc(i))
                j=j+1
            end if
        end do

        descripSurr(1)=minval(vsurrogate)
        descripSurr(2)=maxval(vsurrogate)
        descripSurr(3)=Median(vsurrogate, surrogate(1))

        descripDeces(1)=minval(vdeces)
        descripDeces(2)=maxval(vdeces)
        descripDeces(3)=Median(vdeces, deces(1))


        !====================================================================
        !                         quelques statistiques
        !====================================================================
        if(rang_proc==0) then
        endif
        if(rang_proc==0) then
        endif
        allocate(paGH(ng,nsim_node(7)+1+nsim_node(7) + (nsim_node(7)*(nsim_node(7)-1))/2))
        paGH=0.d0

        select case(typeof) !nombre de paramètes
            case(0)
                np = 2*(nz+2) + nva + nparamfrail + nknots + splines_ord !scl v // qlc
            case(1)
                np = nbintervDC + nbintervR + nva + nparamfrail +1 + nknots
            case(2)
                np = 2*nst + nva + effet + nparamfrail +1 + nknots
        end select

        !allocate(b(np),H_hessOut(np,np),HIHOut(np,np))
        mtaille(1) = mt1 ! vecteur contenant la taille des vecteurs utilises pour ploter lees fonction de risque et de survie
        mtaille(2) = mt2
        mtaille(3) = mt11
        mtaille(4) = mt12

        ! paramete de la weibull a prendre en entree
        shape_weib(1) = paraweib(1)  ! pour le surrogate
        shape_weib(2) = paraweib(2)  ! pour le deces
        scale_weib(1) = paraweib(3)  ! pour le surrogate
        scale_weib(2) = paraweib(4)  ! pour le deces

        paratps = 0 ! indicateur de presence de variable temps-dependente

        MartinGales = 0 ! matrice des residus de martingale pour surrogate, deces, la predistion des fragilits, variance des fragilites predites

        ttU = 0.d0        ! en cas de censure par intervalle, borne superieure de l'intervalle


        if(gener_only.eq.1) then ! ceci voudrait dire qu'on voudrait simplement generer les donnees pas d'estimation
          !if(rang_proc==0) !print*,"suis a la simulation",s_i
          goto 1000
        endif

        !==========Gestion des kappa====================
        !===============================================
        ! ici on considere soit le premier kappa pour toutes les simulation, soit un kappa par simul ou alors un kappa par simul avec possibilite en cas de non convergence
        ! de rentrer utiliser un des 4 premiers kappas ayant permis la convergence, dans ce cas, si le premier kappa ne permet pas la convergence, on prend le suivant
        statut_kappa=0 ! pour dire que c'est un nouveau kappa, s'il permet la convergence on le retient pour les cas de non convergence
        ind_rech=1 ! s'incremente a la sortie du joint si le modele n'a pas converge
        Rech_kappa=0 ! controle la convergence du modele en faisant varie le kappa
        control_kappa=0 ! controle si on a deja divise les kappas dans le cas kappa_use=4
        ind_premier_kappa=0 ! qui est incremente a chaque fois que le kappa considere ne permet pas la convergence du premier jeu de donnee si on teste 5 kappa sur le premier jeu et on n'a pas de resultat, alors on continu avec le jeu de donnee suivant
        control=0 ! pour controler l'acces unique a la modification du kappa en cas de non convergence et lorsque le changement du nbre de point et des kappas issus de la cross validation ne marchent pas
        control2=0! pour s'assurer qu'on ne boucle pas lorsque ni le changement de kappa ni le nbre de point ne permet pas la convergence (en fait on doit relancer le modele une seule fois en changeant simultanement les deux parametre)
        2001 continue
        if(kapa_use.ne.0) then ! on usilise un nouveau kappa pour chaque jeu de donnees
            !if(une_donnee.ne.1 .or. donne_reel .ne.1)read(15,*)ax1,ax2 ! si les deux vallent un alors on utilise les kappas dournis dans le fichiers des parametres: joint_scl_simul
            if(une_donnee.ne.1 .or. donne_reel .ne.1)then
                ax1 = vect_kappa(indice_kapa,1)
                ax2 = vect_kappa(indice_kapa,2)
                indice_kapa = indice_kapa+1 ! si les deux vallent un alors on utilise les kappas dournis dans le fichiers des parametres: joint_scl_simul
            end if
            k0(1)=ax1
            k0(2)=ax2
            k0_save=k0
            if(s_i==1)    k01_save=k0 ! on sauvegarde le premier kappa a utiliser plutard en cas de besoin
            !!print*,"kappa considere vaut:",k01_save
            if(statut_kappa1==1)then ! on utilise le kappa suivant (deuxieme)

                indice_kapa = 1
                !read(15,*)ax1,ax2
                ax1 = vect_kappa(indice_kapa,1)
                ax2 = vect_kappa(indice_kapa,2)
                indice_kapa = indice_kapa + 1
                statut_kappa1=0
            endif
        else
            if((s_i.eq.1).or.((statut_kappa1==0).and.(ind_rech<=n_sim))) then !on considere le premier kappa qui marche pour toute les simul
                !if(une_donnee.ne.1 .or. donne_reel .ne.1) read(15,*)ax1,ax2  ! si les deux vallent un alors on utilise les kappas dournis dans le fichiers des parametres: joint_scl_simul
                if(une_donnee.ne.1 .or. donne_reel .ne.1)then
                    ax1 = vect_kappa(indice_kapa,1)
                    ax2 = vect_kappa(indice_kapa,2)
                    indice_kapa = indice_kapa +1! si les deux vallent un alors on utilise les kappas dournis dans le fichiers des parametres: joint_scl_simul
                 end if
                k0(1)=ax1
                k0(2)=ax2
                k0_save=k0
                !statut_kappa1=1
                ind_rech=ind_rech+1
            endif
        endif



        2000 continue
        2003 continue
        if (((Rech_kappa==1).and.((kapa_use.eq.2).or.(kapa_use.eq.3).or.(kapa_use.eq.4))) .or.(controlgoto==1)) then ! on recherhce parmi les kappas deja utilise s'il ya un ki permet la convergence
            if(((kapa_use.eq.2).or.(kapa_use.eq.4)).and.(controlgoto==0))then
                k0(1)=kappa(ind_rech-1,1)
                k0(2)=kappa(ind_rech-1,2)
            else
                controlgoto=0 ! ceci me permet de mieux gerer le goto 2003
                if((k0(1)-k0(2))>=1000)then
                    k0(1)=k0(1)/10.d0
                else if((k0(1)-k0(2))<=-1000)then
                    k0(2)=k0(2)/10.d0
                else
                    k0(1)=k0(1)/10.d0
                    k0(2)=k0(2)/10.d0
                endif
                !if(rang_proc==0) !print*,"division par 10 de kappa"
                goto 2002
            endif
        endif
        2002 continue
        ! initialisation du vecteur b des parametres a l'aide des parametres de simulation
        b=0.5d0
        if(nsim_node(8)==0)then !model conjoint surrogate classique
            !if(logNormal==1)then
                b(np-2-nparamfrail+indice_eta+indice_theta)=theta !theta2
                !b(np-2-nparamfrail+indice_eta+indice_theta)=0.5d0 !theta2
            !else
            !    b(np-2-nparamfrail+indice_eta+indice_theta)=dsqrt(theta) !theta
            !endif
        else
            if(nsim_node(8)==1)then !effets aleatoires partages
                !b(np-2-nparamfrail+indice_eta+indice_theta)=dsqrt(theta)
                !b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS)=dsqrt(sigma_s)
                !b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT)=dsqrt(sigma_t)
                b(np-2-nknots-splines_ord-nparamfrail+indice_eta+&
                indice_theta)=dsqrt(theta)
                b(np-2-nknots-splines_ord-nparamfrail+indice_eta+&
                indice_theta+indice_varS)=dsqrt(sigma_s)
                b(np-2-nknots-splines_ord-nparamfrail+indice_eta+&
                indice_theta+indice_varS+indice_varT)=dsqrt(sigma_t)
                if(indice_covST==1) then
                    !b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST)=&
                    !rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
                    b(np-2-nknots-splines_ord-nparamfrail+indice_eta+&
                    indice_theta+indice_varS+indice_varT+indice_covST)=&
                    rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
                endif
                if(frailt_base==1) then
                    !if(indice_alpha==1) b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+&
                    !indice_alpha)=alpha_ui
                    !b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+&
                    !indice_gamma)=dsqrt(gamma_ui)
                    if(indice_alpha==1) b(np-2-nknots-splines_ord-&
                    nparamfrail+indice_eta+indice_theta+indice_varS+&
                    indice_varT+indice_covST+indice_alpha)=alpha_ui
                    b(np-2-nknots-splines_ord-nparamfrail+indice_eta+&
                    indice_theta+indice_varS+indice_varT+&
                    indice_covST+indice_alpha+indice_gamma)=dsqrt(gamma_ui)
                endif
            endif

            if(nsim_node(8)==2)then !effets aleatoires correles
                b(np-2-nparamfrail+indice_eta+indice_theta)=dsqrt(theta)
                b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS)=&
                dsqrt(sigma_s)
                b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+&
                indice_varT)=dsqrt(sigma_t)
                b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+&
                indice_varT+indice_covST)&
                =rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
                b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+&
                indice_varT+indice_covST+indice_gamma)&
                      =dsqrt(gamma_ui)
                b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+&
                indice_varT+indice_covST+indice_gamma+&
                      indice_theta_t)=dsqrt(theta2_t)
                b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+&
                indice_varT+indice_covST+indice_gamma+&
                        indice_theta_t+indice_theta_st)=&
                        rsqrt_theta*dsqrt(theta)*dsqrt(theta2_t)
                b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+&
                indice_varT+indice_covST+indice_gamma+&
                        indice_theta_t+indice_theta_st+indice_gamma_t)=&
                        dsqrt(gamma_uit)
                if(indice_gamma_st==1)b(np-2-nparamfrail+indice_eta+&
                indice_theta+indice_varS+indice_varT+&
                indice_covST+indice_gamma+indice_theta_t+indice_theta_st+&
                indice_gamma_t+&
                indice_gamma_st)=rsqrt_gamma_ui*dsqrt(gamma_ui)*dsqrt(gamma_uit)
            endif
            if(nsim_node(8)==3)then !joint frailty-copula
                b(np-nva-nparamfrail+indice_varS)=dsqrt(sigma_s)
                b(np-nva-nparamfrail+indice_varS+indice_varT)=dsqrt(sigma_t)
                if(indice_covST==1) then
                    b(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST)=&
                    rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
                endif
                if(frailt_base==1) then
                    if(indice_alpha==1) b(np-nva-nparamfrail+indice_varS+&
                    indice_varT+indice_covST+&
                        indice_alpha)= alpha_ui
                    b(np-nva-nparamfrail+indice_varS+indice_varT+&
                    indice_covST+indice_alpha+&
                    indice_gamma)=dsqrt(gamma_ui)
                endif
                if(copula_function == 1) b(np-nva) = dlog(thetacopule) ! claton: exp transform
                if(copula_function == 2) b(np-nva) = dsqrt(thetacopule)  ! Gumbel: choleschy transform
                b((np-nva + 1) : (np - nva + nva1)) = vbetas
                b((np-nva2 + 1) : np) = vbetat
            endif
        endif

        !affectation manuelle des parametres initiaux (choleschy obtenue de R)
        if(vrai_val_init==1)then ! on initialise avec les vrai valeurs
            if(theta>=1.d0)then
                theta_init=0.7d0 ! on fait ceci car pour un theta initialiser a 1 le modèle est difficile a estimer
            else
                theta_init=theta
            endif
            theta_init=theta
            sigma_ss_init=sigma_s
            sigma_tt_init=sigma_t
            gamma_init=gamma_ui
            alpha_init=alpha_ui
            if(indice_covST==1) then
                sigma_st_init=rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
            endif
            betas_init=betas
            betat_init=betat
            gamma_init=gamma_ui
            alpha_init=alpha_ui
            zeta_init=eta

            ! pour le modele complet
            thetat_init=theta2_t
            thetast_init=rsqrt_theta*dsqrt(theta)*dsqrt(theta2_t)
            gammat_init=gamma_uit
            gammast_init=rsqrt_gamma_ui*dsqrt(gamma_ui)*dsqrt(gamma_uit)

            ! frailty-copula
            thetacopula_init = thetacopule
            vbetas_intit = vbetas
            vbetat_intit = vbetat
        endif

        if(sigma_ss_init.eq.0.d0) then
            !if(rang_proc==0) !print*,"Attention: revoir la valeur initiale de sigma_SS"
        else
            if((sigma_tt_init-(sigma_st_init**2.d0)/sigma_ss_init)<0)then
            !    if(rang_proc==0) !print*,"Attention: revoir les valeurs initiales de la matrice sigma des effets aleatoire"
            endif
        endif

        if(nsim_node(8)==3)then !joint frailty-copula
            b(np-nva-nparamfrail+indice_varS)=dsqrt(sigma_ss_init)
            b(np-nva-nparamfrail+indice_varS+indice_varT)=dsqrt(sigma_tt_init-&
                    (sigma_st_init**2.d0)/sigma_ss_init)
            if(indice_covST==1) then
                b(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST)=&
                sigma_st_init/dsqrt(sigma_ss_init)
            endif
            if(frailt_base==1) then
                if(indice_alpha==1) b(np-nva-nparamfrail+indice_varS+&
                    indice_varT+indice_covST+indice_alpha)=alpha_init
                    b(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST+&
                    indice_alpha+indice_gamma)=dsqrt(gamma_init)
                endif
            if(copula_function == 1)  b(np-nva) =  dlog(thetacopula_init)! claton: exp transform
            if(copula_function == 2) b(np-nva) = dsqrt(thetacopula_init)  ! Gumbel: choleschy transform
            b((np-nva + 1) : (np - nva + nva1)) = vbetas_intit
            b((np-nva2 + 1) : np) = vbetat_intit
        else
            b(np-nknots-splines_ord-2)=betas_init !b_s
            b(np-nknots-splines_ord-1)=betat_init !b_t
            b((np-nknots-splines_ord+1):np)=0.5d0 !b_gamma
            if(indice_eta==1)then
                b(np-2-nknots-splines_ord-nparamfrail+indice_eta)=zeta_init !eta
            endif
        endif

        ! si modele complet avec effets aleatoires partages
        !if(nsim_node(8)==1)then
        !    b(np-2-nparamfrail+indice_eta+indice_theta)=dsqrt(theta_init)
        !    b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS)=dsqrt(sigma_ss_init)
        !    b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT)=dsqrt(sigma_tt_init-&
        !            (sigma_st_init**2.d0)/sigma_ss_init)
        !    if(indice_covST==1) then
        !        b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST)=sigma_st_init/&
        !            dsqrt(sigma_ss_init)
        !     endif
        !    if(frailt_base==1) then
        !        if(indice_alpha==1) b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+&
        !            indice_alpha)=alpha_init
        !        b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+&
        !        indice_gamma)=dsqrt(gamma_init)
        !    endif
        !endif
        if(nsim_node(8)==1)then
            b(np-2-nknots-splines_ord-nparamfrail+indice_eta+indice_theta)=&
            dsqrt(theta_init)
            b(np-2-nknots-splines_ord-nparamfrail+indice_eta+indice_theta+&
            indice_varS)=dsqrt(sigma_ss_init)
            b(np-2-nknots-splines_ord-nparamfrail+indice_eta+indice_theta+&
            indice_varS+indice_varT)=dsqrt(sigma_tt_init-&
                    (sigma_st_init**2.d0)/sigma_ss_init)
            if(indice_covST==1) then
                b(np-2-nknots-splines_ord-nparamfrail+indice_eta+indice_theta+&
                  indice_varS+indice_varT+indice_covST)=&
                  sigma_st_init/dsqrt(sigma_ss_init)
             endif
            if(frailt_base==1) then
                if(indice_alpha==1) b(np-2-nknots-splines_ord-nparamfrail+&
                  indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+&
                  indice_alpha)=alpha_init
                b(np-2-nknots-splines_ord-nparamfrail+indice_eta+&
                indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+&
                indice_gamma)=dsqrt(gamma_init)
            endif
        endif
        ! si modele complet avec effets aleatoires correles
        if(nsim_node(8)==2)then
            b(np-2-nparamfrail+indice_eta+indice_theta)=dsqrt(theta_init)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS)=&
            dsqrt(sigma_ss_init)
            b(np-2-nparamfrail+indice_eta+indice_theta+&
            indice_varS+indice_varT)&
            =sigma_st_init/dsqrt(sigma_ss_init)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+&
            indice_varT+indice_covST)=dsqrt(sigma_tt_init-&
            (sigma_st_init**2.d0)/sigma_ss_init)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+&
            indice_varT+indice_covST+indice_alpha+&
            indice_gamma)=dsqrt(gamma_init)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS&
            +indice_varT+indice_covST+indice_gamma+&
            indice_theta_t)=thetast_init/dsqrt(theta_init)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+&
            indice_varT+indice_covST+indice_gamma+&
            indice_theta_t+indice_theta_st)=&
            dsqrt(thetat_init-(thetast_init**2.d0)/theta_init)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+&
            indice_varT+indice_covST+indice_gamma+&
            indice_theta_t+indice_theta_st+indice_gamma_t)=&
            gammast_init/dsqrt(gamma_init)
            if(indice_gamma_st==1)b(np-2-nparamfrail+indice_eta+&
            indice_theta+indice_varS+indice_varT+&
            indice_covST+indice_gamma+indice_theta_t+indice_theta_st&
            +indice_gamma_t+&
            indice_gamma_st)=dsqrt(gammat_init-(gammast_init**2.d0)/gamma_init)
        endif
        if((nsim_node(4).ne.0) .and.(nsim_node(4).ne.3))then
            !if(rang_proc==0) !print*,nsim_node(2),np_save
            if((nsim_node(2).ne.np_save).and.(s_i==1))then
                !if(rang_proc==0) !print*,"le nombre de point de quadrature est passé de",np_save,"a:",nsim_node(2)
            endif
        endif


        2005 continue

        EPS =EPS2    ! critere de convergence du modele

        if(rang_proc==0) then

        endif
        call cpu_time(tp1)

        if(debut_exe==1)then
            !if(rang_proc==0) !print*,"Je suis à la simulation",s_i,"kappa=",k0,"EPS=",EPS
        else
            goto 1002 ! le processus courant n'utilise pas ce jeu de donnee
        endif


        k0(1) = k0(1) + ckappa(1)
        k0(2) = k0(2) + ckappa(2)

        if(affiche_itteration == 1) then

        endif
        Call joint_surrogate(nsujet,ng,ntrials,0,nz,nst,k0,tt0,tt1,ic,groupe,trials,pourtrial,nig_Ts,cdc_Ts,0, &
                            tt0dc,tt1dc,icdc,0.d0,0,nva1,vax,nva2,vaxdc,0.d0,noVar1,noVar2,ag,maxiter,np,b,knotsurro,H_hessOut,&
                            HIHOut,resOut,LCV,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out,typeof,equidistant,&
                            nbintervR,nbintervDC,mtaille,ni,cpt,cpt_dc,ier,istop,paraweib,MartinGales,linearpred,&
                            linearpreddc,ziOut,time,timedc,0.d0 , 1 , 0 , ttU , logNormal , paratps , 0 , 0.d0 , 0.d0 , &
                            EPS,nsim_node,indice_esti,indice_covST,0.d0,param_weibull)
        !! compute natural effects
        if(mediation.and.(istop.eq.1))then
          call compute_rt(b,HIHOut,knotsurro,np,(/nz+2,nz+2,nva,nparamfrail,&
                          indice_eta,indice_alpha,nknots,splines_ord,&
                          mtaille(1),mtaille(2)/),ziOut,res_rt(:,1,1),res_rt,rt_boot,&
                          integ_type)
          !deallocate(zmed,zdcmed)
        end if
        if (istop.ne.1) then
            !call intpr("je suis la :", -1, ni, 1)
            ! if(nsim_node(3).ne.1) then ! cas ou on ne fait pas de l'adaptative
                if((nsim_node(4).ne.0) .and.(nsim_node(4).ne.3) .and. (kapa_use.eq.4))then !cas quadrature
                    if ((control_kappa==3) .or. (control_kappa==4)) then ! on emet les vraies valeurs
                        if(np_save==npoint1) nsim_node(2)=npoint2
                        if(np_save==npoint2) nsim_node(2)=npoint1
                    endif

                    if((nsim_node(2).eq.np_save))then ! on change le nombre de point de quadrature avant de continuer
                        if(np_save==npoint1) nsim_node(2)=npoint2
                        if(np_save==npoint2) nsim_node(2)=npoint1
                        goto 2002 ! apres changement du nbre de point on relance l'estimation
                    else
                        if(control_kappa<4)then ! on vient de changer le nombre de point alors on revient
                            if(control_kappa<2)then ! on va aller divise les kappas par 10, en maintenant change le nbre de point
                                !if(rang_proc==0) !print*,"division du kappa par 10 ou par 100"
                                control_kappa=control_kappa+1
                                controlgoto=1
                               ! goto 2003
                            else
                                if(control_kappa==2)then ! on remet les valeurs de kappa pour la suite
                                    k0=k0_save
                                endif
                            endif
                            if(controlgoto==0)then
                                nsim_node(2)=np_save
                                if(control_kappa<4)then ! on relance le modele avec les kappas diviss par 10, sans change le nbre de point
                                    !if(rang_proc==0) !print*,"division du kappa par 10 ou par 100 avec maintien du nombre de", &
                                     !   "points de quadrature"
                                    control_kappa=control_kappa+1
                                    controlgoto=1
                                    !goto 2003
                                endif
                            endif
                        else !on a deja changer aussi bien le nbre de point que le kappa et xa ne marche toujours pas, on va donc changer les deux
                            if(control_kappa==4)then !dans ce cas on change les kappa et les points de quadrature
                                if(np_save==npoint1) nsim_node(2)=npoint2
                                if(np_save==npoint2) nsim_node(2)=npoint1

                                if(control2==0) then !on considere dans ce cas le dernier kappa/100 avec changement aussi du nombre de points de quadrature
                                   ! if(rang_proc==0) !print*,"changement simultané du nbre de point et kappa division par 100"
                                    control2=1 ! pour controler qu'on relance une seule fois le modele
                                    goto 2002 !
                                endif

                                if(s_i>1) then
                                    if(control2==1)then ! on considere dans ce cas le premier kappa qui marche ou le premier kappa avec les points de quadrature de depart
                                        nsim_node(2)=np_save
                                        control2=2
                                   !     if(rang_proc==0) !print*,"changement du kappa (premier kappa qui marche ou tout premier si", &
                                     !       "aucun n'a marché jusqu'ici)"
                                        if(incre_kappa>=2) then ! on a aumoins un kappa qui marche on recupere le tout premier
                                            k0(1)=kappa(1,1)
                                            k0(2)=kappa(1,2)
                                        else ! sinon on recupere le tout premier kappa
                                            k0=k01_save
                                        endif
                                        goto 2002 ! apres changement du nbre de point on relance l'estimation une seule fois
                                    endif
                                endif

                                if(s_i>1) then
                                    if(control2==2)then ! on considere dans ce cas le premier kappa qui marche ou le premier kappa avec changement aussi du nombre de points de quadrature
                                        control2=3
                                       ! if(rang_proc==0) !print*,"changement simultane du nbre de point et kappa (premier kappa", &
                                        !        "qui marche ou tout premier si aucun n'a marche jusqu'ici)"
                                        if(incre_kappa>=2) then ! on a aumoins un kappa qui marche on recupere le tout premier
                                            k0(1)=kappa(1,1)
                                            k0(2)=kappa(1,2)
                                        else ! sinon on recupere le tout premier kappa
                                            k0=k01_save
                                        endif
                                        goto 2002 ! apres changement du nbre de point on relance l'estimation une seule fois
                                    endif
                                endif
                            endif
                        endif
                        if(controlgoto==0)then
                           ! if(rang_proc==0) !print*,"suis en fin là, control_kappa vaut",control_kappa
                        endif
                    endif

                    if(controlgoto==1) goto 2003
                    if(ind_rech==1)then
                        if(incre_kappa>=2) then ! on a aumoins un kappa qui marche on recupere le tout premier
                            k0(1)=kappa(1,1)
                            k0(2)=kappa(1,2)
                            ind_rech=2    ! vu que j'ai utilise le premier kappa qui a permis la convergence alors je passe au suivant
                        else ! sinon on recupere le tout premier kappa
                            k0=k01_save
                        endif
                        !!print*,"k0 pour cette nouvelle relance vaut",k0
                    endif

                    if(control==0) then ! on relance le modèle avec le premier kappa (converge ou pas) en changent les valeurs initiales aux dernieres estimations qui n'ont pas permises la convergence
                        control=1
                        goto 2005    ! les parametres sont initialises avec les dernieres estimations n'ayant pas permis la convergence
                     else
                        if(control==1) then ! ici j'utilise le nombre de points de quadrature de base et je change juste le kappa et les valeurs initiales
                           ! if(rang_proc==0) !print*,"changement du kappa (premier kappa qui marche ou tout", &
                             !       "premier si aucun n'a marché jusqu'ici) ainsi que des valeurs initiales des parametres"
                            nsim_node(2)=np_save
                            control=2
                            goto 2005    ! les parametres sont initialises avec les dernieres estimations n'ayant pas permis la convergence
                         endif
                    endif
                endif
            ! endif

            if((nsim_node(4).eq.3 .or. nsim_node(4).eq.0 .or. ((nsim_node(4).eq.1 .or. nsim_node(4).eq.2).and.&
                nsim_node(3).eq.1)) .and. (kapa_use.eq.4))then !cas Monte carlo et Laplace ou alors pseudo adaptative

                    if(control_kappa<2)then ! on vient de changer le nombre de point alors on revient
                        if(control_kappa<2)then ! on va aller divise les kappas par 10, en maintenant change le nbre de point
                          !  if(rang_proc==0) !print*,"division du kappa par 10 ou par 100"
                            control_kappa=control_kappa+1
                            controlgoto=1
                            !goto 2003
                        else
                            if(control_kappa==2)then ! on remet les valeurs de kappa pour la suite
                                k0=k0_save
                            endif
                        endif
                    else !on a deja changer le kappa et xa ne marche toujours pas
                        if(control_kappa==2)then !dans ce cas on change les kappa et les points de quadrature
                            if(s_i>1) then
                                if(control2==1)then ! on considere dans ce cas le premier kappa qui marche ou le premier kappa avec les points de quadrature de depart
                                    control2=2
                                    ! if(rang_proc==0) !print*,"changement du kappa (premier kappa qui marche ou tout premier si ",&
                                        ! "aucun n'a marché jusqu'ici)"
                                    if(incre_kappa>=2) then ! on a aumoins un kappa qui marche on recupere le tout premier
                                        k0(1)=kappa(1,1)
                                        k0(2)=kappa(1,2)
                                    else ! sinon on recupere le tout premier kappa
                                        k0=k01_save
                                    endif
                                    goto 2002 ! apres changement du kappa on relance l'estimation une seule fois
                                endif
                            endif
                        endif
                    endif
                if(controlgoto==1) goto 2003
                if(ind_rech==1)then
                    if(incre_kappa>=2) then ! on a aumoins un kappa qui marche on recupere le tout premier
                        k0(1)=kappa(1,1)
                        k0(2)=kappa(1,2)
                        ind_rech=2    ! vu que j'ai utilise le premier kappa qui a permis la convergence alors je passe au suivant
                    else ! sinon on recupere le tout premier kappa
                        k0=k01_save
                    endif
                    !!print*,"k0 pour cette nouvelle relance vaut",k0
                endif

                if(control==0) then ! on relance le modèle avec le premier kappa (converge ou pas) en changent les valeurs initiales aux dernieres estimations qui n'ont pas permises la convergence
                    control=1
                    goto 2005    ! les parametres sont initialises avec les dernieres estimations n'ayant pas permis la convergence
                endif
            endif

            ! if(rang_proc==0) !print*,"je vais a la recherche des kappas qui on converges pour n'avoir rien trouver"
            ! if(rang_proc==0) !print*,"ind_rech=",ind_rech,"incre_kappa=",incre_kappa
            if((ind_rech<incre_kappa).and.((kapa_use.eq.2).or.(kapa_use.eq.3).or.(kapa_use.eq.4))) then ! on recherhce parmi les kappas deja utilise s'il ya un qui permet la convergence
                ! if(rang_proc==0) !print*,"test du kappa numero",ind_rech-1,"qui a deja eu à converger"
                Rech_kappa=1
                ind_rech=ind_rech+1
                statut_kappa=1 ! evite qu'on enregistre le kappa trouve s'il permet la convergence

                if(ind_rech==4)then ! on s'arrete et passe a la suivante simulation
                    goto 1000
                endif
                goto 2000
            endif

            ! si le tout premier kappa ne converge pas
            !if((s_i==1).and.(n_sim.ne.1).and.((kapa_use.eq.2).or.(kapa_use.eq.3).or.(kapa_use.eq.4))) then
            if((s_i==1).and.(ind_premier_kappa==4))then

            endif

            if((s_i==1).and.(n_sim.ne.1) .and.(ind_premier_kappa<4)) then
                ! if(rang_proc==0) !print*,"suis la pour recherche kappa: on utilise le suivant"
                ind_premier_kappa=ind_premier_kappa+1
                statut_kappa1=1
                goto 2001 ! on utilise le kappa suivant
            endif
        else
            if(kapa_use==0)then ! si c'est un seul kappa qu'il faut utiliser alors on le maintient et on continu
                statut_kappa1=1
            endif
        endif

        1000 continue
        Rech_kappa=0
        ind_rech=1
        ! extraction des fragilites niveau essai pour evaluation
        k=1
        do i=1,n_essai
            !if(nsim_node(8)==2 .or. nsim_node(8)==3)then
            if(nsim_node(8)==2)then
                donnee_essai(i,:)=don_simulS(k,(/v_s1,v_t1,trialref1,u_i1,u_it/))
            else
                donnee_essai(i,:)=don_simulS(k,(/v_s1,v_t1,trialref1,u_i1/))
            endif
            k=k+tableEssai(i)
            !!print*,donnee_essai(i,4)
        enddo
        nsim_node(2)=np_save ! on remet le nombre de point de quadrature s'il y'a eu modification

        if (istop.ne.1) then
            if(rang_proc==0) then
                code_print = 0
            endif
            ! si pas de convergence on simule un nouveau jeux de donnees
            nbre_rejet=nbre_rejet+1
            tab_rejet(nbre_rejet)=s_i

            ! controle de l'ecriture dans le fichier
            erreur_fichier=0
            996 continue
            !write(8,*,iostat=erreur_fichier)s_i
            if(erreur_fichier .ne.0) then
                goto 996
            endif

            ! parametres empiriques en cas de non convergence
            parametre_empirique_NC(nbre_rejet,1)=sum(don_simul(:,trt1))*100.d0/n_obs! personne traitees
            parametre_empirique_NC(nbre_rejet,2)=variance(don_simul(:,w_ij1))! theta empirique
            parametre_empirique_NC(nbre_rejet,3)=sum(don_simul(:,statusS1))*100.d0/n_obs! progression
            parametre_empirique_NC(nbre_rejet,4)=sum(don_simul(:,statusT1))*100.d0/n_obs! deces

            if(nsim_node(8).ne.0)then !model conjoint surrogate
                parametre_empirique_NC(nbre_rejet,5)=variance(donnee_essai(:,1))! sigma_s empirique
                parametre_empirique_NC(nbre_rejet,6)=variance(donnee_essai(:,2))! sigma_t empirique
                parametre_empirique_NC(nbre_rejet,7)=covariance(donnee_essai(:,1),donnee_essai(:,2))/&
                    (dsqrt(variance(donnee_essai(:,1)))*dsqrt(variance(donnee_essai(:,2))))! correlation empirique: rho(x,y)=cov(x,y)/(sigma_x*sigma_y)
                parametre_empirique_NC(nbre_rejet,8)=variance(donnee_essai(:,4))! gamma S empirique

                if(nsim_node(8)==2)then
                    parametre_empirique_NC(nbre_rejet,9)=variance(don_simul(:,w_ijt))! theta t empirique
                    parametre_empirique_NC(nbre_rejet,10)=covariance(don_simul(:,w_ij1),don_simul(:,w_ijt))/&
                    (dsqrt(variance(don_simul(:,w_ij1)))*dsqrt(variance(don_simul(:,w_ijt))))! correlation empirique wij_ST
                    parametre_empirique_NC(nbre_rejet,11)=variance(donnee_essai(:,5))! gamma T empirique
                    parametre_empirique_NC(nbre_rejet,12)=covariance(donnee_essai(:,4),donnee_essai(:,5))/&
                    (dsqrt(variance(donnee_essai(:,4)))*dsqrt(variance(donnee_essai(:,5))))! correlation empirique ui_ST
                endif
            endif

            !!write(12,*)parametre_empirique_NC(nbre_rejet,:)
            if(gener_only.eq.1) then
                theta_sim=variance(don_simul(:,w_ij1))! variance des effets aleatoires simules
                sigmas_sim=variance(donnee_essai(:,1))! sigma_s empirique
                sigmat_sim=variance(donnee_essai(:,2))! sigma_t empirique
                rho_sim=covariance(donnee_essai(:,1),donnee_essai(:,2))/(dsqrt(variance(donnee_essai(:,1)))*&
                            dsqrt(variance(donnee_essai(:,2))))
                gamma_sim=variance(donnee_essai(:,4)) !gamma_ui empirique
                theta_simt=variance(don_simul(:,w_ijt))! variance des effets aleatoires simules T
                rho_sim_wij=covariance(don_simul(:,w_ij1),don_simul(:,w_ijt))/(dsqrt(variance(don_simul(:,w_ij1)))*&
                            dsqrt(variance(don_simul(:,w_ijt))))! correlation empirique wij_ST
                gamma_simt=variance(donnee_essai(:,5))! gamma T empirique
                rho_sim_ui=covariance(donnee_essai(:,4),donnee_essai(:,5))/(dsqrt(variance(donnee_essai(:,4)))*&
                            dsqrt(variance(donnee_essai(:,5))))! correlation empirique ui_ST
                !!print*,"rho=",rho_sim
                !stop
                tab_var_theta(s_i)=theta_sim
                tab_var_sigma(s_i,1)=sigmas_sim
                tab_var_sigma(s_i,2)=sigmat_sim
                tab_var_sigma(s_i,3)=rho_sim
                tab_var_sigma(s_i,4)=gamma_sim
                tab_var_sigma(s_i,5)=theta_simt
                tab_var_sigma(s_i,6)=rho_sim_wij
                tab_var_sigma(s_i,7)=gamma_simt
                tab_var_sigma(s_i,8)=rho_sim_ui
                moy_trt=moy_trt+(sum(don_simul(:,trt1))*100.d0/n_obs)! proportion des traites
                moy_theta=moy_theta+(theta_sim)  ! on calcule la moyenne empirique des theta S
                moy_thetat=moy_thetat+(theta_simt)  ! on calcule la moyenne empirique des theta T
                moy_rho_wij=moy_rho_wij+(rho_sim_wij)  ! on calcule la moyenne empirique des rho wij
                moy_sigmas=moy_sigmas+(sigmas_sim)  ! on calcule la moyenne empirique des sigmas
                moy_sigmat=moy_sigmat+(sigmat_sim)  ! on calcule la moyenne empirique des sigmat
                moy_rho=moy_rho+(rho_sim)  ! on calcule la moyenne empirique des rho
                moy_gamma=moy_gamma+gamma_sim ! on calcule la moyenne empirique des gamma_ui
                moy_gammat=moy_gammat+gamma_simt ! on calcule la moyenne empirique des gamma_ui
                moy_rho_ui=moy_rho_ui+(rho_sim_ui)  ! on calcule la moyenne empirique des rho wij
                moy_pros=moy_pros+(sum(don_simulS(:,statusS1))*100.d0/n_obs) !progression
                moy_dec=moy_dec+(sum(don_simul(:,statusT1))*100.d0/n_obs) ! deces
            endif
        else
            !on conserve le kappa
            if(((kapa_use.eq.2).or.(kapa_use.eq.3).or.(kapa_use.eq.4)).and.(statut_kappa==0)) then ! on conserve ce kappa juste s'il est nouveau
                !on recupere tous les kappa qui marchent
                kappa(incre_kappa,1)=k0(1)
                kappa(incre_kappa,2)=k0(2)
                incre_kappa=incre_kappa+1
                ! if(rang_proc==0) !print*,"j'incremente le kappa s_i=",s_i,"incre_kappa=",incre_kappa
            endif
            nva=nva1+nva2
            rangparam=np-nva-nparamfrail+indice_eta+indice_theta-nknots-splines_ord
            rangparam_theta=rangparam
            if(indice_eta==1) rangparam_eta=np-nva-nparamfrail+indice_eta-&
                            nknots-splines_ord
            rangparam_sigs=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS-&
                            nknots-splines_ord
            rangparam_sigt=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+&
                            indice_varT-nknots-splines_ord
            rangparam_sigst=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+&
                            indice_varT+indice_covST-nknots-splines_ord
            rangparam_alpha=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+&
                            indice_varT+indice_covST+indice_alpha-nknots-splines_ord
            rangparam_gamma=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+&
                            indice_varT+indice_covST+indice_alpha+indice_gamma-nknots-splines_ord
            rangparam_copula=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+&
                             indice_varT+indice_covST+indice_alpha+indice_gamma + 1-nknots-splines_ord
            if(nsim_node(8)==2)then !effets aleatoires correles
                rangparam_thetat=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+&
                    indice_gamma+indice_theta_t
                rangparam_thetast=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+&
                    indice_gamma+indice_theta_t+indice_theta_st
                rangparam_gammat=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+&
                    indice_gamma+indice_theta_t+indice_theta_st+indice_gamma_t
                rangparam_gammast=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+&
                    indice_gamma+indice_theta_t+indice_theta_st+indice_gamma_t+indice_gamma_st
            endif
            rangparam2=rangparam
            ! elements de calcul des estimations
            theta_sim=variance(don_simul(:,w_ij1))! variance des effets aleatoires simules
            tab_var_theta(s_i-nbre_rejet)=theta_sim
            sigmas_sim=variance(donnee_essai(:,1))! sigma_s empirique
            sigmat_sim=variance(donnee_essai(:,2))! sigma_t empirique
            theta_simt=variance(don_simul(:,w_ijt))! variance des effets aleatoires simules T
            rho_sim_wij=covariance(don_simul(:,w_ij1),don_simul(:,w_ijt))/(dsqrt(variance(don_simul(:,w_ij1)))&
                    *dsqrt(variance(don_simul(:,w_ijt))))! correlation empirique wij_ST
            gamma_simt=variance(donnee_essai(:,5))! gamma T empirique
            rho_sim_ui=covariance(donnee_essai(:,4),donnee_essai(:,5))/(dsqrt(variance(donnee_essai(:,4)))*&
                    dsqrt(variance(donnee_essai(:,5))))! correlation empirique ui_ST
            rho_sim=covariance(donnee_essai(:,1),donnee_essai(:,2))/&
                (dsqrt(variance(donnee_essai(:,1)))*dsqrt(variance(donnee_essai(:,2))))
            gamma_sim=variance(donnee_essai(:,4)) !gamma_ui empirique
            tab_var_sigma(s_i-nbre_rejet,1)=sigmas_sim
            tab_var_sigma(s_i-nbre_rejet,2)=sigmat_sim
            tab_var_sigma(s_i-nbre_rejet,3)=rho_sim
            tab_var_sigma(s_i-nbre_rejet,4)=gamma_sim
            tab_var_sigma(s_i-nbre_rejet,5)=theta_simt
            tab_var_sigma(s_i-nbre_rejet,6)=rho_sim_wij
            tab_var_sigma(s_i-nbre_rejet,7)=gamma_simt
            tab_var_sigma(s_i-nbre_rejet,8)=rho_sim_ui
            moy_trt=moy_trt+(sum(don_simul(:,trt1))*100.d0/n_obs)! proportion des traites
            moy_theta=moy_theta+(theta_sim)  ! on calcule la moyenne empirique des theta
            moy_sigmas=moy_sigmas+(sigmas_sim)  ! on calcule la moyenne empirique des sigmas
            moy_sigmat=moy_sigmat+(sigmat_sim)  ! on calcule la moyenne empirique des sigmat
            moy_rho=moy_rho+(rho_sim)  ! on calcule la moyenne empirique des rho
            moy_gamma=moy_gamma+gamma_sim ! on calcule la moyenne empirique des gamma_ui
            moy_thetat=moy_thetat+(theta_simt)  ! on calcule la moyenne empirique des theta T
            moy_rho_wij=moy_rho_wij+(rho_sim_wij)  ! on calcule la moyenne empirique des rho wij
            moy_gammat=moy_gammat+gamma_simt ! on calcule la moyenne empirique des gamma_ui
            moy_rho_ui=moy_rho_ui+(rho_sim_ui)  ! on calcule la moyenne empirique des rho wij
            !!print*,statusS1,sum(don_simulS(:,statusS1))
            !stop
            moy_pros=moy_pros+(sum(don_simulS(:,statusS1))*100.d0/n_obs) !progression
            moy_dec=moy_dec+(sum(don_simul(:,statusT1))*100.d0/n_obs) ! deces

            moy_theta_est=moy_theta_est+(b(rangparam)**2.d0) ! theta estime
            moy_se_theta=moy_se_theta+(dsqrt(((2.d0*b(rangparam))**2.d0)*H_hessOut(rangparam,rangparam)))
            bi_theta = b(rangparam)**2.d0 - 1.96d0*(dsqrt(((2.d0*b(rangparam))**2.d0)*&
                        H_hessOut(rangparam,rangparam)))
            bs_theta = b(rangparam)**2.d0 + 1.96d0*(dsqrt(((2.d0*b(rangparam))**2.d0)*&
                        H_hessOut(rangparam,rangparam)))
            !taux de couverture
            if(theta>=bi_theta .and. theta<=bs_theta)then ! taux de couverture
                taux_couverture_theta=taux_couverture_theta+1.d0
            endif
            if(nsim_node(8)==3) then
                if(copula_function == 1) then  ! claton: exp transform
                    moy_param_cop = moy_param_cop+dexp(b(rangparam_copula)) ! param copule estime
                    ! par delta methode, SE = sqrt(exp(theta_chap) * H_theta_chap * exp(theta_chap))
                    pour_ic = dexp(b(rangparam_copula)) * dsqrt(H_hessOut(rangparam_copula,rangparam_copula))
                    bi_param_cop = dexp(b(rangparam_copula)) - 1.96d0 * pour_ic
                    bs_param_cop = dexp(b(rangparam_copula)) + 1.96d0 * pour_ic
                endif
                if(copula_function == 2) then  ! Gumbel: choleschy transform
                    moy_param_cop = moy_param_cop+(b(rangparam_copula)**2.d0) ! param copule estime
                    pour_ic = 2.d0*b(rangparam_copula) * dsqrt(H_hessOut(rangparam_copula,rangparam_copula))
                    bi_param_cop = b(rangparam_copula)**2.d0 - 1.96d0 * pour_ic
                    bs_param_cop = b(rangparam_copula)**2.d0 + 1.96d0 * pour_ic
                endif
                moy_se_param_cop=moy_se_param_cop + pour_ic
                !taux de couverture
                if(thetacopule>=bi_param_cop .and. thetacopule<=bs_param_cop)then ! taux de couverture
                    taux_couverture_param_cop = taux_couverture_param_cop+1.d0
                endif
            endif
            if(nsim_node(8)==2) then
                !theta estime. j'utilise les variables au niveau essai juste pour ne pas avoir a declarer de nouvelles et pas pour grand chose
                if(logNormal==1)then
                    thetaS1 = b(rangparam)
                    thetaT1 = b(rangparam_thetat)
                    thetaST1 =b(rangparam_thetast)
                    ! pour eviter d'avoir des matrices de variances-covariances non defini positive, je suppose que c'est la cholesky qui est generee. par consequent sigma=Chol*Chol^T
                    !Chol: matrice triangulaire inferieur
                    Chol=0.d0
                    Chol(1,1)=thetaS1
                    Chol(2,1)=thetaST1
                    !Chol(1,2)=covST1
                    Chol(2,2)=thetaT1
                    theta_st=MATMUL(Chol,TRANSPOSE(Chol))
                    thetaS_es=theta_st(1,1)
                    thetaT_es=theta_st(2,2)
                    thetaST_es=theta_st(1,2)

                    !theta_t
                    moy_thetat_est=moy_thetat_est+thetaT_es ! sigma_t estime
                    moy_se_thetat=moy_se_thetat+2.d0*dsqrt(thetaST1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast)+&
                                2.d0*thetaT1*thetaST1*H_hessOut(rangparam_thetast,rangparam_thetat)+&
                                thetaT1**2.d0*H_hessOut(rangparam_thetat,rangparam_thetat))
                    bi_sigmat = thetaT_es - 1.96d0*2.d0*dsqrt(thetaST1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast)+&
                                2.d0*thetaT1*thetaST1*H_hessOut(rangparam_thetast,rangparam_thetat)+&
                                thetaT1**2.d0*H_hessOut(rangparam_thetat,rangparam_thetat))
                    bs_sigmat = thetaT_es + 1.96d0*2.d0*dsqrt(thetaST1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast)+&
                                2.d0*thetaT1*thetaST1*H_hessOut(rangparam_thetast,rangparam_thetat)+&
                                thetaT1**2.d0*H_hessOut(rangparam_thetat,rangparam_thetat))
                    !taux de couverture

                    if(theta2_t>=bi_sigmat .and. theta2_t<=bs_sigmat)then ! taux de couverture
                        taux_couverture_thetat=taux_couverture_thetat+1.d0
                    endif

                    !sigma_st
                    moy_thetast_est=moy_thetast_est+thetaST_es ! sigma_t estime
                    moy_se_thetast=moy_se_thetast+dsqrt(thetaST1**2.d0*H_hessOut(rangparam,rangparam)+&
                                2.d0*thetaS1*thetaST1*H_hessOut(rangparam,rangparam_thetast)+&
                                thetaS1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast))
                    bi_sigmast = thetaST_es - 1.96d0*dsqrt(thetaST1**2.d0*H_hessOut(rangparam,rangparam)+&
                                2.d0*thetaS1*thetaST1*H_hessOut(rangparam,rangparam_thetast)+&
                                thetaS1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast))
                    bs_sigmast = thetaST_es + 1.96d0*dsqrt(thetaST1**2.d0*H_hessOut(rangparam,rangparam)+&
                                2.d0*thetaS1*thetaST1*H_hessOut(rangparam,rangparam_thetast)+&
                                thetaS1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast))
                    !taux de couverture
                    !sigmast_vrai=rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
                    if(thetast_vrai>=bi_sigmast .and. thetast_vrai<=bs_sigmast)then ! taux de couverture
                        taux_couverture_thetast=taux_couverture_thetast+1.d0
                    endif
                endif
            endif
            moy_ni=moy_ni+ni ! Nombre moyen d'itteration pour la convergence
            if(nsim_node(8).ne.0)then !model conjoint surrogate
                ! gamma estime
                if(frailt_base==1) then
                    ! gamma estime S
                    moy_gamma_est=moy_gamma_est+(b(rangparam_gamma)**2.d0) ! theta estime
                    moy_se_gamma=moy_se_gamma+(dsqrt(((2.d0*b(rangparam_gamma))**2.d0)*&
                                 H_hessOut(rangparam_gamma,rangparam_gamma)))
                    bi_gamma = b(rangparam_gamma)**2.d0 - 1.96d0*(dsqrt(((2.d0*b(rangparam_gamma))**2.d0)*&
                               H_hessOut(rangparam_gamma,rangparam_gamma)))
                    bs_gamma = b(rangparam_gamma)**2.d0 + 1.96d0*(dsqrt(((2.d0*b(rangparam_gamma))**2.d0)*&
                               H_hessOut(rangparam_gamma,rangparam_gamma)))
                    !taux de couverture
                    if(gamma_ui>=bi_gamma .and. gamma_ui<=bs_gamma)then ! taux de couverture
                        taux_couverture_gamma=taux_couverture_gamma+1.d0
                    endif

                    if(nsim_node(8)==2)then
                        gammaS1 = b(rangparam_gamma)
                        gammaT1 = b(rangparam_gammat)
                        gammaST1 =b(rangparam_gammast)
                        ! pour eviter d'avoir des matrices de variances-covariances non defini positive, je suppose que c'est la cholesky qui est generee. par consequent sigma=Chol*Chol^T
                        !Chol: matrice triangulaire inferieur
                        Chol=0.d0
                        Chol(1,1)=gammaS1
                        Chol(2,1)=gammaST1
                        !Chol(1,2)=gammaST1
                        Chol(2,2)=gammaT1
                        gamma_st=MATMUL(Chol,TRANSPOSE(Chol))
                        gammaS_es=gamma_st(1,1)
                        gammaT_es=gamma_st(2,2)
                        gammaST_es=gamma_st(1,2)

                        ! gamma estime T
                        moy_gammat_est=moy_gammat_est+gammaT_es ! sigma_t estime
                        moy_se_gammat=moy_se_gammat+2.d0*dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast)+&
                                    2.d0*gammaT1*gammaST1*H_hessOut(rangparam_gammast,rangparam_gammat)+&
                                    gammaT1**2.d0*H_hessOut(rangparam_gammat,rangparam_gammat))
                        bi_sigmat = gammaT_es - 1.96d0*2.d0*dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast)+&
                                    2.d0*gammaT1*gammaST1*H_hessOut(rangparam_gammast,rangparam_gammat)+&
                                    gammaT1**2.d0*H_hessOut(rangparam_gammat,rangparam_gammat))
                        bs_sigmat = gammaT_es + 1.96d0*2.d0*dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast)+&
                                    2.d0*gammaT1*gammaST1*H_hessOut(rangparam_gammast,rangparam_gammat)+&
                                    gammaT1**2.d0*H_hessOut(rangparam_gammat,rangparam_gammat))
                        !taux de couverture

                        if(gamma_uit>=bi_sigmat .and. gamma_uit<=bs_sigmat)then ! taux de couverture
                            taux_couverture_gammat=taux_couverture_gammat+1.d0
                        endif

                        !gamma_st
                        moy_gammast_est=moy_gammast_est+gammaST_es ! sigma_t estime
                        moy_se_gammast=moy_se_gammast+dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gamma,rangparam_gamma)+&
                                    2.d0*gammaS1*gammaST1*H_hessOut(rangparam_gamma,rangparam_gammast)+&
                                    gammaS1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast))
                        bi_sigmast = gammaST_es - 1.96d0*dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gamma,rangparam_gamma)+&
                                    2.d0*gammaS1*gammaST1*H_hessOut(rangparam_gamma,rangparam_gammast)+&
                                    gammaS1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast))
                        bs_sigmast = gammaST_es + 1.96d0*dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gamma,rangparam_gamma)+&
                                    2.d0*gammaS1*gammaST1*H_hessOut(rangparam_gamma,rangparam_gammast)+&
                                    gammaS1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast))
                        !taux de couverture
                        !sigmast_vrai=rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
                        if(gammast_vrai>=bi_sigmast .and. gammast_vrai<=bs_sigmast)then ! taux de couverture
                            taux_couverture_gammast=taux_couverture_gammast+1.d0
                        endif
                    endif
                endif
                !sigma estimes
                covST_es = 0.d0
                if(logNormal==1)then
                    varS1 = b(rangparam_sigs)
                    varT1 = b(rangparam_sigt)
                    covST1 =b(rangparam_sigst)
                    ! pour eviter d'avoir des matrices de variances-covariances non defini positive, je suppose que c'est la cholesky qui est generee. par consequent sigma=Chol*Chol^T
                    !Chol: matrice triangulaire inferieur
                    Chol=0.d0
                    Chol(1,1)=varS1
                    Chol(2,1)=covST1
                    !Chol(1,2)=covST1
                    Chol(2,2)=varT1
                    sigma_st=MATMUL(Chol,TRANSPOSE(Chol))
                    varS_es=sigma_st(1,1)
                    varT_es=sigma_st(2,2)
                    covST_es=sigma_st(1,2)

                    ! =========Delta methode pour varcov des elements de sigme_v===============
                    ! recherche de la matrice de variance-covariance de (sigma_S,sigma_ST,sigmaT) par la delta methode:
                    ! à partir de la hessienne. voir le raisonnement dans le cahier à la date du 04/01/2019
                    hb(1,:) = (/ 2.d0*Chol(1,1), 0.d0, 0.d0 /)
                    hb(2,:) = (/ 0.d0, 2.d0*Chol(2,2), 2.d0*Chol(2,1) /)
                    hb(3,:) = (/ Chol(2,1), 0.d0, Chol(1,1) /)
                    sigmac(1,:) = (/H_hessOut(rangparam_sigs,rangparam_sigs), H_hessOut(rangparam_sigs,rangparam_sigt), &
                                    H_hessOut(rangparam_sigs,rangparam_sigst)/)
                    sigmac(2,:) = (/H_hessOut(rangparam_sigt,rangparam_sigs), H_hessOut(rangparam_sigt,rangparam_sigt), &
                                    H_hessOut(rangparam_sigt,rangparam_sigst)/)
                    sigmac(3,:) = (/H_hessOut(rangparam_sigst,rangparam_sigs), H_hessOut(rangparam_sigst,rangparam_sigt), &
                                    H_hessOut(rangparam_sigst,rangparam_sigst)/)
                    hb = TRANSPOSE(hb)
                    varcov = MATMUL(TRANSPOSE(hb), sigmac)
                    varcov = MATMUL(varcov, hb)


                    ! ========== Fin delta methode ==================

                    ! ====sauvegarde de la hessienne et du vecteur b des parametres====
                    do i=1,np
                        dataHessian(np*(s_i-nbre_rejet-1) + i,:) = H_hessOut(i,:)
                        dataHessianIH(np*(s_i-nbre_rejet-1) + i,:) = HIHOut(i,:)
                    enddo
                    datab(s_i-nbre_rejet,:) = b
                    R2_trial=(covST1**2)/(covST1**2+varT1**2)
                    se_R2_trial=2.d0*dsqrt((covST1**2 * varT1**4 * H_hessOut(rangparam_sigst,rangparam_sigst)-2.d0*covST1**3 &
                    * varT1**3 *H_hessOut(rangparam_sigst,rangparam_sigt) + varT1**2 * covST1**4 * H_hessOut(rangparam_sigt,&
                    rangparam_sigt)))/(covST1**2 + varT1**2)**2
                    moy_R2_trial=moy_R2_trial+R2_trial
                    moy_se_R2_trial=moy_se_R2_trial+se_R2_trial
                    bi_R2_trial = R2_trial - 1.96d0*se_R2_trial
                    bs_R2_trial = R2_trial + 1.96d0*se_R2_trial
                    moy_bi_R2_trial=moy_bi_R2_trial+bi_R2_trial
                    moy_bs_R2_trial=moy_bs_R2_trial+bs_R2_trial
                    !taux de couverture
                    if((rsqrt**2)>=bi_R2_trial .and. (rsqrt**2)<=bs_R2_trial)then ! taux de couverture
                        taux_couverture_R2_trial=taux_couverture_R2_trial+1.d0
                    endif
                endif
                !sigma_s
                moy_sigmas_est=moy_sigmas_est+varS_es ! sigma_s estime
                moy_se_sigmas=moy_se_sigmas+ dsqrt(varcov(1,1))
                bi_sigmas = varS_es - 1.96d0*dsqrt(varcov(1,1))
                bs_sigmas = varS_es + 1.96d0*dsqrt(varcov(1,1))
                !taux de couverture
                if(sigma_s>=bi_sigmas .and. sigma_s<=bs_sigmas)then ! taux de couverture
                    taux_couverture_sigmas=taux_couverture_sigmas+1.d0
                endif
                !sigma_t
                moy_sigmat_est=moy_sigmat_est+varT_es ! sigma_t estime
                moy_se_sigmat=moy_se_sigmat+dsqrt(varcov(2,2))
                bi_sigmat = varT_es - 1.96d0*dsqrt(varcov(2,2))
                bs_sigmat = varT_es + 1.96d0*dsqrt(varcov(2,2))
                !taux de couverture

                if(sigma_t>=bi_sigmat .and. sigma_t<=bs_sigmat)then ! taux de couverture
                    taux_couverture_sigmat=taux_couverture_sigmat+1.d0
                endif
                !sigma_st
                moy_sigmast_est=moy_sigmast_est+covST_es ! sigma_t estime
                moy_se_sigmast=moy_se_sigmast+dsqrt(varcov(3,3))
                bi_sigmast = covST_es - 1.96d0*dsqrt(varcov(3,3))
                bs_sigmast = covST_es + 1.96d0*dsqrt(varcov(3,3))
                if(sigmast_vrai>=bi_sigmast .and. sigmast_vrai<=bs_sigmast)then ! taux de couverture
                    taux_couverture_sigmast=taux_couverture_sigmast+1.d0
                endif
            endif
            !eta
            if(indice_eta==1)then
                moy_eta=moy_eta+(b(np-nva-nknots-splines_ord-nparamfrail+indice_eta))
                moy_se_eta=moy_se_eta+(dsqrt(H_hessOut(np-nva-nknots-splines_ord-nparamfrail+indice_eta,&
                                       np-nva-nknots-splines_ord-nparamfrail+indice_eta)))
                bi_eta = b(np-nva-nknots-splines_ord-nparamfrail+indice_eta) - 1.96d0*(dsqrt(H_hessOut(np-nva-&
                           nknots-splines_ord-nparamfrail+indice_eta,np-nva-nknots-splines_ord-nparamfrail+indice_eta)))
                bs_eta = b(np-nva-nknots-splines_ord-nparamfrail+indice_eta) + 1.96d0*(dsqrt(H_hessOut(np-nva-&
                            nknots-splines_ord-nparamfrail+indice_eta,np-nva-nknots-splines_ord-nparamfrail+indice_eta)))
                !taux de couverture
                if(eta>=bi_eta .and. eta<=bs_eta)then ! taux de couverture
                    taux_couverture_eta=taux_couverture_eta+1.d0
                endif
            endif
            !alpha
            if(indice_alpha==1)then
                moy_alpha=moy_alpha+(b(rangparam_alpha))
                moy_se_alpha=moy_se_alpha+(dsqrt(H_hessOut(rangparam_alpha,rangparam_alpha)))
                bi_eta = b(rangparam_alpha) - 1.96d0*(dsqrt(H_hessOut(rangparam_alpha,rangparam_alpha)))
                bs_eta = b(rangparam_alpha) + 1.96d0*(dsqrt(H_hessOut(rangparam_alpha,rangparam_alpha)))
                !taux de couverture
                if(alpha_ui>=bi_eta .and. alpha_ui<=bs_eta)then ! taux de couverture
                    taux_couverture_alpha=taux_couverture_alpha+1.d0
                endif
            endif
            if(nsim_node(8).ne.3) then
                do i=1,nva+nknots+splines_ord
                    rangparam=np-nva-nknots-splines_ord+i
                    !Intervalle de confiance
                    bi2 = b(rangparam) - 1.96d0*dsqrt(H_hessOut(rangparam,rangparam))
                    bs2 = b(rangparam) + 1.96d0*dsqrt(H_hessOut(rangparam,rangparam))

                    if(i.eq.1)then
                        ! if(rang_proc==0) !write(*,*)'*** For surrogate ***'
                        moy_betaS(1)=moy_betaS(1)+b(rangparam)
                        moy_betaS_se(1)=moy_betaS_se(1)+(dsqrt(H_hessOut(rangparam,rangparam)))
                        parametre_estimes(s_i-nbre_rejet,5)=b(rangparam) !beta_s_chap
                        parametre_estimes(s_i-nbre_rejet,6)=dsqrt(H_hessOut(rangparam,rangparam)) !sd beta_s_chap
                        if(betas>=bi2 .and. betas<=bs2)then ! taux de couverture
                            taux_couvertureS(1)=taux_couvertureS(1)+1.d0
                        endif
                    endif

                    if(i.eq.nva1+1)then
                        ! if(rang_proc==0) !write(*,*)'*** For true endpoint ***',nva,i
                        moy_betaT(1)=moy_betaT(1)+b(rangparam)
                        moy_betaT_se(1)=moy_betaT_se(1)+(dsqrt(H_hessOut(rangparam,rangparam)))
                        parametre_estimes(s_i-nbre_rejet,7)=b(rangparam) !beta_t_chap
                        parametre_estimes(s_i-nbre_rejet,8)=(dsqrt(H_hessOut(rangparam,rangparam))) !sd beta_t_chap
                        if(betat>=bi2 .and. betat<=bs2)then ! taux de couverture
                            taux_couvertureT(1)=taux_couvertureT(1)+1.d0
                        endif
                    endif
                end do
            else ! copula model
                do i=1,nva1
                    rangparam=np-nva+i
                    !Intervalle de confiance
                    bi2 = b(rangparam) - 1.96d0*dsqrt(H_hessOut(rangparam,rangparam))
                    bs2 = b(rangparam) + 1.96d0*dsqrt(H_hessOut(rangparam,rangparam))
                    moy_betaS(i)=moy_betaS(i)+b(rangparam)
                    moy_betaS_se(i)=moy_betaS_se(i)+(dsqrt(H_hessOut(rangparam,rangparam)))
                    if(i.eq.1) then ! traitement
                        parametre_estimes(s_i-nbre_rejet,5)=b(rangparam) !beta_s_chap
                        parametre_estimes(s_i-nbre_rejet,6)=dsqrt(H_hessOut(rangparam,rangparam)) !sd beta_s_chap
                    else ! sauvegarde des autres covariables
                        parametre_estimes(s_i-nbre_rejet,size(parametre_estimes,2)-2*(nva2-1)-2*nva1+2*(i-1)+1)=b(rangparam)
                        parametre_estimes(s_i-nbre_rejet,size(parametre_estimes,2)-2*(nva2-1)-2*nva1+2*(i-1)+2)=&
                        dsqrt(H_hessOut(rangparam,rangparam))
                    endif
                    if(betas>=bi2 .and. betas<=bs2)then ! taux de couverture
                        taux_couvertureS(i)=taux_couvertureS(i)+1.d0
                    endif
                enddo
                do i=1,nva2
                    rangparam = np - nva + nva1 + i
                    !Intervalle de confiance
                    bi2 = b(rangparam) - 1.96d0*dsqrt(H_hessOut(rangparam,rangparam))
                    bs2 = b(rangparam) + 1.96d0*dsqrt(H_hessOut(rangparam,rangparam))
                    moy_betaT(i)=moy_betaT(i)+b(rangparam)
                    moy_betaT_se(i)=moy_betaT_se(i)+(dsqrt(H_hessOut(rangparam,rangparam)))
                    if(i.eq.1) then ! traitement
                        parametre_estimes(s_i-nbre_rejet,7)=b(rangparam) !beta_t_chap
                        parametre_estimes(s_i-nbre_rejet,8)=dsqrt(H_hessOut(rangparam,rangparam)) !sd beta_t_chap
                    else ! sauvegarde des autres covariables
                        parametre_estimes(s_i-nbre_rejet,size(parametre_estimes,2)-2*nva2+2*(i-1)+1)=b(rangparam)
                        parametre_estimes(s_i-nbre_rejet,size(parametre_estimes,2)-2*nva2+2*(i-1)+2)=&
                        dsqrt(H_hessOut(rangparam,rangparam))
                    endif
                    if(betat>=bi2 .and. betat<=bs2)then ! taux de couverture
                        taux_couvertureT(i)=taux_couvertureT(i)+1.d0
                    endif
                    !endif
                end do
            endif
            ! parametres empiriques
            parametre_empirique(s_i-nbre_rejet,1)=sum(don_simul(:,trt1))*100.d0/n_obs ! personne traitees
            parametre_empirique(s_i-nbre_rejet,2)=variance(don_simul(:,w_ij1)) ! theta empirique
            parametre_empirique(s_i-nbre_rejet,3)=sum(don_simul(:,statusS1))*100.d0/n_obs ! progression
            parametre_empirique(s_i-nbre_rejet,4)=sum(don_simul(:,statusT1))*100.d0/n_obs
            if(nsim_node(8).ne.0)then !model conjoint surrogate
                parametre_empirique(s_i-nbre_rejet,5)=variance(donnee_essai(:,1))! sigma_s empirique
                parametre_empirique(s_i-nbre_rejet,6)=variance(donnee_essai(:,2))! sigma_t empirique
                parametre_empirique(s_i-nbre_rejet,7)=covariance(donnee_essai(:,1),donnee_essai(:,2))/&
                    (dsqrt(variance(donnee_essai(:,1)))*dsqrt(variance(donnee_essai(:,2))))! correlation empirique: rho(x,y)=cov(x,y)/(sigma_x*sigma_y)
                parametre_empirique(s_i-nbre_rejet,8)=variance(donnee_essai(:,4))! gamma empirique
                if(nsim_node(8)==2)then
                    parametre_empirique(s_i-nbre_rejet,9)=variance(don_simul(:,w_ijt))! theta t empirique
                    parametre_empirique(s_i-nbre_rejet,10)=covariance(don_simul(:,w_ij1),don_simul(:,w_ijt))/&
                        (dsqrt(variance(don_simul(:,w_ij1)))*dsqrt(variance(don_simul(:,w_ijt))))! correlation empirique wij_ST
                    parametre_empirique(s_i-nbre_rejet,11)=variance(donnee_essai(:,5))! gamma T empirique
                    parametre_empirique(s_i-nbre_rejet,12)=covariance(donnee_essai(:,4),donnee_essai(:,5))/&
                        (dsqrt(variance(donnee_essai(:,4)))*dsqrt(variance(donnee_essai(:,5))))! correlation empirique ui_ST
                endif
            endif
            if(nsim_node(8).ne.3) then
                parametre_estimes(s_i-nbre_rejet,1)=(b(rangparam2)**2.d0) !theta
                parametre_estimes(s_i-nbre_rejet,2)=(dsqrt(((2.d0*b(rangparam2))**2.d0)*H_hessOut(rangparam2,rangparam2)))! se theta
            else ! copula
                if(copula_function == 1) then ! claton: exp transform
                    parametre_estimes(s_i-nbre_rejet,1) = dexp(b(rangparam_copula)) !theta_copula
                    parametre_estimes(s_i-nbre_rejet,2) = pour_ic
                endif
                if(copula_function == 2) then ! Gumbel
                    parametre_estimes(s_i-nbre_rejet,1) = b(rangparam_copula)**2.d0 !theta_copula
                    parametre_estimes(s_i-nbre_rejet,2) = pour_ic
                endif
            endif
            if(indice_eta==1)then
                parametre_estimes(s_i-nbre_rejet,3)=(b(np-nva-nknots-splines_ord-nparamfrail+indice_eta))!zeta
                parametre_estimes(s_i-nbre_rejet,4)=(dsqrt(H_hessOut(np-nva-nknots-splines_ord-nparamfrail+&
                                  indice_eta,np-nva-nparamfrail+indice_eta))) !se zeta
            else
                parametre_estimes(s_i-nbre_rejet,3)=1.d0!zeta
                parametre_estimes(s_i-nbre_rejet,4)=0.d0 !se zeta
            endif
            if(nsim_node(8).ne.0)then !model conjoint surrogate (se calcule par la delta methode)
                parametre_estimes(s_i-nbre_rejet,9)= varS_es!siqma_s
                parametre_estimes(s_i-nbre_rejet,10)=dsqrt(varcov(1,1))! se sigma_s
                parametre_estimes(s_i-nbre_rejet,11)= varT_es!siqma_t
                parametre_estimes(s_i-nbre_rejet,12)=dsqrt(varcov(2,2))! se sigma_t
                parametre_estimes(s_i-nbre_rejet,13)= covST_es !siqma_st
                parametre_estimes(s_i-nbre_rejet,14)=dsqrt(varcov(3,3))! se sigma_st
                if(frailt_base==1)then
                    parametre_estimes(s_i-nbre_rejet,15)=(b(rangparam_gamma)**2.d0) !gamma
                    parametre_estimes(s_i-nbre_rejet,16)=(dsqrt(((2.d0*b(rangparam_gamma))**2.d0)*H_hessOut(rangparam_gamma,&
                        rangparam_gamma)))! se gamma
                    if(indice_alpha==1)then
                        parametre_estimes(s_i-nbre_rejet,17)=(b(rangparam_alpha))!alpha
                        parametre_estimes(s_i-nbre_rejet,18)=(dsqrt(H_hessOut(rangparam_alpha,rangparam_alpha))) !se alpha
                    else
                        parametre_estimes(s_i-nbre_rejet,17)=1.d0!zeta
                        parametre_estimes(s_i-nbre_rejet,18)=0.d0 !se zeta
                    endif
                endif
                if(nsim_node(8)==1 .or. nsim_node(8)==3)then
                    parametre_estimes(s_i-nbre_rejet,19)=R2_trial     !r2_trial reduite
                    parametre_estimes(s_i-nbre_rejet,20)=se_R2_trial    !se r2_trial reduite
                endif
                if(nsim_node(8)==2)then
                    parametre_estimes(s_i-nbre_rejet,19)= thetaT_es!theta_t
                    parametre_estimes(s_i-nbre_rejet,20)=2.d0*dsqrt(thetaST1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast)+&
                            2.d0*thetaT1*thetaST1*H_hessOut(rangparam_thetast,rangparam_thetat)+&
                            thetaT1**2.d0*H_hessOut(rangparam_thetat,rangparam_thetat))! se theta_t
                    parametre_estimes(s_i-nbre_rejet,21)= thetaST_es !theta_st
                    parametre_estimes(s_i-nbre_rejet,22)=dsqrt(thetaST1**2.d0*H_hessOut(rangparam_theta,rangparam_theta)+&
                            2.d0*thetaS1*thetaST1*H_hessOut(rangparam_theta,rangparam_thetast)+&
                            thetaS1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast))! se theta_st
                    !gamma
                    parametre_estimes(s_i-nbre_rejet,23)= gammaT_es!gamma_t
                    parametre_estimes(s_i-nbre_rejet,24)=2.d0*dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast)+&
                            2.d0*gammaT1*gammaST1*H_hessOut(rangparam_gammast,rangparam_gammat)+&
                            gammaT1**2.d0*H_hessOut(rangparam_gammat,rangparam_gammat))! se gamma_t
                    parametre_estimes(s_i-nbre_rejet,25)= gammaST_es !gamma_st
                    parametre_estimes(s_i-nbre_rejet,26)=dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gamma,rangparam_gamma)+&
                            2.d0*gammaS1*gammaST1*H_hessOut(rangparam_gamma,rangparam_gammast)+&
                            gammaS1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast))! se gamma_st
                    ! pour les R2 trial reduits et complets
                    parametre_estimes(s_i-nbre_rejet,27)=R2_trial     !r2_trial reduite
                    parametre_estimes(s_i-nbre_rejet,28)=se_R2_trial    !se r2_trial reduite
                endif
                ! ======calcul du taux de kendall =================
                call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation pour le calcul du tau de kendall
                if(nsim_node(8)==2)then
                    if(method_int_kendal==4 .or.method_int_kendal==5)then ! 1 seul taux de kendall et modele complet
                        tau_kendal_00=tau_kendall(theta_ST,gamma_st,sigma_st,0,0,method_int_kendal,N_MC_kendall,&
                        parametre_estimes(s_i-nbre_rejet,17),parametre_estimes(s_i-nbre_rejet,3),1)    !tau de kendal des non traites z_11=0,z_21=0
                        moy_kendal_00=moy_kendal_00+tau_kendal_00

                    else    ! 4 taux de kendall (suivant les bras de traitement) et modele complet
                        tau_kendal_11=tau_kendall(theta_ST,gamma_st,sigma_st,1,1,method_int_kendal,N_MC_kendall,&
                            parametre_estimes(s_i-nbre_rejet,17),parametre_estimes(s_i-nbre_rejet,3),1)    !tau de kendal des traites z_11=1,z_21=1
                        tau_kendal_10=tau_kendall(theta_ST,gamma_st,sigma_st,1,0,method_int_kendal,N_MC_kendall,&
                            parametre_estimes(s_i-nbre_rejet,17),parametre_estimes(s_i-nbre_rejet,3),1)    !tau de kendal des 1 traite et l'autre non traite z_11=1,z_21=0
                        tau_kendal_01=tau_kendall(theta_ST,gamma_st,sigma_st,0,1,method_int_kendal,N_MC_kendall,&
                            parametre_estimes(s_i-nbre_rejet,17),parametre_estimes(s_i-nbre_rejet,3),1)    !tau de kendal des traite z_11=0,z_21=1
                        tau_kendal_00=tau_kendall(theta_ST,gamma_st,sigma_st,0,0,method_int_kendal,N_MC_kendall,&
                            parametre_estimes(s_i-nbre_rejet,17),parametre_estimes(s_i-nbre_rejet,3),1)    !tau de kendal des non traites z_11=0,z_21=0
                        moy_kendal_11=moy_kendal_11+tau_kendal_11
                        moy_kendal_10=moy_kendal_10+tau_kendal_10
                        moy_kendal_01=moy_kendal_01+tau_kendal_01
                        moy_kendal_00=moy_kendal_00+tau_kendal_00
                    endif
                endif
                if(nsim_node(8)==1 .or. nsim_node(8)==3)then ! 1 seul taux de kendall et modele reduit
                    allocate(theta_ST_2(1,1),gamma_st_2(1,1))
                    theta_ST_2(1,1)= parametre_estimes(s_i-nbre_rejet,1) !theta estime
                    gamma_st_2(1,1)= parametre_estimes(s_i-nbre_rejet,15) ! gamma estime
                    if(nsim_node(8)==1) then
                        tau_kendal_00=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,method_int_kendal,N_MC_kendall,&
                        parametre_estimes(s_i-nbre_rejet,17),parametre_estimes(s_i-nbre_rejet,3),0)    !tau de kendal des non traites z_11=0,z_21=0
                    else !copula
                        if(copula_function == 1) then! claton
                            tau_kendal_00 = parametre_estimes(s_i-nbre_rejet,1) / (parametre_estimes(s_i-nbre_rejet,1) +2.d0)
                        endif
                        if(copula_function == 2) then ! Gumbel
                            tau_kendal_00 = parametre_estimes(s_i-nbre_rejet,1) / (parametre_estimes(s_i-nbre_rejet,1) +1.d0)
                        endif
                    endif
                    moy_kendal_00=moy_kendal_00+tau_kendal_00
                    !calcul intervalle de confiance du tau de kendall par bootstrap parametrique, a partir de la distribution a posteriorie des parametres
                    if(nsim_node(8).ne.3)then ! si pas modele de copule
                        if(indice_alpha==0 .and. indice_eta==0)then ! on estime ni alpha, ni eta
                            v_chap_kendall=0.d0
                            t_chap_kendall=(/b(rangparam_theta),b(rangparam_gamma)/) ! parametres necessaire: theta, gamma, zeta, alpha
                            v_chap_kendall(1,:)=(/H_hessOut(rangparam_theta,rangparam_theta),&
                            H_hessOut(rangparam_theta,rangparam_gamma)/)
                            v_chap_kendall(2,:)=(/H_hessOut(rangparam_theta,rangparam_gamma),&
                            H_hessOut(rangparam_gamma,rangparam_gamma)/)
                        else ! on estime au moins un des deux
                            if(indice_alpha==1 .and. indice_eta == 1)then !on estime les deux
                                t_chap_kendall=(/b(rangparam_theta),b(rangparam_gamma),b(rangparam_eta),b(rangparam_alpha)/) ! parametres necessaire: theta, gamma, zeta, alpha
                                v_chap_kendall(1,:)=(/H_hessOut(rangparam_theta,rangparam_theta),H_hessOut(rangparam_theta,&
                                                    rangparam_gamma),H_hessOut(rangparam_theta,rangparam_eta),&
                                                    H_hessOut(rangparam_theta,rangparam_alpha)/)
                                v_chap_kendall(2,:)=(/H_hessOut(rangparam_theta,rangparam_gamma),H_hessOut(rangparam_gamma,&
                                                    rangparam_gamma),H_hessOut(rangparam_eta,rangparam_gamma),&
                                                    H_hessOut(rangparam_alpha,rangparam_gamma)/)
                                v_chap_kendall(3,:)=(/H_hessOut(rangparam_theta,rangparam_eta),H_hessOut(rangparam_eta,&
                                                    rangparam_gamma),H_hessOut(rangparam_eta,rangparam_eta),&
                                                    H_hessOut(rangparam_eta,rangparam_alpha)/)
                                v_chap_kendall(4,:)=(/H_hessOut(rangparam_theta,rangparam_alpha),&
                                                    H_hessOut(rangparam_alpha,rangparam_gamma) ,H_hessOut(rangparam_alpha,&
                                                    rangparam_eta),H_hessOut(rangparam_alpha,rangparam_alpha)/)
                            else ! on estime seulement un des deux
                                if(indice_alpha==1)then ! c'est alpha on estime
                                    v_chap_kendall=0.d0
                                    t_chap_kendall=(/b(rangparam_theta),b(rangparam_gamma),b(rangparam_alpha)/) ! parametres necessaire: theta, gamma, alpha
                                    v_chap_kendall(1,:)=(/H_hessOut(rangparam_theta,rangparam_theta),&
                                                         H_hessOut(rangparam_theta,rangparam_gamma)&
                                                        ,H_hessOut(rangparam_theta,rangparam_alpha)/)
                                    v_chap_kendall(2,:)=(/H_hessOut(rangparam_theta,rangparam_gamma),&
                                                        H_hessOut(rangparam_gamma,rangparam_gamma)&
                                                        ,H_hessOut(rangparam_alpha,rangparam_gamma)/)
                                    v_chap_kendall(3,:)=(/H_hessOut(rangparam_theta,rangparam_alpha),&
                                                         H_hessOut(rangparam_gamma,rangparam_alpha)&
                                                        ,H_hessOut(rangparam_alpha,rangparam_alpha)/)
                                else ! c'est eta on estime
                                    t_chap_kendall=(/b(rangparam_theta),b(rangparam_gamma),b(rangparam_eta)/) ! parametres necessaire: theta, gamma, zeta
                                    v_chap_kendall(1,:)=(/H_hessOut(rangparam_theta,rangparam_theta),&
                                                        H_hessOut(rangparam_theta,rangparam_gamma)&
                                                        ,H_hessOut(rangparam_theta,rangparam_eta)/)
                                    v_chap_kendall(2,:)=(/H_hessOut(rangparam_theta,rangparam_gamma),&
                                                        H_hessOut(rangparam_gamma,rangparam_gamma),&
                                                        H_hessOut(rangparam_eta,rangparam_gamma)/)
                                    v_chap_kendall(3,:)=(/H_hessOut(rangparam_theta,rangparam_eta),&
                                                        H_hessOut(rangparam_eta,rangparam_gamma),&
                                                        H_hessOut(rangparam_eta,rangparam_eta)/)
                                endif
                            endif
                        endif
                    endif
                    ! pour R2
                    v_chap_R2=0.d0
                    t_chap_R2=(/b(rangparam_sigs),b(rangparam_sigt),b(rangparam_sigst)/) ! parametres necessaire: sigma_s,sigma_t,sigma_st
                    v_chap_R2(1,:)=(/H_hessOut(rangparam_sigs,rangparam_sigs),H_hessOut(rangparam_sigs,&
                                    rangparam_sigt), H_hessOut(rangparam_sigs,rangparam_sigst)/)
                    v_chap_R2(2,:)=(/H_hessOut(rangparam_sigt,rangparam_sigs),H_hessOut(rangparam_sigt,&
                                    rangparam_sigt), H_hessOut(rangparam_sigt,rangparam_sigst)/)
                    v_chap_R2(3,:)=(/H_hessOut(rangparam_sigst,rangparam_sigs),H_hessOut(rangparam_sigst,&
                                    rangparam_sigt), H_hessOut(rangparam_sigst,rangparam_sigst)/)

                    moy_tau_boots=0.d0
                    moy_R2_boots=0.d0
                    call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation
                    allocate(matrice_generation(1,np))
                    v_chap_copula(1,1) = H_hessOut(rangparam_copula,rangparam_copula)
                    do i=1,nboot_kendal
                        if(nsim_node(8).ne.3) then
                            call rmvnorm(t_chap_kendall,v_chap_kendall,1,0,theta_chap_kendall)
                            theta_ST_2(1,1)= theta_chap_kendall(1,1)**2.d0 !theta simule
                            gamma_st_2(1,1)= theta_chap_kendall(1,2)**2.d0  ! gamma simule
                        else
                            call rmvnorm((/b(rangparam_copula)/),v_chap_copula,1,0,theta_chap_copula)
                        endif
                        if(indice_alpha==0 .and. indice_eta==0)then
                            if(nsim_node(8).ne.3) then
                                vect_kendall_tau(i)=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,&
                                                     method_int_kendal,N_MC_kendall, 1.d0,1.d0,0)    !tau de kendal des non traites z_11=0,z_21=0
                            else
                                if(copula_function == 1) vect_kendall_tau(i) = dexp(theta_chap_copula(1,1)) / &
                                    (dexp(theta_chap_copula(1,1)) +2.d0) ! claton
                                if(copula_function == 2) vect_kendall_tau(i) = theta_chap_copula(1,1)**2.d0 / &
                                    (theta_chap_copula(1,1)**2.d0 +1.d0)  ! Gumbel
                            endif
                        else
                            if(indice_alpha==1 .and. indice_eta==1)then
                                if(nsim_node(8).ne.3) then
                                    vect_kendall_tau(i)=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,&
                                                        method_int_kendal,N_MC_kendall,&
                                                        theta_chap_kendall(1,4),theta_chap_kendall(1,3),0)
                                else
                                    if(copula_function == 1) vect_kendall_tau(i) = dexp(theta_chap_copula(1,1)) / &
                                        (dexp(theta_chap_copula(1,1)) +2.d0) ! claton
                                    if(copula_function == 2) vect_kendall_tau(i) = theta_chap_copula(1,1)**2.d0 / &
                                        (theta_chap_copula(1,1)**2.d0 +1.d0)  ! Gumbel
                                endif
                            else
                                if(indice_alpha==1)then
                                    if(nsim_node(8).ne.3) then
                                        vect_kendall_tau(i)=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,&
                                                        method_int_kendal,N_MC_kendall,&
                                                        theta_chap_kendall(1,4),1.d0,0)
                                    else
                                        if(copula_function == 1) vect_kendall_tau(i) = dexp(theta_chap_copula(1,1)) / &
                                            (dexp(theta_chap_copula(1,1)) +2.d0) ! claton
                                        if(copula_function == 2) vect_kendall_tau(i) = theta_chap_copula(1,1)**2.d0 / &
                                            (theta_chap_copula(1,1)**2.d0 +1.d0)  ! Gumbel
                                    endif
                                else
                                    if(nsim_node(8).ne.3) then
                                        vect_kendall_tau(i)=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,&
                                                        method_int_kendal,N_MC_kendall,&
                                                        1.d0,theta_chap_kendall(1,3),0)
                                    else
                                        if(copula_function == 1) vect_kendall_tau(i) = dexp(theta_chap_copula(1,1)) / &
                                            (dexp(theta_chap_copula(1,1)) +2.d0) ! claton
                                        if(copula_function == 2) vect_kendall_tau(i) = theta_chap_copula(1,1)**2.d0 / &
                                            (theta_chap_copula(1,1)**2.d0 +1.d0)  ! Gumbel
                                    endif
                                endif
                            endif
                        endif
                        moy_tau_boots=moy_tau_boots+vect_kendall_tau(i)
                        call rmvnorm(t_chap_R2,v_chap_R2,1,0,theta_chap_R2)
                        Chol_R2=0.d0
                        Chol_R2(1,1)=theta_chap_R2(1,1) ! associe a sigma S
                        Chol_R2(2,2)=theta_chap_R2(1,2) ! associe a sigma T
                        Chol_R2(2,1)=theta_chap_R2(1,3) ! associe a sigma ST
                        mat_A=MATMUL(Chol_R2,TRANSPOSE(Chol_R2))
                        vect_R2(i)=(Chol_R2(2,1)**2.d0)/(Chol_R2(2,1)**2.d0+Chol_R2(2,2)**2.d0)
                        moy_R2_boots=moy_R2_boots+vect_R2(i)
                        !!write(18,*)vect_R2(i)
                    enddo
                    deallocate(matrice_generation)
                    ! endif
                    call percentile_scl(vect_kendall_tau,nboot_kendal,0.025d0,IC_Inf) !borne inf
                    call percentile_scl(vect_kendall_tau,nboot_kendal,0.975d0,IC_sup) !borne sup
                    !R2
                    call percentile_scl(vect_R2,nboot_kendal,0.025d0,IC_Inf_R2) !borne inf
                    call percentile_scl(vect_R2,nboot_kendal,0.975d0,IC_sup_R2) !borne sup
                    ! sauvegarde des resultats dans un tableau
                    result_bootstrap(s_i,1:3)=(/moy_R2_boots/dble(nboot_kendal),IC_Inf_R2,IC_sup_R2/) ! R2
                    result_bootstrap(s_i,4:6)=(/moy_tau_boots/dble(nboot_kendal),IC_Inf,IC_sup/)    ! Tau de kendall
                    deallocate(gamma_st_2)
                    deallocate(theta_ST_2)
                endif
                if(rang_proc==0) then
                    !write(3,*)moy_tau_boots/dble(nboot_kendal),IC_Inf,IC_sup
                    fichier_kendall(s_i,:)=(/moy_tau_boots/dble(nboot_kendal),IC_Inf,IC_sup/)
                endif
                if(rang_proc==0) then
                    !write(18,*) moy_R2_boots/dble(nboot_kendal),IC_Inf_R2,IC_sup_R2
                    fichier_R2(s_i,:)=(/moy_R2_boots/dble(nboot_kendal),IC_Inf_R2,IC_sup_R2/)
                endif
                if(nsim_node(8)==1)then
                    parametre_estimes(s_i-nbre_rejet,21)=tau_kendal_11    !tau de kendal des traites z_11=1,z_21=1
                    parametre_estimes(s_i-nbre_rejet,22)=tau_kendal_10    !tau de kendal des 1 traite et l'autre non traite z_11=1,z_21=0
                    parametre_estimes(s_i-nbre_rejet,23)=tau_kendal_01    !tau de kendal des traite z_11=0,z_21=1
                    parametre_estimes(s_i-nbre_rejet,24)=tau_kendal_00    !tau de kendal des non traites z_11=0,z_21=0
                endif
                if(nsim_node(8)==2)then
                    parametre_estimes(s_i-nbre_rejet,29)=tau_kendal_11    !tau de kendal des traites z_11=1,z_21=1
                    parametre_estimes(s_i-nbre_rejet,30)=tau_kendal_10    !tau de kendal des 1 traite et l'autre non traite z_11=1,z_21=0
                    parametre_estimes(s_i-nbre_rejet,31)=tau_kendal_01    !tau de kendal des traite z_11=0,z_21=1
                    parametre_estimes(s_i-nbre_rejet,32)=tau_kendal_00    !tau de kendal des non traites z_11=0,z_21=0
                endif
                if(nsim_node(8)==3)then
                    if(nva1 > 1 .or. nva2 >1) then !si plus d'une variable explicative
                        parametre_estimes(s_i-nbre_rejet,22)=tau_kendal_00    !tau de kendal des non traites z_11=0,z_21=0
                    else
                        parametre_estimes(s_i-nbre_rejet,24)=tau_kendal_00    !tau de kendal des non traites z_11=0,z_21=0
                    endif
                    ! SE de kendall_tau par delta method. voir cahier le 18/04/2019 pour demonstration
                    if(copula_function == 1) then
                        if(nva1 > 1 .or. nva2 >1) then !si plus d'une variable explicative
                            parametre_estimes(s_i-nbre_rejet,23) =  2.d0 * dexp(b(rangparam_copula))/&
                                (dexp(b(rangparam_copula)) + 2.d0)**2.d0 * dsqrt(H_hessOut(rangparam_copula,&
                                rangparam_copula))! se_tau_kendal
                        else
                            parametre_estimes(s_i-nbre_rejet,25) =  2.d0 * dexp(b(rangparam_copula))/&
                                (dexp(b(rangparam_copula)) + 2.d0)**2.d0 * dsqrt(H_hessOut(rangparam_copula,&
                                rangparam_copula))! se_tau_kendal
                        endif
                        bi_sigmas = dexp(b(rangparam_copula)) - 1.96d0*parametre_estimes(s_i-nbre_rejet,25)
                        bs_sigmas = dexp(b(rangparam_copula)) + 1.96d0*parametre_estimes(s_i-nbre_rejet,25)
                        !taux de couverture
                        vrai_tau_copula = thetacopule/(thetacopule + 2.d0)
                    endif
                    if(copula_function == 2) then
                        if(nva1 > 1 .or. nva2 >1) then !si plus d'une variable explicative
                            parametre_estimes(s_i-nbre_rejet,23) =   2.d0 * b(rangparam_copula)/&
                                (b(rangparam_copula)**2.d0 + 1.d0)**2.d0 * dsqrt(H_hessOut(rangparam_copula,&
                                rangparam_copula))! se_tau_kendal
                        else
                            parametre_estimes(s_i-nbre_rejet,25) =   2.d0 * b(rangparam_copula)/&
                                (b(rangparam_copula)**2.d0 + 1.d0)**2.d0 * dsqrt(H_hessOut(rangparam_copula,&
                                rangparam_copula))! se_tau_kendal
                        endif
                        bi_sigmas = b(rangparam_copula)**2.d0 - 1.96d0*parametre_estimes(s_i-nbre_rejet,25)
                        bs_sigmas = b(rangparam_copula)**2.d0  + 1.96d0*parametre_estimes(s_i-nbre_rejet,25)
                        !taux de couverture
                        vrai_tau_copula = thetacopule/(thetacopule + 1.d0)
                    endif
                    if(vrai_tau_copula >= bi_sigmas .and. vrai_tau_copula <= bs_sigmas)then ! taux de couverture
                        taux_couverture_tauk=taux_couverture_tauk+1.d0
                    endif
                endif
            endif
            erreur_fichier=0
            if(nb_processus<=1)then ! on ecrit le resultat dans le fichier que si l'on n'a pas plusieurs processus. si oui l'ecriture est faite apres la sortie de la boucle par le processus de rang 0
                997 continue
                !write(11,*,iostat=erreur_fichier)parametre_estimes(s_i-nbre_rejet,:)
                param_estimes(s_i,:) = parametre_estimes(s_i-nbre_rejet,:)
                if(erreur_fichier .ne.0) then
                    goto 997
                endif
            endif
            if(nb_processus>1)then ! on ne le fait que dans on a plus d'un processus
                parametre_estimes_MPI(indice_sim_proc,:)=parametre_estimes(s_i-nbre_rejet,:) ! juste pour le MIP
                indice_sim_proc=indice_sim_proc+1
            endif
        endif
        1002 continue
        s_i=s_i+1
        !deallocate(b)
        !deallocate(H_hessOut,HIHOut)
        deallocate(vdeces,vsurrogate,paGH,don_simulS,tableEssai)
    end do
    EPS2 = EPS ! on restitue les critères de convergence pour le dernier jeux de donnees
    ! ecart-type des parametres simules
    se_theta_sim=dsqrt(variance(tab_var_theta))
    if(nsim_node(8).ne.0)then
        se_sigmas_sim= dsqrt(variance(tab_var_sigma(:,1)))
        se_sigmat_sim= dsqrt(variance(tab_var_sigma(:,2)))
        se_rho_sim= dsqrt(variance(tab_var_sigma(:,3)))
        se_gamma_sim=dsqrt(variance(tab_var_sigma(:,4)))
        if(nsim_node(8)==2)then
            se_theta_simt=dsqrt(variance(tab_var_sigma(:,5)))
            se_rho_sim_wij=dsqrt(variance(tab_var_sigma(:,6)))
            se_gamma_simt=dsqrt(variance(tab_var_sigma(:,7)))
            se_rho_sim_ui=dsqrt(variance(tab_var_sigma(:,8)))
        endif
    endif
    ! ecart-type des parametres estime(SD)
    se_theta_est=dsqrt(variance(parametre_estimes(:,1)))
    if(indice_eta==1) se_eta_est=dsqrt(variance(parametre_estimes(:,3)))
    se_beta_s=dsqrt(variance(parametre_estimes(:,5)))
    se_beta_t=dsqrt(variance(parametre_estimes(:,7)))
    if(nsim_node(8).ne.0)then
        se_sigmas_est= dsqrt(variance(parametre_estimes(:,9)))
        se_sigmat_est= dsqrt(variance(parametre_estimes(:,11)))
        se_cov_est= dsqrt(variance(parametre_estimes(:,13)))
        se_gamma_est=dsqrt(variance(parametre_estimes(:,15)))
        if(indice_alpha==1) se_alpha_est=dsqrt(variance(parametre_estimes(:,17)))
        if(nsim_node(8)==1)then
            se_R2_trial=dsqrt(variance(parametre_estimes(:,19)))
            se_kendal_11=dsqrt(variance(parametre_estimes(:,21)))
            se_kendal_10=dsqrt(variance(parametre_estimes(:,22)))
            se_kendal_01=dsqrt(variance(parametre_estimes(:,23)))
            se_kendal_00=dsqrt(variance(parametre_estimes(:,24)))
        endif
        if(nsim_node(8)==2)then
            se_thetat_est=dsqrt(variance(parametre_estimes(:,19)))
            se_cov_est_wij=dsqrt(variance(parametre_estimes(:,21)))
            se_gamma_estt=dsqrt(variance(parametre_estimes(:,23)))
            se_cov_est_ui=dsqrt(variance(parametre_estimes(:,25)))
            se_R2_trial=dsqrt(variance(parametre_estimes(:,27)))
            se_kendal_11=dsqrt(variance(parametre_estimes(:,29)))
            se_kendal_10=dsqrt(variance(parametre_estimes(:,30)))
            se_kendal_01=dsqrt(variance(parametre_estimes(:,31)))
            se_kendal_00=dsqrt(variance(parametre_estimes(:,32)))
        endif
    endif
    if(nb_processus>1)then ! on ne le fait que dans on a plus d'un processus
        n_sim_exact=dble(n_sim-nbre_rejet)
        if (rang==0) then
            init_i=1 ! ce processus commence a la premiere simulation
        else
            init_i=sum(tableNsim(1:rang))+1 ! ce processus commence a la simulation qui respecte son ordre et doit s'arreter au nombre de simultation dont il a le droit d'executer
        endif
        max_i=init_i+tableNsim(rang+1)-1
        ! ce qui suit pour l'envoie des element de la matrice vecteur par vecteur
        allocate(tampon(tableNsim(rang+1)*size(parametre_estimes_MPI,2)))
        allocate(tampon_all(n_sim_total*size(parametre_estimes_MPI,2)))
        !sauvegarde des element de la matrice ligne par ligne
        k=1
        do i=1,tableNsim(rang+1)
        !do i=1,n_sim_exact
            j=k+size(parametre_estimes_MPI,2)
            tampon(k:j)=parametre_estimes_MPI(i,:)
            k=k+size(parametre_estimes_MPI,2)
        enddo
        if(rang==0)then
            n_sim_exact=n_sim_exact_0
            ! if(rang_proc==0) !print*, "taux de convergence=",n_sim_exact*100.d0/n_sim_total,"(",n_sim_exact,"/",n_sim_total,")"
            !reccuperation des element sauvegardes dans la matrice parametre_estimes_MPI_T ligne par ligne
            k=1
            !do i=1,n_sim_total
            i=1
            do while (i<=n_sim_total)
            !do i=1,max_i ! avec max_i le nombre de simulation qui a converge
                j=k+size(parametre_estimes_MPI,2)
                    parametre_estimes_MPI_T(i,:)=tampon_all(k:j)

                    i=i+1
                !endif
                k=k+size(parametre_estimes_MPI,2)
            enddo
            ! recherche des lignes correspondantes aux simulations qui n'ont pas convergees
            deallocate(parametre_estimes)
            allocate(parametre_estimes(INT(n_sim_exact),size(parametre_estimes_MPI,2)))
            k=1
            do i=1,n_sim_total
                if(sum(parametre_estimes_MPI_T(i,:)).ne.0.d0) then
                    !mise à jour des parametres estimes a partir des valeurs issues de la fusion
                    parametre_estimes(k,:)=parametre_estimes_MPI_T(i,:)
                    k=k+1
                endif
            enddo
            moy_theta_est=moy_theta_est_0
            moy_ni=moy_ni_0
            moy_trt=moy_trt_0
            moy_theta=moy_theta_0
            moy_sigmas=moy_sigmas_0
            moy_sigmat=moy_sigmat_0
            moy_rho=moy_rho_0
            moy_gamma=moy_gamma_0
            moy_pros=moy_pros_0
            moy_dec=moy_dec_0
            tab_var_sigma=tab_var_sigma_0
            moy_se_theta=moy_se_theta_0
            taux_couverture_theta=taux_couverture_theta_0
            moy_gamma_est=moy_gamma_est_0
            moy_se_gamma=moy_se_gamma_0
            taux_couverture_gamma=taux_couverture_gamma_0
            moy_sigmas_est=moy_sigmas_est_0
            moy_se_sigmas=moy_se_sigmas_0
            taux_couverture_sigmas=taux_couverture_sigmas_0
            moy_sigmat_est=moy_sigmat_est_0
            moy_se_sigmat=moy_se_sigmat_0
            taux_couverture_sigmat=taux_couverture_sigmat_0
            moy_sigmast_est=moy_sigmast_est_0
            moy_se_sigmast=moy_se_sigmast_0
            taux_couverture_sigmast=taux_couverture_sigmast_0
            moy_eta=moy_eta_0
            moy_se_eta=moy_se_eta_0
            taux_couverture_eta=taux_couverture_eta_0
            moy_betaS=moy_betaS_0
            moy_betaS_se=moy_betaS_se_0
            taux_couvertureS=taux_couvertureS_0
            moy_betaT=moy_betaT_0
            moy_betaT_se=moy_betaT_se_0
            taux_couvertureT=taux_couvertureT_0
            nbre_rejet=nbre_rejet_0
            se_theta_est=se_theta_est_0
            se_eta_est=se_eta_est_0
            se_beta_s=se_beta_s_0
            se_beta_t=se_beta_t_0
            se_sigmas_est=se_sigmas_est_0
            se_sigmat_est=se_sigmat_est_0
            se_cov_est=se_cov_est_0
            se_gamma_est=se_gamma_est_0
            se_alpha_est=se_alpha_est_0
            se_thetat_est=se_thetat_est_0
            se_cov_est_wij=se_cov_est_wij_0
            se_gamma_estt=se_gamma_estt_0
            se_cov_est_ui=se_cov_est_ui_0
            moy_thetat=moy_thetat_0
            se_theta_simt=se_theta_simt_0
            moy_thetat_est=moy_thetat_est_0
            moy_se_thetat=moy_se_thetat_0
            taux_couverture_thetat=taux_couverture_thetat_0
            moy_rho_wij=moy_rho_wij_0
            se_rho_sim_wij=se_rho_sim_wij_0
            moy_thetast_est=moy_thetast_est_0
            se_cov_est=se_cov_est_0
            moy_se_thetast=moy_se_thetast_0
            taux_couverture_thetast=taux_couverture_thetast_0
            moy_gammat=moy_gammat_0
            se_gamma_estt=se_gamma_estt_0
            se_gamma_simt=se_gamma_simt_0
            moy_gammat_est=moy_gammat_est_0
            moy_se_gammat=moy_se_gammat_0
            taux_couverture_gammat=taux_couverture_gammat_0
            moy_rho_ui=moy_rho_ui_0
            se_rho_sim_ui=se_rho_sim_ui_0
            moy_gammast_est=moy_gammast_est_0
            se_cov_est_ui=se_cov_est_ui_0
            moy_se_gammast=moy_se_gammast_0
            taux_couverture_gammast=taux_couverture_gammast_0
            moy_alpha=moy_alpha_0
            se_alpha_est=se_alpha_est_0
            moy_se_alpha=moy_se_alpha_0
            taux_couverture_alpha=taux_couverture_alpha_0
        endif
        deallocate(tampon,tampon_all)
    endif
    !else ! si l'on a qu'un seul processus alors on suppose que l'on fait que du OpenMP et on fait donc ce qui suit:
   if(rang==0)then
        if(nb_processus<=1) n_sim_exact=dble(n_sim-nbre_rejet)
        if((nsim_node(4).ne.0) .and.(nsim_node(4).ne.3)) then
            if(nsim_node(4).eq.1) then
                if(rang_proc==0) then
                    if(logNormal==1)then
                    else
                    endif
                endif
            endif
            if(nsim_node(4).eq.2) then
                if(rang_proc==0) then
                    if(logNormal==1)then
                    else
                        !print*,"cas impossible pour ce type de modele"
                    endif
                endif
            endif
        else
            if(rang_proc==0) then
                if(logNormal==1)then

                else

                endif
            endif
        endif
        if(rang_proc==0) then
            if(nsim_node(8)==2)then
            endif
            if(nsim_node(8).ne.0)then !model conjoint surrogate
                if(frailt_base==1) then
                endif
                if(nsim_node(8)==2)then
                endif
            endif
            if((nsim_node(4).ne.0) .and.(nsim_node(4).ne.3)) then
                if(nsim_node(4).eq.1) then
                    !print*,"nbre_point:",nsim_node(2)
                endif
                if(nsim_node(4).eq.2) then
                endif
            else
            endif
            ni=int(moy_ni/n_sim_exact)
            call cpu_time(tp2)
        endif
        if(gener_only.eq.1) then
            nbre_rejet=dble(0)
            if(nb_processus<=1)then
                n_sim_exact=dble(n_sim-nbre_rejet)
            else
                n_sim_exact=dble(n_sim_total-nbre_rejet)
            endif
            if(rang_proc==0) then
                if(nsim_node(8)==2)then
                endif
                if(nsim_node(8).ne.0)then !model conjoint surrogate
                    if(frailt_base==1) then
                        if(nsim_node(8)==2)then
                        endif
                    endif
                endif
            endif
            goto 1001
        endif
        if(nb_processus<=1) n_sim_exact=dble(n_sim-nbre_rejet) !nombre de simulations qui ont converges
        if(rang_proc==0) then
            !print*,"moyenne empirique des personnes traitées=",    moy_trt/n_sim_exact
            if(une_donnee==0)then !si on ne genere pas les donnees, pas la peine de les afficher
                    if(nsim_node(8)==2)then
                    endif

                    if(nsim_node(8).ne.0)then !model conjoint surrogate

                        if(frailt_base==1) then
                            if(nsim_node(8)==2)then
                            endif
                        endif
                    endif
            endif
            if(nsim_node(8)==2)then
            endif
            if(indice_eta==1)then
            endif
            if(nsim_node(8).ne.0)then !model conjoint surrogate
                if(frailt_base==1) then

                    if(indice_alpha==1)then
                    endif
                endif

                if(nsim_node(8)==2)then
                endif
            endif
        endif
        ! parametres de validation du surrogate
        if(nsim_node(8).ne.0)then
            if(rang_proc==0) then
            endif
            ! calcul des vrais taux de kendall estimes
            if(nsim_node(8).eq.2)then ! modele complet
                theta_ST0(:,1)= (/theta2,thetast_vrai/)
                theta_ST0(:,2)= (/thetast_vrai,theta2_t/)
                gamma_st0(:,1)= (/gamma_ui,gammast_vrai/)
                gamma_st0(:,2)= (/gammast_vrai,gamma_uit/)
                sigma_st0(:,1)= (/sigma_s,sigmast_vrai/)
                sigma_st0(:,2)= (/sigmast_vrai,sigma_t/)
            endif
            if((nsim_node(8).eq.1) .or. (nsim_node(8).eq.3))then ! modele reduit
                allocate(theta_ST0_2(1,1),gamma_st0_2(1,1))
                theta_ST0_2(1,1)= theta2
                gamma_st0_2(1,1)= gamma_ui
                if(indice_eta==0) eta=1.d0
                if(indice_alpha==0) alpha_ui=1.d0
            endif
            if(nsim_node(8).eq.1)then ! modele reduit
                if(method_int_kendal==4 .or. method_int_kendal==5)then ! 1 seul taux de kendall
                    tau_kendal_00=tau_kendall(theta_ST0_2,gamma_st0_2,sigma_st0,0,0,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,0)!tau de kendal des non traites z_11=0,z_21=0
                else
                    tau_kendal_11=tau_kendall(theta_ST0_2,gamma_st0_2,sigma_st0,1,1,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,0)!tau de kendal des traites z_11=1,z_21=1
                    tau_kendal_10=tau_kendall(theta_ST0_2,gamma_st0_2,sigma_st0,1,0,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,0)!tau de kendal des 1 traite et l'autre non traite z_11=1,z_21=0
                    tau_kendal_01=tau_kendall(theta_ST0_2,gamma_st0_2,sigma_st0,0,1,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,0)!tau de kendal des traite z_11=0,z_21=1
                    tau_kendal_00=tau_kendall(theta_ST0_2,gamma_st0_2,sigma_st0,0,0,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,0)!tau de kendal des non traites z_11=0,z_21=0
                endif
            endif
            if(nsim_node(8).eq.3)then ! modele reduit copule
                tau_kendal_00=vrai_tau_copula
            endif
            deallocate(theta_ST0_2,gamma_st0_2)

            if(nsim_node(8).eq.2)then ! modele complet
                if(method_int_kendal==4)then ! 1 seul taux de kendall
                    tau_kendal_00=tau_kendall(theta_ST0,gamma_st0,sigma_st0,0,0,method_int_kendal,&
                    N_MC_kendall,alpha_ui,eta,1)    !tau de kendal des non traites z_11=0,z_21=0
                else
                    tau_kendal_11=tau_kendall(theta_ST0,gamma_st0,sigma_st0,1,1,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,1)    !tau de kendal des traites z_11=1,z_21=1
                    tau_kendal_10=tau_kendall(theta_ST0,gamma_st0,sigma_st0,1,0,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,1)    !tau de kendal des 1 traite et l'autre non traite z_11=1,z_21=0
                    tau_kendal_01=tau_kendall(theta_ST0,gamma_st0,sigma_st0,0,1,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,1)    !tau de kendal des traite z_11=0,z_21=1
                    tau_kendal_00=tau_kendall(theta_ST0,gamma_st0,sigma_st0,0,0,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,1)    !tau de kendal des non traites z_11=0,z_21=0
                endif
            endif
            if(rang_proc==0) then
                if(method_int_kendal==4 .or. method_int_kendal==5)then ! 1 seul taux de kendall et MC
                else
                endif
                ! recherche des taux de couverture par bootsrap
                CP_R2_boot=0.d0
                CP_ktau_boot=0.d0
                remplnsim=INT(n_sim_exact) ! jute pour avoir un INt pour le Do
                do i=1,remplnsim
                    if((result_bootstrap(i,2)<=rsqrt**2) .and. (result_bootstrap(i,3)>=rsqrt**2)) CP_R2_boot=&
                        CP_R2_boot+1
                    if((result_bootstrap(i,5)<=tau_kendal_00) .and. (result_bootstrap(i,6)>=tau_kendal_00)) &
                        CP_ktau_boot=CP_ktau_boot+1
                enddo
                if(rang_proc==0) then
                endif
            endif
        endif
        if(rang_proc==0) then
                if((nsim_node(4).ne.0) .and.(nsim_node(4).ne.3)) then

            else
                if(logNormal==1)then
                else
                endif
            endif
            if(nsim_node(8)==2)then
            endif
            if((nsim_node(4).ne.0) .and.(nsim_node(4).ne.3)) then
                if(nsim_node(4).eq.1) then
                endif
                if(nsim_node(4).eq.2) then
                endif
            else
            endif
            call cpu_time(tp2)
        endif
        if(gener_only.eq.1) then
            nbre_rejet=dble(0)
            !n_sim_exact=dble(n_sim-nbre_rejet)
            if(nb_processus<=1)then
                n_sim_exact=dble(n_sim-nbre_rejet)
            else
                n_sim_exact=dble(n_sim_total-nbre_rejet)
            endif
            if(rang_proc==0) then
                if(nsim_node(8)==2)then
                endif

                if(nsim_node(8).ne.0)then !model conjoint surrogate
                    if(frailt_base==1) then
                        if(nsim_node(8)==2)then
                        endif
                    endif
                endif
            endif
            goto 1001
        endif
        if(nb_processus<=1) n_sim_exact=dble(n_sim-nbre_rejet) !nombre de simulations qui ont converges
        if(rang_proc==0) then
            if(une_donnee==0)then
                if(nsim_node(8)==2)then
                endif

                if(nsim_node(8).ne.0)then !model conjoint surrogate
                    if(frailt_base==1) then
                        if(nsim_node(8)==2)then
                        endif
                    endif
                endif
            endif

            if(nsim_node(8)==2)then
            endif

            if(indice_eta==1)then
            endif

            if(nsim_node(8).ne.0)then !model conjoint surrogate
                if(frailt_base==1) then
                    if(indice_alpha==1)then
                    endif
                endif

                if(nsim_node(8)==2)then
                endif
            endif
            if(nsim_node(8).ne.0)then
                if(method_int_kendal==4 .or. method_int_kendal==5)then ! 1 seul taux de kendall et MC
                else
                endif
                if(rang_proc==0) then
                endif
            endif
        endif
   endif
   1001 continue
   goto 998
   if (istop.ne.1) then
        !write(*,*)"ERREUR : LE MODELE N'A PAS CONVERGE"
   else
        if (effet.eq.1) then
        endif
        moyrecu=moyrecu/ng
        moyrecu=moyrecu/ng
        if (effet.eq.1)then
            if(AG.eq.1)then
            endif
        endif
        if(nva.gt.0)then
            do i=1,nva
                j=(np-nva-nknots-splines_ord+i)*(np-nva+i+1)/2
                bi = b(np-nva-nknots-splines_ord+i) - 1.96*dsqrt(H_hessOut(np-nva-nknots-splines_ord+i,&
                       np-nva-nknots-splines_ord+i))
                bs = b(np-nva-nknots-splines_ord+i) + 1.96*dsqrt(H_hessOut(np-nva-nknots-splines_ord+i,&
                       np-nva-nknots-splines_ord+i))
                wres=(b(np-nva-nknots-splines_ord+i))/dsqrt(H_hessOut(np-nva-nknots-splines_ord+i,&
                      np-nva-nknots-splines_ord+i))
            end do
        endif
   endif
   998 continue
   call date_and_time(dateamj,heure2,zone,values)
   deallocate(don_simultamp,don_simulStamp)
   deallocate(moy_betaS, moy_betaT,moy_betaS_se, moy_betaT_se,taux_couvertureS, taux_couvertureT, &
               theta_chap_copula, v_chap_copula)
   deallocate(d_S,d_T,vbetas,vbetat)
end subroutine jointsurrogate
    !complilation:
    !mpif90 -fopenmp -O3 -o exe_joint_surr_MPI_OMP  Adonnees.f90 Aparameters.f90 autres_fonctions.f90 Integrant_scl.f90 aaOptim_New_scl.f90 aaOptim_New_scl2.f90 funcpa_laplace.f90 aaOptim.f90 aaOptim_SCL_0.f90 aaOptimres.f90 funcpa_adaptative.f90 Integrale_mult_scl.f90 Pour_Adaptative.f90 aaUseFunction.f90 funcpajsplines_surrogate_scl_1.f90 funcpajsplines_surrogate_scl_2.f90 afuncpasres.f90 aresidusMartingale.f90 distance.f90 joint_surrogate.f90 main_Surr_simulation.f90
    !execution
    !time ./exe_joint_surr_MPI_OMP
    !time mpirun -n 1 ./exe_joint_surr_MPI_OMP
