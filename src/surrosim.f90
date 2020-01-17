! fonction permettant de generer les donnees a l'aide du model joint surrogate
subroutine surrosim(don_simul,don_simulS1,n_obs,n_col,lognormal,affiche_stat,vrai_theta,&
            ng,ver,truealpha,propC,cens_A,gamma1,gamma2,theta2,lambda_S,nu_S,lambda_T,nu_T,betas,&
            betat,n_essai,rsqrt,sigma_s,sigma_t,p,prop_i,gamma,alpha,frailt_base,random_generator0,&
            aleatoire, nbre_sim , graine, nbre_don_non_cons,param_weibull0)
    
    ! nbre_don_non_cons : nombre de jeux de donnees a jeter avant de simuler les donnees a retourner 
    ! ceci permet de se rassurer que les donner geneger lors de simulation tiennent compte du seed fixe a la base
    ! tout comme dans la fonction jointsurrogate. avec cette option, l'on peut etre certain que les donnees 
    ! generees dans jointsurrogate pour les simulations sont les memes que celle generees dans la fonction R 
    ! pour la recherche des kappa par validation croisee
     use var_surrogate, only: random_generator, param_weibull
     use Autres_fonctions

      integer, intent(in)::n_essai,frailt_base,affiche_stat,n_obs,n_col,lognormal,ng,ver,&
                            random_generator0, aleatoire,nbre_sim , graine, param_weibull0
      double precision, intent(in)::truealpha,propC,cens_A,gamma1,gamma2,theta2,gamma,alpha,&
                                    lambda_S,nu_S,lambda_T,nu_T,betas,betat,rsqrt,sigma_s,sigma_t
      double precision,dimension(n_essai),intent(in)::prop_i,p      
      double precision,intent(out)::vrai_theta
      double precision,dimension(n_obs,n_col),intent(out)::don_simulS1,don_simul
      integer :: i
      random_generator = random_generator0
      param_weibull = param_weibull0
      
      if(random_generator==1) then
          if(graine >= 0) call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environement de generation
      endif
      
      if(nbre_don_non_cons > 0) then
         do i = 1, nbre_don_non_cons
            call Generation_surrogate(don_simul,don_simulS1,n_obs,n_col,lognormal,affiche_stat,vrai_theta,&
                 ng,ver,truealpha,propC,cens_A,gamma1,gamma2,theta2,lambda_S,nu_S,lambda_T,nu_T,betas,&
                 betat,n_essai,rsqrt,sigma_s,sigma_t,p,prop_i,gamma,alpha,frailt_base)
         enddo
      endif
      
      call Generation_surrogate(don_simul,don_simulS1,n_obs,n_col,lognormal,affiche_stat,vrai_theta,&
            ng,ver,truealpha,propC,cens_A,gamma1,gamma2,theta2,lambda_S,nu_S,lambda_T,nu_T,betas,&
            betat,n_essai,rsqrt,sigma_s,sigma_t,p,prop_i,gamma,alpha,frailt_base)
     
      ! avec ceci on est certain de retourner le bon jeu de donnees 

endsubroutine surrosim 