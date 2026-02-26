# Defining model
run_MCMC_allcode_1L = function(data, constants, inits, thin = 40, niter = 60000, nburnin = 20000, calc_dist=F) {
  library(nimbleHMC)
  
  code = nimbleCode({ 
    ##### Model #####
    for(i in 1:N) {    
      dist_ij[i, i] <- zero_L
      for(j in J_idx[i,]){
        
        # Dist
        dist_ij[i, j] <- sqrt(sum((xi[i, 1:V] - xi[j, 1:V])^2))
        
        for(p in 1:P){
          d[ i, j, 1, p] <- 1		
          for(k in 2:K){
            d[i, j, k, p] <- exp((k-1)*(theta_S[i] + theta_R[j] - beta[p] 
                                        - lambda[delta[p]+1, p] * dist_ij[i,j]) 
                                 - sum(tau[1:(k-1)]) ) 
          }
          
          cd[i, j, p] <- sum(d[ i, j, 1:K, p]) 
          
          prob[i, j, 1:K, p] <- d[ i, j, 1:K, p]/cd[i, j, p]
          
          y[i, j, p] ~ dcat(prob[i,j,1:K,p]) 
        }
      } 
    }   
    
    # Priors 
    for(p in 1:P){
      ## Beta
      beta_raw[p] ~ dnorm(0, var=4)
      
      ## lambda_ind
      log_lambda[1, p] ~ dnorm(-8,1)
      log_lambda[2, p] ~ dnorm(0.5,1)
      lambda[1:2, p] <- exp(log_lambda[1:2, p])
      omega[p] ~ dbeta(1,1)
      delta[p] ~ dbern(omega[p]) 
    }
    # beta[1:P] <- beta_raw[1:P] - mean(beta_raw[1:P])
    
    beta[1:P] <- beta_raw[1:P] - mean(beta_raw[1:P])
    
    ## tau
    for(k in 1:(K-1)){
      tau_raw[k] ~ dnorm(0, var=4)
    }
    tau[1:(K-1)] <- tau_raw[1:(K-1)] - mean(tau_raw[1:(K-1)])    
 
    ## rho
    Ustar[,] ~ dlkj_corr_cholesky(1,2*L) # upper-triangular
    Sigma[1:(2*L),1:(2*L)] <- t(Ustar[1:(2*L),1:(2*L)]) %*% Ustar[1:(2*L),1:(2*L)]
    
    for(i in 1: N){
      theta[i, 1:2] ~ dmnorm(M[1:(2*L)],  cholesky = Ustar[1:(2*L), 1:(2*L)], prec_param = 0)

      logL_raw[i, i, 1:K, 1:P] <- zero_KP[1:K,1:P]
      
      xi_raw[i,1:V] ~ dmnorm(M_xi[1:V], cov = Sigma_xi[1:V, 1:V])
      
      for(p in 1:P){
        for(j in J_idx[i,]){
          for(k in 1:K){
            logL_raw[i, j, k, p] <- log(prob[i, j, k, p]) * (y[i, j, p] == k)
          } # k
        } # j
      } # p
    } # i
    
    theta_S[1:N] <- theta[1:N, 1]
    theta_R[1:N] <- theta[1:N, 2]
    
    for(v in 1:V){
      xi[1:N, v] <- xi_raw[1:N, v] - mean(xi_raw[1:N, v])
    }
    
  })
  Rmodel <- nimbleModel(code,  data = dta, inits = inits, constants = constants, buildDerivs = TRUE)
  
  ## create default MCMC configuration object 
  monitors = c('theta', 'xi', 'omega', 'lambda', 'beta', 'tau', 'dist_ij', 'Sigma', 'logL_raw')
  if(calc_dist==T){
    conf <- configureMCMC(Rmodel, enableWAIC = F,
                          monitors = monitors) 
  }else{
    conf <- configureMCMC(Rmodel, enableWAIC = F,
                          monitors = monitors[-grep('dist', monitors)]) 
  }

  ## add an HMC sampler operating on b0 and b1 
  addHMC(conf, target = c('theta', 'log_lambda', 'omega', 'xi_raw', 'Ustar', 'beta_raw', 'tau_raw'), replace=T) 
  Rmcmc <- buildMCMC(conf)
  
  Cmodel <- compileNimble(Rmodel) 
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel) 
  samples <- runMCMC(Cmcmc, thin = thin, niter = niter, nburnin = nburnin, nchains = 1,
                     summary = T, WAIC = F, samplesAsCodaMCMC=T)
  
  return(samples)
}
