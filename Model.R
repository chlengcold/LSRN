# Defining model
run_MCMC_allcode = function(data, constants, inits, thin = 40, niter = 60000, nburnin = 20000) {
  library(nimbleHMC)

  code = nimbleCode({ 
    ##### Model #####
    for(i in 1:N) {    
      for(j in J_idx[i,]){
        for(p in 1:P){
          dist[i,j,p] <- sqrt(sum((xi[i,1:V,item_dim[p]] - xi[j,1:V,item_dim[p]])^2))
          eta[i,j,p] <- theta_S[i,item_dim[p]] + theta_R[j,item_dim[p]] - beta[p] - lambda[delta[p]+1,p] * dist[i,j,p]
          
          d[ i, j, 1, p] <- 1	
          for(k in 2:K){
            d[i,j,k,p] <- exp((k-1)*eta[i,j,p] - tau_sum[k, item_dim[p]])
          }
          
          cd[i, j, p] <- sum(d[ i, j, 1:K, p]) 
          
          prob[i, j, 1:K, p] <- d[ i, j, 1:K, p]/cd[i, j, p]
          
          y[i, j, p] ~ dcat(prob[i,j,1:K,p]) 
        }
      } 
    }   
    
    for(l in 1:L){
      tau_sum[1,l] <- 0
      for(k in 2:K){
        tau_sum[k,l] <- sum(tau[1:(k-1),l])
      }
    }
    
    # Priors 
    for(p in 1:P){
      ## Beta
      beta_raw[p] ~ dnorm(0, var=4)
      
      ## lambda
      log_lambda[1, p] ~ dnorm(-8,1)
      log_lambda[2, p] ~ dnorm(0.5,1)
      lambda[1:2, p] <- exp(log_lambda[1:2, p])
      delta[p] ~ dbern(0.5) 
    }
    
    for(l in 1:L){
      beta[(item_seq[2*l-1]):(item_seq[2*l])] <- beta_raw[(item_seq[2*l-1]):(item_seq[2*l])] - mean(beta_raw[(item_seq[2*l-1]):(item_seq[2*l])])
      
      ## tau
      for(k in 1:(K-1)){
        tau_raw[k, l] ~ dnorm(0, var=4)
      }
      tau[1:(K-1), l] <- tau_raw[1:(K-1), l] - mean(tau_raw[1:(K-1), l])    
    }
    
    ## rho
    Ustar[,] ~ dlkj_corr_cholesky(1,2*L) # upper-triangular
    Sigma[1:(2*L),1:(2*L)] <- t(Ustar[1:(2*L),1:(2*L)]) %*% Ustar[1:(2*L),1:(2*L)]
    
    for(i in 1: N){
      theta[i, 1:(2*L)] ~ dmnorm(M[1:(2*L)],  cov = Sigma[1:(2*L),1:(2*L)])
      
      for(l in 1:L){
        xi_raw[i,1:V,l] ~ dmnorm(M_xi[1:V], cov = Sigma_xi[1:V, 1:V])
      }
    } 
    
    theta_S[1:N, 1:L] <- theta[1:N, 1:L]
    theta_R[1:N, 1:L] <- theta[1:N, (L+1):(2*L)]
    
    for(l in 1:L){
      for(v in 1:V){
        xi[1:N, v, l] <- xi_raw[1:N, v, l] - mean(xi_raw[1:N, v, l])
      }
    }
    
  })
  Rmodel <- nimbleModel(code,  data = data, inits = inits, constants = constants, buildDerivs = TRUE)

  monitors = c('theta', 'xi', 'delta', 'lambda', 'beta', 'tau', 'Sigma')
  conf <- configureMCMC(Rmodel, enableWAIC = F, monitors = monitors) 

  addHMC(conf, target = c('theta', 'log_lambda', 'xi_raw', 'Ustar', 'beta_raw', 'tau_raw'), replace=T) 
  Rmcmc <- buildMCMC(conf)
  
  Cmodel <- compileNimble(Rmodel) 
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel) 
  
  samples <- runMCMC(Cmcmc, thin = thin, niter = niter, nburnin = nburnin, nchains = 1,
                     summary = T, WAIC = F, samplesAsCodaMCMC=T)
  
  return(samples)
}
