##### Network Indicies #####
Reci = function(X){
  s_i = apply(X, 1, sum, na.rm=T)
  w_i = matrix(0, nrow(X), nrow(X))
  for(i in 1:nrow(X)){
    for(j in (1:nrow(X))[-i]){
      w_i[i,j] = min(X[i,j],X[j,i],na.rm=T)
    }
  }
  w_i = apply(w_i, 1, sum, na.rm=T)
  S = sum(s_i, na.rm=T)
  W = sum(w_i, na.rm=T)
  r = W/S
  return(r)
}

Clus = function(X,  Max){
  N = ncol(X)
  O = matrix(1, N, N)
  diag(O) = 0
  C_total = matrix(1, Max, N)
  for(i in 0:Max){
    A = 1 * (X >= i)
    AA = (A + t(A))
    AA[is.na(AA)] = 0
    Ci = diag( (AA%*%AA%*%AA) / (2*(AA%*%O%*%AA)) )
    C_total[i,] = Ci
  }
  C_total[is.na(C_total)]=0
  Output = apply(C_total, 2, mean, na.rm=T) %>% mean(., na.rm=T)
  
  return(Output)
}

# Multiplexity (similarity)
multiPlex = function(X,idx){
  N = dim(X)[1]
  sum_wij = 0
  for(i in 1:N){
    for(j in (1:N)[-i]){
      wij = min(X[i,j,idx], na.rm = T)
      wij = ifelse(is.infinite(wij), 0, wij)
      sum_wij = sum_wij + wij
    }
  }
  total_wij = sum(X[,,idx], na.rm = T)
  m = 2*sum_wij/total_wij
  return(m)
} 

s_multiPlex = function(X,idx){
  N = dim(X)[1]
  sum_wi = 0
  for(i in 1:N){
    s = apply(X[i,,idx], 2, sum)
    wi = min(s, na.rm = T)
    wi = ifelse(is.infinite(wi), 0, wi)
    sum_wi = sum_wi + wi
  }
  total_wi = sum(X[,,idx], na.rm = T)
  m = 2*sum_wi/total_wi
  return(m)
} 

r_multiPlex = function(X,idx){
  N = dim(X)[1]
  sum_wj = 0
  for(j in 1:N){
    r = apply(X[,j,idx], 2, sum)
    wj = min(r, na.rm = T)
    wj = ifelse(is.infinite(wj), 0, wj)
    sum_wj = sum_wj + wj
  }
  total_wj = sum(X[,,idx], na.rm = T)
  m = 2*sum_wj/total_wj
  return(m)
} 

# Multireciprocity
multiReci = function(X,idx){
  N = dim(X)[1]
  sum_wij = 0
  for(i in 1:N){
    for(j in (1:N)[-i]){
      wij = min(X[i,j,idx[1]], X[j,i,idx[2]], na.rm = T)
      sum_wij = sum_wij + wij
    }
  }
  total_wij = sum(X[,,idx], na.rm = T)
  m = 2*sum_wij/total_wij
  return(m)
} 

sr_multiReci = function(X,idx){
  N = dim(X)[1]
  sum_wij = 0
  for(i in 1:N){
    s = sum(X[i,,idx[1]])
    r = sum(X[,i,idx[2]])
    wij = min(s, r, na.rm = T)
    sum_wij = sum_wij + wij
  }
  total_wij = sum(X[,,idx], na.rm = T)
  m = 2*sum_wij/total_wij
  return(m)
} 

rs_multiReci = function(X,idx){
  N = dim(X)[1]
  sum_wij = 0
  for(j in 1:N){
    r = sum(X[,j,idx[1]])
    s = sum(X[j,,idx[2]])
    wij = min(s, r, na.rm = T)
    sum_wij = sum_wij + wij
  }
  total_wij = sum(X[,,idx], na.rm = T)
  m = 2*sum_wij/total_wij
  return(m)
} 

##### Procrustes Matching #####
ProcrustesMatching = function(chain_output, nchain=3){
  
  post_samples = NULL
  for(c in 1:nchain){
    post_samples = rbind(post_samples, as.matrix(chain_output[[c]][[1]]$samples))
  }
  n_var = colnames(post_samples)
  sel_Xi= grep('xi', n_var)
  
  len = dim(post_samples)[1]
  
  est_logL = NULL
  for(c in 1:nchain){
    tmp_logL = chain_output[[c]][[1]]$samples[,grep('logL_raw',n_var)] %>% array(dim=c(len/nchain, N,N,K,P))
    est_logL = c(est_logL, rowSums(tmp_logL,dims=1))
  }
  
  ref_point = which.max(est_logL) 
  post_Xi = post_samples[,sel_Xi] %>% array(.,dim=c(len,N,V,L))
  
  Xi_0 = post_Xi[ref_point,,,]

  post_Xi_arr = array(0, dim=c(N, V, L, len))
  
  for(i in 1:len){
    for(l in 1:L){
      proc = MCMCpack::procrustes(post_Xi[i,,,l], Xi_0[,,l], translation = T, dilation = T)
      post_Xi_arr[,,l,i] = proc$X.new # c(N, v, P, len)
    }
  }
  
  est_Xi = rowMeans(post_Xi_arr, dim=3)
  
  return(est_Xi)
}

##### Fit #####
logL = function(y, thetaS, thetaR, beta, tauj, g, w, xi){
  y[is.na(y)] = 0
  logL_outpput = 0
  cats = nrow(tauj)
  N = nrow(thetaS)
  for(p in 1:P){
    for(i in 1: N){
      for(j in (1:N)[-i]){
        # Retrivel from knowledge
        ## Probability of y
        logit = thetaS[i,item_dim[p]] + thetaR[j,item_dim[p]] - beta[p] - 
          g[((w[p]>0.5)+1) ,p] * sqrt(sum((xi[i,,item_dim[p]]- xi[j,,item_dim[p]])^2)) - 
          tauj[,item_dim[p]]
        
        py = exp(cumsum(c(0, logit)))
        py = py / sum(py)
        
        # Reponse of person n
        for(k in 1:(cats+1)){
          log_ijk = (y[i,j,p]==k)*log(py[k])
          logL_outpput = logL_outpput + log_ijk
        }
      }
    }
  }
  
  return(logL_outpput)
}

logL_1L = function(y, thetaS, thetaR, beta, tauj, g, w, xi){
  y[is.na(y)] = 0
  logL_outpput = 0
  cats = length(tauj)
  N = length(thetaS)
  for(p in 1:P){
    for(i in 1: N){
      for(j in (1:N)[-i]){
        # Retrivel from knowledge
        ## Probability of y
        logit = thetaS[i] + thetaR[j] - beta[p] - 
          g[((w[p]>0.5)+1) ,p] * sqrt(sum((xi[i,]- xi[j,])^2)) - 
          tauj
        
        py = exp(cumsum(c(0, logit)))
        py = py / sum(py)
        
        # Reponse of person n
        for(k in 1:(cats+1)){
          log_ijk = (y[i,j,p]==k)*log(py[k])
          logL_outpput = logL_outpput + log_ijk
        }
      }
    }
  }
  
  return(logL_outpput)
}

# DIC
DIC_LSRN = function(y){
  # ppDIC
  DIC = LogL = NULL
  {
    N = nrow(y)
  }
  
  logl = logL(y, est_Theta[,1:L], est_Theta[,(L+1):(2*L)],  est_Beta, est_Tau, est_Lambda, est_Omega, est_Xi)
  
  len = dim(post_samples)[1]
  n_var = colnames(post_samples)
  ppdic = post_samples[,grep('logL', n_var)] %>% as.matrix %>% apply(.,1,sum) %>% sum()
  ppdic = ppdic/len
  
  pdic = 2*(logl-ppdic)
  DIC = c(DIC, -2*logl+2*pdic)
  LogL = c(LogL, -2*logl)
  
  out = list(LogL = LogL, DIC=DIC)
  return(out)
}

# Sender fit
Fit_s = function(y, thetaS, thetaR, beta, tauj, lambda, omega, xiS, xiR, Mod='Dist'){
  
  ProbS = function(thetaS, thetaR, beta, tauj, lambda, omega, xiS, xiR, Mod='Dist'){
    cats = nrow(tauj)
    N = nrow(thetaR)
    Prob = array(0,c(N,P,cats+1))
    for(p in 1:P){
      for(i in 1: N){
        if(Mod=='Dist'){
          Dist_ij = -sqrt(sum((xiS[,item_dim[p]] - xiR[i,,item_dim[p]])^2))
        }
        logit = thetaS[item_dim[p]] + thetaR[i,item_dim[p]] + lambda[(omega[p]>.5)+1,p] * Dist_ij + beta[p] - tauj[,item_dim[p]]

        py = exp(cumsum(c(0, logit)))
        py = py / sum(py)
        
        Prob[i,p,] = py
      }
    }
    
    return(Prob)
  }
  
  K = nrow(tauj)+1
  N = nrow(thetaR) 
  
  Prob_S = ProbS(thetaS, thetaR, beta, tauj, lambda, omega, xiS, xiR, Mod='Dist') 
  m_y = rep(0:(K-1), each=(N)*P) %>% array(dim=c(N,P,K))
  E_S = rowSums(Prob_S * m_y, dims=2)
  m_E_S = array(E_S, dim=c(N,P,K))
  V_S = rowSums(Prob_S * (m_y - m_E_S)^2, dims=2)
  y_S = y-1
  
  S_Fit = sum((y_S - E_S)^2, na.rm=T)/sum(V_S, na.rm=T)
  S_Fit_each = apply((y_S - E_S)^2,2,sum, na.rm=T)/apply(V_S, 2, sum, na.rm=T) 
  
  return(list(S_Fit, S_Fit_each))
}

# Receiver fit
Fit_r = function(y, thetaS, thetaR, beta, tauj, lambda, omega, xiS, xiR, Mod='Dist'){
  
  ProbR = function(thetaS, thetaR, beta, tauj, lambda, omega, xiS, xiR, Mod='Dist'){
    cats = nrow(tauj)
    N = nrow(thetaS)
    Prob = array(0,c(N,P,cats+1))
    for(p in 1:P){
      for(i in 1: N){
        if(Mod=='Dist'){
          Dist_ij = -sqrt(sum((xiS[i,,item_dim[p]] - xiR[,item_dim[p]])^2))
        }
        logit = thetaS[i,item_dim[p]] + thetaR[item_dim[p]] + lambda[(omega[p]>.5)+1,p] * Dist_ij + beta[p] - tauj[,item_dim[p]]

        py = exp(cumsum(c(0, logit)))
        py = py / sum(py)
        
        Prob[i,p,] = py
      }
    }
    
    return(Prob)
  }
  
  K = nrow(tauj)+1
  N = nrow(thetaS) 
  
  Prob_R = ProbR(thetaS, thetaR, beta, tauj, lambda, omega, xiS, xiR, Mod='Dist') 
  m_y = rep(0:(K-1), each=(N)*P) %>% array(dim=c(N,P,K))
  E_R = rowSums(Prob_R * m_y, dims=2)
  m_E_R = array(E_R, dim=c(N,P,K))
  V_R = rowSums(Prob_R * (m_y - m_E_R)^2, dims=2)
  y_R = y-1
  
  R_Fit = sum((y_R - E_R)^2, na.rm=T)/sum(V_R, na.rm=T) 
  R_Fit_each = apply((y_R - E_R)^2,2,sum, na.rm=T)/apply(V_R, 2, sum, na.rm=T) 
  
  return(list(R_Fit, R_Fit_each))
}

data_gen_fit = function(thetaS, thetaR, beta, tauj, lambda, omega, xiS, xiR, SR='s', Mod='Dist'){
  if(SR=='s'){
    cats = nrow(tauj)
    N = nrow(thetaR)
    r = array(0,c(1,N,P))
    for(i in 1:1){
      for(j in 1:N){
        for(p in 1:P){
          if(Mod=='Dist'){
            Dist_ij = -sqrt(sum((xiS[,item_dim[p]] - xiR[j,,item_dim[p]])^2))
          }
          logit = thetaS[item_dim[p]] + thetaR[j,item_dim[p]] + lambda[(omega[p]>.5)+1,p] * Dist_ij + beta[p] - tauj[,item_dim[p]]
   
          py = exp(cumsum(c(0, logit)))
          py = py / sum(py)
          
          x = rcat(1, py)
          
          r[i,j,p] = x
        }
      }
    }
  }else if(SR=='r'){
    cats = nrow(tauj)
    N = nrow(thetaS)
    r = array(0,c(N,1,P))
    for(i in 1:N){
      for(j in 1:1){
        for(p in 1:P){
          if(Mod=='Dist'){
            Dist_ij = -sqrt(sum((xiS[i,,item_dim[p]] - xiR[,item_dim[p]])^2))
          }
          logit = thetaS[i,item_dim[p]] + thetaR[item_dim[p]] + lambda[(omega[p]>.5)+1,p] * Dist_ij +  beta[p] - tauj[,item_dim[p]]

          py = exp(cumsum(c(0, logit)))
          py = py / sum(py)
          
          x = rcat(1, py)
          
          r[i,j,p] = x
        }
      }
    }
  }
  
  return(r)
}

# Network fit
Fit_net = function(y, thetaS, thetaR, beta, tauj, lambda, omega, xi,Mod='Dist'){
  
  Prob = function(thetaS, thetaR, beta, tauj, lambda, omega, xi, Mod='Dist'){
    cats = length(tauj)
    N = nrow(thetaS)
    Prob_Net = array(0,c(N,N,cats+1))
    for(i in 1:N){
      for(j in (1:N)[-i]){
        if(Mod=='Dist'){
          Dist_ij = -sqrt(sum((xi[i,,item_dim[p]] - xi[j,,item_dim[p]])^2))
        }
        logit = thetaS[i,item_dim[p]] + thetaR[j,item_dim[p]] + lambda[(omega>.5)+1] * Dist_ij - beta - tauj

        py = exp(cumsum(c(0, logit)))
        py = py / sum(py)
        
        Prob_Net[i,j,] = py
      }
    }
    
    return(Prob_Net)
  }
  
  K = length(tauj)+1
  N = nrow(thetaS) 
  
  Prob_Net = Prob(thetaS, thetaR, beta, tauj, lambda, omega, xi, Mod='Dist') 
  m_y = rep(0:(K-1), each=N*N) %>% array(dim=c(N,N,K))
  E = rowSums(Prob_Net * m_y, dims=2)
  m_E = array(E, dim=c(N,N,K))
  V = rowSums(Prob_Net * (m_y - m_E)^2, dims=2)
  y_Net = y-1
  
  Net_Fit = sum((y_Net - E)^2, na.rm=T)/sum(V, na.rm=T)
  
  return(Net_Fit)
}

data_gen_fit_p = function(thetaS, thetaR, beta, tauj, lambda, omega, xi, Mod='Dist'){
  cats = length(tauj)
  N = nrow(thetaR)
  r = array(0,c(N,N))
  for(i in 1:N){
    for(j in (1:N)[-i]){
      if(Mod=='Dist'){
        Dist_ij = -sqrt(sum((xi[i,,item_dim[p]] - xi[j,,item_dim[p]])^2))
      }
      logit = thetaS[i,item_dim[p]] + thetaR[j,item_dim[p]] + lambda[(omega>.5)+1] * Dist_ij - beta - tauj
      py = exp(cumsum(c(0, logit)))
      py = py / sum(py)
      
      x = rcat(1, py)
      
      r[i,j] = x
    }
  }
  
  return(r)
}