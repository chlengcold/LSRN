setwd('Working/Noncompensatory Multivariate LSRRM LSM/')

{
  library(nimbleHMC)
  library(dplyr)
  library(doParallel)
  
  library(mcmcr)
  library(Cairo)
  
  source('Sim/Gen_RS.R')
  source('Sim/Utils.R')
  source('Sim/Model_rasch_1L.R')
  
  str_class = function(no, item='q_1'){
    options(expressions = 10000)
    
    toADJ = function(dta, k, end, output){
      edge = dta[k,]
      
      i = edge[1] %>% as.character()
      j = edge[2] %>% as.character()
      Link = edge[3] %>% as.numeric()
      
      output[i,j] = Link
      
      k = k + 1 
      
      if(k <= end){
        toADJ(dta, k, end, output)
      }else{
        return(output)
      }
    }
    
    dta = read.csv('https://www.dropbox.com/s/umr5e99drk0umj9/w1link.csv?dl=1', header=T) %>% filter(., classno==no) 
    classmate = c(dta[,'id.x'], dta[,'id.y'])  %>% unique() %>% na.omit()
    classmate = classmate[order(classmate)] 
    
    nclass = length(classmate)
    adjQ1 = matrix(0, nclass, nclass) 
    row.names(adjQ1) = as.character(classmate)
    colnames(adjQ1) = as.character(classmate)
    
    # q1: 我對他
    dtaQ1 = dta[,c('id.x', 'id.y', item)] %>% na.omit()
    colnames(dtaQ1) = c('from', 'to', 'weight')
    # dtaQ1[,'weight'] = dtaQ1[,'weight'] -1
    # dtaQ1 = filter(dtaQ1, weight!=0)
    
    adjQ1 = toADJ(dtaQ1, 1, nrow(dtaQ1), adjQ1) 
    
    check = apply(adjQ1,1,sum)
    posi = which(check==0)
    if(length(posi)>0){adjQ1 = adjQ1[-posi,-posi]}
    check = apply(adjQ1,2,sum)
    posi = which(check==0)
    if(length(posi)>0){adjQ1 = adjQ1[-posi,-posi]}
    
    for(i in 1:nrow(adjQ1)){
      check = which(adjQ1[i,]==0)
      if(length(check)!=0){adjQ1[i,check]=round(mean(adjQ1[i,c(-i, -check)]))}
    }
    
    diag(adjQ1) = 0
    
    output = structure(list(No = no, studID= row.names(adjQ1), Len=nclass, ADJ=adjQ1))
    
    return(output)
  }
}

no = c(23,48)
for(ClassNo in no){
  # ClassNo = 25 #48  
  doc = paste0('Empirical/', ClassNo, '/1L/')
  
  cat('Class',ClassNo,'\n')
  Class_1 = str_class(ClassNo, 'q_1')
  Class_2 = str_class(ClassNo, 'q_3')
  
  N = Class_1$ADJ %>% nrow
  
  calc_dist = F
  P = 2
  L = 1
  V = 2
  K = 3
  item_dim = rep(1:L, each=P/L)
  item_seq = c(1:2)
  
  y = array(0,c(N,N,2))
  y[,,1] = Class_1$ADJ
  y[,,2] = Class_2$ADJ
  dta = list(y=y) 
  ##### Setting #####
  nchain = 3
  J_idx = NULL
  for(i in 1:N){
    J_idx = rbind(J_idx, (1:N)[-i])
  }
  
  constants = list(V=V, K=K, N=N, P=P, L=L, 
                   zero_L=rep(0,L), zero_KP=matrix(0,K,P), 
                   M=rep(0,2*L), M_xi=rep(0,V), Sigma_xi=diag(1,V), 
                   J_idx=J_idx, item_dim=item_dim, item_seq=item_seq, sigma_idx=rep(1:L,2)) 
  
  inits = list()
  for(i in 1:nchain){
    tmp_xi = array(0,c(N,V))
    for(v in 1:V){
      tmp_xi[,v] = rnorm(N, 0, 1) %>% scale()
    }
    inits[[i]] =
      list(
        theta = cbind(scale(apply(y[,,1:2],1,sum, na.rm=T)),scale(apply(y[,,1:2],2,sum, na.rm=T))),
        Ustar = diag(1,2*L),
        
        beta = rep(0,P),
        tau_raw = (seq(-1.5,1.5,length.out=K-1)), 
        
        xi_raw = tmp_xi,
        
        log_lambda = matrix(c(-5, .5), 2, P),
        omega = rep(.1, P),
        delta = rep(0, P)
      )
  }
  
  # Parallel backend
  my.cluster <- parallel::makeCluster(nchain, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()
  foreach::getDoParWorkers()
  
  chain_output = foreach(c = 1:nchain, .combine = list, .multicombine=TRUE) %dopar% {
    
    tmp = run_MCMC_allcode_1L(data = dta, constants = constants, inits = inits, 
                           thin=10, niter = 15000, nburnin=5000, calc_dist = F)  
    
    out = list(tmp)
    return(out)
  }
  
  stopCluster(my.cluster)
  
  rhat_df = rhat(as.mcmc(
    chain_output[[1]][[1]]$samples), 'parameter', as_df = T)
  write.csv(rhat_df, paste0(doc,'/rhat_df.csv'), row.names = T)
  
  ##### Estimating #####
  Est = NULL
  for(c in 1:nchain){
    Est = rbind(Est, chain_output[[c]][[1]]$summary[,'Mean'])
  }
  Est = apply(Est, 2, mean)
  
  n_var_est = names(Est)
  est_Tau = Est[grep('tau', n_var_est)] %>% matrix(K-1,L)
  est_Sigma = Est[grep('Sigma', n_var_est)] %>% matrix(2*L,2*L)
  est_Theta = Est[grep('theta', n_var_est)] %>% matrix(N, 2*L)
  est_Beta = Est[grep('beta', n_var_est)] 
  
  est_Omega = Est[grep('omega', n_var_est)]
  est_Lambda = Est[grep('lambda', n_var_est)] %>% matrix(2, P)
  
  if(calc_dist){
    est_Dist = Est[grep('dist', n_var_est)] %>% array(c(N, N, L))
  }
  
  ##### Procrustes matching #####
  est_Xi = ProcrustesMatching_1L(chain_output, nchain) 
  
  ##### Standard Deviation #####
  post_samples = NULL
  for(c in 1:nchain){
    post_samples = rbind(post_samples, chain_output[[c]][[1]]$samples)
  }
  n_var = post_samples %>% colnames()
  
  sd_Theta = apply(post_samples[,grep('theta', n_var)],2,sd) %>% matrix(N,2*L)
  sd_Tau = apply(post_samples[,grep('tau', n_var)],2,sd) %>% matrix(K-1, L)
  sd_Sigma = apply(post_samples[,grep('Sigma', n_var)],2,sd) %>% matrix(2*L,2*L)
  sd_Beta =apply(post_samples[,grep('beta', n_var)],2,sd) 
  sd_Omega = apply(post_samples[,grep('omega', n_var)],2,sd) 
  sd_Lambda = apply(post_samples[,grep('lambda', n_var)],2,sd) %>% matrix(2, P)
  
  ##### Saving #####
  write.csv(post_samples, paste0(doc,'/post_samples.csv'), row.names = T)
  
  write.csv(cbind(est_Theta, sd_Theta), paste0(doc,'/Est_Theta.csv'), row.names = T)
  write.csv(cbind(est_Beta, sd_Beta), paste0(doc,'/Est_Beta.csv'), row.names = T)
  write.csv(cbind(as.vector(est_Tau), as.vector(sd_Tau)), paste0(doc,'/Est_Tau.csv'), row.names = T)
  write.csv(cbind((est_Lambda), (sd_Lambda)), paste0(doc,'/Est_Lambda.csv'), row.names = T)
  write.csv(cbind(est_Omega, sd_Omega), paste0(doc,'/Est_Omega.csv'), row.names = T)
  write.csv(est_Xi, paste0(doc,'/Est_Xi.csv'), row.names = T)
  write.csv(cbind(est_Sigma, sd_Sigma), paste0(doc,'/Est_Sigma.csv'), row.names = T)
  
  ##### Interaction plot #####
  CairoPDF(file=paste0(doc,'/Interaction_plot'), onefile=T, width = 8, height=8) 
  minX = c(est_Xi) %>% min() 
  maxX = c(est_Xi, -minX) %>% max()
  plot(est_Xi, col=c("black"), pch=16, xlab='', ylab='', cex=2, cex.lab=3, cex.axis=1.5, xlim=c(-maxX, maxX), ylim=c(-maxX, maxX) )
  dev.off()
}

{
  BIG515 = function(x){
    o = sum(x[13:15])
    c = sum(x[7:9])
    e = sum(x[1:3])
    a = sum(x[4:6])
    n = sum(x[10:12])
    return(c(o,c,e,a,n))
  }
  
  ClassNo = 25
  doc = paste0('Empirical/', ClassNo,'/')
  est_Theta = read.csv(paste0(doc,'Est_Theta.csv'), row.names = 1)[,1:4] %>% as.matrix 
  N = nrow(est_Theta)
  est_Xi = read.csv(paste0(doc,'Est_Xi.csv'), row.names = 1)[,1:4] %>% as.matrix  %>%  array(., dim=c(N, 2, 2)) 
  
  w1_codebook <- read.csv("Empirical/w1_codebook.csv", header=T) %>% filter(school==ClassNo)
  w2_codebook <- read.csv("Empirical/w2_codebook.csv", header=T) %>% filter(school==ClassNo)
  w3_codebook <- read.csv("Empirical/w3_codebook.csv", header=T) %>% filter(school==ClassNo)
  
  Gender = w1_codebook$gender 
  Big5 = w3_codebook %>% select(paste0('ifeel_', 1:15)) %>% apply(.,1,BIG515) %>% t() %>% scale()
  Openness = Big5[,1]
  Conscientiousness = Big5[,2]
  Extraversion = Big5[,3]
  Agreeableness = Big5[,4]
  Neuroticism = Big5[,5]
  Learningmotivation = w1_codebook %>% select(paste0('studystrat_', 1:5)) %>% apply(.,1,sum) %>% scale()
  
  X = cbind(Gender, Learningmotivation, Openness, Conscientiousness, Extraversion, Agreeableness, Neuroticism)
  
  mod = lm(est_Theta[,c(1,3)]~Gender+Learningmotivation+Openness+Conscientiousness+Extraversion+Agreeableness+Neuroticism)
  mod = lm(est_Theta[,c(1,3)]~Gender+Openness+Conscientiousness+Extraversion+Agreeableness+Neuroticism)
  result = summary(mod)
  write.csv(result$`Response V1`$coefficients, 'empirical/result/Valid/mod_1s.csv')
  write.csv(result$`Response V3`$coefficients, 'empirical/result/Valid/mod_1r.csv')
  
  mod2 = lm(est_Theta[,c(2,4)]~Gender+Learningmotivation+Openness+Conscientiousness+Extraversion+Agreeableness+Neuroticism)
  mod2 = lm(est_Theta[,c(2,4)]~Gender+Openness+Conscientiousness+Extraversion+Agreeableness+Neuroticism)
  result2 = summary(mod2)
  write.csv(result2$`Response V2`$coefficients, 'empirical/result/Valid/mod_2s.csv')
  write.csv(result2$`Response V4`$coefficients, 'empirical/result/Valid/mod_2r.csv')
  
  ### gender x variable
  t.test(Openness~Gender)
  Big5 = w3_codebook %>% select(paste0('ifeel_', 1:15)) %>% apply(.,1,BIG515) %>% t()
  Openness = Big5[,1]
  Conscientiousness = Big5[,2]
  Extraversion = Big5[,3]
  Agreeableness = Big5[,4]
  Neuroticism = Big5[,5]
  
  t.test(Openness~Gender)
  t.test(Conscientiousness~Gender)
  t.test(Extraversion~Gender)
  t.test(Agreeableness~Gender)
  t.test(Neuroticism~Gender)
  
  ### Gender ~ xi
  mod = FNN::knn.reg(train=est_Xi[,,1], y=Gender, k=6)
  plot(Gender, mod$pred)
  summary(mod)
  
  ### Gender ~ node
  library(igraph)
  
  ClassNo = 25 #48  
  doc = paste0('Empirical/', ClassNo, '/1L/')
  
  cat('Class',ClassNo,'\n')
  Class_1 = str_class(ClassNo, 'q_1')
  Class_2 = str_class(ClassNo, 'q_3')
  
  N = Class_1$ADJ %>% nrow
  
  y = array(0,c(N,N,2))
  y[,,1] = Class_1$ADJ
  y[,,2] = Class_2$ADJ
  
  #### Familiarity
  p = 1
  g_y = graph_from_adjacency_matrix(y[,,p])
  degree_out = degree(g_y, mode = 'out') 
  degree_in = degree(g_y, mode = 'in')
  trans = transitivity(y[,,p], 3)
  mod = glm(Gender~degree_in + degree_out + trans, 'binomial')
  summary(mod)
  
  t.test(degree_out ~ Gender)
  t.test(degree_in ~ Gender)
  t.test(trans ~ Gender)
  
  #### Liking
  p = 2
  g_y = graph_from_adjacency_matrix(y[,,p])
  degree_out = degree(g_y, mode = 'out') 
  degree_in = degree(g_y, mode = 'in')
  trans = transitivity(y[,,p], 3)
  mod = glm(Gender~degree_in + degree_out + trans, 'binomial')
  summary(mod)
  
  t.test(degree_out ~ Gender)
  t.test(degree_in ~ Gender)
  t.test(trans ~ Gender)
}
