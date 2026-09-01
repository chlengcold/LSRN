# !!!!! Please modify all TODO items !!!!! #

################################
#----- Required Libraries -----#
################################
library(nimbleHMC)
library(dplyr)
library(doParallel)
library(doSNOW)
library(HDInterval)
library(mcmcr)

#########################
#----- LSRN Models -----#
#########################
source('Model.R')

#######################
#----- Utilities -----#
#######################
source('Utils.R')

#####################
#----- Setting -----#
#####################

#!!!!! Important: Data shall be scaled as 1 -- K in a K-point item. !!!!!

Args = list(
  N = 30, # TODO: sample size
  P = 8, # TODO: number of total items
  L = 2, # TODO: number of Layers
  V = 2, # TODO: dimension of \xi 
  K = 5, # TODO: number of categories
  nchain = 3, # TODO: number of chains !!!!!Important: this depends on your temporal memory size!!!!!
  thin = 40, # TODO: thinning number
  niter = 60000, # TODO: number of iterations
  nburnin = 20000, # TODO: number of burn-in
  HPD = .95 # TODO: HPD interval
)

y = array(NA, c(Args[['N']], Args[['N']], Args[['P']]))
for(l in 1:8){
  y[,,l] = read.csv(paste0('toy_dta_',l,'.csv'), row.names = 1, header=T) %>% as.matrix()
}

dta = list(y=y) 

# TODO: items' corresponding dimension
item_dim = rep(1:Args[['L']], each=Args[['P']]/Args[['L']]) 

# TODO: for beta's identifiability: the starting and ending item positions for each dimension, ordered sequentially by dimension
item_seq = c(1, Args[['P']]/Args[['L']], Args[['P']]/Args[['L']]+1, Args[['P']]) 
## e.g., in a two-dimensional case:
## Position 1: the first item of Dimension 1; 
## Position 2: the last item of Dimension 1;
## Position 3: the first item of Dimension 2
## Position 4: the last item of Dimension 2

J_idx = NULL
for(i in 1:Args[['N']]){
  J_idx = rbind(J_idx, (1:Args[['N']])[-i])
}

constants = list(
  V=Args[['V']], 
  K=Args[['K']], 
  N=Args[['N']], 
  P=Args[['P']],
  L=Args[['L']],
  M=rep(0,2*Args[['L']]), 
  M_xi=rep(0,Args[['V']]), Sigma_xi=diag(1,Args[['V']]), 
  J_idx=J_idx, 
  item_dim=item_dim, item_seq=item_seq
) 

inits = list()
for(i in 1:Args[['nchain']]){
  tmp_xi = array(0,c(Args[['N']],Args[['V']],Args[['L']]))
  for(l in 1:Args[['L']]){
    for(v in 1:Args[['V']]){
      tmp_xi[,v,l] = rnorm(Args[['N']], 0, 1) %>% scale()
    }
  }
  inits[[i]] =
    list(
      theta = cbind(scale(apply(y[,,1],1,sum, na.rm=T)), # TODO: sender: dim 1
                    scale(apply(y[,,2],1,sum, na.rm=T)), # TODO: sender: dim 2
                    scale(apply(y[,,1],2,sum, na.rm=T)), # TODO: receiver: dim 1
                    scale(apply(y[,,2],2,sum, na.rm=T))), # TODO: receiver: dim 2
      Ustar = diag(1,2*Args[['L']]),
      
      beta = rep(0,Args[['P']]),
      tau_raw = matrix(seq(-1.5, 1.5, length.out=Args[['K']]-1), Args[['K']]-1, Args[['L']]), 
      
      xi_raw = tmp_xi,
      
      log_lambda = matrix(c(-8, .5), 2, Args[['P']]),
      delta = rep(0, Args[['P']])
    )
}

#######################
#----- Compiling -----#
#######################
my.cluster <- parallel::makeCluster(Args[['nchain']], type = "PSOCK", outfile = 'parallel_log.txt')
doParallel::registerDoParallel(cl = my.cluster)

registerDoSNOW(cl = my.cluster)
pb <- txtProgressBar(max = Args[['nchain']], style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

foreach::getDoParRegistered()
foreach::getDoParWorkers()

chain_output = foreach(chain = 1:Args[['nchain']], .combine = list, .multicombine=TRUE, .options.snow = opts,
                       .packages = c("nimble", "nimbleHMC")
) %dopar% {
  
  current_inits <- inits[[chain]]
  
  tmp = run_MCMC_allcode(data = dta, constants = constants, inits = current_inits,
                         thin=Args[['thin']], niter = Args[['niter']], nburnin = Args[['nburnin']])
  
  out = list(tmp)
  return(out)
}

close(pb)
stopCluster(my.cluster)

gc()

rhat_df_parameter = rhat(as.mcmc(chain_output[[1]][[1]]$samples), 'parameter', as_df = T)
for(chain in 2:Args[['nchain']]){
  rhat_df_parameter = cbind(rhat_df_parameter, rhat2 = rhat(as.mcmc(
    chain_output[[chain]][[1]]$samples), 'parameter', as_df = T)$rhat)
}
rhat_df_parameter = data.frame(rhat_df_parameter$parameter, apply(rhat_df_parameter[,2:(Args[['nchain']]+1)],1,mean))
write.csv(rhat_df_parameter, 'rhat_df_parameter.csv', row.names = T)

########################
#----- Estimating -----#
########################
Est = NULL
for(c in 1:Args[['nchain']]){
  Est = rbind(Est, chain_output[[c]][[1]]$summary[,'Mean'])
}
Est = apply(Est, 2, mean)

n_var_est = names(Est)
est_Tau = Est[grep('tau', n_var_est)] %>% matrix(Args[['K']]-1, Args[['L']])
est_Sigma = Est[grep('Sigma', n_var_est)] %>% matrix(2*Args[['L']],2*Args[['L']])
est_Theta = Est[grep('theta', n_var_est)] %>% matrix(Args[['N']], 2*Args[['L']])
est_Beta = Est[grep('beta', n_var_est)] 

est_Delta = Est[grep('delta', n_var_est)] # \lambda is significant if \delta is greater than 0.5
est_Lambda = Est[grep('lambda', n_var_est)] %>% matrix(2, Args[['P']]) 
# row: 1: insig; 2: sig.
# column: item

#################
#----- HPD -----#
#################
post_samples = NULL
for(c in 1:Args[['nchain']]){
  post_samples = rbind(post_samples, chain_output[[c]][[1]]$samples)
}
n_var = post_samples %>% colnames()

hdi_Theta = post_samples[,grep('theta', n_var)] %>% HDInterval::hdi(.)
hdi_Tau = post_samples[,grep('tau', n_var)] %>% HDInterval::hdi(.)
hdi_Sigma = post_samples[,grep('Sigma', n_var)] %>% HDInterval::hdi(.)
hdi_Beta =post_samples[,grep('beta', n_var)] %>% HDInterval::hdi(.)
hdi_Omega = post_samples[,grep('omega', n_var)] %>% HDInterval::hdi(.)
hdi_Lambda = post_samples[,grep('lambda', n_var)] %>% HDInterval::hdi(.)

#################################
#----- Procrustes matching -----#
#################################
est_Xi = ProcrustesMatching(chain_output, Args[['nchain']]) 

##############################
#----- Model assessment -----#
##############################
# DIC value: LogL: -2logLikelihood; DIC: DIC 
DIC = DIC_LSRN(y)

# Sender Fit & Receiver Fit 
n_var = colnames(post_samples)
Rep_Tau = post_samples[,grep('tau', n_var)] %>% as.matrix() %>% array(dim=c(nrow(post_samples), Args[['K']]-1, Args[['L']]))
Rep_Theta = post_samples[,grep('theta', n_var)] %>% as.matrix() %>% array(dim=c(nrow(post_samples), Args[['N']], 2*Args[['L']]))
Rep_Beta = post_samples[,grep('beta', n_var)] %>% as.matrix()
Rep_Delta = post_samples[,grep('delta', n_var)] %>% as.matrix() %>% array(dim=c(nrow(post_samples), Args[['P']]))
Rep_Lambda = post_samples[,grep('lambda', n_var)] %>% as.matrix() %>% array(dim=c(nrow(post_samples), 2, Args[['P']]))
Rep_Xi = post_samples[,grep('xi', n_var)] %>% as.matrix() %>% array(dim=c(nrow(post_samples), Args[['N']], Args[['V']], Args[['L']]))

s = r = NULL
s_each = r_each = NULL
s_p = r_p = NULL
s_p_each = r_p_each = NULL
for(i in 1:Args[['N']]){
  # preference
  FIT_s = Fit_s(y[i,-i,], est_Theta[i,1:Args[['L']]], est_Theta[-i,(Args[['L']]+1):(2*Args[['L']])], est_Beta, est_Tau, est_Lambda, est_Delta, est_Xi[i,,], est_Xi[-i,,])
  FIT_r = Fit_r(y[-i,i,], est_Theta[-i,1:Args[['L']]], est_Theta[i,(Args[['L']]+1):(2*Args[['L']])], est_Beta, est_Tau, est_Lambda, est_Delta, est_Xi[-i,,], est_Xi[i,,])
  
  tmp_s_p = tmp_r_p = tmp_s_p_each = tmp_r_p_each = 0
  FIT_rep_s_all = FIT_rep_r_all = vector('numeric', nrow(post_samples))
  FIT_rep_s_all_each = FIT_rep_r_all_each = matrix(0, nrow(post_samples), Args[['P']])
  for(reps in 1:nrow(post_samples)){
    
    y_rep_s = data_gen_fit(Rep_Theta[reps,i,1:Args[['L']]], Rep_Theta[reps,-i,(Args[['L']]+1):(2*Args[['L']])], Rep_Beta[reps,], Rep_Tau[reps,,], Rep_Lambda[reps,,], Rep_Delta[reps,], Rep_Xi[reps,i,,], Rep_Xi[reps,-i,,], 's')[1,,]
    FIT_rep_s = Fit_s(y_rep_s, Rep_Theta[reps,i,1:Args[['L']]], Rep_Theta[reps,-i,(Args[['L']]+1):(2*Args[['L']])], Rep_Beta[reps,], Rep_Tau[reps,,], Rep_Lambda[reps,,], Rep_Delta[reps,], Rep_Xi[reps,i,,], Rep_Xi[reps,-i,,], )
    tmp_s_p = tmp_s_p + (FIT_rep_s[[1]] > FIT_s[[1]])
    tmp_s_p_each = tmp_s_p_each + (FIT_rep_s[[2]] > FIT_s[[2]])
    FIT_rep_s_all[reps] = FIT_rep_s[[1]]
    FIT_rep_s_all_each[reps,] = FIT_rep_s[[2]]
    
    y_rep_r = data_gen_fit(Rep_Theta[reps,-i,1:Args[['L']]], Rep_Theta[reps,i,(Args[['L']]+1):(2*Args[['L']])], Rep_Beta[reps,], Rep_Tau[reps,,], Rep_Lambda[reps,,], Rep_Delta[reps,], Rep_Xi[reps,-i,,], Rep_Xi[reps,i,,], 'r')[,1,]
    FIT_rep_r = Fit_r(y_rep_r,Rep_Theta[reps,-i,1:Args[['L']]], Rep_Theta[reps,i,(Args[['L']]+1):(2*Args[['L']])], Rep_Beta[reps,], Rep_Tau[reps,,], Rep_Lambda[reps,,], Rep_Delta[reps,],Rep_Xi[reps,-i,,], Rep_Xi[reps,i,,])
    tmp_r_p = tmp_r_p + (FIT_rep_r[[1]] > FIT_r[[1]])
    tmp_r_p_each = tmp_r_p_each + (FIT_rep_r[[2]] > FIT_r[[2]])
    FIT_rep_r_all[reps] = FIT_rep_r[[1]]
    FIT_rep_r_all_each[reps,] = FIT_rep_r[[2]]
  }
  
  s_p = c(s_p, tmp_s_p / nrow(post_samples))
  r_p = c(r_p, tmp_r_p / nrow(post_samples))
  s_p_each = rbind(s_p_each, tmp_s_p_each / nrow(post_samples))
  r_p_each = rbind(r_p_each, tmp_r_p_each / nrow(post_samples))
  
  s = c(s, FIT_s[[1]])
  r = c(r, FIT_r[[1]])
  s_each = rbind(s_each, FIT_s[[2]])
  r_each = rbind(r_each, FIT_r[[2]])
}

## Sender Fit, p_value
All_result = cbind(SenderFit=s, SenderFit_P=s_p, Receiver_Fit=r, Receiver_Fit_p=r_p)
Each_S_result = matrix(0, Args[['N']], 2*Args[['P']])
colnames(Each_S_result) = paste0(paste0(rep('Dim', Args[['P']]*2), rep(1:Args[['P']],each=2)), rep(c('','_p'), Args[['P']]))
Each_S_result[,2*(1:Args[['P']])-1] = s_each
Each_S_result[,2*(1:Args[['P']])] = s_p_each

## Receiver Fit, p_value
Each_R_result = matrix(0, Args[['N']], 2*Args[['P']])
colnames(Each_R_result) = paste0(paste0(rep('Dim', Args[['P']]*2), rep(1:Args[['P']],each=2)), rep(c('','_p'), Args[['P']]))
Each_R_result[,2*(1:Args[['P']])-1] = r_each
Each_R_result[,2*(1:Args[['P']])] = r_p_each

# Network Fit 
Net = NULL
Net_p = NULL
for(p in 1:Args[['P']]){
  # preference
  FIT_NET = Fit_net(y[,,p],  est_Theta[,1:Args[['L']]], est_Theta[,(Args[['L']]+1):(2*Args[['L']])], est_Beta[p], est_Tau[,item_dim[p]], est_Lambda[,p], est_Delta[p], est_Xi)
  
  tmp_Net_p = 0
  FIT_rep_Net_all = vector('numeric', nrow(post_samples))
  for(reps in 1:nrow(post_samples)){
    
    y_rep_Net = data_gen_fit_p(Rep_Theta[reps,,1:Args[['L']]], Rep_Theta[reps,,(Args[['L']]+1):(2*Args[['L']])], Rep_Beta[reps,p], Rep_Tau[reps,,item_dim[p]], Rep_Lambda[reps,,p], Rep_Delta[reps,p], Rep_Xi[reps,,,])
    FIT_rep_Net = Fit_net(y_rep_Net, Rep_Theta[reps,,1:Args[['L']]], Rep_Theta[reps,,(Args[['L']]+1):(2*Args[['L']])], Rep_Beta[reps,p], Rep_Tau[reps,,item_dim[p]], Rep_Lambda[reps,,p], Rep_Delta[reps,p], Rep_Xi[reps,,,])
    
    tmp_Net_p = tmp_Net_p + (FIT_rep_Net > FIT_NET)
    FIT_rep_Net_all[reps] = FIT_rep_Net
  }
  
  Net_p = c(Net_p, tmp_Net_p / nrow(post_samples))
  Net = c(Net, FIT_NET)
}
## Network Fit, p_value
All_Net_result = cbind(NetworkFit=Net, NetworkFit_P=Net_p)