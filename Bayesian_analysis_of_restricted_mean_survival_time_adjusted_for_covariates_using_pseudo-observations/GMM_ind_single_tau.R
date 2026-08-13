#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#--------- Bayesian Generalized Method of Moments Model with nimble (independence WCM) -------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
########################################################
## Library
########################################################
library(nimble)

########################################################
## Functions
########################################################
GMM_ind_MCMC_fit <- function(data,parM= list(burnin=5000,thin=2,niter=10000),ModelCst=NA,ModelInit=NA){
  ModelCode <- nimbleCode({
    ### Priors  ###
    for (p in 1:np){
      beta[p] ~ dnorm(0, sd = sqrt(1000))
    }
  
    for (i in 1:n) {
      mu[i] <- sum(X[i,1:np] * beta[1:np])
      
      for (p in 1:np) {
        u[i,p] <- (1 / n) * X[i,p] * (Y[i] - mu[i])
      }
    }
    for (p in 1:np) {
      U[p] <- sum(u[1:n,p])
    }
    
    for (i in 1:n) {
      C_temp[i, 1:np, 1:np] <- u[i,1:np] %*% t(u[i,1:np])
    }
    
    for(p in 1:np){
      for(q in 1:np){
        C[p, q] <- sum(C_temp[1:n,p ,q])

      }
    }
    # Calculate Sigma
        Sigma[1:np, 1:np] <- C[1:np, 1:np] - (1 / n) * U[1:np] %*% t(U[1:np])

    
    #### log likelihood
    Sigma_inv[1:np, 1:np] <- inverse(Sigma[1:np, 1:np])
    phi1 <- (t(U[1:np])%*%Sigma_inv[1:np, 1:np]%*%U[1:np])[1,1]  
    phi2 <-  0.5*phi1 + cadj
    zeros1 ~ dpois( phi2 )     # likelihood is exp(-phi[i]) loglik is -phi[i]
  })
    


  ModelConsts <- list(n = data$n,
                      np = data$np,
                      cadj = 1000)
  
  ModelInits <- list(list(beta = c(1, 1)) ,
                     list(beta = c(0.5, 0.5))#,
                     #list(beta = c(2, 1.5))
                     )
  
  ModelData <- list(Y=data$Y,X=data$X,
                    zeros1= 0)
  
  ModelDimensions <-  list(beta = data$np,
                           U = data$np,
                           C_temp = c(data$n, data$np, data$np),
                           C = c(data$np, data$np),
                           Sigma_inv = c(data$np, data$np),
                           Sigma = c(data$np, data$np), 
                           Y = data$n)

  
  mcmc.out <- nimbleMCMC(code = ModelCode, constants = ModelConsts,
                         data = ModelData, inits = ModelInits,
                         nchains = 2, niter = parM$niter,
                         summary = TRUE, WAIC = F,
                         monitors = c("beta") )
  
  mcmc.out$samples$allchain_BT <- data.frame(rbind(cbind(mcmc.out$samples$chain1,chain=rep(1,nrow(mcmc.out$samples$chain1)))[parM$burnin:nrow(mcmc.out$samples$chain1),],
                                                   cbind(mcmc.out$samples$chain2,chain=rep(2,nrow(mcmc.out$samples$chain2)))[parM$burnin:nrow(mcmc.out$samples$chain2),]#,
                                                   #cbind(mcmc.out$samples$chain3,chain=rep(3,nrow(mcmc.out$samples$chain3)))[parM$burnin:nrow(mcmc.out$samples$chain3),]
                                                   )) %>%
    filter(row_number() %% parM$thin == 1) %>% group_by(chain) %>% ggpubr::mutate(sampleid=c(1:n()))
  
  colstopivot  <-   data.frame(name=colnames(mcmc.out$samples$allchain_BT),stringsAsFactors = F) %>% filter(grepl("\\.",name))
  
  # chainwide <- pivot_longer(mcmc.out$samples$allchain_BT,cols =colstopivot$name, names_to = "variable",values_to = "value" ) %>%  
  #   mutate(idpat=as.numeric( gsub("\\.", "", regmatches(x = variable,m = regexpr("\\.[0-9]*",variable, perl=TRUE)) ))) %>%  
  #   mutate(type=  regmatches(x = variable,m = regexpr("[^\\.]*",variable, perl=F))  ) %>% select(.,-"variable")
  #  chainwide <- pivot_wider(chainwide,names_from = "type",values_from = "value")
  
  # mcmc.out$samples$allchainWide <- chainwide
  mcmc.out$samples$colstopivot <-  colstopivot
  mcmc.out$samples$PopulationVariables <-  setdiff(colnames(mcmc.out$samples$allchain_BT), colstopivot$name)
  mcmc.out$data <- data
  mcmc.out$parM <- parM
  mcmc.out$chain <- mcmc.out$samples$allchain_BT
  mcmc.out$prior_parameters <- ModelCst
  return(mcmc.out)
}
