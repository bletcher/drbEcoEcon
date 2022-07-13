tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis"))



# Priors for psi where more likely to stay than move

# states
        
# River      size1  size2   size3
# Main        1       2       3
# Trib        4       5       6
same <- 0.75
shrink <- 0.01
grow1 <- 0.025
grow2 <- 0.001
move <- 0.025

alphaR <- list()
alphaR[[1]] <- alphaR1 <- c(same,   grow1,  grow2, move,   move,    move) 
alphaR[[2]] <- alphaR2 <- c(shrink, same,   grow1, shrink, move,    move)
alphaR[[3]] <- alphaR3 <- c(shrink, shrink, same,  shrink, shrink,  move)
alphaR[[4]] <- alphaR4 <- c(move,   move,   move,  same,   grow1,   grow2)
alphaR[[5]] <- alphaR5 <- c(shrink, move,   move,  shrink, same,    grow1)
alphaR[[6]] <- alphaR6 <- c(shrink, shrink, move,  shrink, shrink,  same)

getDirchPriorsR <- function(nStates, myConstants, alphaR){
  a = array(rep(0, nStates * nStates * (myConstants$T - 1)) , c(nStates, nStates, (myConstants$T - 1) ))
  for(r in 1:nStates){
    for(t in 1:(myConstants$T - 1)){
        dirch <- rdirch(1, alphaR[[r]])
        for (r2 in 1:nStates){
          a[r,r2,t] <- dirch[r2]
        }
      }
    
  }
  return(a)
}

# target_ttt   file name will be same
   # start all with ttt
target_CMR_models_phiT_pT_psiT = 
  tar_plan(
    
    nStates = length(unique(target_eh$data$state)),
    nRivers = nStates, # for now
    first = target_eh$first, #apply(y, 1, function(x) min(which(x !=0)))
    last = target_eh$last,
    
    runData_ttt = list(
      # Updateable model-specific variables 
      nIter = 2500, 
      nBurnin = 1000, 
      nChains = 2,
      thinRate = 5
    ),
    
    inputData_ttt = list(
      y = target_stateMatrix,
      first = target_first, 
      last = target_last, 
      zInits = target_zInits 
    ),
    
    # Proportion of fish in each river on the first observation
    #y1 = tar_read(inputData_ttt)$y,
    y1 = inputData_ttt$y,
    deltaProps = table(y1[y1>0]) / length(y1[y1>0]),
    
    myConstants_ttt = list(
      N = nrow(inputData_ttt$y),
      T = ncol(inputData_ttt$y),
      first = inputData_ttt$first,
      last = inputData_ttt$last,
      
      nRivers = nRivers,
      length = last - first + 1,
      alphaR1 = alphaR1,
      alphaR2 = alphaR2,
      alphaR3 = alphaR3,
      alphaR4 = alphaR4,
      alphaR5 = alphaR5,
      alphaR6 = alphaR6,
      
      deltaProps = deltaProps,
      nStates = nStates
    ),
    
    myData_ttt = list(
      y = inputData_ttt$y + 1
    ),
    
    initialValues_ttt =  
      list(
        betaPhiRiver = array(runif(nRivers, 0, 1), c(nRivers)),
        betaPhi = array(rnorm(nRivers * (myConstants_ttt$T - 1), 0, 1), c(nRivers, (myConstants_ttt$T - 1))),
        
        betaPRiver = array(runif(nRivers, 0, 1), c(nRivers)),
        betaP = array(rnorm(nRivers * (myConstants_ttt$T - 1), 0, 1), c(nRivers, (myConstants_ttt$T - 1))),
  
        psi = getDirchPriorsR(nRivers, myConstants_ttt, alphaR)
      ),
    
    parametersToSave_ttt = c("betaPhi", "betaPhiRiver", 
                             "betaP",   "betaPRiver", 
                             "betaPhiOut", "betaPhiRiverOut",
                             "betaPOut",   "betaPRiverOut", 
                             
                             "psi"),
    
    modelCode_ttt = nimbleCode({
      # Initial distribution among rivers
      delta[1] <- deltaProps[1]                  # Pr(alive t = 1 and in river 1) = 0.4
      delta[2] <- deltaProps[2]
      delta[3] <- deltaProps[3]
      delta[4] <- deltaProps[4]
      delta[5] <- deltaProps[5]
      delta[6] <- deltaProps[6]
      delta[7] <- 0                    # Pr(dead t = 1) = 0

      for (r in 1:nStates){
        betaPhiRiver[r] ~ dnorm(0,1)
        betaPRiver[r] ~ dnorm(0,1)

        betaPhiRiverOut[r] <- ilogit(betaPhiRiver[r])
        betaPRiverOut[r] <- ilogit(betaPRiver[r])

        for (t in 1:(T-1)){
          betaPhi[r,t] ~ dnorm(betaPhiRiver[r],1)
          betaP[r,t] ~ dnorm(betaPRiver[r],1)

          betaPhiOut[r,t] <- ilogit(betaPhi[r,t])
          betaPOut[r,t] <- ilogit(betaP[r,t])

          # move from river 'r' to one of river 1:nStates
          # Nice description of effect of 'alpha' on probabilities:
          # https://stats.stackexchange.com/questions/244917/what-exactly-is-the-alpha-in-the-dirichlet-distribution
          #psi[r,1:nStates,t] ~ ddirch(alpha[1:nStates])
        }
      }
    
      for (t in 1:(T-1)){ # loop over time
          psi[1,1:nRivers,t] ~ ddirch(alphaR1[1:nRivers])
          psi[2,1:nRivers,t] ~ ddirch(alphaR2[1:nRivers])
          psi[3,1:nRivers,t] ~ ddirch(alphaR3[1:nRivers])
          psi[4,1:nRivers,t] ~ ddirch(alphaR4[1:nRivers])
  
      }

      for (t in 1:(T-1)){ # loop over time
        gamma[1,1,t] <- ilogit(betaPhi[1,t]) * psi[1,1,t]
        gamma[1,2,t] <- ilogit(betaPhi[1,t]) * psi[1,2,t]
        gamma[1,3,t] <- ilogit(betaPhi[1,t]) * psi[1,3,t]
        gamma[1,4,t] <- ilogit(betaPhi[1,t]) * psi[1,4,t]
        gamma[1,5,t] <- ilogit(betaPhi[1,t]) * psi[1,5,t]
        gamma[1,6,t] <- ilogit(betaPhi[1,t]) * psi[1,6,t]          
        gamma[1,7,t] <- 1 - ilogit(betaPhi[1,t])
        gamma[2,1,t] <- ilogit(betaPhi[2,t]) * psi[2,1,t]
        gamma[2,2,t] <- ilogit(betaPhi[2,t]) * psi[2,2,t]
        gamma[2,3,t] <- ilogit(betaPhi[2,t]) * psi[2,3,t]
        gamma[2,4,t] <- ilogit(betaPhi[2,t]) * psi[2,4,t]
        gamma[2,5,t] <- ilogit(betaPhi[2,t]) * psi[2,5,t]
        gamma[2,6,t] <- ilogit(betaPhi[2,t]) * psi[2,6,t]
        gamma[2,7,t] <- 1 - ilogit(betaPhi[2,t])
        gamma[3,1,t] <- ilogit(betaPhi[3,t]) * psi[3,1,t]
        gamma[3,2,t] <- ilogit(betaPhi[3,t]) * psi[3,2,t]
        gamma[3,3,t] <- ilogit(betaPhi[3,t]) * psi[3,3,t]
        gamma[3,4,t] <- ilogit(betaPhi[3,t]) * psi[3,4,t]
        gamma[3,5,t] <- ilogit(betaPhi[3,t]) * psi[3,5,t]
        gamma[3,6,t] <- ilogit(betaPhi[3,t]) * psi[3,6,t]
        gamma[3,7,t] <- 1 - ilogit(betaPhi[3,t])
        gamma[4,1,t] <- ilogit(betaPhi[4,t]) * psi[4,1,t]
        gamma[4,2,t] <- ilogit(betaPhi[4,t]) * psi[4,2,t]
        gamma[4,3,t] <- ilogit(betaPhi[4,t]) * psi[4,3,t]
        gamma[4,4,t] <- ilogit(betaPhi[4,t]) * psi[4,4,t]
        gamma[4,5,t] <- ilogit(betaPhi[4,t]) * psi[4,5,t]
        gamma[4,6,t] <- ilogit(betaPhi[4,t]) * psi[4,6,t]
        gamma[4,7,t] <- 1 - ilogit(betaPhi[4,t])
        
        gamma[5,1,t] <- ilogit(betaPhi[5,t]) * psi[5,1,t]
        gamma[5,2,t] <- ilogit(betaPhi[5,t]) * psi[5,2,t]
        gamma[5,3,t] <- ilogit(betaPhi[5,t]) * psi[5,3,t]
        gamma[5,4,t] <- ilogit(betaPhi[5,t]) * psi[5,4,t]
        gamma[5,5,t] <- ilogit(betaPhi[5,t]) * psi[5,5,t]
        gamma[5,6,t] <- ilogit(betaPhi[5,t]) * psi[5,6,t]
        gamma[5,7,t] <- 1 - ilogit(betaPhi[5,t])
        gamma[6,1,t] <- ilogit(betaPhi[6,t]) * psi[6,1,t]
        gamma[6,2,t] <- ilogit(betaPhi[6,t]) * psi[6,2,t]
        gamma[6,3,t] <- ilogit(betaPhi[6,t]) * psi[6,3,t]
        gamma[6,4,t] <- ilogit(betaPhi[6,t]) * psi[6,4,t]
        gamma[6,5,t] <- ilogit(betaPhi[6,t]) * psi[6,5,t]
        gamma[6,6,t] <- ilogit(betaPhi[6,t]) * psi[6,6,t]
        gamma[6,7,t] <- 1 - ilogit(betaPhi[6,t])
        
        gamma[7,1,t] <- 0
        gamma[7,2,t] <- 0
        gamma[7,3,t] <- 0
        gamma[7,4,t] <- 0
        gamma[7,5,t] <- 0
        gamma[7,6,t] <- 0
        gamma[7,7,t] <- 1
        
      }

      # gamma for the last occasion
      for (a in 1:(nStates+1)){
        for (b in 1:nStates){
          gamma[a,b,T] <- 0
        }
        gamma[a,7,T] <- 1
      }

      for (i in 1:N){ # loop over individuals
        # omega for first obs
        omega[1,1,first[i],i] <- 0          # Pr(alive A t -> non-detected t)
        omega[1,2,first[i],i] <- 1          # Pr(alive A t -> detected A t)
        omega[1,3,first[i],i] <- 0          # Pr(alive A t -> detected B t)
        omega[1,4,first[i],i] <- 0          # Pr(alive A t -> detected C t)
        omega[1,5,first[i],i] <- 0          # Pr(alive A t -> detected D t)
        omega[1,6,first[i],i] <- 0
        omega[1,7,first[i],i] <- 0
        omega[2,1,first[i],i] <- 0          # Pr(alive B t -> non-detected t)
        omega[2,2,first[i],i] <- 0          # Pr(alive B t -> detected A t)
        omega[2,3,first[i],i] <- 1          # Pr(alive B t -> detected B t)
        omega[2,4,first[i],i] <- 0          # Pr(alive B t -> detected C t)
        omega[2,5,first[i],i] <- 0          # Pr(alive B t -> detected C t)
        omega[2,6,first[i],i] <- 0
        omega[2,7,first[i],i] <- 0
        omega[3,1,first[i],i] <- 0          # Pr(alive C t -> non-detected t)
        omega[3,2,first[i],i] <- 0          # Pr(alive C t -> detected A t)
        omega[3,3,first[i],i] <- 0          # Pr(alive C t -> detected B t)
        omega[3,4,first[i],i] <- 1          # Pr(alive C t -> detected C t)
        omega[3,5,first[i],i] <- 0          # Pr(alive C t -> detected C t)
        omega[3,6,first[i],i] <- 0
        omega[3,7,first[i],i] <- 0
        omega[4,1,first[i],i] <- 0          # Pr(dead t -> non-detected t)
        omega[4,2,first[i],i] <- 0          # Pr(dead t -> detected A t)
        omega[4,3,first[i],i] <- 0          # Pr(dead t -> detected B t)
        omega[4,4,first[i],i] <- 0          # Pr(dead t -> detected C t)
        omega[4,5,first[i],i] <- 1          # Pr(dead t -> detected C t)
        omega[4,6,first[i],i] <- 0
        omega[4,7,first[i],i] <- 0
        omega[5,1,first[i],i] <- 0          # Pr(dead t -> non-detected t)
        omega[5,2,first[i],i] <- 0          # Pr(dead t -> detected A t)
        omega[5,3,first[i],i] <- 0          # Pr(dead t -> detected B t)
        omega[5,4,first[i],i] <- 0          # Pr(dead t -> detected C t)
        omega[5,5,first[i],i] <- 0          # Pr(dead t -> detected D t)
        omega[5,6,first[i],i] <- 1
        omega[5,7,first[i],i] <- 0
        
        omega[6,1,first[i],i] <- 0          
        omega[6,2,first[i],i] <- 0          
        omega[6,3,first[i],i] <- 0      
        omega[6,4,first[i],i] <- 0        
        omega[6,5,first[i],i] <- 0        
        omega[6,6,first[i],i] <- 0
        omega[6,7,first[i],i] <- 1
        omega[7,1,first[i],i] <- 1          
        omega[7,2,first[i],i] <- 0      
        omega[7,3,first[i],i] <- 0        
        omega[7,4,first[i],i] <- 0          
        omega[7,5,first[i],i] <- 0          
        omega[7,6,first[i],i] <- 0
        omega[7,7,first[i],i] <- 0


        ## DT changes:
        ## time t > first[i]:
        for(t in (first[i]+1):last[i]) {
          logit(pA[t,i]) <- betaP[1,t-1]
          logit(pB[t,i]) <- betaP[2,t-1]
          logit(pC[t,i]) <- betaP[3,t-1]
          logit(pD[t,i]) <- betaP[4,t-1]
          
          logit(pE[t,i]) <- betaP[5,t-1]
          logit(pF[t,i]) <- betaP[6,t-1]

          # probabilities of y(t) given z(t)
          # omega[z, y, t, i]

          # z=1 = alive in River 1, z=2 = alive in River 2...z=5 = dead
          # y=1 = unobserved, y=2 = observed in River 1, y=3 = observed in River 2, etc

          omega[1,1,t,i] <- 1 - pA[t,i]     # Pr(alive A t -> non-detected t)
          omega[1,2,t,i] <- pA[t,i]         # Pr(alive A t -> detected A t)
          omega[1,3,t,i] <- 0               # Pr(alive A t -> detected B t)
          omega[1,4,t,i] <- 0               # Pr(alive A t -> detected C t)
          omega[1,5,t,i] <- 0               # Pr(alive A t -> detected D t)
          omega[1,6,t,i] <- 0
          omega[1,7,t,i] <- 0
          omega[2,1,t,i] <- 1 - pB[t,i]     # Pr(alive B t -> non-detected t)
          omega[2,2,t,i] <- 0               # Pr(alive B t -> detected A t)
          omega[2,3,t,i] <- pB[t,i]         # Pr(alive B t -> detected B t)
          omega[2,4,t,i] <- 0               # Pr(alive B t -> detected C t)
          omega[2,5,t,i] <- 0               # Pr(alive B t -> detected C t)
          omega[2,6,t,i] <- 0
          omega[2,7,t,i] <- 0
          omega[3,1,t,i] <- 1 - pC[t,i]     # Pr(alive C t -> non-detected t)
          omega[3,2,t,i] <- 0               # Pr(alive C t -> detected A t)
          omega[3,3,t,i] <- 0               # Pr(alive C t -> detected B t)
          omega[3,4,t,i] <- pC[t,i]         # Pr(alive C t -> detected C t)
          omega[3,5,t,i] <- 0               # Pr(alive C t -> detected C t)
          omega[3,6,t,i] <- 0
          omega[3,7,t,i] <- 0
          omega[4,1,t,i] <- 1 - pD[t,i]     # Pr(alive D t -> non-detected t))
          omega[4,2,t,i] <- 0               # Pr(dead D t -> detected A t)
          omega[4,3,t,i] <- 0               # Pr(dead D t -> detected B t)
          omega[4,4,t,i] <- 0               # Pr(dead D t -> detected C t)
          omega[4,5,t,i] <- pD[t,i]         # Pr(alive D t -> detected D t)
          omega[4,6,t,i] <- 0
          omega[4,7,t,i] <- 0
          
          omega[5,1,t,i] <- 1 - pE[t,i]               # Pr(dead t -> non-detected t)
          omega[5,2,t,i] <- 0               # Pr(dead t -> detected A t)
          omega[5,3,t,i] <- 0               # Pr(dead t -> detected B t)
          omega[5,4,t,i] <- 0               # Pr(dead t -> detected C t)
          omega[5,5,t,i] <- 0               # Pr(dead t -> detected D t)
          omega[5,6,t,i] <- pE[t,i]
          omega[5,7,t,i] <- 0
          
          omega[6,1,t,i] <- 1 - pF[t,i]               # Pr(dead t -> non-detected t)
          omega[6,2,t,i] <- 0               # Pr(dead t -> detected A t)
          omega[6,3,t,i] <- 0               # Pr(dead t -> detected B t)
          omega[6,4,t,i] <- 0               # Pr(dead t -> detected C t)
          omega[6,5,t,i] <- 0               # Pr(dead t -> detected D t)
          omega[6,6,t,i] <- 0
          omega[6,7,t,i] <- pF[t,i]
          omega[7,1,t,i] <- 1               # Pr(dead t -> non-detected t)
          omega[7,2,t,i] <- 0               # Pr(dead t -> detected A t)
          omega[7,3,t,i] <- 0               # Pr(dead t -> detected B t)
          omega[7,4,t,i] <- 0               # Pr(dead t -> detected C t)
          omega[7,5,t,i] <- 0               # Pr(dead t -> detected D t)
          omega[7,6,t,i] <- 0
          omega[7,7,t,i] <- 0
        }

      } # i loop

      for (i in 1:N){
        y[i,first[i]:last[i]] ~ dDHMMo(init = delta[1:7],
                                       probTrans = gamma[1:7, 1:7, first[i]:last[i]],
                                       probObs =   omega[1:7, 1:7, first[i]:last[i], i],
                                       len = length[i],
                                       checkRowSums = 1)
      }

    }),

    Rmodel_ttt = nimbleModel(
      code = modelCode_ttt,
      constants = myConstants_ttt,
      data = myData_ttt,
      inits = initialValues_ttt,
      calculate = FALSE
    ),
    conf_ttt = configureMCMC(
      Rmodel_ttt,
      monitors = parametersToSave_ttt
    ),

    Rmcmc_ttt = buildMCMC(conf_ttt, useConjugacy = FALSE),
    Cmodel_ttt = compileNimble(Rmodel_ttt),
    Cmcmc_ttt = compileNimble(Rmcmc_ttt, project = Rmodel_ttt),

    model_ttt = runMCMC(
      Cmcmc_ttt,
      niter = runData_ttt$nIter,
      nburnin = runData_ttt$nBurnin,
      thin = runData_ttt$thinRate,
      nchains = runData_ttt$nChains
    ),

    modelOut_ttt =
      list(
        mcmc = model_ttt,
        name = "phiT_pT_psiT",
        myConstants = myConstants_ttt,
        runData = runData_ttt
      )
    
    
    
  ) # tar_plan
#tmp <- tar_read(modelOut_ttt)
#save(tmp, file = paste0('./models/runsOut/mcmc_ttt_', substr(Sys.time(),1,13), '.RData'))

#MCMCplot(object = modelOut_ttt$mcmc, params = "betaPhiOut")
