tar_option_set(packages = c("tidyverse", "nimble"))

# Updateable model-specific variables 
nIter = 250 
nBurnin = 100 
nChains = 2 
thinRate = 5

nStates = 6
#(nStates <- length(unique(tar_read(target_eh)$data$riverN)))# rivers 1:4

target_CMR_models_phiT_pT_psiT = 
  tar_plan(
    
    # inputVars_ttt = list(
    #   nIter = nIter, 
    #   nBurnin = nBurnin, 
    #   nChains = nChains, 
    #   thinRate = thinRate
    # ),   
    
    inputData_ttt = list(
      y = tar_read(target_stateMatrix),
      first = tar_read(target_first), 
      last = tar_read(target_last), 
      zInits = tar_read(target_zInits) 
    ),
    
    myConstants_ttt = list(
      N = nrow(tar_read(target_y)),
      T = ncol(tar_read(target_y)),
      first = tar_read(target_first),
      last = tar_read(target_last)
    ),
    
    myData_ttt = list(
      y = tar_read(target_y) + 1
    ),
    
    initialValues_ttt = function() list(phi = runif(myConstants_ttt$T - 1, 0, 1),
                                        p = runif(myConstants_ttt$T - 1, 0, 1),
                                        z = tar_read(target_zInits)
    ),
    
    parametersToSave_ttt = c("phi", "p", "psi"),
    
    # Proportion of fish in each river on the first observation
    y1 = inputData_ttt$y,
    deltaProps = table(y1[y1>0]) / length(y1[y1>0])
    
    # modelCode_ttt = nimbleCode({
    #   # Initial distribution among rivers
    #   delta[1] <- deltaProps[1]                  # Pr(alive t = 1 and in river 1) = 0.4
    #   delta[2] <- deltaProps[2]
    #   delta[3] <- deltaProps[3]
    #   delta[4] <- deltaProps[4]
    #   delta[5] <- deltaProps[5]
    #   delta[6] <- deltaProps[6]
    #   delta[7] <- 0                    # Pr(dead t = 1) = 0
    #   
    #   for (r in 1:nRivers){
    #     betaPhiRiver[r] ~ dnorm(0,1)
    #     betaPRiver[r] ~ dnorm(0,1)
    #     
    #     betaPhiRiverOut[r] <- ilogit(betaPhiRiver[r])
    #     betaPRiverOut[r] <- ilogit(betaPRiver[r])
    #     
    #     for (c in 1:nCohorts){
    #       betaPhiRiverCohort[r,c] ~ dnorm(betaPhiRiver[r],1)
    #       betaPRiverCohort[r,c] ~ dnorm(betaPRiver[r],1)
    #       
    #       betaPhiRiverCohortOut[r,c] <- ilogit(betaPhiRiverCohort[r,c])
    #       betaPRiverCohortOut[r,c] <- ilogit(betaPRiverCohort[r,c])
    #       for (t in 1:(T-1)){
    #         betaPhi[r,t,c] ~ dnorm(betaPhiRiverCohort[r,c],1)
    #         betaP[r,t,c] ~ dnorm(betaPRiverCohort[r,c],1)
    #         
    #         betaPhiOut[r,t,c] <- ilogit(betaPhi[r,t,c])
    #         betaPOut[r,t,c] <- ilogit(betaP[r,t,c])
    #         
    #         # move from river 'r' to one of river 1:nRivers
    #         # Nice description of effect of 'alpha' on probabilities:
    #         # https://stats.stackexchange.com/questions/244917/what-exactly-is-the-alpha-in-the-dirichlet-distribution
    #         psi[r,1:nRivers,t,c] ~ ddirch(alpha[1:nRivers])
    #       }
    #     }
    #   }
    #   
    #   for (t in 1:(T-1)){ # loop over time
    #     for (c in 1:nCohorts){
    #       
    #       gamma[1,1,t,c] <- ilogit(betaPhi[1,t,c]) * psi[1,1,t,c]
    #       gamma[1,2,t,c] <- ilogit(betaPhi[1,t,c]) * psi[1,2,t,c]
    #       gamma[1,3,t,c] <- ilogit(betaPhi[1,t,c]) * psi[1,3,t,c]
    #       gamma[1,4,t,c] <- ilogit(betaPhi[1,t,c]) * psi[1,4,t,c]
    #       gamma[1,5,t,c] <- 1 - ilogit(betaPhi[1,t,c])
    #       gamma[2,1,t,c] <- ilogit(betaPhi[2,t,c]) * psi[2,1,t,c]
    #       gamma[2,2,t,c] <- ilogit(betaPhi[2,t,c]) * psi[2,2,t,c]
    #       gamma[2,3,t,c] <- ilogit(betaPhi[2,t,c]) * psi[2,3,t,c]
    #       gamma[2,4,t,c] <- ilogit(betaPhi[2,t,c]) * psi[2,4,t,c]
    #       gamma[2,5,t,c] <- 1 - ilogit(betaPhi[2,t,c])
    #       gamma[3,1,t,c] <- ilogit(betaPhi[3,t,c]) * psi[3,1,t,c]
    #       gamma[3,2,t,c] <- ilogit(betaPhi[3,t,c]) * psi[3,2,t,c]
    #       gamma[3,3,t,c] <- ilogit(betaPhi[3,t,c]) * psi[3,3,t,c]
    #       gamma[3,4,t,c] <- ilogit(betaPhi[3,t,c]) * psi[3,4,t,c]
    #       gamma[3,5,t,c] <- 1 - ilogit(betaPhi[3,t,c])
    #       gamma[4,1,t,c] <- ilogit(betaPhi[4,t,c]) * psi[4,1,t,c]
    #       gamma[4,2,t,c] <- ilogit(betaPhi[4,t,c]) * psi[4,2,t,c]
    #       gamma[4,3,t,c] <- ilogit(betaPhi[4,t,c]) * psi[4,3,t,c]
    #       gamma[4,4,t,c] <- ilogit(betaPhi[4,t,c]) * psi[4,4,t,c]
    #       gamma[4,5,t,c] <- 1 - ilogit(betaPhi[4,t,c])
    #       gamma[5,1,t,c] <- 0
    #       gamma[5,2,t,c] <- 0
    #       gamma[5,3,t,c] <- 0
    #       gamma[5,4,t,c] <- 0
    #       gamma[5,5,t,c] <- 1
    #     }
    #   }
    #   
    #   # gamma for the last occasion  
    #   for (c in 1:nCohorts){
    #     for (a in 1:(nRivers+1)){
    #       for (b in 1:nRivers){
    #         gamma[a,b,T,c] <- 0
    #       }  
    #       gamma[a,5,T,c] <- 1
    #     }
    #   }
    #   
    #   
    #   for (i in 1:N){ # loop over individuals
    #     # omega for first obs      
    #     omega[1,1,first[i],i] <- 0          # Pr(alive A t -> non-detected t)
    #     omega[1,2,first[i],i] <- 1          # Pr(alive A t -> detected A t)
    #     omega[1,3,first[i],i] <- 0          # Pr(alive A t -> detected B t)
    #     omega[1,4,first[i],i] <- 0          # Pr(alive A t -> detected C t)
    #     omega[1,5,first[i],i] <- 0          # Pr(alive A t -> detected D t)
    #     omega[2,1,first[i],i] <- 0          # Pr(alive B t -> non-detected t)
    #     omega[2,2,first[i],i] <- 0          # Pr(alive B t -> detected A t)
    #     omega[2,3,first[i],i] <- 1          # Pr(alive B t -> detected B t)
    #     omega[2,4,first[i],i] <- 0          # Pr(alive B t -> detected C t)
    #     omega[2,5,first[i],i] <- 0          # Pr(alive B t -> detected C t)
    #     omega[3,1,first[i],i] <- 0          # Pr(alive C t -> non-detected t)
    #     omega[3,2,first[i],i] <- 0          # Pr(alive C t -> detected A t)
    #     omega[3,3,first[i],i] <- 0          # Pr(alive C t -> detected B t)
    #     omega[3,4,first[i],i] <- 1          # Pr(alive C t -> detected C t)
    #     omega[3,5,first[i],i] <- 0          # Pr(alive C t -> detected C t)
    #     omega[4,1,first[i],i] <- 0          # Pr(dead t -> non-detected t)
    #     omega[4,2,first[i],i] <- 0          # Pr(dead t -> detected A t)
    #     omega[4,3,first[i],i] <- 0          # Pr(dead t -> detected B t)
    #     omega[4,4,first[i],i] <- 0          # Pr(dead t -> detected C t)
    #     omega[4,5,first[i],i] <- 1          # Pr(dead t -> detected C t)
    #     omega[5,1,first[i],i] <- 1          # Pr(dead t -> non-detected t)
    #     omega[5,2,first[i],i] <- 0          # Pr(dead t -> detected A t)
    #     omega[5,3,first[i],i] <- 0          # Pr(dead t -> detected B t)
    #     omega[5,4,first[i],i] <- 0          # Pr(dead t -> detected C t)
    #     omega[5,5,first[i],i] <- 0          # Pr(dead t -> detected D t)
    #     
    #     
    #     ## DT changes:
    #     ## time t > first[i]:
    #     for(t in (first[i]+1):last[i]) {
    #       logit(pA[t,i]) <- betaP[1,t-1,cohort[i]]
    #       logit(pB[t,i]) <- betaP[2,t-1,cohort[i]]
    #       logit(pC[t,i]) <- betaP[3,t-1,cohort[i]]
    #       logit(pD[t,i]) <- betaP[4,t-1,cohort[i]]
    #       
    #       # probabilities of y(t) given z(t)
    #       # omega[z, y, t, i]
    #       
    #       # z=1 = alive in River 1, z=2 = alive in River 2...z=5 = dead
    #       # y=1 = unobserved, y=2 = observed in River 1, y=3 = observed in River 2, etc
    #       
    #       omega[1,1,t,i] <- 1 - pA[t,i]     # Pr(alive A t -> non-detected t)
    #       omega[1,2,t,i] <- pA[t,i]         # Pr(alive A t -> detected A t)
    #       omega[1,3,t,i] <- 0               # Pr(alive A t -> detected B t)
    #       omega[1,4,t,i] <- 0               # Pr(alive A t -> detected C t)
    #       omega[1,5,t,i] <- 0               # Pr(alive A t -> detected D t)
    #       omega[2,1,t,i] <- 1 - pB[t,i]     # Pr(alive B t -> non-detected t)
    #       omega[2,2,t,i] <- 0               # Pr(alive B t -> detected A t)
    #       omega[2,3,t,i] <- pB[t,i]         # Pr(alive B t -> detected B t)
    #       omega[2,4,t,i] <- 0               # Pr(alive B t -> detected C t)
    #       omega[2,5,t,i] <- 0               # Pr(alive B t -> detected C t)
    #       omega[3,1,t,i] <- 1 - pC[t,i]     # Pr(alive C t -> non-detected t)
    #       omega[3,2,t,i] <- 0               # Pr(alive C t -> detected A t)
    #       omega[3,3,t,i] <- 0               # Pr(alive C t -> detected B t)
    #       omega[3,4,t,i] <- pC[t,i]         # Pr(alive C t -> detected C t)
    #       omega[3,5,t,i] <- 0               # Pr(alive C t -> detected C t)
    #       omega[4,1,t,i] <- 1 - pD[t,i]     # Pr(alive D t -> non-detected t))
    #       omega[4,2,t,i] <- 0               # Pr(dead D t -> detected A t)
    #       omega[4,3,t,i] <- 0               # Pr(dead D t -> detected B t)
    #       omega[4,4,t,i] <- 0               # Pr(dead D t -> detected C t)
    #       omega[4,5,t,i] <- pD[t,i]         # Pr(alive D t -> detected D t)
    #       omega[5,1,t,i] <- 1               # Pr(dead t -> non-detected t)
    #       omega[5,2,t,i] <- 0               # Pr(dead t -> detected A t)
    #       omega[5,3,t,i] <- 0               # Pr(dead t -> detected B t)
    #       omega[5,4,t,i] <- 0               # Pr(dead t -> detected C t)
    #       omega[5,5,t,i] <- 0               # Pr(dead t -> detected D t)
    #     }
    #     
    #   } # i loop
    #   
    #   for (i in 1:N){
    #     y[i,first[i]:last[i]] ~ dDHMMo(init = delta[1:5],
    #                                    probTrans = gamma[1:5, 1:5, first[i]:last[i], cohort[i]],
    #                                    probObs =   omega[1:5, 1:5, first[i]:last[i], i],
    #                                    len = length[i],
    #                                    checkRowSums = 1)
    #   }
    #   
    # }),
    # 
    # Rmodel_ttt = nimbleModel(
    #   code = modelCode_ttt,
    #   constants = myConstants_ttt,
    #   data = myData_ttt,
    #   inits = initialValues_ttt(),
    #   calculate = FALSE
    # ),
    # conf_ttt = configureMCMC(
    #   Rmodel_ttt,
    #   monitors = parametersToSave_ttt
    # ),
    # 
    # Rmcmc_ttt = buildMCMC(conf_ttt, useConjugacy = FALSE),
    # Cmodel_ttt = compileNimble(Rmodel_ttt),
    # Cmcmc_ttt = compileNimble(Rmcmc_ttt, project = Rmodel_ttt),
    # 
    # modelOut_ttt = runMCMC(
    #   Cmcmc_ttt,
    #   niter = nIter,
    #   nburnin = nBurnin,
    #   thin = thinRate,
    #   nchains = nChains
    # )
    
  )
    # modelOut_ttt = runModel(
    #   modelCode_ttt,
    #   myConstants_ttt,
    #   myData_ttt,
    #   initialValues_ttt(),
    #   parametersToSave_ttt,
    #   inputVars_ttt
    # )