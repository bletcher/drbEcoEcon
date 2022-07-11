tar_option_set(packages = c("tidyverse", "nimble"))

# Updateable model-specific variables 
nIter = 250 
nBurnin = 100 
nChains = 2 
thinRate = 5

target_CMR_models_phiT_pT_psiT = 
  tar_plan(
    
    # inputVars_ttt = list(
    #   nIter = nIter, 
    #   nBurnin = nBurnin, 
    #   nChains = nChains, 
    #   thinRate = thinRate
    # ),   
    
    inputData_ttt = list(
      y = target_y,
      first = target_first, 
      last = target_last, 
      zInits = target_zInits 
    ),
    
    myConstants_ttt = list(
      N = nrow(target_y),
      T = ncol(target_y),
      first = target_first,
      last = target_last
    ),
    
    myData_ttt = list(
      y = target_y + 1
    ),
    
    initialValues_ttt = function() list(phi = runif(myConstants_ttt$T - 1, 0, 1),
                                        p = runif(myConstants_ttt$T - 1, 0, 1),
                                        z = target_zInits
    ),
    
    parametersToSave_ttt = c("phi", "p"),#, "psi"),
    
    modelCode_ttt = nimbleCode({
      delta[1] <- 1                    # Pr(alive t = 1) = 1
      delta[2] <- 0                    # Pr(dead t = 1) = 0
      for (t in 1:(T-1)){                # loop over time
        phi[t] ~ dunif(0, 1)           # prior survival
        gamma[1,1,t] <- phi[t]         # Pr(alive t -> alive t+1)
        gamma[1,2,t] <- 1 - phi[t]     # Pr(alive t -> dead t+1)
        gamma[2,1,t] <- 0              # Pr(dead t -> alive t+1)
        gamma[2,2,t] <- 1              # Pr(dead t -> dead t+1)
        p[t] ~ dunif(0, 1)             # prior detection
        omega[1,1,t] <- 1 - p[t]       # Pr(alive t -> non-detected t)
        omega[1,2,t] <- p[t]           # Pr(alive t -> detected t)
        omega[2,1,t] <- 1              # Pr(dead t -> non-detected t)
        omega[2,2,t] <- 0              # Pr(dead t -> detected t)
      }
      # likelihood
      for (i in 1:N){
        z[i,first[i]] ~ dcat(delta[1:2])
        for (j in (first[i]+1):last[i]){
          z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
          y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
        }
      }
    }),
    
    Rmodel_ttt = nimbleModel(
      code = modelCode_ttt,
      constants = myConstants_ttt,
      data = myData_ttt,
      inits = initialValues_ttt(),
      calculate = FALSE
    ),
    conf_ttt = configureMCMC(
      Rmodel_ttt,
      monitors = parametersToSave_ttt
    ),
    
    Rmcmc_ttt = buildMCMC(conf_ttt, useConjugacy = FALSE),
    Cmodel_ttt = compileNimble(Rmodel_ttt),
    Cmcmc_ttt = compileNimble(Rmcmc_ttt, project = Rmodel_ttt),
    
    modelOut_ttt = runMCMC(
      Cmcmc_ttt,
      niter = nIter,
      nburnin = nBurnin,
      thin = thinRate,
      nchains = nChains
    )
    
  )
    # modelOut_ttt = runModel(
    #   modelCode_ttt,
    #   myConstants_ttt,
    #   myData_ttt,
    #   initialValues_ttt(),
    #   parametersToSave_ttt,
    #   inputVars_ttt
    # )