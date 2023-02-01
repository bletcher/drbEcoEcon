tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis"))

tt_main = 
  tar_plan(
    
    tt_inputData_main =
      list(
        y = target_eh_main$sizeStateMatrix,
        deltaProps = as.numeric(
          table(target_eh_main$stateMatrix[target_eh_main$stateMatrix > 0]) / 
                           length(target_eh_main$stateMatrix[target_eh_main$stateMatrix > 0])
          ),
        first = target_eh_main$first, 
        last = target_eh_main$last, 
        zInits = getInits_tt_main(target_eh_main$eh),
        
        nStates = length(unique(target_eh_main$data$sizeState)),
        nRivers = length(unique(target_eh_main$data$sizeState)), # for now
        
        dummy = 0 # For resetting targets stream
      ),   
    
    tt_runData_main = list(
      # Updateable model-specific variables 
      nIter = 5000, 
      nBurnin = 2000, 
      nChains = 2,
      thinRate = 5
    ),

    tt_alpha_main = list(
      alphaR1 = c(0.8,0.1,0.1),
      alphaR2 = c(0.1,0.8,0.1),
      alphaR3 = c(0.1,0.1,0.8)
    ),
    
    tt_myConstants_main = list(
      N = nrow(tt_inputData_main$y),
      T = ncol(tt_inputData_main$y),
      first = tt_inputData_main$first,
      last = tt_inputData_main$last,

      nRivers = tt_inputData_main$nRivers,
      length = tt_inputData_main$last - tt_inputData_main$first + 1,

      alphaR1 = tt_alpha_main$alphaR1,
      alphaR2 = tt_alpha_main$alphaR2,
      alphaR3 = tt_alpha_main$alphaR3,

      deltaProps = tt_inputData_main$deltaProps,
      nStates = tt_inputData_main$nStates,
      
     # dateMedianDiffMonth = target_eh_main$dateMedian$dateMedianDiffMonth,
     dateMedianDiffMonthVector = target_eh_main$dateMedianDiffMonthVector
    ),

    tt_myData_main = list(
      y = tt_inputData_main$y + 1
    ),

    tt_modelCode_main = nimbleCode({
      # Initial distribution among rivers
      delta[1] <- deltaProps[1]                  # Pr(alive t = 1 and in river 1) = 0.4
      delta[2] <- deltaProps[2]
      delta[3] <- deltaProps[3]
      delta[4] <- 0                    # Pr(dead t = 1) = 0

      for (r in 1:nStates){
        for (t in 1:(T-1)){
          #betaPhi[r,t] ~ dnorm(0, sd = 1)
          #betaP[r,t] ~ dnorm(0, sd = 1)
          betaPhi[r,t] ~ dunif(0,1)
          betaP[r,t] ~ dunif(0,1)

          betaPhiMonthly[r,t] <- betaPhi[r,t] ^ dateMedianDiffMonthVector[t]
          #betaPhiOut[r,t] <- ilogit(betaPhi[r,t])
          #betaPOut[r,t] <- ilogit(betaP[r,t])
        }
      }

      for (t in 1:(T-1)){ # loop over time
        psi[1,1:nStates,t] ~ ddirch(alphaR1[1:nStates])
        psi[2,1:nStates,t] ~ ddirch(alphaR2[1:nStates])
        psi[3,1:nStates,t] ~ ddirch(alphaR3[1:nStates])
        
        gamma[1,1,t] <- (betaPhi[1,t]) * psi[1,1,t]
        gamma[1,2,t] <- (betaPhi[1,t]) * psi[1,2,t]
        gamma[1,3,t] <- (betaPhi[1,t]) * psi[1,3,t]
        gamma[1,4,t] <- 1 - (betaPhi[1,t])
        
        gamma[2,1,t] <- (betaPhi[2,t]) * psi[2,1,t]
        gamma[2,2,t] <- (betaPhi[2,t]) * psi[2,2,t]
        gamma[2,3,t] <- (betaPhi[2,t]) * psi[2,3,t]
        gamma[2,4,t] <- 1 - (betaPhi[2,t])
        
        gamma[3,1,t] <- (betaPhi[3,t]) * psi[3,1,t]
        gamma[3,2,t] <- (betaPhi[3,t]) * psi[3,2,t]
        gamma[3,3,t] <- (betaPhi[3,t]) * psi[3,3,t]
        gamma[3,4,t] <- 1 - (betaPhi[3,t])
        
        gamma[4,1,t] <- 0
        gamma[4,2,t] <- 0
        gamma[4,3,t] <- 0
        gamma[4,4,t] <- 1

      }

      # gamma for the last occasion
      for (a in 1:(nStates+1)){
        for (b in 1:nStates){
          gamma[a,b,T] <- 0
        }
        gamma[a,4,T] <- 1
      }

      for (i in 1:N){ # loop over individuals
        # omega for first obs
        omega[1,1,first[i],i] <- 0          # Pr(alive A t -> non-detected t)
        omega[1,2,first[i],i] <- 1          # Pr(alive A t -> detected A t)
        omega[1,3,first[i],i] <- 0          # Pr(alive A t -> detected B t)
        omega[1,4,first[i],i] <- 0

        omega[2,1,first[i],i] <- 0          # Pr(alive B t -> non-detected t)
        omega[2,2,first[i],i] <- 0          # Pr(alive B t -> detected A t)
        omega[2,3,first[i],i] <- 1          # Pr(alive B t -> detected B t)
        omega[2,4,first[i],i] <- 0

        omega[3,1,first[i],i] <- 0          # Pr(alive C t -> non-detected t)
        omega[3,2,first[i],i] <- 0          # Pr(alive C t -> detected A t)
        omega[3,3,first[i],i] <- 0          # Pr(alive C t -> detected B t)
        omega[3,4,first[i],i] <- 1

        omega[4,1,first[i],i] <- 1
        omega[4,2,first[i],i] <- 0
        omega[4,3,first[i],i] <- 0
        omega[4,4,first[i],i] <- 0


        ## DT changes:
        ## time t > first[i]:
        for(t in (first[i]+1):last[i]) {
          logit(pA[t,i]) <- betaP[1,t-1]
          logit(pB[t,i]) <- betaP[2,t-1]
          logit(pC[t,i]) <- betaP[3,t-1]

          # probabilities of y(t) given z(t)
          # omega[z, y, t, i]

          # z=1 = alive in River 1, z=2 = alive in River 2...z=5 = dead
          # y=1 = unobserved, y=2 = observed in River 1, y=3 = observed in River 2, etc

          omega[1,1,t,i] <- 1 - pA[t,i]     # Pr(alive A t -> non-detected t)
          omega[1,2,t,i] <- pA[t,i]         # Pr(alive A t -> detected A t)
          omega[1,3,t,i] <- 0               # Pr(alive A t -> detected B t)
          omega[1,4,t,i] <- 0

          omega[2,1,t,i] <- 1 - pB[t,i]     # Pr(alive B t -> non-detected t)
          omega[2,2,t,i] <- 0               # Pr(alive B t -> detected A t)
          omega[2,3,t,i] <- pB[t,i]         # Pr(alive B t -> detected B t)
          omega[2,4,t,i] <- 0

          omega[3,1,t,i] <- 1 - pC[t,i]     # Pr(alive C t -> non-detected t)
          omega[3,2,t,i] <- 0               # Pr(alive C t -> detected A t)
          omega[3,3,t,i] <- 0               # Pr(alive C t -> detected B t)
          omega[3,4,t,i] <- pC[t,i]         # Pr(alive C t -> detected C t)

          omega[4,1,t,i] <- 1               # Pr(dead t -> non-detected t)
          omega[4,2,t,i] <- 0               # Pr(dead t -> detected A t)
          omega[4,3,t,i] <- 0               # Pr(dead t -> detected B t)
          omega[4,4,t,i] <- 0               # Pr(dead t -> detected C t)

        } # t loop
      } # i loop
        
      for (i in 1:N){
        y[i,first[i]:last[i]] ~ dDHMMo(init = delta[1:4],
                                       probTrans = gamma[1:4, 1:4, first[i]:last[i]],
                                       probObs =   omega[1:4, 1:4, first[i]:last[i], i],
                                       len = length[i],
                                       checkRowSums = 1)
      }
      
      # # likelihood 
      # for (i in 1:N){
      #   # latent state at first capture
      #   z[i,first[i]] <- y[i,first[i]] - 1
      #   for (t in (first[i]+1):K){
      #     # z(t) given z(t-1)
      #     z[i,t] ~ dcat(gamma[z[i,t-1],1:4])
      #     # y(t) given z(t)
      #     y[i,t] ~ dcat(omega[z[i,t],1:4])
      #   }
      # }

    }),

    tt_Rmodel_main = nimbleModel(
      code = tt_modelCode_main,
      constants = tt_myConstants_main,
      data = tt_myData_main,
      inits = initialValues_tt_main(tt_inputData_main$nStates, tt_myConstants_main$T, tt_alpha_main),
      calculate = FALSE
    ),

    tt_parametersToSave_main = c("betaPhi", "betaP", "psi", "betaPhiMonthly"),
    
    tt_conf_main = configureMCMC(
      model = tt_Rmodel_main,
      monitors = tt_parametersToSave_main
    ),

    tt_Rmcmc_main = buildMCMC(tt_conf_main, useConjugacy = FALSE),
    tt_Cmodel_main = compileNimble(tt_Rmodel_main),
    tt_Cmcmc_main = compileNimble(tt_Rmcmc_main, project = tt_Rmodel_main),

    tt_model_main = runMCMC(
      tt_Cmcmc_main,
      niter = tt_runData_main$nIter,
      nburnin = tt_runData_main$nBurnin,
      thin = tt_runData_main$thinRate,
      nchains = tt_runData_main$nChains
    ),

    tt_modelOut_main =
      list(
        mcmc = tt_model_main,
        name = "phiT_pT_main",
        modelCode = tt_modelCode_main,
        myConstants = tt_myConstants_main,
        runData = tt_runData_main
      ),

    tt_save_main = saveModelOut_tt_main(tt_modelOut_main)

  ) # tar_plan

############ Functions ###############

saveModelOut_tt_main <- function(d) {
  save(d, file = paste0('./models/runsOut/tt_main_', substr(Sys.time(),1,13), '.RData'))
}

initialValues_tt_main <- function(s, t, a) {
  list(
    #betaPhi = array(rnorm(s * (t - 1), 0, 1), c(s, (t - 1))),
    #betaP =   array(rnorm(s * (t - 1), 0, 1), c(s, (t - 1))),
    
    betaPhi = array(runif(s * (t - 1), 0, 1), c(s, (t - 1))),
    betaP =   array(runif(s * (t - 1), 0, 1), c(s, (t - 1))),
    
    psi = getDirchPriorsR_tt_main(s, t, a)
  )
}

getDirchPriorsR_tt_main <- function(s, tt, aa){
  a <- array(rep(0, s * s * (tt - 1)), c(s, s, (tt - 1)))
  for(t in 1:(tt - 1)){
    a[1,1:s,t] <- rdirch(1, aa$alphaR1)
    a[2,1:s,t] <- rdirch(1, aa$alphaR2)
    a[3,1:s,t] <- rdirch(1, aa$alphaR3)
  }
  return(a)
}

getInits_tt_main <- function(y) {
  zInits <- y + 1 # non-detection -> alive
  zInits[zInits == 2] = 1 # dead -> alive
  return(zInits)
}
