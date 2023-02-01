tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis"))

tt_trib = 
  tar_plan(
    
    tt_inputData_trib =
      list(
        y = target_eh_trib$sizeStateMatrix,
        deltaProps = as.numeric(
          table(target_eh_trib$stateMatrix[target_eh_trib$stateMatrix > 0]) / 
                           length(target_eh_trib$stateMatrix[target_eh_trib$stateMatrix > 0])
          ),
        first = target_eh_trib$first, 
        last = target_eh_trib$last, 
        zInits = getInits_tt_trib(target_eh_trib$eh),
        
        nStates = length(unique(target_eh_trib$data$sizeState)),
        nRivers = length(unique(target_eh_trib$data$sizeState)) # for now
        
        
      ),   
    
    tt_runData_trib = list(
      # Updateable model-specific variables 
      nIter = 5000, 
      nBurnin = 2000, 
      nChains = 2,
      thinRate = 5
    ),

    tt_alpha_trib = list(
      alphaR1 = c(0.8, 0.1, 0.1), 
      alphaR2 = c(0.1, 0.8, 0.1),
      alphaR3 = c(0.1, 0.1, 0.8)
    ),
    
    tt_myConstants_trib = list(
      N = nrow(tt_inputData_trib$y),
      T = ncol(tt_inputData_trib$y),
      first = tt_inputData_trib$first,
      last = tt_inputData_trib$last,

      nRivers = tt_inputData_trib$nRivers,
      length = tt_inputData_trib$last - tt_inputData_trib$first + 1,

      alphaR1 = tt_alpha_trib$alphaR1,
      alphaR2 = tt_alpha_trib$alphaR2,
      alphaR3 = tt_alpha_trib$alphaR3,

      deltaProps = tt_inputData_trib$deltaProps,
      nStates = tt_inputData_trib$nStates,
      
      dateMedianDiffMonth = target_eh_trib$dateMedian$dateMedianDiffMonth
    ),

    tt_myData_trib = list(
      y = tt_inputData_trib$y + 1
    ),

    tt_modelCode_trib = nimbleCode({
      # Initial distribution among rivers
      delta[1] <- deltaProps[1]                  # Pr(alive t = 1 and in river 1) = 0.4
      delta[2] <- deltaProps[2]
      delta[3] <- deltaProps[3]
      delta[4] <- 0                    # Pr(dead t = 1) = 0

      for (s in 1:nStates){
        for (t in 1:(T-1)){
          #betaPhi[r,t] ~ dnorm(0, sd = 1)
          #betaP[r,t] ~ dnorm(0, sd = 1)
          betaPhi[s,t] ~ dunif(0,1)
          betaP[s,t] ~ dunif(0,1)
          
          betaPhiMonthly[r,t] <- betaPhi[r,t] ^ dateMedianDiffMonth[t]
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

    }),

    tt_Rmodel_trib = nimbleModel(
      code = tt_modelCode_trib,
      constants = tt_myConstants_trib,
      data = tt_myData_trib,
      inits = initialValues_tt_trib(tt_inputData_trib$nStates, tt_myConstants_trib$T, tt_alpha_trib),
      calculate = FALSE
    ),

    tt_parametersToSave_trib = c("betaPhi", "betaP", "psi", "betaPhiMonthly"),
    
    tt_conf_trib = configureMCMC(
      tt_Rmodel_trib,
      monitors = tt_parametersToSave_trib
    ),

    tt_Rmcmc_trib = buildMCMC(tt_conf_trib, useConjugacy = FALSE),
    tt_Cmodel_trib = compileNimble(tt_Rmodel_trib),
    tt_Cmcmc_trib = compileNimble(tt_Rmcmc_trib, project = tt_Rmodel_trib),

    tt_model_trib = runMCMC(
      tt_Cmcmc_trib,
      niter = tt_runData_trib$nIter,
      nburnin = tt_runData_trib$nBurnin,
      thin = tt_runData_trib$thinRate,
      nchains = tt_runData_trib$nChains
    ),

    tt_modelOut_trib =
      list(
        mcmc = tt_model_trib, # "Error : invalid nimbleFunction argument"
        name = "phiT_pT_trib",
        modelCode = tt_modelCode_trib,
        myConstants = tt_myConstants_trib,
        runData = tt_runData_trib
      ),

    tt_save_trib = saveModelOut_tt_trib(tt_modelOut_trib)

  ) # tar_plan

############ Functions ###############

saveModelOut_tt_trib <- function(d) {
  save(d, file = paste0('./models/runsOut/tt_trib_', substr(Sys.time(),1,13), '.RData'))
}

initialValues_tt_trib <- function(s, t, a) {
  list(
    betaPhi = array(runif(s * (t - 1), 0, 1), c(s, (t - 1))),
    betaP =   array(runif(s * (t - 1), 0, 1), c(s, (t - 1))),
    
    psi = getDirchPriorsR_tt_trib(s, t, a)
  )
}

getDirchPriorsR_tt_trib <- function(s, tt, aa){
  a <- array(rep(0, s * s * (tt - 1)), c(s, s, (tt - 1)))
  for(t in 1:(tt - 1)){
    a[1,1:s,t] <- rdirch(1, aa$alphaR1)
    a[2,1:s,t] <- rdirch(1, aa$alphaR2)
    a[3,1:s,t] <- rdirch(1, aa$alphaR3)
  }
  return(a)
}

getInits_tt_trib <- function(y) {
  zInits <- y + 1 # non-detection -> alive
  zInits[zInits == 2] = 1 # dead -> alive
  return(zInits)
}
