tar_option_set(packages = c("tidyverse", "nimble", "nimbleEcology", "MCMCvis"))

target_ttt = 
  tar_plan(
    
    ttt_nStates = length(unique(target_eh$data$state)),
    ttt_nRivers = ttt_nStates, # for now
    #ttt_first = target_eh$first, #apply(y, 1, function(x) min(which(x !=0)))
    #ttt_last = target_eh$last,
    
    # Priors for psi where more likely to stay than move
    # states
    
    # River      size1  size2   size3
    # Main        1       2       3
    # Trib        4       5       6
    ttt_a = list(
      same = 0.750,
      shrink = 0.01,
      grow1 = 0.025,
      grow2 = 0.001,
      move = 0.025
    ),
    
    ttt_alpha = list(
      alphaR1 = c(ttt_a$same,   ttt_a$grow1,  ttt_a$grow2, ttt_a$move,   ttt_a$move,    ttt_a$move), 
      alphaR2 = c(ttt_a$shrink, ttt_a$same,   ttt_a$grow1, ttt_a$shrink, ttt_a$move,    ttt_a$move),
      alphaR3 = c(ttt_a$shrink, ttt_a$shrink, ttt_a$same,  ttt_a$shrink, ttt_a$shrink,  ttt_a$move),
      alphaR4 = c(ttt_a$move,   ttt_a$move,   ttt_a$move,  ttt_a$same,   ttt_a$grow1,   ttt_a$grow2),
      alphaR5 = c(ttt_a$shrink, ttt_a$move,   ttt_a$move,  ttt_a$shrink, ttt_a$same,    ttt_a$grow1),
      alphaR6 = c(ttt_a$shrink, ttt_a$shrink, ttt_a$move,  ttt_a$shrink, ttt_a$shrink,  ttt_a$same)
    ),
    
    ttt_inputData = list(
      y = target_eh$stateMatrix,
      first = target_eh$first, 
      last = target_eh$last, 
      zInits = getInits(target_eh$eh),
      
      # Proportion of fish in each river on the first observation
      deltaProps = table(target_eh$stateMatrix[target_eh$stateMatrix > 0]) / length(target_eh$stateMatrix[target_eh$stateMatrix > 0])
    ),   
    
    ttt_runData = list(
      # Updateable model-specific variables 
      nIter = 2500, 
      nBurnin = 1000, 
      nChains = 2,
      thinRate = 5
    ),
    
    ttt_myConstants = list(
      N = nrow(ttt_inputData$y),
      T = ncol(ttt_inputData$y),
      first = target_eh$first,
      last = target_eh$last,
      
      nRivers = ttt_nRivers,
      length = target_eh$last - target_eh$first + 1,
      
      alphaR1 = ttt_alpha$alphaR1,
      alphaR2 = ttt_alpha$alphaR2,
      alphaR3 = ttt_alpha$alphaR3,
      alphaR4 = ttt_alpha$alphaR4,
      alphaR5 = ttt_alpha$alphaR5,
      alphaR6 = ttt_alpha$alphaR6,
      
      deltaProps = ttt_inputData$deltaProps,
      nStates = ttt_nStates
    ),
    
    ttt_myData = list(
      y = ttt_inputData$y + 1
    ),
    
    ttt_parametersToSave = c("betaPhi", "betaPhiRiver", 
                             "betaP",   "betaPRiver", 
                             "betaPhiOut", "betaPhiRiverOut",
                             "betaPOut",   "betaPRiverOut", 
                             
                             "psi"),
    
    ttt_modelCode = nimbleCode({
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
          psi[5,1:nRivers,t] ~ ddirch(alphaR5[1:nRivers])
          psi[6,1:nRivers,t] ~ ddirch(alphaR6[1:nRivers])
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
        } # t loop
      } # i loop

      for (i in 1:N){
        y[i,first[i]:last[i]] ~ dDHMMo(init = delta[1:7],
                                       probTrans = gamma[1:7, 1:7, first[i]:last[i]],
                                       probObs =   omega[1:7, 1:7, first[i]:last[i], i],
                                       len = length[i],
                                       checkRowSums = 1)
      }

    }),

    ttt_Rmodel = nimbleModel(
      code = ttt_modelCode,
      constants = ttt_myConstants,
      data = ttt_myData$y,
      inits = initialValues(ttt_nRivers, ttt_myConstants, ttt_alpha),
      calculate = FALSE
    ),
    
    ttt_conf = configureMCMC(
      ttt_Rmodel,
      monitors = ttt_parametersToSave
    ),

    ttt_Rmcmc = buildMCMC(ttt_conf, useConjugacy = FALSE),
    ttt_Cmodel= compileNimble(ttt_Rmodel),
    ttt_Cmcmc = compileNimble(ttt_Rmcmc, project = ttt_Rmodel),

    ttt_model = runMCMC(
      ttt_Cmcmc,
      niter = ttt_runData$nIter,
      nburnin = ttt_runData$nBurnin,
      thin = ttt_runData$thinRate,
      nchains = ttt_runData$nChains
    ),

    ttt_modelOut =
      list(
        mcmc = ttt_model,
        name = "phiT_pT_psiT",
        myConstants = ttt_myConstants,
        runData = ttt_runData
      )
    
    
    
  ) # tar_plan

############ Functions ###############

initialValues <- function(r, c, aa) {  
  list(
    betaPhiRiver = array(runif(r, 0, 1), c(r)),
    betaPhi = array(rnorm(r * (c$T - 1), 0, 1), c(r, (c$T - 1))),
     
    betaPRiver = array(runif(r, 0, 1), c(r)),
    betaP = array(rnorm(r * (c$T - 1), 0, 1), c(r, (c$T - 1))),
    
    psi = getDirchPriorsR(r, c, aa)
  )
}

getDirchPriorsR <- function(s, c, aa){
  a <- array(rep(0, s * s * (c$T - 1)), c(s, s, (c$T - 1)))
  #for(r in 1:s){
    for(t in 1:(c$T - 1)){
      a[1,1:s,t] <- rdirch(1, aa$alphaR1)
      a[2,1:s,t] <- rdirch(1, aa$alphaR2)
      a[3,1:s,t] <- rdirch(1, aa$alphaR3)
      a[4,1:s,t] <- rdirch(1, aa$alphaR4)
      a[5,1:s,t] <- rdirch(1, aa$alphaR5)
      a[6,1:s,t] <- rdirch(1, aa$alphaR6)
    }
  #}
  return(a)
}

getInits <- function(y) {
  zInits <- y + 1 # non-detection -> alive
  zInits[zInits == 2] = 1 # dead -> alive
  return(zInits)
}
