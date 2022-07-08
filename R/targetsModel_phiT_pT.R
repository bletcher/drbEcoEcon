tar_option_set(packages = c("tidyverse", "nimble"))

target_CMR_models_phiT_pT = 
  tar_plan(
    nIter = 25000, 
    nBurnin = 10000, 
    nChains = 2, 
    thinRate = 5,
    #target_model_phiT_pT = run_CMR_models_phiT_pT(tar_read(target_eh),
    #                                              nIter = nIter, 
    #                                              nBurnin = nBurnin, 
    #                                              nChains = nChains, 
    #                                              thinRate = thinRate)
   
    eh = tar_read(target_eh),
    y = eh$eh,
    first = eh$first, 
    last = eh$last, # check - this should be num Occ for all fish

    hmm.phiT_pT = nimbleCode({
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


    myConstants = list(N = nrow(y),
                       T = ncol(y),
                       first = first,
                       last = last),

    myData = list(y = y + 1),

    zinits = y + 1, # non-detection -> alive
    #zinits[zinits == 2] = 1, # dead -> alive

    initialValues = function() list(phi = runif(myConstants$T - 1, 0, 1),
                                     p = runif(myConstants$T - 1, 0, 1)
                                     #z = zinits
                                    ),

    parametersToSave = c("phi", "p"),

    #start = Sys.time(),

    Rmodel = nimbleModel(
      code = hmm.phiT_pT,
      constants = myConstants,
      data = myData,
      inits = initialValues(),
      calculate = FALSE
    ),
    conf = configureMCMC(
      Rmodel,
      monitors = parametersToSave
    ),

    Rmcmc = buildMCMC(conf, useConjugacy = FALSE),
    Cmodel = compileNimble(Rmodel),
    Cmcmc = compileNimble(Rmcmc, project = Rmodel),

    mcmc.phiT_pT = runMCMC(
      Cmcmc,
      niter = nIter,
      nburnin = nBurnin,
      thin = thinRate,
      nchains = nChains
    ),

    #end = Sys.time(),
    #elapsed_phiT_pT = end - start,

    # to return
    target_model_phiT_pT =
      list(
        mcmc = mcmc.phiT_pT,
       #elapsed = elapsed_phiT_pT,
        name = "phiT_pT",
        myConstants = myConstants,
        nIter = nIter,
        nBurnin = nBurnin,
        thinRate = thinRate,
        #nSeasons = nSeasons,
        #nCohorts = nCohorts,
        nChains = nChains
      )

  )
