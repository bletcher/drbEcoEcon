tar_option_set(packages = c("tidyverse", "nimble"))

target_CMR_models_phiT_pT = 
  tar_plan(
    nIter_ttt = 2500, 
    nBurnin_ttt = 1000, 
    nChains_ttt = 2, 
    thinRate_ttt = 5,
    #target_model_phiT_pT = run_CMR_models_phiT_pT(tar_read(target_eh),
    #                                              nIter = nIter, 
    #                                              nBurnin = nBurnin, 
    #                                              nChains = nChains, 
    #                                              thinRate = thinRate)
   
    eh_ttt = tar_read(target_eh),
    y_ttt = eh_ttt$eh,
    first_ttt = eh_ttt$first, 
    last_ttt = eh_ttt$last, # check - this should be num Occ for all fish

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


    myConstants_ttt = list(N = nrow(y_ttt),
                       T = ncol(y_ttt),
                       first = first_ttt,
                       last = last_ttt),

    myData_ttt = list(y_ttt = y_ttt + 1),

    zinits_ttt = y_ttt + 1, # non-detection -> alive
    #zinits[zinits == 2] = 1, # dead -> alive

    initialValues_ttt = function() list(phi = runif(myConstants_ttt$T - 1, 0, 1),
                                     p = runif(myConstants_ttt$T - 1, 0, 1)
                                     #z = zinits
                                    ),

    parametersToSave_ttt = c("phi", "p", "psi"),

    #start = Sys.time(),

    Rmodel_ttt = nimbleModel(
      code = hmm.phiT_pT,
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

    mcmc.phiT_pT = runMCMC(
      Cmcmc_ttt,
      niter = nIter_ttt,
      nburnin = nBurnin_ttt,
      thin = thinRate_ttt,
      nchains = nChains_ttt
    ),

    #end = Sys.time(),
    #elapsed_phiT_pT = end - start,

    # to return
    target_model_phiT_pT =
      list(
        mcmc = mcmc.phiT_pT,
       #elapsed = elapsed_phiT_pT,
        name = "phiT_pT_psiT",
        myConstants = myConstants_ttt,
        nIter = nIter_ttt,
        nBurnin = nBurnin_ttt,
        thinRate = thinRate_ttt,
        #nSeasons = nSeasons,
        #nCohorts = nCohorts,
        nChains = nChains_ttt
      )

  )
