# tar_option_set(packages = c("tidyverse", "nimble"))
# 
# target_CMR_models_phiT_pT = 
#   tar_plan(
#     nIter_tt = 2500, 
#     nBurnin_tt = 1000, 
#     nChains_tt = 2, 
#     thinRate_tt = 5,
#     #target_model_phiT_pT = run_CMR_models_phiT_pT(tar_read(target_eh),
#     #                                              nIter = nIter, 
#     #                                              nBurnin = nBurnin, 
#     #                                              nChains = nChains, 
#     #                                              thinRate = thinRate)
#    
#     eh_tt = tar_read(target_eh),
#     y_tt = eh_tt$eh,
#     first_tt = eh_tt$first, 
#     last_tt = eh_tt$last, # check - this should be num Occ for all fish
# 
#     hmm.phiT_pT = nimbleCode({
#       delta[1] <- 1                    # Pr(alive t = 1) = 1
#       delta[2] <- 0                    # Pr(dead t = 1) = 0
#       for (t in 1:(T-1)){                # loop over time
#         phi[t] ~ dunif(0, 1)           # prior survival
#         gamma[1,1,t] <- phi[t]         # Pr(alive t -> alive t+1)
#         gamma[1,2,t] <- 1 - phi[t]     # Pr(alive t -> dead t+1)
#         gamma[2,1,t] <- 0              # Pr(dead t -> alive t+1)
#         gamma[2,2,t] <- 1              # Pr(dead t -> dead t+1)
#         p[t] ~ dunif(0, 1)             # prior detection
#         omega[1,1,t] <- 1 - p[t]       # Pr(alive t -> non-detected t)
#         omega[1,2,t] <- p[t]           # Pr(alive t -> detected t)
#         omega[2,1,t] <- 1              # Pr(dead t -> non-detected t)
#         omega[2,2,t] <- 0              # Pr(dead t -> detected t)
#       }
#       # likelihood
#       for (i in 1:N){
#         z[i,first[i]] ~ dcat(delta[1:2])
#         for (j in (first[i]+1):last[i]){
#           z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
#           y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
#         }
#       }
#     }),
# 
# 
#     myConstants_tt = list(N = nrow(y_tt),
#                        T = ncol(y_tt),
#                        first = first_tt,
#                        last = last_tt),
# 
#     myData_tt = list(y_tt = y_tt + 1),
# 
#     zinits_tt = y_tt + 1, # non-detection -> alive
#     #zinits[zinits == 2] = 1, # dead -> alive
# 
#     initialValues_tt = function() list(phi = runif(myConstants_tt$T - 1, 0, 1),
#                                      p = runif(myConstants_tt$T - 1, 0, 1)
#                                      #z = zinits
#                                     ),
# 
#     parametersToSave_tt = c("phi", "p"),
# 
#     #start = Sys.time(),
# 
#     Rmodel_tt = nimbleModel(
#       code = hmm.phiT_pT,
#       constants = myConstants_tt,
#       data = myData_tt,
#       inits = initialValues_tt(),
#       calculate = FALSE
#     ),
#     conf_tt = configureMCMC(
#       Rmodel_tt,
#       monitors = parametersToSave_tt
#     ),
# 
#     Rmcmc_tt = buildMCMC(conf_tt, useConjugacy = FALSE),
#     Cmodel_tt = compileNimble(Rmodel_tt),
#     Cmcmc_tt = compileNimble(Rmcmc_tt, project = Rmodel_tt),
# 
#     mcmc.phiT_pT = runMCMC(
#       Cmcmc_tt,
#       niter = nIter_tt,
#       nburnin = nBurnin_tt,
#       thin = thinRate_tt,
#       nchains = nChains_tt
#     ),
# 
#     #end = Sys.time(),
#     #elapsed_phiT_pT = end - start,
# 
#     # to return
#     target_model_phiT_pT =
#       list(
#         mcmc = mcmc.phiT_pT,
#        #elapsed = elapsed_phiT_pT,
#         name = "phiT_pT",
#         myConstants = myConstants_tt,
#         nIter = nIter_tt,
#         nBurnin = nBurnin_tt,
#         thinRate = thinRate_tt,
#         #nSeasons = nSeasons,
#         #nCohorts = nCohorts,
#         nChains = nChains_tt
#       )
# 
#   )
