##########################################
# Run Complex Biphasic Models on Data
##########################################
library(tidyverse)
library(here)
library(brms)
library(parallelly)
library(cmdstanr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallelly::availableCores())

# load(here("data", "cleaned-data", "2023-07-11_remoxtb-clean.RData"))
load(here("data", "cleaned-data", "2023-08-15_mams-clean.RData"))

m_nonlinear_mams_42 <- brm(bf(log10(dtp_42) | cens(censored_42) ~ alpha + beta1*weeks + beta2*gamma*log((exp((weeks-kappa)/gamma) + exp(-(weeks-kappa)/gamma))/(exp(kappa/gamma) + exp(-kappa/gamma))),
                              alpha + beta1 + beta2 + kappa + gamma ~ 1 + (1|patient.id + Treatm_arm), #specifying individual-level and regimen-level effects
                              nl = TRUE), # specify non-linear
                           family = gaussian(),
                           inits = 0, # Starting here 
                           data = df_analysis_mams,
                           backend = "cmdstanr",
                           normalize = FALSE,
                           prior = c(prior(normal(0,4), nlpar = "alpha"), 
                                     prior(normal(0,4), nlpar = "beta1"),
                                     prior(normal(0,4), nlpar = "beta2"),
                                     set_prior("uniform(2,11)", lb = 2, ub = 11, nlpar = "kappa"), # from Burger + Schall
                                     set_prior("uniform(0.1,2)", lb = 0.1, ub = 2, nlpar = "gamma")), # from Burger + Schall
                           chains = 4,
                           save_pars = save_pars(group = FALSE),
                           control = list(adapt_delta = 0.99),
                           iter = 4000)

save(m_nonlinear_mams_42,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_biphasic-mams-lod-42.RData")))
rm(m_nonlinear_mams_42)

m_nonlinear_mams_30 <- brm(bf(log10(dtp_30) | cens(censored_30) ~ alpha + beta1*weeks + beta2*gamma*log((exp((weeks-kappa)/gamma) + exp(-(weeks-kappa)/gamma))/(exp(kappa/gamma) + exp(-kappa/gamma))),
                              alpha + beta1 + beta2 + kappa + gamma ~ 1 + (1|patient.id + Treatm_arm), #specifying individual-level and regimen-level effects
                              nl = TRUE), # specify non-linear
                           family = gaussian(),
                           inits = "0", # Starting here 
                           data = filter(df_analysis_mams, DV != -99),
                           backend = "cmdstanr",
                           normalize = FALSE,
                           prior = c(prior(normal(0,4), nlpar = "alpha"), 
                                     prior(normal(0,4), nlpar = "beta1"),
                                     prior(normal(0,4), nlpar = "beta2"),
                                     set_prior("uniform(2,11)", lb = 2, ub = 11, nlpar = "kappa"), # from Burger + Schall
                                     set_prior("uniform(0.1,2)", lb = 0.1, ub = 2, nlpar = "gamma")), # from Burger + Schall
                           chains = 4,
                           save_pars = save_pars(group = FALSE),
                           control = list(adapt_delta = 0.99),
                           iter = 6000)
save(m_nonlinear_mams_30,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_biphasic-mams-lod-30.RData")))

################  
# REMoxTB
################  
load(here("data", "cleaned-data", "2023-07-31_remoxtb-clean.RData"))

m_nonlinear_remox_42 <- brm(bf(log10(dtp_42) | cens(censored_42) ~ alpha + beta1*weeks + beta2*gamma*log((exp((weeks-kappa)/gamma) + exp(-(weeks-kappa)/gamma))/(exp(kappa/gamma) + exp(-kappa/gamma))),
                              alpha + beta1 + beta2 + kappa + gamma ~ 1 + (1|trial_no + treat), #specifying individual-level and regimen-level effects
                              nl = TRUE), # specify non-linear
                           family = gaussian(),
                           inits = "0", # Starting here 
                           data = df_analysis_remox,
                           backend = "cmdstanr",
                           normalize = FALSE,
                           prior = c(prior(normal(0,4), nlpar = "alpha"), 
                                     prior(normal(0,4), nlpar = "beta1"),
                                     prior(normal(0,4), nlpar = "beta2"),
                                     set_prior("uniform(2,11)", lb = 2, ub = 11, nlpar = "kappa"), # from Burger + Schall
                                     set_prior("uniform(0.1,2)", lb = 0.1, ub = 2, nlpar = "gamma")), # from Burger + Schall
                           chains = 4,
                           save_pars = save_pars(group = FALSE),
                           control = list(adapt_delta = 0.99),
                           iter = 8000)

save(m_nonlinear_remox_42,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_biphasic-remox-lod-42.RData")))
rm(m_nonlinear_remox_42)

m_nonlinear_remox_30 <- brm(bf(log10(dtp_30) | cens(censored_30) ~ alpha + beta1*weeks + beta2*gamma*log((exp((weeks-kappa)/gamma) + exp(-(weeks-kappa)/gamma))/(exp(kappa/gamma) + exp(-kappa/gamma))),
                              alpha + beta1 + beta2 + kappa + gamma ~ 1 + (1|trial_no + treat), #specifying individual-level and regimen-level effects
                              nl = TRUE), # specify non-linear
                           family = gaussian(),
                           inits = "0", # Starting here 
                           data = df_analysis_remox,
                           backend = "cmdstanr",
                           normalize = FALSE,
                           prior = c(prior(normal(0,4), nlpar = "alpha"), 
                                     prior(normal(0,4), nlpar = "beta1"),
                                     prior(normal(0,4), nlpar = "beta2"),
                                     set_prior("uniform(2,11)", lb = 2, ub = 11, nlpar = "kappa"), # from Burger + Schall
                                     set_prior("uniform(0.1,2)", lb = 0.1, ub = 2, nlpar = "gamma")), # from Burger + Schall
                           chains = 4,
                           save_pars = save_pars(group = FALSE),
                           control = list(adapt_delta = 0.99),
                           iter = 8000)
save(m_nonlinear_remox_30,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_biphasic-remox-lod-30.RData")))
