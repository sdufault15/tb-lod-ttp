##########################################
# Run Linear Models on Data
##########################################
library(tidyverse)
library(here)
library(brms)
library(parallelly)
library(cmdstanr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallelly::availableCores())

###########################
# Run the models on MAMS-TB
###########################
load(here("data", "cleaned-data", "2023-07-20_mams-clean.RData"))
m_linear_mams_42 <- brm(log10(dtp_42) | cens(censored_42) ~ weeks + (1 + weeks | patient.id + Treatm_arm), # run the model
                        data = df_analysis_mams,
                        inits = 0,
                        chains = 4,
                        control = list(adapt_delta = 0.99),
                        iter = 2000,
                        backend = "cmdstanr", # attempt to improve convergence speed
                        save_pars=save_pars(group=FALSE), # attempt to decrease file size
                        normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                        prior = prior(normal(0,4), class = "Intercept"))
save(m_linear_mams_42,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-mams-lod-42.RData")))
rm(m_linear_mams_42)

m_linear_mams_30 <- brm(log10(dtp_30) | cens(censored_30) ~ weeks + (1 + weeks | patient.id + Treatm_arm), # run the model
                        data = df_analysis_mams,
                        inits = 0,
                        chains = 4,
                        control = list(adapt_delta = 0.99),
                        iter = 2000,
                        backend = "cmdstanr", # attempt to improve convergence speed
                        save_pars=save_pars(group=FALSE), # attempt to decrease file size
                        normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                        prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_mams_30,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-mams-lod-30.RData")))
rm(m_linear_mams_30, df_analysis_mams)
beepr::beep()

###########################
# Run the models on REMox-TB
###########################
load(here("data", "cleaned-data", "2023-07-20_remoxtb-clean.RData"))
# Model took 4.4 hours to run locally
m_linear_remox_42 <- brm(log10(dtp_42) | cens(censored_42) ~ weeks + (1 + weeks | trial_no + treat), # run the model
                         data = df_analysis_remox,
                         inits = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         save_pars=save_pars(group=FALSE), # attempt to decrease file size
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))
save(m_linear_remox_42,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-remox-lod-42.RData")))
rm(m_linear_remox_42)
beepr::beep()

m_linear_remox_30 <- brm(log10(dtp_30) | cens(censored_30) ~ weeks + (1 + weeks | trial_no + treat), # run the model
                         data = df_analysis_remox,
                         inits = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         save_pars=save_pars(group=FALSE), # attempt to decrease file size
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_remox_30,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-remox-lod-30.RData")))
beepr::beep()
