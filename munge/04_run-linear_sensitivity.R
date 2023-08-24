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
# Run the models on MAMS-TB - REMOVING BASELINE MEASUREMENT
###########################
load(here("data", "cleaned-data", "2023-08-15_mams-clean.RData"))
m_linear_mams_42s <- brm(log10(dtp_42) | cens(censored_42) ~ weeks + (1 + weeks | patient.id + Treatm_arm), # run the model
                        data = filter(df_analysis_mams, weeks != 0),
                        inits = 0,
                        chains = 4,
                        control = list(adapt_delta = 0.99),
                        iter = 2000,
                        backend = "cmdstanr", # attempt to improve convergence speed
                        save_pars=save_pars(group=FALSE), # attempt to decrease file size
                        normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                        prior = prior(normal(0,4), class = "Intercept"))
save(m_linear_mams_42s,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-mams-lod-42_sensitivity.RData")))
rm(m_linear_mams_42s)

m_linear_mams_30s <- brm(log10(dtp_30) | cens(censored_30) ~ weeks + (1 + weeks | patient.id + Treatm_arm), # run the model
                        data = filter(df_analysis_mams, weeks != 0),
                        inits = 0,
                        chains = 4,
                        control = list(adapt_delta = 0.99),
                        iter = 2000,
                        backend = "cmdstanr", # attempt to improve convergence speed
                        save_pars=save_pars(group=FALSE), # attempt to decrease file size
                        normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                        prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_mams_30s,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-mams-lod-30_sensitivity.RData")))
rm(m_linear_mams_30s)

m_linear_mams_25s <- brm(log10(dtp_25) | cens(censored_25) ~ weeks + (1 + weeks | patient.id + Treatm_arm), # run the model
                        data = filter(df_analysis_mams, weeks != 0),
                        inits = 0,
                        chains = 4,
                        control = list(adapt_delta = 0.99),
                        iter = 2000,
                        backend = "cmdstanr", # attempt to improve convergence speed
                        save_pars=save_pars(group=FALSE), # attempt to decrease file size
                        normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                        prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_mams_25s,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-mams-lod-25_sensitivity.RData")))
rm(m_linear_mams_25s, df_analysis_mams)

###########################
# Run the models on REMox-TB
###########################
load(here("data", "cleaned-data", "2023-07-31_remoxtb-clean.RData"))
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
