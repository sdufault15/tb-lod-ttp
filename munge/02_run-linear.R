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
                        normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                        prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_mams_30,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-mams-lod-30.RData")))
rm(m_linear_mams_30, df_analysis_mams)
beepr::beep()

load(here("data", "cleaned-data", "2024-02-16_mams-clean.RData"))
m_linear_mams_25 <- brm(log10(dtp_25) | cens(censored_25) ~ weeks + (1 + weeks | patient.id + Treatm_arm), # run the model
                        data = df_analysis_mams,
                        inits = 0,
                        chains = 4,
                        control = list(adapt_delta = 0.99),
                        iter = 4000,
                        backend = "cmdstanr", # attempt to improve convergence speed
                        normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                        prior = prior(normal(0,4), class = "b"))

save(m_linear_mams_25,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-mams-lod-25.RData")))
rm(m_linear_mams_25, df_analysis_mams)
beepr::beep()

###########################
# Run the models on REMox-TB
###########################
load(here("data", "cleaned-data", "2023-08-15_remoxtb-clean.RData"))
# Model took 4.4 hours to run locally
m_linear_remox_42 <- brm(log10(dtp_42) | cens(censored_42) ~ weeks + (1 + weeks | trial_no + treat), # run the model
                         data = df_analysis_remox,
                         inits = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
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
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_remox_30,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-remox-lod-30.RData")))

m_linear_remox_25 <- brm(log10(dtp_25) | cens(censored_25) ~ weeks + (1 + weeks | trial_no + treat), # run the model
                         data = df_analysis_remox,
                         inits = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_remox_25,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-remox-lod-25.RData")))

###########################
# Run the models on NC-002
###########################
load(here("data", "cleaned-data", "2024-03-26_NC-002-clean.RData"))
# Model took 4.4 hours to run locally
m_linear_nc002_42 <- brm(log10(dtp_42) | cens(censored_42) ~ weeks + (1 + weeks | USUBJID + ACTARM), # run the model
                         data = df_analysis_nc002,
                         init = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))
save(m_linear_nc002_42,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-nc002-lod-42.RData")))
rm(m_linear_nc002_42)
beepr::beep()

m_linear_nc002_30 <- brm(log10(dtp_30) | cens(censored_30) ~ weeks + (1 + weeks | USUBJID + ACTARM), # run the model
                         data = df_analysis_nc002,
                         init = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_nc002_30,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-nc002-lod-30.RData")))

m_linear_nc002_25 <- brm(log10(dtp_25) | cens(censored_25) ~ weeks + (1 + weeks | USUBJID + ACTARM), # run the model
                         data = df_analysis_nc002,
                         init = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_nc002_25,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-nc002-lod-25.RData")))

###########################
# Run the models on NC-005
###########################
load(here("data", "cleaned-data", "2024-03-26_NC-005-clean.RData"))
m_linear_nc005_42 <- brm(log10(dtp_42) | cens(censored_42) ~ weeks + (1 + weeks | patient.id + ACTARM), # run the model
                         data = df_analysis_nc005,
                         init = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))
save(m_linear_nc005_42,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-nc005-lod-42.RData")))
rm(m_linear_nc005_42)

m_linear_nc005_30 <- brm(log10(dtp_30) | cens(censored_30) ~ weeks + (1 + weeks | patient.id + ACTARM), # run the model
                         data = df_analysis_nc005,
                         init = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_nc005_30,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-nc005-lod-30.RData")))

m_linear_nc005_25 <- brm(log10(dtp_25) | cens(censored_25) ~ weeks + (1 + weeks | patient.id + ACTARM), # run the model
                         data = df_analysis_nc005,
                         init = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_nc005_25,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-nc005-lod-25.RData")))


###########################
# Run the models on NC-006
###########################
load(here("data", "cleaned-data", "2024-04-15_NC-006-clean.RData"))
m_linear_nc006_42 <- brm(log10(dtp_42) | cens(censored_42) ~ weeks + (1 + weeks | patient.id + ACTARM), # run the model
                         data = df_analysis_nc006,
                         init = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))
save(m_linear_nc006_42,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-nc006-lod-42.RData")))
rm(m_linear_nc006_42)

m_linear_nc006_30 <- brm(log10(dtp_30) | cens(censored_30) ~ weeks + (1 + weeks | patient.id + ACTARM), # run the model
                         data = df_analysis_nc006,
                         init = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_nc006_30,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-nc006-lod-30.RData")))

m_linear_nc006_25 <- brm(log10(dtp_25) | cens(censored_25) ~ weeks + (1 + weeks | patient.id + ACTARM), # run the model
                         data = df_analysis_nc006,
                         init = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_nc006_25,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-nc006-lod-25.RData")))

###########################
# Run the models on Study 29
###########################
load(here("data", "cleaned-data", "2024-03-26_TBTC-S29-clean.RData"))
# Model took 4.4 hours to run locally
m_linear_s29_42 <- brm(log10(dtp_42) | cens(censored_42) ~ weeks + (1 + weeks | patient.id + ACTARM), # run the model
                         data = df_analysis_s29,
                         init = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))
save(m_linear_s29_42,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-s29-lod-42.RData")))
rm(m_linear_s29_42)
beepr::beep()

m_linear_s29_30 <- brm(log10(dtp_30) | cens(censored_30) ~ weeks + (1 + weeks | patient.id + ACTARM), # run the model
                         data = df_analysis_s29,
                         init = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_s29_30,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-s29-lod-30.RData")))

m_linear_s29_25 <- brm(log10(dtp_25) | cens(censored_25) ~ weeks + (1 + weeks | patient.id + ACTARM), # run the model
                         data = df_analysis_s29,
                         init = 0,
                         chains = 4,
                         control = list(adapt_delta = 0.99),
                         iter = 4000,
                         backend = "cmdstanr", # attempt to improve convergence speed
                         normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                         prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_s29_25,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-s29-lod-25.RData")))

###########################
# Run the models on Study 29x
###########################
load(here("data", "cleaned-data", "2024-03-26_TBTC-S29x-clean.RData"))
# Model took 4.4 hours to run locally
m_linear_s29x_42 <- brm(log10(dtp_42) | cens(censored_42) ~ weeks + (1 + weeks | patient.id + ACTARM), # run the model
                       data = df_analysis_s29x,
                       init = 0,
                       chains = 4,
                       control = list(adapt_delta = 0.99),
                       iter = 4000,
                       backend = "cmdstanr", # attempt to improve convergence speed
                       normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                       prior = prior(normal(0,4), class = "Intercept"))
save(m_linear_s29x_42,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-s29x-lod-42.RData")))
rm(m_linear_s29x_42)
beepr::beep()

m_linear_s29x_30 <- brm(log10(dtp_30) | cens(censored_30) ~ weeks + (1 + weeks | patient.id + ACTARM), # run the model
                       data = df_analysis_s29x,
                       init = 0,
                       chains = 4,
                       control = list(adapt_delta = 0.99),
                       iter = 4000,
                       backend = "cmdstanr", # attempt to improve convergence speed
                       normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                       prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_s29x_30,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-s29x-lod-30.RData")))

m_linear_s29x_25 <- brm(log10(dtp_25) | cens(censored_25) ~ weeks + (1 + weeks | patient.id + ACTARM), # run the model
                       data = df_analysis_s29x,
                       init = 0,
                       chains = 4,
                       control = list(adapt_delta = 0.99),
                       iter = 4000,
                       backend = "cmdstanr", # attempt to improve convergence speed
                       normalize = FALSE, # attempt to improve speed of convergence (https://discourse.mc-stan.org/t/faster-convergence/21532)
                       prior = prior(normal(0,4), class = "Intercept"))

save(m_linear_s29x_25,
     file = here("data", "model-generated", 
                 paste0(Sys.Date(), "_linear-s29x-lod-25.RData")))
