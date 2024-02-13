#########################
# Suzanne M Dufault
# Data Cleaning
# Origin: 2023-06-16
# Revision: 2023-06-16
#########################


#########################
# REMoxTB
#########################
dta_patient <- read_dta(here("data", "cl_patientdata4sunita.dta"))
dta_culture <- read_dta(here("data","cl_culturedata4sunita.dta"))

dta_full_remox <- full_join(dta_patient,
                            dta_culture)

df_analysis_remox <- dta_full_remox %>% 
  filter(bact == 1, weeks <= 8, result < 6, !box_mitt %in% c(10,20,30)) %>% 
  # Add a censoring variable (if result is negative and is.na(DTP), then censored) 
  mutate(dtp_42 = ifelse(is.na(dtp), 42, dtp),
         censored_42 = ifelse(result == 0 & is.na(dtp), "right", "none")) %>% 
  # Adding in 30 day censoring (for modeling)
  mutate(dtp_30 = ifelse(dtp_42 >= 30, 30, dtp_42),
         censored_30 = ifelse(dtp_42 >= 30, "right", "none")) %>% 
  # Adding in 25 day censoring (for modeling)
  mutate(dtp_25 = ifelse(dtp_42 >= 25, 25, dtp_42),
         censored_25 = ifelse(dtp_42 >= 25, "right", "none")) 

save(df_analysis_remox,
     file = here("data", "cleaned-data", 
                 paste0(Sys.Date(), "_remoxtb-clean.RData")))

#########################
# PanACEA MAMS-TB
#########################
library(readr)
df_mams <- read_csv("data/MGITreg_MAMS_20170228_red.csv")
df_analysis_mams <- df_mams %>% 
  dplyr::select(arm = ARM, 
                patient.id = ID, 
                Treatm_arm,
                weeks = WEEK, 
                date = DATE,
                DV,
                POS) %>% 
  # CHECK WIH ELIN REGARDING DV = -99
  filter(weeks <= 8, !DV %in% c("-99", 
                                ".")) %>% 
  # FOLLOWING ELIN
  mutate(DV = ifelse(DV == "0", 5, as.numeric(DV))) %>%
  mutate(dtp = DV/24) %>% 
  mutate(dtp_42 = ifelse(POS == "0", 42, dtp),
         censored_42 = ifelse(POS == "0", "right", "none")) %>% 
  # Adding in 30 day censoring (for modeling)
  mutate(dtp_30 = ifelse(dtp_42 >= 30, 30, dtp_42),
         censored_30 = ifelse(dtp_42 >= 30, "right", "none")) %>% 
  # Adding in 25 day censoring (for modeling)
  mutate(dtp_25 = ifelse(dtp_42 >= 25, 25, dtp_42),
         censored_25 = ifelse(dtp_42 >= 25, "right", "none")) %>% 
  filter(!is.na(dtp))

save(df_analysis_mams,
     file = here("data", "cleaned-data", 
                 paste0(Sys.Date(), "_mams-clean.RData")))


#########################
# NC-002 (M-Pa-Z)
#########################

load("~/Documents/ucsf-research-git/tb-pacts-general/generated-datasets/2024-02-06_tb-pacts-data-with-ttp-measures.RData")
list2env(tb_pacts, envir = environment())

df_analysis_nc002 <- df_NC_002 %>% 
  filter(MBDY >= 0) %>% 
  dplyr::select(USUBJID, MBDY, MBSTRESN, MBSTRESC, ACTARM) %>% 
  distinct() %>% 
  # Taking the first 8 weeks of observation
  filter(MBDY <= 7*8 & MBDY >= 0,
         # There are 2 odd observations where MCSTRESC = POSITIVE or is missing with no corresponding MBSTRESN, we'll remove
         MBSTRESC != "POSITIVE",
         !is.na(MBSTRESC)) %>% 
  # Setting NA values to 42 for plotting purposes
  mutate(weeks_rounded = floor(MBDY/7),
         # Oddly, there are values of MBSTRESN > 42, we will reset these to 42
         dtp_42 = ifelse(is.na(MBSTRESN) | MBSTRESN > 42, 42, MBSTRESN),
         censored_42 = ifelse(MBSTRESC == ">42", "right", "none")) %>% 
  # Adding in 30 day censoring (for modeling)
  mutate(dtp_30 = ifelse(dtp_42 >= 30, 30, dtp_42),
         censored_30 = ifelse(dtp_42 >= 30, "right", "none")) %>% 
  # Adding in 25 day censoring (for modeling)
  mutate(dtp_25 = ifelse(dtp_42 >= 25, 25, dtp_42),
         censored_25 = ifelse(dtp_42 >= 25, "right", "none")) %>% 
  dplyr::select(-MBSTRESN, -MBSTRESC) %>% 
  mutate(weeks = MBDY/7)

save(df_analysis_nc002,
     file = here("data", "cleaned-data", 
                 paste0(Sys.Date(), "_NC-002-clean.RData")))

#########################
# NC-005 (M-Pa-Z-B)
#########################

df_analysis_nc005 <- df_NC_005 %>% 
  filter(MBDY >= 0) %>% 
  dplyr::select(patient.id = USUBJID, MBDY, MBSTRESN, MBSTRESC, ACTARM) %>% 
  distinct() %>% 
  # Taking the first 8 weeks of observation
  filter(MBDY <= 7*8 & MBDY >= 0) %>%
  mutate(weeks_rounded = floor(MBDY/7),
         dtp_42 = ifelse(MBSTRESC == "NEGATIVE" | MBSTRESN > 42, 42, MBSTRESN),
         censored_42 = ifelse(MBSTRESC == "NEGATIVE" | MBSTRESN > 42, "right", "none")) %>% 
  # Adding in 30 day censoring (for modeling)
  mutate(dtp_30 = ifelse(dtp_42 >= 30, 30, dtp_42),
         censored_30 = ifelse(dtp_42 >= 30, "right", "none")) %>% 
  # Adding in 25 day censoring (for modeling)
  mutate(dtp_25 = ifelse(dtp_42 >= 25, 25, dtp_42),
         censored_25 = ifelse(dtp_42 >= 25, "right", "none")) %>% 
  dplyr::select(-MBSTRESN, -MBSTRESC) %>% 
  mutate(weeks = MBDY/7)

save(df_analysis_nc005,
     file = here("data", "cleaned-data", 
                 paste0(Sys.Date(), "_NC-005-clean.RData")))

#########################
# Study 29 and 29X 
#########################

df_analysis_s29 <- df_TBTC_S29 %>% 
  filter(MBSTRESC == "NO GROWTH" | !is.na(MBSTRESN)) %>% 
  dplyr::select(patient.id = USUBJID, MBDY, MBSTRESN, MBSTRESC, MBTSTDTL, VISIT) %>% 
  distinct() %>% 
  # Taking the first 8 weeks of observation
  filter(MBDY < 7*9 & MBDY >= 0 | VISIT == "SCREENING") %>% 
  mutate(weeks_rounded = ifelse(MBDY > 0,floor(MBDY/7), 0),
         dtp_42 = ifelse(MBSTRESC == "NO GROWTH", 42, MBSTRESN),
         censored_42 = ifelse(MBSTRESC == "NO GROWTH", "right", "none")) %>% 
  # Adding in 30 day censoring (for modeling)
  mutate(dtp_30 = ifelse(dtp_42 >= 30, 30, dtp_42),
         censored_30 = ifelse(dtp_42 >= 30, "right", "none")) %>% 
  # Adding in 25 day censoring (for modeling)
  mutate(dtp_25 = ifelse(dtp_42 >= 25, 25, dtp_42),
         censored_25 = ifelse(dtp_42 >= 25, "right", "none")) %>% 
  dplyr::select(-MBSTRESN, -MBSTRESC) %>% 
  mutate(weeks = MBDY/7)

save(df_analysis_s29,
     file = here("data", "cleaned-data", 
                 paste0(Sys.Date(), "_TBTC-S29-clean.RData")))

#########################
# Study 29 and 29X 
#########################

df_analysis_s29x <- df_TBTC_S29x %>% 
  filter(MBSTRESC == "NO GROWTH" | !is.na(MBSTRESN)) %>% 
  dplyr::select(patient.id = USUBJID, MBDY, MBSTRESN, MBSTRESC, MBTSTDTL, VISIT) %>% 
  distinct() %>% 
  # Taking the first 8 weeks of observation
  filter(MBDY < 7*9 & MBDY >= 0 | VISIT == "SCREENING") %>% 
  mutate(weeks_rounded = ifelse(MBDY > 0,floor(MBDY/7), 0),
         dtp_42 = ifelse(MBSTRESC == "NO GROWTH", 42, MBSTRESN),
         censored_42 = ifelse(MBSTRESC == "NO GROWTH", "right", "none")) %>% 
  # Adding in 30 day censoring (for modeling)
  mutate(dtp_30 = ifelse(dtp_42 >= 30, 30, dtp_42),
         censored_30 = ifelse(dtp_42 >= 30, "right", "none")) %>% 
  # Adding in 25 day censoring (for modeling)
  mutate(dtp_25 = ifelse(dtp_42 >= 25, 25, dtp_42),
         censored_25 = ifelse(dtp_42 >= 25, "right", "none")) %>% 
  dplyr::select(-MBSTRESN, -MBSTRESC) %>% 
  mutate(weeks = MBDY/7)

save(df_analysis_s29,
     file = here("data", "cleaned-data", 
                 paste0(Sys.Date(), "_TBTC-S29x-clean.RData")))
