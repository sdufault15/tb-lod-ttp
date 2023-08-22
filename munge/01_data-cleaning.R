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
  filter(weeks <= 8, !DV %in% c("-99", ".")) %>% 
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
