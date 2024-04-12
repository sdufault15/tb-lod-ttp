# REPLICATE CLEANING

#########################
# MAMS-TB
#########################
df_mams_replicates <- df_analysis_mams %>% 
  dplyr::select(arm, patient.id, weeks, dtp_42) %>% 
  group_by(patient.id, weeks) %>% 
  mutate(replicate = row_number()) %>% 
  arrange(patient.id, weeks) %>% 
  pivot_wider(values_from = dtp_42,
              names_from = replicate) %>%
  # Keep replicates only (no singlets)
  filter(!is.na(`1`) & !is.na(`2`)) %>% 
  rename(rep1 = `1`, rep2 = `2`)

#########################
# TB-PACTS: NC_002
#########################
load(here("..", "tb-pacts-general", "generated-datasets","2024-03-26_tb-pacts-data-with-ttp-measures.RData"))
list2env(tb_pacts, environment())
rm(tb_pacts)

df_NC_002_replicates <- df_NC_002 %>% 
  # Smaller dataset to work with
  dplyr::select(USUBJID, MBDY, MBSTRESN, ACTARM) %>% 
  distinct() %>% 
  # Taking the first 8 weeks of observation
  filter(MBDY <= 7*8 & MBDY >= 0) %>% 
  # Setting NA values to 42 for plotting purposes
  mutate(dtp_42 = ifelse(is.na(MBSTRESN), 42, MBSTRESN)) %>% #,
  # censored_42 = ifelse(is.na(MBSTRESN), "right", "none")) %>% 
  dplyr::select(-MBSTRESN) %>% 
  # Identify replicates
  group_by(USUBJID, MBDY) %>% 
  mutate(replicate = row_number()) %>% 
  pivot_wider(values_from = dtp_42,
              names_from = replicate) %>%
  # Keep replicates only (no singlets)
  filter(!is.na(`1`) & !is.na(`2`)) %>% 
  rename(rep1 = `1`, rep2 = `2`)

df_replicates <- list(df_mams_replicates = df_mams_replicates,
                      df_NC_002_replicates = df_NC_002_replicates)
save(df_replicates,
     file = here("data", "cleaned-data", paste0(Sys.Date(), "_replicate-datasets.RData")))

#########################
# TB-PACTS: NC_005
#########################
# Doesn't appear to have duplicates (except for individual 10101 at week 8 and 15)
df_NC_005_replicates <- df_NC_005 %>% 
  # Smaller dataset to work with
  dplyr::select(USUBJID, MBDY, MBSTRESN, MBSTRESC, ACTARM) %>% 
  distinct() %>% 
  # Taking the first 8 weeks of observation
  filter(MBDY <= 7*8 & MBDY >= 0) %>% 
  # Setting NA values to 42 for plotting purposes
  mutate(dtp_42 = ifelse(is.na(MBSTRESN) & MBSTRESC == "NEGATIVE", 42, MBSTRESN)) %>% #arrange(USUBJID, MBDY)#,
  # censored_42 = ifelse(is.na(MBSTRESN), "right", "none")) %>% 
  dplyr::select(-MBSTRESN, -MBSTRESC) %>% 
  # Identify replicates
  group_by(USUBJID, MBDY) %>% 
  mutate(replicate = row_number()) %>% 
  pivot_wider(values_from = dtp_42,
              names_from = replicate) %>% 
  filter(!is.na(`1`) & !is.na(`2`)) %>% 
  arrange(USUBJID, MBDY)

#########################
# TB-PACTS: TBTC Study 29(x)
#########################
# Also doesn't seem to have replicates (except patient 0096 and 0444 at week 8)
df_TBTC_S29 %>% 
  # Smaller dataset to work with
  dplyr::select(USUBJID, MBDY, MBSTRESN, MBSTRESC, ACTARM) %>% 
  distinct() %>% 
  # Taking the first 8 weeks of observation
  filter(MBDY <= 7*8 & MBDY >= 0) %>% 
  # Setting NA values to 42 for plotting purposes
  mutate(dtp_42 = ifelse(is.na(MBSTRESN) & MBSTRESC == "NO GROWTH", 42, MBSTRESN)) %>% 
  dplyr::select(-MBSTRESN, -MBSTRESC) %>% 
  # Identify replicates
  group_by(USUBJID, MBDY) %>% 
  mutate(replicate = row_number()) %>% 
  pivot_wider(values_from = dtp_42,
              names_from = replicate) %>% 
  filter(!is.na(`1`) & !is.na(`2`)) %>% arrange(USUBJID, MBDY)

df_TBTC_S29x %>% 
  # Also only 2 individuals 
  dplyr::select(USUBJID, MBDY, MBSTRESN, MBSTRESC, ACTARM) %>% 
  distinct() %>% 
  # Taking the first 8 weeks of observation
  filter(MBDY <= 7*8 & MBDY >= 0) %>% 
  # Setting NA values to 42 for plotting purposes
  mutate(dtp_42 = ifelse(is.na(MBSTRESN) & MBSTRESC == "NO GROWTH", 42, MBSTRESN)) %>% 
  dplyr::select(-MBSTRESN, -MBSTRESC) %>% 
  # Identify replicates
  group_by(USUBJID, MBDY) %>% 
  mutate(replicate = row_number()) %>% 
  pivot_wider(values_from = dtp_42,
              names_from = replicate) %>% 
  filter(!is.na(`1`) & !is.na(`2`)) %>% arrange(USUBJID, MBDY)
