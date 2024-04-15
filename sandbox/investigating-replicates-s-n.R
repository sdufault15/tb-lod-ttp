library(olsrr)
library(tidyverse)
people_with_replicates <- df_analysis_remox %>% 
  group_by(trial_no, weeks) %>% 
  summarise(n = n_distinct(dtp)) %>% 
  filter(n > 1)

df_subset <- df_analysis_remox %>% 
  left_join(people_with_replicates) %>% 
  filter(!is.na(n))

head(df_subset)

df_subset <- df_subset %>% 
  group_by(trial_no, weeks) %>% 
  mutate(replicate = row_number()) %>% 
  dplyr::select(trial_no, weeks, dtp, replicate) %>% 
  pivot_wider(values_from = dtp,
              names_from = replicate) 

df_subset %>% 
  ggplot(aes(x = `1`, y = `2`)) + 
  coord_fixed(xlim = c(0,42),
              ylim = c(0,42)) + 
  geom_abline(slope = 1, 
              intercept = 0,
              col = "red") +
  geom_point(alpha = 0.3) + 
  labs(x = "DTP measure 1",
       y = "DTP measure 2")

model <- lm(`1` ~ `2`, data = df_subset)
ols_test_breusch_pagan(model)  

df_subset %>% 
  ggplot(aes(x = `1`, y = `2`-`1`)) +
  geom_hline(yintercept = 0,
             col = "red") + 
  geom_point(alpha = 0.3)
  