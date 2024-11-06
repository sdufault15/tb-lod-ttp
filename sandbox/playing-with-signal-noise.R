library(tidyverse)
library(latex2exp)
library(here)
source(here("data", "cleaned-data", "2023-08-15_remoxtb-clean.RData"))
theme_set(ggpubr::theme_pubr())

df_replicates <- df_analysis_remox %>% 
  group_by(trial_no, weeks) %>% 
  mutate(replicate = row_number()) %>% 
  mutate(replicate = max(replicate)) %>% 
  filter(replicate > 1)

df_replicates %>% 
  dplyr::select(trial_no, weeks, treat, dtp, replicate) %>% 
  mutate(replicate = row_number()) %>% 
  pivot_wider(names_from = "replicate",
              values_from = "dtp") %>% 
  ggplot(aes(x = `1`)) +
  geom_abline(intercept = 0, slope = 1, col = 'red') +
  geom_point(aes(y = `2`),
             alpha = 0.8) +
  geom_point(aes(y = `3`),
             alpha = 0.8) +
  coord_equal() + 
  labs(x = TeX("DTP$_1$"), y = TeX("DTP$_k$, $k = \\{2,3\\}"))

df_replicates %>% 
  dplyr::select(trial_no, weeks, treat, dtp, replicate) %>% 
  mutate(replicate = row_number()) %>% 
  mutate(abs.diff = abs(dtp - dtp[replicate == 1])) %>% 
  filter(replicate != 1) %>% 
  ggplot(aes(x = dtp, y = abs.diff)) +
  geom_hline(yintercept = 0,
             col = "red") +
  geom_point() +
  coord_equal() + 
  labs(x = TeX("DTP$_1$"), y = TeX("|DTP$_k$ - DTP$_1$| \\ , when $k=\\{2,3\\}$"))

temp <- df_replicates %>% 
  dplyr::select(trial_no, weeks, treat, dtp, replicate) %>% 
  mutate(replicate = row_number()) %>% 
  pivot_wider(names_from = replicate,
              values_from = dtp) %>%
  mutate(dtp_cat = cut(`1`,breaks = c(0,15,20,25,30,42))) 

df_15 <- temp %>% 
  filter(dtp_cat == "(0,15]")
df_20 <- temp %>% 
  filter(dtp_cat == "(15,20]")
df_25 <- temp %>% 
  filter(dtp_cat == "(20,25]")
df_30 <- temp %>% 
  filter(dtp_cat == "(25,30]")
df_42 <- temp %>% 
  filter(dtp_cat == "(30,42]")

cor(x = df_15$`1`, 
    y = df_15$`2`,
    use = "complete.obs")

cor(x = df_20$`1`, 
    y = df_20$`2`,
    use = "complete.obs")

cor(x = df_25$`1`, 
    y = df_25$`2`,
    use = "complete.obs")

cor(x = df_30$`1`, 
    y = df_30$`2`,
    use = "complete.obs")

cor(x = df_42$`1`, 
    y = df_42$`2`,
    use = "pairwise.complete.obs")

m1 <- df_replicates %>% 
  dplyr::select(trial_no, weeks, treat, dtp, replicate) %>% 
  mutate(replicate = row_number()) %>% 
  pivot_wider(names_from = replicate,
              values_from = dtp) %>%
  lm(`1` ~ `2`,
         data = .)

preds <- predict(m1)

df_replicates_wide <- df_replicates %>% 
  dplyr::select(trial_no, weeks, treat, dtp, replicate) %>% 
  mutate(replicate = row_number()) %>% 
  pivot_wider(names_from = replicate,
              values_from = dtp)  

df_replicates_wide %>%   
  bind_cols(data.frame(pred.y = predict(m1, newdata = df_replicates_wide))) %>% 
  mutate(pred.error = pred.y - `1`) %>% 
  ggplot(aes(x = `1`, y = pred.error)) + 
  geom_point()


df_analysis_remox %>% 
  dplyr::select(trial_no, weeks, treat, dtp) %>% 
  arrange(trial_no, weeks) %>% 
  group_by(trial_no) %>% 
  slice_head(n= 5) %>% 
  ggplot(aes(x = weeks, y = dtp)) + 
  facet_wrap(~trial_no) +
  geom_point()
