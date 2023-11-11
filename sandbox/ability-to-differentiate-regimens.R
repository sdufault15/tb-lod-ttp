library(tidyverse)
library(here)
library(latex2exp)
load(here("data", "model-generated", "2023-07-20_linear-mams-lod-42.RData"))
load(here("data", "model-generated", "2023-07-20_linear-mams-lod-30.RData"))
load(here("data", "model-generated", "2023-07-20_linear-remox-lod-42.RData"))
load(here("data", "model-generated", "2023-07-31_linear-remox-lod-30.RData"))
source(here("lib", "ucsf-color-palette.R"))

# Relative change in slope (Relative to HRZE)

theme_set(ggpubr::theme_pubr())

relative_function_mams <- function(model_obj, loq){
  df_coef <- coef(model_obj, summary = F)$Treatm_arm %>% 
    as_tibble()
  
  df_relative <- df_coef %>% 
    mutate(credible_est_no = row_number()) %>% 
    pivot_longer(cols = HR20ZM.Intercept:HRZQ.weeks,
                 names_to = "regimen",
                 values_to = "coef") %>% 
    filter(str_detect(regimen, "weeks")) %>%
    group_by(credible_est_no) %>% 
    mutate(rel_slope = coef/coef[regimen=="HRZE.weeks"]) 
  
  p_relative_slopes <- df_relative %>% 
    filter(regimen != "HRZE.weeks") %>% 
    mutate(regimen = str_remove(regimen, ".weeks")) %>% 
    ggplot(aes(x = rel_slope, fill = regimen)) +
    facet_wrap(~regimen, ncol = 1) + 
    geom_vline(xintercept = 1, 
               lty = 2) +
    geom_density(alpha = 0.3,
                 col = "white") +
    labs(x = TeX("$\\beta_k / \\beta_1, k \\neq 1$"),
         y = "Posterior Distribution",
         subtitle = paste0("LOQ = ",loq," days")) + 
    theme(legend.title = element_blank())
  
  p_percent_change_slopes <- df_relative %>% 
    filter(regimen != "HRZE.weeks") %>% 
    mutate(regimen = str_remove(regimen, ".weeks")) %>% 
    mutate(perc_change_slope = (rel_slope-1)*100) %>% 
    ggplot(aes(x = perc_change_slope, fill = regimen)) +
    facet_wrap(~regimen, ncol = 1) + 
    geom_vline(xintercept = 0, 
               lty = 2) +
    geom_density(alpha = 0.3,
                 col = "white") +
    labs(x = "% change in log10(TTP) slope",
         y = "Posterior Distribution",
         subtitle = paste0("LOQ = ",loq," days")) + 
    theme(legend.title = element_blank())
  
  df_posterior_distribution <- df_relative %>% 
    # Under the assumption that we're ONLY interested in LARGER slopes
    filter(regimen != "HRZE.weeks") %>% 
    mutate(regimen = str_remove(regimen, ".weeks")) %>% 
    mutate(perc_change_slope = (rel_slope-1)*100) %>% 
    group_by(regimen) %>% 
    summarise(p.change.geq.00 = mean(perc_change_slope >= 0),
              p.change.geq.05 = mean(perc_change_slope >= 5),
              p.change.geq.10 = mean(perc_change_slope >= 10),
              p.change.geq.15 = mean(perc_change_slope >= 15),
              p.change.geq.20 = mean(perc_change_slope >= 20)) %>% 
    pivot_longer(cols = p.change.geq.00:p.change.geq.20,
                 values_to = "posterior.prob",
                 names_to = "percent.change") %>% 
    mutate(`Percent Change` = case_when(str_detect(percent.change, ".00") ~ 0,
                                        str_detect(percent.change, ".05") ~ 5,
                                        str_detect(percent.change, ".10") ~ 10,
                                        str_detect(percent.change, ".15") ~ 15,
                                        str_detect(percent.change, ".20") ~ 20))
  
  p_posterior_tile <- df_posterior_distribution %>% 
    ggplot(aes(x = `Percent Change`, y = regimen, fill = posterior.prob)) +
    # geom_line(aes(group = regimen)) + 
    geom_tile() + 
    geom_text(aes(label = round(posterior.prob, 2)),
              col = "white") + 
    labs(subtitle = paste0("LOQ: ",loq," days")) +
    theme(legend.position = "none",
          axis.title.y = element_blank()) 
  
  output <- list(plots = list(p_relative_slopes = p_relative_slopes,
                              p_percent_change_slopes = p_percent_change_slopes,
                              p_posterior_tile = p_posterior_tile),
                 output = df_posterior_distribution)
  return(output)
}
relative_function_remox <- function(model_obj, loq){
  df_coef <- coef(model_obj, summary = F)$treat %>% 
    as_tibble()
  
  df_relative <- df_coef %>% 
    mutate(credible_est_no = row_number()) %>% 
    pivot_longer(cols = `1. 2EHRZ/4HR.Intercept`:`3. 2EMRZ/2MR.weeks`,
                 names_to = "regimen",
                 values_to = "coef") %>% 
    filter(str_detect(regimen, "weeks")) %>%
    group_by(credible_est_no) %>% 
    mutate(rel_slope = coef/coef[regimen=="1. 2EHRZ/4HR.weeks"]) 
  
  p_relative_slopes <- df_relative %>% 
    filter(!str_detect(regimen, "EHRZ")) %>% 
    mutate(regimen = str_remove(regimen, ".weeks")) %>% 
    ggplot(aes(x = rel_slope, fill = regimen)) +
    facet_wrap(~regimen, ncol = 1) + 
    geom_vline(xintercept = 1, 
               lty = 2) +
    geom_density(alpha = 0.3,
                 col = "white") +
    labs(x = TeX("$\\beta_k / \\beta_1, k \\neq 1$"),
         y = "Posterior Distribution",
         subtitle = paste0("LOQ = ",loq," days")) + 
    theme(legend.title = element_blank())
  
  p_percent_change_slopes <- df_relative %>% 
    filter(!str_detect(regimen, "EHRZ")) %>% 
    mutate(regimen = str_remove(regimen, ".weeks")) %>% 
    mutate(perc_change_slope = (rel_slope-1)*100) %>% 
    ggplot(aes(x = perc_change_slope, fill = regimen)) +
    facet_wrap(~regimen, ncol = 1) + 
    geom_vline(xintercept = 0, 
               lty = 2) +
    geom_density(alpha = 0.3,
                 col = "white") +
    labs(x = "% change in log10(TTP) slope",
         y = "Posterior Distribution",
         subtitle = paste0("LOQ = ",loq," days")) + 
    theme(legend.title = element_blank())
  
  df_posterior_distribution <- df_relative %>% 
    # Under the assumption that we're ONLY interested in LARGER slopes
    filter(!str_detect(regimen, "EHRZ")) %>% 
    mutate(regimen = str_remove(regimen, ".weeks")) %>% 
    mutate(perc_change_slope = (rel_slope-1)*100) %>% 
    group_by(regimen) %>% 
    summarise(p.change.geq.00 = mean(perc_change_slope >= 0),
              p.change.geq.05 = mean(perc_change_slope >= 5),
              p.change.geq.10 = mean(perc_change_slope >= 10),
              p.change.geq.15 = mean(perc_change_slope >= 15),
              p.change.geq.20 = mean(perc_change_slope >= 20)) %>% 
    pivot_longer(cols = p.change.geq.00:p.change.geq.20,
                 values_to = "posterior.prob",
                 names_to = "percent.change") %>% 
    mutate(`Percent Change` = case_when(str_detect(percent.change, ".00") ~ 0,
                                        str_detect(percent.change, ".05") ~ 5,
                                        str_detect(percent.change, ".10") ~ 10,
                                        str_detect(percent.change, ".15") ~ 15,
                                        str_detect(percent.change, ".20") ~ 20))
  
  p_posterior_tile <- df_posterior_distribution %>% 
    ggplot(aes(x = `Percent Change`, y = regimen, fill = posterior.prob)) +
    # geom_line(aes(group = regimen)) + 
    geom_tile() + 
    geom_text(aes(label = round(posterior.prob, 2)),
              col = "white") + 
    labs(subtitle = paste0("LOQ: ",loq," days")) +
    theme(legend.position = "none",
          axis.title.y = element_blank()) 
  
  output <- list(plots = list(p_relative_slopes = p_relative_slopes,
                              p_percent_change_slopes = p_percent_change_slopes,
                              p_posterior_tile = p_posterior_tile),
                 output = df_posterior_distribution)
  return(output)
}

df_42m <- relative_function_mams(m_linear_mams_42, loq = 42)
df_30m <- relative_function_mams(m_linear_mams_30, loq = 30)
df_42r <- relative_function_remox(m_linear_remox_42, loq = 42)
df_30r <- relative_function_remox(m_linear_remox_30, loq = 30)

library(patchwork)
df_42$plots$p_relative_slopes + df_30$plots$p_relative_slopes
df_42$plots$p_posterior_tile + df_30$plots$p_posterior_tile

bind_rows(mutate(df_42m$output, LOQ = "42 days"),
          mutate(df_30m$output, LOQ = "30 days")) %>% 
  ggplot(aes(x = `Percent Change`, y = posterior.prob, col = LOQ)) + 
  facet_wrap(~regimen) +
  geom_point() + 
  scale_color_manual(values = ucsf_col_pal[c(2,8)]) + 
  labs(y = TeX("Pr(\\theta \\geq x)"),
       x = "x, Relative change in log10(TTP) slope as compared to HRZE")

bind_rows(mutate(df_42r$output, LOQ = "42 days"),
          mutate(df_30r$output, LOQ = "30 days")) %>% 
  ggplot(aes(x = `Percent Change`, y = posterior.prob, col = LOQ)) + 
  facet_wrap(~regimen,
             ncol = 1) +
  geom_point() + 
  scale_color_manual(values = ucsf_col_pal[c(2,8)]) + 
  labs(y = TeX("Pr(\\theta \\geq x)"),
       x = "x, Relative change in log10(TTP) slope as compared to HRZE")
