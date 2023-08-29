# Example: https://discourse.mc-stan.org/t/non-linear-specification-variable-is-not-a-valid-distributional-or-non-linear-parameter/23521/9
Data <- expand.grid(r = factor(c(1:5)),
                    t = factor(c(1:8))
) %>%
  mutate(C = round(rnorm(nrow(.),500, 150)),
         D = round(rnorm(nrow(.),10000, 1000)),
         R = C/D,
         Srt = rnorm(nrow(.),0.5, 0.15))

bform <- bf(
  Srt ~ I,
  I ~ (1|r) + (1|t),
  nlf(I ~ C * B),
  nlf(B ~ a * R^b + c),
  a ~ 1,
  b ~ 1,
  c ~ 1,
  nl = TRUE
)

prior_list <- get_prior(bform,
                        data = Data,
                        family = Beta())


nl_priors <- c(
  set_prior("normal(0, 10)",
            class = "b",
            nlpar = "a"),
  set_prior("normal(0, 10)",
            class = "b",
            nlpar = "b"),
  set_prior("normal(0, 10)",
            class = "b",
            nlpar = "c")
)

NL_Mod <- brm(bform,
              Data,
              family = Beta(), 
              prior = nl_priors,
              inits = "random",
              iter = 2000,
              warmup = 1000,
              chains = 4,
              cores = ncores,
              backend = "cmdstan",
              normalize = FALSE,
              save_pars = save_pars(all = TRUE),
              control = list(adapt_delta = 0.9,
                             max_treedepth = 12)
)

posterior_epred(NL_Mod, nlpar = "I")
