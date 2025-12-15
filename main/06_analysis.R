if (!require("pacman")) install.packages("pacman")

pkgs <- c("brms", "here", "tidyverse")

pacman::p_load(pkgs, character.only = T)

hh_data <- read_rds(here("data", "hh_data.rds"))
i_data <- read_rds(here("data", "i_data.rds"))
serodata <- read_rds(here("data", "serology.rds"))

village_state = tibble(abbreviation = c("zug", "dye", "iky", "oki", "oga", "ofo", "eze", "eny", "off"),
                       village = c("Zugu", "Dyegh", "Ikyogbakpev", "Okimbongha", "Ogamanna", "Ofonekom", "Ezeakataka", "Enyandulogu", "Offianka"),
                       lga = c("Vande Ikya", "Vande Ikya", "Vande Ikya", "Obubra", "Obubra", "Obubra", "Izzi", "Izzi", "Izzi"),
                       state = c("Benue", "Benue", "Benue", "Cross River", "Cross River", "Cross River", "Ebonyi", "Ebonyi", "Ebonyi"))

pop_estimates <- read_rds(here("project_wide_data", "village_pop_estimates.rds"))

# Village level prevalence ------------------------------------------------

village_agg <- serodata |>
  mutate(abbreviation = str_split(id, pattern = "-", simplify = TRUE)[, 1],
         y = case_when(interpretation == "Positive" ~ 1,
                       TRUE ~ 0)) |>
  left_join(village_state, by = "abbreviation") |>
  group_by(village) |>
  summarise(positive = sum(y == 1),
            total = n(),
            .groups = "drop")

bayes_prev_model <- brm(
  positive | trials(total) ~ 0 + village,  # No global intercept, one parameter per village
  data = village_agg,
  family = binomial(link = "logit"),
  prior = prior(normal(0, 5), class = "b"),  # weakly informative prior
  chains = 4, cores = 4, iter = 4000, seed = 2025,
  control = list(adapt_delta = 0.95)
)

village_levels <- village_state$village

village_preds <- posterior_summary(bayes_prev_model, probs = c(0.025, 0.975)) |>
  as.data.frame() |>
  rownames_to_column(var = "village") |>
  filter(!village %in% c("lprior", "lp__")) |>  # remove unwanted rows
  mutate(
    village = str_remove_all(village, "b_village"),  # clean names
    village = factor(village, levels = village_levels),  # match factor levels
    prevalence = round(plogis(Estimate) * 100, 2),
    lower = round(plogis(Q2.5) * 100, 2),
    upper = round(plogis(Q97.5) * 100, 2),
    CrI = paste0(lower, "-", upper)
  ) |>
  arrange(village) |>
  left_join(serodata |>
              group_by(village) |>
              summarise(N = n(),
                        N_positive = sum(y == 1)),
            by = "village") |>
  select(Village = village, `N sampled` = N, `N positive` = N_positive, `Prevalence (%)` = prevalence, `Credible Interval (2.5% - 97.5%)` = CrI)


# State level prevalence --------------------------------------------------

state_agg <- serodata |>
  mutate(abbreviation = str_split(id, pattern = "-", simplify = TRUE)[, 1],
         y = case_when(interpretation == "Positive" ~ 1,
                       TRUE ~ 0)) |>
  left_join(village_state, by = "abbreviation") |>
  group_by(state) |>
  summarise(positive = sum(y == 1),
            total = n(),
            .groups = "drop") |>
  mutate(state = factor(state, levels = unique(village_state$state)))

bayes_prev_model <- brm(
  positive | trials(total) ~ 0 + state,  # No global intercept, one parameter per state
  data = state_agg,
  family = binomial(link = "logit"),
  prior = prior(normal(0, 5), class = "b"),  # weakly informative prior
  chains = 4, cores = 4, iter = 4000, seed = 2025,
  control = list(adapt_delta = 0.95)
)

state_preds <- posterior_summary(bayes_prev_model, probs = c(0.025, 0.975)) |>
  as.data.frame() |>
  rownames_to_column(var = "state") |>
  filter(!state %in% c("lprior", "lp__")) |>  # exclude log-posterior, etc.
  mutate(
    state = case_when(str_detect(state, "Benue") ~ "Benue",
                      str_detect(state, "Cross") ~ "Cross River",
                      str_detect(state, "Ebonyi") ~ "Ebonyi"),  # clean up names
    prevalence = round(plogis(Estimate) * 100, 2),
    lower = round(plogis(Q2.5) * 100, 2),
    upper = round(plogis(Q97.5) * 100, 2),
    CrI = paste0(lower, "-", upper)
  ) |>
  left_join(state_agg, by = "state") |>
  select(State = state, `N sampled` = total, `Prevalence (%)` = prevalence, `Credible Interval (2.5% - 97.5%)` = CrI)

formula <- bf(
  y ~ hh_size + age + time_in_village + sex + (1 | village/household_id),
  family = bernoulli(link = "logit")
)

priors <- c(
  set_prior("normal(0, 2.5)", class = "b"),           # For fixed effects
  set_prior("normal(0, 1)", class = "Intercept"),     # For the intercept
  set_prior("exponential(1)", class = "sd")           # For random effects (household)
)

fit <- brm(
  formula = formula,
  data = serodata,
  prior = priors,
  chains = 4,              # Number of MCMC chains
  iter = 4000,             # Total iterations per chain
  warmup = 1000,           # Number of warmup (burn-in) samples
  cores = 4,               # Number of CPU cores to use
  seed = 1234,             # Reproducibility
  control = list(adapt_delta = 0.95)  # Improves stability for hierarchical models
)

summary(fit)

plot(fit)

posterior_summary <- fixef(fit, probs = c(0.05, 0.95)) %>% 
  as.data.frame() %>%
  rownames_to_column("term")

ggplot(posterior_summary, aes(x = term, y = Estimate)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Q5, ymax = Q95), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Posterior Estimates for Fixed Effects",
       y = "Log-odds (with 90% credible interval)",
       x = "") +
  theme_minimal()

# Extract village-level random intercepts
village_ranef <- ranef(fit)$village[, , "Intercept"] %>%
  as.data.frame() %>%
  rownames_to_column("village") %>%
  rename(
    est = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )

# Plot
ggplot(village_ranef, aes(x = reorder(village, est), y = est)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(
    title = "Village-level random effects",
    x = "Village",
    y = "Log-odds of seropositivity (random effect)"
  ) +
  coord_flip() +
  theme_minimal()

# Create a dataset for marginal predictions by village
newdata <- serodata %>%
  group_by(village) %>%
  summarise(
    hh_size = mean(hh_size, na.rm = TRUE),
    age = mean(age, na.rm = TRUE),
    time_in_village = mean(time_in_village, na.rm = TRUE),
    sex = "Male"  # Can repeat for "Female" if needed
  )

# Predict probability with uncertainty (marginal effects)
preds <- posterior_epred(fit, newdata = newdata, re_formula = ~(1 | village))  # Only include village-level effects

# Summarise predictions
village_preds <- tibble(
  village = newdata$village,
  mean_prob = colMeans(preds),
  lower = apply(preds, 2, quantile, 0.025),
  upper = apply(preds, 2, quantile, 0.975)
) %>%
  mutate(village = fct(as.character(village), levels = village_state$village)) %>%
  arrange(village)

# Plot
ggplot(village_preds, aes(x = reorder(village, mean_prob), y = mean_prob)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(
    title = "Estimated prevalence of seropositivity by village",
    x = "Village",
    y = "Predicted prevalence"
  ) +
  coord_flip() +
  theme_minimal()

vc <- VarCorr(fit)
print(vc, digits = 2)

# ICCs (intra-class correlation coefficients)
village_sd <- vc$village$sd[1]
hh_sd <- vc$`village:household_id`$sd[1]
resid_var <- pi^2 / 3  # Logistic distribution residual variance

icc_village <- village_sd^2 / (village_sd^2 + hh_sd^2 + resid_var)
icc_household <- hh_sd^2 / (village_sd^2 + hh_sd^2 + resid_var)
