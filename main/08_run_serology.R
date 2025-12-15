# --- 1. SETUP ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(brms, tidyverse, here)

# Create an output directory if it doesn't exist
output_dir <- here("output", "models")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Load data

hh_data <- read_rds(here("data", "hh_data.rds"))
i_data <- read_rds(here("data", "i_data.rds"))
serodata <- read_rds(here("data", "serology.rds"))

pop_estimates <- read_rds(here("project_wide_data", "village_pop_estimates.rds"))

# Define village/state helper table
village_state <- tibble(
  abbreviation = c("zug", "dye", "iky", "oki", "oga", "ofo", "eze", "eny", "off"),
  village = c("Zugu", "Dyegh", "Ikyogbakpev", "Okimbongha", "Ogamanna", "Ofonekom", "Ezeakataka", "Enyandulogu", "Offianka"),
  state = c("Benue", "Benue", "Benue", "Cross River", "Cross River", "Cross River", "Ebonyi", "Ebonyi", "Ebonyi")
)


# --- 2. DATA PREPARATION ---
serodata <- i_data |>
  select(id, age, sex, village_residence) |>
  distinct() |>
  left_join(serodata |>
              mutate(abbreviation = str_split(id, pattern = "-", simplify = TRUE)[, 1],
                     y = case_when(interpretation == "Positive" ~ 1,
                                   TRUE ~ 0)),
            by = c("id")) |>
  mutate(
    y = if_else(interpretation == "Positive", 1, 0),
    village = factor(village_residence, levels = village_state$village),
    age = round(age, 0)
  ) |>
  select(y, age, village)

# --- 3. STATE AND VILLAGE PREVALENCE MODELS ---
# Aggregate data for binomial model
village_agg <- serodata |>
  group_by(village) |>
  summarise(positive = sum(y, na.rm = TRUE), total = n(), .groups = "drop")

state_agg <- serodata |>
  left_join(village_state |> select(village, state), by = "village") |>
  group_by(state) |>
  summarise(positive = sum(y, na.rm = TRUE), total = n(), .groups = "drop")

# Fit models
village_model <- brm(
  positive | trials(total) ~ 0 + village, data = village_agg,
  family = binomial(), prior = prior(normal(0, 5), class = "b"),
  chains = 4, cores = 4, iter = 4000, seed = 1111, control = list(adapt_delta = 0.95)
)

state_model <- brm(
  positive | trials(total) ~ 0 + state, data = state_agg,
  family = binomial(), prior = prior(normal(0, 5), class = "b"),
  chains = 4, cores = 4, iter = 4000, seed = 1111, control = list(adapt_delta = 0.95)
)

# Process and save results
village_prevalence_summary <- as_draws_df(village_model) |>
  select(starts_with("b_")) |>
  mutate(across(everything(), ~ plogis(.))) |> # Convert to probability scale
  pivot_longer(everything(), names_to = "village", values_to = "prevalence") |>
  mutate(village = str_remove(village, "b_village"))
saveRDS(village_prevalence_summary, file = here(output_dir, "village_prevalence_draws.rds"))


state_prevalence_summary <- as_draws_df(state_model) |>
  select(starts_with("b_")) |>
  mutate(across(everything(), ~ plogis(.))) |>
  pivot_longer(everything(), names_to = "state", values_to = "prevalence") |>
  mutate(state = str_remove(state, "b_state"))
saveRDS(state_prevalence_summary, file = here(output_dir, "state_prevalence_draws.rds"))


# --- 4. AGE-SEROPREVALENCE GAM MODEL ---
age_village_model <- brm(
  y ~ 0 + village + s(age, by = village),
  family = bernoulli(), data = serodata,
  prior = c(prior(normal(0, 5), class = "b"), prior(exponential(1), class = "sds")),
  chains = 4, cores = 4, iter = 4000, seed = 1111, control = list(adapt_delta = 0.95)
)

library(tidybayes)

newdata <- expand.grid(
  age = seq(0, max(serodata$age, na.rm = TRUE), by = 1),
  village = levels(serodata$village)
)

age_pred_draws <- newdata |>
  add_epred_draws(age_village_model)


# Save prediction data
saveRDS(age_pred_draws, file = here(output_dir, "age_seroprevalence_data.rds"))

# --- 5. OVERALL SEROPREVALENCE MODEL ---
overall_agg <- serodata |>
  summarise(
    positive = sum(y, na.rm = TRUE),
    total = n()
  )

overall_model <- brm(
  positive | trials(total) ~ 1,
  data = overall_agg,
  family = binomial(),
  chains = 4, cores = 4, seed = 1111
)

overall_prevalence_draws <- as_draws_df(overall_model) |>
  mutate(overall_prevalence = plogis(b_Intercept))
saveRDS(overall_prevalence_draws, file = here(output_dir, "overall_prevalence_draws.rds"))
