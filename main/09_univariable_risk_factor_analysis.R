# --- 09_univariable_risk_factor_analysis.R ---
# This script runs all univariable risk factor models and saves the outputs.

# --- 1. SETUP ---
message("Loading packages and data...")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(brms, tidyverse, here, tidybayes, ggridges)

# Create an output directory
output_dir <- here("manuscript", "cross_sectional_study_outputs", "models")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


run_univariable_models <- function(data, outcome, vars, group_var = NULL,
                                   iter = 2000, chains = 4, cores = 4, seed = 1111,
                                   return_models = TRUE) {
  fit_one_model <- function(var) {
    rhs <- if (!is.null(group_var)) {
      paste0(var, " + (1 | ", group_var, ")")
    } else {
      var
    }
    
    formula <- as.formula(paste0(outcome, " ~ ", rhs))
    
    fit <- brm(
      formula = formula,
      data = data,
      family = bernoulli(),
      iter = iter,
      chains = chains,
      cores = cores,
      seed = seed,
      silent = TRUE,
      refresh = 0
    )
    
    return(fit)
  }
  
  models <- set_names(map(vars, fit_one_model), vars)
  
  if (return_models) {
    return(models)
  } else {
    # Extract summary statistics only
    tidy_summary <- map_dfr(models, function(fit) {
      tidy(fit, effects = "fixed", conf.level = 0.95) %>%
        filter(term != "(Intercept)") %>%
        mutate(
          OR = exp(estimate),
          lower = exp(conf.low),
          upper = exp(conf.high),
          excludes_1 = lower > 1 | upper < 1
        )
    }, .id = "predictor") %>%
      select(predictor, OR, lower, upper, excludes_1)
    
    return(tidy_summary)
  }
}

plot_posterior_ORs <- function(model_list, return_table = TRUE, max_log10_width = 6) {
  
  # 1. Extract and summarise posterior draws
  posterior_draws <- map_dfr(names(model_list), function(var) {
    model <- model_list[[var]]
    
    model |>
      gather_draws(`b_.*`, regex = TRUE) |>
      filter(.variable != "b_Intercept") |>
      mutate(
        predictor = var,
        term = str_remove(.variable, "^b_"),
        OR = exp(.value)
      )
  })
  
  summary_table <- posterior_draws |>
    group_by(predictor, term) |>
    summarise(
      median_OR = median(OR, na.rm = TRUE),
      lower = quantile(OR, 0.025, na.rm = TRUE),
      upper = quantile(OR, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 2. Add N and N_positive to summary_table
  Ns_table <- map_dfr(names(model_list), function(var) {
    model <- model_list[[var]]
    data <- model$data
    vars <- all.vars(model$formula$formula)
    outcome <- vars[1]
    predictor_var <- vars[2]
    
    if (is.factor(data[[predictor_var]]) || is.character(data[[predictor_var]])) {
      levels_vec <- levels(factor(data[[predictor_var]]))
      
      # Skip first level (reference), handled separately
      map_dfr(levels_vec[-1], function(level) {
        idx <- which(data[[predictor_var]] == level)
        tibble(
          predictor = var,
          term = paste0(predictor_var, level),
          N = length(idx),
          N_positive = sum(data[[outcome]][idx] == 1, na.rm = TRUE)
        )
      })
    } else {
      # For continuous predictors like age
      tibble(
        predictor = var,
        term = predictor_var,
        N = sum(!is.na(data[[predictor_var]])),
        N_positive = sum(data[[outcome]][!is.na(data[[predictor_var]])] == 1, na.rm = TRUE)
      )
    }
  })
  
  summary_table <- left_join(summary_table, Ns_table, by = c("predictor", "term"))
  
  # 3. Add reference level OR = 1 rows
  ref_rows <- map_dfr(names(model_list), function(var) {
    model <- model_list[[var]]
    data <- model$data
    
    vars <- all.vars(model$formula$formula)
    outcome <- vars[1]
    predictor_var <- vars[2]
    
    if (is.factor(data[[predictor_var]]) || is.character(data[[predictor_var]])) {
      ref_level <- levels(factor(data[[predictor_var]]))[1]  # assumes first level is reference
      idx <- which(data[[predictor_var]] == ref_level)
      tibble(
        predictor = var,
        term = paste0(predictor_var, ref_level),  # match brms naming like sexMale
        median_OR = 1,
        lower = 1,
        upper = 1,
        N = length(idx),
        N_positive = sum(data[[outcome]][idx] == 1, na.rm = TRUE)
      )
    } else {
      NULL  # No ref row for numeric predictors
    }
  })
  
  summary_table <- bind_rows(summary_table, ref_rows)
  
  # 4. Clean up and label terms for plotting
  valid_terms <- summary_table$term
  
  posterior_draws <- posterior_draws |>
    filter(term %in% valid_terms) |>
    mutate(
      label = if_else(
        str_starts(term, predictor),
        str_remove(term, paste0("^", predictor)),
        term
      ),
      label = if_else(label == "", predictor, paste(predictor, label, sep = ": "))
    )
  
  label_order <- posterior_draws |>
    group_by(label) |>
    summarise(median_OR = median(OR), .groups = "drop") |>
    arrange(median_OR) |>
    pull(label)
  
  posterior_draws <- posterior_draws |>
    mutate(label = factor(label, levels = label_order))
  
  # 5. Plotting
  or_plot <- posterior_draws |>
    ggplot(aes(x = OR, y = label, fill = label)) +
    geom_density_ridges(
      scale = 1.2,
      rel_min_height = 0.01,
      alpha = 0.6,
      color = NA,
      na.rm = TRUE
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_x_log10(limits = c(10^-max_log10_width, 10^max_log10_width)) +
    scale_fill_viridis_d(option = "C", guide = "none") +
    labs(
      x = "Odds Ratio (log scale)",
      y = NULL,
      title = "Posterior Distributions of Odds Ratios (Univariable Models)"
    ) +
    theme_minimal(base_size = 14)
  
  # 6. Return output
  if (return_table) {
    return(list(plot = or_plot, summary_table = summary_table))
  } else {
    return(or_plot)
  }
}


i_data_models <- tbl_2_data_processed |> 
  mutate(seropositive = case_when(interpretation == "Positive" ~ TRUE,
                                  TRUE ~ FALSE),
         household_id = str_remove(id, "-\\d+$"),
         # Collapse non-christian religions
         religion = fct_collapse(religion,
                                 "Christian" = "Christian",
                                 "Non-christian" = c("Muslim", "Traditionalist", "Other")),
         # Collapse educational levels
         education = case_when(str_detect(education, "None") ~ "None",
                                                    str_detect(education, "Primary") ~ "Primary",
                                                    str_detect(education, "secondary") ~ "Secondary",
                                                    TRUE ~ "Other"),
         education = factor(education, levels = c("None", "Primary", "Secondary", "Other")))

hh_data_models <- tbl_1_data |>
  mutate(positive_hh = case_when(n_positive == 0 ~ FALSE,
                                 n_positive >= 1 ~ TRUE))

# Only need to run if not previously produced
if(!file.exists(here("manuscript", "cross_sectional_study_outputs", "models", "demographic_models_results.rds"))) {
  # Demographic models - individual level
  # age
  # sex
  # education - grouped as none, primary, junior secondary, senior secondary and above
  # religion
  
  demographic_models_uni <- run_univariable_models(
    data = i_data_models,
    outcome = "seropositive",
    vars = c("age", "sex", "education", "religion"),
    group_var = "village_residence",
    return_models = TRUE
  )
  
  demographic_models_results <- plot_posterior_ORs(demographic_models_uni)
  
  write_rds(demographic_models_results, here("household_questionnaire", "manuscript", "cross_sectional_study_outputs", "models", "demographic_models_results.rds"))
  
  # Environment models - hh level
  # number of people in household
  # number of buildings
  # n of single room buildings
  # n of multi room buildings
  # bush in the proximity of household buildings
  # farms in the proximity of household buildings
  # toilet use - field defecation vs. other
  # rodents enter
  
  environmental_hh_models_uni <- run_univariable_models(
    data = hh_data_models,
    outcome = "positive_hh",
    vars = c("n_people", "n_buildings", "n_single_room", "n_multi_room", "proximity_bush", "proximity_farm", "toilet_use", "rodents_enter"),
    group_var = "village",
    return_models = TRUE
  )
  
  environmental_hh_models_results <- plot_posterior_ORs(environmental_hh_models_uni)
  
  write_rds(environmental_hh_models_results, here("household_questionnaire", "manuscript", "cross_sectional_study_outputs", "models", "environmental_hh_results.rds"))
  
  # Environment models - individual level
  # forest entry
  # field entry
  # rodent excreta cleaning
  
  environmental_models_uni <- run_univariable_models(
    data = i_data_models,
    outcome = "seropositive",
    vars = c("forest_entry", "field_entry", "excreta_cleaning"),
    group_var = "village_residence",
    return_models = TRUE
  )
  
  environmental_models_results <- plot_posterior_ORs(environmental_models_uni)
  
  write_rds(environmental_models_results, here("household_questionnaire", "manuscript", "cross_sectional_study_outputs", "models", "environmental_uni_results.rds"))
  
  # Behavioural models - household level
  # rodent control - animals, burrows, structural changes, conact
  # rodent removal method - animals, poison, sticks, traps
  # rodent carcass disposal - dispose, eat, feed to animals, sell
  
  behavioural_hh_models_uni <- run_univariable_models(
    data = hh_data_models,
    outcome = "positive_hh",
    vars = c("rodent_control_animals", "rodent_control_burrows", "rodent_control_structure", "rodent_control_containers", "rodent_method_animal", "rodent_method_poison", "rodent_method_sticks", "rodent_method_trap",
             "rodent_remove_dispose", "rodent_remove_eat", "rodent_remove_feed", "rodent_remove_sell"),
    group_var = "village",
    return_models = TRUE)
  
  behavioural_hh_models_results <- plot_posterior_ORs(behavioural_hh_models_uni)
  
  write_rds(behavioural_hh_models_results, here("household_questionnaire", "manuscript", "cross_sectional_study_outputs", "models", "behavioural_hh_results.rds"))
  
  # Behavioural models - individual
  # rodent consumption - current consumption, past consumption
  # rodent eat reason - taste, availability, cheap, nutrition, cultural, hunger, other.y
  
  behavioural_models_uni <- run_univariable_models(
    data = i_data_models,
    outcome = "seropositive",
    vars = c("current_rodent_consumption", "past_rodent_consumption", "taste", "availability", "cheap", "nutrition", "cultural", "hunger", "other.y"),
    group_var = "village_residence",
    return_models = TRUE
  )
  
  behavioural_models_results <- plot_posterior_ORs(behavioural_models_uni)
  
  write_rds(behavioural_models_results, here("household_questionnaire", "manuscript", "cross_sectional_study_outputs", "models", "behavioural_uni_results.rds"))
  
  # Occupational models - individual
  # months resident
  # community pob
  # occupation - farming, trader, driver, teacher, student, animal, hunter_trapper, artisan, ag_worker, government_worker, pensioner, clergy, ntfps, fishing, timer, other.x
  
  occupational_models_uni <- run_univariable_models(
    data = i_data_models,
    outcome = "seropositive",
    vars = c("months_residence_in_year", "community_pob", "farming", "ag_worker", "trader", "driver", "teacher", "student", "animal", "hunter_trapper", "artisan", "government_worker", "pensioner", "clergy", "ntfps",
             "fishing", "timber", "other.x"),
    group_var = "village_residence",
    return_models = TRUE
  )
  
  occupational_models_results <- plot_posterior_ORs(occupational_models_uni)
  
  write_rds(occupational_models_results, here("household_questionnaire", "manuscript", "cross_sectional_study_outputs", "models", "occupational_uni_results.rds"))
  
} else {
  
  demographic_models_results <- read_rds(here("household_questionnaire", "manuscript", "cross_sectional_study_outputs", "models", "demographic_models_results.rds"))
  environmental_hh_models_results <- read_rds(here("household_questionnaire", "manuscript", "cross_sectional_study_outputs", "models", "environmental_hh_results.rds"))
  environmental_models_results <- read_rds(here("household_questionnaire", "manuscript", "cross_sectional_study_outputs", "models", "environmental_uni_results.rds"))
  behavioural_hh_models_results <- read_rds(here("household_questionnaire", "manuscript", "cross_sectional_study_outputs", "models", "behavioural_hh_results.rds"))
  behavioural_models_results <- read_rds(here("household_questionnaire", "manuscript", "cross_sectional_study_outputs", "models", "behavioural_uni_results.rds"))
  occupational_models_results <- read_rds(here("household_questionnaire", "manuscript", "cross_sectional_study_outputs", "models", "occupational_uni_results.rds"))
  
}