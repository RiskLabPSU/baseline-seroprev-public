if (!require("pacman")) install.packages("pacman")

pkgs <- c("brms",
          "here",
          "janitor",
          "readxl",
          "tidybayes",
          "tidyverse",
          "tidyterra")


pacman::p_load(pkgs, character.only = T)

hh_data <- read_rds(here("data", "hh_data.rds"))
i_data <- read_rds(here("data", "i_data.rds"))
serodata <- read_rds(here("data", "serology.rds"))

pop_estimates <- read_rds(here("project_wide_data", "village_pop_estimates.rds"))

village_state = tibble(abbreviation = c("zug", "dye", "iky", "oki", "oga", "ofo", "eze", "eny", "off"),
                       village = c("Zugu", "Dyegh", "Ikyogbakpev", "Okimbongha", "Ogamanna", "Ofonekom", "Ezeakataka", "Enyandulogu", "Offianka"),
                       lga = c("Vande Ikya", "Vande Ikya", "Vande Ikya", "Obubra", "Obubra", "Obubra", "Izzi", "Izzi", "Izzi"),
                       state = c("Benue", "Benue", "Benue", "Cross River", "Cross River", "Cross River", "Ebonyi", "Ebonyi", "Ebonyi"))

compare_groups_hierarchical <- function(data, outcome_var, group_var) {
  
  outcome_quo <- enquo(outcome_var)
  group_quo <- enquo(group_var)
  
  if (is.numeric(pull(data, !!outcome_quo)) && !is.factor(pull(data, !!outcome_quo))) {
    model_family <- gaussian()
  } else {
    outcome_vec <- pull(data, !!outcome_quo)
    n_unique_vals <- length(unique(na.omit(outcome_vec)))
    if (n_unique_vals <= 1) stop("Outcome variable has only one level.")
    model_family <- if (n_unique_vals == 2) bernoulli() else categorical()
  }
  
  formula <- bf(paste(quo_name(outcome_quo), "~ 1 + (1 |", quo_name(group_quo), ")"))
  
  message(paste("Fitting Hierarchical", model_family$family, "model for:", quo_name(outcome_quo)))
  
  fit <- brm(
    formula = formula,
    data = data,
    family = model_family,
    chains = 2, iter = 2000, warmup = 500,
    silent = 2, refresh = 0, seed = 123,
    control = list(adapt_delta = 0.95) # Higher adapt_delta is good for hierarchical models
  )
  
  # Find all random effect SD parameters for the grouping variable
  sd_param_names <- tidybayes::get_variables(fit) %>%
    stringr::str_subset(paste0("^sd_", quo_name(group_quo), "__"))
  
  # If no SD parameters are found (e.g., model failed), return FALSE
  if (length(sd_param_names) == 0) {
    warning("No random effect SD parameters found for this model.")
    return(list(evidence_of_difference = FALSE, model_fit = fit))
  }
  
  posterior_summary <- posterior_summary(fit, variable = sd_param_names)
  
  # Evidence of difference if the lower 95% CI for ANY SD parameter is > 0.1
  lower_ci_bounds <- posterior_summary[, "Q2.5"]
  evidence <- any(lower_ci_bounds > 0.1)
  
  return(list(
    evidence_of_difference = evidence, 
    model_summary = posterior_summary, # Return the full summary table
    model_fit = fit
  ))
}

# Table 1 analysis --------------------------------------------------------

tbl_1_data <- hh_data |>
  select(`_index`, household_id, n_people, compound, n_single_room, n_multi_room, rodent_removal_home_method, rodent_mitigation_method, rodent_remove_use, livestock, toilet, storage_cooked_containers, storage_packaged_containers, storage_seed_crop_containers_c, storage_seed_crop_containers_g, latitude, longitude, surroundings_household, surroundings_outside_other) |>
  mutate(rodents_enter = case_when(is.na(rodent_removal_home_method) ~ "No",
                                   TRUE ~ "Yes"),
         proximity_bush = case_when(str_detect(surroundings_household, "bush") ~ "Yes",
                                    TRUE ~ "No"),
         proximity_farm = case_when(str_detect(surroundings_household, "farm") ~ "Yes",
                                    TRUE ~ "No"),
         toilet_use = fct(case_when(str_detect(toilet, "toilet_plumbing") ~ "Plumbed toilet",
                                    str_detect(toilet, "pit_latrine") ~ "Pit latrine",
                                    str_detect(toilet, "open_system") ~ "Open system",
                                    str_detect(toilet, "field_defecation") ~ "Field defecation",
                                    TRUE ~ "Other"),
                          levels = c("Plumbed toilet", "Pit latrine", "Open system", "Field defecation", "Other")),
         rodent_control_animals = case_when(str_detect(rodent_mitigation_method, "cats|dogs") ~ "Yes",
                                            TRUE ~ "No"),
         rodent_control_containers = case_when(str_detect(rodent_mitigation_method, "containers") ~ "Yes",
                                               TRUE ~ "No"),
         rodent_control_burrows = case_when(str_detect(rodent_mitigation_method, "seal_burrows") ~ "Yes",
                                            TRUE ~ "No"),
         rodent_control_structure = case_when(str_detect(rodent_mitigation_method, "structure|convert") ~ "Yes",
                                              TRUE ~ "No"),
         rodent_remove_dispose = case_when(str_detect(rodent_remove_use, "dispose_them") ~ "Yes",
                                           TRUE ~ "No"),
         rodent_remove_eat = case_when(str_detect(rodent_remove_use, "eat_them") ~ "Yes",
                                       TRUE ~ "No"),
         rodent_remove_feed = case_when(str_detect(rodent_remove_use, "feed_them_to_animals") ~ "Yes",
                                        TRUE ~ "No"),
         rodent_remove_sell = case_when(str_detect(rodent_remove_use, "sell_them") ~ "Yes",
                                        TRUE ~ "No"),
         rodent_method_animal = case_when(str_detect(rodent_removal_home_method, "cat|dog") ~ "Yes",
                                          TRUE ~ "No"),
         rodent_method_trap = case_when(str_detect(rodent_removal_home_method, "gumtrap|traps") ~ "Yes",
                                        TRUE ~ "No"),
         rodent_method_poison = case_when(str_detect(rodent_removal_home_method, "poison") ~ "Yes",
                                          TRUE ~ "No"),
         rodent_method_sticks = case_when(str_detect(rodent_removal_home_method, "sticks") ~ "Yes",
                                          TRUE ~ "No")
  )  |>
  rowwise() |>
  mutate(n_buildings = sum(n_single_room, n_multi_room, na.rm = TRUE)) |>
  select(`_index`, household_id, n_people, n_buildings, n_single_room, n_multi_room, proximity_bush, proximity_farm, toilet_use, rodents_enter, rodent_control_animals, rodent_control_burrows, rodent_control_structure, rodent_control_containers,
         rodent_method_animal, rodent_method_poison, rodent_method_sticks, rodent_method_trap, rodent_remove_dispose, rodent_remove_eat, rodent_remove_feed, rodent_remove_sell) %>%
  left_join(serodata |>
              mutate(household_id = str_remove(as.character(id), "-\\d+$"),
                     is_positive = interpretation == "Positive") |>
              group_by(household_id) |>
              summarise(n_positive = sum(is_positive, na.rm = TRUE)),
            by = "household_id") |>
  mutate(abbreviation = str_split(household_id, "-", simplify = TRUE)[, 1]) |>
  left_join(village_state, by = "abbreviation") |>
  mutate(village = fct(village, levels = c(village_state$village)),
         n_single_room = as.integer(n_single_room),
         n_multi_room = as.integer(n_multi_room))


# Descriptive statistics

hh_size_comparison <- compare_groups_hierarchical(tbl_1_data, n_people, village)

n_buildings_comp <- compare_groups_hierarchical(tbl_1_data, n_buildings, village)

n_single_room_comp <- compare_groups_hierarchical(tbl_1_data, n_single_room, village)

proximity_bush_comp <- compare_groups_hierarchical(tbl_1_data, proximity_bush, village)

proximity_farm_comp <- compare_groups_hierarchical(tbl_1_data, proximity_farm, village)

toilet_use_comp <- compare_groups_hierarchical(tbl_1_data |>
                                              mutate(toilet_use_two_levels = case_when(toilet_use == "Field defecation" ~ "Field defecation",
                                                                                       TRUE ~ "Other"),
                                                     toilet_use_two_levels = fct(toilet_use_two_levels, levels = c("Field defecation", "Other"))), toilet_use_two_levels, village)

rodents_enter_comp <- compare_groups_hierarchical(tbl_1_data, rodents_enter, village)

removal_animal_comp <- compare_groups_hierarchical(tbl_1_data, rodent_method_animal, village)

removal_poison_comp <- compare_groups_hierarchical(tbl_1_data, rodent_method_poison, village)

removal_sticks_comp <- compare_groups_hierarchical(tbl_1_data, rodent_method_sticks, village)

removal_trap_comp <- compare_groups_hierarchical(tbl_1_data, rodent_method_trap, village)

use_dispose_comp <- compare_groups_hierarchical(tbl_1_data, rodent_remove_dispose, village)

use_eat_comp <- compare_groups_hierarchical(tbl_1_data, rodent_remove_eat, village)

use_feed_comp <- compare_groups_hierarchical(tbl_1_data, rodent_remove_feed, village)

#use_sell_comp <- compare_groups_hierarchical(tbl_1_data, rodent_remove_sell, village)

bayesian_results <- list(
  n_people = TRUE,
  n_buildings = TRUE,
  n_single_room = TRUE,
  proximity_bush = FALSE,
  proximity_farm = FALSE,
  toilet_use = FALSE,
  rodents_enter = FALSE,
  rodent_method_animal = FALSE,
  rodent_method_poison = FALSE,
  rodent_method_sticks = FALSE,
  rodent_method_trap = FALSE,
  rodent_remove_dispose = FALSE,
  rodent_remove_eat = FALSE,
  rodent_remove_feed = FALSE,
  rodent_remove_sell = FALSE
)

# Table 2 analysis --------------------------------------------------------

tbl_2_data <- i_data |>
  select(id, age, sex, ethnicity, religion, education, village_residence, months_residence_in_year, community_pob, income_individual, income, field_entry, field_work_freq, forest_entry, forest_entry_freq, current_rodent_consumption, past_rodent_consumption,
         rodent_eat_reason, excreta_cleaning, rodent_diseases) |>
  left_join(serodata,
            by = "id")
tbl_2_data$sex[tbl_2_data$id == "dye-6-1"] <- "Male"

income_dummies <- tbl_2_data |>
  select(id, income) |>
  separate_rows(income, sep = " ") |>
  mutate(income = str_trim(income)) |>
  dplyr::distinct(id, income) |>
  mutate(present = 1) |>
  pivot_wider(names_from = income, values_from = present, values_fill = 0) |>
  mutate(ag_worker = case_when(ag_worker == 0 & ag_worker_elsewhere == 1 ~ 1,
                               TRUE ~ ag_worker))|>
  select(-`NA`, -ag_worker_elsewhere)

rodent_eat_reason_dummies <- tbl_2_data |>
  select(id, rodent_eat_reason) |>
  separate_rows(rodent_eat_reason, sep = " ") |>
  mutate(rodent_eat_reason = str_trim(rodent_eat_reason)) |>
  dplyr::distinct(id, rodent_eat_reason) |>
  mutate(present = 1) |>
  pivot_wider(names_from = rodent_eat_reason, values_from = present, values_fill = 0) |>
  select(-`NA`)

tbl_2_data_processed <- tbl_2_data %>%
  left_join(income_dummies, by = "id") %>%
  left_join(rodent_eat_reason_dummies, by = "id") %>%
  mutate(village_residence = factor(village_residence, levels = c(village_state$village)),
         past_rodent_consumption = case_when(current_rodent_consumption == "Yes" ~ NA,
                                             TRUE ~ past_rodent_consumption),
         education = factor(education, levels = c("None/No formal schooling",
                                                  "Primary education (up until age 12)",
                                                  "Junior secondary education (up until age 16)",
                                                  "Senior secondary education (up until age 19)",
                                                  "University level education",
                                                  "Master's, Doctoral level degree",
                                                  "Other")))


age_comp <- compare_groups_hierarchical(tbl_2_data_processed, age, village_residence)

sex_comp <- compare_groups_hierarchical(tbl_2_data_processed, sex, village_residence)

education_comp <- compare_groups_hierarchical(tbl_2_data_processed |>
                                        mutate(education_four_levels = case_when(str_detect(education, "None") ~ "None",
                                                                                 str_detect(education, "Primary") ~ "Primary",
                                                                                 str_detect(education, "secondary") ~ "Secondary",
                                                                                 TRUE ~ "Other"),
                                               education_four_levels = factor(education_four_levels, levels = c("None", "Primary", "Secondary", "Other"))), education_four_levels, village_residence)

months_residence_comp <-  compare_groups_hierarchical(tbl_2_data_processed, months_residence_in_year, village_residence)

pob_comp <-  compare_groups_hierarchical(tbl_2_data_processed, community_pob, village_residence)

#income_comp <-  compare_groups_hierarchical(tbl_2_data_processed, income, village_residence)

field_entry_comp <-  compare_groups_hierarchical(tbl_2_data_processed, field_entry, village_residence)

forest_entry_comp <-  compare_groups_hierarchical(tbl_2_data_processed, forest_entry, village_residence)

current_rodent_consump_comp <-  compare_groups_hierarchical(tbl_2_data_processed, current_rodent_consumption, village_residence)

past_rodent_consump_comp <-  compare_groups_hierarchical(tbl_2_data_processed, past_rodent_consumption, village_residence)

excreta_comp <-  compare_groups_hierarchical(tbl_2_data_processed, excreta_cleaning, village_residence)

farming_comp <- compare_groups_hierarchical(tbl_2_data_processed, farming, village_residence)

ag_worker_comp <- compare_groups_hierarchical(tbl_2_data_processed, ag_worker, village_residence)

taste_comp <- compare_groups_hierarchical(tbl_2_data_processed, taste, village_residence)

availability_comp <- compare_groups_hierarchical(tbl_2_data_processed, availability, village_residence)

cheap_comp <- compare_groups_hierarchical(tbl_2_data_processed, cheap, village_residence)

nutrition_comp <- compare_groups_hierarchical(tbl_2_data_processed, nutrition, village_residence)

cultural_comp <- compare_groups_hierarchical(tbl_2_data_processed, cultural, village_residence)

# --- Store the final TRUE/FALSE results in a list ---
bayesian_results_individual <- list(
  # Demographics & Residency
  age = FALSE,
  sex = FALSE,
  education = TRUE,
  months_residence_in_year = FALSE,
  community_pob = TRUE,
  
  # Occupations
  farming = FALSE,
  ag_worker = FALSE,
  
  # Behaviours
  field_entry = FALSE,
  forest_entry = TRUE,
  current_rodent_consumption = TRUE,
  past_rodent_consumption = TRUE,
  excreta_cleaning = FALSE,
  
  # Reasons for Rodent Consumption
  taste = FALSE,
  availability = FALSE,
  cheap = FALSE,
  nutrition = FALSE,
  cultural = FALSE
)


# Supplementary -----------------------------------------------------------
aggregate_bayesian_comparison_results <- function(...) {
  results_list <- list(...)
  all_summaries <- list()
  
  for (name in names(results_list)) {
    result <- results_list[[name]]
    
    if (inherits(result, "list") && !is.null(result$model_summary)) {
      summary_df <- result$model_summary |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "variable") |>
        as_tibble() |>
        mutate(
          # 1. Map the internal variable names to the Table Labels
          Variable = case_match(
            name,
            "hh_size" ~ "Household size (individuals)",
            "n_buildings" ~ "Number of buildings",
            "n_single_room" ~ "Single-room buildings",
            "proximity_bush" ~ "Proximity to bush",
            "proximity_farm" ~ "Proximity to farm",
            "toilet_use" ~ "Type of toilet",
            "rodents_enter" ~ "Rodents enter home",
            "removal_animal" ~ "Rodent removal: Animals (Cat/Dog)",
            "removal_poison" ~ "Rodent removal: Poison",
            "removal_sticks" ~ "Rodent removal: Sticks",
            "removal_trap" ~ "Rodent removal: Traps",
            "use_dispose" ~ "Rodent use: Dispose",
            "use_eat" ~ "Rodent use: Eat",
            "use_feed" ~ "Rodent use: Feed to animals",
            "age" ~ "Age (years)",
            "sex" ~ "Sex",
            "education" ~ "Education",
            "months_residence" ~ "Months residence in year",
            "pob" ~ "Born in study village",
            "field_entry" ~ "Field entry",
            "forest_entry" ~ "Forest entry",
            "current_rodent_consump" ~ "Rodent consumption (Current)",
            "past_rodent_consump" ~ "Past rodent consumption only",
            "excreta" ~ "Cleaned rodent excreta",
            "farming" ~ "Occupation: Farming (own land)",
            "ag_worker" ~ "Occupation: Agricultural work (hired)",
            "taste" ~ "Reason: Taste",
            "availability" ~ "Reason: Availability",
            "cheap" ~ "Reason: Cheap",
            "nutrition" ~ "Reason: Nutrition",
            "cultural" ~ "Reason: Cultural",
            .default = name
          ),
          # 2. Clean up Parameter names (handling brms sd and mu prefixes)
          Parameter = variable |>
            stringr::str_remove("^sd_\\w+__") |>
            stringr::str_remove("^mu") |>
            stringr::str_replace("Other_Intercept", "Post-secondary (Intercept)") |>
            stringr::str_replace("_Intercept", " (Intercept)"),
          Estimate = round(Estimate, 3),
          `95% CrI` = paste0("[", round(Q2.5, 3), ", ", round(Q97.5, 3), "]"),
          Evidence_of_Difference = if_else(Q2.5 > 0.1, "Strong Evidence", "Weak/No Evidence")
        ) |>
        select(Variable, Parameter, Estimate, `95% CrI`, Evidence_of_Difference)
      
      all_summaries[[name]] <- summary_df
    }
  }
  
  bind_rows(all_summaries)
}

supplementary_table_comparison <- aggregate_bayesian_comparison_results(
  # Table 1 results
  hh_size = hh_size_comparison,
  n_buildings = n_buildings_comp,
  n_single_room = n_single_room_comp,
  proximity_bush = proximity_bush_comp,
  proximity_farm = proximity_farm_comp,
  toilet_use = toilet_use_comp,
  rodents_enter = rodents_enter_comp,
  removal_animal = removal_animal_comp,
  removal_poison = removal_poison_comp,
  removal_sticks = removal_sticks_comp,
  removal_trap = removal_trap_comp,
  use_dispose = use_dispose_comp,
  use_eat = use_eat_comp,
  use_feed = use_feed_comp,
  # Table 2 results
  age = age_comp,
  sex = sex_comp,
  education = education_comp,
  months_residence = months_residence_comp,
  pob = pob_comp,
  field_entry = field_entry_comp,
  forest_entry = forest_entry_comp,
  current_rodent_consump = current_rodent_consump_comp,
  past_rodent_consump = past_rodent_consump_comp,
  excreta = excreta_comp,
  farming = farming_comp,
  ag_worker = ag_worker_comp,
  taste = taste_comp,
  availability = availability_comp,
  cheap = cheap_comp,
  nutrition = nutrition_comp,
  cultural = cultural_comp
)

print(supplementary_table_comparison, n = Inf)

write_csv(supplementary_table_comparison, here("household_questionnaire", "output", "supplementary_table_bayesian_comparison.csv"))
