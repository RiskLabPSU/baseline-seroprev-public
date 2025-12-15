if (!require("pacman")) install.packages("pacman")

pkgs <- c("brms", "here", "tidyverse", "gt", "gtsummary")

pacman::p_load(pkgs, character.only = T)

hh_data <- read_rds(here("data", "hh_data.rds"))
i_data <- read_rds(here("data", "i_data.rds"))
serology_data <- read_rds(here("data", "serology.rds"))

village_state = tibble(abbreviation = c("zug", "dye", "iky", "oki", "oga", "ofo", "eze", "eny", "off"),
                       village = c("Zugu", "Dyegh", "Ikyogbakpev", "Okimbongha", "Ogamanna", "Ofonekom", "Ezeakataka", "Enyandulogu", "Offianka"),
                       lga = c("Vande Ikya", "Vande Ikya", "Vande Ikya", "Obubra", "Obubra", "Obubra", "Izzi", "Izzi", "Izzi"),
                       state = c("Benue", "Benue", "Benue", "Cross River", "Cross River", "Cross River", "Ebonyi", "Ebonyi", "Ebonyi"))

pop_estimates <- read_rds(here("project_wide_data", "village_pop_estimates.rds"))

# Now integrated directly in QMD