if (!require("pacman")) install.packages("pacman")

pkgs <- c("cowplot",
          "crayon",
          "crul",
          "dm",
          "DT",
          "flextable",
          "forcats",
          "geodata",
          "ggrepel",
          "glue",
          "googledrive",
          "gtsummary",
          "haven",
          "here",
          "httr",
          "janitor",
          "kableExtra",
          "knitr",
          "Microsoft365R",
          "progress",
          "RcppSimdJson",
          "readxl",
          "terra",
          "tidyverse",
          "tidyterra",
          "viridisLite",
          "widyr",
          "zip")

pacman::p_load(pkgs, character.only = T)

project_CRS = "EPSG:4326"
utm_nigeria_CRS = "EPSG:26332"

village_state = tibble(abbreviation = c("zug", "dye", "iky", "oki", "oga", "ofo", "eze", "eny", "off"),
                       village = c("Zugu", "Dyegh", "Ikyogbakpev", "Okimbongha", "Ogamanna", "Ofonekom", "Ezeakataka", "Enyandulogu", "Offianka"),
                       lga = c("Vande Ikya", "Vande Ikya", "Vande Ikya", "Obubra", "Obubra", "Obubra", "Izzi", "Izzi", "Izzi"),
                       state = c("Benue", "Benue", "Benue", "Cross River", "Cross River", "Cross River", "Ebonyi", "Ebonyi", "Ebonyi"))

# Number of states and villages per state
num_states <- 3
villages_per_state <- 3

# Create a vector of distinct states
states <- rep(1:num_states, each = villages_per_state)

# Define the colour palette for states
state_colours <- viridis(num_states)

# Create a function to generate perceptually distinct colours for villages within each state
generate_village_colours <- function(num_colours) {
  village_colours <- viridis(num_colours)
  return(village_colours)
}

# Create a colour table with distinct but related colours for villages within each state
colour_table <- tibble(
  state_colour = rep(state_colours, each = villages_per_state),
  village_colour = generate_village_colours(num_states * villages_per_state)
)

village_state <- tibble(village_state,
                        village_colour = colour_table$village_colour,
                        state_colour = colour_table$state_colour)

dir.create(here("project_wide_data"))
write_rds(village_state, here("project_wide_data", "village_state.rds"))

village_colours <- setNames(village_state$village_colour, village_state$village)
write_rds(village_colours, here("project_wide_data", "village_colours.rds"))

occupation_matching <- data.frame(
  "Individual questionnaire" = c("Farming", "Assist with agricultural work (in household fields)", "Assist with agricultural work (in other households fields)", 
                                 "Hunter/Trapper", "Fishing", "Timber", "Collect forest goods (NTFPs)", 
                                 "Animal husbandry", "Trader", "Artisan/Handiwork/Carpenter", "Driver", 
                                 "Teacher", "Clergy (minister, pastor)", "Student", "Government worker", 
                                 "Pensioner", "Other"),
  "SCAPES code" = c("farming", "ag_worker", "ag_worker_elsewhere", 
                    "hunter_trapper", "fishing", "timber", "ntfps", "animal",
                    "trader", "artisan", "driver", "teacher", "clergy", "student",
                    "government_worker", "pensioner", "other"),
  "DHS 2018" = c("Agriculture", "Agriculture",
                 "Agriculture", "Agriculture",
                 "Agriculture", "Agriculture",
                 "Agriculture", "Agriculture",
                 "Sales", "Skilled manual", "Other", "Professional/technical/managerial", 
                 "Professional/technical/managerial", "Not working/NA/Did not work", 
                 "Clerical/Professional/technical/managerial", "Not working/NA/Did not work", "Unskilled manual/Services"),
  "ISCO-08" = c("Field Crop and Vegetable Growers", "Subsistence Crop Farmers", "Subsistence Crop Farmers", 
                "Subsistence Farmers, Fishers, Hunters and Gatherers", "Subsistence Farmers, Fishers, Hunters and Gatherers", 
                "Subsistence Farmers, Fishers, Hunters and Gatherers", "Subsistence Fishers, Hunters, Trappers and Gatherers", 
                "Subsistence Livestock Farmers", "Street and Related Sales and Service Workers", 
                "Craft and Related Trades Workers", "Car, Van and Motorcycle Drivers", 
                "Teaching Professionals", "Social and Religious Professionals", "N/A", 
                "Business and Administration Associate Professionals", "N/A", "N/A"),
  "ISCO-08 code" = c(6111, 631, 631, 6340, 6340, 6340, 6340, 6320, 95, 7, 832, 23, 263, NA, 33, NA, NA)
)

nigeria_states <- c("Abia", "Adamawa", "Akwa Ibom", "Anambra", "Bauchi", "Bayelsa", "Benue", "Borno", 
                    "Cross-River", "Delta", "Ebonyi", "Edo", "Ekiti", "Enugu", "Gombe", "Imo", "Jigawa", 
                    "Kaduna", "Kano", "Katsina", "Kebbi", "Kogi", "Kwara", "Lagos", "Nasarawa", "Niger", 
                    "Ogun", "Ondo", "Osun", "Oyo", "Plateau", "Rivers", "Sokoto", "Taraba", "Yobe", "Zamfara",
                    "Federal Capital Territory")

Nigeria <- gadm(country = "NGA", level = 0, path = here("project_wide_data", "spatial"))
lga <- gadm(country = "NGA", level = 2, path = here("project_wide_data", "spatial"))
