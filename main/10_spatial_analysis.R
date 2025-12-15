# --- 1. SETUP ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(sf, spdep, tidyverse, here, broom, cowplot)

# Create an output directory
output_dir <- here("output", "spatial_analysis")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load data
hh_data <- read_rds(here("data", "hh_data.rds"))
i_data <- read_rds(here("data", "i_data.rds"))
serodata <- read_rds(here("data", "serology.rds"))

pop_estimates <- read_rds(here("project_wide_data", "village_pop_estimates.rds"))
village_state = tibble(abbreviation = c("zug", "dye", "iky", "oki", "oga", "ofo", "eze", "eny", "off"),
                       village = c("Zugu", "Dyegh", "Ikyogbakpev", "Okimbongha", "Ogamanna", "Ofonekom", "Ezeakataka", "Enyandulogu", "Offianka"),
                       lga = c("Vande Ikya", "Vande Ikya", "Vande Ikya", "Obubra", "Obubra", "Obubra", "Izzi", "Izzi", "Izzi"),
                       state = c("Benue", "Benue", "Benue", "Cross River", "Cross River", "Cross River", "Ebonyi", "Ebonyi", "Ebonyi"))
project_CRS = "EPSG:4326"
utm_nigeria_CRS = "EPSG:26332"

# --- 2. DATA PREPARATION ---
hh_sf <- serodata |>
  mutate(household_id = str_remove(id, "\\-\\d+$"),
         y = case_when(interpretation == "Positive" ~ 1,
                       TRUE ~ 0)) |>
  group_by(household_id) |>
  summarise(
    n_tested = n(),
    n_positive = sum(y),
    any_positive = as.integer(any(y == 1)),
    prop_positive = mean(y)) |>
  ungroup() |>
  left_join(hh_data |> select(household_id, longitude, latitude), by = "household_id") |>
  mutate(abbreviation = str_split(household_id, "-", simplify = TRUE)[, 1]) |>
  left_join(village_state |> select(abbreviation, village), by = "abbreviation") |>
  drop_na(longitude, latitude) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = project_CRS) |>
  st_transform(crs = utm_nigeria_CRS) |>
  group_by(household_id) |>
  mutate(n_count = n(), # Create a suffix (e.g., "-a", "-b") only if the ID is duplicated
         suffix = if_else(n_count > 1, paste0("-", letters[row_number()]), "")) |>
  ungroup() |>
  # Create the final, unique ID
  mutate(household_id = paste0(household_id, suffix)) |>
  select(-n_count, -suffix)

# --- 3. SPATIAL ANALYSIS (Village by Village) ---
dist_threshold <- 300

spatial_results <- purrr::map(unique(hh_sf$village), function(v) {
  village_sf <- hh_sf |> filter(village == v)
  if (nrow(village_sf) < 5) return(NULL)
  
  nb <- dnearneigh(village_sf, 0, dist_threshold)
  lw <- nb2listw(nb, style = "B", zero.policy = TRUE)
  
  # UPDATED: Run tests on prop_positive
  moran_result <- moran.test(village_sf$prop_positive, lw, zero.policy = TRUE)
  gi_result <- localG(village_sf$prop_positive, lw)
  
  list(
    moran_summary = tibble(village = v, p_value = moran_result$p.value),
    gi_summary = tibble(
      household_id = village_sf$household_id,
      gi_star_z = as.numeric(scale(gi_result))
    )
  )
})

# --- 3b. ASSESS IMPACT OF dist_threshold on ANALYSIS (Village by Village) ---
# run_village_analysis <- function(village_sf, dist_threshold) {
#   # Define neighbours
#   nb <- dnearneigh(village_sf, 0, dist_threshold)
#   lw <- nb2listw(nb, style = "B", zero.policy = TRUE)
#   # Calculate subgraphs (the number of disconnected components)
#   isolated_count <- sum(spdep::card(nb) == 0)
#   # Check if the network is so sparse that no connections exist
#   if (sum(spdep::card(nb)) == 0) {
#     return(tibble(
#       Distance = dist_threshold, 
#       Village = unique(village_sf$village), 
#       Isolated_Households = isolated_count, 
#       Moran_P = NA, 
#       Status = "Fully Disconnected"
#     ))
#   }
#   moran_result <- moran.test(village_sf$prop_positive, lw, zero.policy = TRUE)
#   return(tibble(
#     Distance = dist_threshold,
#     Village = unique(village_sf$village),
#     Isolated_Households = isolated_count,
#     Moran_P = moran_result$p.value,
#     Status = if (isolated_count == 0) "Fully Connected" else "Partially Connected"
#   ))
# }
# 
# all_villages_sf_list <- hh_sf %>% group_split(village)
# distances_to_test <- c(150, 300)
# 
# comparison_results <- purrr::map_dfr(distances_to_test, function(d) {
#   purrr::map_dfr(all_villages_sf_list, function(v_sf) {
#     if (nrow(v_sf) >= 5) {
#       run_village_analysis(v_sf, d)
#     } else {
#       NULL
#     }
#   })
# })
# 
# comparison_results %>%
#   arrange(Village, Distance) %>%
#   select(Distance, Village, Isolated_Households, Moran_P, Status) %>%
#   print(n = Inf)

# --- 4. PROCESS & SAVE RESULTS ---
moran_summary <- map_dfr(spatial_results, "moran_summary")
saveRDS(moran_summary, file = here(output_dir, "moran_summary.rds"))

gi_summary <- map_dfr(spatial_results, "gi_summary") |>
  mutate(
    cluster_type = case_when(
      gi_star_z > 1.96  ~ "Hotspot",
      gi_star_z < -1.96 ~ "Coldspot",
      TRUE              ~ "Not significant"
    ),
    cluster_type = factor(cluster_type, levels = c("Hotspot", "Coldspot", "Not significant"))
  )

hh_sf_with_clusters <- hh_sf |> left_join(gi_summary, by = "household_id")

post_hoc_glm <- glm(I(n_positive > 0) ~ I(cluster_type != "Not significant"), 
                    data = hh_sf_with_clusters, 
                    family = "binomial")
glm_summary <- tidy(post_hoc_glm)
saveRDS(glm_summary, file = here(output_dir, "glm_summary.rds"))

# --- 5. CREATE & SAVE FIGURE ---
scale_positions <- c("tl", "tr", "tl", "tl", "tl", "tr", "br", "tl", "tr")

plot_list <- hh_sf_with_clusters %>%
  left_join(village_state, by = "village") %>%
  arrange(state, village) %>%
  mutate(
    village = fct_inorder(village),
    pos = if_else(n_positive >= 1, "Seropositive", "Seronegative")
  ) %>%
  group_split(village) %>%
  purrr::imap(function(data_subset, index) {
    ggplot(data = data_subset) +
      geom_sf(aes(colour = cluster_type, shape = pos, size = pos)) +
      ggspatial::annotation_scale(location = scale_positions[index], width_hint = 0.3) +
      scale_shape_manual(name = "Household Serostatus", values = c("Seropositive" = 19, "Seronegative" = 17)) +
      scale_size_manual(name = "Household Serostatus", values = c("Seropositive" = 2, "Seronegative" = 1.2)) +
      scale_color_manual(
        name = "Local Gi* Cluster",
        values = c("Hotspot" = "#E41A1C", "Coldspot" = "#377EB8", "Not significant" = "grey50"),
        drop = FALSE
      ) +
      guides(size = "none") +
      labs(title = unique(data_subset$village)) +
      theme_bw() +
      theme(axis.text = element_blank())
  })

# Assemble the final plot with a shared legend
legend <- get_legend(plot_list[[5]] + theme(legend.position = "bottom", legend.box = "horizontal"))
plots_no_legend <- map(plot_list, ~ .x + theme(legend.position = "none"))
map_grid <- plot_grid(plotlist = plots_no_legend, ncol = 3)
final_spatial_plot <- plot_grid(map_grid, legend, ncol = 1, rel_heights = c(1, 0.1))

# Save the final plot
saveRDS(final_spatial_plot, file = here(output_dir, "spatial_plot.rds"))

# --- 6. CREATE & SAVE FIGURE ---
cluster_counts <- hh_sf_with_clusters |>
  as_tibble() |>
  count(cluster_type)

print(cluster_counts)

pos_in_hotspot <- hh_sf_with_clusters |>
  as_tibble() |>
  filter(n_positive > 0) |>
  count(cluster_type)

print(pos_in_hotspot)
2/546 * 100
