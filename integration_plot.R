library(ggplot2)
library(dplyr)
library(patchwork)

# ==========================
# 1. Data Loading and Preparation
# ==========================

# Load chromosome length information
chr_info <- read.table("sequence_report.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(Role == "assembled-molecule", Chromosome.name %in% c(1:19, "X", "Y")) %>%
  select(Chromosome.name, Seq.length) %>%
  mutate(
    Chromosome = factor(Chromosome.name, levels = c(1:19, "X", "Y")),
    length = as.numeric(Seq.length)
  )

# Load integration site data
data <- read.csv("integration_sites.csv", stringsAsFactors = FALSE) %>%
  select(Mouse, Host_Chromosome, Host_Start) %>%
  filter(!(Mouse != "Control" & (is.na(Host_Chromosome) | is.na(Host_Start))))

# HARDCODED MODIFICATION: Create specific labels for all samples
data <- data %>%
  mutate(
    Host_Chromosome_clean = gsub("^Chr", "", Host_Chromosome),
    # Create proper labels for all samples with explicit hardcoding for Cas9
    Mouse = case_when(
      Mouse == "SN-20d-103" ~ "SN-20d-103: 5HA-Fancc-3HA",
      Mouse == "SN-20d-113" ~ "SN-20d-113: 5HA-Fancc-3HA",
      # Hardcode the three Cas9-sg6-158 samples based on their specific positions
      Mouse == "Cas9-sg6-158" & Host_Start == 63550724 ~ "Cas9-sg6-158: Fancc-3HA-ITR",  # Fancc on Chr13
      Mouse == "Cas9-sg6-158" & Host_Start == 63550004 ~ "Cas9-sg6-158: Fancc-3HA-ITR",  # 3HA-ITR on Chr13  
      Mouse == "Cas9-sg6-158" & Host_Start == 26412211 ~ "Cas9-sg6-158: sk-polyA", # sk-polyA on Chr3
      TRUE ~ Mouse
    )
  ) %>%
  left_join(chr_info, by = c("Host_Chromosome_clean" = "Chromosome.name"))

# Define a function to spread out close integration sites for visualization
spread_close_positions <- function(data, threshold, spread_factor) {
  data %>%
    # Important to arrange by position before calculating differences
    arrange(Host_Chromosome_clean, Host_Start) %>%
    group_by(Host_Chromosome_clean) %>%
    # Identify clusters of points closer than the threshold
    mutate(cluster_id = cumsum(c(0, diff(Host_Start) > threshold))) %>%
    group_by(Host_Chromosome_clean, cluster_id) %>%
    # For each cluster, calculate new visual positions
    mutate(
      n_in_cluster = n(),
      cluster_rank = row_number(),
      # Only spread points if there's more than one in the cluster
      Host_Start_Visual = if_else(
        n_in_cluster > 1,
        # Create evenly spaced positions around the mean, with consistent spacing
        mean(Host_Start) + (cluster_rank - (n_in_cluster + 1) / 2) * spread_factor,
        as.numeric(Host_Start) # Otherwise, use the original position
      )
    ) %>%
    ungroup() %>%
    # Additional step: ensure all sites within a chromosome have consistent minimal spacing
    group_by(Host_Chromosome_clean) %>%
    arrange(Host_Start_Visual) %>%
    mutate(
      # Calculate the difference between consecutive visual positions
      pos_diff = c(Inf, diff(Host_Start_Visual)),
      # If the difference is less than spread_factor, adjust to maintain consistent spacing
      Host_Start_Visual = cumsum(c(first(Host_Start_Visual), 
                                   pmax(diff(Host_Start_Visual), spread_factor)[-length(Host_Start_Visual)]))
    ) %>%
    select(-cluster_id, -n_in_cluster, -cluster_rank, -pos_diff) %>% # Clean up helper columns
    ungroup()
}

# Apply the spreading logic to the dataset
# A threshold of 50kb defines "close" sites.
# A spread_factor of 480kb determines the visual spacing between these sites.
data <- spread_close_positions(data, threshold = 50000, spread_factor = 100000)

# HARDCODED FIX: For SN samples, manually adjust positions to ensure even spacing
# Create evenly spaced positions for all 6 integration sites
data <- data %>%
  group_by(Mouse, Host_Chromosome_clean) %>%
  mutate(
    Host_Start_Visual = case_when(
      # For SN samples on Chr13, create evenly spaced positions
      Mouse == "SN-20d-103: 5HA-Fancc-3HA" & Host_Chromosome_clean == "13" ~ 
        63000000 + (row_number() - 1) * 400000,  # 2 sites: 63000000, 63400000
      Mouse == "SN-20d-113: 5HA-Fancc-3HA" & Host_Chromosome_clean == "13" ~ 
        63800000 + (row_number() - 1) * 400000,  # 4 sites: 63800000, 64200000, 64600000, 65000000
      # Keep all other positions as they were
      TRUE ~ Host_Start_Visual
    )
  ) %>%
  ungroup()

# ==========================
# 2. Plotting Setup
# ==========================

# Define all sample colors - NOW INCLUDING THE NEW LABEL
sample_colors <- c(
  "Control" = "black",
  "SN-20d-103: 5HA-Fancc-3HA" = "#E41A1C",
  "SN-20d-113: 5HA-Fancc-3HA" = "#377EB8",
  "Cas9-sg6-158: Fancc-3HA-ITR" = "#2D7F2D",           # Original green
  "Cas9-sg6-158: sk-polyA" = "#4DAF4A"           # Darker green
)

# Create dummy control data for legend
dummy_control <- data.frame(
  x = 0, y = 0, Mouse = "Control"
)

# Function to create a subplot for a single chromosome
create_chr_subplot <- function(chr_data, chr_name, colors_to_use, show_legend = FALSE) {
  min_pos <- min(chr_data$Host_Start_Visual, na.rm = TRUE) - 50000000
  max_pos <- max(chr_data$Host_Start_Visual, na.rm = TRUE) + 50000000
  
  # Ensure we don't go below 0
  min_pos <- max(0, min_pos)
  
  p <- ggplot() +
    # Chromosome bar showing only the relevant region
    geom_rect(aes(xmin = min_pos/1000, xmax = max_pos/1000, 
                  ymin = 0.5, ymax = 1.5),
              fill = "grey80") +
    
    # Invisible dummy point for Control legend
    geom_point(data = dummy_control,
               aes(x = x, y = y, color = Mouse),
               size = 0, alpha = 0, show.legend = show_legend) +
    
    # Actual integration site segments using the VISUAL position
    geom_segment(data = chr_data,
                 aes(x = Host_Start_Visual/1000, xend = Host_Start_Visual/1000,
                     y = 0.5, yend = 1.5,
                     color = Mouse),
                 # Increased size for thicker bars
                 size = 1.5, alpha = 0.9, show.legend = show_legend) +
    
    # Manual colors and legend
    scale_color_manual(
      values = colors_to_use,
      breaks = names(colors_to_use),
      limits = names(colors_to_use),
      drop = FALSE,
      name = "Mouse"
    ) +
    
    # Axis scales
    scale_x_continuous(name = "Position (kb)", 
                       limits = c(min_pos/1000, max_pos/1000),
                       expand = c(0.01, 0),
                       labels = function(x) format(as.integer(x), big.mark = ",")) +
    scale_y_continuous(limits = c(0, 2), 
                       breaks = 1,
                       labels = paste0("Chr", chr_name),
                       expand = c(0.1, 0.1)) +
    
    # Theme
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 9, face = "bold"),
          axis.title.y = element_blank(),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          plot.margin = margin(5, 5, 5, 5))
  
  if (show_legend) {
    p <- p + guides(color = guide_legend(
      override.aes = list(
        linetype = 1,
        # Make legend key thicker to match the plot
        size = 2, 
        alpha = 1
      ),
      keywidth  = unit(0.8, "cm"),
      keyheight = unit(0.4, "cm")
    ))
  }
  
  return(p)
}

# Function to create an empty placeholder subplot
create_empty_subplot <- function() {
  ggplot() +
    theme_void() +
    theme(plot.margin = margin(5, 5, 5, 5))
}

# ==========================
# 3. PLOT 1: Cas9-sg6-158 and Control
# ==========================

cat("Creating Plot 1: Cas9-sg6-158 samples and Control\n")

# Filter data for Cas9-sg6-158 and Control
data_cas9 <- data %>%
  filter(Mouse %in% c("Control", "Cas9-sg6-158: Fancc-3HA-ITR", "Cas9-sg6-158: sk-polyA"))

# Define colors for this plot
colors_cas9 <- sample_colors[c("Control", "Cas9-sg6-158: Fancc-3HA-ITR", "Cas9-sg6-158: sk-polyA")]

# Get chromosomes with integrations for Cas9
chr_with_integrations_cas9 <- unique(data_cas9$Host_Chromosome_clean[!is.na(data_cas9$Host_Chromosome_clean)])
chr_with_integrations_cas9 <- chr_with_integrations_cas9[order(as.numeric(gsub("[^0-9]", "99", chr_with_integrations_cas9)))]

# Create subplots for each chromosome
plot_list_cas9 <- list()
for (i in seq_along(chr_with_integrations_cas9)) {
  chr <- chr_with_integrations_cas9[i]
  chr_data <- data_cas9 %>% filter(Host_Chromosome_clean == chr)
  
  # Only show legend on the first subplot
  show_legend <- (i == 1)
  
  plot_list_cas9[[i]] <- create_chr_subplot(chr_data, chr, colors_cas9, show_legend)
}

# Store the number of chromosomes in Cas9 plot for reference
n_chr_cas9 <- length(plot_list_cas9)

# Combine plots for Cas9
if (length(plot_list_cas9) > 0) {
  p_cas9 <- wrap_plots(plot_list_cas9, ncol = 1) + 
    plot_annotation(
      title = "",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  print(p_cas9)
  
  # Adjust height based on number of chromosomes
  plot_height_cas9 <- 2 + length(plot_list_cas9) * 1.5
  
  ggsave("aav_integration_sites_cas9.png", plot = p_cas9, 
         width = 12, height = plot_height_cas9, dpi = 300, 
         bg = "white")
  cat("Cas9 plot saved as 'aav_integration_sites_cas9.png'\n\n")
}

# ==========================
# 4. PLOT 2: SN-20d-103, SN-20d-113, and Control
# ==========================

cat("Creating Plot 2: SN-20d samples and Control\n")

# Filter data for SN samples and Control
data_sn <- data %>%
  filter(Mouse %in% c("Control", "SN-20d-103: 5HA-Fancc-3HA", "SN-20d-113: 5HA-Fancc-3HA"))

# Define colors for this plot
colors_sn <- sample_colors[c("Control", "SN-20d-103: 5HA-Fancc-3HA", "SN-20d-113: 5HA-Fancc-3HA")]

# Get chromosomes with integrations for SN samples
chr_with_integrations_sn <- unique(data_sn$Host_Chromosome_clean[!is.na(data_sn$Host_Chromosome_clean)])
chr_with_integrations_sn <- chr_with_integrations_sn[order(as.numeric(gsub("[^0-9]", "99", chr_with_integrations_sn)))]

# Create subplots for each chromosome
plot_list_sn <- list()
for (i in seq_along(chr_with_integrations_sn)) {
  chr <- chr_with_integrations_sn[i]
  chr_data <- data_sn %>% filter(Host_Chromosome_clean == chr)
  
  # Only show legend on the first subplot
  show_legend <- (i == 1)
  
  plot_list_sn[[i]] <- create_chr_subplot(chr_data, chr, colors_sn, show_legend)
}

# Add empty placeholder subplots to match the number of subplots in Cas9 plot
n_chr_sn <- length(plot_list_sn)
if (n_chr_sn < n_chr_cas9) {
  n_empty_plots <- n_chr_cas9 - n_chr_sn
  for (j in 1:n_empty_plots) {
    plot_list_sn[[n_chr_sn + j]] <- create_empty_subplot()
  }
}

# Combine plots for SN samples
if (length(plot_list_sn) > 0) {
  p_sn <- wrap_plots(plot_list_sn, ncol = 1) + 
    plot_annotation(
      title = "",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  print(p_sn)
  
  # Use the same height calculation as Cas9 plot for consistency
  plot_height_sn <- 2 + n_chr_cas9 * 1.5
  
  ggsave("aav_integration_sites_sn.png", plot = p_sn, 
         width = 12, height = plot_height_sn, dpi = 300, 
         bg = "white")
  cat("SN samples plot saved as 'aav_integration_sites_sn.png'\n\n")
}

# ==========================
# 5. Print summaries for both plots
# ==========================

cat("=== Summary of Cas9-sg6-158 Integration Sites ===\n")
data_cas9 %>%
  group_by(Mouse) %>%
  summarise(
    Count = n(),
    Positions = paste(Host_Start, collapse = ", ")
  ) %>%
  print()

cat("\n=== Summary of SN-20d Integration Sites ===\n")
data_sn %>%
  group_by(Mouse) %>%
  summarise(
    Count = n(),
    Positions = paste(Host_Start, collapse = ", ")
  ) %>%
  print()