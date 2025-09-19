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
    # Create proper labels for all samples
    Mouse = case_when(
      Mouse == "SN-20d-103" ~ "SN-20d-103: 5HA-Fancc-3HA",
      Mouse == "SN-20d-113" ~ "SN-20d-113: 5HA-Fancc-3HA",
      Mouse == "Cas9-sg6-158" & Host_Start >= 127500 & Host_Start <= 127800 ~ "Cas9-sg6-158: 5HA-Fancc",
      Mouse == "Cas9-sg6-158" ~ "Cas9-sg6-158: ITR-3HA",
      TRUE ~ Mouse
    ),
    Host_Chromosome_clean = gsub("^Chr", "", Host_Chromosome)
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
        # Center the spread points around the mean position of the original cluster
        mean(Host_Start) + (cluster_rank - (n_in_cluster + 1) / 2) * spread_factor,
        as.numeric(Host_Start) # Otherwise, use the original position
      )
    ) %>%
    ungroup() %>%
    select(-cluster_id, -n_in_cluster, -cluster_rank) # Clean up helper columns
}

# Apply the spreading logic to the dataset
# A threshold of 50kb defines "close" sites.
# A spread_factor of 25kb determines the visual spacing between these sites.
data <- spread_close_positions(data, threshold = 50000, spread_factor = 480000)

# ==========================
# 2. Plotting Setup
# ==========================

# Define all sample colors - NOW INCLUDING THE NEW LABEL
sample_colors <- c(
  "Control" = "black",
  "SN-20d-103: 5HA-Fancc-3HA" = "#E41A1C",
  "SN-20d-113: 5HA-Fancc-3HA" = "#377EB8",
  "Cas9-sg6-158: ITR-3HA" = "#4DAF4A",           # Original green
  "Cas9-sg6-158: 5HA-Fancc" = "#2D7F2D"          # Darker green
)

# Create dummy control data for legend
dummy_control <- data.frame(
  x = 0, y = 0, Mouse = "Control"
)

# Function to create a subplot for a single chromosome
create_chr_subplot <- function(chr_data, chr_name, colors_to_use, show_legend = FALSE) {
  min_pos <- min(chr_data$Host_Start_Visual, na.rm = TRUE) - 127515
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


# ==========================
# 3. COMBINED PLOT: All samples together
# ==========================

# Filter data for all samples we want to include
# Note: The Mouse column now has the modified labels for Cas9-sg6-158
data_combined <- data %>%
  filter(Mouse %in% names(sample_colors))

# Get chromosomes with integrations
chr_with_integrations <- unique(data_combined$Host_Chromosome_clean[!is.na(data_combined$Host_Chromosome_clean)])
chr_with_integrations <- chr_with_integrations[order(as.numeric(gsub("[^0-9]", "99", chr_with_integrations)))]

# Use all colors
colors_combined <- sample_colors

# Create subplots for each chromosome
plot_list <- list()
for (i in seq_along(chr_with_integrations)) {
  chr <- chr_with_integrations[i]
  chr_data <- data_combined %>% filter(Host_Chromosome_clean == chr)
  
  # Only show legend on the first subplot
  show_legend <- (i == 1)
  
  plot_list[[i]] <- create_chr_subplot(chr_data, chr, colors_combined, show_legend)
}

# Combine plots using patchwork
if (length(plot_list) > 0) {
  p_combined <- wrap_plots(plot_list, ncol = 1) + 
    plot_annotation(
      # title = "AAV Integration Sites - All Samples",
      title = "",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  print(p_combined)
  
  # Adjust height based on number of chromosomes
  plot_height <- 2 + length(plot_list) * 1.5
  
  ggsave("aav_integration_sites_all_samples.png", plot = p_combined, 
         width = 12, height = plot_height, dpi = 300, 
         bg = "white")
  cat("Combined plot saved as 'aav_integration_sites_all_samples.png'\n")
}

# Print summary to verify the labeling
cat("\n=== Summary of Integration Sites ===\n")
data_combined %>%
  group_by(Mouse) %>%
  summarise(
    Count = n(),
    Positions = paste(Host_Start, collapse = ", ")
  ) %>%
  print()