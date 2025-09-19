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
  filter(!(Mouse != "Control" & (is.na(Host_Chromosome) | is.na(Host_Start)))) %>%
  mutate(Host_Chromosome_clean = gsub("^Chr", "", Host_Chromosome)) %>%
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
data <- spread_close_positions(data, threshold = 50000, spread_factor = 400)

# ==========================
# 2. Plotting Setup
# ==========================

# Define all sample colors
sample_colors <- c(
  "Control" = "black",
  "SN-20d-103" = "#E41A1C",
  "SN-20d-113" = "#377EB8",
  "Cas9-sg6-158" = "#4DAF4A",
  "SN-16d-350" = "#984EA3"
)

# Create dummy control data for legend
dummy_control <- data.frame(
  x = 0, y = 0, Mouse = "Control"
)

# Function to create a subplot for a single chromosome
create_chr_subplot <- function(chr_data, chr_name, colors_to_use, show_legend = FALSE) {
  # Calculate the range based on VISUAL positions with a 50kb buffer
  min_pos <- min(chr_data$Host_Start_Visual, na.rm = TRUE) - 50000
  max_pos <- max(chr_data$Host_Start_Visual, na.rm = TRUE) + 50000
  
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
                 # MODIFICATION: Increased size for thicker bars
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
# 3. PLOT 1: Control + SN samples
# ==========================
# Filter data for Plot 1
data_plot1 <- data %>%
  filter(Mouse %in% c("Control", "SN-20d-103", "SN-20d-113", "SN-16d-350"))

# Get chromosomes with integrations in THIS plot's data
chr_with_integrations_plot1 <- unique(data_plot1$Host_Chromosome_clean[!is.na(data_plot1$Host_Chromosome_clean)])
chr_with_integrations_plot1 <- chr_with_integrations_plot1[order(as.numeric(gsub("[^0-9]", "99", chr_with_integrations_plot1)))]

# Colors for Plot 1
colors_plot1 <- sample_colors[c("Control", "SN-20d-103", "SN-20d-113", "SN-16d-350")]

# Create subplots for each chromosome
plot_list1 <- list()
for (i in seq_along(chr_with_integrations_plot1)) {
  chr <- chr_with_integrations_plot1[i]
  chr_data <- data_plot1 %>% filter(Host_Chromosome_clean == chr)
  
  # Only show legend on the first subplot
  show_legend <- (i == 1)
  
  plot_list1[[i]] <- create_chr_subplot(chr_data, chr, colors_plot1, show_legend)
}

# Combine plots using patchwork
if (length(plot_list1) > 0) {
  p1 <- wrap_plots(plot_list1, ncol = 1) + 
    plot_annotation(
      title = "",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  print(p1)
  
  # Adjust height based on number of chromosomes
  plot_height <- 2 + length(plot_list1) * 1.5
  
  ggsave("aav_integration_sites_SN_samples.png", plot = p1, 
         width = 12, height = plot_height, dpi = 300, 
         bg = "white")
  cat("Plot 1 saved as 'aav_integration_sites_SN_samples.png'\n")
}


# ==========================
# 4. PLOT 2: Control + Cas9 sample
# ==========================
# Filter data for Plot 2
data_plot2 <- data %>%
  filter(Mouse %in% c("Control", "Cas9-sg6-158"))

# Get chromosomes with integrations in THIS plot's data
chr_with_integrations_plot2 <- unique(data_plot2$Host_Chromosome_clean[!is.na(data_plot2$Host_Chromosome_clean)])
chr_with_integrations_plot2 <- chr_with_integrations_plot2[order(as.numeric(gsub("[^0-9]", "99", chr_with_integrations_plot2)))]

# Colors for Plot 2
colors_plot2 <- sample_colors[c("Control", "Cas9-sg6-158")]

# Create subplots for each chromosome
plot_list2 <- list()
for (i in seq_along(chr_with_integrations_plot2)) {
  chr <- chr_with_integrations_plot2[i]
  chr_data <- data_plot2 %>% filter(Host_Chromosome_clean == chr)
  
  # Only show legend on the first subplot
  show_legend <- (i == 1)
  
  plot_list2[[i]] <- create_chr_subplot(chr_data, chr, colors_plot2, show_legend)
}

# Combine plots using patchwork
if (length(plot_list2) > 0) {
  p2 <- wrap_plots(plot_list2, ncol = 1) + 
    plot_annotation(
      title = "",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  print(p2)
  
  # Adjust height based on number of chromosomes
  plot_height <- 2 + length(plot_list2) * 1.5
  
  ggsave("aav_integration_sites_Cas9_sample.png", plot = p2, 
         width = 12, height = plot_height, dpi = 300, 
         bg = "white")
  cat("Plot 2 saved as 'aav_integration_sites_Cas9_sample.png'\n")
}