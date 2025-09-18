library(ggplot2)
library(dplyr)

chr_info <- read.table("sequence_report.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(Role == "assembled-molecule", Chromosome.name %in% c(1:19, "X", "Y")) %>%
  select(Chromosome.name, Seq.length) %>%
  mutate(
    Chromosome = factor(Chromosome.name, levels = c(1:19, "X", "Y")),
    length = as.numeric(Seq.length),
    y_pos = as.numeric(Chromosome)
  )

data <- read.csv("integration_sites.csv", stringsAsFactors = FALSE) %>%
  select(Mouse, Host_Chromosome, Host_Start) %>%
  filter(!(Mouse != "Control" & (is.na(Host_Chromosome) | is.na(Host_Start)))) %>%
  mutate(Host_Chromosome_clean = gsub("^Chr", "", Host_Chromosome)) %>%
  left_join(chr_info, by = c("Host_Chromosome_clean" = "Chromosome.name"))

spread_close_positions <- function(data, threshold = 2e6, spread_distance = 1e6) {
  valid_data <- data %>% filter(!is.na(y_pos))
  
  spread_data <- valid_data %>%
    group_by(Host_Chromosome_clean) %>%
    arrange(Host_Start) %>%
    mutate(
      group = cumsum(c(1, diff(Host_Start) > threshold)),
      n_in_group = n(),
      Host_Start_spread = if_else(
        n_in_group > 1 & Host_Chromosome_clean == "13" & 
          Host_Start > 63500000 & Host_Start < 63560000,
        Host_Start + seq_along(Host_Start) * spread_distance - 
          (n() + 1) * spread_distance / 2,
        as.numeric(Host_Start)
      )
    ) %>%
    group_by(Host_Chromosome_clean, group) %>%
    mutate(
      group_size = n(),
      group_idx = row_number(),
      Host_Start_visual = if_else(
        Host_Chromosome_clean == "13" & 
          Host_Start > 63500000 & Host_Start < 63560000 & group_size > 1,
        mean(Host_Start) + (group_idx - (group_size + 1)/2) * spread_distance,
        as.numeric(Host_Start)
      )
    ) %>%
    ungroup()
  
  return(spread_data)
}

# Apply spreading
data_spread <- spread_close_positions(data, threshold = 5000, spread_distance = 4e5)

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

# ==========================
# PLOT 1: Control + SN samples
# ==========================

# Filter data for Plot 1
data_plot1 <- data_spread %>%
  filter(Mouse %in% c("Control", "SN-20d-103", "SN-20d-113", "SN-16d-350"))

# Colors for Plot 1
colors_plot1 <- sample_colors[c("Control", "SN-20d-103", "SN-20d-113", "SN-16d-350")]

p1 <- ggplot() +
  # Chromosomes as thick gray bars
  geom_segment(data = chr_info,
               aes(x = 0, xend = length/1e6, y = y_pos, yend = y_pos),
               size = 7, color = "grey80") +
  
  # Invisible dummy point for Control legend
  geom_point(data = dummy_control,
             aes(x = x, y = y, color = Mouse),
             size = 0, alpha = 0, show.legend = TRUE) +
  
  # Actual integration site segments
  geom_segment(data = data_plot1,
               aes(x = Host_Start_visual/1e6, xend = Host_Start_visual/1e6,
                   y = y_pos - 0.32, yend = y_pos + 0.32,
                   color = Mouse),
               size = 0.8, alpha = 0.9) +
  
  # Manual colors and legend
  scale_color_manual(
    values = colors_plot1,
    breaks = names(colors_plot1),
    limits = names(colors_plot1),
    drop = FALSE,
    name = "Mouse"
  ) +
  guides(color = guide_legend(
    override.aes = list(
      linetype = 1,
      size = 2,
      alpha = 1
    ),
    keywidth  = unit(0.8, "cm"),
    keyheight = unit(0.4, "cm")
  )) +
  
  # Axis scales
  scale_x_continuous(name = "Position (Mbp)", expand = c(0.01, 0)) +
  scale_y_continuous(breaks = chr_info$y_pos,
                     labels = paste0("Chr", chr_info$Chromosome),
                     expand = c(0.05, 0)) +
  
  # Theme
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9),
        legend.position = c(0.9, 0.7),
        legend.background = element_rect(fill = NA, color = NA)) +
  
  labs(title = "AAV Integration Sites in Mouse Genome - SN Samples",
       y = "Chromosome")

print(p1)

ggsave("aav_integration_sites_SN_samples.png", plot = p1, 
       width = 12, height = 8, dpi = 300, 
       bg = "white")

cat("Plot 1 saved as 'aav_integration_sites_SN_samples.png'\n")

# ==========================
# PLOT 2: Control + Cas9 sample
# ==========================

# Filter data for Plot 2
data_plot2 <- data_spread %>%
  filter(Mouse %in% c("Control", "Cas9-sg6-158"))

# Colors for Plot 2
colors_plot2 <- sample_colors[c("Control", "Cas9-sg6-158")]

p2 <- ggplot() +
  # Chromosomes as thick gray bars
  geom_segment(data = chr_info,
               aes(x = 0, xend = length/1e6, y = y_pos, yend = y_pos),
               size = 7, color = "grey80") +
  
  # Invisible dummy point for Control legend
  geom_point(data = dummy_control,
             aes(x = x, y = y, color = Mouse),
             size = 0, alpha = 0, show.legend = TRUE) +
  
  # Actual integration site segments
  geom_segment(data = data_plot2,
               aes(x = Host_Start_visual/1e6, xend = Host_Start_visual/1e6,
                   y = y_pos - 0.32, yend = y_pos + 0.32,
                   color = Mouse),
               size = 0.8, alpha = 0.9) +
  
  # Manual colors and legend
  scale_color_manual(
    values = colors_plot2,
    breaks = names(colors_plot2),
    limits = names(colors_plot2),
    drop = FALSE,
    name = "Mouse"
  ) +
  guides(color = guide_legend(
    override.aes = list(
      linetype = 1,
      size = 2,
      alpha = 1
    ),
    keywidth  = unit(0.8, "cm"),
    keyheight = unit(0.4, "cm")
  )) +
  
  # Axis scales
  scale_x_continuous(name = "Position (Mbp)", expand = c(0.01, 0)) +
  scale_y_continuous(breaks = chr_info$y_pos,
                     labels = paste0("Chr", chr_info$Chromosome),
                     expand = c(0.05, 0)) +
  
  # Theme
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9),
        legend.position = c(0.9, 0.7),
        legend.background = element_rect(fill = NA, color = NA)) +
  
  labs(title = "AAV Integration Sites in Mouse Genome - Cas9 Sample",
       y = "Chromosome")

print(p2)

ggsave("aav_integration_sites_Cas9_sample.png", plot = p2, 
       width = 12, height = 8, dpi = 300, 
       bg = "white")

cat("Plot 2 saved as 'aav_integration_sites_Cas9_sample.png'\n")

