#### Plot 1: Heavy drinking at each time poinit

library(tidyverse)
library(ggpubr)
library(grid)
library(extrafont)


here::i_am("analysis/post_hoc/raw_data_plots.R")


sbirt <- readRDS("analysis/data/sbirt_clean.rds")


# Number of zeros:

zero_counts <- sbirt |> group_by(visit) |> summarize(zero = sum(heavy == 0, na.rm = TRUE))

# Plots of raw data at each time point

plot_list <- list()

visits <- c("Baseline", "3 Month", "6 Month", "12 Month")

for (i in seq_along(visits)){
  plot_list[[i]] <- sbirt |>
    filter(!is.na(heavy), visit == i) |>
    ggplot() +
    geom_histogram(aes(x = heavy), alpha = .6, position = 'identity', binwidth = 1, color="white") +
    coord_cartesian(ylim = c(0, 85)) +
    theme(legend.position = "none",
        panel.background = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_line(colour = "grey90"),
        title = element_text(size = 12),
        plot.title = element_text(hjust = .5),
        panel.border = element_rect(colour = NA, fill = NA, linewidth = .5)) +
    labs(x = "Days of Heavy Drinking", y = "Count", title = visits[i]) +
    annotate(geom = "text", x = 50, y = 50, label = paste("Zeros =", zero_counts[i, 2]), size = 4) +
    scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
}

inflation_all_plot <- ggarrange(
  plot_list[[1]] + rremove("ylab") + rremove("xlab"),
  plot_list[[2]] + rremove("ylab") + rremove("xlab"),
  plot_list[[3]] + rremove("ylab") + rremove("xlab"),
  plot_list[[4]] + rremove("ylab") + rremove("xlab"),
  nrow = 2
)

annotate_figure(inflation_all_plot, left = textGrob("Count", rot = 90, vjust = 1, gp = gpar(cex = 1.3, fontsize=10)),
                bottom = textGrob("Days of Heavy Drinking", gp = gpar(cex = 1.3, fontsize = 10)))

ggsave("manuscript/figures/inflation_heavy_plot.pdf", width = 6, height = 3 )




heavy_labs <- tibble(visit_string = factor(c("Baseline", "3 Month", "6 Month", "12 Month"), levels = c("Baseline", "3 Month", "6 Month", "12 Month")),
                     x = rep(25,4),
                     y = rep(50,4),
                     label = paste0("Zeros = ", zero_counts$zero))


inflationplot_heavy_all <-   sbirt %>% filter(!is.na(heavy)) %>% 
  mutate(visit_string = case_when(
    visit == 1 ~ "Baseline",
    visit == 2 ~ "3 Month",
    visit == 3 ~ "6 Month",
    visit == 4 ~ "12 Month"
  )) %>%
  mutate(visit_string = factor(visit_string, levels = c("Baseline", "3 Month", "6 Month", "12 Month"))) %>%
  ggplot() +
  geom_histogram(mapping=aes(x=heavy), alpha=.6, position = 'identity', binwidth = 1, fill = "#1f77b4") +
  coord_cartesian(ylim = c(0, 85)) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        text = element_text(family = "sans"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey90"),
        #title = element_text(size = 10),
        #plot.title = element_text(hjust = .5),
        #plot.subtitle = element_text(hjust = .5),
        panel.border = element_rect(colour = NA, fill = NA, linewidth = .5)) +
  labs(x = "Days of Heavy Drinking", y = "Count") +
  #annotate(geom="text",x=50, y=50, label="Zeros = 270", size=5)+
  geom_label(data=heavy_labs, aes(label = label, x = x, y = y), label.size = NA, size = 3) +
  scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~visit_string, nrow = 2) +
  theme(
    strip.text = element_text(size = rel(1.2)),
    strip.background = element_rect(fill = "white", colour = "white", size = 1)
  )


ggsave("manuscript/figures/inflation_heavy_plot.pdf", plot = inflationplot_heavy_all, width = 6, height = 3 )
