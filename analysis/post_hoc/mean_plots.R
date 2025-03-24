# Processing of posterior mean summaries.

# To simulate from posterior predictive distribution
library(dplyr)
library(ggplot2)
library(extrafont)
library(ggpubr)

### Directory organization

here::i_am("analysis/post_hoc/mean_plots.R")


post_summary_df <- readRDS("analysis/data/processed/post_summary_means.rds")

post_summary_df_in_samp <- readRDS("analysis/data/processed/post_summary_means_insample.rds")


sbirt <- readRDS("analysis/data/sbirt_clean.rds")


# Let's turn this data wrangling into a funciton that will work on both the posterior summary and in-sample data
# Add code so that the SBIRT group has a baseline observation as well
# Then the plots should have one row for each outcome and one column for each treatment group

#Posterior checks

mean_baseline <- sbirt |> filter(visit == 1) |> pull(heavy) |> mean() 

mean_followup <- sbirt |> filter(visit > 1) |> group_by(group, visit) |> summarize(mean = mean(heavy, na.rm=TRUE)) |> rename(Group = group)



base_df <- data.frame(Group = c(0, 1), visit = c(1, 1), mean = c(mean_baseline, mean_baseline))

raw_df <- rbind(base_df, mean_followup) |>
  mutate(Group = ifelse(Group == 0, "Control", "SBIRT")) |>
  mutate(Month = case_when(visit == 1 ~ 0,
                           visit == 2 ~ 3,
                           visit == 3 ~ 6,
                           visit == 4 ~ 12)) |>
  mutate(Month = as.numeric(Month))

# First plot, let's plot the means with upper and lower error bars for heavy
# Let's do it with three sets of lines and error bars

# Add months column for plotting
post_summary_df <- post_summary_df |>
                   mutate(Month = case_when(visit == "visit1" ~ 0,
                          visit == "visit2" ~ 3,
                          visit == "visit3" ~ 6,
                          visit == "visit4" ~ 12,
                          visit == "visit2trt" ~ 3,
                          visit == "visit3trt" ~ 6,
                          visit == "visit4trt" ~ 12)) |>
                          mutate(Month = as.numeric(Month)) |>
                          mutate(Group = ifelse(visit %in% c("visit1", "visit2", "visit3", "visit4"), "Control", "SBIRT"))

# Add baseline for post_summary_df Treatment group

baseline <- post_summary_df |> filter(visit == "visit1")
baseline_trt <- baseline
baseline_trt$Group <- "SBIRT"

post_summary_df <- rbind(post_summary_df, baseline_trt)



#Some plotting values
pd <- position_dodge(.7)

models_of_interest <- c("un", "indcv", "ri", "ad")
outcome_of_interest <- "heavy"

best_models <- post_summary_df |>
  filter(outcome == outcome_of_interest, model %in% models_of_interest) |>
  rename(Model = model) |>
  mutate(Model = case_when(Model == "un" ~ "UN",
                           Model == "indcv" ~ "INDcv",
                           Model == "ri" ~ "RI",
                           Model == "ad" ~ "AD")) |>
  mutate(Model = factor(Model, levels = c("UN", "AD", "INDcv", "RI")))

other_heavy_means <- post_summary_df |> filter(outcome == "heavy", est == "avg", !model %in% models_of_interest)


# Add rows to raw_df and heavy_means to connect SBIRT group to baseline
# Fix legend
# Facet by group

mean_plot <- filter(best_models, est == "avg") |>
  ggplot(aes(x = Month, y = avg)) +
  geom_line(aes(color = Model, group = Model), position = pd, linewidth = .5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = Model, group = Model) , position = pd, linewidth = .5, width = .5) +
  scale_color_brewer(type = "qual", aesthetics = c("colour", "fill"), palette = 2) +
  geom_line(data = raw_df, aes(x = Month, y = mean), linewidth = .5, color = "black", linetype = "dashed") +
  geom_point(data = raw_df, aes(x = Month, y = mean), color = "black", shape = 15) +
  geom_point(aes(color = Model, group = Model), size = 1, shape = 21, fill = "white", position = pd, alpha = .5) +
  labs(x = "Month", y = "Heavy Drinking Days") +
  ggtitle("Zero-Inflated Models") +
  theme_bw(base_size = 12) +
  facet_wrap(~Group) +
  theme(legend.title = element_blank(),
  plot.title = element_text(hjust = 0.5),
  legend.position = "inside",  # New way to specify that the legend is inside
  legend.position.inside =  c(0.85, 0.7),
  legend.background = element_rect(fill = alpha("white", 0.6)))

ggsave("manuscript/figures/mean_plot.pdf", plot = mean_plot, width = 6, height = 3 )


# Now let's do the same for the zero and count models data

zero_plot <- filter(best_models, est == "theta") |>
  ggplot(aes(x = Month, y = avg)) +
  geom_line(aes(color = Model, group = Model), position = pd, linewidth = .5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = Model, group = Model) , position = pd, linewidth = .5, width = .5) +
  scale_color_brewer(type = "qual", aesthetics = c("colour", "fill"), palette = 2) +
  geom_point(aes(color = Model, group = Model), size = 1, shape = 21, fill = "white", position = pd, alpha = .5) +
  labs(x = "Month", y = "Proportion Heavy Drinkers") +
  theme_bw(base_size = 12) +
  ggtitle("Zero Models") +
  facet_wrap(~Group) +
  theme(legend.title = element_blank(),
  plot.title = element_text(hjust = 0.5),
  legend.position = "inside",  # New way to specify that the legend is inside
  legend.position.inside =  c(0.87, 0.72),
  legend.background = element_rect(fill = alpha("white", 0.6)))

ggsave("manuscript/figures/zero_plot.pdf", plot = zero_plot, width = 6, height = 3 )

## Count model

count_plot <- filter(best_models, est == "pi") |>
  mutate(avg = avg * 90, upper = upper * 90, lower = lower * 90) |>
  ggplot(aes(x = Month, y = avg)) +
  geom_line(aes(color = Model, group = Model), position = pd, linewidth = .5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = Model, group = Model) , position = pd, linewidth = .5, width = .5) +
  scale_color_brewer(type = "qual", aesthetics = c("colour", "fill"), palette = 2) +
  geom_point(aes(color = Model, group = Model), size = 1, shape = 21, fill = "white", position = pd, alpha = .5) +
  labs(x = "Month", y = "Heavy Drinking Days") +
  ggtitle("Count Models") +
  theme_bw(base_size = 12) +
  facet_wrap(~Group) +
  theme(legend.title = element_blank(),
  plot.title = element_text(hjust = 0.5),
  legend.position = "inside",  # New way to specify that the legend is inside
  legend.position.inside =  c(0.85, 0.7),
  legend.background = element_rect(fill = alpha("white", 0.6)))

ggsave("manuscript/figures/count_plot.pdf", plot = count_plot, width = 6, height = 3 )

strp_plot <- function(plot) {
  stripped_plot <- plot + theme(legend.title = element_blank(),
  legend.position = "none",
  axis.title.x = element_blank(),
  axis.text.x = element_blank()
  )

  return(stripped_plot)
}

mean_all_plot <- ggarrange(
    strp_plot(mean_plot), strp_plot(count_plot), zero_plot, nrow = 3
  )
mean_all_plot

ggsave("manuscript/figures/mean_all_plot.pdf", plot = mean_all_plot, width = 6, height = 9 )


#### Now let's make a Difference of Differences plot

# 1) I want posterior distributions of difference of differences for each time point

'
Steps:
1) Get 500 samples of the difference of differences for each time point
'

setwd("analysis/post_hoc/")
source("post_functions.r")


post <- post_means("heavy", "indcv")

DoD_indcv <- array(data = NA, c(dim(post$mu_est)[1], 3))
DoD_indcv[, 1] <- apply(post$mu_est[, , 5] - post$mu_est[, , 2], 1, mean)
DoD_indcv[, 2] <- apply(post$mu_est[, , 6] - post$mu_est[, , 3], 1, mean)
DoD_indcv[, 3] <- apply(post$mu_est[, , 7] - post$mu_est[, , 4], 1, mean)

rm(post)

post <- post_means("heavy", "ri")

DoD_ri <- array(data = NA, c(dim(post$mu_est)[1], 3))
DoD_ri[, 1] <- apply(post$mu_est[, , 5] - post$mu_est[, , 2], 1, mean)
DoD_ri[, 2] <- apply(post$mu_est[, , 6] - post$mu_est[, , 3], 1, mean)
DoD_ri[, 3] <- apply(post$mu_est[, , 7] - post$mu_est[, , 4], 1, mean)

rm(post)

setwd("../../")

bayes_p_ri <- rep(NA, 3)
bayes_p_indcv <- rep(NA, 3)
for (i in 1:3) {
  bayes_p_ri[i] <- sum(DoD_ri[, i] > 0) / length(DoD_ri[, i])
  bayes_p_indcv[i] <- sum(DoD_indcv[, i] > 0) / length(DoD_indcv[, i])
}

DoD <- data.frame(DoD = c(DoD_ri[, 1]), Model = rep("RI", length(DoD_ri[, 1])), Visit = rep("3 Month", length(DoD_ri[, 1])))

DoD <- rbind(DoD, data.frame(DoD = c(DoD_indcv[, 1]), Model = rep("INDcv", length(DoD_indcv[, 1])), Visit = rep("3 Month", length(DoD_indcv[, 1]))))
DoD <- rbind(DoD, data.frame(DoD = c(DoD_ri[, 2]), Model = rep("RI", length(DoD_ri[, 2])), Visit = rep("6 Month", length(DoD_ri[, 2]))))
DoD <- rbind(DoD, data.frame(DoD = c(DoD_indcv[, 2]), Model = rep("INDcv", length(DoD_indcv[, 2])), Visit = rep("6 Month", length(DoD_indcv[, 2]))))
DoD <- rbind(DoD, data.frame(DoD = c(DoD_ri[, 3]), Model = rep("RI", length(DoD_ri[, 3])), Visit = rep("12 Month", length(DoD_ri[, 3]))))
DoD <- rbind(DoD, data.frame(DoD = c(DoD_indcv[, 3]), Model = rep("INDcv", length(DoD_indcv[, 3])), Visit = rep("12 Month", length(DoD_indcv[, 3]))))


DoD <- DoD |> mutate(Visit = factor(Visit, levels = c("3 Month", "6 Month", "12 Month")))

DoD_plot <- ggplot(data = DoD) +
  geom_density(aes(x = DoD, fill = Model), alpha = .5, linewidth = .4) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  scale_fill_manual(values = c("blue", "orange")) +
  labs(x = "DoD Days of Heavy Drinking", y = "Density") +
  theme_bw() +
  facet_grid(Visit ~ .) +
  theme(legend.text = element_text(size = 8),
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        legend.position = "inside",
        legend.position.inside =  c(0.85, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.6)))


ggsave("manuscript/figures/DoD_heavy_plot.pdf", plot = DoD_plot, width = 6, height = 3 )

