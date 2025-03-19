# Processing of posterior mean summaries.

# To simulate from posterior predictive distribution
library(dplyr)
library(ggplot2)
library(extrafont)

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
  mutate(Group = ifelse(Group == 0, "Control", "Treatment")) |>
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
                          mutate(Group = ifelse(visit %in% c("visit1", "visit2", "visit3", "visit4"), "Control", "Treatment"))

# Add baseline for post_summary_df Treatment group

baseline <- post_summary_df |> filter(visit == "visit1")
baseline_trt <- baseline
baseline_trt$Group <- "Treatment"

post_summary_df <- rbind(post_summary_df, baseline_trt)



#Some plotting values
pd <- position_dodge(.7)

models_of_interest <- c("un", "indcv", "ri")
heavy_means <- post_summary_df |>
  filter(outcome == "heavy", est == "avg", model %in% models_of_interest) |>
  rename(Model = model) |>
  mutate(Model = case_when(Model == "un" ~ "UN",
                           Model == "indcv" ~ "INDcv",
                           Model == "ri" ~ "RI"))

other_heavy_means <- post_summary_df |> filter(outcome == "heavy", est == "avg", !model %in% models_of_interest)


# Add rows to raw_df and heavy_means to connect SBIRT group to baseline
# Fix legend
# Facet by group

mean_plot <- ggplot(data = heavy_means, aes(x = Month, y = avg)) +
  geom_line(aes(color = Model, group = Model), position = pd, linewidth = .5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = Model, group = Model) , position = pd, linewidth = .5, width = .5) +
  geom_point(aes(color = Model, group = Model), size = 1, shape = 21, fill = "white", position = pd) +
  scale_color_manual(values = c("blue", "orange", "green")) +
  scale_fill_manual(values = c("blue", "orange", "green")) +
  geom_line(data = raw_df, aes(x = Month, y = mean), linewidth = .5, color = "grey", linetype = "dashed") +
  geom_point(data = raw_df, aes(x = Month, y = mean), color = "grey", size = .5) +
  theme_bw() +
  labs(x = "Month", y = "Average Heavy Drinking Days") +
  facet_wrap(~Group) +
  theme(text = element_text(family = "sans", size = 12),
  legend.title = element_blank(),
  legend.position = "inside",  # New way to specify that the legend is inside
  legend.position.inside =  c(0.85, 0.75),
  legend.background = element_rect(fill = alpha("white", 0.6)))

ggsave("manuscript/figures/mean_plot.pdf", plot = mean_plot, width = 6, height = 3 )



#### Now let's make a Difference of Differences plot

# 1) I want posterior distributions of difference of differences for each time point

'
Steps:
1) Get 500 samples of the difference of differences for each time point
'

source("post_functions.r")
