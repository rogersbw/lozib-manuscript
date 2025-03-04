# Processing of posterior mean summaries.

# To simulate from posterior predictive distribution
library(dplyr)
library(ggplot2)

### Directory organization

here::i_am("analysis/post_hoc/mean_plots.R")


post_summary_df <- readRDS("analysis/data/processed/post_summary_means.rds")


sbirt <- readRDS("analysis/data/sbirt_clean.rds")


mean_baseline <- sbirt |> filter(visit==1) |> pull(heavy) |> mean() 

mean_followup <- sbirt |> filter(visit>1) |> group_by(group, visit) |> summarize(mean = mean(heavy, na.rm=TRUE)) 

mean_baseline
mean_followup

base_df <- data.frame(group = 0, visit = 1, mean = mean_baseline)

raw_df <- rbind(base_df, mean_followup) |>
    mutate(group = ifelse(group == 0, "Control", "Treatment")) |>
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

#Some plotting values
pd <- position_dodge(.7)

models_of_interest <- c("un", "ri")
heavy_means <- post_summary_df |> filter(outcome == "heavy", est == "avg", model %in% models_of_interest)
other_heavy_means <- post_summary_df |> filter(outcome == "heavy", est == "avg", !model %in% models_of_interest)


# Add rows to raw_df and heavy_means to connect SBIRT group to baseline
# Fix legend


mean_plot <- ggplot(data = heavy_means, aes(x = Month, y = avg)) +
  geom_line(aes(color = model, group=interaction(model, Group), linetype=Group), position = pd, linewidth = 1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = model, group=interaction(model, Group), linetype=Group) , position = pd, linewidth = 1.5, width = .5) +
  geom_point(aes(color = model, group=interaction(model, Group)), size = 5, shape = 21, fill = "white", position = pd) +
  scale_color_manual(values = c("blue", "orange")) +
  scale_fill_manual(values = c("blue", "orange")) +
  geom_line(data = raw_df, aes(x = Month, y = mean, linetype = group), linewidth = 1.5) +
  geom_point(data = raw_df, aes(x = Month, y = mean, shape = group), size = 5) +
  theme_bw()

mean_plot

  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_line(data = other_heavy_means, aes(x = visit, y = avg), color = "gray") +
  geom_point(data = other_heavy_means, aes(x = visit, y = avg),  color = "gray") +
  geom_errorbar(data = other_heavy_means, aes(x = visit, ymin = lower, ymax = upper), width = 0.2, color = "gray") +
  labs(x = "Visit", y = "Average", title = "Model Comparisons with Error Bars") +
  theme_bw()

mean_plot <- post_summary_df |> filter(outcome == "heavy", est == "avg") |>
  ggplot(aes(x = model, y = avg, group = model,
             color = ifelse(model == "UN", "UN", "Other"),
             alpha = ifelse(model == "UN", "UN", "Other"))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  scale_color_manual(values = c("UN" = "orange", "Other" = "grey")) +
  scale_alpha_manual(values = c("UN" = 1, "Other" = 0.7)) +
  theme_minimal() +
  labs(x = "Visit", y = "Average", title = "Model Comparisons with Error Bars") +
  theme(legend.title = element_blank())

mean_visit_df |> filter(outcome == "heavy", est == "avg") |> ggplot(aes(x = outcome, y = avg, fill = model)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = lower, ymax = upper), position = pd, width = 0.25) +
    labs(title = "Posterior Mean Estimates for Each Outcome and Covariate Model",
         x = "Outcome", y = "Mean Estimate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set1")
