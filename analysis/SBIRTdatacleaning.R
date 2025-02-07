# Script to read in and clean SBIRT data

library(tidyverse)

here::i_am("analysis/SBIRTdatacleaning.R")


# Read in and tidy data
sbirt <- read.csv("analysis/data/sbirt.csv") |>
  rename(id = ID1, visit = VISITNUM, alc = SUM_ALCofdrinkingdays,
         heavy = Heavy_Drinkingdays, thc = SUM_THC, stim = sum_Stim,
         group = Group) |>
  mutate(group = group - 1) |>
  group_by(id) |>
  mutate(id = cur_group_id()) |>
  mutate(visit0 = ifelse(visit == 0, 1, 0),
         visit1 = ifelse(visit == 1, 1, 0),
         visit2 = ifelse(visit == 2, 1, 0),
         visit3 = ifelse(visit == 3, 1, 0),
         visit1xtrt = visit1 * group,
         visit2xtrt = visit2 * group,
         visit3xtrt = visit3 * group) |>
  select(id, visit, alc, heavy, thc, stim, group, starts_with("visit")) |>
  mutate(visit = visit + 1)

# Save to file
saveRDS(sbirt, "data/sbirt_clean.rds")
