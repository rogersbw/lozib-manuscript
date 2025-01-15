library(tidyverse)

setwd("analysis/")
library(tidyverse)

sbirt <- read.csv("data/sbirt.csv") %>% 
  rename(id = ID1, visit = VISITNUM, alc = SUM_ALCofdrinkingdays, heavy= Heavy_Drinkingdays, thc = SUM_THC, stim = sum_Stim, group=Group) %>% 
  mutate(group = group - 1) %>% 
  group_by(id) %>% mutate(id=cur_group_id()) %>% #Renumbers the ids starting at 1
  mutate( visit0 = ifelse(visit==0, 1, 0),
          visit1 = ifelse(visit==1, 1, 0),
          visit2 = ifelse(visit==2, 1, 0),
          visit3 = ifelse(visit==3, 1, 0),
          visit1xtrt = visit1*group,
          visit2xtrt = visit2*group,
          visit3xtrt = visit3*group) %>% 
  select(id, visit, alc, heavy, thc, stim, group, starts_with("visit")) %>% 
  mutate(visit = visit + 1)


# Plot of the data


inflationplot_all <- sbirt %>% 
  pivot_longer(c("alc", "heavy", "thc", "stim"), names_to = "drug", values_to = "count") %>% 
  ggplot() +
  geom_histogram(mapping=aes(x=count), bins = 90) +
  coord_cartesian(ylim = c(0, 50)) +
  facet_grid(drug ~ visit)




inflationplot_heavy <- sbirt %>% 
  ggplot() +
  geom_histogram(mapping=aes(x=heavy), bins = 90) +
  coord_cartesian(ylim = c(0, 50)) +
  facet_grid( ~ visit)


saveRDS(sbirt, "data/sbirt_clean.rds")
