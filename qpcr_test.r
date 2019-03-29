#### Stepwise calculation for testing ####
# Justin stay away!

library(tidyverse)
library(stats)

HKG = "RP4"
GOI = levels(lipid$Target)[levels(lipid$Target) != HKG]

means =
  lipid %>%
  group_by(Location, Date, Target, Beetle) %>%
  summarise(Mean = mean(Cq, na.rm = T))
means

gms =
  means %>%
  filter(Location == "VV") %>%
  summarise(GM = mean(Mean, na.rm = T)) %>%
  ungroup() %>%
  mutate(Location = NULL)
gms


# join primer efficiency, compute pfaffle value
joined =
  left_join(means, gms) %>%
  mutate(Delta = GM - Mean) %>%
  left_join(primers) %>%
  mutate(Pfaffl = Efficiency ^ Delta)
joined

pfaffl =
  joined %>%
  select(c("Location", "Date", "Target", "Beetle", "Pfaffl")) %>%
  spread(key = Target, value = Pfaffl)
pfaffl

rge =
  pfaffl %>%
  mutate_at(GOI, funs(. / get(HKG)))
rge

meanrge.wide =
  rge %>%
  summarise_all(list("mean", "sd"), na.rm = T)
meanrge.wide

targets = c("ELO", "PLA2", "POX", "RP4")
meanrge.long =
  rge %>%
  gather(key = "Target", value = "Ratio", targets) %>%
  group_by(Target, add = TRUE) %>%
  summarise(Mean = mean(Ratio),
            SD = sd(Ratio))

meanrge.long
