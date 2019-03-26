library(tidyverse)
library(stats)


# read data
qpcr.in =
  read.csv("qpcr.csv") %>%
  mutate(Beetle = as.factor(Beetle),
         Rep = as.factor(Rep))

primers =
  read.csv("primers.csv")



# Cq technical rep means and grand means
qpcr.means =
  qpcr.in %>%
  group_by(., Location, Date, Target, Beetle) %>%
  summarise(Mean = mean(Cq, na.rm = T)) %>%
  left_join(summarise(., GM = mean(Mean, na.rm = T))) %>%
  mutate(Delta = GM - Mean)

# join primer efficiency, compute pfaffle value
qpcr.joined =
  qpcr.means %>%
  left_join(primers) %>%
  mutate(Pfaffl = Efficiency ^ Delta)

# spread on target genes
qpcr.pfaffl =
  qpcr.joined %>%
  select(c("Location", "Date", "Target", "Beetle", "Pfaffl")) %>%
  spread(key = Target, value = Pfaffl)

# divide efficiency-weighted gene expressions by housekeeping gene
qpcr.rge =
  qpcr.pfaffl %>%
  mutate_at(c("ABC", "CYP", "GST"), funs(. / RP4)) %>%
  select(-"RP4")

# mean and sd of relative gene efficiences per beetle
qpcr.rgemeans =
  qpcr.rge %>%
  select(-"Beetle") %>%
  summarise_all(list("mean", "sd"), na.rm = T)

