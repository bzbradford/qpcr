library(tidyverse)
library(stats)


# read data
qpcr =
  read.csv("qpcr.csv") %>%
  mutate(Beetle = as.factor(Beetle),
         Rep = as.factor(Rep))

primers =
  read.csv("primers.csv")



# pfaffl procedure
qpcr.means =
  qpcr %>%
  group_by(Location, Date, Target, Beetle) %>%
  summarise(mean = mean(Cq))

# generate grand means  
qpcr.gms =
  qpcr.means %>%
  group_by(Location, Date, Target) %>%
  summarise(gm = mean(mean, na.rm = T))

# join grand mean back to main data, join primer efficiency, compute pfaffle value
qpcr.joined =
  qpcr.means %>%
  left_join(qpcr.gms) %>%
  mutate(delta = gm - mean) %>%
  left_join(primers) %>%
  mutate(pfaffl = Efficiency ^ delta)

# spread 
qpcr.pfaffl =
  qpcr.joined %>%
  select(c("Location", "Date", "Target", "Beetle", "pfaffl")) %>%
  spread(key = Target, value = pfaffl)

# means
qpcr.pfaffl %>%
  group_by(Location, Date) %>%
  mutate(ABC = ABC/RP4,
         CYP = CYP/RP4,
         GST = GST/RP4,
         RP4 = RP4/RP4) %>%
  summarise_at(c("ABC", "CYP", "GST", "RP4"), mean, na.rm = TRUE)

# sd
qpcr.pfaffl %>%
  group_by(Location, Date) %>%
  summarise_at(c("ABC", "CYP", "GST", "RP4"), sd, na.rm = TRUE)

