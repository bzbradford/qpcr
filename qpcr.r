#### read packages ####
library(tidyverse)
library(stats)


#### define function ####

pfaffl =
  function(df, refpop, HKG) {
    require(tidyverse)
    require(stats)
    options(warn = -1)
    datname = deparse(substitute(df))
    
    # determine GOIs
    GOI = levels(df$Target)[levels(df$Target) != HKG]
    
    # generate technical rep means
    means =
      df %>%
      group_by(Location, Date, Target, Beetle) %>%
      summarise(Mean = mean(Cq, na.rm = T))
    
    gms =
      means %>%
      filter(Location == refpop) %>%
      summarise(GM = mean(Mean, na.rm = T)) %>%
      ungroup() %>%
      mutate(Location = NULL)
      
    # join primer efficiency, compute pfaffle value
    joined =
      left_join(means, gms) %>%
      mutate(Delta = GM - Mean) %>%
      left_join(primers) %>%
      mutate(Pfaffl = Efficiency ^ Delta)
    
    # spread on target genes
    pfaffl =
      joined %>%
      select(c("Location", "Date", "Target", "Beetle", "Pfaffl")) %>%
      spread(key = Target, value = Pfaffl)
    
    # divide efficiency-weighted gene expressions by housekeeping gene
    ratios =
      pfaffl %>%
      mutate_at(GOI, funs(. / get(HKG))) %>%
      select(-HKG)
    
    # mean and sd of relative gene efficiences per beetle
    meanratios =
      ratios %>%
      select(-"Beetle") %>%
      summarise_all(list("mean", "sd"), na.rm = T) %>%
      tibble(Data = rep(datname, nrow(.)), .)
    
    return(meanratios)
  }


#### read data ####

lipid = read.csv("lipid_qpcr.csv")
detox = read.csv("detox_qpcr.csv")
primers = read.csv("primers.csv")
HKG = "RP4"
GOI = levels(lipid$Target)[levels(lipid$Target) != HKG]

pfaffl(lipid, "VV", "RP4")



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


#### Stepwise calculation ####

# Cq technical rep means and grand means
qpcr.means =
  lipid %>%
  group_by(., Location, Date, Target, Beetle) %>%
  summarise(Mean = mean(Cq, na.rm = T)) %>%
  left_join(summarise(., GM = mean(Mean, na.rm = T))) %>%
  mutate(Delta = GM - Mean)
qpcr.means

# join primer efficiency, compute pfaffle value
qpcr.joined =
  qpcr.means %>%
  left_join(primers) %>%
  mutate(Pfaffl = Efficiency ^ Delta)
qpcr.joined

# spread on target genes
qpcr.pfaffl =
  qpcr.joined %>%
  select(c("Location", "Date", "Target", "Beetle", "Pfaffl")) %>%
  spread(key = Target, value = Pfaffl)
qpcr.pfaffl

# divide efficiency-weighted gene expressions by housekeeping gene
qpcr.rge =
  qpcr.pfaffl %>%
  mutate_at(GOI, funs(. / get(HKG))) %>%
  select(-HKG)

# mean and sd of relative gene efficiences per beetle
qpcr.rgemeans =
  qpcr.rge %>%
  select(-"Beetle") %>%
  summarise_all(list("mean", "sd"), na.rm = T)

qpcr.rgemeans %>% write.csv("pfaffl_out.csv")


