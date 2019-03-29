#### read packages ####

library(tidyverse)
library(stats)


#### read data ####
# Justin do this

lipid = read.csv("lipid_qpcr.csv")
detox = read.csv("detox_qpcr.csv")
detox1.1 = read.csv("detox1v2.csv")
detox1.2 = read.csv("detox2.csv")
primers = read.csv("primers.csv")


#### run computations ####
# Justin do this

pfaffl(lipid, "VV", "RP4")
pfaffl(detox, "VV", "RP4")
pfaffl(detox1.1, "VV", "RP4")
pfaffl(detox1.2, "VV", "RP4")

pfaffl2(detox1.2, "VV", "RP4")



#### define function ####

pfaffl =
  function(df, refpop, HKG, style = "long", write = TRUE) {
    require(tidyverse)
    require(stats)
    options(warn = -1)
    datname = deparse(substitute(df))
    
    # determine genes of interest
    targets = levels(df$Target)
    GOI = targets[targets != HKG]
    
    # generate technical rep means
    means =
      df %>%
      group_by(Location, Date, Target, Beetle) %>%
      na.omit() %>%
      summarise(Mean = mean(Cq))
    
    # generate grand means of reference population
    gms =
      means %>%
      filter(Location == refpop) %>%
      summarise(GM = mean(Mean)) %>%
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
      select(c("Location", "Date", "Target", "Pfaffl")) %>%
      spread(key = Target, value = Pfaffl)
    
    # divide efficiency-weighted gene expressions by housekeeping gene
    ratios =
      pfaffl %>%
      mutate_at(GOI, funs(. / get(HKG)))
    
    # mean and sd of relative gene efficiences per beetle
    meanratios.wide =
      ratios %>%
      select(-"Beetle") %>%
      summarise_all(list("mean", "sd")) %>%
      cbind(Data = rep(datname, nrow(.)), .)
    
    meanratios.long =
      ratios %>%
      select(-"Beetle") %>%
      gather(key = "Target", value = "Ratio", targets)
    
    write.csv(meanratios, paste0("out/", datname, " out.csv"))
    
    return(meanratios)
  }


pfaffl2 =
  function(df, refpop, HKG, style = "wide", write = TRUE) {
    require(tidyverse)
    require(stats)
    options(warn = -1)
    datname = deparse(substitute(df))
    
    # determine genes of interest
    targets = levels(df$Target)
    GOI = targets[targets != HKG]
    
    # generate technical rep means
    means =
      df %>%
      group_by(Location, Date, Target, Beetle) %>%
      na.omit() %>%
      summarise(Mean = mean(Cq))
    
    # generate grand means of reference population
    gms =
      means %>%
      filter(Location == refpop) %>%
      summarise(GM = mean(Mean)) %>%
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
      mutate_at(GOI, funs(. / get(HKG)))
    
    # mean and sd of relative gene efficiences per beetle (wide format)
    meanratios.wide =
      ratios %>%
      select(-"Beetle") %>%
      summarise_all(list("mean", "sd")) %>%
      cbind(Data = rep(datname, nrow(.)), .)
    
    # mean and sd of relative gene efficiences per beetle (long format)
    meanratios.long =
      ratios %>%
      gather(key = "Target", value = "Ratio", targets) %>%
      group_by(Target, add = TRUE) %>%
      summarise(Mean = mean(Ratio), SD = sd(Ratio)) %>%
      cbind(Data = rep(datname, nrow(.)), .)
    
    if (write == TRUE) {
      write.csv(meanratios.long, paste0("out/", datname, " out (long).csv"))
      write.csv(meanratios.wide, paste0("out/", datname, " out (wide).csv"))
    }
    
    if (style == "long") {
      return(meanratios.long)
    } else {
      return(meanratios.wide)
    }
  }







#### Stepwise calculation for testing ####
# Justin stay away!

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
