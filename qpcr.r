#### read packages ####

library(tidyverse)
library(stats)


#### define function ####

pfaffl =
  function(df, refpop, HKG, style = "long", write = T) {
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
    print(gms)
    
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
      summarise_all(list("mean", "sd"), na.rm = T) %>%
      cbind(Data = rep(datname, nrow(.)), .)
    
    # mean and sd of relative gene efficiences per beetle (long format)
    meanratios.long =
      ratios %>%
      gather(key = "Target", value = "Ratio", targets) %>%
      group_by(Target, add = T) %>%
      summarise(Mean = mean(Ratio, na.rm = T),
                SD = sd(Ratio, na.rm = T)) %>%
      cbind(Data = rep(datname, nrow(.)), .)
    
    if (write == TRUE) {
      write.csv(meanratios.long, paste0("out/", datname, " out (long).csv"))
      write.csv(meanratios.wide, paste0("out/", datname, " out (wide).csv"))
    }
    
    if (style == "wide") {
      return(meanratios.wide)
    } else {
      return(meanratios.long)
    }
  }


#### Do the work ####
# Justin do this

# read datasets
# primer names must start with a letter
primers = read.csv("primers.csv")
lipid = read.csv("lipid_qpcr.csv")
detox = read.csv("detox_qpcr.csv")
detox1.1 = read.csv("detox1v2.csv")
detox1.2 = read.csv("detox2.csv")


# run computations
pfaffl(lipid, "VV", "RP4")
pfaffl(detox1.2, "VV", "RP4")



# detox results
detox.results = pfaffl(detox1.2, "VV", "RP4", write = F)
dodge = position_dodge(width = .5)
detox.results %>%
  ggplot(aes(x = Location,
             y = Mean,
             ymin = Mean,
             ymax = Mean + SD,
             fill = Target)) +
  geom_col(position = dodge, color = "black", width = .5) +
  geom_errorbar(position = dodge, width = .25) +
  facet_wrap(~ Date) +
  labs(y = "Mean gene expression ratio",
       title = "Detox qPCR (Reference population = VV")


# lipid results
lipid.results = pfaffl(lipid, "VV", "RP4", write = F)
lipid.results
dodge = position_dodge(width = .5)
lipid.results %>%
  ggplot(aes(x = Location,
             y = Mean,
             ymin = Mean,
             ymax = Mean + SD,
             fill = Target)) +
  geom_col(position = dodge, color = "black", width = .5) +
  geom_errorbar(position = dodge, width = .25) +
  facet_wrap(~ Date) +
  labs(y = "Mean gene expression ratio",
       title = "Lipid qPCR (Reference population = VV")


# this is wrong
anova(lm(Mean ~ Target + Date + Location, data = lipid.results))
