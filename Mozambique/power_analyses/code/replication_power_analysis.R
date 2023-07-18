"
Authors: Arran J. Davis
Emails: arran.davis@anthro.ox.ac.uk | davis.arran@gmail.com
Affiliation: Social Body Lab, Institute of Human Sciences, University of Oxford
Date: 22 September 2022
"

library(extrafont)
library(faux)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(sjPlot)
library(tidyr)
library(broom.mixed)
library(purrr)
library(censReg)
library(plm)
library(reshape2)
library(data.table)
library(effsize)

#clean environment
rm(list = ls())

#set working directory
setwd(getSrcDirectory()[1])
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load plot themes
source("plot_theme.R")

################################################################################################################################################

### DATA SIMULATION FOR POWER ANALYSIS ###

#create a function for power simulations (based on above data creation procedure); `participant_n` is the number of participants per neighbourhood
power_simulation = function(participant_n = 40,
                            pozzuoli_soc_stm = 1.260  + (8 * 0.479),
                            pozzuoli_soc_wm = 0.285 + (8 * 0.489),
                            high_stress_soc_stm_effect = 0.841,
                            high_stress_soc_wm_effect = 0.499,
                            pozzuoli_stm_sd = 1,
                            pozzuoli_wm_sd = 0.75,
                            high_stress_stm_sd = 1,
                            high_stress_wm_sd = 1,
                            test_score_cors_pozzuoli = 0.7,
                            test_score_cors_high_stress = 0.6,
                            ...) {
  
  #create the neighbourhood variable 
  between = list(neighbourhood = c(Scampia = "Scampia", 
                                   Pozzuoli = "Pozzuoli "))
  
  #create the test variable
  within = list(test = c("Visuospatial social short-term memory",
                         "Visuospatial social working memory"))
  
  within_short = list(test = c("social-short-term",
                               "social-working"))
  
  #create scores for each neighbourhood for social short-term memory
  scampia_soc_stm = pozzuoli_soc_stm + high_stress_soc_stm_effect
  pozzuoli_soc_stm = pozzuoli_soc_stm
  
  #create scores for each neighbourhood for verbal working memory
  scampia_soc_wm = pozzuoli_soc_wm + high_stress_soc_wm_effect
  pozzuoli_soc_wm = pozzuoli_soc_wm
  
  #create the mean scores for each neighbourhood on each memory test 
  hood_means = list(Scampia = c(soc_stm = scampia_soc_stm, 
                                soc_wm = scampia_soc_wm),
                    Pozzuoli = c(soc_stm = pozzuoli_soc_stm,
                                 soc_wm = pozzuoli_soc_wm))
  
  #create the standard deviation scores for each neighbourhood on each memory test 
  hood_sds = list(Scampia = c(soc_stm = high_stress_stm_sd, 
                              soc_wm = high_stress_wm_sd),
                  Pozzuoli = c(soc_stm = pozzuoli_stm_sd,
                               soc_wm = pozzuoli_wm_sd))
  
  #set the correlation of test score results for each neighbourhood
  hood_cors = list(Scampia = test_score_cors_high_stress, Pozzuoli = test_score_cors_pozzuoli)
  
  #create the dataframe
  dat = sim_design(within_short, between, n = participant_n, 
                   mu = hood_means, sd = hood_sds, r = hood_cors,
                   empirical = FALSE, plot = FALSE)
  
  #make the data long format
  dat_long = melt(setDT(dat), id.vars = c("id","neighbourhood"), variable.name = "test")
  
  #ensure no values are below 2 or above 8
  dat_long$value = ifelse(dat_long$value < 2, 2,
                          ifelse(dat_long$value > 8, 8, dat_long$value))
  
  #set contrasts of neighbourhood so that "Pozzuoli" is the reference
  contrasts(dat_long$neighbourhood) = contr.treatment(2, base = 2)
  
  #round the test scores to a whole number
  dat_long$value = round(dat_long$value, 0)  
  
  ### ### ###

  #create a dataset for each test type
  soc_stm_df = droplevels(subset(dat_long, dat_long$test == "social-short-term"))
  soc_wm_df = droplevels(subset(dat_long, dat_long$test == "social-working"))

  ### ### ###
  
  #set the contrast for neighbourhood for each dataset
  contrasts(soc_stm_df$neighbourhood) = contr.treatment(2, base = 2)
  contrasts(soc_wm_df$neighbourhood) = contr.treatment(2, base = 2)
  
  #run model for each test type
  soc_stm_mod = lm(value ~ neighbourhood, data = soc_stm_df)
  soc_wm_mod = lm(value ~ neighbourhood, data = soc_wm_df)
  
  ### ### ###
  
  #create groups for one of the high-stress neighbourhoods (they are drawn from the same population) and calculate effect sizes for each outcome
  scampia_soc_stm = subset(soc_stm_df, soc_stm_df$neighbourhood == "Scampia")
  scampia_soc_stm_list = as.numeric(scampia_soc_stm$value)
  pozzuoli_soc_stm = subset(soc_stm_df, soc_stm_df$neighbourhood == "Pozzuoli")
  pozzuoli_soc_stm_list = as.numeric(pozzuoli_soc_stm$value)
  
  scampia_wm = subset(soc_wm_df, soc_wm_df$neighbourhood == "Scampia")
  scampia_soc_wm_list = as.numeric(scampia_wm$value)
  pozzuoli_soc_wm = subset(soc_wm_df, soc_wm_df$neighbourhood == "Pozzuoli")
  pozzuoli_soc_wm_list = as.numeric(pozzuoli_soc_wm$value)
  
  
  #get Cohen's d for the comparisons between groups
  effect_size_results_soc_stm = effsize::cohen.d(scampia_soc_stm_list, pozzuoli_soc_stm_list)
  effect_size_soc_stm = as.numeric(effect_size_results_soc_stm$estimate)
  effect_size_soc_stm_lower = as.numeric(effect_size_results_soc_stm$conf.int[1])
  effect_size_soc_stm_upper = as.numeric(effect_size_results_soc_stm$conf.int[2])
  
  effect_size_results_soc_wm = effsize::cohen.d(scampia_soc_wm_list, pozzuoli_soc_wm_list)
  effect_size_soc_wm = as.numeric(effect_size_results_soc_wm$estimate)
  effect_size_soc_wm_lower = as.numeric(effect_size_results_soc_wm$conf.int[1])
  effect_size_soc_wm_upper = as.numeric(effect_size_results_soc_wm$conf.int[2])

  ### ### ###

  #return a dataframe of the model results 
  soc_stm_results = broom.mixed::tidy(soc_stm_mod)
  soc_stm_results$outcome = "Social short-term memory"
  soc_stm_results$high_stress_effect = c(NA, effect_size_soc_stm)
  soc_stm_results$high_stress_effect_lower = c(NA, effect_size_soc_stm_lower)
  soc_stm_results$high_stress_effect_upper = c(NA, effect_size_soc_stm_upper)
  
  soc_wm_results = broom.mixed::tidy(soc_wm_mod)
  soc_wm_results$outcome = "Social working memory"
  soc_wm_results$high_stress_effect = c(NA, effect_size_soc_wm)
  soc_wm_results$high_stress_effect_lower = c(NA, effect_size_soc_wm_lower)
  soc_wm_results$high_stress_effect_upper = c(NA, effect_size_soc_wm_upper)
  
  do.call("rbind", list(soc_stm_results, soc_wm_results))

}

#run repeated simulations with dataset variants for each effect size
simulations = crossing(replications = 1:1000,
                       participant_n = c(15, 20, 25, 30, 35),
                       pozzuoli_soc_stm = 1.260  + (8 * 0.479),
                       pozzuoli_soc_wm = 0.285 + (8 * 0.489),
                       high_stress_soc_stm_effect = c(0.2, 0.4, 0.841, 1, 1.2),
                       high_stress_soc_wm_effect = c(0.25, 0.499, 0.75,  1, 1.25),
                       pozzuoli_stm_sd = 1,
                       pozzuoli_wm_sd = 0.75,
                       high_stress_stm_sd = 1,
                       high_stress_wm_sd = 1,
                       test_score_cors_pozzuoli = 0.7,
                       test_score_cors_high_stress = 0.6) %>% mutate(analysis = pmap(., power_simulation)) %>% unnest(analysis)

#save the dataframes
write.csv(simulations, "../data/study2_replication_social_memory_data_simulations.csv")

################################################################################################################################################

### POWER ANALYSIS ###

#subset the data to each outcome type
stm_sims = subset(simulations, simulations$outcome == "Social short-term memory")
wm_sims = subset(simulations, simulations$outcome == "Social working memory")

#create dataset to plot power analysis for neighbourhood main effects (both high-stress neighbourhoods were estimated to have the same effect)
simulation_results_neighbourhood_stm = filter(stm_sims, term == "neighbourhood1") %>%
                                              group_by(high_stress_soc_stm_effect, participant_n) %>% 
                                              summarise(power = mean(p.value < .05), .groups = "drop")


simulation_results_neighbourhood_wm = filter(wm_sims, term == "neighbourhood1") %>%
                                           group_by(high_stress_soc_wm_effect, participant_n) %>% 
                                           summarise(power = mean(p.value < .05), .groups = "drop")

### ### ###

#plot power analyse for neighbourhood effect on social short-term memory
stm_pa = ggplot(aes(as.character(high_stress_soc_stm_effect), participant_n, fill = power), 
                data = simulation_results_neighbourhood_stm) +
          geom_tile() +
          geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
          scale_y_continuous(breaks = c(15, 20, 25, 30, 35)) +
          scale_fill_viridis_c(name = "Power",
                               limits = c(0, 1), 
                               breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                               labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
          xlab("High-stress environment effect on social short-term memory") + 
          ylab("Participant sample size (per environment)") +
          avenir_theme

ggsave("../plots/power_analysis_social_short_term_memory.jpg", stm_pa, width = 10, height = 5)

### ### ###

#plot power analyse for neighbourhood effect on social working memory
wm_pa = ggplot(aes(as.character(high_stress_soc_wm_effect), participant_n, fill = power), 
               data = simulation_results_neighbourhood_wm) +
           geom_tile() +
           geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
           scale_y_continuous(breaks = c(15, 20, 25, 30, 35)) +
           scale_fill_viridis_c(name = "Power",
                                limits = c(0, 1), 
                                breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
           xlab("High-stress environment effect on social working memory") + 
           ylab("Participant sample size (per neighbourhood)") +
           avenir_theme

ggsave("../plots/power_analysis_social_working_memory.jpg", verb_pa, width = 10, height = 5)
