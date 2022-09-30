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

#clean environment
rm(list = ls())

#set working directory
setwd(getSrcDirectory()[1])
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load plot themes
source("plot_theme.R")

################################################################################################################################################

### CREATE DATA ###

#create the neighbourhood variable 
between = list(neighbourhood = c(Scampia = "Scampia", 
                                 Roma = "Roma camp",
                                 Pozzuoli = "Pozzuoli "))

#create the test variable
within = list(test = c("Visuospatial social short-term memory",
                       "Visuospatial social working memory"))

within_short = list(test = c("social-short-term",
                             "social-working"))

#create scores for each neighbourhood for visuospatial social short-term memory
scampia_vss_stm = 5.15
roma_vss_stm = 5.15
pozzuoli_vss_stm = 4.65

#create scores for each neighbourhood for visuospatial social working memory
scampia_vss_wm = 4.15
roma_vss_wm = 4.15
pozzuoli_vss_wm = 3.65

#create the mean scores for each neighbourhood on each memory test 
hood_means = list(Scampia = c(vss_stm = scampia_vss_stm, 
                              vss_wm = scampia_vss_wm),
                  Roma = c(vss_stm = roma_vss_stm,
                           vss_wm = roma_vss_wm),
                  Pozzuoli = c(vss_stm = pozzuoli_vss_stm,
                               vss_wm = pozzuoli_vss_wm))

#create the standard deviation scores for each neighbourhood on each memory test 
hood_sds = list(Scampia = c(vss_stm = 1, 
                            vss_wm = 1),
                Roma = c(vss_stm = 1,
                         vss_wm = 1),
                Pozzuoli = c(vss_stm = 1,
                             vss_wm = 0.75))

#set the correlation of test score results for each neighbourhood
hood_cors = list(Scampia = .6, Roma = .6, Pozzuoli = .7)

#create the dataframe
dat = sim_design(within_short, between, n = 50, 
                 mu = hood_means, sd = hood_sds, r = hood_cors,
                 empirical = FALSE, plot = FALSE)

#make the data long format
dat_long = melt(setDT(dat), id.vars = c("id","neighbourhood"), variable.name = "test")

#ensure no values are below 2 or above 8
dat_long$value = ifelse(dat_long$value < 2, 2,
                        ifelse(dat_long$value > 8, 8, dat_long$value))

#set contrasts of neighbourhood so that "Pozzuoli" is the reference
contrasts(dat_long$neighbourhood) = contr.treatment(3, base = 3)

#round the test scores to a whole number
dat_long$value = round(dat_long$value, 0)

### ### ###

#plot the data and save it
test_labels = c("Visuospatial social\nshort-term memory",
                "Visuospatial social\nworking memory")

hood_colors = c("#999999", "#009E73", "#F0E442")

scores_by_hood = ggplot(dat_long, aes(x = test, y = value)) + 
                  geom_violin(aes(fill = neighbourhood), trim = TRUE, position = position_dodge(0.9)) +
                  geom_boxplot(aes(fill = neighbourhood), width = 0.15, position = position_dodge(0.9)) +
                  scale_x_discrete(labels = test_labels) +
                  scale_y_continuous(limits = c(0,8), breaks = c(seq(0,8,1))) +
                  scale_fill_manual(values = hood_colors) +
                  ylab("Test score") +
                  xlab("Test type") +
                  labs(fill = "Neighbourhood") +
                  aes(ymin = 0) +
                  avenir_theme

ggsave("../plots/plotted_simulated_data_study_2.jpg", scores_by_hood, width = 10, height = 5)

################################################################################################################################################

### ANALYSES ###

#create a dataset for each test type
soc_stm_df = droplevels(subset(dat_long, dat_long$test == "social-short-term"))
soc_wm_df = droplevels(subset(dat_long, dat_long$test == "social-working"))

### ### ###

#run models for each test type (first set the contrast for neighbourhood)
contrasts(soc_stm_df$neighbourhood) = contr.treatment(3, base = 3)
soc_stm_mod = lm(value ~ neighbourhood, data = soc_stm_df)
summary(soc_stm_mod)

contrasts(soc_wm_df$neighbourhood) = contr.treatment(3, base = 3)
soc_wm_mod = lm(value ~ neighbourhood, data = soc_wm_df)
summary(soc_wm_mod)

### ### ###

library(effsize)

#create groups for each neighbourhood and test
roma_soc_stm = subset(soc_stm_df, soc_stm_df$neighbourhood == "Roma")
roma_soc_stm_list = as.numeric(roma_soc_stm$value)

pozzuoli_soc_stm = subset(soc_stm_df, soc_stm_df$neighbourhood == "Pozzuoli")
pozzuoli_soc_stm_list = as.numeric(pozzuoli_soc_stm$value)

#get Cohen's d for the comparisons between groups
cohen.d(roma_soc_stm_list, pozzuoli_soc_stm_list)

################################################################################################################################################

### DATA SIMULATION FOR POWER ANALYSIS ###

#create a function for power simulations (based on above data creation procedure); `participant_n` is the number of participants per neighbourhood
power_simulation = function(participant_n = 40,
                            pozzuoli_soc_stm = 4.65,
                            pozzuoli_soc_wm = 3.65,
                            high_stress_soc_effect = 0.5,
                            pozzuoli_stm_sd = 1,
                            pozzuoli_wm_sd = 0.75,
                            high_stress_stm_sd = 1,
                            high_stress_wm_sd = 1,
                            test_score_cors_pozzuoli = 0.7,
                            test_score_cors_high_stress = 0.6,
                            ...) {
  
  #create the neighbourhood variable 
  between = list(neighbourhood = c(Scampia = "Scampia", 
                                   Roma = "Roma camp",
                                   Pozzuoli = "Pozzuoli "))
  
  #create the test variable
  within = list(test = c("Visuospatial social short-term memory",
                         "Visuospatial social working memory"))
  
  within_short = list(test = c("social-short-term",
                               "social-working"))
  
  #create scores for each neighbourhood for social short-term memory
  scampia_soc_stm = pozzuoli_soc_stm + high_stress_soc_effect
  roma_soc_stm = pozzuoli_soc_stm + high_stress_soc_effect
  pozzuoli_soc_stm = pozzuoli_soc_stm
  
  #create scores for each neighbourhood for verbal working memory
  scampia_soc_wm = pozzuoli_soc_wm + high_stress_soc_effect
  roma_soc_wm = pozzuoli_soc_wm + high_stress_soc_effect
  pozzuoli_soc_wm = pozzuoli_soc_wm
  
  #create the mean scores for each neighbourhood on each memory test 
  hood_means = list(Scampia = c(soc_stm = scampia_soc_stm, 
                                soc_wm = scampia_soc_wm),
                    Roma = c(soc_stm = roma_soc_stm,
                             soc_wm = roma_soc_wm),
                    Pozzuoli = c(soc_stm = pozzuoli_soc_stm,
                                 soc_wm = pozzuoli_soc_wm))
  
  #create the standard deviation scores for each neighbourhood on each memory test 
  hood_sds = list(Scampia = c(soc_stm = high_stress_stm_sd, 
                              soc_wm = high_stress_wm_sd),
                  Roma = c(soc_stm = high_stress_stm_sd,
                           soc_wm = high_stress_wm_sd),
                  Pozzuoli = c(soc_stm = pozzuoli_stm_sd,
                               soc_wm = pozzuoli_wm_sd))
  
  #set the correlation of test score results for each neighbourhood
  hood_cors = list(Scampia = test_score_cors_high_stress, Roma = test_score_cors_high_stress, Pozzuoli = test_score_cors_pozzuoli)
  
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
  contrasts(dat_long$neighbourhood) = contr.treatment(3, base = 3)
  
  #round the test scores to a whole number
  dat_long$value = round(dat_long$value, 0)  
  
  ### ### ###

  #create a dataset for each test type
  soc_stm_df = droplevels(subset(dat_long, dat_long$test == "social-short-term"))
  soc_wm_df = droplevels(subset(dat_long, dat_long$test == "social-working"))

  ### ### ###
  
  #set the contrast for neighbourhood for each dataset
  contrasts(soc_stm_df$neighbourhood) = contr.treatment(3, base = 3)
  contrasts(soc_wm_df$neighbourhood) = contr.treatment(3, base = 3)
  
  #run model for each test type
  soc_stm_mod = lm(value ~ neighbourhood, data = soc_stm_df)
  soc_wm_mod = lm(value ~ neighbourhood, data = soc_wm_df)
  
  ### ### ###
  
  #create groups for one of the high-stress neighbourhoods (they are drawn from the same population) and calculate effect sizes for each outcome
  roma_soc_stm = subset(soc_stm_df, soc_stm_df$neighbourhood == "Roma")
  roma_soc_stm_list = as.numeric(roma_soc_stm$value)
  pozzuoli_soc_stm = subset(soc_stm_df, soc_stm_df$neighbourhood == "Pozzuoli")
  pozzuoli_soc_stm_list = as.numeric(pozzuoli_soc_stm$value)
  
  roma_soc_wm = subset(soc_wm_df, soc_wm_df$neighbourhood == "Roma")
  roma_soc_wm_list = as.numeric(roma_soc_wm$value)
  pozzuoli_soc_wm = subset(soc_wm_df, soc_wm_df$neighbourhood == "Pozzuoli")
  pozzuoli_soc_wm_list = as.numeric(pozzuoli_soc_wm$value)
  
  #get Cohen's d for the comparisons between groups
  effect_size_results_soc_stm = cohen.d(roma_soc_stm_list, pozzuoli_soc_stm_list)
  effect_size_soc_stm = as.numeric(effect_size_results_soc_stm$estimate)
  effect_size_soc_stm_lower = as.numeric(effect_size_results_soc_stm$conf.int[1])
  effect_size_soc_stm_upper = as.numeric(effect_size_results_soc_stm$conf.int[2])
  
  effect_size_results_soc_wm = cohen.d(roma_soc_wm_list, pozzuoli_soc_wm_list)
  effect_size_soc_wm = as.numeric(effect_size_results_soc_wm$estimate)
  effect_size_soc_wm_lower = as.numeric(effect_size_results_soc_wm$conf.int[1])
  effect_size_soc_wm_upper = as.numeric(effect_size_results_soc_wm$conf.int[2])
  
  ### ### ###

  #return a dataframe of the model results 
  soc_stm_results = broom.mixed::tidy(soc_stm_mod)
  soc_stm_results$outcome = "Visuospatial social short-term memory"
  soc_stm_results$high_stress_effect = c(NA, effect_size_soc_stm, NA)
  soc_stm_results$high_stress_effect_lower = c(NA, effect_size_soc_stm_lower, NA)
  soc_stm_results$high_stress_effect_upper = c(NA, effect_size_soc_stm_upper, NA)
  
  soc_wm_results = broom.mixed::tidy(soc_wm_mod)
  soc_wm_results$outcome = "Visuospatial social working memory"
  soc_wm_results$high_stress_effect = c(NA, effect_size_soc_wm, NA)
  soc_wm_results$high_stress_effect_lower = c(NA, effect_size_soc_wm_lower, NA)
  soc_wm_results$high_stress_effect_upper = c(NA, effect_size_soc_wm_upper, NA)
  
  do.call("rbind", list(soc_stm_results, soc_wm_results))

}

#run repeated simulations with dataset variants for each effect size
simulations = crossing(replications = 1:1000,
                       participant_n = c(50, 75, 100, 125, 150),
                       pozzuoli_soc_stm = 4.5,
                       pozzuoli_soc_wm = 3.5,
                       high_stress_soc_effect = c(0.1, 0.2, 0.3, 0.4, 0.5),
                       pozzuoli_stm_sd = 1,
                       pozzuoli_wm_sd = 0.75,
                       high_stress_stm_sd = 1,
                       high_stress_wm_sd = 1,
                       test_score_cors_pozzuoli = 0.7,
                       test_score_cors_high_stress = 0.6) %>% mutate(analysis = pmap(., power_simulation)) %>% unnest(analysis)

#save the dataframes
write.csv(simulations, "../data/study2_social_memory_data_simulations.csv")

################################################################################################################################################

### POWER ANALYSIS ###

#subset the data to each outcome type
stm_sims = subset(simulations, simulations$outcome == "Visuospatial social short-term memory")
wm_sims = subset(simulations, simulations$outcome == "Visuospatial social working memory")

#create dataset to plot power analysis for neighbourhood main effects (both high-stress neighbourhoods were estimated to have the same effect)
simulation_results_neighbourhood_stm = filter(stm_sims, term == "neighbourhood1") %>%
                                            group_by(high_stress_soc_effect, participant_n) %>% 
                                            summarise(power = mean(p.value < .05), .groups = "drop")


simulation_results_neighbourhood_wm = filter(wm_sims, term == "neighbourhood1") %>%
                                           group_by(high_stress_soc_effect, participant_n) %>% 
                                           summarise(power = mean(p.value < .05), .groups = "drop")

### ### ###

#plot power analyse for neighbourhood effect on visuospatial social short-term memory
stm_pa = ggplot(aes(as.character(high_stress_soc_effect), participant_n, fill = power), 
                data = simulation_results_neighbourhood_stm) +
          geom_tile() +
          geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
          scale_y_continuous(breaks = c(50, 75, 100, 125, 150)) +
          scale_fill_viridis_c(name = "Power",
                               limits = c(0, 1), 
                               breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                               labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
          xlab("High-stress neighbourhood effect on visuospatial social short-term memory") + 
          ylab("Participant sample size (per neighbourhood)") +
          avenir_theme

ggsave("../plots/power_analysis_visuospatial_social_short_term_memory.jpg", stm_pa, width = 10, height = 5)

### ### ###

#plot power analyse for neighbourhood effect on visuospatial social working memory
verb_pa = ggplot(aes(as.character(high_stress_soc_effect), participant_n, fill = power), 
                 data = simulation_results_neighbourhood_wm) +
           geom_tile() +
           geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
           scale_y_continuous(breaks = c(50, 75, 100, 125, 150)) +
           scale_fill_viridis_c(name = "Power",
                                limits = c(0, 1), 
                                breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
           xlab("High-stress neighbourhood effect on visuospatial social working memory") + 
           ylab("Participant sample size (per neighbourhood)") +
           avenir_theme

ggsave("../plots/power_analysis_visuospatial_social_working_memory.jpg", verb_pa, width = 10, height = 5)

################################################################################################################################################

### EFFECT SIZE ESTIMATES ###

#create datasets
simulation_results_effect_sizes_stm = stm_sims[complete.cases(stm_sims), ]
simulation_results_effect_sizes_wm = wm_sims[complete.cases(wm_sims), ]

#select a participant sample size (per neighbourhood)
neighbourhood_n = 150
effect = 0.4

#get the mean effect size and CI for the assumed effect and a sample size of 150 for short-term memory (to estimate comparing effect sizes)
assumed_stm_effect = subset(simulation_results_effect_sizes_stm,
                            simulation_results_effect_sizes_stm$participant_n == neighbourhood_n & 
                            simulation_results_effect_sizes_stm$high_stress_soc_effect == effect)

assumed_wm_effect = subset(simulation_results_effect_sizes_wm,
                           simulation_results_effect_sizes_wm$participant_n == neighbourhood_n & 
                           simulation_results_effect_sizes_wm$high_stress_soc_effect == effect)


combined = data.frame(stm_lower = assumed_stm_effect$high_stress_effect_lower,
                      stm_effect = assumed_stm_effect$high_stress_effect,
                      stm_upper = assumed_stm_effect$high_stress_effect_upper,
                      wm_lower = assumed_wm_effect$high_stress_effect_lower,
                      wm_effect = assumed_wm_effect$high_stress_effect,
                      wm_upper = assumed_wm_effect$high_stress_effect_upper)

#add a variable that is whether or not the confidence intervals overlap
combined$overlapping_CI = ifelse(combined$stm_lower <= combined$stm_upper, TRUE, FALSE)
print(paste0("PERCENTAGE OF SIGNIFICANT EFFECT SIZE DIFFERENCES: ",
      table(combined$overlapping_CI)[1] / sum(table(combined$overlapping_CI)) * 100, "%"))

print(paste0("MEAN EFEFCT SIZES: ", mean(assumed_stm_effect$high_stress_effect)))
print(paste0("MEAN EFEFCT SIZES: ", mean(assumed_wm_effect$high_stress_effect)))

################################################################################################################################################

### COMPARISON WITH STUDY 1 TEST SCORES (RELATIVE CHANGE) - EXAMPLE DATA FOR SECONDARY ANALYSIS ###

#set expected means and SDs (data will be scaled within test type) for each group and test
pozzuoli_standard_mean = 0.6
pozzuoli_standard_sd = 1
roma_standard_mean = -0.3
roma_standard_sd = 1
scampia_standard_mean = -0.3
scampia_standard_sd = 1

pozzuoli_verbal_mean = 0.6
pozzuoli_verbal_sd = 1
roma_verbal_mean = -0.3
roma_verbal_sd = 1
scampia_verbal_mean = -0.3
scampia_verbal_sd = 1

pozzuoli_social_mean = -0.6
pozzuoli_social_sd = 1
roma_social_mean = 0.3
roma_social_sd = 1
scampia_social_mean = 0.3
scampia_social_sd = 1

#create a dataframe of just 20 participants, some with two and one measures (for illustrative purposes)
data = data.frame(test_score = c(rnorm(20, pozzuoli_standard_mean, pozzuoli_standard_sd),
                                 rnorm(19, roma_standard_mean, pozzuoli_standard_sd),
                                 rnorm(18, scampia_standard_mean, scampia_standard_sd),
                                 rnorm(20, pozzuoli_verbal_mean, pozzuoli_verbal_sd),
                                 rnorm(19, roma_verbal_mean, roma_verbal_sd),
                                 rnorm(18, scampia_verbal_mean, scampia_verbal_sd),
                                 rnorm(17, pozzuoli_social_mean, pozzuoli_social_sd),
                                 rnorm(18, roma_social_mean, roma_social_sd),
                                 rnorm(16, scampia_social_mean, scampia_social_sd)),
                  test_type = c(rep("standard", 20),
                                rep("standard", 19),
                                rep("standard", 18),
                                rep("verbal", 20),
                                rep("verbal", 19),
                                rep("verbal", 18),
                                rep("social", 17),
                                rep("social", 18),
                                rep("social", 16)),
                  hood =  c(rep("Pozzuoli", 20),
                            rep("Roma", 19),
                            rep("Scampia", 18),
                            rep("Pozzuoli", 20),
                            rep("Roma", 19),
                            rep("Scampia", 18),
                            rep("Pozzuoli", 17),
                            rep("Roma", 18),
                            rep("Scampia", 16)),
                  participant_id = c(seq(1, 20, 1),
                                   seq(21, 24, 1), seq(26, 40, 1),
                                   41, 42, 43, seq(45, 56, 1), 58, 59, 60,
                                   seq(1, 20, 1),
                                   seq(21, 24, 1), seq(26, 40, 1),
                                   41, 42, 43, seq(45, 56, 1), 58, 59, 60,
                                   seq(1, 7, 1), seq(9, 13, 1), seq(15, 19, 1),
                                   21, 22, seq(24, 37, 1), 39, 40,
                                   42, 44, seq(45, 53, 1), 55, 56, 57, 59, 60)) 

#make relevant variables factors
data$hood = factor(data$hood, levels = c("Scampia", "Roma", "Pozzuoli"))
data$test_type = factor(data$test_type, levels = c("standard", "verbal", "social"))

### ### ###

#plot the data
scores_by_hood_and_test = ggplot(data, aes(x = test_type, y = test_score)) + 
                            geom_violin(aes(fill = hood), trim = TRUE, position = position_dodge(0.9)) +
                            geom_boxplot(aes(fill = hood), width = 0.15, position = position_dodge(0.9)) +
                            scale_fill_manual(values = hood_colors) +
                            scale_x_discrete(labels = c("Visuospatial\nshort-term memory", "Verbal\nshort-term memory", "Visuospatial social\nshort-term memory")) +
                            ylab("Test score") +
                            xlab("Test type") +
                            labs(fill = "Neighbourhood") +
                            aes(ymin = 0) +
                            avenir_theme

ggsave("../plots/plotted_simulated_data_for_test_type_comparisons.jpg", scores_by_hood_and_test, width = 10, height = 5)

#run the model (first set contrasts)
contrasts(data$hood) = contr.treatment(3, base = 3)
contrasts(data$test_type) = contr.treatment(3, base = 3)

comparison_model = lmer(test_score ~ hood*test_type + (1 | participant_id), data = data)
summary(comparison_model)

### ### ###

library(emmeans)
  
#compare the contrasts
contrasts = as.data.frame(summary(emmeans(comparison_model, data = data, pairwise ~ test_type | hood)))  

#report results
print(paste0("VISUOSPATIAL SHORT-TERM MEMORY TEST SCORE - VISUOSPATIAL SOCIAL SHORT-TERM MEMORY TEST SCORE FOR SCAMPIA:", 
             " estimate = ", round(contrasts$contrasts.estimate[2], 3),
             ", SE = ", round(contrasts$contrasts.SE[2], 3),
             ", t = ", round(contrasts$contrasts.t.ratio[2], 3),
             ", p = ", round(contrasts$contrasts.p.value[2], 4),
             ", 95% CI = [", round(contrasts$emmeans.lower.CL[2], 3), ", ", round(contrasts$emmeans.upper.CL[1], 3), "]"))

print(paste0("VERBAL SHORT-TERM MEMORY TEST SCORE - VISUOSPATIAL SOCIAL SHORT-TERM MEMORY TEST SCORE FOR SCAMPIA:", 
             " estimate = ", round(contrasts$contrasts.estimate[3], 3),
             ", SE = ", round(contrasts$contrasts.SE[3], 3),
             ", t = ", round(contrasts$contrasts.t.ratio[3], 3),
             ", p = ", round(contrasts$contrasts.p.value[3], 4),
             ", 95% CI = [", round(contrasts$emmeans.lower.CL[3], 3), ", ", round(contrasts$emmeans.upper.CL[1], 3), "]"))

print(paste0("VISUOSPATIAL SHORT-TERM MEMORY TEST SCORE - VISUOSPATIAL SOCIAL SHORT-TERM MEMORY TEST SCORE FOR THE ROMA CAMP:", 
             " estimate = ", round(contrasts$contrasts.estimate[5], 3),
             ", SE = ", round(contrasts$contrasts.SE[5], 3),
             ", t = ", round(contrasts$contrasts.t.ratio[5], 3),
             ", p = ", round(contrasts$contrasts.p.value[5], 4),
             ", 95% CI = [", round(contrasts$emmeans.lower.CL[5], 3), ", ", round(contrasts$emmeans.upper.CL[1], 3), "]"))

print(paste0("VERBAL SHORT-TERM MEMORY TEST SCORE - VISUOSPATIAL SOCIAL SHORT-TERM MEMORY TEST SCORE FOR THE ROMA CAMP:", 
             " estimate = ", round(contrasts$contrasts.estimate[6], 3),
             ", SE = ", round(contrasts$contrasts.SE[6], 3),
             ", t = ", round(contrasts$contrasts.t.ratio[6], 3),
             ", p = ", round(contrasts$contrasts.p.value[6], 4),
             ", 95% CI = [", round(contrasts$emmeans.lower.CL[6], 3), ", ", round(contrasts$emmeans.upper.CL[1], 3), "]"))

print(paste0("VISUOSPATIAL SHORT-TERM MEMORY TEST SCORE - VISUOSPATIAL SOCIAL SHORT-TERM MEMORY TEST SCORE FOR POZZUOLI:", 
             " estimate = ", round(contrasts$contrasts.estimate[8], 3),
             ", SE = ", round(contrasts$contrasts.SE[8], 3),
             ", t = ", round(contrasts$contrasts.t.ratio[8], 3),
             ", p = ", round(contrasts$contrasts.p.value[8], 4),
             ", 95% CI = [", round(contrasts$emmeans.lower.CL[8], 3), ", ", round(contrasts$emmeans.upper.CL[1], 3), "]"))

print(paste0("VERBAL SHORT-TERM MEMORY TEST SCORE - VISUOSPATIAL SOCIAL SHORT-TERM MEMORY TEST SCORE FOR POZZUOLI:", 
             " estimate = ", round(contrasts$contrasts.estimate[9], 3),
             ", SE = ", round(contrasts$contrasts.SE[9], 3),
             ", t = ", round(contrasts$contrasts.t.ratio[9], 3),
             ", p = ", round(contrasts$contrasts.p.value[9], 4),
             ", 95% CI = [", round(contrasts$emmeans.lower.CL[9], 3), ", ", round(contrasts$emmeans.upper.CL[1], 3), "]"))
