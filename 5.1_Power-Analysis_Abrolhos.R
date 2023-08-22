###
# Project: Methods Comparison
# Data:    Abrolhos BOSS & BRUV fish
# Task:    Power analysis for length groups Abrolhos
# Author:  Charlotte 
# Date:    July 2023
##

rm(list=ls())

# libraries----
library(devtools)
# devtools::install_github("beckyfisher/FSSgam_package") # Run once
library(FSSgam)
library(tidyverse)
library(mgcv)
library(MuMIn)
library(car)
library(doBy)
library(doSNOW)
# devtools::install_github("UWAMEGFisheries/GlobalArchive") # Run once
library(GlobalArchive)
library(ggplot2)
library(corrr)
library(pwr)
library(effectsize)
library(broom)

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")
out_dir <- paste(working.dir, "fssgam_output", sep="/")

study <- "2021-05_Abrolhos_BOSS-BRUV" 
name <- study

##### READ IN FORMATTED DATA ####
setwd(data_dir)

dat <-  readRDS(paste0(name, sep="_", "dat_length.rds")) %>% 
  glimpse()

# Set the predictors for modeling
pred.vars <- c("depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended") 


#* GREATER LESS THAN ALL SPECIES ####
dat.response <- dat %>% 
  filter(status %in% "No-take") %>% 
  group_by(method, sample, scientific, Maturity2) %>% 
  summarise(Abundance = length(Maturity2)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Maturity2", values_from="Abundance", id_cols=c("method", "sample", "scientific")) %>% 
  mutate(less_mat = ifelse(is.na(less_mat), 0, less_mat),
         greater_mat = ifelse(is.na(greater_mat), 0, greater_mat)) %>% 
  pivot_longer(cols= c(less_mat, greater_mat), names_to="Maturity2", values_to="Abundance") %>% 
  mutate(Maturity2 = fct_recode(Maturity2,  "< Length Maturity" = "less_mat", "> Length Maturity"="greater_mat" )) %>% 
  group_by(method, sample, Maturity2) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  glimpse()

dat.preds <- dat %>% 
  dplyr::select("method","sample","depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended") %>% 
  distinct()

dat.greater.less <- dat.response %>% 
  inner_join(., dat.preds, by=c("sample", "method")) %>% 
  mutate(method = as.factor(method),
         sample = as.factor(sample))

#* GREATER LESS THAN BY SPECIES ####

dat.response <- dat %>% 
  filter(status %in% "No-take") %>% 
  group_by(method, sample, scientific, Maturity2) %>% 
  summarise(Abundance = length(Maturity2)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Maturity2", values_from="Abundance", id_cols=c("method", "sample", "scientific")) %>% 
  mutate(less_mat = ifelse(is.na(less_mat), 0, less_mat),
         greater_mat = ifelse(is.na(greater_mat), 0, greater_mat)) %>% 
  pivot_longer(cols= c(less_mat, greater_mat), names_to="Maturity2", values_to="Abundance") %>% 
  mutate(Maturity2 = fct_recode(Maturity2,  "< Length Maturity" = "less_mat", "> Length Maturity"="greater_mat" )) %>% 
  glimpse()

dat.preds <- dat %>% 
  dplyr::select("method","sample","depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended") %>% 
  distinct()

dat.greater.less.species <- dat.response %>% 
  inner_join(., dat.preds, by=c("sample", "method")) %>% 
  mutate(method = as.factor(method),
         sample = as.factor(sample)) %>% 
  mutate(Mat.Species = paste0(Maturity2, scientific))

#* ALL SPECIES IN FOUR GROUPS ####

dat.response <- dat %>% 
  filter(status %in% "No-take") %>% 
  group_by(method, sample, scientific, Maturity) %>% 
  summarise(Abundance = length(Maturity)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Maturity", values_from="Abundance", id_cols=c("method", "sample", "scientific")) %>% 
  mutate(greater_mat_125 = ifelse(is.na(greater_mat_125), 0, greater_mat_125),
         greater_mat_less_125 = ifelse(is.na(greater_mat_less_125), 0, greater_mat_less_125),
         greater_50_less_mat = ifelse(is.na(greater_50_less_mat), 0, greater_50_less_mat)) %>% 
  pivot_longer(cols= c(greater_mat_125, greater_mat_less_125, greater_50_less_mat), names_to="Maturity", values_to="Abundance") %>% 
  mutate(Maturity = factor(Maturity, levels = c("less_50","greater_50_less_mat", "greater_mat_less_125", "greater_mat_125"))) %>% 
  mutate(Maturity = fct_recode(Maturity, "< 50 Length Maturity" = "less_50", ">50 Maturity but < Maturity"="greater_50_less_mat",
                               "> Length Maturity but < 1.25x Maturity" = "greater_mat_less_125",
                               "> 1.25x Maturity"="greater_mat_125")) %>% 
  group_by(method, sample, Maturity) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  glimpse()

dat.preds <- dat %>% 
  dplyr::select("method","sample","depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended") %>% 
  distinct()

dat.groups <- dat.response %>% 
  inner_join(., dat.preds, by=c("sample", "method")) %>% 
  mutate(method = as.factor(method),
         sample = as.factor(sample))

#* FOUR GROUPS BY SPECIES ####
dat.response <- dat %>% 
  filter(status %in% "No-take") %>% 
  group_by(method, sample, scientific, Maturity) %>% 
  summarise(Abundance = length(Maturity)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Maturity", values_from="Abundance", id_cols=c("method", "sample", "scientific")) %>% 
  mutate(greater_mat_125 = ifelse(is.na(greater_mat_125), 0, greater_mat_125),
         greater_mat_less_125 = ifelse(is.na(greater_mat_less_125), 0, greater_mat_less_125),
         greater_50_less_mat = ifelse(is.na(greater_50_less_mat), 0, greater_50_less_mat)) %>% 
  pivot_longer(cols= c(greater_mat_125, greater_mat_less_125, greater_50_less_mat), names_to="Maturity", values_to="Abundance") %>% 
  mutate(Maturity = factor(Maturity, levels = c("less_50","greater_50_less_mat", "greater_mat_less_125", "greater_mat_125"))) %>% 
  mutate(Maturity = fct_recode(Maturity, "< 50 Length Maturity" = "less_50", ">50 Maturity but < Maturity"="greater_50_less_mat",
                               "> Length Maturity but < 1.25x Maturity" = "greater_mat_less_125",
                               "> 1.25x Maturity"="greater_mat_125")) %>%  
  glimpse()

dat.preds <- dat %>% 
  dplyr::select("method","sample","depth", "macroalgae", "biog", "mean.relief","tpi", "roughness", "detrended") %>% 
  distinct()

dat.groups.species <- dat.response %>% 
  inner_join(., dat.preds, by=c("sample", "method")) %>% 
  mutate(method = as.factor(method),
         sample = as.factor(sample)) %>% 
  mutate(Mat.Species = paste0(Maturity, scientific))

#### CALCULATE POWER FOR EACH GROUP ####
#* Greater or less than maturity all species ####
use.dat <- dat.greater.less %>% 
  filter(Maturity2 %in% c("> Length Maturity"))

boss.bruv.mean <- use.dat %>% 
  group_by(method) %>% 
  summarise(Mean = mean(Abundance))

boss.og <- rpois(65, 3.06) # draw 65 numbers from poisson distribution with a mean of 3.06, we have 65 BOSS samples (Poisson distributions have the same mean and variance)
bruv.og <- rpois(42, 5)

# Effect size of 25%
boss.25 <- rpois(65, (3.06*1.25))
bruv.25 <- rpois(42, (5*1.25))

# Effect size of 50%
boss.50 <- rpois(65, (3.06*1.5))
bruv.50 <- rpois(42, (5*1.5))
  
  
  
  
