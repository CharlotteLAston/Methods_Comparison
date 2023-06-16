###
# Project: Parks - Abrolhos, Ningaloo, South West Coast 
# Data:    Indicator Species 
# Task:    Exploratory
# author:  Charlotte
# date:    Nov-Dec 2021
##

rm(list=ls())

# libraries----
# detach("package:plyr", unload=TRUE)#will error - don't worry
library(tidyr)
library(dplyr)
options(dplyr.width = Inf) #enables head() to display all coloums
library(mgcv)
library(MuMIn)
library(car)
library(doBy)
library(gplots)
library(RColorBrewer)
# library(doParallel) #this can removed?
library(doSNOW)
library(gamm4)
library(RCurl) #needed to download data from GitHub
library(FSSgam)
library(GlobalArchive)
library(ggplot2)
library(forcats)

#### SET DIRECTORIES AND READ IN DATA ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir <- paste(working.dir, "tidy_data", sep="/")
fig_dir <- paste(working.dir, "figures", sep="/")
out_dir <- paste(working.dir, "fssgam_output", sep="/")

length.groups <- c("less_50", "greater_50_less_mat", "less_mat", "greater_mat", "greater_mat_less_125", "greater_mat_125")

#### ABROLHOS PLOTS ####
setwd(data_dir)
study <- "2021-05_Abrolhos_BOSS-BRUV" 
name <- study

ab.indicators <- c("Sparidae Chrysophrys auratus", "Labridae Choerodon rubescens", "Lethrinidae Lethrinus miniatus",
                   "Glaucosomatidae Glaucosoma hebraicum", "Berycidae Centroberyx gerrardi")

ab.length <-  readRDS(paste0(name, sep="_", "dat_length.rds"))

ab.length.all <- ab.length %>% 
  filter(status %in% "No-take") %>% 
  filter(!is.na(length)) %>% 
  filter(scientific %in% ab.indicators) %>% # Make sure it's only Abrolhos indicator species 
  group_by(method, sample, Maturity) %>% 
  summarise(Abundance = length(Maturity)) %>% 
  ungroup() %>% 
  mutate(Maturity = fct_recode(Maturity, "< 50% Length Maturity" = "less_50", ">50% Maturity but < Maturity"="greater_50_less_mat",
                               "< Length Maturity" = "less_mat", "> Length Maturity" ,"> Length Maturity but < 1.25x Maturity" = "greater_mat_less_125",
                               "> 1.25x Maturity"="greater_mat_125" ))

  

ab.length.species <- ab.length %>% 
  filter(status %in% "No-take") %>% 
  filter(!is.na(length)) %>% 
  filter(scientific %in% ab.indicators) %>% # Make sure it's only Abrolhos indicator species 
  group_by(method, sample, scientific, Maturity) %>% 
  summarise(Abundance = length(Maturity)) %>% 
  ungroup() %>% 
  mutate(Maturity = fct_recode(Maturity, "< 50% Length Maturity" = "less_50", ">50% Maturity but < Maturity"="greater_50_less_mat",
                               "< Length Maturity" = "less_mat", "> Length Maturity" ,"> Length Maturity but < 1.25x Maturity" = "greater_mat_less_125",
                               "> 1.25x Maturity"="greater_mat_125" ))

  

setwd(fig_dir)

for(G in 1:6){
  
  ## Species Levels
  for(S in 1:5){
    
    species <- ab.indicators[S] 
    
    ## Species by 4
    dat <- ab.length %>% 
      filter(status %in% "No-take") %>% 
      filter(!is.na(length)) %>% 
      filter(scientific %in% ab.indicators) %>% # Make sure it's only Abrolhos indicator species 
      group_by(method, sample, scientific, Maturity) %>% 
      summarise(Abundance = length(Maturity)) %>% 
      ungroup() %>% 
      mutate(Maturity = factor(Maturity, levels = c("less_50","greater_50_less_mat", "greater_mat_less_125", "greater_mat_125"))) %>% 
      mutate(Maturity = fct_recode(Maturity, "< 50% Length Maturity" = "less_50", ">50% Maturity but < Maturity"="greater_50_less_mat",
                                   "> Length Maturity but\n< 1.25x Maturity" = "greater_mat_less_125",
                                   "> 1.25x Maturity"="greater_mat_125")) %>% 
      filter(scientific %in% species)
    

      length.plot <- dat %>% 
        mutate(method = as.factor(method)) %>% 
        ggplot(.)+
        geom_boxplot(aes(x = method, y=Abundance, fill=method), alpha=0.3)+
        geom_jitter(aes(x = method, y=Abundance, colour=method), size=1.25, alpha=0.5, width=0.1)+
        facet_grid(~Maturity, drop=FALSE)+
        stat_summary(aes(x = method, y=Abundance), fun.y=mean, geom="point", shape=15, size=2, color="#C77CFF", fill="red")+
        guides(fill=guide_legend(title="Method"))+
        theme_classic()+
        theme(legend.position = "none")+
        ylab("Abundance")+
        xlab("Method")
      
      filename <- paste0(name, sep="_",species, sep="_", "All_Groups", ".png")
      ggsave(filename = filename, plot=length.plot, dpi=300, width = 7,height=7, unit="in") 
      
      ## Species by 2
      dat <- ab.length %>% 
        filter(status %in% "No-take") %>% 
        filter(!is.na(length)) %>% 
        filter(scientific %in% ab.indicators) %>% # Make sure it's only Abrolhos indicator species 
        group_by(method, sample, scientific, Maturity2) %>% 
        summarise(Abundance = length(Maturity2)) %>% 
        ungroup() %>% 
        mutate(Maturity2 = fct_recode(Maturity2,  "< Length Maturity" = "less_mat", "> Length Maturity"="greater_mat" )) %>% 
        filter(scientific %in% species)
      
      length.plot <- dat %>% 
        mutate(method = as.factor(method)) %>% 
        ggplot(.)+
        geom_boxplot(aes(x = method, y=Abundance, fill=method), alpha=0.3)+
        geom_jitter(aes(x = method, y=Abundance, colour=method), size=1.25, alpha=0.5, width=0.1)+
        facet_grid(~Maturity2, drop=FALSE)+
        stat_summary(aes(x = method, y=Abundance), fun.y=mean, geom="point", shape=15, size=2, color="#C77CFF", fill="red")+
        guides(fill=guide_legend(title="Method"))+
        theme_classic()+
        theme(legend.position = "none")+
        ylab("Abundance")+
        xlab("Method")
      
      filename <- paste0(name, sep="_",species, sep="_", "Greater_Less", ".png")
      ggsave(filename = filename, plot=length.plot, dpi=300, width = 7,height=7, unit="in") 

    
  }
  # All species by 4
  dat <- ab.length %>% 
    filter(status %in% "No-take") %>% 
    filter(!is.na(length)) %>% 
    filter(scientific %in% ab.indicators) %>% # Make sure it's only Abrolhos indicator species 
    group_by(method, sample, Maturity) %>% 
    summarise(Abundance = length(Maturity)) %>% 
    ungroup() %>% 
    mutate(Maturity = factor(Maturity, levels = c("less_50","greater_50_less_mat", "greater_mat_less_125", "greater_mat_125"))) %>% 
    mutate(Maturity = fct_recode(Maturity, "< 50% Length Maturity" = "less_50", ">50% Maturity but < Maturity"="greater_50_less_mat",
                                 "> Length Maturity but\n< 1.25x Maturity" = "greater_mat_less_125",
                                 "> 1.25x Maturity"="greater_mat_125"))
  
  length.plot <- dat %>% 
    mutate(method = as.factor(method)) %>% 
    ggplot(.)+
    geom_boxplot(aes(x = method, y=Abundance, fill=method), alpha=0.3)+
    geom_jitter(aes(x = method, y=Abundance, colour=method), size=1.25, alpha=0.5, width=0.1)+
    facet_grid(~Maturity, drop=FALSE)+
    stat_summary(aes(x = method, y=Abundance), fun.y=mean, geom="point", shape=15, size=2, color="#C77CFF", fill="red")+
    guides(fill=guide_legend(title="Method"))+
    theme_classic()+
    theme(legend.position = "none")+
    ylab("Abundance")+
    xlab("Method")
    
    filename <- paste0(name, sep="_", "All", sep="_" ,"All_Groups", ".png")
    ggsave(filename = filename, plot=length.plot, dpi=300, width = 7,height=7, unit="in") 
    
    ## All species by 2
    dat <- ab.length %>% 
      filter(status %in% "No-take") %>% 
      filter(!is.na(length)) %>% 
      filter(scientific %in% ab.indicators) %>% # Make sure it's only Abrolhos indicator species 
      group_by(method, sample, scientific, Maturity2) %>% 
      summarise(Abundance = length(Maturity2)) %>% 
      ungroup() %>% 
      mutate(Maturity2 = fct_recode(Maturity2,  "< Length Maturity" = "less_mat", "> Length Maturity"="greater_mat" ))
    
    length.plot <- dat %>% 
      mutate(method = as.factor(method)) %>% 
      ggplot(.)+
      geom_boxplot(aes(x = method, y=Abundance, fill=method), alpha=0.3)+
      geom_jitter(aes(x = method, y=Abundance, colour=method), size=1.25, alpha=0.5, width=0.1)+
      facet_grid(~Maturity2, drop=FALSE)+
      stat_summary(aes(x = method, y=Abundance), fun.y=mean, geom="point", shape=15, size=2, color="#C77CFF", fill="red")+
      guides(fill=guide_legend(title="Method"))+
      theme_classic()+
      theme(legend.position = "none")+
      ylab("Abundance")+
      xlab("Method")
    
    filename <- paste0(name, sep="_", "All", sep="_", "Greater_Less", ".png")
    ggsave(filename = filename, plot=length.plot, dpi=300, width = 7,height=7, unit="in")  

}

#### NINGALOO PT CLOATES PLOTS ####
study <- "Ningaloo_PtCloates_BOSS-BRUV" 
name <- study
setwd(data_dir)

ni.indicators <- c("Sparidae Chrysophrys auratus", "Lethrinidae Lethrinus nebulosus",
                   "Lutjanidae Pristipomoides multidens")

ni.length <-  readRDS(paste0(name, sep="_", "dat_length.rds"))

# ni.length.all <- ni.length %>%
#   filter(!status %in% "Fished") %>%
#   filter(!is.na(length)) %>%
#   filter(scientific %in% ni.indicators) %>% 
#   pivot_longer(cols=c(less_50, greater_50_less_mat, less_mat, greater_mat, greater_mat_less_125, greater_mat_125), names_to="Maturity", values_to = "Number") %>% 
#   filter(Number == "Y") %>% 
#   group_by(method, sample, Maturity) %>% 
#   count(Number, .drop=FALSE)
# 
# ni.length.species <- ni.length %>% 
#   filter(status %in% "No-take") %>% 
#   filter(!is.na(length)) %>% 
#   filter(scientific %in% ni.indicators) %>% # Make sure it's only Abrolhos indicator species 
#   pivot_longer(cols=c(less_50, greater_50_less_mat, less_mat, greater_mat, greater_mat_less_125, greater_mat_125), names_to="Maturity", values_to = "Number") %>% 
#   filter(Number == "Y") %>% 
#   group_by(method, sample, scientific, Maturity) %>% 
#   count(Number, .drop=FALSE)

setwd(fig_dir)

for(G in 1:6){
  
  ## Species Levels
  for(S in 1:5){
    
    species <- ni.indicators[S] 
    
    ## Species by 4
    dat <- ni.length %>% 
      filter(status %in% "No-take") %>% 
      filter(!is.na(length)) %>% 
      filter(scientific %in% ni.indicators) %>% # Make sure it's only Abrolhos indicator species 
      group_by(method, sample, scientific, Maturity) %>% 
      summarise(Abundance = length(Maturity)) %>% 
      ungroup() %>% 
      mutate(Maturity = factor(Maturity, levels = c("less_50","greater_50_less_mat", "greater_mat_less_125", "greater_mat_125"))) %>% 
      mutate(Maturity = fct_recode(Maturity, "< 50% Length Maturity" = "less_50", ">50% Maturity but < Maturity"="greater_50_less_mat",
                                   "> Length Maturity but\n< 1.25x Maturity" = "greater_mat_less_125",
                                   "> 1.25x Maturity"="greater_mat_125")) %>% 
      filter(scientific %in% species)
    
    
    length.plot <- dat %>% 
      mutate(method = as.factor(method)) %>% 
      ggplot(.)+
      geom_boxplot(aes(x = method, y=Abundance, fill=method), alpha=0.3)+
      geom_jitter(aes(x = method, y=Abundance, colour=method), size=1.25, alpha=0.5, width=0.1)+
      facet_grid(~Maturity, drop=FALSE)+
      stat_summary(aes(x = method, y=Abundance), fun.y=mean, geom="point", shape=15, size=2, color="#C77CFF", fill="red")+
      guides(fill=guide_legend(title="Method"))+
      theme_classic()+
      theme(legend.position = "none")+
      ylab("Abundance")+
      xlab("Method")
    
    filename <- paste0(name, sep="_",species, sep="_", "All_Groups", ".png")
    ggsave(filename = filename, plot=length.plot, dpi=300, width = 7,height=7, unit="in") 
    
    ## Species by 2
    dat <- ni.length %>% 
      filter(status %in% "No-take") %>% 
      filter(!is.na(length)) %>% 
      filter(scientific %in% ni.indicators) %>% # Make sure it's only Abrolhos indicator species 
      group_by(method, sample, scientific, Maturity2) %>% 
      summarise(Abundance = length(Maturity2)) %>% 
      ungroup() %>% 
      mutate(Maturity2 = fct_recode(Maturity2,  "< Length Maturity" = "less_mat", "> Length Maturity"="greater_mat" )) %>% 
      filter(scientific %in% species)
    
    length.plot <- dat %>% 
      mutate(method = as.factor(method)) %>% 
      ggplot(.)+
      geom_boxplot(aes(x = method, y=Abundance, fill=method), alpha=0.3)+
      geom_jitter(aes(x = method, y=Abundance, colour=method), size=1.25, alpha=0.5, width=0.1)+
      facet_grid(~Maturity2, drop=FALSE)+
      stat_summary(aes(x = method, y=Abundance), fun.y=mean, geom="point", shape=15, size=2, color="#C77CFF", fill="red")+
      guides(fill=guide_legend(title="Method"))+
      theme_classic()+
      theme(legend.position = "none")+
      ylab("Abundance")+
      xlab("Method")
    
    filename <- paste0(name, sep="_",species, sep="_", "Greater_Less", ".png")
    ggsave(filename = filename, plot=length.plot, dpi=300, width = 7,height=7, unit="in") 
    
    
  }
  # All species by 4
  dat <- ni.length %>% 
    filter(status %in% "No-take") %>% 
    filter(!is.na(length)) %>% 
    filter(scientific %in% ni.indicators) %>% # Make sure it's only Abrolhos indicator species 
    group_by(method, sample, Maturity) %>% 
    summarise(Abundance = length(Maturity)) %>% 
    ungroup() %>% 
    mutate(Maturity = factor(Maturity, levels = c("less_50","greater_50_less_mat", "greater_mat_less_125", "greater_mat_125"))) %>% 
    mutate(Maturity = fct_recode(Maturity, "< 50% Length Maturity" = "less_50", ">50% Maturity but < Maturity"="greater_50_less_mat",
                                 "> Length Maturity but\n< 1.25x Maturity" = "greater_mat_less_125",
                                 "> 1.25x Maturity"="greater_mat_125"))
  
  length.plot <- dat %>% 
    mutate(method = as.factor(method)) %>% 
    ggplot(.)+
    geom_boxplot(aes(x = method, y=Abundance, fill=method), alpha=0.3)+
    geom_jitter(aes(x = method, y=Abundance, colour=method), size=1.25, alpha=0.5, width=0.1)+
    facet_grid(~Maturity, drop=FALSE)+
    stat_summary(aes(x = method, y=Abundance), fun.y=mean, geom="point", shape=15, size=2, color="#C77CFF", fill="red")+
    guides(fill=guide_legend(title="Method"))+
    theme_classic()+
    theme(legend.position = "none")+
    ylab("Abundance")+
    xlab("Method")
  
  filename <- paste0(name, sep="_", "All", sep="_" ,"All_Groups", ".png")
  ggsave(filename = filename, plot=length.plot, dpi=300, width = 7,height=7, unit="in") 
  
  ## All species by 2
  dat <- ni.length %>% 
    filter(status %in% "No-take") %>% 
    filter(!is.na(length)) %>% 
    filter(scientific %in% ni.indicators) %>% # Make sure it's only Abrolhos indicator species 
    group_by(method, sample, scientific, Maturity2) %>% 
    summarise(Abundance = length(Maturity2)) %>% 
    ungroup() %>% 
    mutate(Maturity2 = fct_recode(Maturity2,  "< Length Maturity" = "less_mat", "> Length Maturity"="greater_mat" ))
  
  length.plot <- dat %>% 
    mutate(method = as.factor(method)) %>% 
    ggplot(.)+
    geom_boxplot(aes(x = method, y=Abundance, fill=method), alpha=0.3)+
    geom_jitter(aes(x = method, y=Abundance, colour=method), size=1.25, alpha=0.5, width=0.1)+
    facet_grid(~Maturity2, drop=FALSE)+
    stat_summary(aes(x = method, y=Abundance), fun.y=mean, geom="point", shape=15, size=2, color="#C77CFF", fill="red")+
    guides(fill=guide_legend(title="Method"))+
    theme_classic()+
    theme(legend.position = "none")+
    ylab("Abundance")+
    xlab("Method")
  
  filename <- paste0(name, sep="_", "All", sep="_", "Greater_Less", ".png")
  ggsave(filename = filename, plot=length.plot, dpi=300, width = 7,height=7, unit="in") 
  
}
















