library(tidyverse)
library(readxl)
library (writexl)
library(ggpubr)
library(tableone)
library(rms)
library(survival)
library(survminer)
library(oddsratio)
library(gridExtra)
library(nnet)
library(grid)

# loading files 

setwd("/Users/Zuzanna_Bien/Data/Files")

blood <- read.csv("covidaki_bloods_01dec20.csv", stringsAsFactors = FALSE)
patients<- read.csv("covidaki_pt_04feb21.csv", stringsAsFactors = FALSE)

# only use creatinine values 
bloods <- blood %>% filter (ClinicalEvent == "Creatinine Serum")

# filter out patients with ESRF and those who did not have any creatinine values recorded on admission
ESRF_ID <- patients %>% filter (esrf == 1) %>% select (ID) %>% unlist (use.names = FALSE)
no_SCr <- setdiff(blood$ID, bloods$ID)
excluded_pt <- union(no_SCr, ESRF_ID)

#final data frame after excluding the above groups 
patients_1873 <- patients %>% filter (!ID %in%excluded_pt)

