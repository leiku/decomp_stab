#Workflow for the manuscript "Biodiversity stabilizes ecological communities through statistical-averaging effects
# rather than compensatory dynamics"
#
#Coded by Lei Zhao (lei.zhao@cau.edu.cn) and Daniel Reuman (d294r143@ku.edu)
#
#There are additional authors that contributed to the manuscript itself but not to the codes here 
#


#####  Preparing  #####
rm(list=ls())

#folder for data, results, and figures,
dir.create("Data", showWarnings = FALSE)   #create a folder to store the cleaned data
dir.create("Figs", showWarnings = FALSE)   #create a folder to store the figs
dir.create("Results", showWarnings = FALSE)   #create a folder to store the results


#Call packages

library(ggplot2)
library(patchwork)
library(ggsci)
library(ggpmisc)  #for stat_poly_eq to show fitting equation 

library(ggrepel)  #to make labels for each point in ggplot

library(dplyr)
library(tidyr)

########################################################
#################### Data cleaning #####################
########################################################

#### For Survey set:

#Aim: convert the original data frame into a list of population matrix
#Input: "LTER-grasslands-master-2012.csv", downloaded from 
#http://pasta-s.lternet.edu/package/data/eml/edi/358/3/5a9c0695f6641848bc17fd8a43b6ec74
#Output: "Data_Survey.RDS" in Data folder

source("Data_prep_Survey.R")



#### For Cedar: BigBio:

#Aim: convert the original data frame into a list of population matrix; clean data; distinguish NA from 0
#Input: "e120_Plant aboveground biomass data.csv", downloaded from 
#https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-cdr.273.10
#Output: "Data_BigBio.Rds" in Data folder

source("Data_prep_BigBio.R")    #by Reuman



#### For Cedar: BioCON:

#Aim: convert the original data frame into a list of population matrix; clean data; distinguish NA from 0
#Input: "e141_Plant aboveground biomass data.txt", downloaded from 
#https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-cdr.302.12
#Output: "Data_BioCON.Rds" in Data folder

source("Data_prep_BioCON.R")


#### For Cedar:
#Aim: get rid of some years, and some plots, because of NAs; combine the two experiments together
#Input: "Data_BigBio.Rds", and "Data_BioCON.Rds" in Data folder
#Output: "Data_Cedar_clean.Rds" in Data folder

source("Data_Cedar_dataclean.R")



##################################################################
################### Decompose asynchrony   #######################
##################################################################

###### Survey data
#Aim: decompose community stability into 3 parts: Spop * SAE * CPE
#Call: "Fn_decomposition.R"  and  "Fn_myplot.R"
#Input: "Data_Survey.Rds" in Data folder
#Output: "result_decomposition_Survey.csv" in Results folder
#        "Fig_Survey_decomposition.tif" in Figs folder
#        "Fig_Survey_richness.tif" in Figs folder
#        "Fig_Survey_shannon.tif" in Figs folder
#        "Fig_Survey_inv.simpson.tif" in Figs folder
source("script_Survey_3parts.R", echo=TRUE)


###### Cedar data
#Aim: decompose community stability into 3 parts: Spop * SAE * CPE
#Call: "Fn_decomposition.R"  and  "Fn_myplot.R"
#Input: "Data_Cedar_clean.Rds" in Data folder
#Output: "result_decomposition_Cedar.csv" in Results folder
#        "Fig_Cedar_decomposition.tif" in Figs folder
#        "Fig_Cedar_richness.tif" in Figs folder
#        "Fig_Cedar_shannon.tif" in Figs folder
#        "Fig_Cedar_inv.simpson.tif" in Figs folder

source("script_Cedar_3parts.R", echo=TRUE)




##################################################################
################### Further decompose CPE  #######################
##################################################################


#################### Constructing fake community using monoculture data 

#Aim: decompose community stability into 4 parts: Spop * SAE * CPEenv * CPEint
#Call: "Fn_decomposition_monoculture.R"  and  "Fn_myplot.R"
#Input: "Data_Cedar_clean.Rds" in Data folder
#Output: "result_decomposition_Cedar_monoculture.csv" in Results folder
#        "Fig_SI_Cedar_4parts.tif" in Figs folder

source("script_Cedar_4parts.R", echo=TRUE)



#################### Constructing surrogate community using populations from different plots

#Aim: decompose community stability into 4 parts: Spop * SAE * CPEenv * CPEint
#Call: "Fn_decomposition_surrogate_covariance.R"  and  "Fn_myplot.R"
#Input: "Data_Survey.Rds" in Data folder
#Output: "result_decomposition_Survey_surrogate.csv" in Results folder
#        "Fig_SI_Survey_4parts.tif" in Figs folder

## Notice:  time consuming!!!  Half an hours!!!
source("script_Survey_4parts.R", echo=TRUE)



#################### make plots for SAE * CPEenv

#Aim: show relationship between community stability (and richness) with SAE * CPEenv
#Call:  "Fn_myplot.R"
#Input: "result_decomposition_Survey_surrogate.csv" in Data folder
#        "result_decomposition_Cedar_monoculture.csv" in Data folder
#Output: "Fig_SI_SAE_prod_CPEenv.tif" in Figs folder

source("script_SAE_prod_CPEenv.R", echo=TRUE)




##################################################################
################### Get significance for CPE   ###################
##################################################################

#Aim: test significance of CPE using surrogate method
#Call: "Fn_decomposition.R" 
#Input: "Data_Survey.Rds" in Data folder
#Output: directly print out two tables (one for Survey data and the other for Cedar)

## Notice:  time consuming!

source("script_AAFT.R", echo=TRUE)


