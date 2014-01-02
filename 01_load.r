setwd("/Users/jillian/R_projects/Allometry")
#load(".RData")
library("ggplot2")
library('plyr')
library('gridExtra')
library(sig)
library(smatr)
library(quantreg)

# Can use load Rdata to restore workspace when sublime or REPL has been closed
#load(".RData")

# Data for all specimens where a species has more than 
# 3 samples with gape and SL data
# This file was used for my thesis:
#fish <- read.csv('Allometry Data/fish_gape_cleaned_more_than_3_jan_24.csv', 
                 # header = TRUE),

# This file will be used for the manuscript
fish <- read.csv('Allometry Data/gape_data_for_MS_Jul11_cleaned.csv', 
	header=TRUE)

# Data for prey size analysis:
prey <- 
	read.csv(
		'Allometry Data/dissection_data_for_prey_size_2013_03_10_cleaned.csv',
         header = TRUE)

# Data for prey counts:
freq <- 
	read.csv(
	'Allometry Data/dissection_data_for_prey_size_2013_03_10_invert_fish_freq.csv',
     header = TRUE)

# Fishbase maximum sizes for selected species:

#spp <- c("CA.MELA", "AP.FURC", "LU.BOHA", "LU.KASM", "CE.ARGU",
         #"CE.UROD", "VA.LOUT", "PA.ARCA", "MO.GRAN", "PA.INSU",
         #"AC.NIGR", "AC.OLIV", "CE.FLAV", "CH.SORD", "SC.FREN",
         #"SC.RUBR", "CA.TERE", "PT.TILE", "CH.VAND", "PS.BART",
         #"PS.DISP", "PS.OLIV", "CH.ORNA")
max_len <- c(1170, 700, 900, 400, 600, 280, 830, 200, 600, 300, 
             213, 350, 140, 400, 470, 700, 400, 300, 45,  90,
             95,  88, 200)
len <- c("FL", rep("TL", 17), "SL", "TL", "SL", "SL", "TL")

j_fg <- c('Pi', 'BI', 'ZP', 'He', 'C')
spp_fg <- c( rep((j_fg[1]), 7), rep(j_fg[2], 3), rep(j_fg[3], 6), rep(j_fg[4], 6), 
			 rep(j_fg[5], 1) )

j_fg <- factor(j_fg, levels=c('Pi', 'BI', 'ZP', 'He', 'C'))
spp_fg <- factor(spp_fg, levels=j_fg)

pento_order <- c("CA.MELA", "AP.FURC", "LU.BOHA", "LU.KASM", "CE.ARGU",
                 "CE.UROD", "VA.LOUT", "PA.ARCA", "MO.GRAN", "PA.INSU",
                 "AC.NIGR", "AC.OLIV", "CE.FLAV", "CH.SORD", "SC.FREN",
                 "SC.RUBR", "CA.TERE", "PT.TILE", "CH.VAND", "PS.BART",
                 "PS.DISP", "PS.OLIV", "CH.ORNA")

pento_order <- factor(pento_order, levels = pento_order)

maxL <- data.frame(SpeciesCode=pento_order, j_fg=spp_fg, maxL=max_len, length=len)


