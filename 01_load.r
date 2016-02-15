library(ggplot2)
library(plyr)
library(gridExtra)
library(grid)
library(doMC)
library(smatr)
library(quantreg)

# Data for the allometric analysis:
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

# Data taken mannually from fishbase maximum sizes for selected species:
max_len <- c(1170, 700, 900, 400, 600, 280, 830, 200, 600, 300, 
             213, 350, 140, 400, 470, 700, 400, 300, 45,  90,
             95,  88, 200)
len <- c("FL", rep("TL", 17), "SL", "TL", "SL", "SL", "TL")

j_fg <- c('Pi', 'BI', 'ZP', 'He', 'C')
spp_fg <- c( rep((j_fg[1]), 7), rep(j_fg[2], 3), rep(j_fg[3], 6), rep(j_fg[4], 6), 
			 rep(j_fg[5], 1) )

# Setting up factor levels
j_fg <- factor(j_fg, levels=c('Pi', 'BI', 'ZP', 'He', 'C'))
spp_fg <- factor(spp_fg, levels=j_fg)

pento_order <- c("CA.MELA", "AP.FURC", "LU.BOHA", "LU.KASM", "CE.ARGU",
                 "CE.UROD", "VA.LOUT", "PA.ARCA", "MO.GRAN", "PA.INSU",
                 "AC.NIGR", "AC.OLIV", "CE.FLAV", "CH.SORD", "SC.FREN",
                 "SC.RUBR", "CA.TERE", "PT.TILE", "CH.VAND", "PS.BART",
                 "PS.DISP", "PS.OLIV", "CH.ORNA")

pento_order <- factor(pento_order, levels = pento_order)

maxL <- data.frame(SpeciesCode=pento_order, j_fg=spp_fg, maxL=max_len, length=len)


