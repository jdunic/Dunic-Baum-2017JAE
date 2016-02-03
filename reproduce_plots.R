#Script order to reproduce ms outputs

if (getwd() != '/Users/jillian/R_projects/Allometry') setwd('Allometry')

source('01_load.r')
source('02_clean.r')
source('03_func.r')

#source('KI_map.R')
#source(bootstrap_for_pgls.r)

source('11_pgls_with_yerror.R')

source('plots_for_ms_species_gape_heights.R')
source('plots_for_ms_species_gape_widths.R')

# Gets the maximum lengths from fishbase using rfishbase. These values were 
# manually included in Table 1.
source('get_max_lengths_from_fishbase.R')
