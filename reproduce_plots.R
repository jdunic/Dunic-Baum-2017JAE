# Load and clean data
source('01_load.r')
source('02_clean.r')

# Load functions
source('03_func.r')

# Produce map of KI (Figure 1)
source('KI_map.R')

# Produce 10000 bootstrapped samples for gape height and widths
# commented out because this is slow
# bootstrapped samples saved and loaded in '05_pgls_with_yerror.R'
#source('04_bootstrap_for_pgls.r')

# Comparison of allometric coefficient estimates across functional groups
# includes analysis that accounts for phylogenetic non-independence (Figure 2)
# and the analysis that does not account for phylogeny (Figure S1)
source('05_pgls_with_yerror.R')

# Species - level allometric analyses for gape heights and widths
# includes relative gape size plots (Figs 3 and S2)
source('06_plots_for_ms_species_gape_heights.R')
source('07_plots_for_ms_species_gape_widths.R')

# Gets the maximum lengths from fishbase using rfishbase. These values were 
# manually included in Table 1.
source('get_max_lengths_from_fishbase.R')

# Check random effects S7 - 
source('09_gape_ms_SOM_plots.R')

beepr::beep()
