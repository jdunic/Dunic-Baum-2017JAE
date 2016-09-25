source('01_load.r')
source('02_clean.r')
source('03_func.r')

#-------------------------------------------------------------------------------
# Bootstrap SMA allometric coefficient estimates for each species
#-------------------------------------------------------------------------------
# Resample species with replacement
resample_spp <- function(data) {
    len <- length(data[[1]])
    N <- 1:len
    samp <- data[sample(N, size = len, replace = TRUE, prob = NULL), ]
    return(samp)
}

# Get 100 resamples - for now. This should be bumped up to 10 000
n <- 10000
spp_resamples <- rlply(.n = n, ddply(pento, .(SpeciesCode), resample_spp), .progress = 'text') 

# Species bootstrap commented out so that I don't accidentally do it again. 
# Takes a long time.

# Run an SMA for each species
# Learned how to use tryCatch because there are cases where model convergence 
# fails and these need to be passed over.
# http://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r
#boot_spp_SMA_gh <- llply(spp_resamples, function(x) {
#  spp_sma <- tryCatch(
#        sma(gh ~ SL * SpeciesCode, data = x, log = "xy", 
#                       method = "SMA", robust = T, slope.test = 2 ), 
#        error = function(e) e, 
#        warning = function(w) w
#    )
#    if (inherits(spp_sma, 'try-error')) { next }
#    return(spp_sma)
#   }, .progress = 'text')

#boot_spp_SMA_gw <- llply(spp_resamples, function(x) {
#  spp_sma <- tryCatch(
#        sma(gw ~ SL * SpeciesCode, data = x, log = "xy", 
#                       method = "SMA", robust = T, slope.test = 1), 
#        error = function(e) e, 
#        warning = function(w) w
#    )
#    if (inherits(spp_sma, 'try-error')) { next }
#    return(spp_sma)
#   }, .progress = 'text')
#beep()

#boot_spp_SMA_gh_mass <- llply(spp_resamples, function(x) {
#  spp_sma <- tryCatch(
#        sma(gh ~ wt * SpeciesCode, data = x, log = "xy", 
#                       method = "SMA", robust = T, slope.test = 3), 
#        error = function(e) e, 
#        warning = function(w) w
#    )
#    if (inherits(spp_sma, 'try-error')) { next }
#    return(spp_sma)
#   }, .progress = 'text')

#boot_spp_SMA_gw_mass <- llply(spp_resamples, function(x) {
#  spp_sma <- tryCatch(
#        sma(gw ~ wt * SpeciesCode, data = x, log = "xy", 
#                       method = "SMA", robust = T, slope.test = 3), 
#        error = function(e) e, 
#        warning = function(w) w
#    )
#    if (inherits(spp_sma, 'try-error')) { next }
#    return(spp_sma)
#   }, .progress = 'text')

#save(boot_spp_SMA_gh, file = "boot_spp_SMA_gh_10000.RData")
#save(boot_spp_SMA_gw, file = "boot_spp_SMA_gw_10000.RData")


#load('boot_spp_SMA_gw_10000.RData')

# Get the summary dataframe
boot_spp_summ_list <- list()
for (i in 1:length(boot_spp_SMA)) {
    summ <- boot_spp_SMA[[i]]$groupsummary
    if (is.null(summ)) {
        next
    }
    boot_spp_summ_list[[i]] <- summ
    
}


boot_spp_summ_list_gh <- list()
for (i in 1:length(boot_spp_SMA_gh)) {
    summ <- boot_spp_SMA_gh[[i]]$groupsummary
    if (is.null(summ)) {
        next
    }
    boot_spp_summ_list_gh[[i]] <- summ
    
}

boot_spp_summ_list_gw <- list()
for (i in 1:length(boot_spp_SMA_gw)) {
    summ <- boot_spp_SMA_gw[[i]]$groupsummary
    if (is.null(summ)) {
        next
    }
    boot_spp_summ_list_gw[[i]] <- summ  
}

boot_spp_summary <- rbind.fill(boot_spp_summ_list)

boot_spp_summary_gh <- rbind.fill(boot_spp_summ_list_gh)
boot_spp_summary_gw <- rbind.fill(boot_spp_summ_list_gw)

#write.csv(boot_spp_summary, 'species_bootstrapped_coefficients_10000.csv')
#write.csv(boot_spp_summary_gh, 'gape_height_species_bootstrapped_coefficients_10000.csv')
#write.csv(boot_spp_summary_gw, 'gape_width_species_bootstrapped_coefficients_10000.csv')

boot_spp_summ_list_gh_mass <- list()
for (i in 1:length(boot_spp_SMA_gh_mass)) {
    summ <- boot_spp_SMA_gh_mass[[i]]$groupsummary
    if (is.null(summ)) {
        next
    }
    boot_spp_summ_list_gh_mass[[i]] <- summ    
}

boot_spp_summ_list_gw_mass <- list()
for (i in 1:length(boot_spp_SMA_gw_mass)) {
    summ <- boot_spp_SMA_gw_mass[[i]]$groupsummary
    if (is.null(summ)) {
        next
    }
    boot_spp_summ_list_gw_mass[[i]] <- summ  
}

boot_spp_summary_gh_mass <- rbind.fill(boot_spp_summ_list_gh_mass)
boot_spp_summary_gw_mass <- rbind.fill(boot_spp_summ_list_gw_mass)
#write.csv(boot_spp_summary_gh_mass, 'gape_height_mass_species_bootstrapped_coefficients_10000.csv')
#write.csv(boot_spp_summary_gw_mass, 'gape_width_mass_species_bootstrapped_coefficients_10000.csv')
