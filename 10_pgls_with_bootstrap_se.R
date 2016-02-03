# I think the phylogenetic independent contrasts don't work for the same reason
# that we are having problems with this analysis. It requires that you use the
# mean values for standard length and gape size for each species, and calculate
# the differences, but it means that we don't have enough data points for the
# different species?
# No, the pics would have to be used to calculate the sma regression estimates
# for each functional group by itself. 

# Okay it has been confirmed - pics don't work because there are just too few 
# data points and the sma regression doesn't converge...


library(plyr)

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

# Run an SMA for each species
# Learned how to use tryCatch because there are cases where model convergence 
# fails and these need to be passed over.
# http://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r
boot_spp_SMA <- llply(spp_resamples, function(x) {
  spp_sma <- tryCatch(
        sma(ga ~ SL * SpeciesCode, data = x, log = "xy", 
                       method = "SMA", robust = T, slope.test = 2 ), 
        error = function(e) e, 
        warning = function(w) w
    )
    if (inherits(spp_sma, 'try-error')) { next }
    return(spp_sma)
   }, .progress = 'text')
beep()


# Get the summary dataframe
boot_spp_summ_list <- list()
for (i in 1:length(boot_spp_SMA)) {
    summ <- boot_spp_SMA[[i]]$groupsummary
    if (is.null(summ)) {
        next
    }
    boot_spp_summ_list[[i]] <- summ
    
}
boot_spp_summary <- rbind.fill(boot_spp_summ_list)

write.csv(boot_spp_summary, 'species_bootstrapped_coefficients_10000.csv')

#-------------------------------------------------------------------------------
# Read in the bootstrapped data that was saved from above.
#-------------------------------------------------------------------------------
boot_spp_summary <- read.csv('species_bootstrapped_coefficients.csv')

#-------------------------------------------------------------------------------
# Run pgls incorporating sampling error in the y-measurement
# SE from the bootstrapped estimates of the allometric coefficients
#-------------------------------------------------------------------------------
library(ape)
library(nlme)
library(multcomp)
library(visreg)

source('01_load.r')
source('02_clean.r')
source('03_func.r')

dev.new()
dev.new()

# Read in Jon's tree
fishTree <- read.tree('Jillian tree.txt')

# drop outgroup
fishTree_no_out <- drop.tip(fishTree, "OUTGROUP_Myxine_glutinosa")

spp_names <- fishTree_no_out$tip.label
spp_codes <- c('AC.NIGR', 'AC.OLIV', 'AP.FURC', 'LU.BOHA', 'CA.TERE', 
               'PT.TILE', 'PS.DISP', 'PS.OLIV', 'CE.ARGU', 'CE.UROD', 
               'VA.LOUT', 'CH.VAND', 'MO.GRAN', 'CH.SORD', 'SC.FREN', 
               'SC.RUBR', 'PA.ARCA', 'CA.MELA', 'PA.INSU', 'CE.FLAV', 
               'CH.ORNA')

# Create lookup dataframes
fg_lookup <- unique(data.frame('codes' = pento$Species, 'fg' = pento$j_fg))
spp_lookup_df <- data.frame('species' = spp_names, 'codes' = spp_codes)
fam_lookup <- ddply(fish, .(SpeciesCode), summarise, 'Family' = unique(Family), 
                    'Order' = unique(Order))
spp_lookup_df <- merge(spp_lookup_df, fg_lookup)

all_spp_summ_df <- merge(fam_lookup, boot_spp_summary, 
                         by.x = 'SpeciesCode', by.y = 'group')
all_spp_summ_df <- merge(spp_lookup_df, all_spp_summ_df, by.x = 'codes', 
                         by.y = 'SpeciesCode')

mean_spp_summ <- ddply(all_spp_summ_df, .(codes), summarise, 
                       'slope' = mean(Slope), 
                       'fg' = unique(fg), 
                       'species' = unique(species), 
                       'se' = sd(Slope))
row.names(mean_spp_summ) <- mean_spp_summ$species

#-------------------------------------------------------------------------------
# Prepare phylogenetic information (excluding PS.BART for now)
#-------------------------------------------------------------------------------

# Make tree ultrametric
chronos_tree1 <- chronos(fishTree_no_out, lambda = 1)
#chronos_tree0 <- chronos(fishTree_no_out, lambda = 0)

# Get the correlated error structure from the phylogeny assuming Brownian 
# motion
bm_fish <- corBrownian(phy = chronos_tree1)


# Plain GLS
#-------------------------------------------------------------------------------
# Functional group comparison - not including phylogeny
gls_spp_err <- gls(Slope ~ fg - 1, data = all_spp_summ_df, method = "ML")
gls_no_err <- gls(slope ~ fg - 1, data = mean_spp_summ, method = "ML")
dev.set(2)
visreg(gls_spp_err)
dev.set(3)
visreg(gls_no_err)

# Multiple hypothesis test for isometry
# Including bootstrap error
spp_err_iso <- glht(gls_spp_err, linfct = c('fgPi = 2', 
                           'fgBI = 2', 'fgZP = 2', 'fgHe = 2', 'fgC = 2')
)
summary(spp_err_iso)

# Not including bootstrap error
no_err_iso <- glht(gls_no_err, linfct = c('fgPi = 2', 
                           'fgBI = 2', 'fgZP = 2', 'fgHe = 2', 'fgC = 2')
)
summary(no_err_iso)


# Phylogenetic GLS
#-------------------------------------------------------------------------------
# Functional group comparison - including phylogeny
bm_gls <- gls(slope ~ fg - 1, correlation = bm_fish, data = mean_spp_summ, 
              method = "ML")

visreg(bm_gls)

pgls_iso <- glht(bm_gls, linfct = c('fgPi = 2', 'fgBI = 2', 'fgZP = 2', '
                fgHe = 2', 'fgC = 2'))
summary(pgls_iso)


# Functional group comparison - including phylogeny AND allometric coefficient 
# estimate errors
#-------------------------------------------------------------------------------

# Custom function to calculate the error that includes phylogenetic correlation
# and measurement error (from bootstrap).
# Measurement error is set to NULL as the default, such that using pgls using 
# the covariance matrix calculated in this function or using just the plain
# phylogenetic correlation structure should give the same answer.
lk <- function (sig2, y, X, C, v=NULL, opt=TRUE) {
    n <- nrow(C)

    if (is.null(v)) { v <- rep(0, n) }
    
    V <- sig2 * C + diag(v)
    beta <- solve(t(X) %*% solve(V) %*% X) %*% (t(X) %*% 
        solve(V) %*% y)
    logL <- -(1/2) * t(y - X %*% beta) %*% solve(V) %*% 
        (y - X %*% beta) - (1/2) * determinant(V, 
        logarithm = TRUE)$modulus - (n / 2) * log(2 * pi)
    
    if (opt == TRUE) { 
        return(-logL)
        } else { 
            return(list(beta = beta, sig2e = sig2, logL = logL))
        }
}

#-------------------------------------------------------------------------------
# Testing method
#-------------------------------------------------------------------------------
# Calculate design matrix for mean values
design_mat <- model.matrix(~ fg - 1, data = mean_spp_summ)
design_mat <- cbind(mean_spp_summ$slope, design_mat)

# Optimise for sig2
fit_lik <- optimize(lk, interval = c(0, 1000), y = design_mat[, 1], 
                    X = design_mat[, 2:6], C = vcv(chronos_tree1), v = NULL, 
                    opt = TRUE)

# Calculate fitted estimates
fitted_test <- lk(sig2 = fit_lik$minimum, y = design_mat[, 1], 
             X = design_mat[, 2:6], C = vcv(chronos_tree1), v = NULL, opt = FALSE)


# Compare these estimates to the gls with only the phylogenetic correlations
fitted_test

pgls_err_test <- 
    gls(slope ~ fg - 1, correlation = corBrownian(1, chronos_tree1), 
        data = mean_spp_summ, method = "ML")
pgls_err_test

#-------------------------------------------------------------------------------
# Now adding the standard errors from the bootstrapping.
#-------------------------------------------------------------------------------
# Calculate design matrix for mean values
design_mat <- model.matrix(~ fg - 1, data = mean_spp_summ)
slopes <- mean_spp_summ$slope
ses <- mean_spp_summ$se

# Optimise for sig2
fit_lik <- optimize(lk, interval = c(0, 1000), y = slopes, 
                    X = design_mat, C = vcv(chronos_tree1), v = ses, 
                    opt = TRUE)

# Calculate fitted estimates
fitted_test <- lk(sig2 = fit_lik$minimum, y = slopes, 
             X = design_mat, C = vcv(chronos_tree1), v = NULL, opt = FALSE)


# Values after 
fitted_test


pgls_err_test <- 
    gls(slope ~ fg - 1, correlation = corBrownian(1, chronos_tree1), 
        data = mean_spp_summ, method = "ML")
pgls_err_test

# What about using just an nlme to approximate and account for family groupings?

# Get dataframe with family groupings


boot_spp_summary

test1 <- lme(Slope ~ 1, random = ~ 1 | factor(fg), data = all_spp_summ_df)
test2 <- lme(Slope ~ fg - 1, random = ~ 1 | factor(fg), data = all_spp_summ_df)
summary(test)
visreg(test)

anova(test)

#-------------------------------------------------------------------------------
# Simply use pgls and a t-test with the bootstrapped SEs
#-------------------------------------------------------------------------------
pgls <- 
    gls(slope ~ fg - 1, correlation = corBrownian(1, chronos_tree1), 
        data = mean_spp_summ, method = "ML")
pgls

lm()

intervals(pgls)

#load("bootSMA_10000.RData")
#boot_summ  <- ldply(bootSMA, mk_boot_summ_df)
#boot_means <- ddply(boot_summ, .(j_fg), get_mean_slope_int_r2)
boot_means

t_test <- function(xbar, mu, se) {
    t <- (xbar - mu) / se
    pval <- pnorm(t) * 2
    return(list(t_stat = t, pval = pval))
}

xbar = coef(pgls)[['fgPi']]
mu = 2
se = boot_means[which(boot_means$j_fg == 'Pi'), 'slope_se']
p_test <- t_test(xbar, mu, se)

xbar = coef(pgls)[['fgBI']]
mu = 2
se = boot_means[which(boot_means$j_fg == 'BI'), 'slope_se']
b_test <- t_test(xbar, mu, se)

xbar = coef(pgls)[['fgHe']]
mu = 2
se = boot_means[which(boot_means$j_fg == 'He'), 'slope_se']
h_test <- t_test(xbar, mu, se)
h_test

xbar = coef(pgls)[['fgZP']]
mu = 2
se = boot_means[which(boot_means$j_fg == 'ZP'), 'slope_se']
zp_test <- t_test(xbar, mu, se)
zp_test

p_test
b_test
h_test
zp_test


