library(plyr)
library(doMC)

source('01_load.r')
source('02_clean.r')

# Bootstrap slopes for each species
# The results from this will be used to calculate a pgls and test whether
# slopes are different from 2, and whether slopes differ from each other.

#-------------------------------------------------------------------------------
# Resample each species
#-------------------------------------------------------------------------------
# Why am I bootstrapping here? And is it the correct thing to do?
# I think that I am boostrapping because I would like to calculate the 
# allometric coefficient with for each species AND the variance of that 
# estimate so that I can use that variance when performing the gls with 
# functional groups. Including this variance increases my power... does it? 
# Is this acceptable to do?

resample_spp <- function(data) {
    len <- length(data[[1]])
    N <- 1:len
    samp <- data[sample(N, size = len, replace = TRUE, prob = NULL), ]
    return(samp)
}

# Get 100 resamples - for now. This should be bumped up to 10 000
n <- 100
spp_resamples <- rlply(.n = n, ddply(pento, .(SpeciesCode), resample_spp), .progress = 'text') 

# Run an SMA for each species
# Learned how to use tryCatch
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
# I don't know whats throwing the errors.
# The warnings are just failure to converge.

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


#-------------------------------------------------------------------------------
# jacknifing vs bootstrap? when do use these techniques
# nah, SO seems to think that bootstrapping is generally better
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Read in the bootstrapped data that was saved from above.
#-------------------------------------------------------------------------------
# write the summary object to a csv so that I don't have to repeat all this
#write.csv(boot_spp_summary, 'species_bootstrapped_coefficients')
# read in this summary object so that I can figure stuff with Liam
boot_spp_summary <- read.csv('species_bootstrapped_coefficients.csv')
#-------------------------------------------------------------------------------


# Assumption that the variance of the bootstrap

# Assumption that the variance of the bootstrap




#hist(boot_spp_summary[which(boot_spp_summary$group == 'AC.NIGR'), 'Slope'])

#-------------------------------------------------------------------------------
# Read in the bootstrapped data that was saved from above.
#-------------------------------------------------------------------------------

library(ape)
library(nlme)
library(multcomp)
library(visreg)

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
spp_lookup_df <- merge(spp_lookup_df, fg_lookup)

all_spp_summ_df <- merge(spp_lookup_df, boot_spp_summary, by.x = 'codes', 
                         by.y = 'group')

mean_spp_summ <- ddply(all_spp_summ_df, .(codes), summarise, 
                       'Slope' = mean(Slope), 
                       'fg' = unique(fg), 
                       'species' = unique(species))
row.names(mean_spp_summ) <- mean_spp_summ$species


#-------------------------------------------------------------------------------
# Prepare phylogenetic information
#-------------------------------------------------------------------------------



# Make tree ultrametric
chronos_tree1 <- chronos(fishTree_no_out, lambda = 1)

# Remove PS.BART for now:
no_psbart <- all_spp_summ_df[which(all_spp_summ_df$codes != 'PS.BART'), ]
no_psbart2 <- all_spp_summ_df[which(mean_spp_summ$codes != 'PS.BART'), ]

# Functional group comparison - not including phylogeny
no_cor_gls <- gls(Slope ~ fg - 1, data = no_psbart, method = "ML")
no_cor_gls <- gls(Slope ~ fg - 1, data = mean_spp_summ, method = "ML")
visreg(no_cor_gls)


# Brownian motion - phylogenetic generalised least squares regression
bm_fish <- corBrownian(phy = chronos_tree1)
bm_gls <- gls(Slope ~ fg - 1, correlation = bm_fish, data = mean_spp_summ, 
              method = "ML")

bm_gls <- gls(Slope ~ fg - 1, correlation = bm_fish, data = no_psbart, 
              method = "ML")

visreg(bm_gls)

anova(no_cor_gls)
anova(bm_gls)

test <- glht(no_cor_gls, linfct = c('fgPi = 2', 'fgBI = 2', 'fgZP = 2', 'fgHe = 2', 
    'fgC = 2'))


test2 <- glht(bm_gls, linfct = c('fgPi = 2', 'fgBI = 2', 'fgZP = 2', 'fgHe = 2', 
    'fgC = 2'))

summary(test)
summary(test2)



# standard estimating equation for least squares regression (multiple linear regression)
XX transpose, XY - cross covariance matrix between x and Y

vector of regression coefficients

m + 1, are the variables 

for anova, make a design matrix where the columns are the categorical variable coding

there are equations that let you get the variance

- generalised least squares ()


# Reproducible example
library(phytools)
library(nlme)

# Create data matrix
spp <- c("Acanthurus_nigricans", "Acanthurus_olivaceus", "Aphareus_furca", 
         "Caranx_melampygus", "Caesio_teres", "Cephalopholis_argus", 
         "Centropyge_flavissima", "Cephalopholis_urodeta", 
         "Chaetodon_ornatissimus", "Chlorurus_sordidus", "Chromis_vanderbilti", 
         "Lutjanus_bohar", "Monotaxis_grandoculis", "Paracirrhites_arcatus",
         "Parupeneus_insularis", "Pseudanthias_dispar", 
         "Pseudanthias_olivaceus", "Pterocaesio_tile", "Scarus_frenatus", 
         "Scarus_rubroviolaceus", "Variola_louti")
slope <- c(2.219673, 2.199110, 1.648774, 1.961538, 1.880211, 1.852918, 2.215627, 
           1.960175, 3.757264, 2.431778, 3.156724, 1.660703, 2.723232, 1.491198, 
           1.947656, 2.206506, 1.927678, 2.536812, 2.576764, 2.796249, 1.669269)
fg <- c('He', 'He', 'Pi', 'Pi', 'ZP', 'Pi', 'He', 'Pi', 'C', 'He', 'ZP', 'Pi', 
        'BI', 'BI', 'BI', 'ZP', 'ZP', 'ZP', 'He', 'He', 'Pi')
data <- data.frame('slope' = slope, 'fg' = fg)
row.names(data) <- spp

design_mat <- model.matrix(~ fg -1, data = data)
design_mat <- cbind(slope, design_mat)

tree <- pbtree(n = 21, tip.label = spp, scale = 1)

# gls
fit.gls <- gls(slope ~ fg - 1, data = data, correlation=corBrownian(1, tree), method="ML")

# custom function
lk <- function (sig2, y, X, C, v=NULL, opt=TRUE) {
    n <- nrow(C)
    if(is.null(v)) v <- rep(0,n)
    V <- sig2 * C + diag(v)
    beta <- solve(t(X) %*% solve(V) %*% X) %*% (t(X) %*% 
        solve(V) %*% y)
    logL <- -(1/2) * t(y - X %*% beta) %*% solve(V) %*% 
        (y - X %*% beta) - (1/2) * determinant(V, 
        logarithm = TRUE)$modulus - (n/2) * log(2 * pi)
    if(opt) -logL else list(beta=beta,sig2e=sig2,logL=logL)
}

fit_lik <- optimize(lk, interval = c(0, 1000), y = design_mat[, 1], 
                    X = design_mat[, 2:6], C = vcv(tree), v = NULL, 
                    opt = TRUE)

fitted <- lk(sig2 = fit_lik$minimum, y = design_mat[, 1], X = design_mat[, 2:6], 
             C = vcv(tree), v = NULL, opt = FALSE)

fit.gls
fitted

# Okay, this reproducible example works... what is going on with my stuff?



bm_fish <- corBrownian(value = 1, phy = chronos_tree1)

# Get best unbiased linear predictor for the sig2 (proportionality constant)
fit_lik <- optimize(lk, interval = c(0, 1000), y = design_mat[, 1], 
                    X = design_mat[, 2:6], C = vcv(chronos_tree1), v = NULL, 
                    opt = TRUE)

# Plug the BLUP for sig2 back into fit_lik to get the fitted values.
fitted <- lk(sig2 = fit_lik$minimum, y = design_mat[, 1], X = design_mat[, 2:6], 
             C = vcv(chronos_tree1), v = NULL, opt = FALSE)
fitted

bm_gls <- gls(Slope ~ fg - 1, correlation = bm_fish, data = mean_spp_summ, 
              method = "ML", x = TRUE)
bm_gls


# create the design matrix for the anova
# each row will be a species, each column will be a functional group
design_mat <- model.matrix(~ fg - 1, data = mean_spp_summ)
design_mat <- cbind(mean_spp_summ$Slope, design_mat)
dimnames(design_mat)[[2]][1] <- 'Slope'

# Liam's function to incorporate intraspecies variation.
lk <- function (sig2, y, X, C, v=NULL, opt=TRUE) {
    n <- nrow(C)
    if(is.null(v)) v <- rep(0,n)
    V <- sig2 * C + diag(v)
    beta <- solve(t(X) %*% solve(V) %*% X) %*% (t(X) %*% 
        solve(V) %*% y)
    logL <- -(1/2) * t(y - X %*% beta) %*% solve(V) %*% 
        (y - X %*% beta) - (1/2) * determinant(V, 
        logarithm = TRUE)$modulus - (n/2) * log(2 * pi)
    if(opt) -logL else list(beta=beta,sig2e=sig2,logL=logL)
}

no_cor_gls <- gls(Slope ~ fg, data = no_psbart, method = "ML")
test <- glht(no_cor_gls, linfct = c('fgPi = 2', 'fgBI = 2', 'fgZP = 2', 'fgHe = 2', 
    'fgC = 2'))

# Brownian motion - phylogenetic generalised least squares regression
bm_fish <- corBrownian(phy = chronos_tree1)

# Get best unbiased linear predictor for the sig2 (proportionality constant)
fit_lik <- optimize(lk, interval = c(0, 1000), y = design_mat[, 1], 
                    X = design_mat[, 2:6], C = vcv(chronos_tree1), v = NULL, 
                    opt = TRUE)

# Plug the BLUP for sig2 back into fit_lik to get the fitted values.
fitted <- lk(sig2 = fit_lik$minimum, y = design_mat[, 1], X = design_mat[, 2:6], 
             C = vcv(chronos_tree1), v = NULL, opt = FALSE)
fitted

bm_gls <- gls(Slope ~ fg - 1, correlation = bm_fish, data = mean_spp_summ, 
              method = "ML", x = TRUE)
bm_gls


library(phytools)
library(nlme)
tree <- pbtree(n=26, tip.label=letters, scale=1)
X <- fastBM(tree,nsim=2)
colnames(X) <- c("x1", "x2")
y <- cbind(rep(1, Ntip(tree)), X) %*% c(1, 2, 3) + fastBM(tree)
fit.gls <- gls(y ~ x1 + x2, data = data.frame(y, X), correlation=corBrownian(1, tree), method="ML")
fit.lk <- optimize(lk, c(0,1000), y=y, X=cbind(rep(1, Ntip(tree)), X), C=vcv(tree))
fitted <- lk(fit.lk$minimum, y=y, X=cbind(rep(1, Ntip(tree)), X), C=vcv(tree),
    opt=FALSE)
v <- setNames(rexp(n=26), tree$tip.label)
fit.lk <- optimize(lk, c(0,1000), y=y, X=cbind(rep(1, Ntip(tree)), X),C=vcv(tree),
    v=v)
fitted <-lk(fit.lk$minimum, y=y, X=cbind(rep(1, Ntip(tree)), X),C=vcv(tree),
    v=v, opt=FALSE)




# Using a gls that accounts for the variance associated with each bootstrapped value
# For now, ignoring the 
no_cor_gls <- gls(Slope ~ fg - 1, data = no_psbart, method = "ML")
test <- glht(no_cor_gls, linfct = c('fgPi = 2', 'fgBI = 2', 'fgZP = 2', 'fgHe = 2', 
    'fgC = 2'))
summary(test)

# Without using bootstrapped values
source('03_func.r')
ga_sma <- function(df) {
  sma(ga~SL, data=df, log="xy", method="SMA", robust=T, slope.test=2)
}

all_sppGA <- dlply(pento, .(SpeciesCode), ga_sma)
all_sppGA_summ <- cbind(attributes(all_sppGA)$split_labels, mk_spp_summary(all_sppGA, 22))
spp_names <- fishTree_no_out$tip.label
spp_codes <- c('AC.NIGR', 'AC.OLIV', 'AP.FURC', 'LU.BOHA', 'CA.TERE', 
               'PT.TILE', 'PS.DISP', 'PS.OLIV', 'CE.ARGU', 'CE.UROD', 
               'VA.LOUT', 'CH.VAND', 'MO.GRAN', 'CH.SORD', 'SC.FREN', 
               'SC.RUBR', 'PA.ARCA', 'CA.MELA', 'PA.INSU', 'CE.FLAV', 
               'CH.ORNA')

# Create lookup dataframes
fg_lookup <- unique(data.frame('codes' = pento$Species, 'fg' = pento$j_fg))
spp_lookup_df <- data.frame('species' = spp_names, 'codes' = spp_codes)
spp_lookup_df <- merge(spp_lookup_df, fg_lookup)

all_sppGA_summ <- merge(spp_lookup_df, all_sppGA_summ, by.x = 'codes', 
                         by.y = 'SpeciesCode')
row.names(all_sppGA_summ) <- all_sppGA_summ$species


# GLS
no_cor_gls <- gls(slope ~ fg - 1, data = all_sppGA_summ, method = "ML")

# PGLS
bm_fish <- corBrownian(phy = chronos_tree1)
bm_gls <- gls(slope ~ fg - 1, correlation = bm_fish, data = all_sppGA_summ, 
              method = "ML")

test <- glht(no_cor_gls, linfct = c('fgPi = 2', 'fgBI = 2', 'fgZP = 2', 'fgHe = 2', 
    'fgC = 2'))
summary(test)

test2 <- glht(bm_gls, linfct = c('fgPi = 2', 'fgBI = 2', 'fgZP = 2', 'fgHe = 2', 
    'fgC = 2'))
summary(test2)

#-------------------------------------------------------------------------------
# Using ga ~ mass

ga_sma <- function(df) {
  sma(ga~wt, data=df, log="xy", method="SMA", robust=T, slope.test=3)
}

all_sppGA <- dlply(pento, .(SpeciesCode), ga_sma)
all_sppGA_summ <- cbind(attributes(all_sppGA)$split_labels, mk_spp_summary(all_sppGA, 22))
spp_names <- fishTree_no_out$tip.label
spp_codes <- c('AC.NIGR', 'AC.OLIV', 'AP.FURC', 'LU.BOHA', 'CA.TERE', 
               'PT.TILE', 'PS.DISP', 'PS.OLIV', 'CE.ARGU', 'CE.UROD', 
               'VA.LOUT', 'CH.VAND', 'MO.GRAN', 'CH.SORD', 'SC.FREN', 
               'SC.RUBR', 'PA.ARCA', 'CA.MELA', 'PA.INSU', 'CE.FLAV', 
               'CH.ORNA')

# Create lookup dataframes
fg_lookup <- unique(data.frame('codes' = pento$Species, 'fg' = pento$j_fg))
spp_lookup_df <- data.frame('species' = spp_names, 'codes' = spp_codes)
spp_lookup_df <- merge(spp_lookup_df, fg_lookup)

all_sppGA_summ <- merge(spp_lookup_df, all_sppGA_summ, by.x = 'codes', 
                         by.y = 'SpeciesCode')
row.names(all_sppGA_summ) <- all_sppGA_summ$species


# GLS
no_cor_gls <- gls(slope ~ fg - 1, data = all_sppGA_summ, method = "ML")

# PGLS
bm_fish <- corBrownian(phy = chronos_tree1)
bm_gls <- gls(slope ~ fg - 1, correlation = bm_fish, data = all_sppGA_summ, 
              method = "ML")

test <- glht(no_cor_gls, linfct = c('fgPi = 0.666', 'fgBI = 0.666', 'fgZP = 0.666', 'fgHe = 0.666', 
    'fgC = 0.666'))
summary(test)

test2 <- glht(bm_gls, linfct = c('fgPi = 0.666', 'fgBI = 0.666', 'fgZP = 0.666', 'fgHe = 0.666', 
    'fgC = 0.666'))
summary(test2)

ga_sma <- function(df) {
  sma(ga~wt, data=df, log="xy", method="SMA", robust=T, slope.test=)
}

ga_wt <- sma(ga~wt*j_fg, data=pento, log="xy", method="SMA", robust=T, slope.test=2/3, multcomp = TRUE, multcompmethod = 'adjusted')


#### Not sure that I need this (August 20, 2015)

library(plyr)

#-------------------------------------------------------------------------------
# Simply use pgls and a t-test with the bootstrapped SEs
#-------------------------------------------------------------------------------
pgls <- 
    gls(slope ~ fg - 1, correlation = corBrownian(1, chronos_tree1), 
        data = mean_spp_summ, method = "ML")
pgls

#load("bootSMA_10000.RData")

boot_summ  <- ldply(bootSMA, mk_boot_summ_df)
boot_means <- ddply(boot_summ, .(j_fg), get_mean_slope_int_r2)
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

rbind.all(list(p_test, b_test, h_test, zp_test))
