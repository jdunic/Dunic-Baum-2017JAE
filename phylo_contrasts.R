# Consideration of phylogenetic non-independence
source('01_load.r')
source('02_clean.r')

library(ape)
library(geiger)
#library(nlme)

# Read in Jon's tree
fishTree <- read.tree('Jillian tree.txt')

# drop outgroup
fishTree_no_out <- drop.tip(fishTree, "OUTGROUP_Myxine_glutinosa")


# Prepare lookup tables to match up the tree and data.
# ----------------------------------------------------------------------------
# get list of species names from the phylogeny (excluding outgroup)
spp_names <- fishTree_no_out$tip.label

# get the species codes that match the names from the phylogeny
spp_codes <- c('AC.NIGR', 'AC.OLIV', 'AP.FURC', 'LU.BOHA', 'CA.TERE', 
               'PT.TILE', 'PS.DISP', 'PS.OLIV', 'CE.ARGU', 'CE.UROD', 
               'VA.LOUT', 'CH.VAND', 'MO.GRAN', 'CH.SORD', 'SC.FREN', 
               'SC.RUBR', 'PA.ARCA', 'CA.MELA', 'PA.INSU', 'CE.FLAV', 
               'CH.ORNA')

spp_lookup_df <- data.frame('species' = spp_names, 'codes' = spp_codes)

# create lookup table for species codes and their functional group 
# classification
fg_lookup <- unique(data.frame('codes' = fish$Species, 'fg' = fish$j_fg))


# lookup table containing speceis codes, functional group, and species names 
# matching the names found in the fishTree
spp_lookup_df <- merge(spp_lookup_df, fg_lookup)


# ----------------------------------------------------------------------------
# Prepare phylogenetic independent contrasts for gape area and standard length
# ----------------------------------------------------------------------------
# Working through: http://phylo.wikidot.com/phylogenetically-independent-contrasts
# Also see: http://www.clarku.edu/faculty/pbergmann/biostats/Biol206%20-%20Lab10-%20PICs.pdf

# get mean gape and body sizes to calculate pics
mean_sizes <- ddply(fish, .(SpeciesCode), summarise, 'mean_gh' = mean(log10(gh)), 
                    'mean_gw' = mean(log10(gw)), 'mean_ga' = mean(log10(ga)), 
                    'mean_sl' = mean(log10(SL)), 'mean_wt' = mean(log10(wt))
              )

# make df to get pics without P. bartlettorum
mean_sizes2 <- mean_sizes[which(mean_sizes$SpeciesCode != 'PS.BART'), ]

# Are our data normal and linearly related?
par(mfrow = c(1, 3))
plot(mean_sizes$mean_gh ~ mean_sizes$mean_sl)
plot(mean_sizes$mean_gw ~ mean_sizes$mean_sl)
plot(mean_sizes$mean_ga ~ mean_sizes$mean_sl)

plot(mean_sizes$mean_gh ~ mean_sizes$mean_wt)
plot(mean_sizes$mean_gw ~ mean_sizes$mean_wt)
plot(mean_sizes$mean_ga ~ mean_sizes$mean_wt)

# Yes they are.
qqnorm(mean_sizes$mean_gh)
qqline(mean_sizes$mean_gh)
qqnorm(mean_sizes$mean_gw)
qqline(mean_sizes$mean_gw)
qqnorm(mean_sizes$mean_ga)
qqline(mean_sizes$mean_ga)

par(mfrow = c(1, 2))
qqnorm(mean_sizes$mean_sl)
qqline(mean_sizes$mean_sl)
qqnorm(mean_sizes$mean_wt)
qqline(mean_sizes$mean_wt)


# break to work through: http://www.clarku.edu/faculty/pbergmann/biostats/Biol206%20-%20Lab10-%20PICs.pdf


# Transform trees to calculate pics
chronos_tree0 <- chronos(fishTree_no_out, lambda = 0)
chronos_tree0.1 <- chronos(fishTree_no_out, lambda = 0.1)
chronos_tree0.5 <- chronos(fishTree_no_out, lambda = 0.5)
chronos_tree0.9 <- chronos(fishTree_no_out, lambda = 0.9)
chronos_tree1 <- chronos(fishTree_no_out, lambda = 1)
chronos_tree_strict <- chronos(fishTree_no_out, lambda = 0, model = 'discrete')
chronos_tree_relaxed <- chronos(fishTree_no_out, lambda = 0, model = 'relaxed')

# Using a log transformation to get properly standardised pics
log_tree <- chronos_tree0
log_tree$edge.length <- log10(chronos_tree0$edge.length + 1)

# Using a sqrt transformation to get properly standardised pics
sqrt_tree <- chronos_tree0
sqrt_tree$edge.length <- sqrt(chronos_tree0$edge.length)


# Use the diagnostic plot script written by David Ackerly:
# http://ib.berkeley.edu/courses/ib200b/scripts/diagnostics_v3.R
source('phylo_diagnostics.R')

# From: http://phylo.wikidot.com/phylogenetically-independent-contrasts
# So, we have our contrasts, but just like with our linear correlations, we need 
# to make sure that our data fit the assumptions of the model.  Luckily David 
# Ackerly has provided us with a function that makes this all very easy. 

# Six plots will appear.  The two on the left will test whether the positivized 
# PICs fit a half-normal distribution centered at 0.  The top left is just a 
# histogram of the PICs and it should look like just half of a normal 
# distribution.  The plot on the bottom right is a QQ-plot of the PICs against a
# half normal.  It should fit approximately a line.  Both those plots look good.
#
# You want the next three plots to show no significant relationship.  The top
# middle is the most important; it plots the standardized contrasts against 
# their branch lengths. If there is a significant negative relationship here, 
# then that implies your branch lengths are bad.  Because contrasts are 
# standardized by dividing by branch lengths, long arbitrary branch lengths can 
# make the contrasts too small.  This problem can be alleviated by replacing all
# your branch lengths with their logs.  There is a negative slope here; it's not
# significant, so we wouldn't normally worry about it, but just for fun let's 
# log convert the branches:

diagnostics(mean_sizes2$mean_gh, chronos_tree0)  # okay
diagnostics(mean_sizes2$mean_gh, chronos_tree1)  # bad
diagnostics(mean_sizes2$mean_gh, chronos_tree_strict)  # bad 
diagnostics(mean_sizes2$mean_gh, fishTree_log)  # v. sig abs val ~ sd contrast
diagnostics(mean_sizes2$mean_gh, sqrt_tree)  # best
diagnostics(mean_sizes2$mean_gh, log_tree) # okay
diagnostics(mean_sizes2$mean_gh, chronos_tree_relaxed) # okay
# Let's use sqrt_tree

diagnostics(mean_sizes2$mean_gw, chronos_tree0)
diagnostics(mean_sizes2$mean_gw, chronos_tree1)
diagnostics(mean_sizes2$mean_gw, chronos_tree_strict)
diagnostics(mean_sizes2$mean_gw, log_tree)
diagnostics(mean_sizes2$mean_gw, sqrt_tree)
# Let's use chronos_tree0

diagnostics(mean_sizes2$mean_ga, chronos_tree0)
diagnostics(mean_sizes2$mean_ga, chronos_tree1)
diagnostics(mean_sizes2$mean_ga, chronos_tree_strict)
diagnostics(mean_sizes2$mean_ga, log_tree)
diagnostics(mean_sizes2$mean_ga, sqrt_tree)
# Let's use chronos_tree0

diagnostics(mean_sizes2$mean_sl, chronos_tree1)
diagnostics(mean_sizes2$mean_sl, log_tree)
diagnostics(mean_sizes2$mean_sl, sqrt_tree)
# Best behaviour comes from sqrt transformed ultrametric tree

diagnostics(mean_sizes2$mean_wt, chronos_tree1)
diagnostics(mean_sizes2$mean_wt, log_tree)
diagnostics(mean_sizes2$mean_sl, sqrt_tree)
# Best behaviour comes from sqrt transformed ultrametric tree


#
mean_sizes2 <- merge(spp_lookup_df, mean_sizes, by.x = 'codes', by.y = 'SpeciesCode')
row.names(mean_sizes2) <- mean_sizes2$species


gape_pics <- apply(mean_sizes2[, c('mean_gh', 'mean_gw', 'mean_ga')], 2, 
                   function(x) { 
                     pic(x, sqrt_tree)
                   })
gape_pics <- as.data.frame(gape_pics)
names(gape_pics) <- c('gh_pics', 'gw_pics', 'ga_pics')

body_pics <- apply(mean_sizes2[, c('mean_sl', 'mean_wt')], 2, 
                   function(x) { 
                     pic(x, sqrt_tree)
                   })
body_pics <- as.data.frame(body_pics)
names(body_pics) <- c('sl_pics', 'wt_pics')

pics <- cbind(mean_sizes2[, c('codes', 'species', 'fg')], gape_pics)

pics <- cbind(mean_sizes2[1], gape_pics)

pics <- merge(rel_len, rel_gape)
pics_n_more <- merge(pics, all_spp_summ_df)





# Include Pseudanthias bartlettorum using phytools add species to genus
# Do it multiple times with PS.BART added randomly to double check that the
# results are robust, but I think that it shouldn't be a problem.
# See Liam's blog post here: http://blog.phytools.org/2013/11/new-function-to-add-species-to-genus-in.html

library(phytools)

# Best to use an ultrametric tree.
# Let's do the same thing for both trees.
chronos_tree1
fishTree_sqrt


species <- 'Pseudanthias_bartlettorum'
tree1 <- add.species.to.genus(chronos_tree1, species)
set.seed(1)
tree_ran1 <- add.species.to.genus(chronos_tree1, species, where = 'random')

> t1<-add.species.to.genus(tree,species)
> plotTree(t1,ftype="i")


# This might just need deleting
#pics_fish <- merge(fish, spp_lookup_df, by.x = 'SpeciesCode', by.y = 'codes')

all_spp_summ_df <- merge(spp_lookup_df, all_spp_summ_df, by.x = 'codes', by.y = 'species')

# Create summary df for pics
pics_fish <- fish[which(fish$SpeciesCode != 'PS.BART'), ]

rel_gape <- ddply(pics_fish, .(SpeciesCode), summarise, 'mean_ga' = mean(log(ga))) 
rel_len  <- ddply(pics_fish, .(SpeciesCode), summarise, 'mean_sl' = mean(log(SL)))

# Making tree ultrametric
chronos_tree0 <- chronos(fishTree_no_out, lambda = 0)
chronos_tree1 <- chronos(fishTree_no_out, lambda = 1)
chronos_tree_strict <- chronos(fishTree_no_out, lambda = 0, model = 'discrete')

pic_log_ga <- pic(rel_gape$mean_ga, chronos_tree0)
pic_log_sl <- pic(rel_len$mean_sl, chronos_tree0)

plot(sma(pic_log_ga ~ pic_log_sl))

cor.test(pic_log_ga, pic_log_sl)




# Working through: http://phylo.wikidot.com/phylogenetically-independent-contrasts
# Also see: http://www.clarku.edu/faculty/pbergmann/biostats/Biol206%20-%20Lab10-%20PICs.pdf
plot(rel_ga$mean_ga, rel_len$mean_sl)
cor.test(rel_gape$mean_ga, rel_len$mean_sl, alternative = 'g')

qqnorm(rel_gape$mean_ga)
qqline(rel_gape$mean_ga)

qqnorm(rel_len$mean_sl)
qqline(rel_len$mean_sl)


diagnostics(rel_gape$mean_ga, chronos_tree0)
diagnostics(rel_gape$mean_ga, chronos_tree1)
diagnostics(rel_gape$mean_ga, chronos_tree_strict)

diagnostics(rel_len$mean_sl, chronos_tree0)


chronos_tree0.1 <- chronos(fishTree_no_out, lambda = 0.1)
chronos_tree0.2 <- chronos(fishTree_no_out, lambda = 0.2)
chronos_tree0.3 <- chronos(fishTree_no_out, lambda = 0.3)
chronos_tree0.4 <- chronos(fishTree_no_out, lambda = 0.4)
chronos_tree0.5 <- chronos(fishTree_no_out, lambda = 0.5)
chronos_tree0.6 <- chronos(fishTree_no_out, lambda = 0.6)
chronos_tree0.7 <- chronos(fishTree_no_out, lambda = 0.7)
chronos_tree0.8 <- chronos(fishTree_no_out, lambda = 0.8)
chronos_tree0.9 <- chronos(fishTree_no_out, lambda = 0.9)
chronos_tree1 <- chronos(fishTree_no_out, lambda = 1)
chronos_tree_strict <- chronos(fishTree_no_out, lambda = 0, model = 'discrete')


fishTree_log = chronos_tree1
fishTree_log$edge.length = log(chronos_tree1$edge.length + 1)

fishTree_sqrt = chronos_tree1
fishTree_sqrt$edge.length = sqrt(chronos_tree1$edge.length + 1)

pic_log_sl0.1 <- pic(rel_len$mean_sl, chronos_tree0.1)
pic_log_sl1 <- pic(rel_len$mean_sl, chronos_tree1)
pic_sqrt_sl <- pic(rel_len$mean_sl, fishTree_sqrt)

diagnostics(rel_len$mean_sl, chronos_tree0.1)
diagnostics(rel_len$mean_sl, chronos_tree0.2)
diagnostics(rel_len$mean_sl, chronos_tree0.3)
diagnostics(rel_len$mean_sl, chronos_tree0.4)
diagnostics(rel_len$mean_sl, chronos_tree0.5)
diagnostics(rel_len$mean_sl, chronos_tree0.6)
diagnostics(rel_len$mean_sl, chronos_tree0.7)
diagnostics(rel_len$mean_sl, chronos_tree0.8)
diagnostics(rel_len$mean_sl, chronos_tree0.9)
diagnostics(rel_len$mean_sl, chronos_tree1)
diagnostics(rel_len$mean_sl, chronos_tree_strict)
diagnostics(rel_len$mean_sl, fishTree_log)
diagnostics(rel_len$mean_sl, fishTree_sqrt)

pics <- merge(rel_len, rel_gape)
pics_n_more <- merge(pics, all_spp_summ_df)



pics_sma <- sma(mean_ga ~ mean_sl * fg, data = pics_n_more, log = "", method = "SMA", slope.test = 2)


pics_sma <- sma(mean_ga ~ mean_sl * fg - 1, data = pics_n_more, log = "", method = "SMA", slope.test = 2)

all_sma <- sma(ga ~ SL * j_fg, data = fish, log = "xy", method = "SMA", robust = T, slope.test = 2)


# Species level analysis for phylogenetic generalized least squares regression
library(ape)
library(nlme)
library(visreg)
library(multcomp)

source('03_func.r')

# Set up bootstrap values for comparison once we've done the pgls.
load('bootSMA_10000.RData')
slp_check_df <- data.frame(j_fg = c("Pi", "BI", "ZP", "He", "C"), iso = 0, 
                           neg = 0, pos = 0)

# Create dataframe to determine whether sma regression was neg, iso, or pos
# Checking allometry
check_df <- ldply(bootSMA, function(x) slp_check(x))
allo_check_df <- ddply(check_df, .(j_fg), summarise, 
             iso = sum(iso), neg = sum(neg), pos = sum(pos)
             )
allo_check_df$j_fg <- factor(allo_check_df$j_fg, levels = c("Pi", "BI", "ZP", "He", "C"))
allo_check_df <- arrange(allo_check_df, j_fg)

#-------------------------------------------------------------------------------
# Multiple Comparison Check
# For each bootSMA and the multiple pairwise comparison between functional groups
# count if there was a significant difference
multcomp_all_df <- ldply(bootSMA, function(x) multcomp_check(x))
multcomp_df <- ddply(multcomp_all_df, .(j_fg_1, j_fg_2), summarise, sig_diff = sum(sig_diff))
multcomp_df$prop <- multcomp_df$sig_diff / 10000
multcomp_df

#-------------------------------------------------------------------------------
# Average functional group slopes, intercepts, and r2
boot_summ  <- ldply(bootSMA, mk_boot_summ_df)
boot_means <- ddply(boot_summ, .(j_fg), get_mean_slope_int_r2)
boot_means

# boot_means is what I will ultimately compare the results of the pgls with.


fishTree <- read.tree('Jillian tree.txt')
# drop outgroup
fishTree_no_out <- drop.tip(fishTree, "OUTGROUP_Myxine_glutinosa")

ga_sma <- function(df) {
  sma(ga~SL, data=df, log="xy", method="SMA", robust=T, slope.test=2)
}

all_sppGA <- dlply(fish, .(SpeciesCode), ga_sma)
all_sppGA_summ <- cbind(attributes(all_sppGA)$split_labels, mk_spp_summary(all_sppGA, 22))


# Need to deail with PS.BART which is missing from the phylogeny. 
spp_names <- fishTree_no_out$tip.label
spp_codes <- c('AC.NIGR', 'AC.OLIV', 'AP.FURC', 'LU.BOHA', 'CA.TERE', 
               'PT.TILE', 'PS.DISP', 'PS.OLIV', 'CE.ARGU', 'CE.UROD', 
               'VA.LOUT', 'CH.VAND', 'MO.GRAN', 'CH.SORD', 'SC.FREN', 
               'SC.RUBR', 'PA.ARCA', 'CA.MELA', 'PA.INSU', 'CE.FLAV', 
               'CH.ORNA')

# Create lookup dataframes
fg_lookup <- unique(data.frame('codes' = fish$Species, 'fg' = fish$j_fg))
spp_lookup_df <- data.frame('species' = spp_names, 'codes' = spp_codes)
spp_lookup_df <- merge(spp_lookup_df, fg_lookup)

all_spp_summ_df <- merge(spp_lookup_df, all_sppGA_summ, by.x = 'codes', by.y = 'SpeciesCode')

allom_coef <- all_spp_summ_df$slope
fg <- all_spp_summ_df$fg

fish_df <- data.frame(allom_coef, fg, row.names = all_spp_summ_df$species)
fish_df[fishTree_no_out$tip.label, ]
fish_df


no_cor_gls <- gls(allom_coef ~ fg - 1, data = fish_df[which(fish_df$fg != 'C'), ], method = "ML")

# Brownian motion:
# First make tree ultrametric
chronos_tree1 <- chronos(fishTree_no_out, lambda = 1)
bm_fish <- corBrownian(phy = fishTree_no_out)
bm_fish <- corBrownian(phy = chronos_tree1)
bm_gls <- gls(allom_coef ~ fg - 1, correlation = bm_fish, data = fish_df, method = "ML")

# Now to see whether this phylogenetically-corrected model supports the same 
# results that we found using the SMA, non-phylogenetically corrected model.
summary(no_cor_gls)
summary(bm_gls)

# Let's look back at boot_means and see if we still have the same differences
boot_no_cor_gls <- gls(slope ~ j_fg - 1, data = boot_means, method = "ML")
boot_means

summary(no_cor_gls)
anova1 <- anova(no_cor_gls)
anova1
tukey1 <- glht(no_cor_gls, lincfct = mcp(fg = "Tukey"))
tukey1
summary(bm_gls)
anova2 <- anova(bm_gls)
anova2
tukey2 <- glht(bm_gls, lincfct = mcp(fg = "Tukey"))
tukey2

postHoc <- TukeyHSD(anova1)

par(mfrow = c(1, 2))
visreg(no_cor_gls)
visreg(bm_gls)









# Using phylogenetically independent contrasts
# http://www.eve.ucdavis.edu/~wainwrightlab/Intro2Phylo_S5.R

## phylogenetic independent contrasts:

# When we have species-level data on any set of traits, we are often interested in the correlation between traits. such correlation (given a model) can help us understand the causes and effects. for example, we can test the correlation between body size and brain size, assuming that larger body size provides the resources needed to evolve large brains. However, species are not independent of each other, and this dependence (summeriesed by the phylogeny) is validating one of the core assumptions of regression/correlation. Felsenstein (1985. Phylogenies and the comparative method. American Naturalist 125: 1â€“15) proposed a method that fully takes phylogeny into account in the analysis of comparative data. The idea behind the â€œcontrastsâ€ method is that, if we assume that a continuous trait evolves randomly in any direction (i.e., the Brownian motion model), then the â€œcontrastâ€ between two species is expected to have a distribution centered on zero and a variance proportional to the time since divergence. If the contrasts are scaled with the latter, then they have a variance equal to one. This can be done also for internal nodes because under the assumptions of the Brownian model the ancestral state of the variable can be calculated; a rescaling of the internal branches eventually occurs. 
#In this formulation, the tree needs to be binary (fully dichotomous), and a contrast is computed for each node. Thus for n species, n âˆ’ 1 contrasts will be computed. The contrasts are independent with respect to the phylogeny (unlike the original values of x), and standard statistical methods for continuous variables can be used.

pic.log_gape=pic(mysorteddata$log_gape,mytree)
#pic's are calculated for each trait separately. they are point estimates for the rate of change in trait values at the node, so the mean of pic.log_gape should be similar to the brownian rate estimated from the fitContinous model
mean(pic.log_gape)
log_gape.BM$Trait1$beta

#pics are used in comparative studies, typically when comparing two traits in a regression-like analysis. let's add a second variable
pic.log_protrusion=pic(mysorteddata$log_protrusion,mytree)

# Garland et al. (1993, Phylogenetic analysis of covariance by computer simulation. Systematic Biology 42: 265â€“292) recommended that linear regressions with PICs should be done through the origin (i.e. the intercept is set to zero). this can be done using linear regression model (GLM) or using major axis regression using the smatr module in R

# linear regression analysis of the raw data. slope is NOT forced through the origin. the order of variables is important as residuals are calculated on the Y dimension only
lm.raw=lm(mysorteddata$log_gape~mysorteddata$log_protrusion) 
summary(lm.raw)

# linear regression analysis of the PICs. slope is forced through the origin (by adding -1 in the formula). the order of variables is important as residuals are calculated on the Y dimension only
lm.pic=lm(pic.log_gape~pic.log_protrusion-1) 
summary(lm.pic)

#an alternative test is the major axis regression, recommended when there is no cause and effect. residuals are calculated based on X and Y distance from the regression line. intercept is forced through the origin (intercept=F). Y variable comes first. the method by which regression line is fitted can be 'OLS' or 0 for linear regression, 'SMA' or 1 for standardized major axis (this is the default), and 'MA' or 2 for major axis

library(smatr)
ma.pic=slope.test(pic.log_protrusion,pic.log_gape, test.value = 0, intercept = F,method = 2 )
ma.pic

# positive correlation between the variables would be interpreted as a correlated evolution between the two traits such that evolutionary increases in one trait are associated with evolutionary increases in another. here, we see no such correlation, indicating that the observed correlation is an artifact of the phylogeny. however often with contrasts, very short branches turn out to be outliers. you can try to remove contrasts for very short branches
short.pic.gape=pic.log_gape[-which((pic.log_gape)==min(pic.log_gape))]
short.pic.prot=pic.log_protrusion[-which((pic.log_gape)==min(pic.log_gape))]
ma.pic.short=slope.test(short.pic.prot,short.pic.gape, test.value = 0, intercept = F,method = 2 )

par(mfrow=c(2,2))#remember from day1 session2 that this sets the graphical parameters so that the plotting device has 2 row and 2 columns, so we can now plot the different contrast plots.
plot(pic.log_gape~pic.log_protrusion) 
#we can plot the regression line from the lm (full line)
abline(a=0,b=summary(lm.pic)$coefficients[1])
#and plot the regression line from the MA (dashed line)
abline(a=0,b=ma.pic$b,lty = 2)

plot(mysorteddata$log_gape~mysorteddata$log_protrusion)
 #we can plot the regression line
abline(a=summary(lm.raw)$coefficients[1,1], b=summary(lm.raw)$coefficients[2,1])

plot(short.pic.prot,short.pic.gape)
 #we can plot the regression line
abline(a=0, b=ma.pic.short$b)

# Making tree ultrametric
chronos_tree0 <- chronos(fishTree_no_out, lambda = 0)
chronos_tree1 <- chronos(fishTree_no_out, lambda = 0)
chronos_tree_strict <- chronos(fishTree_no_out, lambda = 0, model = 'discrete')

require("ape")
tr = rtree(10) # generate a random, non-ultrametric tree
pdf("clock-trees_from_chronos.pdf")
par(mfcol = c(4,2))
plot(tr, main = "original tree") # original tree
plot(chronos(tr, lambda = 0, model = "correlated"), main = "correlated, l=0") # correlated clock model, lambda parameter set to 0
axisPhylo()
plot(chronos(tr, lambda = 0.1, model = "correlated"), main = "correlated, l=0.1") # correlated clock model, lambda parameter set to 0.1
axisPhylo()
plot(chronos(tr, lambda = 1.0, model = "correlated"), main = "correlated, l=1.0") # correlated clock model, lambda parameter set to 1.0
axisPhylo()
plot(chronos(tr, model = "discrete", control = chronos.control(nb.rate.cat = 1)), main = "strict") # strict clock model
axisPhylo()
plot(chronos(tr, model = "relaxed", lambda = 0), main = "relaxed, l=0") # relaxed clock model, lambda parameter set to 0
axisPhylo()
plot(chronos(tr, model = "relaxed", lambda = 0.1), main = "relaxed, l=0.1") # relaxed clock model, lambda parameter set to 0.1
axisPhylo()
plot(chronos(tr, model = "relaxed", lambda = 1.0), main = "relaxed, l=1.0") # relaxed clock model, lambda parameter set to 1.0
axisPhylo()
dev.off()


########
# Trying the PICs with smatr, doing it for each group
########
# Okay, here I have exhausted the attempt to use PICs + smatr. We have too few 
# data.
########

drop.tip(sqrt_tree, sqrt_tree$tip.label[!(sqrt_tree$tip.label %in% row.names(mean_sizes2[which(mean_sizes2$fg == 'Pi'), ]))])

library(sma)

row.names(mean_sizes2) <- mean_sizes2$species

# calculate pics for each species
get_pic <- function(data, x, y, tree, xpic_name, ypic_name) {
  #browser()
  tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% data$species)])
  x_pic <- pic(data[[x]], tree)
  y_pic <- pic(data[[y]], tree)
  fg <- data$fg[-1]
  ret <- data.frame(fg = fg, xpic = x_pic, ypic = y_pic)
  names(ret) <- c('fg', xpic_name, ypic_name)
  return(ret)
}

# ignoring corallivore
fg_pic_list <- dlply(mean_sizes2[-which(mean_sizes2$fg == 'C'), ], .(fg), 
                function(x) get_pic(data = x, x = 'mean_sl', y = 'mean_ga', tree = sqrt_tree, 'sl_pic', 'ga_pic'))

fg_pic_df <- rbind.fill(fg_pic_list)

test <- sma(ga_pic ~ sl_pic * fg, data = fg_pic_df, log = "xy", method = "SMA", 
            robust = T, slope.test = 2)


ga_sma <- function(df) {
  sma(ga_pic ~ sl_pic, data=df, log="xy", method="SMA", robust=T, slope.test=2)
}

ga_sma(df = fg_pic_list$Pi)

pi_pics <- pic(mysorteddata$log_gape, sqrt_tree)

pi_tree <- drop.tip(sqrt_tree, sqrt_tree$tip.label[!(sqrt_tree$tip.label %in% mean_sizes2[which(mean_sizes2$fg == 'Pi'), 'species'])])

pic(mean_sizes2[which(mean_sizes2$fg == 'Pi'), 'mean_sl'], pi_tree)

pi_GA_SL_sma <- ga_sma(df = mean_sizes2[which(mean_sizes2$fg == 'Pi'), ])
bi_GA_SL_sma <- ga_sma(df = mean_sizes2[which(mean_sizes2$fg == 'BI'), ])
he_GA_SL_sma <- ga_sma(df = mean_sizes2[which(mean_sizes2$fg == 'He'), ])
zp_GA_SL_sma <- ga_sma(df = mean_sizes2[which(mean_sizes2$fg == 'ZP'), ])