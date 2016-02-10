if (getwd() != '/Users/jillian/R_projects/Allometry') setwd('Allometry')

library(ape)
library(nlme)
library(multcomp)
library(beepr)
library(visreg)
library(dplyr)

# Prepare phylogenetic tree and bootstrapped species data for analysis
#-------------------------------------------------------------------------------
boot_spp_summary <- read.csv('species_bootstrapped_coefficients_10000.csv')

boot_spp_summary_gh <- 
  read.csv('gape_height_species_bootstrapped_coefficients_10000.csv')

boot_spp_summary_gw <- 
  read.csv('gape_width_species_bootstrapped_coefficients_10000.csv')

# Read in Jon's tree
fishTree <- read.tree('Jillian tree.txt')

# drop outgroup
fishTree_no_out <- drop.tip(fishTree, "OUTGROUP_Myxine_glutinosa")
fishTree_no_out <- drop.tip(fishTree_no_out, "Chaetodon_ornatissimus")

spp_names <- fishTree_no_out$tip.label
spp_codes <- c('AC.NIGR', 'AC.OLIV', 'AP.FURC', 'LU.BOHA', 'CA.TERE', 
               'PT.TILE', 'PS.DISP', 'PS.OLIV', 'CE.ARGU', 'CE.UROD', 
               'VA.LOUT', 'CH.VAND', 'MO.GRAN', 'CH.SORD', 'SC.FREN', 
               'SC.RUBR', 'PA.ARCA', 'CA.MELA', 'PA.INSU', 'CE.FLAV')#, 
               #'CH.ORNA')

# Create lookup dataframes
fg_lookup <- unique(data.frame('codes' = pento$Species, 'fg' = pento$j_fg))
spp_lookup_df <- data.frame('species' = spp_names, 'codes' = spp_codes)
fam_lookup <- ddply(fish, .(SpeciesCode), summarise, 'Family' = unique(Family), 
                    'Order' = unique(Order))
spp_lookup_df <- merge(spp_lookup_df, fg_lookup)

# Gape area
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

# Gape height
all_spp_summ_df_gh <- merge(fam_lookup, boot_spp_summary_gh, 
                         by.x = 'SpeciesCode', by.y = 'group')
all_spp_summ_df_gh <- merge(spp_lookup_df, all_spp_summ_df_gh, by.x = 'codes', 
                         by.y = 'SpeciesCode')

mean_spp_summ_gh <- ddply(all_spp_summ_df_gh, .(codes), summarise, 
                       'slope' = mean(Slope), 
                       'fg' = unique(fg), 
                       'species' = unique(species), 
                       'se' = sd(Slope))
row.names(mean_spp_summ_gh) <- mean_spp_summ_gh$species

# Gape height
all_spp_summ_df_gw <- merge(fam_lookup, boot_spp_summary_gw, 
                         by.x = 'SpeciesCode', by.y = 'group')
all_spp_summ_df_gw <- merge(spp_lookup_df, all_spp_summ_df_gw, by.x = 'codes', 
                         by.y = 'SpeciesCode')

mean_spp_summ_gw <- ddply(all_spp_summ_df_gw, .(codes), summarise, 
                       'slope' = mean(Slope), 
                       'fg' = unique(fg), 
                       'species' = unique(species), 
                       'se' = sd(Slope))
row.names(mean_spp_summ_gw) <- mean_spp_summ_gw$species

#-------------------------------------------------------------------------------
# Prepare phylogenetic information (excluding PS.BART for now)
#-------------------------------------------------------------------------------

# Make tree ultrametric
chronos_tree1 <- chronos(fishTree_no_out, lambda = 1)
#chronos_tree0 <- chronos(fishTree_no_out, lambda = 0)

#-------------------------------------------------------------------------------
# Fit pgls and account for intraspecies variation
#-------------------------------------------------------------------------------

# Need to filter out the CH.ORNA data otherwise the calculation of the likelihood
# solution fails to converge. Must be removed from both the mean_spp_summ and 
# the design matrix.
mean_spp_summ <- dplyr::filter(mean_spp_summ, codes != 'CH.ORNA')
# Gape area
design_mat <- model.matrix(~ fg - 1, data = dplyr::filter(mean_spp_summ, codes != 'CH.ORNA'))[, -5]
slopes <- setNames(mean_spp_summ$slope, mean_spp_summ$species)
ses <- setNames(mean_spp_summ$se, mean_spp_summ$species)

# Gape height
mean_spp_summ_gh <- dplyr::filter(mean_spp_summ_gh, codes != 'CH.ORNA')
design_mat_gh <- model.matrix(~ fg - 1, data = dplyr::filter(mean_spp_summ_gh, codes != 'CH.ORNA'))[, -5]
slopes_gh <- setNames(mean_spp_summ_gh$slope, mean_spp_summ_gh$species)
ses_gh <- setNames(mean_spp_summ_gh$se, mean_spp_summ_gh$species)

# Gape height
mean_spp_summ_gw <- dplyr::filter(mean_spp_summ_gw, codes != 'CH.ORNA')
design_mat_gw <- model.matrix(~ fg - 1, data = dplyr::filter(mean_spp_summ_gw, codes != 'CH.ORNA'))[, -5]
slopes_gw <- setNames(mean_spp_summ_gw$slope, mean_spp_summ_gw$species)
ses_gw <- setNames(mean_spp_summ_gw$se, mean_spp_summ_gw$species)


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

# Gape area fit model
fit.lk <- optimize(lk, c(0, 1000), y = slopes, X = design_mat, C = vcv(chronos_tree1), v = ses)
fitted <- lk(fit.lk$minimum, y = slopes, X = design_mat, C = vcv(chronos_tree1), v = ses, opt = FALSE)
fitted

# Gape height fit model
fit.lk_gh <- optimize(lk, c(0, 1000), y = slopes_gh, X = design_mat_gh, C = vcv(chronos_tree1), v = ses_gh)
fitted_gh <- lk(fit.lk_gh$minimum, y = slopes_gh, X = design_mat_gh, C = vcv(chronos_tree1), v = ses_gh, opt = FALSE)
fitted_gh

# Gape width fit model
fit.lk_gw <- optimize(lk, c(0, 1000), y = slopes_gw, X = design_mat_gw, C = vcv(chronos_tree1), v = ses_gw)
fitted_gw <- lk(fit.lk_gw$minimum, y = slopes_gw, X = design_mat_gw, C = vcv(chronos_tree1), v = ses_gw, opt = FALSE)
fitted_gw



# From Liam's blog post: 
# http://blog.phytools.org/2015/05/pgls-with-measurement-or-sampling-error.html
# "Now, it also turns out that after we have optimized sig2e we can actually 
# coerce gls into giving us the correct fitted model & likelihood. Here, I do 
# this by distorting the edge lengths of our tree to take into account the 
# fitted sig2e and within-species errors in y."

# Gape area Coerce tree
tt <- chronos_tree1
tt$edge.length <- tt$edge.length * fitted$sig2e
for(i in 1:length(ses)){
    tip <- which(tt$tip.label == names(ses)[i])
    ii <- which(tt$edge[,2] == tip)
    tt$edge.length[ii] <- tt$edge.length[ii] + ses[i]
}
vv <- diag(vcv(tt))
w <- varFixed(~vv)
fit.gls <- gls(slopes ~ fg - 1, data = dplyr::filter(mean_spp_summ, codes != 'CH.ORNA'), correlation = corBrownian(1, tt), method = "ML", weights = w)
fit.gls
summary(fit.gls)

# Gape height Coerce tree
tt <- chronos_tree1
tt$edge.length <- tt$edge.length * fitted$sig2e
for(i in 1:length(ses_gh)){
    tip <- which(tt$tip.label == names(ses_gh)[i])
    ii <- which(tt$edge[,2] == tip)
    tt$edge.length[ii] <- tt$edge.length[ii] + ses_gh[i]
}
vv <- diag(vcv(tt))
w <- varFixed(~vv)
fit.gls_gh <- gls(slopes_gh ~ fg - 1, data = mean_spp_summ_gh, correlation = corBrownian(1, tt), method = "ML", weights = w)
fit.gls_gh
summary(fit.gls_gh)

# Gape width Coerce tree
tt <- chronos_tree1
tt$edge.length <- tt$edge.length * fitted$sig2e
for(i in 1:length(ses_gw)){
    tip <- which(tt$tip.label == names(ses_gw)[i])
    ii <- which(tt$edge[,2] == tip)
    tt$edge.length[ii] <- tt$edge.length[ii] + ses_gw[i]
}
vv <- diag(vcv(tt))
w <- varFixed(~vv)
fit.gls_gw <- gls(slopes_gw ~ fg - 1, data = mean_spp_summ_gw, correlation = corBrownian(1, tt), method = "ML", weights = w)
fit.gls_gw
summary(fit.gls_gw)


# Test for isometry
#-------------------------------------------------------------------------------
spp_err_iso <- glht(fit.gls, linfct = c('fgPi = 2', 
                                        'fgBI = 2', 
                                        'fgZP = 2', 
                                        'fgHe = 2'))#, 
                                        #'fgC  = 2'))
summary(spp_err_iso)

spp_err_iso_gh <- glht(fit.gls_gh, linfct = c('fgPi = 1', 
                                              'fgBI = 1', 
                                              'fgZP = 1', 
                                              'fgHe = 1'))#, 
                                              #'fgC  = 1'))
summary(spp_err_iso_gh)

spp_err_iso_gw <- glht(fit.gls_gw, linfct = c('fgPi = 1', 
                                              'fgBI = 1', 
                                              'fgZP = 1', 
                                              'fgHe = 1'))#, 
                                              #'fgC  = 1'))
summary(spp_err_iso_gw)


# Multiple comparisons
#-------------------------------------------------------------------------------
# Functions from: 
# http://rstudio-pubs-static.s3.amazonaws.com/13472_0daab9a778f24d3dbf38d808952455ce.html

#model.matrix.gls <- function(object, ...) {
#    model.matrix(terms(object), data = getData(object), ...)
#}
#model.frame.gls <- function(object, ...) {
#    model.frame(formula(object), data = getData(object), ...)
#}
#terms.gls <- function(object, ...) {
#    terms(model.frame(object), ...)
#}

# The line below gives pairwise comparisons now.  Note that the above
# performs t-tests for all pairwise differences.
#multCompTukey <- glht(fit.gls, linfct = mcp(fg = "Tukey"))

#multCompTukey <- glht(fit.gls_gh, mcp(fg = "Tukey"))

#summary(multCompTukey)


#multCompTukey <- glht(fit.gls_gw, mcp(fg = "Tukey"))

#summary(multCompTukey)


#-------------------------------------------------------------------------------
# Double check that our modified pgls with y-error matches under the test case
# that the y-error = 0.
#-------------------------------------------------------------------------------
#design_mat <- model.matrix(~ fg - 1, data = mean_spp_summ)
#slopes <- setNames(mean_spp_summ$slope, mean_spp_summ$species)
#ses <- setNames(mean_spp_summ$se, mean_spp_summ$species)
#ses <- setNames(rep(0, 21), mean_spp_summ$species)

#fit.lk <- optimize(lk, c(0, 1000),y = slopes, X = design_mat, C = vcv(chronos_tree1), v = NULL)
#fitted <- lk(fit.lk$minimum, y = slopes, X = design_mat, C = vcv(chronos_tree1), v = NULL, opt = FALSE)
#fitted

# Coerce tree:
#tt <- chronos_tree1
#tt$edge.length <- tt$edge.length * fitted$sig2e
#for(i in 1:length(ses)){
#    tip <- which(tt$tip.label == names(ses)[i])
#    ii <- which(tt$edge[,2] == tip)
#    tt$edge.length[ii] <- tt$edge.length[ii] + ses[i]
#}
#vv <- diag(vcv(tt))
#w <- varFixed(~vv)

#fit.gls <- gls(slopes ~ fg - 1, data = mean_spp_summ, correlation = corBrownian(1, tt), method = "ML", weights = w)
#fit.gls
#summary(fit.gls)

# Regular pgls using the chronos_tree1
#fit.gls.test <- gls(slopes ~ fg - 1, data = mean_spp_summ, correlation = corBrownian(1, chronos_tree1), method = "ML")
#fit.gls.test

# Woohoo!
#all.equal(coef(fit.gls), coef(fit.gls.test))
#all.equal(fit.gls$varBeta, fit.gls.test$varBeta)

# Slope comparison plots for MS
fit.gls$data <- dplyr::filter(mean_spp_summ, codes != 'CH.ORNA')
v <- visreg(fit.gls)

v$fit$visregLwr <- as.data.frame(confint(fit.gls, param.CI=0.95))[[1]]
v$fit$visregUpr <- as.data.frame(confint(fit.gls, param.CI=0.95))[[2]]


# Gape height
fit.gls_gh$data <- dplyr::filter(mean_spp_summ_gh, codes != 'CH.ORNA')
vgh <- visreg(fit.gls_gh)

vgh$fit$visregLwr <- as.data.frame(confint(fit.gls_gh, param.CI=0.95))[[1]]
vgh$fit$visregUpr <- as.data.frame(confint(fit.gls_gh, param.CI=0.95))[[2]]

# Gape width
fit.gls_gw$data <- dplyr::filter(mean_spp_summ_gw, codes != 'CH.ORNA')

vgw <- visreg(fit.gls_gw)

vgw$fit$visregLwr <- as.data.frame(confint(fit.gls_gw, param.CI=0.95))[[1]]
vgw$fit$visregUpr <- as.data.frame(confint(fit.gls_gw, param.CI=0.95))[[2]]

#------------------------------------------------------------------------------
# Phylogenetically corrected comparisons of slopes
#------------------------------------------------------------------------------
mean_spp_summ_gh2 <- left_join(mean_spp_summ_gh, vgh$fit) %>% 
  mutate(id = as.numeric(codes)) %>% 
  mutate(fg = factor(fg, levels = levels(mean_spp_summ_gh$fg))) %>% 
  arrange(., fg)

spacing <- c(seq(0, 1, length.out = 6), seq(0, 1, length.out = 3), 
             seq(0, 1, length.out = 5), seq(0, 1, length.out = 6))
mean_spp_summ_gh2$spacing <- spacing

line_fits_gh <- 
  mean_spp_summ_gh2 %>% group_by(fg, visregLwr, visregUpr, visregFit) %>% summarise(count = length(fg)) %>% as.data.frame

gh_plot <- 
ggplot(mean_spp_summ_gh2) +
  geom_rect(data = line_fits_gh, aes(xmin = -0.05, xmax = 1.05, ymin = visregLwr, ymax = visregUpr), fill = "gray85") +
  geom_point(aes(x = -0.2, y = 1), alpha = 0) +
  geom_point(aes(x = spacing, y = slope), size = 1, colour = "gray50") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.margin = unit(0, "cm")) +
  theme(panel.border = element_blank()) +
  theme(plot.margin = unit(c(1.5,0.5,0.5,0.5), "lines")) +
  theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(axis.line = element_line()) +
  theme(strip.text.x = element_blank()) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(margin = margin(0, 6, 0, 0)), 
        #axis.ticks.x=element_blank(), 
        axis.title.x=element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks.length = unit(-0.1, "cm")) +
  geom_segment(data=line_fits_gh, aes(y = visregFit, yend = visregFit, x = -0.05, xend = 1.05), colour = '#099DFFFF', size = 1, lineend = 'round') +
  geom_hline(aes(yintercept = 1), colour = 'red', linetype = 2) +
  scale_x_continuous(breaks = 0.5) +
  scale_y_continuous(breaks = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), lim =c(0.6, 2)) +
  facet_grid(. ~ fg) 

mean_spp_summ_gw2 <- left_join(mean_spp_summ_gw, vgw$fit) %>% 
  mutate(id = as.numeric(codes)) %>% 
  mutate(fg = factor(fg, levels = levels(mean_spp_summ_gw$fg))) %>% 
  arrange(., fg)

spacing <- c(seq(0, 1, length.out = 6), seq(0, 1, length.out = 3), 
             seq(0, 1, length.out = 5), seq(0, 1, length.out = 6))
mean_spp_summ_gw2$spacing <- spacing

line_fits_gw <- 
  mean_spp_summ_gw2 %>% group_by(fg, visregLwr, visregUpr, visregFit) %>% summarise(count = length(fg)) %>% as.data.frame

gw_plot <- 
ggplot(mean_spp_summ_gw2) +
  geom_rect(data = line_fits_gw, aes(xmin = -0.05, xmax = 1.05, ymin = visregLwr, ymax = visregUpr), fill = "gray85") +
  geom_point(aes(x = -0.2, y = 1), alpha = 0) +
  geom_point(aes(x = spacing, y = slope), size = 1, colour = "gray50") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.margin = unit(0, "cm")) +
  theme(panel.border = element_blank()) +
  theme(plot.margin = unit(c(1.5,0.5,0.5,0.5), "lines")) +
  theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(axis.line = element_line()) +
  theme(strip.text.x = element_blank()) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(margin = margin(0, 6, 0, 0)), 
        #axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.length = unit(-0.1, "cm")) +
  geom_segment(data=line_fits_gw, aes(y = visregFit, yend = visregFit, x = -0.05, xend = 1.05), colour = '#099DFFFF', size = 1, lineend = 'round') +
  geom_hline(aes(yintercept = 1), colour = 'red', linetype = 2) +
  scale_x_continuous(breaks = 0.5) +
  scale_y_continuous(breaks = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), lim =c(0.6, 2)) +
  #ylim(c(0.6, 2)) +
  facet_grid(. ~ fg) 

dev.new(height = 6.6, width = 4.8)

master_layout <- 
grid.layout(nrow = 4, ncol = 6, 
            widths = unit(c(0.2, 0.3, 1, 1, 1, 1), "null"),
            heights = unit(c(1, 1, 0.2, 0.1), "null"))
grid.newpage()
pushViewport(viewport(layout=master_layout))
top <- viewport(layout.pos.row = 1, layout.pos.col = 2:6)
bottom <- viewport(layout.pos.row = 2, layout.pos.col = 2:6)
print(gh_plot, vp = top)
print(gw_plot, vp = bottom)
grid.text(
  "a)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1), 
  gp = gpar(fontsize = 9), vjust = -12
  )
grid.text(
  "b)", vp = viewport(layout.pos.row = 2, layout.pos.col = 1), 
  gp = gpar(fontsize = 9), vjust = -12
  )
grid.text(
  "Gape height slope estimate",
  vp = viewport(layout.pos.row = 1, layout.pos.col = 1),
  rot = 90, gp = gpar(fontsize = 9), vjust = 1, hjust = 0.6
    )
grid.text(
  "Gape width slope estimate",
  vp = viewport(layout.pos.row = 2, layout.pos.col = 1),
  rot = 90, gp = gpar(fontsize = 9), vjust = 1, hjust = 0.6
    )
grid.text(
    "Piscivore",
    vp = viewport(layout.pos.row = 3, layout.pos.col = 3),
    vjust = -1.2, gp = gpar(fontsize = 9), hjust = 0.45
    )
grid.text(
    expression("  Benthic \nInvertivore"),
    vp = viewport(layout.pos.row = 3, layout.pos.col = 4),
    vjust = 0.5, gp = gpar(fontsize = 9), hjust = 0.45
    )
grid.text(
    "Zooplanktivore",
    vp = viewport(layout.pos.row = 3, layout.pos.col = 5),
    vjust = -1.2, gp = gpar(fontsize = 9), hjust = 0.5
    )
grid.text(
    "Herbivore",
    vp = viewport(layout.pos.row = 3, layout.pos.col = 6),
    vjust = -1.2, gp = gpar(fontsize = 9), hjust = 0.55
    )
grid.text("Functional group", vp = viewport(layout.pos.row = 4, layout.pos.col = 3:6),
    gp = gpar(fontsize = 9), vjust = -1, hjust = 0.5)
grid.text("Figure 2", vp = viewport(layout.pos.row = 4, layout.pos.col = 1:3),
    gp = gpar(fontsize = 9), vjust = 0, hjust = 1.5)

dev.copy2eps(device = quartz, file = "panel_plots/phylogenetically_corrected_allometrice_slope_estimates.eps")