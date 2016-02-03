################################################################################
#############                     Bootstrapping                    #############
################################################################################
# bootstrapping sma()
ddply(data = p, .(SpeciesCode), fun(x) sample(x, size = 12, replace = TRUE, prob = NULL))

# getting smallest sample within each functional group to write mk_all_spp_samp_df()
x <- ddply(.data = pento, .(j_fg, SpeciesCode), summarise, len = length(SpeciesCode))
ddply(x, .(j_fg), summarise, min(len))

#-------------------------------------------------------------------------------
# The bootstrapping
# The actual bootstrapping code has been commented out so that if this file is
# called with source, it will not spend 10 hours rerunning the code, it will 
# just load the SMA bootstrapped data from bootSMA_10000.RData located in the
# ~/R_projects/Allometry directory.

#m <- 10000  ## bootstrap replicates to use for offical MS values
#bootX <- rlply(m, mk_all_spp_samp_df(p, b, z, h, c), .progress = "text")

#registerDoMC()
#start = proc.time()
#bootSMA <- llply(bootX, function(x) {
#	sma(ga ~ SL * j_fg, data = x, log = "xy", method = "SMA", robust = T, 
#		slope.test = 2, multcomp = T, multcompmethod = "adjusted")
#	}, .parallel = TRUE)
#end = proc.time()
#total = (end[3]-start[3])
#cat(paste("Time Elapsed: ", total, "\n"))

#save(bootX, bootSMA, file="Allometry/bootSMA_10000.RData")
load("bootSMA_10000.RData")
load("boot_spp_SMA_gh_10000.RData")
load("boot_spp_SMA_gw_10000.RData")

#-------------------------------------------------------------------------------
# Allometry Check
# for each bootSMA and each functiona group in that bootSMA determine whether the 
# sma regression was negatively allometric, isometric, or positively allometric
slp_check_df <- data.frame(j_fg = c("Pi", "BI", "ZP", "He", "C"), iso = 0, 
						   neg = 0, pos = 0)

# Create dataframe to determine whether sma regression was neg, iso, or pos
# Checking allometry
check_df <- ldply(bootSMA, function(x) slp_check(x, slope_value = 2))
allo_check_df <- ddply(check_df, .(j_fg), summarise, 
             iso = sum(iso), neg = sum(neg), pos = sum(pos)
             )
allo_check_df$j_fg <- factor(allo_check_df$j_fg, levels = c("Pi", "BI", "ZP", "He", "C"))
allo_check_df <- arrange(allo_check_df, j_fg)

# Check allometry for gape height
check_df <- ldply(boot_spp_SMA_gh_10000, function(x) slp_check(x, slope_value = 1))
allo_check_df <- ddply(check_df, .(j_fg), summarise, 
             iso = sum(iso), neg = sum(neg), pos = sum(pos)
             )
allo_check_df$j_fg <- factor(allo_check_df$j_fg, levels = c("Pi", "BI", "ZP", "He", "C"))
allo_check_df <- arrange(allo_check_df, j_fg)

# Check allometry for gape width
check_df <- ldply(boot_spp_SMA_gw_10000, function(x) slp_check(x, slope_value = 1))
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
boot_summ$j_fg <- factor(boot_summ$j_fg, levels = c("Pi", "BI", "ZP", "He", "C"))
boot_means <- ddply(boot_summ, .(j_fg), get_mean_slope_int_r2)

# Calculate confidence intervals using studentized-t
boot_means$lowerCI <- boot_means$slope - (1.96 * boot_means$slope_se)
boot_means$upperCI <- boot_means$slope + (1.96 * boot_means$slope_se)
boot_means
#  j_fg    slope        int         r2   slope_se  lowerCI  upperCI
#1   Pi 1.531990 -0.3067777 0.67674821 0.09725793 1.341365 1.722616
#2   BI 1.521871 -0.6038637 0.82719876 0.08753882 1.350295 1.693447
#3   ZP 2.042337 -2.0824801 0.85913588 0.08884415 1.868203 2.216472
#4   He 2.449285 -3.2007591 0.89168539 0.10539008 2.242721 2.655850
#5    C 4.003380 -6.3401771 0.07056937 0.00000000 4.003380 4.003380

# Use regular SMA confidence interval for Corallivore lower and upper CI
cspp <- sma(ga~SL, data = c, log = "xy", method = "SMA", robust = T, slope.test = 2)
c_lowerCI <- cspp$coef[[1]]$lower[2]
c_upperCI <- cspp$coef[[1]]$upper[2]

# Replacing Corallivore CIs in boot_means
boot_means$lowerCI[5] <- c_lowerCI
boot_means$upperCI[5] <- c_upperCI
boot_means

# Getting plus-minus
boot_means$PlusMinus <- boot_means$slope - boot_means$lowerCI
boot_means


# Checking that my histograms are symmetric 
ggplot(data = boot_summ, aes(x = slope)) +
  geom_histogram() +
  geom_density() +
  facet_wrap(~ j_fg, scales = "free")

# Use the overall functional group to find the midpoint of the cluster of points
allGA <- sma(ga~SL, data = pento, log = "xy", method = "SMA", robust = T, slope.test = 2)
#check_assump(allGA, "Pento Gape Area All")
allGA_summ <- mk_sma_summary(allGA, 1)
allGA_graph_df <- mk_sma_graph_df(allGA_summ, 1, "j_fg")
names(allGA_graph_df)

# SMA regression for each functional group
all_fg_GA <- sma(ga~SL * j_fg, data = pento, log = "xy", method = "SMA", robust = T, 
  slope.test = 2, multcomp = T, multcompmethod = "adjusted")
all_fg_GA_summ <- mk_spp_summary(all_fg_GA, 5, grouping=TRUE)
all_fg_GA_graph_df <- mk_smaSPP_graph_df(all_fg_GA_summ, 5, "j_fg")

# Get x and y to and from values to plot geom_segment for each functional group

# Get x to and from values
to_from_n_df <- ddply(pento, .(j_fg), summarise, from = min(SL), to = max(SL), 
  n = length(j_fg))
boot_graph_df <- merge(x = to_from_n_df, y = boot_means, by = "j_fg")
boot_graph_df <- arrange(boot_graph_df, j_fg)

# Get y to and from values
boot_y_FromTo <- ddply(boot_graph_df, .(j_fg), get_y_FromTo)
boot_graph_df <- merge(x = boot_y_FromTo, y = boot_graph_df, by = "j_fg")
boot_graph_df <- arrange(boot_graph_df, j_fg)
boot_graph_df

# Get bootSMA equations for the graph
bootSMA_eqns <- ddply(boot_graph_df, .(j_fg), function(x) write_bootSMA_eqn(x))
bootSMA_eqns

allGA <- sma(ga~SL, data = pento, log = "xy", method = "SMA", robust = T, slope.test = 2)
#check_assump(allGA, "Pento Gape Area All")
allGA_summ <- mk_sma_summary(allGA, 1)
pento_line <- mk_sma_graph_df(allGA_summ, 1, "j_fg")

#===============================================================================
# Plot generation for gape size-body size MS

fg_sma <- sma(ga~SL * j_fg, data = pento, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")

fg_sma <- sma(ga~SL + j_fg, data = pento, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")


# Species differences in slopes within functional groups
pspp <- sma(ga~SL * SpeciesCode, data = p, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
bspp <- sma(ga~SL * SpeciesCode, data = b, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
zspp <- sma(ga~SL * SpeciesCode, data = z, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
hspp <- sma(ga~SL * SpeciesCode, data = h, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")

cspp <- sma(ga~SL, data = c, log = "xy", method = "SMA", robust = T, slope.test = 2)

# Species differences in elevation within functional groups
pspp_elev <- sma(ga~SL + SpeciesCode, data = p, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
bspp_elev <- sma(ga~SL + SpeciesCode, data = b, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
zspp_elev <- sma(ga~SL + SpeciesCode, data = z, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
hspp_elev <- sma(ga~SL + SpeciesCode, data = h, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")

#===============================================================================
# Multipanel plot with functional group and species
#===============================================================================
# SMA regression values for entire community
#allGA <- sma(ga~SL, data = pento, log = "xy", method = "SMA", robust = T, slope.test = 2)
#check_assump(allGA, "Pento Gape Area All")
#allGA_summ <- mk_sma_summary(allGA, 1)
#allGA_graph_df <- mk_sma_graph_df(allGA_summ, 1, "j_fg")

# SMA regression for each functional group
#all_fg_GA <- sma(ga~SL * j_fg, data = pento, log = "xy", method = "SMA", robust = T, 
#	slope.test = 2, multcomp = T, multcompmethod = "adjusted")
#check_assump(allGA, "Pento Gape Area All")
#all_fg_GA_summ <- mk_spp_summary(all_fg_GA, 5, grouping=TRUE)
#all_fg_GA_graph_df <- mk_smaSPP_graph_df(all_fg_GA_summ, 5, "j_fg")

# The above SMA regressions for each functional group have been superceded by 
# the output from the bootstrapped values for each functional group produced
# using the code in <bootstrap_for_ms.r> Below is the loading of the summary
# graph dataframe produced through this method. 
# This object is called: boot_graph_df
load("~/R_projects/Allometry/boot_graph_df_10000.RData")


# SMA regression for each species
all_spp_GA <- sma(ga~SL * SpeciesCode, data = pento, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
#check_assump(allGA, "Species Gape Area All")
allGA_bySPP_summ <- mk_spp_summary(all_spp_GA, 22, grouping=TRUE)
allGA_bySPP_graph_df <- mk_smaSPP_graph_df(allGA_bySPP_summ, 22, "SpeciesCode")

#-------------------------------------------------------------------------------
# Setting up values to plot the lines for the entire community, functional group
# level, and individual species level
#-------------------------------------------------------------------------------
#pento_line <- allGA_graph_df
fg_lines <- all_fg_GA_graph_df
spp_lines <- allGA_bySPP_graph_df

#-------------------------------------------------------------------------------
# Setting up equation, r^2, and n values that will be written on the graphs
#-------------------------------------------------------------------------------
#fg_sma_eqns <- write_group_sma_eqn(all_fg_GA_summ, all_fg_GA_summ$group)
#names(fg_sma_eqns) <- c("j_fg", "eqn_r2", "eqn", "r2", "n")

spp_sma_eqns <- write_group_sma_eqn(allGA_bySPP_summ, allGA_bySPP_summ$group)
names(spp_sma_eqns) <- c("SpeciesCode", "eqn_r2", "eqn", "r2", "n")

################################################################################
############  					Panel Plots 						############
################################################################################
# These plots have been superceded by the bootstrapped plots from 
# <bootstrap_for_ms.r>

# Functional group panel plot
# Notes: in the mk_multipanel_plots2 functions the following text settings were
# used. It may be best to add these as options in the function, but for now I'll
# just leave things like this. 

pisc_plot <- 
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = p, 
	spp_line_df_row = boot_graph_df[1, ], eqn_df = bootSMA_eqns[1, ], 
	eqn_x = 700, eqn_y = 1.55, r2_x = 700, r2_y = 4, n_x = 700, n_y = 8.4,
	fg_line_intercept = pento_line$ref_intercept, 
	x_axis_labels = FALSE, x_axis_text = TRUE, y_axis_labels = FALSE, 
	y_axis_text = TRUE, plot_title = "Piscivores") +
	scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
	scale_x_log10(breaks = c(50, 100, 500))
	#annotation_logticks(base = 10, sides = "b")
benth_plot <- 
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = b, 
	spp_line_df_row = boot_graph_df[2, ], eqn_df = bootSMA_eqns[2, ], 
	eqn_x = 700, eqn_y = 1.55, r2_x = 700, r2_y = 4, n_x = 700, n_y = 8.4,
	fg_line_intercept = pento_line$ref_intercept, 
	x_axis_labels = FALSE, x_axis_text = TRUE, y_axis_labels = FALSE, 
	y_axis_text = FALSE, plot_title = "Benthic Invertivores") +
	scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
	scale_x_log10(breaks = c(50, 100, 500))
zoop_plot <- 
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = z, 
	spp_line_df_row = boot_graph_df[3, ], eqn_df = bootSMA_eqns[3, ], 
	eqn_x = 700, eqn_y = 1.55, r2_x = 700, r2_y = 4, n_x = 700, n_y = 8.4,
	fg_line_intercept = pento_line$ref_intercept, 
	x_axis_labels = FALSE, x_axis_text = TRUE, y_axis_labels = FALSE, 
	y_axis_text = FALSE, plot_title = "Zooplanktivores") +
	scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
	scale_x_log10(breaks = c(50, 100, 500))
herb_plot <- 
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = h, 
	spp_line_df_row = boot_graph_df[4, ], eqn_df = bootSMA_eqns[4, ], 
	eqn_x = 700, eqn_y = 1.55, r2_x = 700, r2_y = 4, n_x = 700, n_y = 8.4,
	fg_line_intercept = pento_line$ref_intercept, 
	x_axis_labels = FALSE, x_axis_text = TRUE, y_axis_labels = FALSE, 
	y_axis_text = FALSE, plot_title = "Herbivores") +
	scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
	scale_x_log10(breaks = c(50, 100, 500))
coral_plot <- 
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = c, 
	spp_line_df_row = boot_graph_df[5, ], eqn_df = bootSMA_eqns[5, ], 
	eqn_x = 700, eqn_y = 1.55, r2_x = 700, r2_y = 4, n_x = 700, n_y = 8.4,
	fg_line_intercept = pento_line$ref_intercept, 
	x_axis_labels = FALSE, x_axis_text = TRUE, y_axis_labels = FALSE, 
	y_axis_text = FALSE, plot_title = "Corallivore" ) +
	scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
	scale_x_log10(breaks = c(50, 100, 500))
dev.new(height = 2.5, width = 9)
master_layout <- 
grid.layout(nrow = 2, ncol = 6, 
			widths = unit(c(0.1, 1, 0.9, 0.9, 0.9, 0.9), "null"), 
			heights = unit(c(1, 0.1), "null")
			)

# With "Figure" label
dev.new(height = 2.5, width = 9)
master_layout <- 
grid.layout(nrow = 2, ncol = 6, 
			widths = unit(c(0.1, 1, 0.9, 0.9, 0.9, 0.9), "null"), 
			heights = unit(c(1, 0.1), "null")
			)

grid.newpage()
pushViewport(viewport(layout = master_layout))
print(pisc_plot, vp = set_vp(1, 2))
print(benth_plot, vp = set_vp(1, 3))
print(zoop_plot, vp = set_vp(1, 4))
print(herb_plot, vp = set_vp(1, 5))
print(coral_plot, vp = set_vp(1, 6))

# Figure label
grid.text("Figure 2", vp = viewport(layout.pos.row = 2, layout.pos.col = 1), 
	gp = gpar(fontsize = 9), hjust = -1, vjust = 1)

grid.text(
	expression( paste("gape area (", mm^2, ")", sep = "") ), 
	vp = viewport(layout.pos.row = 1, layout.pos.col = 1),
	rot = 90, gp = gpar(fontsize = 9), 
	vjust = 1
	)
grid.text(
	"standard length (mm)",
	vp = viewport(layout.pos.row = 2, layout.pos.col = 4),
	vjust = -0.3, gp = gpar(fontsize = 9)
	)

dev.copy2eps(device = quartz, file = "panel_plots/fg_plot_figure_label.eps")

#-------------------------------------------------------------------------------
# Multipanel comparison of 3 predatory functional groups
#-------------------------------------------------------------------------------
apfurc <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$AP.FURC, 
	spp_line_df_row = spp_lines[2, ], #ref_intercept_row = spp_lines$ref_intercept[2], 
	eqn_df = spp_sma_eqns[2, ], eqn_x = 700, eqn_y = 130, r2_x = 700, r2_y = 250, 
	n_x = 700, n_y = 400, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[1], 
	x_axis_text = FALSE, y_axis_text = TRUE, plot_title = "Aphareus furca") +
	geom_abline(data = boot_graph_df[1, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250, 500)) +
	geom_point(aes(x = 100, y = 95), alpha = 0)
luboha <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$LU.BOHA, 
	spp_line_df_row = spp_lines[3, ], #ref_intercept_row = spp_lines$ref_intercept[3],
	eqn_df = spp_sma_eqns[3, ], eqn_x = 700, eqn_y = 130, r2_x = 700, r2_y = 250, 
	n_x = 700, n_y = 400, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[1], 
	x_axis_text = FALSE, y_axis_text = FALSE, plot_title = "Lutjanus bohar") +
	geom_abline(data = boot_graph_df[1, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250, 500)) +
	geom_point(aes(x = 100, y = 95), alpha = 0)
valout <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$VA.LOUT, 
	spp_line_df_row = spp_lines[6, ], #ref_intercept_row = spp_lines$ref_intercept[7],
	eqn_df = spp_sma_eqns[6, ], eqn_x = 700, eqn_y = 130, r2_x = 700, r2_y = 250, 
	n_x = 700, n_y = 400, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[1], 
	x_axis_text = FALSE, y_axis_text = FALSE, plot_title = "Variola louti") +
	geom_abline(data = boot_graph_df[1, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250, 500)) +
	geom_point(aes(x = 100, y = 95), alpha = 0)
ceargu <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CE.ARGU, 
	spp_line_df_row = spp_lines[4, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[4, ], eqn_x = 700, eqn_y = 130, r2_x = 700, r2_y = 250, 
	n_x = 700, n_y = 400, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[1], 
	x_axis_text = TRUE, y_axis_text = TRUE, plot_title = "Cephalopholis argus") +
	geom_abline(data = boot_graph_df[1, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250, 500)) +
	geom_point(aes(x = 100, y = 95), alpha = 0)
ceurod <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CE.UROD, 
	spp_line_df_row = spp_lines[5, ], #ref_intercept_row = spp_lines$ref_intercept[6], 
	eqn_df = spp_sma_eqns[5, ], eqn_x = 700, eqn_y = 130, r2_x = 700, r2_y = 250, 
	n_x = 700, n_y = 400, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[1], 
	x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Cephalopholis urodeta") +
	geom_abline(data = boot_graph_df[1, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250, 500)) +
	geom_point(aes(x = 100, y = 95), alpha = 0)
camela <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CA.MELA, 
	spp_line_df_row = spp_lines[1, ], #ref_intercept_row = spp_lines$ref_intercept[1],
	eqn_df = spp_sma_eqns[1, ], eqn_x = 700, eqn_y = 130, r2_x = 700, r2_y = 250, 
	n_x = 700, n_y = 400, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[1], 
	x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Caranx melampygus") +
	geom_abline(data = boot_graph_df[1, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250, 500)) +
	geom_point(aes(x = 100, y = 95), alpha = 0)
# Benthic invertivores
paarca <-
mk_multipanel_plots2(fg_point_df = b, spp_point_df = b_spp_dfs$PA.ARCA, 
	spp_line_df_row = spp_lines[7, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[7, ], eqn_x = 350, eqn_y = 10, r2_x = 350, r2_y = 20, 
	n_x = 350, n_y = 35, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[2], x_axis_text = TRUE, 
	y_axis_text = TRUE, plot_title = "Paracirrhites arcatus") +
	geom_abline(data = boot_graph_df[2, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250)) +
	geom_point(aes(x = 100, y = 7), alpha = 0)
painsu <-
mk_multipanel_plots2(fg_point_df = b, spp_point_df = b_spp_dfs$PA.INSU, 
	spp_line_df_row = spp_lines[9, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[9, ], eqn_x = 350, eqn_y = 10, r2_x = 350, r2_y = 20, 
	n_x = 350, n_y = 35, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[2], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Parupeneus insularis") +
	geom_abline(data = boot_graph_df[2, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250)) +
	geom_point(aes(x = 100, y = 7), alpha = 0)
mogran <-
mk_multipanel_plots2(fg_point_df = b, spp_point_df = b_spp_dfs$MO.GRAN, 
	spp_line_df_row = spp_lines[8, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[8, ], eqn_x = 350, eqn_y = 10, r2_x = 350, r2_y = 20, 
	n_x = 350, n_y = 35, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[2], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Monotaxis grandoculis") +
	geom_abline(data = boot_graph_df[2, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250)) +
	geom_point(aes(x = 100, y = 7), alpha = 0)
# Zooplanktivores
psbart <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.BART, 
	spp_line_df_row = spp_lines[13, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[13, ], eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 4.2, 
	n_x = 290, n_y = 8.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[3], x_axis_text = FALSE, 
	y_axis_text = TRUE, plot_title = "Pseudanthias bartlettorum") +
	geom_abline(data = boot_graph_df[3, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250))
catere <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$CA.TERE, 
	spp_line_df_row = spp_lines[10, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[10, ], eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 4.2, 
	n_x = 290, n_y = 8.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[3], x_axis_text = FALSE, 
	y_axis_text = FALSE, plot_title = "Caesio teres") +
	geom_abline(data = boot_graph_df[3, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250))
psdisp <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.DISP, 
	spp_line_df_row = spp_lines[14, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[14, ], eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 4.2, 
	n_x = 290, n_y = 8.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[3], x_axis_text = FALSE, 
	y_axis_text = FALSE, plot_title = "Pseudanthias dispar") +
	geom_abline(data = boot_graph_df[3, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250))
psoliv <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.OLIV, 
	spp_line_df_row = spp_lines[15, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[15, ], eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 4.2, 
	n_x = 290, n_y = 8.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[3], x_axis_text = TRUE, 
	y_axis_text = TRUE, plot_title = "Pseudanthias olivaceus") +
	geom_abline(data = boot_graph_df[3, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250))
pttile <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$PT.TILE, 
	spp_line_df_row = spp_lines[11, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[11, ], eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 4.2, 
	n_x = 290, n_y = 8.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[3], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Pterocaesio tile") +
	geom_abline(data = boot_graph_df[3, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250))
chvand <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$CH.VAND, 
	spp_line_df_row = spp_lines[12, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[12, ], eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 4.2, 
	n_x = 290, n_y = 8.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[3], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Chromis vanderbilti") +
	geom_abline(data = boot_graph_df[3, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250))

# Plotting multipanel piscivores and benthic invertivores	
dev.new(height = 10, width = 7)
master_layout <- 
grid.layout(nrow = 8, ncol = 4, 
			widths = unit(c(0.2, 1, 0.9, 0.9), "null"),
			heights = unit(c(1, 1, 0.2, 1, 0.2, 1, 1, 0.2), "null"))
grid.newpage()

# With "Figure" labelled
master_layout <- 
grid.layout(nrow = 8, ncol = 4, 
			widths = unit(c(0.2, 1, 0.9, 0.9), "null"),
			heights = unit(c(1, 1, 0.2, 1, 0.2, 1, 1, 0.2), "null"))
dev.new(height = 10, width = 7)

pushViewport(viewport(layout = master_layout))
# piscs
print(apfurc, vp = set_vp(1, 2))
print(luboha, vp = set_vp(1, 3))
print(valout, vp = set_vp(1, 4))
print(ceargu, vp = set_vp(2, 2))
print(ceurod, vp = set_vp(2, 3))
print(camela, vp = set_vp(2, 4))
# benths
print(paarca, vp = set_vp(4, 2))
print(painsu, vp = set_vp(4, 3))
print(mogran, vp = set_vp(4, 4))
# zoops
print(psbart, vp = set_vp(6, 2))
print(catere, vp = set_vp(6, 3))
print(psdisp, vp = set_vp(6, 4))
print(psoliv, vp = set_vp(7, 2))
print(pttile, vp = set_vp(7, 3))
print(chvand, vp = set_vp(7, 4))

# Figure label
grid.text("Figure 4", vp = viewport(layout.pos.row = 8, layout.pos.col = 1),
	gp = gpar(fontsize = 9), hjust = -1, vjust = 1)

# piscs
grid.text(
	"a)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1), 
	gp = gpar(fontsize = 9), vjust = -8
	)
grid.text(
	expression( paste("gape area (", mm^2, ")", sep = "") ), 
	vp = viewport(layout.pos.row = 1:2, layout.pos.col = 1),
	rot = 90, gp = gpar(fontsize = 9), 
	vjust = 1
	)
grid.text(
	"standard length (mm)",
	vp = viewport(layout.pos.row = 3, layout.pos.col = 3),
	vjust = -1.55, gp = gpar(fontsize = 9)
	)
# benths
grid.text(
	"b)", vp = viewport(layout.pos.row = 4, layout.pos.col = 1), 
	gp = gpar(fontsize = 9), vjust = -10
	)
grid.text(
	expression( paste("gape area (", mm^2, ")", sep = "") ), 
	vp = viewport(layout.pos.row = 4, layout.pos.col = 1),
	rot = 90, gp = gpar(fontsize = 9), 
	vjust = 1
	)
grid.text(
	"standard length (mm)",
	vp = viewport(layout.pos.row = 5, layout.pos.col = 3),
	vjust = -1.55, gp = gpar(fontsize = 9)
	)
# zoops
grid.text(
	"c)", vp = viewport(layout.pos.row = 6, layout.pos.col = 1), 
	gp = gpar(fontsize = 9), vjust = -8
	)
grid.text(
	expression( paste("gape area (", mm^2, ")", sep = "") ), 
	vp = viewport(layout.pos.row = 6:7, layout.pos.col = 1),
	rot = 90, gp = gpar(fontsize = 9), 
	vjust = 1
	)
grid.text(
	"standard length (mm)",
	vp = viewport(layout.pos.row = 8, layout.pos.col = 3),
	vjust = -1.55, gp = gpar(fontsize = 9)
	)

dev.copy2eps(device = quartz, file = "panel_plots/benth_pisc_zoop_panel_figure_label.eps")

#-------------------------------------------------------------------------------
# Multipanel herbivores
#-------------------------------------------------------------------------------
acnigr <-
mk_multipanel_plots2(fg_point_df = h, spp_point_df = h_spp_dfs$AC.NIGR, 
	spp_line_df_row = spp_lines[16, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[16, ], eqn_x = 440, eqn_y = 7, r2_x = 440, r2_y =16,
	n_x = 440, n_y = 30, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[4], x_axis_text = FALSE, 
	y_axis_text = TRUE, plot_title = "Acanthurus nigricans") +
	geom_abline(data = boot_graph_df[4, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250)) + 
	geom_point(aes(x = 40, y = 5), alpha = 0)
ceflav <-
mk_multipanel_plots2(fg_point_df = h, spp_point_df = h_spp_dfs$CE.FLAV, 
	spp_line_df_row = spp_lines[18, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[18, ], eqn_x = 440, eqn_y = 7, r2_x = 440, r2_y =16,
	n_x = 440, n_y = 30, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[4], x_axis_text = FALSE, 
	y_axis_text = FALSE, plot_title = "Centropyge flavissima") +
	geom_abline(data = boot_graph_df[4, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250)) + 
	geom_point(aes(x = 40, y = 5), alpha = 0)
chsord <-
mk_multipanel_plots2(fg_point_df  = h, spp_point_df  = h_spp_dfs$CH.SORD, 
	spp_line_df_row = spp_lines[19, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[19, ], eqn_x = 440, eqn_y = 7, r2_x = 440, r2_y =16,
	n_x = 440, n_y = 30, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[4], x_axis_text = FALSE, 
	y_axis_text = FALSE, plot_title = "Chlororus sordidus") +
	geom_abline(data = boot_graph_df[4, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250)) + 
	geom_point(aes(x = 40, y = 5), alpha = 0)
scrubr <-
mk_multipanel_plots2(fg_point_df  = h, spp_point_df  = h_spp_dfs$SC.RUBR, 
	spp_line_df_row = spp_lines[21, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[21, ], eqn_x = 440, eqn_y = 7, r2_x = 440, r2_y =16,
	n_x = 440, n_y = 30, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[4], x_axis_text = TRUE, 
	y_axis_text = TRUE, plot_title = "Scarus rubroviolaceus") +
	geom_abline(data = boot_graph_df[4, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250)) + 
	geom_point(aes(x = 40, y = 5), alpha = 0)
acoliv <-
mk_multipanel_plots2(fg_point_df = h, spp_point_df = h_spp_dfs$AC.OLIV, 
	spp_line_df_row = spp_lines[17, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[17, ], eqn_x = 440, eqn_y = 7, r2_x = 440, r2_y =16,
	n_x = 440, n_y = 30, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[4], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Acanthurus olivaceus") +
	geom_abline(data = boot_graph_df[4, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250)) + 
	geom_point(aes(x = 40, y = 5), alpha = 0)
scfren <-
mk_multipanel_plots2(fg_point_df  = h, spp_point_df  = h_spp_dfs$SC.FREN, 
	spp_line_df_row = spp_lines[20, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[20, ], eqn_x = 440, eqn_y = 7, r2_x = 440, r2_y =16,
	n_x = 440, n_y = 30, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = boot_graph_df$ref_intercept[4], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Scarus frenatus") +
	geom_abline(data = boot_graph_df[4, ], aes(slope = slope, intercept = int), 
				linetype = 2) +
	scale_x_log10(breaks = c(50, 100, 250)) + 
	geom_point(aes(x = 40, y = 5), alpha = 0)
dev.new(height = 4, width = 7)
master_layout <- 
grid.layout(nrow = 3, ncol = 4, 
			widths = unit(c(0.2, 1, 0.9, 0.9), "null"),
			heights = unit(c(1, 1, 0.2), "null"))
grid.newpage()

# With "Figure" label
dev.new(height = 4, width = 7)
master_layout <- 
grid.layout(nrow = 3, ncol = 4, 
			widths = unit(c(0.2, 1, 0.9, 0.9), "null"),
			heights = unit(c(1, 1, 0.2), "null"))
grid.newpage()


pushViewport(viewport(layout = master_layout))
print(acnigr, vp = set_vp(1, 2))
print(ceflav, vp = set_vp(1, 3))
print(chsord, vp = set_vp(1, 4))
print(scrubr, vp = set_vp(2, 2))
print(acoliv, vp = set_vp(2, 3))
print(scfren, vp = set_vp(2, 4))

	# Figure label
grid.text("Figure 5", vp = viewport(layout.pos.row = 3, layout.pos.col = 1), 
	gp = gpar(fontsize = 9), hjust = -1, vjust = 1)

grid.text(
	expression( paste("gape area (", mm^2, ")", sep = "") ), 
	vp = viewport(layout.pos.row = 1:2, layout.pos.col = 1),
	rot = 90, gp = gpar(fontsize = 9), 
	vjust = 1
	)
grid.text(
	"standard length (mm)",
	vp = viewport(layout.pos.row = 3, layout.pos.col = 3),
	vjust = -1.2, gp = gpar(fontsize = 9)
	)

dev.copy2eps(device = quartz, file = "panel_plots/herb_panel_figure_label.eps")

#===============================================================================
# Relative gape size
#===============================================================================

# Pento factored by functional group then slope
SpeciesCode <- c("AP.FURC", "LU.BOHA", "VA.LOUT", "CA.MELA", "CE.ARGU", "CE.UROD",
				 "PA.ARCA", "PA.INSU", "MO.GRAN",
				 "PS.BART", "CA.TERE", "PS.DISP", "PS.OLIV", "PT.TILE", "CH.VAND",
				 "AC.NIGR", "AC.OLIV", "CE.FLAV", "CH.SORD", "SC.FREN", "SC.RUBR", 
				 "CH.ORNA"
				 )

sp_name_by_slope <- 
	c("Aphareus furca", "Lutjanus bohar", "Variola louti", "Caranx melampygus", 
	  "Cephalopholis argus", "Cephalopholis urodeta", 
	  "Paracirrhites arcatus", "Parupeneus insularis", "Monotaxis grandoculis", 
	  "Pseudanthias bartlettorum", "Caesio teres", "Pseudanthias dispar", 
	  "Pseudanthias olivaceus", "Pterocaesio tile", "Chromis vanderbilti", 
	  "Acanthurus nigricans", "Acanthurus olivaceus", "Centropyge flavissima", 
	  "Chlororus sordidus", "Scarus frenatus", "Scarus rubroviolaceus", 
	  "Chaetodon ornatissimus"
	  )

spp_key <- data.frame(SpeciesCode, sp_name_by_slope)
pento_by_slope <- merge(x = pento, y = spp_key, all.x = TRUE, all.y = FALSE)
#pento_by_slope$SpeciesCode <- factor(pento_by_slope$SpeciesCode, levels = SpeciesCode)
pento_by_slope$sp_name_by_slope <- factor(pento_by_slope$sp_name_by_slope, levels = sp_name_by_slope)

#-------------------------------------------------------------------------------
# Removing ceargu_out (CE.ARGU outlier)
#-------------------------------------------------------------------------------
# Finding outlier in boxplot Cephalopholis argus
ceargu_df <- fish[(which(fish$SpeciesCode == "CE.ARGU")), ]

ceargu_out <- which(ceargu_df$ga_ratio == max(ceargu_df$ga_ratio))
ceargu_df[ceargu_out, ]
#    SpecimenID     Family       Order         Genus SpeciesCode j_fg Site
#633  KIF12_171 Serranidae Perciformes Cephalopholis     CE.ARGU   Pi   40
#    Region    TL    SL FL  wt length.cm. a..cm. b..cm. calc_wt   gh   gw
#633  HP.MF 138.6 134.2 NA 380         NA     NA     NA      NA 63.9 66.5
#    dissected_by  stomach_contents prey_size coll_notes dis_notes       ga
#633           MW shrimp; see photo    30.1mm                      3337.432
#    gh_ratio  gw_ratio  ga_ratio observer_id
#633 0.476155 0.4955291 0.1853136          13

bp_outlier <- which(pento_by_slope$SpecimenID == "KIF12_171")
pento_by_slope <- pento_by_slope[-bp_outlier, ]


# With facets:
fgs <- list('Pi' = "PI",
            'BI' = "BI", 
            'ZP' = "ZP", 
            'He' = "HE", 
            'C'  = "CO")

fg_labeller <- function(variable, value){
  return(fgs[value])
}

rel_ga <- ggplot(pento_by_slope, aes(.id, value, dodge = j_fg)) +
  geom_boxplot(aes(x=sp_name_by_slope, y=ga/(SL^2))) +
  xlab("Species") +
  ylab("Relative gape area") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle=45, colour="black", hjust=1, size = 9)) +
  theme(axis.text.y = element_text(colour="black")) +
  theme(axis.title.x = element_text(vjust = -0.2, size = 9)) +
  theme(axis.title.y = element_text(vjust = 0.15, size = 9)) +
  theme(legend.position = c(0.64, 0.9)) +
  theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
  theme(legend.text = element_text(size = 8)) +
  #theme(legend.key.size = unit(0.5, "cm")) +
  #theme(legend.key = element_blank()) + 
  #theme(legend.direction = "horizontal") +
  #guides(fill = guide_legend(override.aes = list(size = 0))) +
  #scale_fill_discrete(name = element_blank(),
  #					  labels=c("Piscivore", "Benthic \ninvertivore", 
  #					  		   "Zooplanktivore", "Herbivore", "Corallivore"))
  #scale_fill_discrete(name = "Functional \n Group")
  theme(panel.margin = unit(0.5, "cm")) +
  theme(panel.border = element_blank()) +
  theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(axis.line = element_line()) +
  facet_grid(. ~ j_fg, space = "free", scales = "free", labeller=fg_labeller)
#dev.off()
dev.new(height = 5, width = 8)
rel_ga

# With "Figure" label
master_layout <- 
grid.layout(nrow = 2, ncol = 1, 
			widths = unit(c(1), "null"),
			heights = unit(c(1, 0.05), "null"))
grid.newpage()
pushViewport(viewport(layout = master_layout))
print(rel_ga, vp = set_vp(1, 1))

# Figure label
grid.text("Figure 3", vp = viewport(layout.pos.row = 2, layout.pos.col = 1), 
	gp = gpar(fontsize = 9), hjust = 8, vjust = -1)
dev.copy2eps(device = quartz, file = "panel_plots/rel_ga_figure_label.eps")


# Absolute gape area
abs_ga <- ggplot(pento_by_slope, aes(.id, value, dodge = j_fg)) +
  geom_boxplot(aes(x=sp_name_by_slope, y=ga)) +
  xlab("Species") +
  ylab("Relative gape area") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle=45, colour="black", hjust=1, size = 9)) +
  theme(axis.text.y = element_text(colour="black")) +
  theme(axis.title.x = element_text(vjust = -0.2, size = 9)) +
  theme(axis.title.y = element_text(vjust = 0.15, size = 9)) +
  theme(legend.position = c(0.64, 0.9)) +
  theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
  theme(legend.text = element_text(size = 8)) +
  #theme(legend.key.size = unit(0.5, "cm")) +
  #theme(legend.key = element_blank()) + 
  #theme(legend.direction = "horizontal") +
  #guides(fill = guide_legend(override.aes = list(size = 0))) +
  #scale_fill_discrete(name = element_blank(),
  #					  labels=c("Piscivore", "Benthic \ninvertivore", 
  #					  		   "Zooplanktivore", "Herbivore", "Corallivore"))
  #scale_fill_discrete(name = "Functional \n Group")
  theme(panel.margin = unit(0.5, "cm")) +
  theme(panel.border = element_blank()) +
  theme(strip.background = element_rect(fill = "white", colour = "white")) +
  facet_grid(. ~ j_fg, space = "free", scales = "free", labeller=fg_labeller)

dev.off()
dev.new(height = 5, width = 8)
abs_ga


# WITHOUT FACETS
rel_ga <- ggplot(pento_by_slope, aes(.id, value, dodge = j_fg)) +
  geom_boxplot(aes(x=sp_name_by_slope, y=ga/(SL^2))) +
  xlab("Species") +
  ylab("Relative gape area") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle=45, colour="black", hjust=1, size = 9)) +
  theme(axis.text.y = element_text(colour="black")) +
  theme(axis.title.x = element_text(vjust = -0.2, size = 9)) +
  theme(axis.title.y = element_text(vjust = 0.15, size = 9)) +
  theme(legend.position = c(0.64, 0.9)) +
  theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
  theme(legend.text = element_text(size = 8))
  #theme(legend.key.size = unit(0.5, "cm")) +
  #theme(legend.key = element_blank()) + 
  #theme(legend.direction = "horizontal") +
  #guides(fill = guide_legend(override.aes = list(size = 0))) +
  #scale_fill_discrete(name = element_blank(),
  #					  labels=c("Piscivore", "Benthic \ninvertivore", 
  #					  		   "Zooplanktivore", "Herbivore", "Corallivore"))
  #scale_fill_discrete(name = "Functional \n Group")

dev.off()
dev.new(height = 5, width = 8)
rel_ga

dev.copy2eps(device = quartz, file = "panel_plots/rel_ga.eps")


rel_gh <- ggplot(pento_by_slope, aes(.id, value, dodge = j_fg)) +
  geom_boxplot(aes(x=sp_name_by_slope, y=gh^3/wt)) +
  xlab("Species") +
  ylab("Relative gape height") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle=45, colour="black", hjust=1, size = 9)) +
  theme(axis.text.y = element_text(colour="black")) +
  theme(axis.title.x = element_text(vjust = -0.2, size = 9)) +
  theme(axis.title.y = element_text(vjust = 0.15, size = 9)) +
  theme(legend.position = c(0.64, 0.9)) +
  theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
  theme(legend.text = element_text(size = 8))

rel_gw <- ggplot(pento_by_slope, aes(.id, value, dodge = j_fg)) +
  geom_boxplot(aes(x=sp_name_by_slope, y=gw^3/wt)) +
  xlab("Species") +
  ylab("Relative gape width") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle=45, colour="black", hjust=1, size = 9)) +
  theme(axis.text.y = element_text(colour="black")) +
  theme(axis.title.x = element_text(vjust = -0.2, size = 9)) +
  theme(axis.title.y = element_text(vjust = 0.15, size = 9)) +
  theme(legend.position = c(0.64, 0.9)) +
  theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
  theme(legend.text = element_text(size = 8))

  theme_bw() +
  theme( plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank()
  ) +
    theme(axis.line = element_line(color = 'black'))


ggplot(data = b, aes(x = SL, y = ga/(SL^2), colour = SpeciesCode)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Standard length, mm") +
  ylab("Relative gape area")


#===============================================================================
# Predator prey size graphs
#===============================================================================

#===============================================================================
# STANDARD LENGTH PRED - PREY
# Graphs for piscivores prey size - predator STANDARD LENGTH
P_prey_SL <- prey_PiBI[which(prey_PiBI$j_fg == 'Pi'), ]
df.n <- ddply(.data=P_prey_SL, .(j_fg), summarize, n=paste("n ==", length(SL)))
pisc_prey <-
ggplot(data = P_prey_SL, aes(x = SL, y = pSize)) +
  geom_point(aes(shape = pType)) +
  scale_shape_manual(values=c(1, 19)) +
  geom_text(data = df.n, aes(x = 90, y = 230, label = n), parse = TRUE, 
            size = 3, hjust = 0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  theme(legend.position = "none") +
  theme(axis.title = element_blank()) +
  #geom_smooth(method = "lm") +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.50, 0.90), method = "rq", 
                colour = "black") +
  theme(axis.ticks.length = unit(-0.2, "cm")) +
  theme(axis.ticks.margin = unit(0.3, "cm"))

# Benthic invertivore stomach contents for predator STANDARD LENGTH
B_prey_SL <- prey_PiBI[which(prey_PiBI$j_fg == 'BI'), ]
# Graph of just benthic invertivore stomach contents:
df.n <- ddply(.data=B_prey_SL, .(j_fg), summarize, n=paste("n ==", length(j_fg)))
benth_prey <- 
ggplot(data = B_prey_SL, aes(x = SL, y = pSize)) +
  geom_point(aes(shape = pType)) +
  geom_point(aes(x = 510, y = 100), alpha = 0.0) +
  geom_point(aes(x = 49, y = 258), alpha = 0.0) +
  scale_shape_manual(values=c(1, 19)) +
  geom_text(data = df.n, aes(x = 90, y = 230, label = n), parse = TRUE, 
            size = 3, hjust = 0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  theme(legend.position = "none") +
  theme(axis.title = element_blank()) +
  #geom_smooth(method = "lm") +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.50, 0.90), method = "rq", 
                colour = "black") +
  theme(axis.ticks.length = unit(-0.2, "cm")) +
  theme(axis.ticks.margin = unit(0.3, "cm"))

#camela_prey_SL <- prey_PiBI[which(prey_PiBI$SpeciesCode == 'CA.MELA'), ]
#ggplot(data = camela_prey_SL, aes(x = SL, y = pSize)) +
 # geom_point(aes(shape = pType))

dev.new(height = 3.2, width = 7.5)

master_layout <- 
grid.layout(nrow = 2, ncol = 4, 
      widths = unit(c(0.1, 1, 0.1, 1), "null"),
      heights = unit(c(1, 0.15), "null"))
grid.newpage()
pushViewport(viewport(layout = master_layout))
print(pisc_prey, vp = set_vp(1, 2))
print(benth_prey, vp = set_vp(1, 4))
grid.text(
  expression( paste("Standard length (", mm, ")", sep = "") ), 
  vp = viewport(layout.pos.row = 2, layout.pos.col = 2:4), 
  gp = gpar(fontsize = 10), vjust = -0.25
  )
grid.text(
  "Prey total length (mm)",
  vp = viewport(layout.pos.row = 1, layout.pos.col = 1), 
  gp = gpar(fontsize = 10), rot = 90, vjust = 2
  )
grid.text(
  "a)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1), 
  gp = gpar(fontsize = 9), vjust = -13
  )
grid.text(
  "b)", vp = viewport(layout.pos.row = 1, layout.pos.col = 3), 
  gp = gpar(fontsize = 9), vjust = -13
  )

# Figure label
grid.text("Figure 6", vp = viewport(layout.pos.row = 2, layout.pos.col = 1), 
	gp = gpar(fontsize = 9), hjust = -1, vjust = 1)

dev.copy2eps(device = quartz, file = "panel_plots/pred_prey_SL_figure_label.eps")

#===============================================================================
# Quantile Regressions
#===============================================================================
# Quantile regression for prey size ~ prey SL relationships
p_10_SL <- rq(data = P_prey_SL, pSize~SL, tau = 0.1)
p_50_SL <- rq(data = P_prey_SL, pSize~SL, tau = 0.5)
p_90_SL <- rq(data = P_prey_SL, pSize~SL, tau = 0.9)

summary(p_10_SL)
#Call: rq(formula = pSize ~ SL, tau = 0.1, data = P_prey_SL)
#tau: [1] 0.1
#Coefficients:
#            coefficients lower bd  upper bd 
#(Intercept) -11.87179    -17.07483   7.80187
#SL            0.10256     -0.00543   0.12278
# SL CI contains zero - not significant relationship
summary(p_50_SL)
#Call: rq(formula = pSize ~ SL, tau = 0.5, data = P_prey_SL)
#tau: [1] 0.5
#Coefficients:
#            coefficients lower bd  upper bd 
#(Intercept) -11.39259    -23.02733   7.37066
#SL            0.22963      0.07650   0.29805
# ga CI does not contain zero - significant relationship
#summary(p_90_SL)
#Call: rq(formula = pSize ~ SL, tau = 0.9, data = P_prey_SL)
#tau: [1] 0.9
#Coefficients:
#            coefficients lower bd  upper bd 
#(Intercept)  10.26230    -15.00095  26.65528
#SL            0.32787      0.24871   0.51036
# ga CI does not contain zero - significant relationship


# Quantile regression for Bethic Invertivores Gape Area
b_10_SL <- rq(data = B_prey_SL, pSize~SL, tau = 0.1)
b_50_SL <- rq(data = B_prey_SL, pSize~SL, tau = 0.5)
b_90_SL <- rq(data = B_prey_SL, pSize~SL, tau = 0.9)

summary(b_10_SL)
#Call: rq(formula = pSize ~ SL, tau = 0.1, data = B_prey_SL)
#tau: [1] 0.1
#Coefficients:
#            coefficients lower bd  upper bd 
#(Intercept)  -1.03306    -15.93614   5.41412
#SL            0.02479     -0.00727   0.09538
# ga CI contains zero - not significant relationship
summary(b_50_SL)
#Call: rq(formula = pSize ~ ga, tau = 0.5, data = B_prey_SL)
#tau: [1] 0.5
#Coefficients:
#            coefficients lower bd upper bd
#(Intercept)  9.36709      3.24418 28.11570
#ga           0.00323     -0.01891  0.01193
# ga CI contains zero - not significant relationship
summary(b_90_SL)
#Call: rq(formula = pSize ~ SL, tau = 0.9, data = B_prey_SL)
#tau: [1] 0.9
#Coefficients:
#            coefficients lower bd upper bd
#(Intercept) 41.77647     22.21822 60.66797
#SL          -0.01176     -0.12021  0.06571
# ga CI contains zero - not significant relationship

# Quantile regression for prey size ~ prey GA relationships
p_10_ga <- rq(data = P_prey_GA, pSize~ga, tau = 0.1)
p_50_ga <- rq(data = P_prey_GA, pSize~ga, tau = 0.5)
p_90_ga <- rq(data = P_prey_GA, pSize~ga, tau = 0.9)

summary(p_10_ga)
#Call: rq(formula = pSize ~ ga, tau = 0.1, data = P_prey_GA)
#tau: [1] 0.1
#Coefficients:
#            coefficients lower bd  upper bd 
#(Intercept)  -1.19193    -15.23673   0.30247
#ga            0.00815      0.00505   0.00950
# ga CI does not contain zero - significant relationship (though marginally)
summary(p_50_ga)
#Call: rq(formula = pSize ~ ga, tau = 0.5, data = P_prey_GA)
#tau: [1] 0.5
#Coefficients:
#            coefficients lower bd upper bd
#(Intercept)  1.72908     -0.30403  9.37498
#ga           0.01724      0.01059  0.01956
# ga CI does not contain zero - significant relationship
summary(p_90_ga)
#Call: rq(formula = pSize ~ ga, tau = 0.9, data = P_prey_GA)
#tau: [1] 0.9
#Coefficients:
#            coefficients lower bd upper bd
#(Intercept) 34.26744      8.39304 86.76394
#ga           0.03420      0.01374  0.04846
# ga CI does not contain zero - significant relationship


# Quantile regression for Bethic Invertivores Gape Area
b_10_ga <- rq(data = B_prey_GA, pSize~ga, tau = 0.1)
b_50_ga <- rq(data = B_prey_GA, pSize~ga, tau = 0.5)
b_90_ga <- rq(data = B_prey_GA, pSize~ga, tau = 0.9)

summary(b_10_ga)
#Call: rq(formula = pSize ~ ga, tau = 0.1, data = B_prey_GA)
#tau: [1] 0.1
#Coefficients:
#            coefficients lower bd upper bd
#(Intercept)  0.55458     -4.87786  3.71942
#ga           0.00591     -0.00213  0.01307
# ga CI contains zero - not significant relationship
summary(b_50_ga)
#Call: rq(formula = pSize ~ ga, tau = 0.5, data = B_prey_GA)
#tau: [1] 0.5
#Coefficients:
#            coefficients lower bd upper bd
#(Intercept)  9.36709      3.24418 28.11570
#ga           0.00323     -0.01891  0.01193
# ga CI contains zero - not significant relationship
summary(b_90_ga)
#Call: rq(formula = pSize ~ ga, tau = 0.9, data = B_prey_GA)
#tau: [1] 0.9
#Coefficients:
#            coefficients lower bd upper bd
#(Intercept) 50.53799     36.24255 93.15310
#ga          -0.01379     -0.03299  0.00143
# ga CI contains zero - not significant relationship

p_10_SL <- rq(data = prey3[(which(prey3$fg=='Pi')), ], psize~sl, tau = 0.1)
p_50_SL <- rq(data = prey3[(which(prey3$fg=='Pi')), ], psize~sl, tau = 0.5)
p_90_SL <- rq(data = prey3[(which(prey3$fg=='Pi')), ], psize~sl, tau = 0.9)


p_90_SL <- rq(data = P_prey_GA, pSize~SL, tau = 0.9)

P_prey_GA$

#===============================================================================
# GAPE AREA PRED - PREY
# Piscivores stomach contents for predator GAPE AREA
P_prey_GA <- prey_ga[which(prey_ga$j_fg == "Pi"), ]
# Graph of just piscivores stomach contents:
df.n <- ddply(.data=P_prey_GA, .(j_fg), summarize, n=paste("n ==", length(j_fg)))
pisc_prey <-
ggplot(data = P_prey_GA, aes(x = ga, y = pSize)) +
  geom_point(aes(shape = pType)) +
  scale_shape_manual(values=c(1, 19)) +
  geom_text(data = df.n, aes(x = 800, y = 230, label = n), parse = TRUE, 
            size = 3, hjust = 0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  theme(legend.position = "none") +
  theme(axis.title = element_blank()) +
  #geom_smooth(method = "lm") +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.50, 0.90), method = "rq", 
                colour = "black") +
  theme(axis.ticks.length = unit(-0.2, "cm")) +
  theme(axis.ticks.margin = unit(0.3, "cm"))

# Benthic invertivore stomach contents for predator GAPE AREA
B_prey_GA <- prey_ga[which(prey_ga$j_fg == "BI"), ]
# Graph of just benthic invertivore stomach contents:
df.n <- ddply(.data=B_prey_GA, .(j_fg), summarize, n=paste("n ==", length(j_fg)))
benth_prey <- 
ggplot(data = B_prey_GA, aes(x = ga, y = pSize)) +
  geom_point(aes(shape = pType)) +
  geom_point(aes(x = 293.7, y = 258), alpha = 0.0) +
  scale_shape_manual(values=c(1, 19)) +
  geom_text(data = df.n, aes(x = 400, y = 230, label = n), parse = TRUE, 
            size = 3, hjust = 0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  theme(legend.position = "none") +
  theme(axis.title = element_blank()) +
  #geom_smooth(method = "lm") +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.50, 0.90), method = "rq", 
                colour = "black") +
  theme(axis.ticks.length = unit(-0.2, "cm")) +
  theme(axis.ticks.margin = unit(0.3, "cm")) +
  scale_x_continuous(limits = c(0, 3000)) 


dev.new(height = 3.2, width = 7.5)

master_layout <- 
grid.layout(nrow = 2, ncol = 4, 
      widths = unit(c(0.1, 1, 0.1, 1), "null"),
      heights = unit(c(1, 0.15), "null"))
grid.newpage()
pushViewport(viewport(layout = master_layout))
print(pisc_prey, vp = set_vp(1, 2))
print(benth_prey, vp = set_vp(1, 4))
grid.text(
  expression( paste("Gape area (", mm^2, ")", sep = "") ), 
  vp = viewport(layout.pos.row = 2, layout.pos.col = 2:4), 
  gp = gpar(fontsize = 10), vjust = -0.25
  )
grid.text(
  "Prey total length (mm)",
  vp = viewport(layout.pos.row = 1, layout.pos.col = 1), 
  gp = gpar(fontsize = 10), rot = 90, vjust = 2
  )
grid.text(
  "a)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1), 
  gp = gpar(fontsize = 9), vjust = -13
  )
grid.text(
  "b)", vp = viewport(layout.pos.row = 1, layout.pos.col = 3), 
  gp = gpar(fontsize = 9), vjust = -13
  )

dev.copy2eps(device = quartz, file = "panel_plots/pred_prey_GA.eps")

# Quantile regression for prey size ~ prey GA relationships
b_10 <- rq(formula = pSize ~ ga, tau = 0.1, data = B_prey_GA)
b_50 <- rq(formula = pSize ~ ga, tau = 0.5, data = B_prey_GA)
b_90 <- rq(formula = pSize ~ ga, tau = 0.9, data = B_prey_GA)

anova(b_10, b_50, b_90)

p_10 <- rq(formula = pSize ~ ga, tau = 0.1, data = P_prey_GA)
p_50 <- rq(formula = pSize ~ ga, tau = 0.5, data = P_prey_GA)
p_90 <- rq(formula = pSize ~ ga, tau = 0.9, data = P_prey_GA)

anova(p_10, p_50, p_90)


#===============================================================================
# Old versions:

# Graph of just piscivores stomach contents:
ggplot(na.omit(prey5), aes(x = pi*((gh/2)+ (gw/2)), y = psize)) +
  geom_point(aes(shape = ptype)) +
  scale_shape_manual(values=c(1, 19)) +
  #geom_text(data = df.n, aes(x = 100, y = 200, label = n), parse = TRUE) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  #geom_smooth(method = "lm") +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.50, 0.90), method = "rq", 
                colour = "black") +
  theme(axis.ticks.length = unit(-0.2, "cm")) +
  theme(axis.ticks.margin = unit(0.3, "cm")) +
  xlab(expression(paste("gape area (", mm^2, ")", sep= ""))) +
  ylab("prey total length (mm)")

# Graph of benthic invertivore and piscivores:
df.n <- ddply(.data=na.omit(prey5), .(fg), summarize, n=paste("n ==", length(fg)))
pisc_prey <-
ggplot(na.omit(prey5), aes(x = pi*((gh/2)+ (gw/2)), y = psize)) +
  geom_point(aes(shape = ptype)) +
  scale_shape_manual(values=c(1, 19)) +
  geom_text(data = df.n, aes(x = 90, y = 230, label = n), parse = TRUE, 
            size = 3, hjust = 0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  theme(legend.position = "none") +
  theme(axis.title = element_blank()) +
  #geom_smooth(method = "lm") +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.50, 0.90), method = "rq", 
                colour = "black") +
  theme(axis.ticks.length = unit(-0.2, "cm")) +
  theme(axis.ticks.margin = unit(0.3, "cm")) +
  xlab(expression(paste("gape area (", mm^2, ")", sep= ""))) +
  ylab("prey total length (mm)")


df.n <- ddply(.data=prey6, .(fg), summarize, n=paste("n ==", length(fg)))
benth_prey <-
ggplot(na.omit(prey6), aes(x = pi*((gh/2)+ (gw/2)), y = psize)) +
  geom_point(aes(shape = ptype)) +
  geom_point(aes(x = 293.7, y = 258), alpha = 0.0) +
  scale_shape_manual(values=c(1, 19)) +
  geom_text(data = df.n, aes(x = 50, y = 230, label = n), parse = TRUE, 
            size = 3, hjust = 0) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  theme(legend.position = "none") +
  theme(axis.title = element_blank()) +
  #geom_smooth(method = "lm") +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.50, 0.90), method = "rq", 
                colour = "black") +
  theme(axis.ticks.length = unit(-0.2, "cm")) +
  theme(axis.ticks.margin = unit(0.3, "cm")) +
  xlab(expression(paste("gape area (", mm^2, ")", sep= ""))) #+
  #ylab("prey total length (mm)")

