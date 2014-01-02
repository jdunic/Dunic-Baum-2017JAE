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
allGA <- sma(ga~SL, data = pento, log = "xy", method = "SMA", robust = T, slope.test = 2)
#check_assump(allGA, "Pento Gape Area All")
allGA_summ <- mk_sma_summary(allGA, 1)
allGA_graph_df <- mk_sma_graph_df(allGA_summ, 1, "j_fg")

# SMA regression for each functional group
all_fg_GA <- sma(ga~SL * j_fg, data = pento, log = "xy", method = "SMA", robust = T, 
	slope.test = 2, multcomp = T, multcompmethod = "adjusted")
#check_assump(allGA, "Pento Gape Area All")
all_fg_GA_summ <- mk_spp_summary(all_fg_GA, 5, grouping=TRUE)
all_fg_GA_graph_df <- mk_smaSPP_graph_df(all_fg_GA_summ, 5, "j_fg")

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
pento_line <- allGA_graph_df
fg_lines <- all_fg_GA_graph_df
spp_lines <- allGA_bySPP_graph_df

#-------------------------------------------------------------------------------
# Setting up equation, r^2, and n values that will be written on the graphs
#-------------------------------------------------------------------------------
fg_sma_eqns <- write_group_sma_eqn(all_fg_GA_summ, all_fg_GA_summ$group)
names(fg_sma_eqns) <- c("j_fg", "eqn_r2", "eqn", "r2", "n")

spp_sma_eqns <- write_group_sma_eqn(allGA_bySPP_summ, allGA_bySPP_summ$group)
names(spp_sma_eqns) <- c("SpeciesCode", "eqn_r2", "eqn", "r2", "n")

################################################################################
############  					Panel Plots 						############
################################################################################

# Functional group panel plot
# Notes: in the mk_multipanel_plots2 functions the following text settings were
# used. It may be best to add these as options in the function, but for now I'll
# just leave things like this. 

pisc_plot <- 
	mk_multipanel_plots2(fg_point_df = pento, spp_point_df = p, 
		spp_line_df_row = fg_lines[1, ], eqn_df = fg_sma_eqns[1, ], eqn_x = 540, 
		eqn_y = 1.25, r2_x = 540, r2_y = 2.6, n_x = 540, n_y = 4.9,
		fg_line_intercept = pento_line$ref_intercept, 
		x_axis_labels = FALSE, x_axis_text = TRUE, y_axis_labels = FALSE, 
		y_axis_text = TRUE, plot_title = "Piscivores") #+ 
		#geom_point(data = p, aes(colour = dissected_by))
benth_plot <- 
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = b, 
	spp_line_df_row = fg_lines[2, ], eqn_df = fg_sma_eqns[2, ], eqn_x = 540, 
	eqn_y = 1.25, r2_x = 540, r2_y = 2.6,  n_x = 540, n_y = 4.9,
	fg_line_intercept = pento_line$ref_intercept, 
	x_axis_labels = FALSE, x_axis_text = TRUE, y_axis_labels = FALSE, 
	y_axis_text = FALSE, plot_title = "Benthic Invertivores")
zoop_plot <- 
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = z, 
	spp_line_df_row = fg_lines[3, ], eqn_df = fg_sma_eqns[3, ], eqn_x = 540, 
	eqn_y = 1.25, r2_x = 540, r2_y = 2.6,  n_x = 540, n_y = 4.9,
	fg_line_intercept = pento_line$ref_intercept, 
	x_axis_labels = FALSE, x_axis_text = TRUE, y_axis_labels = FALSE, 
	y_axis_text = FALSE, plot_title = "Zooplanktivores")
herb_plot <- 
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = h, 
	spp_line_df_row = fg_lines[4, ], eqn_df = fg_sma_eqns[4, ], eqn_x = 540, 
	eqn_y = 1.25, r2_x = 540, r2_y = 2.6,  n_x = 540, n_y = 4.9,
	fg_line_intercept = pento_line$ref_intercept, 
	x_axis_labels = FALSE, x_axis_text = TRUE, y_axis_labels = FALSE, 
	y_axis_text = FALSE, plot_title = "Herbivores")
coral_plot <- 
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = c, 
	spp_line_df_row = fg_lines[5, ], eqn_df = fg_sma_eqns[5, ], eqn_x = 540, 
	eqn_y = 1.25, r2_x = 540, r2_y = 2.6,  n_x = 540, n_y = 4.9,
	fg_line_intercept = pento_line$ref_intercept, 
	x_axis_labels = FALSE, x_axis_text = TRUE, y_axis_labels = FALSE, 
	y_axis_text = FALSE, plot_title = "Corallivore" )
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
grid.text(
	expression( paste("log(gape area, ", mm^2, ")", sep = "") ), 
	vp = viewport(layout.pos.row = 1, layout.pos.col = 1),
	rot = 90, gp = gpar(fontsize = 9), 
	vjust = 1
	)
grid.text(
	"log(standard length, mm)",
	vp = viewport(layout.pos.row = 2, layout.pos.col = 4),
	vjust = -0.3, gp = gpar(fontsize = 9)
	)

dev.copy2eps(device = quartz, file = "panel_plots/fg_plot.eps")

#-------------------------------------------------------------------------------
# Multipanel comparison of 3 predatory functional groups
#-------------------------------------------------------------------------------
apfurc <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$AP.FURC, 
	spp_line_df_row = spp_lines[2, ], #ref_intercept_row = spp_lines$ref_intercept[2], 
	eqn_df = spp_sma_eqns[2, ], eqn_x = 700, eqn_y = 50, r2_x = 700, r2_y = 97, 
	n_x = 700, n_y = 170, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[1], 
	x_axis_text = FALSE, y_axis_text = TRUE, plot_title = "Aphareus furca") #+
	#geom_point(data = p_spp_dfs$AP.FURC, aes(colour = dissected_by), size = 3)
luboha <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$LU.BOHA, 
	spp_line_df_row = spp_lines[3, ], #ref_intercept_row = spp_lines$ref_intercept[3],
	eqn_df = spp_sma_eqns[3, ], eqn_x = 700, eqn_y = 50, r2_x = 700, r2_y = 97, 
	n_x = 700, n_y = 170, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[1], 
	x_axis_text = FALSE, y_axis_text = FALSE, plot_title = "Lutjanus bohar") #+
	#geom_point(data = p_spp_dfs$LU.BOHA, aes(colour = dissected_by), size = 3)
valout <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$VA.LOUT, 
	spp_line_df_row = spp_lines[6, ], #ref_intercept_row = spp_lines$ref_intercept[7],
	eqn_df = spp_sma_eqns[6, ], eqn_x = 700, eqn_y = 50, r2_x = 700, r2_y = 97, 
	n_x = 700, n_y = 170, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[1], 
	x_axis_text = FALSE, y_axis_text = FALSE, plot_title = "Variola louti") #+
	#geom_point(data = p_spp_dfs$VA.LOUT, aes(colour = dissected_by), size = 3)
ceargu <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CE.ARGU, 
	spp_line_df_row = spp_lines[4, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[4, ], eqn_x = 700, eqn_y = 50, r2_x = 700, r2_y = 97, 
	n_x = 700, n_y = 170, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[1], 
	x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Cephalopholis argus") #+
	#geom_point(data = p_spp_dfs$CE.ARGU, aes(colour = dissected_by), size = 3)
ceurod <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CE.UROD, 
	spp_line_df_row = spp_lines[5, ], #ref_intercept_row = spp_lines$ref_intercept[6], 
	eqn_df = spp_sma_eqns[5, ], eqn_x = 700, eqn_y = 50, r2_x = 700, r2_y = 97, 
	n_x = 700, n_y = 170, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[1], 
	x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Cephalopholis urodeta") #+
	#geom_point(data = p_spp_dfs$CE.UROD, aes(colour = dissected_by), size = 3)
camela <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CA.MELA, 
	spp_line_df_row = spp_lines[1, ], #ref_intercept_row = spp_lines$ref_intercept[1],
	eqn_df = spp_sma_eqns[1, ], eqn_x = 700, eqn_y = 50, r2_x = 700, r2_y = 97, 
	n_x = 700, n_y = 170, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[1], 
	x_axis_text = TRUE, y_axis_text = TRUE, plot_title = "Caranx melampygus") #+
	#geom_point(data = p_spp_dfs$CA.MELA, aes(colour = dissected_by), size = 3)
# Benthic invertivores
paarca <-
mk_multipanel_plots2(fg_point_df = b, spp_point_df = b_spp_dfs$PA.ARCA, 
	spp_line_df_row = spp_lines[7, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[7, ], eqn_x = 350, eqn_y = 10, r2_x = 350, r2_y = 20, 
	n_x = 350, n_y = 35, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[2], x_axis_text = TRUE, 
	y_axis_text = TRUE, plot_title = "Paracirrhites arcatus") #+ 
	#geom_point(data = b_spp_dfs$PA.ARCA, aes(colour = dissected_by))
painsu <-
mk_multipanel_plots2(fg_point_df = b, spp_point_df = b_spp_dfs$PA.INSU, 
	spp_line_df_row = spp_lines[9, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[9, ], eqn_x = 350, eqn_y = 10, r2_x = 350, r2_y = 20, 
	n_x = 350, n_y = 35, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[2], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Parupeneus insularis") #+ 
	#geom_point(data = b_spp_dfs$PA.INSU, aes(colour = dissected_by))
mogran <-
mk_multipanel_plots2(fg_point_df = b, spp_point_df = b_spp_dfs$MO.GRAN, 
	spp_line_df_row = spp_lines[8, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[8, ], eqn_x = 350, eqn_y = 10, r2_x = 350, r2_y = 20, 
	n_x = 350, n_y = 35, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[2], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Monotaxis grandoculis") #+ 
	#geom_point(data = b_spp_dfs$MO.GRAN, aes(colour = dissected_by))
# Zooplanktivores
psbart <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.BART, 
	spp_line_df_row = spp_lines[13, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[13, ], eqn_x = 250, eqn_y = 2, r2_x = 250, r2_y = 4.5, 
	n_x = 250, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[4], x_axis_text = FALSE, 
	y_axis_text = TRUE, plot_title = "Pseudanthias bartlettorum")
catere <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$CA.TERE, 
	spp_line_df_row = spp_lines[10, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[10, ], eqn_x = 250, eqn_y = 2, r2_x = 250, r2_y = 4.5, 
	n_x = 250, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[4], x_axis_text = FALSE, 
	y_axis_text = FALSE, plot_title = "Caesio teres")
psdisp <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.DISP, 
	spp_line_df_row = spp_lines[14, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[14, ], eqn_x = 250, eqn_y = 2, r2_x = 250, r2_y = 4.5, 
	n_x = 250, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[4], x_axis_text = FALSE, 
	y_axis_text = FALSE, plot_title = "Pseudanthias dispar")
psoliv <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.OLIV, 
	spp_line_df_row = spp_lines[15, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[15, ], eqn_x = 250, eqn_y = 2, r2_x = 250, r2_y = 4.5, 
	n_x = 250, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[4], x_axis_text = TRUE, 
	y_axis_text = TRUE, plot_title = "Pseudanthias olivaceus")
pttile <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$PT.TILE, 
	spp_line_df_row = spp_lines[11, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[11, ], eqn_x = 250, eqn_y = 2, r2_x = 250, r2_y = 4.5, 
	n_x = 250, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[4], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Pterocaesio tile")
chvand <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$CH.VAND, 
	spp_line_df_row = spp_lines[12, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[12, ], eqn_x = 250, eqn_y = 2, r2_x = 250, r2_y = 4.5, 
	n_x = 250, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[4], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Chromis vanderbilti")


# Plotting multipanel piscivores and benthic invertivores	
dev.new(height = 10, width = 7)
master_layout <- 
grid.layout(nrow = 8, ncol = 4, 
			widths = unit(c(0.2, 1, 0.9, 0.9), "null"),
			heights = unit(c(1, 1, 0.2, 1, 0.2, 1, 1, 0.2), "null"))
grid.newpage()
pushViewport(viewport(layout = master_layout))
# piscs
print(apfurc, vp = set_vp(1, 2))
print(luboha, vp = set_vp(1, 3))
print(valout, vp = set_vp(1, 4))
print(camela, vp = set_vp(2, 2))
print(ceargu, vp = set_vp(2, 3))
print(ceurod, vp = set_vp(2, 4))
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



# piscs
grid.text(
	"a)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1), 
	gp = gpar(fontsize = 9), vjust = -8
	)
grid.text(
	expression( paste("log(gape area, ", mm^2, ")", sep = "") ), 
	vp = viewport(layout.pos.row = 1:2, layout.pos.col = 1),
	rot = 90, gp = gpar(fontsize = 9), 
	vjust = 1
	)
grid.text(
	"log(standard length, mm)",
	vp = viewport(layout.pos.row = 3, layout.pos.col = 3),
	vjust = -1.55, gp = gpar(fontsize = 9)
	)
# benths
grid.text(
	"b)", vp = viewport(layout.pos.row = 4, layout.pos.col = 1), 
	gp = gpar(fontsize = 9), vjust = -10
	)
grid.text(
	expression( paste("log(gape area, ", mm^2, ")", sep = "") ), 
	vp = viewport(layout.pos.row = 4, layout.pos.col = 1),
	rot = 90, gp = gpar(fontsize = 9), 
	vjust = 1
	)
grid.text(
	"log(standard length, mm)",
	vp = viewport(layout.pos.row = 5, layout.pos.col = 3),
	vjust = -1.55, gp = gpar(fontsize = 9)
	)
# zoops
grid.text(
	"c)", vp = viewport(layout.pos.row = 6, layout.pos.col = 1), 
	gp = gpar(fontsize = 9), vjust = -8
	)
grid.text(
	expression( paste("log(gape area, ", mm^2, ")", sep = "") ), 
	vp = viewport(layout.pos.row = 6:7, layout.pos.col = 1),
	rot = 90, gp = gpar(fontsize = 9), 
	vjust = 1
	)
grid.text(
	"log(standard length, mm)",
	vp = viewport(layout.pos.row = 8, layout.pos.col = 3),
	vjust = -1.55, gp = gpar(fontsize = 9)
	)

dev.copy2eps(device = quartz, file = "panel_plots/benth_pisc_zoop_panel.eps")

#-------------------------------------------------------------------------------
# Multipanel herbivores
#-------------------------------------------------------------------------------
acnigr <-
mk_multipanel_plots2(fg_point_df = h, spp_point_df = h_spp_dfs$AC.NIGR, 
	spp_line_df_row = spp_lines[16, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[16, ], eqn_x = 400, eqn_y = 7, r2_x = 400, r2_y =16,
	n_x = 400, n_y = 30, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[3], x_axis_text = FALSE, 
	y_axis_text = TRUE, plot_title = "Acanthurus nigricans")
acoliv <-
mk_multipanel_plots2(fg_point_df = h, spp_point_df = h_spp_dfs$AC.OLIV, 
	spp_line_df_row = spp_lines[17, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[17, ], eqn_x = 400, eqn_y = 7, r2_x = 400, r2_y =16,
	n_x = 400, n_y = 30, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[3], x_axis_text = FALSE, 
	y_axis_text = FALSE, plot_title = "Acanthurus olivaceus")
ceflav <-
mk_multipanel_plots2(fg_point_df = h, spp_point_df = h_spp_dfs$CE.FLAV, 
	spp_line_df_row = spp_lines[18, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[18, ], eqn_x = 400, eqn_y = 7, r2_x = 400, r2_y =16,
	n_x = 400, n_y = 30, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[3], x_axis_text = FALSE, 
	y_axis_text = FALSE, plot_title = "Centropyge flavissima")
chsord <-
mk_multipanel_plots2(fg_point_df  = h, spp_point_df  = h_spp_dfs$CH.SORD, 
	spp_line_df_row = spp_lines[19, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[19, ], eqn_x = 400, eqn_y = 7, r2_x = 400, r2_y =16,
	n_x = 400, n_y = 30, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[3], x_axis_text = TRUE, 
	y_axis_text = TRUE, plot_title = "Chlororus sordidus")
scfren <-
mk_multipanel_plots2(fg_point_df  = h, spp_point_df  = h_spp_dfs$SC.FREN, 
	spp_line_df_row = spp_lines[20, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[20, ], eqn_x = 400, eqn_y = 7, r2_x = 400, r2_y =16,
	n_x = 400, n_y = 30, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[3], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Scarus frenatus")
scrubr <-
mk_multipanel_plots2(fg_point_df  = h, spp_point_df  = h_spp_dfs$SC.RUBR, 
	spp_line_df_row = spp_lines[21, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[21, ], eqn_x = 400, eqn_y = 7, r2_x = 400, r2_y =16,
	n_x = 400, n_y = 30, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[3], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Scarus rubroviolaceus")
dev.new(height = 4, width = 7)
master_layout <- 
grid.layout(nrow = 3, ncol = 4, 
			widths = unit(c(0.2, 1, 0.9, 0.9), "null"),
			heights = unit(c(1, 1, 0.2), "null"))
grid.newpage()
pushViewport(viewport(layout = master_layout))
print(acnigr, vp = set_vp(1, 2))
print(acoliv, vp = set_vp(1, 3))
print(ceflav, vp = set_vp(1, 4))
print(chsord, vp = set_vp(2, 2))
print(scfren, vp = set_vp(2, 3))
print(scrubr, vp = set_vp(2, 4))
grid.text(
	expression( paste("log(gape area, ", mm^2, ")", sep = "") ), 
	vp = viewport(layout.pos.row = 1:2, layout.pos.col = 1),
	rot = 90, gp = gpar(fontsize = 9), 
	vjust = 1
	)
grid.text(
	"log(standard length, mm)",
	vp = viewport(layout.pos.row = 3, layout.pos.col = 3),
	vjust = -1.2, gp = gpar(fontsize = 9)
	)

dev.copy2eps(device = quartz, file = "panel_plots/herb_panel.eps")

#===============================================================================
# Relative gape size
#===============================================================================
# Pento factored by functional group then slope
spp_by_slope <- c("AP.FURC", "LU.BOHA", "VA.LOUT", "CA.MELA", "CE.ARGU", "CE.UROD",
				 "PA.ARCA", "PA.INSU", "MO.GRAN",
				 "PS.BART", "CA.TERE", "PS.DISP", "PS.OLIV", "PT.TILE", "CH.VAND",
				 "AC.NIGR", "AC.NIGR", "AC.OLIV", "CE.FLAV", "CH.SORD", "SC.FREN", "SC.RUBR",
				 "CH.ORNA"
				 )
pento_by_slope <- pento
pento_by_slope$SpeciesCode <- factor(pento_by_slope$SpeciesCode, levels = spp_by_slope)

rel_ga <- ggplot(pento_by_slope, aes(x=SpeciesCode, y=ga/(SL^2), fill=j_fg)) +
  geom_boxplot() +
  xlab("Species") +
  ylab("Relative gape area") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle=45, colour="black", vjust=0.5, size = 8)) +
  theme(axis.text.y = element_text(colour="black")) +
  theme(axis.title.x = element_text(vjust = -0.2, size = 10)) +
  theme(axis.title.y = element_text(vjust = 0.15, size = 10)) +
  theme(legend.position = c(0.95, 0.75)) +
  theme(legend.title = element_text(size = 8)) +
  theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
  scale_fill_discrete(name = "Functional \n Group")
  #scale_fill_discrete(name = "Functional \n Group")
rel_ga

dev.copy2eps(device = quartz, file = "panel_plots/rel_ga.eps")


rel_gh <- ggplot(pento_by_slope, aes(x=SpeciesCode, y=gh_ratio, fill=j_fg)) +
  geom_boxplot() +
  xlab("Species") +
  ylab("Relative gape height") +
  theme(axis.text.x = element_text(angle=45, colour="black", vjust=0.5)) +
  theme(axis.text.y = element_text(colour="black")) +
  theme(axis.title.y = element_text(vjust = 0.3)) +
  scale_fill_discrete(name = "Functional \n Group")

rel_gw <- ggplot(pento_by_slope, aes(x=SpeciesCode, y=gw_ratio, fill=j_fg)) +
  geom_boxplot() +
  xlab("Species") +
  ylab("Relative gape width") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=45, colour="black", vjust=0.5)) +
  theme(axis.text.y = element_text(colour="black")) +
  theme(axis.title.y = element_text(vjust = 0.3)) +
  scale_fill_discrete(name = "Functional \n Group")


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