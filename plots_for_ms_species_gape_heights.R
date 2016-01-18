#===============================================================================
# Gape height multipanel plot with functional group and species
#===============================================================================
# SMA regression for each species
all_spp_GH <- sma(gh ~ SL * SpeciesCode, data = pento, log = "xy", method = "SMA", 
    robust = T, slope.test = 1, multcomp = T, multcompmethod = "adjusted")
#check_assump(allGA, "Species Gape Area All")
allGH_bySPP_summ <- mk_spp_summary(all_spp_GH, 22, grouping=TRUE)
allGH_bySPP_graph_df <- mk_smaSPP_graph_df(allGH_bySPP_summ, 22, "SpeciesCode")

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
############                    Panel Plots                         ############
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
