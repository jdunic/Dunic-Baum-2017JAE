#===============================================================================
# Multipanel plot with functional group and species
#===============================================================================

# Get community level midpoint
allGW <- sma(gw ~ SL, data = pento, log = "xy", method = "SMA", robust = T, slope.test = 1)
#check_assump(allGW, "Pento Gape Area All")
allGW_summ <- mk_sma_summary(allGW, 1)
allGW_graph_df <- mk_sma_graph_df(allGW_summ, 1, "j_fg")
names(allGW_graph_df)

# Get functional group level metrics: midpoints of data
# Using a single sma() iteration because this does not alter results, nor the 
# visual display of the data.
pento_slopes_gw <- coef(fit.gls_gw)
names(pento_slopes_gw) <- c('Pi', 'Bi', 'Zp', 'He')

all_fg_GW <- sma(gw ~ SL * j_fg, data = pento, log = "xy", method = "SMA", 
                 robust = T, slope.test = 1, multcomp = F)
all_fg_GW_summ <- mk_spp_summary(all_fg_GW, 5, grouping=TRUE)
all_fg_GW_graph_df <- mk_smaSPP_graph_df(all_fg_GW_summ, 5, "j_fg")
all_fg_GW_graph_df$boot_slope <- c(pento_slopes_gh, NA)

all_fg_GW_graph_df <- 
  all_fg_GW_graph_df %>% 
    dplyr::mutate(boot_ref_int = log10(midpoint_y / (midpoint_x ^ boot_slope)))

# SMA regression for each species
all_spp_GW <- sma(gw ~ SL * SpeciesCode, data = pento, log = "xy", 
                  method = "SMA", robust = T, slope.test = 1,
                  multcomp = F, multcompmethod = "adjusted")
#check_assump(allGA, "Species Gape Area All")
allGW_bySPP_summ <- mk_spp_summary(all_spp_GW, 22, grouping=TRUE)
allGW_bySPP_graph_df <- mk_smaSPP_graph_df(allGW_bySPP_summ, 22, "SpeciesCode")

rbind_all(list(mutate(allGH_bySPP_summ, gape = 'GH'), dplyr::mutate(allGW_bySPP_summ, gape = 'GW'))) %>% 
  readr::write_csv(., 'gape_height_width_species_summary.csv')




#-------------------------------------------------------------------------------
# Setting up values to plot the lines at the species level
spp_lines <- allGW_bySPP_graph_df

# Setting up equation, r^2, and n values that will be written on the graphs
spp_sma_eqns <- write_group_sma_eqn(allGW_bySPP_summ, allGW_bySPP_summ$group)
names(spp_sma_eqns) <- c("SpeciesCode", "eqn_r2", "eqn", "r2", "n")


spp_lines[1, ]

ggplot() + 
  theme_bw() +
  geom_point(data = p, aes(x = SL, y = gh), colour = 'grey', shape = 1) +
  geom_point(data = p_spp_dfs$CA.MELA, aes(x = SL, y = gh), colour = 'black') + 
  geom_segment(aes(x = 178.3708, xend = 577.9654, y = 10^(1.0391357*log10(178.3708) + -0.9399020), yend = 10^(1.0391357*log10(577.9654) + -0.9399020))) + 
  scale_x_log10() + 
  scale_y_log10()

#-------------------------------------------------------------------------------
# Multipanel comparison of 3 predatory functional groups
#-------------------------------------------------------------------------------
apfurc <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$AP.FURC, 
    spp_line_df_row = spp_lines[2, ], eqn_df = spp_sma_eqns[2, ], 
    eqn_x = 700, eqn_y = 20, r2_x = 700, r2_y = 25.8,
    n_x = 700, n_y = 31, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[1], 
    x_axis_text = FALSE, y_axis_text = TRUE, 
    plot_title = "Aphareus furca", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[1, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250, 500)) +
    geom_point(aes(x = 100, y = 95), alpha = 0)
luboha <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$LU.BOHA, 
    spp_line_df_row = spp_lines[3, ], eqn_df = spp_sma_eqns[3, ], 
    eqn_x = 700, eqn_y = 20, r2_x = 700, r2_y = 25.8, 
    n_x = 700, n_y = 31, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept[1], 
    x_axis_text = FALSE, y_axis_text = FALSE, 
    plot_title = "Lutjanus bohar", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[1, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250, 500)) +
    geom_point(aes(x = 100, y = 95), alpha = 0)
valout <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$VA.LOUT, 
    spp_line_df_row = spp_lines[6, ], eqn_df = spp_sma_eqns[6, ], 
    eqn_x = 700, eqn_y = 20, r2_x = 700, r2_y = 25.8, 
    n_x = 700, n_y = 31, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept[1], 
    x_axis_text = FALSE, y_axis_text = FALSE, 
    plot_title = "Variola louti", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[1, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250, 500)) +
    geom_point(aes(x = 100, y = 95), alpha = 0)
ceargu <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CE.ARGU, 
    spp_line_df_row = spp_lines[4, ], eqn_df = spp_sma_eqns[4, ], 
    eqn_x = 700, eqn_y = 20, r2_x = 700, r2_y = 25.8, 
    n_x = 700, n_y = 31, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept[1], 
    x_axis_text = TRUE, y_axis_text = TRUE, 
    plot_title = "Cephalopholis argus", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[1, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250, 500)) +
    geom_point(aes(x = 100, y = 95), alpha = 0)
ceurod <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CE.UROD, 
    spp_line_df_row = spp_lines[5, ], eqn_df = spp_sma_eqns[5, ], 
    eqn_x = 700, eqn_y = 20, r2_x = 700, r2_y = 25.8, 
    n_x = 700, n_y = 31, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept[1], 
    x_axis_text = TRUE, y_axis_text = FALSE, 
    plot_title = "Cephalopholis urodeta", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[1, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250, 500)) +
    geom_point(aes(x = 100, y = 95), alpha = 0)
camela <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CA.MELA, 
    spp_line_df_row = spp_lines[1, ], eqn_df = spp_sma_eqns[1, ], 
    eqn_x = 700, eqn_y = 20, r2_x = 700, r2_y = 25.8, 
    n_x = 700, n_y = 31, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept[1], 
    x_axis_text = TRUE, y_axis_text = FALSE, 
    plot_title = "Caranx melampygus", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[1, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250, 500)) +
    geom_point(aes(x = 100, y = 95), alpha = 0)
# Benthic invertivores
paarca <-
mk_multipanel_plots2(fg_point_df = b, spp_point_df = b_spp_dfs$PA.ARCA, 
    spp_line_df_row = spp_lines[7, ], eqn_df = spp_sma_eqns[7, ],
    eqn_x = 350, eqn_y = 4.6, r2_x = 350, r2_y = 6.8, 
    n_x = 350, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[2], 
    x_axis_text = TRUE, y_axis_text = TRUE,
    plot_title = "Paracirrhites arcatus", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[2, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) +
    geom_point(aes(x = 100, y = 7), alpha = 0)
painsu <-
mk_multipanel_plots2(fg_point_df = b, spp_point_df = b_spp_dfs$PA.INSU, 
    spp_line_df_row = spp_lines[9, ], eqn_df = spp_sma_eqns[9, ], 
    eqn_x = 350, eqn_y = 4.6, r2_x = 350, r2_y = 6.8, 
    n_x = 350, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept[2], x_axis_text = TRUE, 
    y_axis_text = FALSE, 
    plot_title = "Parupeneus insularis", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[2, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) +
    geom_point(aes(x = 100, y = 7), alpha = 0)
mogran <-
mk_multipanel_plots2(fg_point_df = b, spp_point_df = b_spp_dfs$MO.GRAN, 
    spp_line_df_row = spp_lines[8, ], eqn_df = spp_sma_eqns[8, ], 
    eqn_x = 350, eqn_y = 4.6, r2_x = 350, r2_y = 6.8, 
    n_x = 350, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept[2], 
    x_axis_text = TRUE, y_axis_text = FALSE, 
    plot_title = "Monotaxis grandoculis", gape_dim = 'gw' ) +
    geom_abline(data = all_fg_GW_graph_df[2, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) +
    geom_point(aes(x = 100, y = 7), alpha = 0)
# Zooplanktivores
psbart <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.BART, 
    spp_line_df_row = spp_lines[13, ], eqn_df = spp_sma_eqns[13, ], 
    eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 2.625, 
    n_x = 290, n_y = 3.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[3], 
    x_axis_text = FALSE, y_axis_text = TRUE, 
    plot_title = "Pseudanthias bartlettorum", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[3, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250))
catere <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$CA.TERE, 
    spp_line_df_row = spp_lines[10, ], eqn_df = spp_sma_eqns[10, ], 
    eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 2.625, 
    n_x = 290, n_y = 3.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[3], 
    x_axis_text = FALSE, y_axis_text = FALSE, 
    plot_title = "Caesio teres", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[3, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250))
psdisp <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.DISP, 
    spp_line_df_row = spp_lines[14, ], eqn_df = spp_sma_eqns[14, ], 
    eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 2.625, 
    n_x = 290, n_y = 3.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[3], 
    x_axis_text = TRUE, y_axis_text = TRUE, 
    plot_title = "Pseudanthias dispar", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[3, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250))
psoliv <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.OLIV, 
    spp_line_df_row = spp_lines[15, ], eqn_df = spp_sma_eqns[15, ], 
    eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 2.625, 
    n_x = 290, n_y = 3.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[3], 
    x_axis_text = FALSE, y_axis_text = FALSE, 
    plot_title = "Pseudanthias olivaceus", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[3, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250))
pttile <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$PT.TILE, 
    spp_line_df_row = spp_lines[11, ], eqn_df = spp_sma_eqns[11, ], 
    eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 2.625, 
    n_x = 290, n_y = 3.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[3], 
    x_axis_text = TRUE, y_axis_text = FALSE, 
    plot_title = "Pterocaesio tile", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[3, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250))
chvand <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$CH.VAND, 
    spp_line_df_row = spp_lines[12, ], eqn_df = spp_sma_eqns[12, ], 
    eqn_x = 290, eqn_y = 1.2, r2_x = 290, r2_y = 1.75, 
    n_x = 290, n_y = 2.4, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[3], 
    x_axis_text = TRUE, y_axis_text = FALSE, 
    plot_title = "Chromis vanderbilti", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[3, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250))

# Plotting multipanel piscivores and benthic invertivores   
#dev.new(height = 10, width = 7)
#master_layout <- 
#grid.layout(nrow = 8, ncol = 4, 
#            widths = unit(c(0.2, 1, 0.9, 0.9), "null"),
#            heights = unit(c(1, 1, 0.1, 1, 0.1, 1, 1, 0.1), "null"))
#grid.newpage()

# With "Figure" labelled
master_layout <- 
grid.layout(nrow = 8, ncol = 4, 
            widths = unit(c(0.2, 1, 0.9, 0.9), "null"),
            heights = unit(c(1, 1.1, 0.2, 1.1, 0.2, 1, 1.1, 0.2), "null"))
dev.new(height = 10, width = 7)

grid.newpage()
pushViewport(viewport(layout = master_layout))
# piscs
print(apfurc, vp = set_vp(1, 2))
print(luboha, vp = set_vp(1, 3))
print(valout, vp = set_vp(1, 4))
print(ceargu, vp = set_vp(2, 2))
print(camela, vp = set_vp(2, 3))
print(ceurod, vp = set_vp(2, 4))
# benths
print(paarca, vp = set_vp(4, 2))
print(painsu, vp = set_vp(4, 3))
print(mogran, vp = set_vp(4, 4))
# zoops
print(psbart, vp = set_vp(6, 2))
print(catere, vp = set_vp(6, 3))
print(psoliv, vp = set_vp(6, 4))
print(psdisp, vp = set_vp(7, 2))
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
    expression( paste("Gape width (", mm, ")", sep = "") ), 
    vp = viewport(layout.pos.row = 1:2, layout.pos.col = 1),
    rot = 90, gp = gpar(fontsize = 9), 
    vjust = 1
    )
grid.text(
    "Standard length (mm)",
    vp = viewport(layout.pos.row = 3, layout.pos.col = 3),
    vjust = -1.2, gp = gpar(fontsize = 9)
    )
# benths
grid.text(
    "b)", vp = viewport(layout.pos.row = 4, layout.pos.col = 1), 
    gp = gpar(fontsize = 9), vjust = -10
    )
grid.text(
    expression( paste("Gape width (", mm, ")", sep = "") ), 
    vp = viewport(layout.pos.row = 4, layout.pos.col = 1),
    rot = 90, gp = gpar(fontsize = 9), 
    vjust = 1
    )
grid.text(
    "Standard length (mm)",
    vp = viewport(layout.pos.row = 5, layout.pos.col = 3),
    vjust = -1.2, gp = gpar(fontsize = 9)
    )
# zoops
grid.text(
    "c)", vp = viewport(layout.pos.row = 6, layout.pos.col = 1), 
    gp = gpar(fontsize = 9), vjust = -8
    )
grid.text(
    expression( paste("Gape width (", mm, ")", sep = "") ), 
    vp = viewport(layout.pos.row = 6:7, layout.pos.col = 1),
    rot = 90, gp = gpar(fontsize = 9), 
    vjust = 1
    )
grid.text(
    "Standard length (mm)",
    vp = viewport(layout.pos.row = 8, layout.pos.col = 3),
    vjust = -1.2, gp = gpar(fontsize = 9)
    )

dev.copy2eps(device = quartz, file = "panel_plots/gw_benth_pisc_zoop_panel_figure_label.eps")

#-------------------------------------------------------------------------------
# Multipanel herbivores
#-------------------------------------------------------------------------------
acnigr <-
mk_multipanel_plots2(fg_point_df = h, spp_point_df = h_spp_dfs$AC.NIGR, 
    spp_line_df_row = spp_lines[16, ], eqn_df = spp_sma_eqns[16, ], 
    eqn_x = 440, eqn_y = 2.5, r2_x = 440, r2_y = 3.8,
    n_x = 440, n_y = 5.1, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[4], 
    x_axis_text = FALSE, y_axis_text = TRUE, 
    plot_title = "Acanthurus nigricans", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[4, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    geom_point(aes(x = 40, y = 5), alpha = 0)
ceflav <-
mk_multipanel_plots2(fg_point_df = h, spp_point_df = h_spp_dfs$CE.FLAV, 
    spp_line_df_row = spp_lines[18, ], eqn_df = spp_sma_eqns[18, ], 
    eqn_x = 440, eqn_y = 2.5, r2_x = 440, r2_y = 3.8,
    n_x = 440, n_y = 5.1, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[4], 
    x_axis_text = TRUE, y_axis_text = FALSE, 
    plot_title = "Centropyge flavissima", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[4, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    geom_point(aes(x = 40, y = 5), alpha = 0)
chsord <-
mk_multipanel_plots2(fg_point_df  = h, spp_point_df  = h_spp_dfs$CH.SORD, 
    spp_line_df_row = spp_lines[19, ], eqn_df = spp_sma_eqns[19, ], 
    eqn_x = 440, eqn_y = 2.5, r2_x = 440, r2_y = 3.8,
    n_x = 440, n_y = 5.1, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[4], 
    x_axis_text = FALSE, y_axis_text = FALSE, 
    plot_title = "Chlororus sordidus", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[4, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    geom_point(aes(x = 40, y = 5), alpha = 0)
scrubr <-
mk_multipanel_plots2(fg_point_df  = h, spp_point_df  = h_spp_dfs$SC.RUBR, 
    spp_line_df_row = spp_lines[21, ], eqn_df = spp_sma_eqns[21, ], 
    eqn_x = 440, eqn_y = 2.5, r2_x = 440, r2_y = 3.8,
    n_x = 440, n_y = 5.1, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[4], 
    x_axis_text = TRUE, y_axis_text = TRUE, 
    plot_title = "Scarus rubroviolaceus", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[4, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    geom_point(aes(x = 40, y = 5), alpha = 0)
acoliv <-
mk_multipanel_plots2(fg_point_df = h, spp_point_df = h_spp_dfs$AC.OLIV, 
    spp_line_df_row = spp_lines[17, ], eqn_df = spp_sma_eqns[17, ], 
    eqn_x = 440, eqn_y = 2.5, r2_x = 440, r2_y = 3.8,
    n_x = 440, n_y = 5.1, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[4], 
    x_axis_text = FALSE, y_axis_text = FALSE, 
    plot_title = "Acanthurus olivaceus", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[4, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    geom_point(aes(x = 40, y = 5), alpha = 0)
scfren <-
mk_multipanel_plots2(fg_point_df  = h, spp_point_df  = h_spp_dfs$SC.FREN, 
    spp_line_df_row = spp_lines[20, ], eqn_df = spp_sma_eqns[20, ], 
    eqn_x = 440, eqn_y = 2.5, r2_x = 440, r2_y = 3.8,
    n_x = 440, n_y = 5.1, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GW_graph_df$ref_intercept_iso[4], 
    x_axis_text = TRUE, y_axis_text = FALSE, 
    plot_title = "Scarus frenatus", gape_dim = 'gw') +
    geom_abline(data = all_fg_GW_graph_df[4, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    geom_point(aes(x = 40, y = 5), alpha = 0)

# With "Figure" label
dev.new(height = 4, width = 7)
master_layout <- 
grid.layout(nrow = 3, ncol = 4, 
            widths = unit(c(0.2, 1, 0.9, 0.9), "null"),
            heights = unit(c(1, 1.1, 0.2), "null"))
grid.newpage()

pushViewport(viewport(layout = master_layout))
print(acnigr, vp = set_vp(1, 2))
print(acoliv, vp = set_vp(1, 3))
print(chsord, vp = set_vp(1, 4))
print(scrubr, vp = set_vp(2, 2))
print(scfren, vp = set_vp(2, 3))
print(ceflav, vp = set_vp(2, 4))

# Figure label
grid.text("Figure 5", vp = viewport(layout.pos.row = 3, layout.pos.col = 1), 
    gp = gpar(fontsize = 9), hjust = -1, vjust = 1)

grid.text(
    expression( paste("Gape width (", mm, ")", sep = "") ), 
    vp = viewport(layout.pos.row = 1:2, layout.pos.col = 1),
    rot = 90, gp = gpar(fontsize = 9), 
    vjust = 2
    )
grid.text(
    "Standard length (mm)",
    vp = viewport(layout.pos.row = 3, layout.pos.col = 3),
    vjust = -1.2, gp = gpar(fontsize = 9)
    )

dev.copy2eps(device = quartz, file = "panel_plots/gw_herb_panel_figure_label.eps")

#===============================================================================
# Relative gape size
#===============================================================================

# Pento factored by functional group then slope
SpeciesCode <- c("AP.FURC", "LU.BOHA", "VA.LOUT", "CA.MELA", "CE.ARGU", "CE.UROD",
                 "PA.ARCA", "PA.INSU", "MO.GRAN",
                 "PS.BART", "CA.TERE", "PS.OLIV", "PS.DISP", "PT.TILE", "CH.VAND",
                 "AC.NIGR", "AC.OLIV",  "CH.SORD", "SC.RUBR", "SC.FREN", "CE.FLAV",
                 "CH.ORNA"
                 )

sp_name_by_slope <- 
    c("Aphareus furca", "Lutjanus bohar", "Variola louti", "Caranx melampygus", 
      "Cephalopholis argus", "Cephalopholis urodeta", 
      "Paracirrhites arcatus", "Parupeneus insularis", "Monotaxis grandoculis", 
      "Pseudanthias bartlettorum", "Caesio teres", "Pseudanthias olivaceus", "Pseudanthias dispar", "Pterocaesio tile", "Chromis vanderbilti", 
      "Acanthurus nigricans", "Acanthurus olivaceus", "Chlororus sordidus", "Scarus rubroviolaceus", "Scarus frenatus", "Centropyge flavissima", 
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

ceargu_out <- which(ceargu_df$gw_ratio == max(ceargu_df$gw_ratio))
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

rel_gw <- 
ggplot(pento_by_slope, aes(.id, value, dodge = j_fg)) +
  geom_boxplot(aes(x=sp_name_by_slope, y=gw/(SL))) +
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
  theme(legend.text = element_text(size = 8)) +
  theme(panel.margin = unit(0.5, "cm")) +
  theme(panel.border = element_blank()) +
  theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(axis.line = element_line()) +
  facet_grid(. ~ j_fg, space = "free", scales = "free", labeller = labeller(fgs = fg_labeller))
dev.new(height = 5, width = 8)
rel_gw

# With "Figure" label
master_layout <- 
grid.layout(nrow = 2, ncol = 1, 
            widths = unit(c(1), "null"),
            heights = unit(c(1, 0.05), "null"))
grid.newpage()
pushViewport(viewport(layout = master_layout))
print(rel_gw, vp = set_vp(1, 1))

# Figure label
#grid.text("Figure 3", vp = viewport(layout.pos.row = 2, layout.pos.col = 1), 
 #   gp = gpar(fontsize = 9), hjust = 8, vjust = -1)
dev.copy2eps(device = quartz, file = "panel_plots/rel_gw_figure_label.eps")


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
  theme(axis.ticks.length = unit(-0.2, "cm"),
        axis.text.y = element_text(margin = margin(0, 8, 0, 0)), 
        axis.text.x = element_text(margin = margin(8, 0, 0, 0), vjust = 1))

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
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.50, 0.90), method = "rq", 
                colour = "black") +
  theme(axis.ticks.length = unit(-0.2, "cm"), 
        axis.text.y = element_text(margin = margin(0, 8, 0, 0)), 
        axis.text.x = element_text(margin = margin(8, 0, 0, 0), vjust = 1))

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
  gp = gpar(fontsize = 10), rot = 90, vjust = 1.7
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

