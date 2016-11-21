#===============================================================================
# Multipanel plot with functional group and species
#===============================================================================

# Get community level midpoint
allGH <- sma(gh ~ SL, data = pento, log = "xy", method = "SMA", robust = T, slope.test = 2)
#check_assump(allGH, "Pento Gape Area All")
allGH_summ <- mk_sma_summary(allGH, 1)
allGH_graph_df <- mk_sma_graph_df(allGH_summ, 1, "j_fg")
names(allGH_graph_df)


# Get functional group level metrics: midpoints of data
# Using a single sma() iteration because this does not alter results, nor the 
# visual display of the data.
pento_slopes_gh <- coef(fit.gls_gh)
names(pento_slopes_gh) <- c('Pi', 'Bi', 'Zp', 'He')

all_fg_GH <- sma(gh ~ SL * j_fg, data = pento, log = "xy", method = "SMA", 
                 robust = T, slope.test = 1, multcomp = F)
all_fg_GH_summ <- mk_spp_summary(all_fg_GH, 5, grouping=TRUE)
all_fg_GH_graph_df <- mk_smaSPP_graph_df(all_fg_GH_summ, 5, "j_fg")
all_fg_GH_graph_df$boot_slope <- c(pento_slopes_gh, all_fg_GH_graph_df$slp[5])

all_fg_GH_graph_df <- 
  all_fg_GH_graph_df %>% 
    dplyr::mutate(boot_ref_int = log10(midpoint_y / (midpoint_x ^ boot_slope)))

# SMA regression for each species
all_spp_GH <- sma(gh ~ SL * SpeciesCode, data = pento, log = "xy", 
                  method = "SMA", robust = T, slope.test = 1,
                  multcomp = F, multcompmethod = "adjusted")
#check_assump(allGA, "Species Gape Area All")
allGH_bySPP_summ <- mk_spp_summary(all_spp_GH, 22, grouping=TRUE)

allometry <- NA
for (i in seq_along(allGH_bySPP_summ[[1]])) {
  allometry[i] <- get_allometry(slope = allGH_bySPP_summ$slope[i], 
                                p_val = allGH_bySPP_summ$slp_p_val[i])
}

sig <- NA
for (i in seq_along(allGH_bySPP_summ[[1]])) {
  sig[i] <- get_sig(p_val = allGH_bySPP_summ$slp_p_val[i])
}

allGH_bySPP_summ$allometry <- allometry
allGH_bySPP_summ$sig <- '*'
allGH_bySPP_summ <- 
  mutate(allGH_bySPP_summ, allometry = replace(allometry, group == 'CH.ORNA', ''), 
         sig = replace(sig, group == 'CH.ORNA', ''))

allGH_bySPP_graph_df <- mk_smaSPP_graph_df(allGH_bySPP_summ, 22, "SpeciesCode")
allGH_bySPP_graph_df$allometry <- allometry


#-------------------------------------------------------------------------------
# Setting up values to plot the lines at the species level
spp_lines <- allGH_bySPP_graph_df

# Setting up equation, r^2, and n values that will be written on the graphs
spp_sma_eqns <- write_group_sma_eqn(allGH_bySPP_summ, allGH_bySPP_summ$group)
names(spp_sma_eqns) <- c("SpeciesCode", "eqn_r2", "eqn", "r2", "n")



#-------------------------------------------------------------------------------
# Multipanel comparison of 3 predatory functional groups
#-------------------------------------------------------------------------------
apfurc <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$AP.FURC, 
    spp_line_df_row = spp_lines[2, ], eqn_df = spp_sma_eqns[2, ], 
    eqn_x = 700, eqn_y = 20, r2_x = 700, r2_y = 25.8,
    n_x = 700, n_y = 31, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[1], 
    x_axis_text = FALSE, y_axis_text = FALSE, plot_title = "Aphareus furca") +
    geom_abline(data = all_fg_GH_graph_df[1, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250, 500)) +
    scale_y_log10(breaks = c(20, 50, 100)) +
    geom_point(aes(x = 100, y = 95), alpha = 0) + 
    geom_text(data = allGH_bySPP_graph_df[2, ], aes_string(x = 80, y = 100, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
luboha <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$LU.BOHA, 
    spp_line_df_row = spp_lines[3, ], eqn_df = spp_sma_eqns[3, ], 
    eqn_x = 700, eqn_y = 20, r2_x = 700, r2_y = 25.8, 
    n_x = 700, n_y = 31, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept[1], 
    x_axis_text = FALSE, y_axis_text = FALSE, plot_title = "Lutjanus bohar") +
    geom_abline(data = all_fg_GH_graph_df[1, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250, 500)) +
    scale_y_log10(breaks = c(20, 50, 100)) +
    geom_point(aes(x = 100, y = 95), alpha = 0) +
    geom_text(data = allGH_bySPP_graph_df[3, ], aes_string(x = 80, y = 100, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
valout <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$VA.LOUT, 
    spp_line_df_row = spp_lines[6, ], eqn_df = spp_sma_eqns[6, ], 
    eqn_x = 700, eqn_y = 20, r2_x = 700, r2_y = 25.8, 
    n_x = 700, n_y = 31, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept[1], 
    x_axis_text = FALSE, y_axis_text = TRUE, plot_title = "Variola louti") +
    geom_abline(data = all_fg_GH_graph_df[1, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250, 500)) +
    scale_y_log10(breaks = c(20, 50, 100)) +
    geom_point(aes(x = 100, y = 95), alpha = 0) + 
    geom_text(data = allGH_bySPP_graph_df[6, ], aes_string(x = 80, y = 100, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
ceargu <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CE.ARGU, 
    spp_line_df_row = spp_lines[4, ], eqn_df = spp_sma_eqns[4, ], 
    eqn_x = 700, eqn_y = 20, r2_x = 700, r2_y = 25.8, 
    n_x = 700, n_y = 31, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept[1], 
    x_axis_text = TRUE, y_axis_text = TRUE, plot_title = "Cephalopholis argus") +
    geom_abline(data = all_fg_GH_graph_df[1, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250, 500)) +
    scale_y_log10(breaks = c(20, 50, 100)) +
    geom_point(aes(x = 100, y = 95), alpha = 0) + 
    geom_text(data = allGH_bySPP_graph_df[4, ], aes_string(x = 80, y = 100, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
ceurod <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CE.UROD, 
    spp_line_df_row = spp_lines[5, ], eqn_df = spp_sma_eqns[5, ], 
    eqn_x = 700, eqn_y = 20, r2_x = 700, r2_y = 25.8, 
    n_x = 700, n_y = 31, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept[1], 
    x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Cephalopholis urodeta") +
    geom_abline(data = all_fg_GH_graph_df[1, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250, 500)) +
    scale_y_log10(breaks = c(20, 50, 100)) +
    geom_point(aes(x = 100, y = 95), alpha = 0) + 
    geom_text(data = allGH_bySPP_graph_df[5, ], aes_string(x = 80, y = 100, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
camela <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CA.MELA, 
    spp_line_df_row = spp_lines[1, ], eqn_df = spp_sma_eqns[1, ], 
    eqn_x = 700, eqn_y = 20, r2_x = 700, r2_y = 25.8, 
    n_x = 700, n_y = 31, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept[1], 
    x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Caranx melampygus") +
    geom_abline(data = all_fg_GH_graph_df[1, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250, 500)) +
    scale_y_log10(breaks = c(20, 50, 100)) +
    geom_point(aes(x = 100, y = 95), alpha = 0) + 
    geom_text(data = allGH_bySPP_graph_df[1, ], aes_string(x = 80, y = 100, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
# Benthic invertivores
paarca <-
mk_multipanel_plots2(fg_point_df = b, spp_point_df = b_spp_dfs$PA.ARCA, 
    spp_line_df_row = spp_lines[7, ], eqn_df = spp_sma_eqns[7, ],
    eqn_x = 350, eqn_y = 4.6, r2_x = 350, r2_y = 6.8, 
    n_x = 350, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[2], 
    x_axis_text = TRUE, y_axis_text = TRUE, plot_title = "Paracirrhites arcatus") +
    geom_abline(data = all_fg_GH_graph_df[2, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) +
    scale_y_log10(breaks = c(5, 10, 50)) +
    geom_point(aes(x = 100, y = 7), alpha = 0) + 
    geom_text(data = allGH_bySPP_graph_df[7, ], aes_string(x = 50, y = 60, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
painsu <-
mk_multipanel_plots2(fg_point_df = b, spp_point_df = b_spp_dfs$PA.INSU, 
    spp_line_df_row = spp_lines[9, ], eqn_df = spp_sma_eqns[9, ], 
    eqn_x = 350, eqn_y = 4.6, r2_x = 350, r2_y = 6.8, 
    n_x = 350, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept[2], x_axis_text = TRUE, 
    y_axis_text = FALSE, plot_title = "Parupeneus insularis") +
    geom_abline(data = all_fg_GH_graph_df[2, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) +
    scale_y_log10(breaks = c(5, 10, 50)) +
    geom_point(aes(x = 100, y = 7), alpha = 0) + 
    geom_text(data = allGH_bySPP_graph_df[9, ], aes_string(x = 50, y = 60, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
mogran <-
mk_multipanel_plots2(fg_point_df = b, spp_point_df = b_spp_dfs$MO.GRAN, 
    spp_line_df_row = spp_lines[8, ], eqn_df = spp_sma_eqns[8, ], 
    eqn_x = 350, eqn_y = 4.6, r2_x = 350, r2_y = 6.8, 
    n_x = 350, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept[2], x_axis_text = TRUE, 
    y_axis_text = FALSE, plot_title = "Monotaxis grandoculis") +
    geom_abline(data = all_fg_GH_graph_df[2, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) +
    scale_y_log10(breaks = c(5, 10, 50)) +
    geom_point(aes(x = 100, y = 7), alpha = 0) + 
    geom_text(data = allGH_bySPP_graph_df[8, ], aes_string(x = 50, y = 60, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
# Zooplanktivores
psbart <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.BART, 
    spp_line_df_row = spp_lines[13, ], eqn_df = spp_sma_eqns[13, ], 
    eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 2.625, 
    n_x = 290, n_y = 3.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[3], 
    x_axis_text = FALSE, y_axis_text = FALSE, plot_title = "Pseudanthias bartlettorum") +
    geom_abline(data = all_fg_GH_graph_df[3, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    scale_y_log10(breaks = c(2, 5, 20)) + 
    geom_text(data = allGH_bySPP_graph_df[13, ], aes_string(x = 30, y = 25, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
catere <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$CA.TERE, 
    spp_line_df_row = spp_lines[10, ], eqn_df = spp_sma_eqns[10, ], 
    eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 2.625, 
    n_x = 290, n_y = 3.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[3], x_axis_text = FALSE, 
    y_axis_text = FALSE, plot_title = "Caesio teres") +
    geom_abline(data = all_fg_GH_graph_df[3, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    scale_y_log10(breaks = c(2, 5, 20)) + 
    geom_text(data = allGH_bySPP_graph_df[10, ], aes_string(x = 30, y = 25, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
psdisp <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.DISP, 
    spp_line_df_row = spp_lines[14, ], eqn_df = spp_sma_eqns[14, ], 
    eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 2.625, 
    n_x = 290, n_y = 3.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[3], x_axis_text = FALSE, 
    y_axis_text = TRUE, plot_title = "Pseudanthias dispar") +
    geom_abline(data = all_fg_GH_graph_df[3, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    scale_y_log10(breaks = c(2, 5, 20)) + 
    geom_text(data = allGH_bySPP_graph_df[14, ], aes_string(x = 30, y = 25, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
psoliv <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.OLIV, 
    spp_line_df_row = spp_lines[15, ], eqn_df = spp_sma_eqns[15, ], 
    eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 2.625, 
    n_x = 290, n_y = 3.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[3], x_axis_text = TRUE, 
    y_axis_text = TRUE, plot_title = "Pseudanthias olivaceus") +
    geom_abline(data = all_fg_GH_graph_df[3, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    scale_y_log10(breaks = c(2, 5, 20)) + 
    geom_text(data = allGH_bySPP_graph_df[15, ], aes_string(x = 30, y = 25, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
pttile <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$PT.TILE, 
    spp_line_df_row = spp_lines[11, ], eqn_df = spp_sma_eqns[11, ], 
    eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 2.625, 
    n_x = 290, n_y = 3.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[3], x_axis_text = TRUE, 
    y_axis_text = FALSE, plot_title = "Pterocaesio tile") +
    geom_abline(data = all_fg_GH_graph_df[3, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    scale_y_log10(breaks = c(2, 5, 20)) + 
    geom_text(data = allGH_bySPP_graph_df[11, ], aes_string(x = 30, y = 25, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
chvand <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$CH.VAND, 
    spp_line_df_row = spp_lines[12, ], eqn_df = spp_sma_eqns[12, ], 
    eqn_x = 290, eqn_y = 1.8, r2_x = 290, r2_y = 2.625, 
    n_x = 290, n_y = 3.5, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[3], x_axis_text = TRUE, 
    y_axis_text = FALSE, plot_title = "Chromis vanderbilti") +
    geom_abline(data = all_fg_GH_graph_df[3, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    scale_y_log10(breaks = c(2, 5, 20)) + 
    geom_text(data = allGH_bySPP_graph_df[12, ], aes_string(x = 30, y = 25, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)

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
print(valout, vp = set_vp(1, 2))
print(luboha, vp = set_vp(1, 3))
print(apfurc, vp = set_vp(1, 4))
print(ceargu, vp = set_vp(2, 2))
print(ceurod, vp = set_vp(2, 3))
print(camela, vp = set_vp(2, 4))
# benths
print(paarca, vp = set_vp(4, 2))
print(painsu, vp = set_vp(4, 3))
print(mogran, vp = set_vp(4, 4))
# zoops
print(psdisp, vp = set_vp(6, 2))
print(catere, vp = set_vp(6, 3))
print(psbart, vp = set_vp(6, 4))
print(psoliv, vp = set_vp(7, 2))
print(pttile, vp = set_vp(7, 3))
print(chvand, vp = set_vp(7, 4))

# Figure label
grid.text("Figure 5", vp = viewport(layout.pos.row = 8, layout.pos.col = 1),
    gp = gpar(fontsize = 9), hjust = -1, vjust = 1)

# piscs
grid.text(
    "a)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1), 
    gp = gpar(fontsize = 9), vjust = -8
    )
grid.text(
    expression( paste("Gape height (", mm, ")", sep = "") ), 
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
    expression( paste("Gape height (", mm, ")", sep = "") ), 
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
    expression( paste("Gape height (", mm, ")", sep = "") ), 
    vp = viewport(layout.pos.row = 6:7, layout.pos.col = 1),
    rot = 90, gp = gpar(fontsize = 9), 
    vjust = 1
    )
grid.text(
    "Standard length (mm)",
    vp = viewport(layout.pos.row = 8, layout.pos.col = 3),
    vjust = -1.2, gp = gpar(fontsize = 9)
    )

dev.copy2eps(device = quartz, file = "panel_plots/DunicBaum_f5.eps")

#-------------------------------------------------------------------------------
# Multipanel herbivores
#-------------------------------------------------------------------------------
acnigr <-
mk_multipanel_plots2(fg_point_df = h, spp_point_df = h_spp_dfs$AC.NIGR, 
    spp_line_df_row = spp_lines[16, ], eqn_df = spp_sma_eqns[16, ], 
    eqn_x = 440, eqn_y = 3.9, r2_x = 440, r2_y = 5.3,
    n_x = 440, n_y = 6.8, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[4], 
    x_axis_text = TRUE, y_axis_text = TRUE, plot_title = "Acanthurus nigricans") +
    geom_abline(data = all_fg_GH_graph_df[4, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    scale_y_log10(breaks = c(5, 10, 40)) +
    geom_point(aes(x = 40, y = 5), alpha = 0) + 
    geom_text(data = allGH_bySPP_graph_df[16, ], aes_string(x = 50, y = 40, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
ceflav <-
mk_multipanel_plots2(fg_point_df = h, spp_point_df = h_spp_dfs$CE.FLAV, 
    spp_line_df_row = spp_lines[18, ], eqn_df = spp_sma_eqns[18, ], 
    eqn_x = 440, eqn_y = 3.5, r2_x = 440, r2_y = 5.2,
    n_x = 440, n_y = 6.8, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[4], 
    x_axis_text = FALSE, y_axis_text = TRUE, plot_title = "Centropyge flavissima") +
    geom_abline(data = all_fg_GH_graph_df[4, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    scale_y_log10(breaks = c(5, 10, 40)) +
    geom_point(aes(x = 40, y = 5), alpha = 0) +
    geom_text(data = allGH_bySPP_graph_df[18, ], aes_string(x = 50, y = 40, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
chsord <-
mk_multipanel_plots2(fg_point_df  = h, spp_point_df  = h_spp_dfs$CH.SORD, 
    spp_line_df_row = spp_lines[19, ], eqn_df = spp_sma_eqns[19, ], 
    eqn_x = 440, eqn_y = 3.9, r2_x = 440, r2_y = 5.3,
    n_x = 440, n_y = 6.8, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[4], 
    x_axis_text = FALSE, y_axis_text = FALSE, plot_title = "Chlorurus sordidus") +
    geom_abline(data = all_fg_GH_graph_df[4, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    scale_y_log10(breaks = c(5, 10, 40)) +
    geom_point(aes(x = 40, y = 5), alpha = 0) + 
    geom_text(data = allGH_bySPP_graph_df[19, ], aes_string(x = 50, y = 40, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
scrubr <-
mk_multipanel_plots2(fg_point_df  = h, spp_point_df  = h_spp_dfs$SC.RUBR, 
    spp_line_df_row = spp_lines[21, ], eqn_df = spp_sma_eqns[21, ], 
    eqn_x = 440, eqn_y = 3.9, r2_x = 440, r2_y = 5.3,
    n_x = 440, n_y = 6.8, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[4], 
    x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Scarus rubroviolaceus") +
    geom_abline(data = all_fg_GH_graph_df[4, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    scale_y_log10(breaks = c(5, 10, 40)) +
    geom_point(aes(x = 40, y = 5), alpha = 0) + 
    geom_text(data = allGH_bySPP_graph_df[21, ], aes_string(x = 50, y = 40, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
acoliv <-
mk_multipanel_plots2(fg_point_df = h, spp_point_df = h_spp_dfs$AC.OLIV, 
    spp_line_df_row = spp_lines[17, ], eqn_df = spp_sma_eqns[17, ], 
    eqn_x = 440, eqn_y = 3.9, r2_x = 440, r2_y = 5.3,
    n_x = 440, n_y = 6.8, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[4], 
    x_axis_text = FALSE, y_axis_text = FALSE, plot_title = "Acanthurus olivaceus") +
    geom_abline(data = all_fg_GH_graph_df[4, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    scale_y_log10(breaks = c(5, 10, 40)) +
    geom_point(aes(x = 40, y = 5), alpha = 0) + 
    geom_text(data = allGH_bySPP_graph_df[17, ], aes_string(x = 50, y = 40, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)
scfren <-
mk_multipanel_plots2(fg_point_df  = h, spp_point_df  = h_spp_dfs$SC.FREN, 
    spp_line_df_row = spp_lines[20, ], eqn_df = spp_sma_eqns[20, ], 
    eqn_x = 440, eqn_y = 3.9, r2_x = 440, r2_y = 5.3,
    n_x = 440, n_y = 6.8, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[4], 
    x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Scarus frenatus") +
    geom_abline(data = all_fg_GH_graph_df[4, ], 
                aes(slope = boot_slope, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    scale_y_log10(breaks = c(5, 10, 40)) +
    geom_point(aes(x = 40, y = 5), alpha = 0) + 
    geom_text(data = allGH_bySPP_graph_df[20, ], aes_string(x = 50, y = 40, 
        label = 'allometry'), parse = TRUE, size = 3, hjust = 1)

#dev.new(height = 4, width = 7)
#master_layout <- 
#grid.layout(nrow = 3, ncol = 4, 
#            widths = unit(c(0.2, 1, 0.9, 0.9), "null"),
#            heights = unit(c(1, 1, 0.2), "null"))
#grid.newpage()

# With "Figure" label
dev.new(height = 4, width = 7)
master_layout <- 
grid.layout(nrow = 3, ncol = 4, 
            widths = unit(c(0.2, 1, 0.9, 0.9), "null"),
            heights = unit(c(1, 1.05, 0.2), "null"))
grid.newpage()


pushViewport(viewport(layout = master_layout))
print(ceflav, vp = set_vp(1, 2))
print(chsord, vp = set_vp(1, 3))
print(acoliv, vp = set_vp(1, 4))
print(acnigr, vp = set_vp(2, 2))
print(scfren, vp = set_vp(2, 3))
print(scrubr, vp = set_vp(2, 4))

# Figure label
grid.text("Figure 6", vp = viewport(layout.pos.row = 3, layout.pos.col = 1), 
    gp = gpar(fontsize = 9), hjust = -1, vjust = 1)
grid.text(
    expression( paste("Gape height (", mm, ")", sep = "") ), 
    vp = viewport(layout.pos.row = 1:2, layout.pos.col = 1),
    rot = 90, gp = gpar(fontsize = 9), 
    vjust = 2
    )
grid.text(
    "Standard length (mm)",
    vp = viewport(layout.pos.row = 3, layout.pos.col = 3),
    vjust = -1.2, gp = gpar(fontsize = 9)
    )

dev.copy2eps(device = quartz, file = "panel_plots/DunicBaum_f6.eps")


#-------------------------------------------------------------------------------
# Corallivore plot
#-------------------------------------------------------------------------------
dev.new(height = 2.4, width = 2.7)
master_layout <- 
grid.layout(nrow = 2, ncol = 2, 
            widths = unit(c(0.2, 1), "null"),
            heights = unit(c(1, 0.2), "null"))

grid.newpage()
chorna <-
mk_corallivore_plot(fg_point_df = c, spp_point_df = c, 
    eqn_df = spp_sma_eqns[22, ], 
    eqn_x = 440, eqn_y = 5.4, r2_x = 440, r2_y = 6.8,
    n_x = 440, n_y = 8.2, x_axis_labels = FALSE, y_axis_labels = FALSE, 
    fg_line_intercept = all_fg_GH_graph_df$ref_intercept_iso[5], 
    x_axis_text = TRUE, y_axis_text = TRUE, plot_title = "Chaetodon ornatissimus") +
    geom_abline(data = all_fg_GH_graph_df[5, ], 
                aes(slope = slp, intercept = boot_ref_int), linetype = 2) +
    scale_x_log10(breaks = c(50, 100, 250)) + 
    scale_y_log10(breaks = c(5, 10, 40)) +
    geom_point(aes(x = 40, y = 5), alpha = 0) +
    geom_point(aes(x = 40, y = 40), alpha = 0) 

pushViewport(viewport(layout = master_layout))
print(chorna, vp = set_vp(1, 2))
grid.text("Figure 3", vp = viewport(layout.pos.row = 2, layout.pos.col = 1), 
    gp = gpar(fontsize = 9), hjust = 0, vjust = 1)
grid.text(
    expression( paste("Gape height (", mm, ")", sep = "") ), 
    vp = viewport(layout.pos.row = 1, layout.pos.col = 1),
    rot = 90, gp = gpar(fontsize = 9), 
    vjust = 2
    )
grid.text(
    "Standard length (mm)",
    vp = viewport(layout.pos.row = 2, layout.pos.col = 2),
    vjust = -1.2, gp = gpar(fontsize = 9)
    )

dev.copy2eps(device = quartz, file = "panel_plots/DunicBaum_f3.eps")


#===============================================================================
# Relative gape size
#===============================================================================
# Pento factored by functional group then slope
SpeciesCode <- c("VA.LOUT", "LU.BOHA", "AP.FURC", "CE.ARGU", "CE.UROD", "CA.MELA",
                 "PA.ARCA", "PA.INSU", "MO.GRAN",
                 "PS.DISP", "CA.TERE", "PS.BART", "PS.OLIV", "PT.TILE", "CH.VAND",
                 "CE.FLAV", "CH.SORD", "AC.OLIV", "AC.NIGR", "SC.FREN", "SC.RUBR", 
                 "CH.ORNA"
                 )

sp_name_by_slope <- 
    c("Variola louti", "Lutjanus bohar", "Aphareus furca", "Cephalopholis argus", "Cephalopholis urodeta", "Caranx melampygus", 
      "Paracirrhites arcatus", "Parupeneus insularis", "Monotaxis grandoculis", 
      "Pseudanthias dispar", "Caesio teres", "Pseudanthias bartlettorum", 
      "Pseudanthias olivaceus", "Pterocaesio tile", "Chromis vanderbilti", 
      "Centropyge flavissima", "Chlorurus sordidus", "Acanthurus olivaceus", "Acanthurus nigricans", "Scarus frenatus", "Scarus rubroviolaceus", 
      "Chaetodon ornatissimus"
      )

spp_key <- data.frame(SpeciesCode, sp_name_by_slope)
pento_by_slope <- merge(x = pento, y = spp_key, all.x = TRUE, all.y = FALSE)
#pento_by_slope$SpeciesCode <- factor(pento_by_slope$SpeciesCode, levels = SpeciesCode)
pento_by_slope$sp_name_by_slope <- factor(pento_by_slope$sp_name_by_slope, levels = rev(sp_name_by_slope))

#-------------------------------------------------------------------------------
# Removing ceargu_out (CE.ARGU outlier)
#-------------------------------------------------------------------------------
# Finding outlier in boxplot Cephalopholis argus
ceargu_df <- fish[(which(fish$SpeciesCode == "CE.ARGU")), ]

ceargu_out <- which(ceargu_df$gh_ratio == max(ceargu_df$gh_ratio))
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

rel_gh <- 
ggplot(pento_by_slope, aes(.id, value, dodge = j_fg)) +
  geom_boxplot(aes(x=sp_name_by_slope, y=gh/(SL))) +
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
  theme(legend.text = element_text(size = 8)) +
  theme(panel.margin = unit(0.5, "cm")) +
  theme(panel.border = element_blank()) +
  theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(axis.line = element_line()) +
  facet_grid(. ~ j_fg, space = "free", scales = "free", labeller = labeller(fgs = fg_labeller))
dev.new(height = 5, width = 8)
rel_gh

# With "Figure" label
master_layout <- 
grid.layout(nrow = 2, ncol = 2, 
            widths = unit(c(0.03, 1), "null"),
            heights = unit(c(1, 0.05), "null"))
grid.newpage()
pushViewport(viewport(layout = master_layout))
print(rel_gh, vp = set_vp(1, 2))

# Figure label
grid.text("Figure 4", vp = viewport(layout.pos.row = 2, layout.pos.col = 1), 
    gp = gpar(fontsize = 9), hjust = 0, vjust = 0)

dev.copy2eps(device = quartz, file = "panel_plots/DunicBaum_f4.eps")


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
grid.text("Figure 7", vp = viewport(layout.pos.row = 2, layout.pos.col = 1), 
    gp = gpar(fontsize = 9), hjust = -1, vjust = 1)

dev.copy2eps(device = quartz, file = "panel_plots/DunicBaum_f7.eps")

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
