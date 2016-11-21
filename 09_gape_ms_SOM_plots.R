#data = p
#site = random effect
#observer = random effect
#gape = predicted
#standard_length = fixed

#===============================================================================
# Region and observer effect visualizations
#===============================================================================
# Regions representing the gradiant in fishing pressure and oceanographic 
# productivity were used instead of site as the groupings because of limited 
# sample sizes making effects hard to test given low sample sizes at some sites.

# Note: levels/factors are dropped when sample size < 3 using sma()

#-------------------------------------------------------------------------------
# Piscivore REGION effects
p_site <- sma(gh~SL * Region, data = p, log = "xy", method = "SMA", robust = T,
			  slope.test = 2, multcomp = T, multcompmethod = "adjusted")

p_site_summ <- mk_spp_summary(p_site, length(p_site$groups), grouping = TRUE)
p_site_graph_df <- mk_smaSPP_graph_df(p_site_summ, length(p_site$groups), "Region")

p_reg_plot <- 
	mk_SMAplot(df_points = p, df_lines = p_site_graph_df, facets = FALSE, 
			   x = "SL", gapeType = "gh", grouping = "Region", labels = "None", 
			   axis_labels = FALSE) +
	theme(axis.title = element_text(size = 10)) + 
	theme(axis.title.x = element_text(vjust = -0.5)) +
	theme(axis.title.y = element_text(vjust = 0.4))


# A uniroot error occurs with the following code:
#> p_obs <- sma(ga~SL * dissected_by, data = p, log = "xy", method = "SMA", robust = T,
#+   slope.test = 2, multcomp = T, multcompmethod = "adjusted")
#Warning: dropped level of grouping variable (sample size < 3) : dissected_by  =  GAJ SC/MW 
#Error in uniroot(lr.b.com, c(b, b.p), tol = 1e-04, arguments = arguments) : 
#  f() values at end points not of opposite sign
# I think this is because there are too many observers that do not have enough
# observations. After testing different sets of p_fewer_dis, I have found that 
# sma seems to produce a uniroot error when there is more than one observer with 
# only 3 observations. Therefore I have removed all observers with < 4 
# observations prior to running the sma().
ddply(p, .(dissected_by), summarise, length(dissected_by))

p_fewer_dis <- 
	p[p$dissected_by %in% c("", "AB", "LW", "RT", "RT/SC", "SC"), ]

p_obs <- sma(gh~SL * dissected_by, data = p_fewer_dis, log = "xy", method = "SMA", robust = T,
			 slope.test = 2, multcomp = T, multcompmethod = "adjusted")

p_obs_summ <- mk_spp_summary(p_obs, length(p_obs$groups), grouping = TRUE)
p_obs_graph_df <- mk_smaSPP_graph_df(p_obs_summ, length(p_obs$groups), "dissected_by")

p_obs_plot <- 
	mk_SMAplot(df_points = p_fewer_dis, df_lines = p_obs_graph_df, facets = FALSE, 
			   x = "SL", gapeType = "gh", grouping = "dissected_by", labels = "None", 
			   axis_labels = FALSE) +
	scale_colour_discrete("Observer", 
		labels = paste(" ", as.character(seq(1:length(p_obs_graph_df$dissected_by)) )), 
		guide = guide_legend(label.hjust = 1)
	) +
	theme(axis.title = element_text(size = 10)) + 
	theme(axis.title.x = element_text(vjust = -0.5)) +
	theme(axis.title.y = element_text(vjust = 0.4)) + 
	scale_x_log10(breaks=number_ticks(3)) +
  scale_y_log10(breaks=number_ticks(3))

dev.new(height = 3.5, width = 8)
master_layout <- 
grid.layout(nrow = 3, ncol = 5, 
			widths = unit(c(0.05, 1.1, 0.05, 0.05, 1.1), "null"), 
			heights = unit(c(0.1, 1, 0.1), "null")
			)
grid.newpage()
pushViewport(viewport(layout = master_layout))
print(p_reg_plot, vp = set_vp(2, 2))
print(p_obs_plot, vp = set_vp(2, 5))
grid.text(
	"a)", 
	vp = viewport(layout.pos.row = 2, layout.pos.col = 1), 
	vjust = -14, gp = gpar(fontsize = 10))
grid.text(
	"b)", 
	vp = viewport(layout.pos.row = 2, layout.pos.col = 4), 
	vjust = -14, gp = gpar(fontsize = 10))


dev.copy2eps(device = quartz, file = "panel_plots/SOM_pisc.eps")
#-------------------------------------------------------------------------------
# Benthic Invertivore
# REGION effects
b_site <- sma(gh~SL * Region, data = b, log = "xy", method = "SMA", robust = T,
			  slope.test = 2, multcomp = T, multcompmethod = "adjusted")
b_site_summ <- mk_spp_summary(b_site, length(b_site$groups), grouping = TRUE)
b_site_graph_df <- mk_smaSPP_graph_df(b_site_summ, length(b_site$groups), "Region")

b_reg_plot <- 
	mk_SMAplot(df_points = b, df_lines = b_site_graph_df, facets = FALSE, 
			   x = "SL", gapeType = "gh", grouping = "Region", labels = "None", 
			   axis_labels = FALSE) +
	theme(axis.title = element_text(size = 10)) + 
	theme(axis.title.x = element_text(vjust = -0.5)) +
	theme(axis.title.y = element_text(vjust = 0.4))

# OBSERVER effects
# include only values that were included in the SMA, so I have to manually 
# remove the points with < 3 observations.
b_fewer_dis <- 
	b[b$dissected_by %in% c("", "AB", "JB", "JD", "KP", "LW", "MW", "RT", 
							"RT/SC", "SC"), ]

b_obs <- sma(gh~SL * dissected_by, data = b_fewer_dis, log = "xy", method = "SMA", 
			 robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
b_obs_summ <- mk_spp_summary(b_obs, length(b_obs$groups), grouping = TRUE)
b_obs_graph_df <- mk_smaSPP_graph_df(b_obs_summ, length(b_obs$groups), "dissected_by")

b_obs_plot <- 
	mk_SMAplot(df_points = b_fewer_dis, df_lines = b_obs_graph_df, facets = FALSE, 
			   x = "SL", gapeType = "gh", grouping = "dissected_by", labels = "None", 
			   axis_labels = FALSE) +
	scale_colour_discrete("Observer", 
		labels = as.character(seq(1:length(b_obs_graph_df$dissected_by)))
	) + 
	theme(axis.title = element_text(size = 10)) + 
	theme(axis.title.x = element_text(vjust = -0.5)) +
	theme(axis.title.y = element_text(vjust = 0.4)) +
	guides(col=guide_legend(ncol=2))

dev.new(height = 3.5, width = 8)
master_layout <- 
grid.layout(nrow = 3, ncol = 5, 
			widths = unit(c(0.05, 1.1, 0.05, 0.05, 1.1), "null"), 
			heights = unit(c(0.1, 1, 0.1), "null")
			)
grid.newpage()
pushViewport(viewport(layout = master_layout))
print(b_reg_plot, vp = set_vp(2, 2))
print(b_obs_plot, vp = set_vp(2, 5))
grid.text(
	"a)", 
	vp = viewport(layout.pos.row = 2, layout.pos.col = 1), 
	vjust = -14, gp = gpar(fontsize = 10))
grid.text(
	"b)", 
	vp = viewport(layout.pos.row = 2, layout.pos.col = 4), 
	vjust = -14, gp = gpar(fontsize = 10))


dev.copy2eps(device = quartz, file = "panel_plots/SOM_benth.eps")

#-------------------------------------------------------------------------------
# Zooplanktivore
# REGION effects
ddply(z, .(dissected_by), summarise, length(dissected_by))
z_site <- sma(gh~SL * Region, data = z, log = "xy", method = "SMA", robust = T,
			  slope.test = 2, multcomp = T, multcompmethod = "adjusted")
z_site_summ <- mk_spp_summary(z_site, length(z_site$groups), grouping = TRUE)
z_site_graph_df <- mk_smaSPP_graph_df(z_site_summ, length(z_site$groups), "Region")

z_reg_plot <- 
	mk_SMAplot(df_points = z, df_lines = z_site_graph_df, facets = FALSE, 
			   x = "SL", gapeType = "gh", grouping = "Region", labels = "None", 
			   axis_labels = FALSE) +
	theme(axis.title = element_text(size = 10)) + 
	theme(axis.title.x = element_text(vjust = -0.5)) +
	theme(axis.title.y = element_text(vjust = 0.4))

# OBSERVER effects
# same uniroot error as was found with the piscivores, again have removed values with < 6 observers
z_fewer_dis <- 
	z[z$dissected_by %in% c("", "AO/KP", "KP", "KP/AO", 
						    "RT/AO", "RT/Adrian", "RT/SC", 
						    "rowan.angeleen"), ]

z_obs <- sma(gh~SL * dissected_by, data = z_fewer_dis, log = "xy", method = "SMA", 
			  robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
z_obs_summ <- mk_spp_summary(z_obs, length(z_obs$groups), grouping = TRUE)
z_obs_graph_df <- mk_smaSPP_graph_df(z_obs_summ, length(z_obs$groups), "dissected_by")

z_obs_plot <- 
	mk_SMAplot(df_points = z_fewer_dis, df_lines = z_obs_graph_df, facets = FALSE, 
			   x = "SL", gapeType = "gh", grouping = "dissected_by", labels = "None", 
			   axis_labels = FALSE) +
	scale_colour_discrete("Observer", 
		labels = as.character(seq(1:length(z_obs_graph_df$dissected_by)))
	) + 
	theme(axis.title = element_text(size = 10)) + 
	theme(axis.title.x = element_text(vjust = -0.5)) +
	theme(axis.title.y = element_text(vjust = 0.4)) +
	guides(col=guide_legend(ncol=2))

dev.new(height = 3.5, width = 8)
master_layout <- 
grid.layout(nrow = 3, ncol = 5, 
			widths = unit(c(0.05, 1.1, 0.05, 0.05, 1.1), "null"), 
			heights = unit(c(0.1, 1, 0.1), "null")
			)
grid.newpage()
pushViewport(viewport(layout = master_layout))
print(z_reg_plot, vp = set_vp(2, 2))
print(z_obs_plot, vp = set_vp(2, 5))
grid.text(
	"a)", 
	vp = viewport(layout.pos.row = 2, layout.pos.col = 1), 
	vjust = -14, gp = gpar(fontsize = 10))
grid.text(
	"b)", 
	vp = viewport(layout.pos.row = 2, layout.pos.col = 4), 
	vjust = -14, gp = gpar(fontsize = 10))

dev.copy2eps(device = quartz, file = "panel_plots/SOM_zoop.eps")



#-------------------------------------------------------------------------------
# Herbivore
# REGION effects
h_site <- sma(gh~SL * Region, data = h, log = "xy", method = "SMA", robust = T,
			  slope.test = 2, multcomp = T, multcompmethod = "adjusted")
h_site_summ <- mk_spp_summary(h_site, length(h_site$groups), grouping = TRUE)
h_site_graph_df <- mk_smaSPP_graph_df(h_site_summ, length(h_site$groups), "Region")

h_reg_plot <- 
	mk_SMAplot(df_points = h, df_lines = h_site_graph_df, facets = FALSE, 
			   x = "SL", gapeType = "gh", grouping = "Region", labels = "None", 
			   axis_labels = FALSE) +
	theme(axis.title = element_text(size = 10)) + 
	theme(axis.title.x = element_text(vjust = -0.5)) +
	theme(axis.title.y = element_text(vjust = 0.4))

# OBSERVER effects
ddply(h, .(dissected_by), summarise, length(dissected_by))

h_fewer_dis <- 
	h[h$dissected_by %in% c("", "AB", "GAJ", "LW", "RT", "RT/SC", "SC"), ] 

h_obs <- sma(gh~SL * dissected_by, data = h_fewer_dis, log = "xy", method = "SMA", 
			  robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
h_obs_summ <- mk_spp_summary(h_obs, length(h_obs$groups), grouping = TRUE)
h_obs_graph_df <- mk_smaSPP_graph_df(h_obs_summ, length(h_obs$groups), "dissected_by")

h_obs_plot <- 
	mk_SMAplot(df_points = h_fewer_dis, df_lines = h_obs_graph_df, facets = FALSE, 
			   x = "SL", gapeType = "gh", grouping = "dissected_by", labels = "None", 
			   axis_labels = FALSE) +
	scale_colour_discrete("Observer", 
		labels = as.character(seq(1:length(h_obs_graph_df$dissected_by)))
	) + 
	theme(axis.title = element_text(size = 10)) + 
	theme(axis.title.x = element_text(vjust = -0.5)) +
	theme(axis.title.y = element_text(vjust = 0.4)) +
	guides(col=guide_legend(ncol=2))

dev.new(height = 3.5, width = 8)
master_layout <- 
grid.layout(nrow = 3, ncol = 5, 
			widths = unit(c(0.05, 1.1, 0.05, 0.05, 1.1), "null"), 
			heights = unit(c(0.1, 1, 0.1), "null")
			)
grid.newpage()
pushViewport(viewport(layout = master_layout))
print(h_reg_plot, vp = set_vp(2, 2))
print(h_obs_plot, vp = set_vp(2, 5))
grid.text(
	"a)", 
	vp = viewport(layout.pos.row = 2, layout.pos.col = 1), 
	vjust = -14, gp = gpar(fontsize = 10))
grid.text(
	"b)", 
	vp = viewport(layout.pos.row = 2, layout.pos.col = 4), 
	vjust = -14, gp = gpar(fontsize = 10))

dev.copy2eps(device = quartz, file = "panel_plots/SOM_herb.eps")


# Corallivore REGION effects
c_site <- sma(gh~SL * Region, data = c, log = "xy", method = "SMA", robust = T,
			  slope.test = 2, multcomp = T, multcompmethod = "adjusted")
c_site_summ <- mk_spp_summary(c_site, length(c_site$groups), grouping = TRUE)
c_site_graph_df <- mk_smaSPP_graph_df(c_site_summ, length(c_site$groups), "Region")

c_reg_plot <- 
	mk_SMAplot(df_points = c, df_lines = c_site_graph_df, facets = FALSE, 
			   x = "SL", gapeType = "gh", grouping = "Region", labels = "None", 
			   axis_labels = FALSE) +
	theme(axis.title = element_text(size = 10)) + 
	theme(axis.title.x = element_text(vjust = -0.5)) +
	theme(axis.title.y = element_text(vjust = 0.4))

ddply(c, .(dissected_by), summarise, length(dissected_by))


c_fewer_dis <- 
	c[c$dissected_by %in% c("AB", "LW", "MW", "RT", "SC"), ] 

c_obs <- sma(gh~SL * dissected_by, data = c_fewer_dis, log = "xy", method = "SMA", 
			 robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")

c_obs_summ <- mk_spp_summary(c_obs, length(c_obs$groups), grouping = TRUE)
c_obs_graph_df <- mk_smaSPP_graph_df(c_obs_summ, length(c_obs$groups), "dissected_by")

c_obs_plot <- 
	mk_SMAplot(df_points = c_fewer_dis, df_lines = c_obs_graph_df, facets = FALSE, 
			   x = "SL", gapeType = "gh", grouping = "dissected_by", labels = "None", 
			   axis_labels = FALSE) +
	scale_colour_discrete("Observer", 
		labels = as.character(seq(1:length(c_obs_graph_df$dissected_by)))
	) + 
	theme(axis.title = element_text(size = 10)) + 
	theme(axis.title.x = element_text(vjust = -0.5)) +
	theme(axis.title.y = element_text(vjust = 0.4)) +
	guides(col=guide_legend(ncol=2))

dev.new(height = 3.5, width = 8)
master_layout <- 
grid.layout(nrow = 3, ncol = 5, 
			widths = unit(c(0.05, 1.1, 0.05, 0.05, 1.1), "null"), 
			heights = unit(c(0.1, 1, 0.1), "null")
			)
grid.newpage()
pushViewport(viewport(layout = master_layout))
print(c_reg_plot, vp = set_vp(2, 2))
print(c_obs_plot, vp = set_vp(2, 5))
grid.text(
	"a)", 
	vp = viewport(layout.pos.row = 2, layout.pos.col = 1), 
	vjust = -14, gp = gpar(fontsize = 10))
grid.text(
	"b)", 
	vp = viewport(layout.pos.row = 2, layout.pos.col = 4), 
	vjust = -14, gp = gpar(fontsize = 10))

dev.copy2eps(device = quartz, file = "panel_plots/SOM_coral.eps")