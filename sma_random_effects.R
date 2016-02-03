# Corallivores
cGA <- sma(ga~SL, data=c, log="xy", method="SMA", robust=T, slope.test=2)
check_assump(cGA, "Co Gape Area All")
cGA_summ <- mk_sma_summary(cGA, 1)
cGA_graph_df <- mk_sma_graph_df(cGA_summ, 1, "j_fg")
cGA_graph_df[, 1] <- j_fg[5]
cGA_plot <- 
mk_SMAplot(df_points = c, df_lines = cGA_graph_df, gapeType = "ga", 
	point_colour = "dissected_by", labels = "dissected_by")


cGA <- sma(ga~SL*Region, data=c, log="xy", method="SMA", robust=T, slope.test=2,
	multcomp = TRUE, multcompmethod = "adjusted")
check_assump(cGA, "Co Gape Area by Regions")
cGA_regSumm <- mk_spp_summary(cGA, grouping=T)
cGA_reg_graphing <- mk_smaSPP_graph_df(cGA_regSumm, 5)
names(cGA_reg_graphing)[1] <- "Region"
cReg_by_dis <- 
mk_SMAplot(df_points = c, df_lines = cGA_reg_graphing, gapeType = "ga", 
	point_colour = "Region", line_colour = "Region", labels = "dissected_by")
mk_SMAfacets(df_points = c, df_lines = cGA_reg_graphing, gapeType = "ga",
	point_colour = "dissected_by", labels = "dissected_by", facetting = "Region")


cGA <- sma(ga~SL*dissected_by, data=c, log="xy", method="SMA", robust=T, 
	slope.test = 2, multcomp = TRUE, multcompmethod = "adjusted")
check_assump(cGA, "Co Gape Area dissected_by")
cGA_disSumm <- mk_spp_summary(cGA, grouping=T)
cGA_dis_graphing <- mk_smaSPP_graph_df(cGA_disSumm, length(unique(c$dissected_by)))
names(cGA_dis_graphing)[1] <- "dissected_by"
cDis_by_dis <- 
mk_SMAplot(df_points = c, df_lines = cGA_dis_graphing, gapeType = "ga", 
	point_colour = "dissected_by", line_colour = "dissected_by", 
	labels = "None")
mk_SMAfacets(df_points = c, df_lines = cGA_dis_graphing, gapeType = "ga",
	point_colour = "Region", labels = "None", facetting = "dissected_by")
	

pdf(file = "random_effects_plots/cGA.pdf", width=7.8, height = 6.2)
cReg_by_dis
cDis_by_reg
dev.off()

dev.copy2pdf(device=quartz, file = "random_effects_plots/cGA_reg.pdf", 
	width=7.8, height=6.2)

#-------------------------------------------------------------------------------
# Piscivores
pGA <- sma(ga~SL, data=p, log="xy", method="SMA", robust=T, slope.test=2)
check_assump(pGA, "Pi Gape Area All")
pGA_summ <- mk_sma_summary(pGA, 1)
pGA_graph_df <- mk_sma_graph_df(pGA_summ, 1, "j_fg")
pGA_graph_df[, 1] <- j_fg[1]
pGA_all <- 

mk_SMAplot(df_points = p, df_lines = pGA_graph_df, gapeType = "ga", 
	point_colour = "SpeciesCode", line_colour = "j_fg", labels = "None")


pGA <- sma(ga~SL*Region, data = p_w_sites, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = TRUE, multcompmethod = "adjusted")
check_assump(pGA, "Pi Gape Area by Regions")
pGA_regSumm <- mk_spp_summary(pGA, grouping=T)
pGA_reg_graphing <- mk_smaSPP_graph_df(pGA_regSumm, 6)
names(pGA_reg_graphing)[1] <- "Region"
pReg_by_dis <- 
mk_SMAplot( df_points = p, df_lines = pGA_reg_graphing, gapeType = "ga", 
	point_colour = "Region", line_colour = "Region", labels = "None")

mk_SMAfacets(df_points = p, df_lines = pGA_reg_graphing, gapeType = "ga", 
	point_colour = "SpeciesCode", labels = "dissected_by", facetting = "Region")

p_fewer_dis <- p[p$dissected_by %in% c("", "AB", "AO", "AO/KP", "LW", "RT",
	"RT/SC", "SC", "angeleen"), ]

pGA <- sma(ga~SL + dissected_by, data = p_fewer_dis, log = "xy", method = "SMA", 
	robust = T, elev.com = T, multcomp = TRUE, multcompmethod = "adjusted")

pGA <- sma(ga~SL*dissected_by, data = p_fewer_dis, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = TRUE, multcompmethod = "adjusted")
check_assump(pGA, "Pi Gape Area dissected_by")
pGA_disSumm <- mk_spp_summary(pGA, grouping=T)
pGA_dis_graphing <- mk_smaSPP_graph_df(pGA_disSumm, length(unique(pGA$groups)),
	"dissected_by")
pDis_by_reg <- 
mk_SMAplot(df_points = p, df_lines = pGA_dis_graphing, gapeType = "ga", 
	point_colour = "dissected_by", line_colour = "dissected_by", labels = "None")
pDis_by_spp <- 
mk_SMAfacets(df_points = p_fewer_dis, df_lines = pGA_dis_graphing, gapeType = "ga",
	point_colour = "SpeciesCode", labels = "None", facetting = "dissected_by")

dev.copy2pdf(device=quartz, file = "random_effects_plots/pGA_disby_spp.pdf", 
	width=7.8, height=6.2)

pdf(file = "random_effects_plots/pGA.pdf", width=7.8, height = 6.2)
pReg_by_dis
pDis_by_reg
pDis_by_spp
dev.off()

#-------------------------------------------------------------------------------
# Makes list of Pi species dataframes (subsets of fish)
p_spp_dfs <- split(p, p$SpeciesCode, drop=TRUE)

camela <- sma(ga~SL, data = p_spp_dfs$CA.MELA, log = "xy", method = "SMA", 
	robust = T, slope.test = 2)
check_assump(camela, "CA.MELA Gape Area")
camela_summ <- mk_sma_summary(camela, 1)
camela_graph_df <- mk_sma_graph_df(camela_summ, 1, "SpeciesCode")
camela_graph_df[, 1] <- pento_order[1]
camela_plot <- 
mk_SMAplot(df_points = p_spp_dfs$CA.MELA, df_lines = camela_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode",
	labels = "dissected_by")

camela <- sma(ga~SL*Region, data=p_spp_dfs$CA.MELA, log="xy", method="SMA", 
	robust = T, slope.test=2)
check_assump(camela, "CA.MELA Gape Area by Region")
camela_summ <- mk_spp_summary(camela, grouping=T)
camela_graph_df <- mk_smaSPP_graph_df(camela_summ, 2, "Region")
camela_byReg <- 
mk_SMAplot(df_points = p_spp_dfs$CA.MELA, df_lines = camela_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "dissected_by")
camela_byReg <- 
mk_SMAfacets(df_points = p_spp_dfs$CA.MELA, df_lines = camela_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "dissected_by", 
	facetting = "Region")

camela <- sma(ga~SL*dissected_by, data=p_spp_dfs$CA.MELA, log="xy", method="SMA", 
	robust = T, slope.test=2)
check_assump(camela, "CA.MELA Gape Area by Dis")
camela_summ <- mk_spp_summary(camela, 4, grouping=T)
camela_graph_df <- mk_smaSPP_graph_df(camela_summ, 4, "dissected_by")
camela_byDis <- 
mk_SMAplot(df_points = p_spp_dfs$CA.MELA, df_lines = camela_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")

mk_SMAfacets(df_points = p_spp_dfs$CA.MELA, df_lines = camela_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "Region", 
	facetting = "dissected_by")

pdf(file = "random_effects_plots/CaMela.pdf", width=7.8, height = 6.2)
camela_plot
camela_byReg
camela_byDis
dev.off()

dev.copy2pdf(device=quartz, file = "random_effects_plots/ca_mela_residuals.pdf", 
	width=6.7, height=10.8)

#-------------------------------------------------------------------------------
apfurc <- sma(ga~SL, data=p_spp_dfs$AP.FURC, log="xy", method="SMA", robust=T,
	slope.test=2)
check_assump(apfurc, "AP.FURC Gape Area")
apfurc_summ <- mk_sma_summary(apfurc, 1)
apfurc_graph_df <- mk_sma_graph_df(apfurc_summ, 1, "SpeciesCode")
apfurc_graph_df[, 1] <- pento_order[2]
apfurc_plot <- 
mk_SMAplot(df_points = p_spp_dfs$AP.FURC, df_lines = apfurc_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "None")

apfurc <- sma(ga~SL + Region, data=p_spp_dfs$AP.FURC, log="xy", method="SMA", 
	robust = T, elev.com=T, multcomp = T, multcompmethod = "adjusted")


apfurc <- sma(ga~SL*Region, data=p_spp_dfs$AP.FURC, log="xy", method="SMA", 
	robust = T, slope.test=2, multcomp = T, multcompmethod = "adjusted")
check_assump(apfurc, "AP.FURC Gape Area by Region")
apfurc_summ <- mk_spp_summary(apfurc, 4, grouping=T)
apfurc_graph_df <- mk_smaSPP_graph_df(apfurc_summ, 4, "Region")
mk_SMAplot(df_points = p_spp_dfs$AP.FURC, df_lines = apfurc_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "dissected_by")
apfurc_byReg <- 
mk_SMAfacets(df_points = p_spp_dfs$AP.FURC, df_lines = apfurc_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "dissected_by", 
	facetting = "Region")

apfurc <- sma(ga~SL + dissected_by, data=p_spp_dfs$AP.FURC, log="xy", 
	method="SMA", robust = T, elev.com=T, multcomp = T, 
	multcompmethod = "adjusted")


apfurc <- sma(ga~SL*dissected_by, data=p_spp_dfs$AP.FURC, log="xy", method="SMA", 
	robust = T, slope.test=2, multcomp = T, multcompmethod = "adjusted")
check_assump(apfurc, "AP.FURC Gape Area by Dis")
apfurc_summ <- mk_spp_summary(apfurc, 4, grouping=T)
apfurc_graph_df <- mk_smaSPP_graph_df(apfurc_summ, 4, "dissected_by")
mk_SMAplot(df_points = p_spp_dfs$AP.FURC, df_lines = apfurc_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "Region")
apfurc_byDis <-
mk_SMAfacets(df_points = p_spp_dfs$AP.FURC, df_lines = apfurc_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "Region", 
	facetting = "dissected_by")


dev.copy2pdf(device=quartz, file = "random_effects_plots/ap_furc_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/ApFurc.pdf", width=7.8, height = 6.2)
apfurc_plot
apfurc_byReg
apfurc_byDis
dev.off()

#-------------------------------------------------------------------------------
# To remove outlier and see how results are affected
row_kif12_130 <- which(p_spp_dfs$LU.BOHA$SpecimenID == 'KIF12_130')

luboha <- sma(ga~SL, data =p_spp_dfs$LU.BOHA)
luboha_summ <- mk_sma_summary(luboha, 1)
luboha_graph_df <- mk_sma_graph_df(luboha_summ, 1, "SpeciesCode")
luboha_graph_df[, 1] <- pento_order[3]
mk_SMAplot(df_points = p_spp_dfs$LU.BOHA, df_lines = luboha_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "SpecimenID")


luboha <- sma(ga~SL, data=p_spp_dfs$LU.BOHA[-row_kif12_130, ], log="xy", 
	method="SMA", robust=T, slope.test=2)
check_assump(luboha, "LU.BOHA Gape Area")
luboha_summ <- mk_sma_summary(luboha, 1)
luboha_graph_df <- mk_sma_graph_df(luboha_summ, 1, "SpeciesCode")
luboha_graph_df[, 1] <- pento_order[3]
luboha_plot <- 
mk_SMAplot(df_points = p_spp_dfs$LU.BOHA, df_lines = luboha_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "None")

luboha <- sma(ga~SL + Region, data=p_spp_dfs$LU.BOHA[-row_kif12_130, ], 
	log="xy", method="SMA", robust = T, elev.com=T, multcomp = T, 
	multcompmethod = "adjusted")


luboha <- sma(ga~SL*Region, data=p_spp_dfs$LU.BOHA[-row_kif12_130, ], log="xy", 
	method="SMA", robust = T, slope.test=2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(luboha, "LU.BOHA Gape Area by Region")
luboha_summ <- mk_spp_summary(luboha, 5, grouping=T)
luboha_graph_df <- mk_smaSPP_graph_df(luboha_summ, 5, "Region")
luboha_byReg <- 
mk_SMAplot(df_points = p_spp_dfs$LU.BOHA, df_lines = luboha_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "None")
luboha_byReg <- 
mk_SMAfacets(df_points = p_spp_dfs$LU.BOHA, df_lines = luboha_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "dissected_by", 
	facetting = "Region")

luboha <- sma(ga~SL + dissected_by, data=p_spp_dfs$LU.BOHA[-row_kif12_130, ], 
	log="xy", method="SMA", robust = T, elev.com=T, multcomp = T, 
	multcompmethod = "adjusted")


luboha <- sma(ga~SL*dissected_by, data=p_spp_dfs$LU.BOHA[-row_kif12_130, ], 
	log="xy", method="SMA", robust = T, slope.test=2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(luboha, "LU.BOHA Gape Area by Dis")
luboha_summ <- mk_spp_summary(luboha, 5, grouping=T)
luboha_graph_df <- mk_smaSPP_graph_df(luboha_summ, 5, "dissected_by")
luboha_byDis <-
mk_SMAplot(df_points = p_spp_dfs$LU.BOHA, df_lines = luboha_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")
luboha_byDis <-
mk_SMAfacets(df_points = p_spp_dfs$LU.BOHA, df_lines = luboha_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")


dev.copy2pdf(device=quartz, file = "random_effects_plots/lu_boha_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/LuBoha.pdf", width=7.8, height = 6.2)
luboha_plot
luboha_byReg
luboha_byDis
dev.off()

#-------------------------------------------------------------------------------
# lukasm was excluded from the species specific analyses
#-------------------------------------------------------------------------------
row_kif12_171 <- which(p_spp_dfs$CE.ARGU$SpecimenID == 'KIF12_171')

ceargu <- sma(ga~SL, data=p_spp_dfs$CE.ARGU, log="xy", method="SMA", robust=T,
	slope.test=2)
check_assump(ceargu, "CE.ARGU Gape Area")
ceargu_summ <- mk_sma_summary(ceargu, 1)
ceargu_graph_df <- mk_sma_graph_df(ceargu_summ, 1, "SpeciesCode")
ceargu_graph_df[, 1] <- pento_order[5]
ceargu_plot <- 
mk_SMAplot(df_points = p_spp_dfs$CE.ARGU, df_lines = ceargu_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "None")

ceargu <- sma(ga~SL + Region, data=p_spp_dfs$CE.ARGU[-row_kif12_171, ], log="xy", method="SMA", 
	robust = T, elev.com=T, multcomp = T, multcompmethod = "adjusted")


ceargu <- sma(ga~SL*Region, data=p_spp_dfs$CE.ARGU[-row_kif12_171, ], 
	log="xy", method="SMA", robust = T, slope.test=2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(ceargu, "CE.ARGU Gape Area by Region")
ceargu_summ <- mk_spp_summary(ceargu, 4, grouping=T)
ceargu_graph_df <- mk_smaSPP_graph_df(ceargu_summ, 4, "Region")
ceargu_byReg <- 
mk_SMAplot(df_points = p_spp_dfs$CE.ARGU, df_lines = ceargu_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "dissected_by")
ceargu_byReg <- 
mk_SMAfacets(df_points = p_spp_dfs$CE.ARGU, df_lines = ceargu_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "Region", 
	facetting = "Region")

ceargu <- sma(ga~SL + dissected_by, data=p_spp_dfs$CE.ARGU, log="xy", 
	method="SMA", robust = T, elev.com=T, multcomp = T, 
	multcompmethod = "adjusted")

ceargu <- sma(ga~SL*dissected_by, data=p_spp_dfs$CE.ARGU, log="xy", 
	method="SMA", robust = T, slope.test = 2)
check_assump(ceargu, "CE.ARGU Gape Area by Dis")
ceargu_summ <- mk_spp_summary(ceargu, 8, grouping=T)
ceargu_graph_df <- mk_smaSPP_graph_df(ceargu_summ, 8, "dissected_by")
ceargu_byDis <-
mk_SMAplot(df_points = p_spp_dfs$CE.ARGU, df_lines = ceargu_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")
ceargu_byDis <-
mk_SMAfacets(df_points = p_spp_dfs$CE.ARGU, df_lines = ceargu_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")

ddply(p_spp_dfs$CE.ARGU, .(dissected_by), summarise, length(SpecimenID))
ddply(argus_fewer_dis, .(dissected_by), summarise, length(SpecimenID))

dev.copy2pdf(device=quartz, file = "random_effects_plots/ce_argu_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/CeArgu.pdf", width=7.8, height = 6.2)
ceargu_plot
ceargu_byReg
ceargu_byDis
dev.off()

#-------------------------------------------------------------------------------
ceurod <- sma(ao_corrected~SL, data=p_spp_dfs$CE.UROD, log="xy", method="SMA", 
	robust=T, slope.test=2)
plot(ceurod, which = "residual")
abline(h = 0, col = "red")
plot(ceurod)

ceurod <- sma(ga~SL * dissected_by, data=p_spp_dfs$CE.UROD, log="xy", 
	method="SMA", robust=T, slope.test=2)
plot(ceurod, which = "residual")
abline(h = 0, col = "red")
plot(ceurod)
ceurod_summ <- mk_spp_summary(ceurod, 6, grouping=T)
ceurod_graph_df <- mk_smaSPP_graph_df(ceurod_summ, 6, "dissected_by")
mk_SMAplot(df_points = p_spp_dfs$CE.UROD, df_lines = ceurod_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by", 
	labels = "SpecimenID")


check_assump(ceurod, "CE.UROD Gape Area")
ceurod_summ <- mk_sma_summary(ceurod, 1)
ceurod_graph_df <- mk_sma_graph_df(ceurod_summ, 1, "SpeciesCode")
ceurod_graph_df[, 1] <- pento_order[6]
ceurod_plot <- 
mk_SMAplot(df_points = p_spp_dfs$CE.UROD, df_lines = ceurod_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "dissected_by")

ceurod <- sma(ga~SL + Region, data=p_spp_dfs$CE.UROD, log="xy", method="SMA", 
	robust = T, elev.com=T, multcomp = T, multcompmethod = "adjusted")


ceurod <- sma(ga~SL*Region, data=p_spp_dfs$CE.UROD, log="xy", method="SMA", 
	robust = T, slope.test=2, multcomp = T, multcompmethod = "adjusted")
check_assump(ceurod, "CE.UROD Gape Area by Region")
ceurod_summ <- mk_spp_summary(ceurod, 5, grouping=T)
ceurod_graph_df <- mk_smaSPP_graph_df(ceurod_summ, 5, "Region")
ceurod_byReg <- 
mk_SMAplot(df_points = p_spp_dfs$CE.UROD, df_lines = ceurod_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "dissected_by")
ceurod_byReg <- 
mk_SMAfacets(df_points = p_spp_dfs$CE.UROD, df_lines = ceurod_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "dissected_by", 
	facetting = "Region")

ceurod <- sma(ga~SL + dissected_by, data=p_spp_dfs$CE.UROD, log="xy", 
	method="SMA", robust = T, elev.com=T, multcomp = T, 
	multcompmethod = "adjusted")

ceurod <- sma(ga~SL*observer_id, data=p_spp_dfs$CE.UROD, log="xy", 
	method="SMA", robust = T, slope.test = 2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(ceurod, "CE.UROD Gape Area by Dis")
ceurod_summ <- mk_spp_summary(ceurod, 6, grouping=T)
ceurod_graph_df <- mk_smaSPP_graph_df(ceurod_summ, 6, "observer_id")
ceurod_byDis <-
mk_SMAplot(df_points = p_spp_dfs$CE.UROD, df_lines = ceurod_graph_df,
	gapeType = "ga", point_colour = "observer_id", line_colour = "observer_id",
	labels = "Region") +
	
ceurod_byDis <-
mk_SMAfacets(df_points = p_spp_dfs$CE.UROD, df_lines = ceurod_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "Region", 
	facetting = "dissected_by")

ddply(p_spp_dfs$CE.UROD, .(dissected_by), summarise, length(SpecimenID))
ddply(argus_fewer_dis, .(dissected_by), summarise, length(SpecimenID))

dev.copy2pdf(device=quartz, file = "random_effects_plots/ce_urod_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/CeUrod.pdf", width=7.8, height = 6.2)
ceurod_plot
ceurod_byReg
ceurod_byDis
dev.off()

# Taking CEUORD farther -- correcting for large and significantly lower 
# elevations for measurements made by AO/Angeleen

observer_grouping <- c(rep(1, 4), 2, rep(1, 5))
dissected_by <- unique(p_spp_dfs$CE.UROD$dissected_by)
observer_grouping <- data.frame(dissected_by, observer_grouping)
ceurod_ao_sep <- merge(p_spp_dfs$CE.UROD, observer_grouping, by = "dissected_by")

ceurod_2dis <- sma(ga~SL*observer_grouping, data = ceurod_ao_sep, log = "xy",
	method = "SMA", robust = T, slope.test = 2)

ceurod_2dis <- sma(ga~SL + observer_grouping, data = ceurod_ao_sep, log = "xy",
	method = "SMA", robust = T, slope.test = 2)

#-------------------------------------------------------------------------------
valout <- sma(ga~SL, data=p_spp_dfs$VA.LOUT, log="xy", method="SMA", robust=T,
	slope.test=2)
check_assump(valout, "VA.LOUT Gape Area")
valout_summ <- mk_sma_summary(valout, 1)
valout_graph_df <- mk_sma_graph_df(valout_summ, 1, "SpeciesCode")
valout_graph_df[, 1] <- pento_order[7]
valout_plot1 <- 
mk_SMAplot(df_points = p_spp_dfs$VA.LOUT, df_lines = valout_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "Region")

valout_plot2 <- 
mk_SMAplot(df_points = p_spp_dfs$VA.LOUT, df_lines = valout_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "dissected_by")


valout <- sma(ga~SL + Region, data=p_spp_dfs$VA.LOUT, log="xy", method="SMA", 
	robust = T, elev.com=T, multcomp = T, multcompmethod = "adjusted")


valout <- sma(ga~SL*Region, data=p_spp_dfs$VA.LOUT, log="xy", method="SMA", 
	robust = T, slope.test=2, multcomp = T, multcompmethod = "adjusted")
check_assump(valout, "VA.LOUT Gape Area by Region")
valout_summ <- mk_spp_summary(valout, 4, grouping=T)
valout_graph_df <- mk_smaSPP_graph_df(valout_summ, 4, "Region")
valout_byReg <- 
mk_SMAplot(df_points = p_spp_dfs$VA.LOUT, df_lines = valout_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "dissected_by")
valout_byReg <- 
mk_SMAfacets(df_points = p_spp_dfs$VA.LOUT, df_lines = valout_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "dissected_by", 
	facetting = "Region")

valout <- sma(ga~SL + dissected_by, data=p_spp_dfs$VA.LOUT, log="xy", 
	method="SMA", robust = T, elev.com=T)

valout <- sma(ga~SL*dissected_by, data=p_spp_dfs$VA.LOUT, log="xy", 
	method="SMA", robust = T, slope.test = 2)
check_assump(valout, "VA.LOUT Gape Area by Dis")
valout_summ <- mk_spp_summary(valout, 6, grouping=T)
valout_graph_df <- mk_smaSPP_graph_df(valout_summ, 6, "dissected_by")
valout_byDis <-
mk_SMAplot(df_points = p_spp_dfs$VA.LOUT, df_lines = valout_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "Region")
valout_byDis <-
mk_SMAfacets(df_points = p_spp_dfs$VA.LOUT, df_lines = valout_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "Region", 
	facetting = "dissected_by")

dev.copy2pdf(device=quartz, file = "random_effects_plots/va_lout_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/VaLout.pdf", width=7.8, height = 6.2)
valout_plot1
valout_plot2
dev.off()


#===============================================================================
# Family and species level varation in the Piscivores
#===============================================================================
# All
pGA <- sma(ga~SL, data=p, log="xy", method="SMA", robust=T, slope.test=2)
check_assump(pGA, "Pi Gape Area All")
pGA_summ <- mk_sma_summary(pGA, 1)
pGA_graph_df <- mk_sma_graph_df(pGA_summ, 1, "j_fg")
pGA_graph_df[, 1] <- j_fg[1]
pGA_all <- 
mk_SMAplot(df_points = p, df_lines = pGA_graph_df, gapeType = "ga", 
	point_colour = "SpeciesCode", line_colour = "j_fg", labels = "None")


pfam <- sma(ga~SL + Family, data = p, log = "xy", method = "SMA", robust = T, 
	slope.test = 2, multcomp = T, multcompmethod = "adjusted")

pfam <- sma(ga~SL * Family, data = p, log = "xy", method = "SMA", robust = T, 
	slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(pfam, "Pi Gape Area by Family")
pfam_summ <- mk_spp_summary(pfam, 3, grouping = T)
pfam_graph_df <- mk_smaSPP_graph_df(pfam_summ, 3, "Family")

pfam_plot <- 
mk_SMAplot(df_points = p, df_lines = pfam_graph_df, gapeType = "ga",
 point_colour = "Family", line_colour = "Family", labels = "None")

pfam_facets <- 
mk_SMAfacets(df_points = p, df_lines = pfam_graph_df, gapeType = "ga",
	point_colour = "SpeciesCode", labels = "None", facetting = "Family")

ddply(p, .(Family), summarise, length(SpecimenID))
ddply(p, .(SpeciesCode), summarise, length(SpecimenID))

pspp <- sma(ga~SL * SpeciesCode, data = p, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")

pspp <- sma(ga~SL * SpeciesCode, data = p, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(pspp, "Pi Gape Area by Species")
pspp_summ <- mk_spp_summary(pspp, 7, grouping = T)
pspp_graph_df <- mk_smaSPP_graph_df(pspp_summ, 7, "SpeciesCode")
pspp_plot <- 
mk_SMAplot(df_points = p, df_lines = pspp_graph_df, gapeType = "ga",
 point_colour = "SpeciesCode", line_colour = "SpeciesCode", labels = "None")

# Removing CE.UROD for testing differences in slopes and elevation between
# remaining species:

ceurod_rows <- which(p$SpeciesCode == 'CE.UROD')
no_ceurod <- p[-ceurod_rows, ]
pspp <- sma(ga~SL * SpeciesCode, data = no_ceurod, log = "xy", method = "SMA",
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")

pspp <- sma(ga~SL + SpeciesCode, data = no_ceurod, log = "xy", method = "SMA",
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")

pfam <- sma(ga~SL * Family, data = no_ceurod, log = "xy", method = "SMA", robust = T, 
	slope.test = 2, multcomp = T, multcompmethod = "adjusted")

# Removing AO observations 
ao_rows <- which(p$SpeciesCode == 'CE.UROD' & p$dissected_by == 'AO')
no_ao <- p[-ao_rows, ]
no_ao_slope <- with(no_ao, slope.com(log(ga), log(SL), SpeciesCode, method = 'SMA', robust = T))

pspp <- sma(ga~SL + SpeciesCode, data = no_ao, log = "xy", method = "SMA",
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")


pspp <- sma(ga~SL * SpeciesCode, data = no_ao, log = "xy", method = "SMA",
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")

# Obtaining overall group slope estimate:
with(p, slope.com(log(ga), log(SL), SpeciesCode, method = 'SMA', robust = T))
with(no_ceurod, slope.com(log(ga), log(SL), Family, method = 'SMA', robust = T))

with(no_ceurod, slope.com(log(ga), log(SL), SpeciesCode, method = 'SMA'))


dev.copy2pdf(device=quartz, file = "random_effects_plots/pfam_pspp_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/pfam_pspp.pdf", width=7.8, height = 6.2)
pfam_plot
pfam_facets
pspp_plot
dev.off()

#-------------------------------------------------------------------------------
# Benthic Invertivores

bGA <- sma(ga~SL, data = b= b, log="xy", method="SMA", robust=T, slope.test=2)
check_assump(bGA, "BI Gape Area All")
bGA_summ <- mk_sma_summary(bGA, "j_fg")
bGA_graph_df <- mk_sma_graph_df(bGA_summ, 1, "j_fg")
bGA_graph_df[, 1] <- j_fg[2]
bGA_plot <- 
mk_SMAplot(df_points = b, df_lines = bGA_graph_df, gapeType = "ga", 
	point_colour = "j_fg", line_colour = "j_fg", labels = "None")

with(b, slope.com(log(ga), log(SL), Region, method = 'SMA', robust = T, 
	slope.test = 2))

bGA <- sma(ga~SL * Region, data = b, log = "xy", method = "SMA", robust = T,
	slope.test = 2, multcomp = TRUE, multcompmethod = "adjusted")
check_assump(bGA, "BI Gape Area by Regions")
bGA_regSumm <- mk_spp_summary(bGA, grouping=T)
bGA_reg_graphing <- mk_smaSPP_graph_df(bGA_regSumm, 5, "Region")
bReg_by_dis <- 
mk_SMAplot(df_points = b, df_lines = bGA_reg_graphing, gapeType = "ga", 
	point_colour = "Region", line_colour = "Region", labels = "None")
mk_SMAfacets(df_points = b, df_lines = bGA_reg_graphing, gapeType = "ga",
	point_colour = "SpeciesCode", labels = "dissected_by", facetting = "Region")

aokp_gaj_rows <- which(b$dissected_by %in% c('AO/KP', 'GAJ', 'AO'))
no_aokp_gaj <- b[-aokp_gaj_rows, ]

bGA <- sma(ga~SL + dissected_by, data = no_aokp_gaj, log = "xy", method = "SMA", robust = T, 
	slope.test = 2, multcomp = TRUE, multcompmethod = "adjusted")

with(no_aokp_gaj, slope.com(log(ga), log(SL), dissected_by, method = 'SMA', robust = T, 
	slope.test = 2))

bGA <- sma(ga~SL * dissected_by, data = no_aokp_gaj, log = "xy", method = "SMA", robust = T, 
	slope.test = 2, multcomp = TRUE, multcompmethod = "adjusted")
check_assump(bGA, "BI Gape Area dissected_by")
bGA_disSumm <- mk_spp_summary(bGA, grouping=T)
bGA_dis_graphing <- mk_smaSPP_graph_df(bGA_disSumm, length(unique(b$dissected_by)),
	"dissected_by")
bDis_by_spp <- 
mk_SMAplot(df_points = b, df_lines = bGA_dis_graphing, gapeType = "ga", 
	point_colour = "dissected_by", line_colour = "dissected_by", 
	labels = "None")
mk_SMAfacets(df_points = b, df_lines = bGA_dis_graphing, gapeType = "ga",
	point_colour = "SpeciesCode", labels = "None", facetting = "dissected_by")
	
ddply(no_aokp_gaj, .(dissected_by), summarise, (length(SpecimenID)))

pdf(file = "random_effects_plots/bGA.pdf", width=7.8, height = 6.2)
bGA_plot
bReg_by_dis
bDis_by_spp
dev.off()

dev.copy2pdf(device=quartz, file = "random_effects_plots/bGA_residuals.pdf", 
	width=6.7, height=10.8)

dev.copy2pdf(device=quartz, file = "random_effects_plots/bGA_AO_remov_residuals.pdf", 
	width=6.7, height=4)

dev.copy2pdf(device=quartz, file = "random_effects_plots/bGA_reg.pdf", 
	width=7.8, height=6.2)

#-------------------------------------------------------------------------------
# Makes list of BI species dataframes (subsets of fish)
b_spp_dfs <- split(b, b$SpeciesCode, drop=TRUE)

# Monotaxis grandoculis
mogran <- sma(ga~SL, data=b_spp_dfs$MO.GRAN, log="xy", method="SMA", robust=T,
	slope.test=2)
check_assump(mogran, "MO.GRAN Gape Area")
mogran_summ <- mk_sma_summary(mogran, "SpeciesCode")
mogran_graph_df <- mk_sma_graph_df(mogran_summ, 1, "SpeciesCode")
mogran_graph_df[, 1] <- pento_order[9]
mogran_plot <- 
mk_SMAplot(df_points = b_spp_dfs$MO.GRAN, df_lines = mogran_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "None")

mogran <- sma(ga~SL + Region, data=b_spp_dfs$MO.GRAN, log="xy", method="SMA", 
	robust = T, elev.com=T, multcomp = T, multcompmethod = "adjusted")

extra_row <- which(b_spp_dfs$MO.GRAN$Region == 'EXTRA')
kif11_107 <- which(b_spp_dfs$MO.GRAN$SpecimenID == 'KIF11_107')



with(b_spp_dfs$MO.GRAN[-c(kif11_107, extra_row), ], slope.com(log(ga), log(SL), Region, method = 'SMA', 
	robust = T, slope.test = 2))

test <- b_spp_dfs$MO.GRAN[-kif11_107, ]
lplf <- which(b_spp_dfs$MO.GRAN$Region == 'LP.LF')
no_lplf <- b_spp_dfs$MO.GRAN[-lplf, ]

mogran <- sma(ga~SL * Region, data = test, log = "xy", 
	method = "SMA", robust = T, slope.test = 2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(mogran, "MO.GRAN Gape Area by Region")
mogran_summ <- mk_spp_summary(mogran, 5, grouping=T)
mogran_graph_df <- mk_smaSPP_graph_df(mogran_summ, 5, "Region")
mogran_byReg <- 
mk_SMAplot(df_points = b_spp_dfs$MO.GRAN, df_lines = mogran_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "None")
mogran_byReg <- 
mk_SMAfacets(df_points = b_spp_dfs$MO.GRAN, df_lines = mogran_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "None", 
	facetting = "Region")

mogran <- sma(ga~SL + dissected_by, data = b_spp_dfs$MO.GRAN, log = "xy", 
	method = "SMA", robust = T, elev.com = T)

mogran <- sma(ga~SL * dissected_by, data = b_spp_dfs$MO.GRAN, log = "xy", 
	method = "SMA", robust = T, slope.test = 2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(mogran, "MO.GRAN Gape Area by Dis")
mogran_summ <- mk_spp_summary(mogran, 6, grouping = T)
mogran_graph_df <- mk_smaSPP_graph_df(mogran_summ, 6, "dissected_by")
mogran_byDis <-
mk_SMAplot(df_points = b_spp_dfs$MO.GRAN, df_lines = mogran_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")
mogran_byDis <-
mk_SMAfacets(df_points = b_spp_dfs$MO.GRAN, df_lines = mogran_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")

dev.copy2pdf(device=quartz, file = "random_effects_plots/mo_gran_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/MoGran.pdf", width=7.8, height = 6.2)
mogran_plot
mogran_byReg
mogran_byDis
dev.off()

#-------------------------------------------------------------------------------
paarca <- sma(ga~SL, data = b_spp_dfs$PA.ARCA, log = "xy", method = "SMA", 
	robust = T, slope.test = 2)
check_assump(paarca, "PA.ARCA Gape Area")
paarca_summ <- mk_sma_summary(paarca, "SpeciesCode")
paarca_graph_df <- mk_sma_graph_df(paarca_summ, 1, "SpeciesCode")
paarca_graph_df[, 1] <- pento_order[8]
paarca_plot <- 
mk_SMAplot(df_points = b_spp_dfs$PA.ARCA, df_lines = paarca_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "None")

paarca <- sma(ga~SL + Region, data=b_spp_dfs$PA.ARCA, log="xy", method="SMA", 
	robust = T, elev.com=T, multcomp = T, multcompmethod = "adjusted")


paarca <- sma(ga~SL * Region, data = b_spp_dfs$PA.ARCA, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(paarca, "PA.ARCA Gape Area by Region")
paarca_summ <- mk_spp_summary(paarca, 4, grouping=T)
paarca_graph_df <- mk_smaSPP_graph_df(paarca_summ, 3, "Region")
paarca_byReg <- 
mk_SMAplot(df_points = b_spp_dfs$PA.ARCA, df_lines = paarca_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "dissected_by")
paarca_byReg <- 
mk_SMAfacets(df_points = b_spp_dfs$PA.ARCA, df_lines = paarca_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "dissected_by", 
	facetting = "Region")

paarca <- sma(ga~SL + dissected_by, data = b_spp_dfs$PA.ARCA, log = "xy", 
	method = "SMA", robust = T, elev.com = T)

paarca <- sma(ga~SL * dissected_by, data = b_spp_dfs$PA.ARCA, log = "xy", 
	method = "SMA", robust = T, slope.test = 2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(paarca, "PA.ARCA Gape Area by Dis")
paarca_summ <- mk_spp_summary(paarca, 9, grouping=T)
paarca_graph_df <- mk_smaSPP_graph_df(paarca_summ, 9, "dissected_by")
paarca_byDis <-
mk_SMAplot(df_points = b_spp_dfs$PA.ARCA, df_lines = paarca_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "Region")
paarca_byDis <-
mk_SMAfacets(df_points = b_spp_dfs$PA.ARCA, df_lines = paarca_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")

dev.copy2pdf(device=quartz, file = "random_effects_plots/pa_arca_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/PaArca.pdf", width=7.8, height = 6.2)
paarca_plot
paarca_byReg
paarca_byDis
dev.off()

#-------------------------------------------------------------------------------
painsu <- sma(ga~SL, data = b_spp_dfs$PA.INSU, log = "xy", method = "SMA", 
	robust = T, slope.test = 2)
check_assump(painsu, "PA.INSU Gape Area")
painsu_summ <- mk_sma_summary(painsu, "SpeciesCode")
painsu_graph_df <- mk_sma_graph_df(painsu_summ, 1, "SpeciesCode")
painsu_graph_df[, 1] <- pento_order[10]
painsu_plot <- 
mk_SMAplot(df_points = b_spp_dfs$PA.INSU, df_lines = painsu_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "None")

painsu <- sma(ga~SL + Region, data=b_spp_dfs$PA.INSU, log="xy", method="SMA", 
	robust = T, elev.com=T, multcomp = T, multcompmethod = "adjusted")


painsu <- sma(ga~SL * Region, data = b_spp_dfs$PA.INSU, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(painsu, "PA.INSU Gape Area by Region")
painsu_summ <- mk_spp_summary(painsu, 5, grouping=T)
painsu_graph_df <- mk_smaSPP_graph_df(painsu_summ, 5, "Region")
painsu_byReg <- 
mk_SMAplot(df_points = b_spp_dfs$PA.INSU, df_lines = painsu_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "None")
painsu_byReg <- 
mk_SMAfacets(df_points = b_spp_dfs$PA.INSU, df_lines = painsu_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "None", 
	facetting = "Region")

painsu <- sma(ga~SL + dissected_by, data = b_spp_dfs$PA.INSU, log = "xy", 
	method = "SMA", robust = T, elev.com = T, multcomp = T, 
	multcompmethod = "adjusted")

painsu <- sma(ga~SL * dissected_by, data = b_spp_dfs$PA.INSU, log = "xy", 
	method = "SMA", robust = T, slope.test = 2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(painsu, "PA.INSU Gape Area by Dis")
painsu_summ <- mk_spp_summary(painsu, 14, grouping=T)
painsu_graph_df <- mk_smaSPP_graph_df(painsu_summ, 14, "dissected_by")
painsu_byDis <-
mk_SMAplot(df_points = b_spp_dfs$PA.INSU, df_lines = painsu_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")
painsu_byDis <-
mk_SMAfacets(df_points = b_spp_dfs$PA.INSU, df_lines = painsu_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")

dev.copy2pdf(device=quartz, file = "random_effects_plots/pa_insu_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/PaInsu.pdf", width=7.8, height = 6.2)
painsu_plot
painsu_byReg
painsu_byDis
dev.off()


#===============================================================================
# Species level varation in the Benthic Invertivores
#===============================================================================
# All
bGA <- sma(ga~SL, data = b, log = "xy", method = "SMA", robust = T, 
	slope.test = 2)
check_assump(bGA, "BI Gape Area All")
bGA_summ <- mk_sma_summary(bGA, 1)
bGA_graph_df <- mk_sma_graph_df(bGA_summ, 1, "j_fg")
bGA_graph_df[, 1] <- j_fg[1]
bGA_all <- 
mk_SMAplot(df_points = b, df_lines = bGA_graph_df, gapeType = "ga", 
	point_colour = "SpeciesCode", line_colour = "j_fg", labels = "None")

ddply(b, .(SpeciesCode), summarise, length(SpecimenID))

# Species
bspp <- sma(ga~SL + SpeciesCode, data = b, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")


bspp <- sma(ga~SL * SpeciesCode, data = b, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(bspp, "BI Gape Area by Species")
bspp_summ <- mk_spp_summary(bspp, 3, grouping = T)
bspp_graph_df <- mk_smaSPP_graph_df(bspp_summ, 3, "SpeciesCode")
bspp_plot <- 
mk_SMAplot(df_points = b, df_lines = bspp_graph_df, gapeType = "ga",
 point_colour = "SpeciesCode", line_colour = "SpeciesCode", labels = "None")


# Obtaining overall group slope estimate:
with(b, slope.com(log(ga), log(SL), SpeciesCode, method = 'SMA', robust = T))

dev.copy2pdf(device=quartz, file = "random_effects_plots/bspp_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/bspp.pdf", width=7.8, height = 6.2)
bspp_plot
dev.off()


#-------------------------------------------------------------------------------
# Zooplanktivores

zGA <- sma(ga~SL, data = z, log = "xy", method = "SMA", robust = T, slope.test = 2)
check_assump(zGA, "ZP Gape Area All")
zGA_summ <- mk_sma_summary(zGA, "j_fg")
zGA_graph_df <- mk_sma_graph_df(zGA_summ, 1, "j_fg")
zGA_graph_df[, 1] <- j_fg[3]
zGA_plot <- 
mk_SMAplot(df_points = z, df_lines = zGA_graph_df, gapeType = "ga", 
	point_colour = "j_fg", line_colour = "j_fg", labels = "None")

with(z, slope.com(log(ga), log(SL), Region, method = 'SMA', robust = T, 
	slope.test = 2))

zGA <- sma(ga~SL * Region, data = z, log = "xy", method = "SMA", robust = T,
	slope.test = 2, multcomp = TRUE, multcompmethod = "adjusted")
check_assump(zGA, "ZP Gape Area by Regions")
zGA_regSumm <- mk_spp_summary(zGA, 6, grouping = T)
zGA_reg_graphing <- mk_smaSPP_graph_df(zGA_regSumm, 5, "Region")
zReg_by_spp <- 
mk_SMAplot(df_points = z, df_lines = zGA_reg_graphing, gapeType = "ga", 
	point_colour = "Region", line_colour = "Region", labels = "None")
mk_SMAfacets(df_points = z, df_lines = zGA_reg_graphing, gapeType = "ga",
	point_colour = "SpeciesCode", labels = "None", facetting = "Region")


zGA <- sma(ga~SL + dissected_by, data = no_aokp_gaj, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = TRUE, multcompmethod = "adjusted")

ao_adrian <- which(z$dissected_by == 'AO/Adrian')
no_ao_adrian <- z[-ao_adrian, ]

zGA <- sma(ga~SL * dissected_by, data = no_ao_adrian, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = TRUE, multcompmethod = "adjusted")
check_assump(zGA, "ZP Gape Area dissected_by")
zGA_disSumm <- mk_spp_summary(zGA, 20, grouping = T)
zGA_dis_graphing <- mk_smaSPP_graph_df(zGA_disSumm, 16, "dissected_by")
zDis_by_spp <- 
mk_SMAplot(df_points = z, df_lines = zGA_dis_graphing, gapeType = "ga", 
	point_colour = "dissected_by", line_colour = "dissected_by", 
	labels = "None")
mk_SMAfacets(df_points = z, df_lines = zGA_dis_graphing, gapeType = "ga",
	point_colour = "SpeciesCode", labels = "None", facetting = "dissected_by")
	

pdf(file = "random_effects_plots/zGA.pdf", width=7.8, height = 6.2)
zGA_plot
zReg_by_spp
zDis_by_spp
dev.off()

dev.copy2pdf(device=quartz, file = "random_effects_plots/zGA_residuals.pdf", 
	width=6.7, height=10.8)

#-------------------------------------------------------------------------------
# Makes list of ZP species dataframes (subsets of fish)
z_spp_dfs <- split(z, z$SpeciesCode, drop=TRUE)

pttile <- sma(ga~SL, data = z_spp_dfs$PT.TILE, log = "xy", method = "SMA", 
	robust = T, slope.test = 2)
check_assump(pttile, "PT.TILE Gape Area")
pttile_summ <- mk_sma_summary(pttile, "SpeciesCode")
pttile_graph_df <- mk_sma_graph_df(pttile_summ, 1, "SpeciesCode")
pttile_graph_df[, 1] <- pento_order[18]
pttile_plot <- 
mk_SMAplot(df_points = z_spp_dfs$PT.TILE, df_lines = pttile_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "None")

pttile <- sma(ga~SL + Region, data=b_spp_dfs$PT.TILE, log="xy", method="SMA", 
	robust = T, elev.com=T, multcomp = T, multcompmethod = "adjusted")


pttile <- sma(ga~SL * Region, data = z_spp_dfs$PT.TILE, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(pttile, "PT.TILE Gape Area by Region")
pttile_summ <- mk_spp_summary(pttile, 5, grouping=T)
pttile_graph_df <- mk_smaSPP_graph_df(pttile_summ, 5, "Region")
pttile_byReg <- 
mk_SMAplot(df_points = z_spp_dfs$PT.TILE, df_lines = pttile_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "None")
pttile_byReg <- 
mk_SMAfacets(df_points = z_spp_dfs$PT.TILE, df_lines = pttile_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "dissected_by", 
	facetting = "Region")

pttile <- sma(ga~SL + dissected_by, data = z_spp_dfs$PT.TILE, log = "xy", 
	method = "SMA", robust = T, elev.com = T, multcomp = T, 
	multcompmethod = "adjusted")

pttile <- sma(ga~SL * dissected_by, data = z_spp_dfs$PT.TILE, log = "xy", 
	method = "SMA", robust = T, slope.test = 2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(pttile, "PT.TILE Gape Area by Dis")
pttile_summ <- mk_spp_summary(pttile, 14, grouping=T)
pttile_graph_df <- mk_smaSPP_graph_df(pttile_summ, 5, "dissected_by")
pttile_byDis <-
mk_SMAplot(df_points = z_spp_dfs$PT.TILE, df_lines = pttile_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")
pttile_byDis <-
mk_SMAfacets(df_points = z_spp_dfs$PT.TILE, df_lines = pttile_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")

dev.copy2pdf(device=quartz, file = "random_effects_plots/pt_tile_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/PtTile.pdf", width=7.8, height = 6.2)
pttile_plot
pttile_byReg
pttile_byDis
dev.off()

#-------------------------------------------------------------------------------
# Makes list of ZP species dataframes (subsets of fish)
z_spp_dfs <- split(z, z$SpeciesCode, drop=TRUE)

catere <- sma(ga~SL, data = z_spp_dfs$CA.TERE, log = "xy", method = "SMA", 
	robust = T, slope.test = 2)
check_assump(catere, "CA.TERE Gape Area")
catere_summ <- mk_sma_summary(catere, "SpeciesCode")
catere_graph_df <- mk_sma_graph_df(catere_summ, 1, "SpeciesCode")
catere_graph_df[, 1] <- pento_order[18]
catere_plot <- 
mk_SMAplot(df_points = z_spp_dfs$CA.TERE, df_lines = catere_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "None")

catere <- sma(ga~SL * Region, data = z_spp_dfs$CA.TERE, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(catere, "CA.TERE Gape Area by Region")
catere_summ <- mk_spp_summary(catere, 5, grouping=T)
catere_graph_df <- mk_smaSPP_graph_df(catere_summ, 5, "Region")
catere_byReg <- 
mk_SMAplot(df_points = z_spp_dfs$CA.TERE, df_lines = catere_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "None")
catere_byReg <- 
mk_SMAfacets(df_points = z_spp_dfs$CA.TERE, df_lines = catere_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "dissected_by", 
	facetting = "Region")

catere <- sma(ga~SL * dissected_by, data = z_spp_dfs$CA.TERE, log = "xy", 
	method = "SMA", robust = T, slope.test = 2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(catere, "CA.TERE Gape Area by Dis")
catere_summ <- mk_spp_summary(catere, 14, grouping=T)
catere_graph_df <- mk_smaSPP_graph_df(catere_summ, 5, "dissected_by")
catere_byDis <-
mk_SMAplot(df_points = z_spp_dfs$CA.TERE, df_lines = catere_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")
catere_byDis <-
mk_SMAfacets(df_points = z_spp_dfs$CA.TERE, df_lines = catere_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")

dev.copy2pdf(device=quartz, file = "random_effects_plots/ca_tere_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/CaTere.pdf", width=7.8, height = 6.2)
catere_plot
catere_byReg
catere_byDis
dev.off()

#-------------------------------------------------------------------------------
chvand <- sma(ga~SL, data = z_spp_dfs$CH.VAND, log = "xy", method = "SMA", 
	robust = T, slope.test = 2)
check_assump(chvand, "CH.VAND Gape Area")
chvand_summ <- mk_sma_summary(chvand, "SpeciesCode")
chvand_graph_df <- mk_sma_graph_df(chvand_summ, 1, "SpeciesCode")
chvand_graph_df[, 1] <- pento_order[19]
chvand_plot <- 
mk_SMAplot(df_points = z_spp_dfs$CH.VAND, df_lines = chvand_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "None")

chvand <- sma(ga~SL + Region, data=b_spp_dfs$CH.VAND, log="xy", method="SMA", 
	robust = T, elev.com=T, multcomp = T, multcompmethod = "adjusted")


chvand <- sma(ga~SL * Region, data = z_spp_dfs$CH.VAND, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(chvand, "CH.VAND Gape Area by Region")
chvand_summ <- mk_spp_summary(chvand, 4, grouping=T)
chvand_graph_df <- mk_smaSPP_graph_df(chvand_summ, 4, "Region")
chvand_byReg <- 
mk_SMAplot(df_points = z_spp_dfs$CH.VAND, df_lines = chvand_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "None")
chvand_byReg <- 
mk_SMAfacets(df_points = z_spp_dfs$CH.VAND, df_lines = chvand_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "None", 
	facetting = "Region")

chvand <- sma(ga~SL + dissected_by, data = z_spp_dfs$CH.VAND, log = "xy", 
	method = "SMA", robust = T, elev.com = T, multcomp = T, 
	multcompmethod = "adjusted")

chvand <- sma(ga~SL * dissected_by, data = z_spp_dfs$CH.VAND, log = "xy", 
	method = "SMA", robust = T, slope.test = 2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(chvand, "CH.VAND Gape Area by Dis")
chvand_summ <- mk_spp_summary(chvand, 5, grouping=T)
chvand_graph_df <- mk_smaSPP_graph_df(chvand_summ, 4, "dissected_by")
chvand_byDis <-
mk_SMAplot(df_points = z_spp_dfs$CH.VAND, df_lines = chvand_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")
chvand_byDis <-
mk_SMAfacets(df_points = z_spp_dfs$CH.VAND, df_lines = chvand_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")

dev.copy2pdf(device=quartz, file = "random_effects_plots/ch_vand_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/ChVand.pdf", width=7.8, height = 6.2)
chvand_plot
chvand_byReg
chvand_byDis
dev.off()

#-------------------------------------------------------------------------------
psdisp <- sma(ga~SL, data = z_spp_dfs$PS.DISP, log = "xy", method = "SMA", 
	robust = T, slope.test = 2)
check_assump(psdisp, "PS.DISP Gape Area")
psdisp_summ <- mk_sma_summary(psdisp, "SpeciesCode")
psdisp_graph_df <- mk_sma_graph_df(psdisp_summ, 1, "SpeciesCode")
psdisp_graph_df[, 1] <- pento_order[21]
psdisp_plot <- 
mk_SMAplot(df_points = z_spp_dfs$PS.DISP, df_lines = psdisp_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "None")

psdisp <- sma(ga~SL + Region, data=b_spp_dfs$PS.DISP, log="xy", method="SMA", 
	robust = T, elev.com=T, multcomp = T, multcompmethod = "adjusted")


psdisp <- sma(ga~SL * Region, data = z_spp_dfs$PS.DISP, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(psdisp, "PS.DISP Gape Area by Region")
psdisp_summ <- mk_spp_summary(psdisp, 4, grouping=T)
psdisp_graph_df <- mk_smaSPP_graph_df(psdisp_summ, 4, "Region")
psdisp_byReg <- 
mk_SMAplot(df_points = z_spp_dfs$PS.DISP, df_lines = psdisp_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "None")
psdisp_byReg <- 
mk_SMAfacets(df_points = z_spp_dfs$PS.DISP, df_lines = psdisp_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "None", 
	facetting = "Region")

psdisp <- sma(ga~SL + dissected_by, data = z_spp_dfs$PS.DISP, log = "xy", 
	method = "SMA", robust = T, elev.com = T, multcomp = T, 
	multcompmethod = "adjusted")

psdisp <- sma(ga~SL * dissected_by, data = z_spp_dfs$PS.DISP, log = "xy", 
	method = "SMA", robust = T, slope.test = 2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(psdisp, "PS.DISP Gape Area by Dis")
psdisp_summ <- mk_spp_summary(psdisp, 4, grouping=T)
psdisp_graph_df <- mk_smaSPP_graph_df(psdisp_summ, 4, "dissected_by")
psdisp_byDis <-
mk_SMAplot(df_points = z_spp_dfs$PS.DISP, df_lines = psdisp_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")
psdisp_byDis <-
mk_SMAfacets(df_points = z_spp_dfs$PS.DISP, df_lines = psdisp_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")

dev.copy2pdf(device=quartz, file = "random_effects_plots/ps_disp_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/PsDisp.pdf", width=7.8, height = 6.2)
psdisp_plot
psdisp_byReg
psdisp_byDis
dev.off()

#-------------------------------------------------------------------------------
psbart <- sma(ga~SL, data = z_spp_dfs$PS.BART, log = "xy", method = "SMA", 
	robust = T, slope.test = 2)
#-------------------------------------------------------------------------------
psoliv <- sma(ga~SL, data = z_spp_dfs$PS.OLIV, log = "xy", method = "SMA", 
	robust = T, slope.test = 2)
check_assump(psoliv, "PS.OLIV Gape Area")
psoliv_summ <- mk_sma_summary(psoliv, "SpeciesCode")
psoliv_graph_df <- mk_sma_graph_df(psoliv_summ, 1, "SpeciesCode")
psoliv_graph_df[, 1] <- pento_order[22]
psoliv_plot <- 
mk_SMAplot(df_points = z_spp_dfs$PS.OLIV, df_lines = psoliv_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "None")

psoliv <- sma(ga~SL + Region, data=b_spp_dfs$PS.OLIV, log="xy", method="SMA", 
	robust = T, elev.com=T, multcomp = T, multcompmethod = "adjusted")


psoliv <- sma(ga~SL * Region, data = z_spp_dfs$PS.OLIV, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(psoliv, "PS.OLIV Gape Area by Region")
psoliv_summ <- mk_spp_summary(psoliv, 6, grouping=T)
psoliv_graph_df <- mk_smaSPP_graph_df(psoliv_summ, 5, "Region")
psoliv_byReg <- 
mk_SMAplot(df_points = z_spp_dfs$PS.OLIV, df_lines = psoliv_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "None")
psoliv_byReg <- 
mk_SMAfacets(df_points = z_spp_dfs$PS.OLIV, df_lines = psoliv_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "None", 
	facetting = "Region")

psoliv <- sma(ga~SL + dissected_by, data = z_spp_dfs$PS.OLIV, log = "xy", 
	method = "SMA", robust = T, elev.com = T, multcomp = T, 
	multcompmethod = "adjusted")

ao_adrian <- which(z_spp_dfs$PS.OLIV$dissected_by == 'AO/Adrian')
no_ao_adrian <- z_spp_dfs$PS.OLIV[-ao_adrian, ]

psoliv <- sma(ga~SL * dissected_by, data = no_ao_adrian, log = "xy", 
	method = "SMA", robust = T, slope.test = 2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(psoliv, "PS.OLIV Gape Area by Dis")
psoliv_summ <- mk_spp_summary(psoliv, 10, grouping=T)
psoliv_graph_df <- mk_smaSPP_graph_df(psoliv_summ, 8, "dissected_by")
psoliv_byDis <-
mk_SMAplot(df_points = z_spp_dfs$PS.OLIV, df_lines = psoliv_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")
psoliv_byDis <-
mk_SMAfacets(df_points = z_spp_dfs$PS.OLIV, df_lines = psoliv_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")


dev.copy2pdf(device=quartz, file = "random_effects_plots/ps_oliv_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/PsOliv.pdf", width=7.8, height = 6.2)
psoliv_plot
psoliv_byReg
psoliv_byDis
dev.off()

#===============================================================================
# Family and species level varation in the Herbivores
#===============================================================================
# All
zGA <- sma(ga~SL, data = z, log = "xy", method = "SMA", robust = T, slope.test = 2)
check_assump(zGA, "ZP Gape Area All")
zGA_summ <- mk_sma_summary(zGA, 1)
zGA_graph_df <- mk_sma_graph_df(zGA_summ, 1, "j_fg")
zGA_graph_df[, 1] <- j_fg[3]
zGA_all <- 
mk_SMAplot(df_points = z, df_lines = zGA_graph_df, gapeType = "ga", 
	point_colour = "j_fg", line_colour = "j_fg", labels = "None")

mk_SMAplot(df_points = z, df_lines = zGA_graph_df, gapeType = "ga", 
	point_colour = "j_fg", line_colour = "j_fg", labels = "SpecimenID")


zfam <- sma(ga~SL + Family, data = z, log = "xy", method = "SMA", robust = T, 
	slope.test = 2, multcomp = T, multcompmethod = "adjusted")

zfam <- sma(ga~SL * Family, data = z, log = "xy", method = "SMA", robust = T, 
	slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(zfam, "ZP Gape Area by Family")
zfam_summ <- mk_spp_summary(zfam, 3, grouping = T)
zfam_graph_df <- mk_smaSPP_graph_df(zfam_summ, 3, "Family")
zfam_plot <- 
mk_SMAplot(df_points = z, df_lines = zfam_graph_df, gapeType = "ga",
 point_colour = "Family", line_colour = "Family", labels = "None")
zfam_facets <- 
mk_SMAfacets(df_points = z, df_lines = zfam_graph_df, gapeType = "ga",
	point_colour = "SpeciesCode", labels = "None", facetting = "Family")

ddply(z, .(Family), summarise, length(SpecimenID))
ddply(z, .(SpeciesCode), summarise, length(SpecimenID))

zspp <- sma(ga~SL * SpeciesCode, data = z, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(zspp, "ZP Gape Area by Species")
zspp_summ <- mk_spp_summary(zspp, 6, grouping = T)
zspp_graph_df <- mk_smaSPP_graph_df(zspp_summ, 6, "SpeciesCode")
zspp_plot <- 
mk_SMAplot(df_points = z, df_lines = zspp_graph_df, gapeType = "ga",
 point_colour = "SpeciesCode", line_colour = "SpeciesCode", labels = "None")
zspp_facets <- 
mk_SMAfacets(df_points = z, df_lines = zfam_graph_df, gapeType = "ga",
	point_colour = "SpeciesCode", labels = "None", facetting = "Family")

# Obtaining overall group slope estimate:
with(z, slope.com(log(ga), log(SL), SpeciesCode, method = 'SMA', robust = T))

dev.copy2pdf(device=quartz, file = "random_effects_plots/zfam_zspp_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/zfam_zspp.pdf", width=7.8, height = 6.2)
zfam_plot
zfam_facets
zspp_plot
dev.off()

#-------------------------------------------------------------------------------
# Herbivores

hGA <- sma(ga~SL, data = h, log = "xy", method = "SMA", robust = T, slope.test = 2)
check_assump(hGA, "He Gape Area All")
hGA_summ <- mk_sma_summary(hGA, "j_fg")
hGA_graph_df <- mk_sma_graph_df(hGA_summ, 1, "j_fg")
hGA_graph_df[, 1] <- j_fg[4]
hGA_plot <- 
mk_SMAplot(df_points = h, df_lines = hGA_graph_df, gapeType = "ga", 
	point_colour = "j_fg", line_colour = "j_fg", labels = "None")

with(h, slope.com(log(ga), log(SL), Region, method = 'SMA', robust = T, 
	slope.test = 2))

hGA <- sma(ga~SL * Region, data = h, log = "xy", method = "SMA", robust = T,
	slope.test = 2, multcomp = TRUE, multcompmethod = "adjusted")
check_assump(hGA, "He Gape Area by Regions")
hGA_regSumm <- mk_spp_summary(hGA, 6, grouping = T)
hGA_reg_graphing <- mk_smaSPP_graph_df(hGA_regSumm, 5, "Region")
hReg_by_spp <- 
mk_SMAplot(df_points = h, df_lines = hGA_reg_graphing, gapeType = "ga", 
	point_colour = "Region", line_colour = "Region", labels = "None")
mk_SMAfacets(df_points = h, df_lines = hGA_reg_graphing, gapeType = "ga",
	point_colour = "SpeciesCode", labels = "None", facetting = "Region")

jb <- which(h$dissected_by == 'JB')
no_jb <- h[-jb, ]

hGA <- sma(ga~SL * dissected_by, data = no_jb, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = TRUE, multcompmethod = "adjusted")
check_assump(hGA, "He Gape Area dissected_by")
hGA_disSumm <- mk_spp_summary(hGA, 11, grouping = T)
hGA_dis_graphing <- mk_smaSPP_graph_df(hGA_disSumm, 11, "dissected_by")
hDis_by_spp <- 
mk_SMAplot(df_points = h, df_lines = hGA_dis_graphing, gapeType = "ga", 
	point_colour = "dissected_by", line_colour = "dissected_by", 
	labels = "None")
mk_SMAfacets(df_points = h, df_lines = hGA_dis_graphing, gapeType = "ga",
	point_colour = "SpeciesCode", labels = "None", facetting = "dissected_by")
	

pdf(file = "random_effects_plots/hGA.pdf", width=7.8, height = 6.2)
hGA_plot
hReg_by_spp
hDis_by_spp
dev.off()

dev.copy2pdf(device=quartz, file = "random_effects_plots/hGA_residuals.pdf", 
	width=6.7, height=10.8)

#-------------------------------------------------------------------------------
# Makes list of He species dataframes (subsets of fish)
h_spp_dfs <- split(h, h$SpeciesCode, drop=TRUE)

acnigr <- sma(ga~SL, data = h_spp_dfs$AC.NIGR, log = "xy", method = "SMA", 
	robust = T, slope.test = 2)
check_assump(acnigr, "AC.NIGR Gape Area")
acnigr_summ <- mk_sma_summary(acnigr, "SpeciesCode")
acnigr_graph_df <- mk_sma_graph_df(acnigr_summ, 1, "SpeciesCode")
acnigr_graph_df[, 1] <- pento_order[11]
acnigr_plot <- 
mk_SMAplot(df_points = h_spp_dfs$AC.NIGR, df_lines = acnigr_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "None")

hphf <- which(h_spp_dfs$AC.NIGR$Region == 'HP.HF')
no_hphf <- h_spp_dfs$AC.NIGR[-hphf, ]

acnigr <- sma(ga~SL * Region, data = no_hphf, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(acnigr, "AC.NIGR Gape Area by Region")
acnigr_summ <- mk_spp_summary(acnigr, 5, grouping=T)
acnigr_graph_df <- mk_smaSPP_graph_df(acnigr_summ, 3, "Region")
acnigr_byReg <- 
mk_SMAplot(df_points = h_spp_dfs$AC.NIGR, df_lines = acnigr_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "None")
acnigr_byReg <- 
mk_SMAfacets(df_points = h_spp_dfs$AC.NIGR, df_lines = acnigr_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "dissected_by", 
	facetting = "Region")

acnigr <- sma(ga~SL * dissected_by, data = h_spp_dfs$AC.NIGR, log = "xy", 
	method = "SMA", robust = T, slope.test = 2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(acnigr, "AC.NIGR Gape Area by Dis")
acnigr_summ <- mk_spp_summary(acnigr, 7, grouping=T)
acnigr_graph_df <- mk_smaSPP_graph_df(acnigr_summ, 4, "dissected_by")
acnigr_byDis <-
mk_SMAplot(df_points = h_spp_dfs$AC.NIGR, df_lines = acnigr_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")
acnigr_byDis <-
mk_SMAfacets(df_points = h_spp_dfs$AC.NIGR, df_lines = acnigr_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")

dev.copy2pdf(device=quartz, file = "random_effects_plots/ac_nigr_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/AcNigr.pdf", width=7.8, height = 6.2)
acnigr_plot
acnigr_byReg
acnigr_byDis
dev.off()

#-------------------------------------------------------------------------------
# AC.OlIV not analysed for differences across region or observer
#-------------------------------------------------------------------------------
ceflav <- sma(ga~SL, data = h_spp_dfs$CE.FLAV, log = "xy", method = "SMA", 
	robust = T, slope.test = 2)
check_assump(ceflav, "CE.FLAV Gape Area")
ceflav_summ <- mk_sma_summary(ceflav, "SpeciesCode")
ceflav_graph_df <- mk_sma_graph_df(ceflav_summ, 1, "SpeciesCode")
ceflav_graph_df[, 1] <- pento_order[13]
ceflav_plot <- 
mk_SMAplot(df_points = h_spp_dfs$CE.FLAV, df_lines = ceflav_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "dissected_by")

ceflav <- sma(ga~SL * Region, data = h_spp_dfs$CE.FLAV, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(ceflav, "CE.FLAV Gape Area by Region")
ceflav_summ <- mk_spp_summary(ceflav, 4, grouping=T)
ceflav_graph_df <- mk_smaSPP_graph_df(ceflav_summ, 4, "Region")
ceflav_byReg <- 
mk_SMAplot(df_points = h_spp_dfs$CE.FLAV, df_lines = ceflav_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "dissected_by")
ceflav_byReg <- 
mk_SMAfacets(df_points = h_spp_dfs$CE.FLAV, df_lines = ceflav_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "dissected_by", 
	facetting = "Region")

ceflav <- sma(ga~SL * dissected_by, data = h_spp_dfs$CE.FLAV, log = "xy", 
	method = "SMA", robust = T, slope.test = 2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(ceflav, "CE.FLAV Gape Area by Dis")
ceflav_summ <- mk_spp_summary(ceflav, 7, grouping=T)
ceflav_graph_df <- mk_smaSPP_graph_df(ceflav_summ, 3, "dissected_by")
ceflav_byDis <-
mk_SMAplot(df_points = h_spp_dfs$CE.FLAV, df_lines = ceflav_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")
ceflav_byDis <-
mk_SMAfacets(df_points = h_spp_dfs$CE.FLAV, df_lines = ceflav_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")

dev.copy2pdf(device=quartz, file = "random_effects_plots/ce_flav_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/CeFlav.pdf", width=7.8, height = 6.2)
ceflav_plot
ceflav_byReg
ceflav_byDis
dev.off()

#-------------------------------------------------------------------------------
chsord <- sma(ga~SL, data = h_spp_dfs$CH.SORD, log = "xy", method = "SMA", 
	robust = T, slope.test = 2)
check_assump(chsord, "CH.SORD Gape Area")
chsord_summ <- mk_sma_summary(chsord, "SpeciesCode")
chsord_graph_df <- mk_sma_graph_df(chsord_summ, 1, "SpeciesCode")
chsord_graph_df[, 1] <- pento_order[14]
chsord_plot <- 
mk_SMAplot(df_points = h_spp_dfs$CH.SORD, df_lines = chsord_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "None")

chsord <- sma(ga~SL * Region, data = h_spp_dfs$CH.SORD, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(chsord, "CH.SORD Gape Area by Region")
chsord_summ <- mk_spp_summary(chsord, 5, grouping=T)
chsord_graph_df <- mk_smaSPP_graph_df(chsord_summ, 5, "Region")
chsord_byReg <- 
mk_SMAplot(df_points = h_spp_dfs$CH.SORD, df_lines = chsord_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "None")
chsord_byReg <- 
mk_SMAfacets(df_points = h_spp_dfs$CH.SORD, df_lines = chsord_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "dissected_by", 
	facetting = "Region")

chsord <- sma(ga~SL * dissected_by, data = h_spp_dfs$CH.SORD, log = "xy", 
	method = "SMA", robust = T, slope.test = 2, multcomp = T, 
	multcompmethod = "adjusted")
check_assump(chsord, "CH.SORD Gape Area by Dis")
chsord_summ <- mk_spp_summary(chsord, 6, grouping=T)
chsord_graph_df <- mk_smaSPP_graph_df(chsord_summ, 6, "dissected_by")
chsord_byDis <-
mk_SMAplot(df_points = h_spp_dfs$CH.SORD, df_lines = chsord_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")
chsord_byDis <-
mk_SMAfacets(df_points = h_spp_dfs$CH.SORD, df_lines = chsord_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")

dev.copy2pdf(device=quartz, file = "random_effects_plots/ch_sord_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/ChSord.pdf", width=7.8, height = 6.2)
chsord_plot
chsord_byReg
chsord_byDis
dev.off()

#-------------------------------------------------------------------------------
scfren <- sma(ga~SL, data = h_spp_dfs$SC.FREN, log = "xy", method = "SMA", 
	robust = T, slope.test = 2)
check_assump(scfren, "SC.FREN Gape Area")
scfren_summ <- mk_sma_summary(scfren, "SpeciesCode")
scfren_graph_df <- mk_sma_graph_df(scfren_summ, 1, "SpeciesCode")
scfren_graph_df[, 1] <- pento_order[15]
scfren_plot <- 
mk_SMAplot(df_points = h_spp_dfs$SC.FREN, df_lines = scfren_graph_df,
	gapeType = "ga", point_colour = "SpeciesCode", line_colour = "SpeciesCode", 
	labels = "dissected_by")

kif12_048 <- which(h_spp_dfs$SC.FREN$SpecimenID == 'KIF12_048')
no_kif12_048 <- h_spp_dfs$SC.FREN[-kif12_048, ]

scfren <- sma(ga~SL * Region, data = no_kif12_048, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(scfren, "SC.FREN Gape Area by Region")
scfren_summ <- mk_spp_summary(scfren, 5, grouping=T)
scfren_graph_df <- mk_smaSPP_graph_df(scfren_summ, 4, "Region")
scfren_byReg <- 
mk_SMAplot(df_points = h_spp_dfs$SC.FREN, df_lines = scfren_graph_df,
	gapeType = "ga", point_colour = "Region", line_colour = "Region",
	labels = "None")
scfren_byReg <- 
mk_SMAfacets(df_points = h_spp_dfs$SC.FREN, df_lines = scfren_graph_df,
	gapeType = "ga", point_colour = "dissected_by", labels = "dissected_by", 
	facetting = "Region")

scfren <- sma(ga~SL * dissected_by, data = h_spp_dfs$SC.FREN, log = "xy", 
	method = "SMA", robust = T, slope.test = 2), multcomp = T, 
	multcompmethod = "adjusted")
check_assump(scfren, "SC.FREN Gape Area by Dis")
scfren_summ <- mk_spp_summary(scfren, 6, grouping=T)
scfren_graph_df <- mk_smaSPP_graph_df(scfren_summ, 6, "dissected_by")
scfren_byDis <-
mk_SMAplot(df_points = h_spp_dfs$SC.FREN, df_lines = scfren_graph_df,
	gapeType = "ga", point_colour = "dissected_by", line_colour = "dissected_by",
	labels = "None")
scfren_byDis <-
mk_SMAfacets(df_points = h_spp_dfs$SC.FREN, df_lines = scfren_graph_df,
	gapeType = "ga", point_colour = "Region", labels = "None", 
	facetting = "dissected_by")

dev.copy2pdf(device=quartz, file = "random_effects_plots/sc_fren_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/ScFren.pdf", width=7.8, height = 6.2)
scfren_plot
scfren_byReg
scfren_byDis
dev.off()


#-------------------------------------------------------------------------------
# Too few data points to analyse region and observer effects in SC.RUBR
#-------------------------------------------------------------------------------

#===============================================================================
# Family and species level varation in the Herbivores
#===============================================================================
# All
hGA <- sma(ga~SL, data = h, log = "xy", method = "SMA", robust = T, slope.test = 2)
check_assump(hGA, "He Gape Area All")
hGA_summ <- mk_sma_summary(hGA, 1)
hGA_graph_df <- mk_sma_graph_df(hGA_summ, 1, "j_fg")
hGA_graph_df[, 1] <- j_fg[4]
hGA_all <- 
mk_SMAplot(df_points = h, df_lines = hGA_graph_df, gapeType = "ga", 
	point_colour = "j_fg", line_colour = "j_fg", labels = "None")


hfam <- sma(ga~SL + Family, data = h, log = "xy", method = "SMA", robust = T, 
	slope.test = 2, multcomp = T, multcompmethod = "adjusted")

hfam <- sma(ga~SL * Family, data = h, log = "xy", method = "SMA", robust = T, 
	slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(hfam, "He Gape Area by Family")
hfam_summ <- mk_spp_summary(hfam, 3, grouping = T)
hfam_graph_df <- mk_smaSPP_graph_df(hfam_summ, 3, "Family")
hfam_plot <- 
mk_SMAplot(df_points = h, df_lines = hfam_graph_df, gapeType = "ga",
 point_colour = "Family", line_colour = "Family", labels = "None")
hfam_facets <- 
mk_SMAfacets(df_points = h, df_lines = hfam_graph_df, gapeType = "ga",
	point_colour = "SpeciesCode", labels = "None", facetting = "Family")

ddply(h, .(Family), summarise, length(SpecimenID))
ddply(h, .(SpeciesCode), summarise, length(SpecimenID))

with(h, slope.com(log(ga), log(SL), groups = SpeciesCode, method = "SMA", 
	robust = TRUE, slope.test = 2, ci = TRUE))

hspp <- sma(ga~SL + SpeciesCode, data = h, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")


hspp <- sma(ga~SL * SpeciesCode, data = h, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(hspp, "He Gape Area by Species")
hspp_summ <- mk_spp_summary(hspp, 6, grouping = T)
hspp_graph_df <- mk_smaSPP_graph_df(hspp_summ, 6, "SpeciesCode")
hspp_plot <- 
mk_SMAplot(df_points = h, df_lines = hspp_graph_df, gapeType = "ga",
 point_colour = "SpeciesCode", line_colour = "SpeciesCode", labels = "None")
hspp_facets <- 
mk_SMAfacets(df_points = h, df_lines = hfam_graph_df, gapeType = "ga",
	point_colour = "SpeciesCode", labels = "None", facetting = "Family")

# Obtaining overall group slope estimate:
with(h, slope.com(log(ga), log(SL), SpeciesCode, method = 'SMA', robust = T))

with(p, slope.com(log(ga), log(SL), SpeciesCode, method = 'SMA', robust = T))


dev.copy2pdf(device=quartz, file = "random_effects_plots/hfam_hspp_residuals.pdf", 
	width=6.7, height=10.8)

pdf(file = "random_effects_plots/hfam_hspp.pdf", width=7.8, height = 6.2)
hfam_plot
hfam_facets
hspp_plot
dev.off()


#===============================================================================
# Multipanel plot with functional group and species
#===============================================================================
with(pento, slope.com(log(ga), log(SL), j_fg, method = "SMA", robust = TRUE))

with(p, slope.com(log(ga), log(SL), SpeciesCode, method = "SMA", robust = TRUE))
with(b, slope.com(log(ga), log(SL), SpeciesCode, method = "SMA", robust = TRUE))
with(z, slope.com(log(ga), log(SL), SpeciesCode, method = "SMA", robust = TRUE))
with(h, slope.com(log(ga), log(SL), SpeciesCode, method = "SMA", robust = TRUE))
sma(ga ~ SL, data = c, log = "xy", method = "SMA", robust = T, slope.test = 2)

# SMA regression for entire community
allGA <- sma(ga~SL, data = pento, log = "xy", method = "SMA", robust = T, slope.test = 2)
check_assump(allGA, "Pento Gape Area All")
allGA_summ <- mk_sma_summary(allGA, 1)
allGA_graph_df <- mk_sma_graph_df(allGA_summ, 1, "j_fg")
allGA_all <- 
mk_SMAplot(df_points = pento, df_lines = allGA_graph_df, gapeType = "ga", 
	point_colour = "j_fg", line_colour = "j_fg", labels = "None")

# SMA regression for each functional group
all_fg_GA <- sma(ga~SL * j_fg, data = pento, log = "xy", method = "SMA", robust = T, 
	slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(allGA, "Pento Gape Area All")
all_fg_GA_summ <- mk_spp_summary(all_fg_GA, 5, grouping=TRUE)
all_fg_GA_graph_df <- mk_smaSPP_graph_df(all_fg_GA_summ, 5, "j_fg")
all_fg_GA_all <- 
mk_SMAplot(df_points = pento, df_lines = all_fg_GA_graph_df, gapeType = "ga", 
	point_colour = "j_fg", line_colour = "j_fg", labels = "None")

# SMA regression for each species
all_spp_GA <- sma(ga~SL * SpeciesCode, data = pento, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")
check_assump(allGA, "Species Gape Area All")
allGA_bySPP_summ <- mk_spp_summary(all_spp_GA, 22, grouping=TRUE)
allGA_bySPP_graph_df <- mk_smaSPP_graph_df(allGA_bySPP_summ, 22, "SpeciesCode")
allGA_all <- 
mk_SMAplot(df_points = pento, df_lines = allGA_graph_df, gapeType = "ga", 
	point_colour = "j_fg", line_colour = "j_fg", labels = "None")


pento_line <- allGA_graph_df
fg_lines <- allGA_graph_df
spp_lines <- allGA_bySPP_graph_df
spp_lines$colour <- 
	c("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7",
	  "#F8766D", "#C49A00", "#53B400",
	  "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF",
	  "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF",
	  "#F8766D"
	  )

# Not actually sure what this code was for anymore...
names(allGA_bySPP_summ)[1] = "SpeciesCode"
spp_sma_maxL$SpeciesCode <- factor(spp_sma_maxL$SpeciesCode, levels=pento_order)

spp_sma_maxL <- merge(allGA_bySPP_summ, maxL, by = "SpeciesCode")
sort_by_maxL <- spp_sma_maxL[order(spp_sma_maxL$maxL), ]
sort_by_maxL$y <- c(1:22)

ggplot(data = sort_by_maxL) +
	geom_point(aes(x = slope, y = y, colour = j_fg)) +
	geom_segment(aes(x = upper, xend = lower, y = y, yend = y)) +
	geom_vline(xintercept = 2) +
	geom_text()

#-------------------------------------------------------------------------------
# Functional groups
#-------------------------------------------------------------------------------
fg_sma_eqns <- write_group_sma_eqn(all_fg_GA_summ, all_fg_GA_summ$group)
names(fg_sma_eqns) <- c("j_fg", "eqn_r2", "eqn", "r2", "n")

spp_sma_eqns <- write_group_sma_eqn(allGA_bySPP_summ, allGA_bySPP_summ$group)
names(spp_sma_eqns) <- c("SpeciesCode", "eqn_r2", "eqn", "r2", "n")

#-------------------------------------------------------------------------------
# Within functional groups
#-------------------------------------------------------------------------------

p_spp_GA <- sma(ga~SL * SpeciesCode, data = p, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")

z_spp_GA <- sma(ga~SL * SpeciesCode, data = z, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")

h_spp_GA <- sma(ga~SL * SpeciesCode, data = h, log = "xy", method = "SMA", 
	robust = T, slope.test = 2, multcomp = T, multcompmethod = "adjusted")


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


dev.copy2pdf(device = quartz, file = "panel_plots/fg_plot.pdf")

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
# Thinking about leaving out... but then should drop from main pisc_plot
#lukasm <- 
#mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$LU.KASM, 
#	spp_line_df_row = spp_lines[4, ], ref_intercept_row = spp_lines$ref_intercept[4], 
#	eqn_df = spp_sma_eqns[3, ], eqn_x = 700, eqn_y = 50, r2_x = 700, r2_y = 97 ,
#n_x =7000, n_y = 670, #	x_axis_labels = FALSE, 	y_axis_labels = FALSE, 
#	fg_line_intercept = fg_lines$ref_intercept[1], 
#   x_axis_text = FALSE, y_axis_text = FALSE, plot_title = )
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
catere <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$CA.TERE, 
	spp_line_df_row = spp_lines[10, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[10, ], eqn_x = 250, eqn_y = 2, r2_x = 250, r2_y = 4.5, 
	n_x = 250, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[4], x_axis_text = FALSE, 
	y_axis_text = TRUE, plot_title = "Caesio teres")
pttile <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$PT.TILE, 
	spp_line_df_row = spp_lines[11, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[11, ], eqn_x = 250, eqn_y = 2, r2_x = 250, r2_y = 4.5, 
	n_x = 250, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[4], x_axis_text = FALSE, 
	y_axis_text = FALSE, plot_title = "Pterocaesio tile")
chvand <-
mk_multipanel_plots2(fg_point_df = z, spp_point_df = z_spp_dfs$CH.VAND, 
	spp_line_df_row = spp_lines[12, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[12, ], eqn_x = 250, eqn_y = 2, r2_x = 250, r2_y = 4.5, 
	n_x = 250, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[4], x_axis_text = FALSE, 
	y_axis_text = FALSE, plot_title = "Chromis vanderbilti")
psbart <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.BART, 
	spp_line_df_row = spp_lines[13, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[13, ], eqn_x = 250, eqn_y = 2, r2_x = 250, r2_y = 4.5, 
	n_x = 250, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[4], x_axis_text = TRUE, 
	y_axis_text = TRUE, plot_title = "Pseudanthias bartlettorum")
psdisp <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.DISP, 
	spp_line_df_row = spp_lines[14, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[14, ], eqn_x = 250, eqn_y = 2, r2_x = 250, r2_y = 4.5, 
	n_x = 250, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[4], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Pseudanthias dispar")
psoliv <-
mk_multipanel_plots2(fg_point_df  = z, spp_point_df  = z_spp_dfs$PS.OLIV, 
	spp_line_df_row = spp_lines[15, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[15, ], eqn_x = 250, eqn_y = 2, r2_x = 250, r2_y = 4.5, 
	n_x = 250, n_y = 9, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[4], x_axis_text = TRUE, 
	y_axis_text = FALSE, plot_title = "Pseudanthias olivaceus")

# Plotting multipanel piscivores and benthic invertivores	
#dev.new(height = 10, width = 7)
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
print(catere, vp = set_vp(6, 2))
print(pttile, vp = set_vp(6, 3))
print(chvand, vp = set_vp(6, 4))
print(psbart, vp = set_vp(7, 2))
print(psdisp, vp = set_vp(7, 3))
print(psoliv, vp = set_vp(7, 4))
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

dev.copy2pdf(device = quartz, file = "panel_plots/pisc_panel.pdf")


# Herbivores
#dev.copy2pdf(device = quartz, file = "panel_plots/zoop_panel.pdf")
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
#dev.new(height = 4, width = 7)
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


chorna <-
mk_multipanel_plots2(fg_point_df  = c, spp_point_df  = c, 
	spp_line_df_row = spp_lines[22, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[22, ], eqn_x = 145, eqn_y = 18, r2_x = 145, r2_y = 21,
	x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[5], x_axis_text = TRUE, 
	y_axis_text = TRUE, plot_title = "C. ornatissimus")

dev.copy2pdf(device = quartz, file = "panel_plots/coral_panel.pdf")


camela <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CA.MELA, 
	spp_line_df_row = spp_lines[1, ], #ref_intercept_row = spp_lines$ref_intercept[1],
	eqn_df = spp_sma_eqns[1, ], eqn_x = 700, eqn_y = 50, r2_x = 700, r2_y = 97, 
	n_x = 700, n_y = 170, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[1], 
	x_axis_text = FALSE, y_axis_text = TRUE, plot_title = "Caranx melampygus") +
	geom_point(data = p_spp_dfs$CA.MELA, aes(colour = dissected_by), size = 3) +
	scale_fill_discrete("Observer", labels = p_spp_dfs$dissected_by)
apfurc <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$AP.FURC, 
	spp_line_df_row = spp_lines[2, ], #ref_intercept_row = spp_lines$ref_intercept[2], 
	eqn_df = spp_sma_eqns[2, ], eqn_x = 700, eqn_y = 50, r2_x = 700, r2_y = 97, 
	n_x = 700, n_y = 170, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[1], 
	x_axis_text = FALSE, y_axis_text = FALSE, plot_title = "Aphareus furca") +
	geom_point(data = p_spp_dfs$AP.FURC, aes(colour = dissected_by), size = 3)
luboha <- 
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$LU.BOHA, 
	spp_line_df_row = spp_lines[3, ], #ref_intercept_row = spp_lines$ref_intercept[3],
	eqn_df = spp_sma_eqns[3, ], eqn_x = 700, eqn_y = 50, r2_x = 700, r2_y = 97, 
	n_x = 700, n_y = 170, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[1], 
	x_axis_text = FALSE, y_axis_text = FALSE, plot_title = "Lutjanus bohar") +
	geom_point(data = p_spp_dfs$LU.BOHA, aes(colour = dissected_by), size = 3)
ceargu <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CE.ARGU, 
	spp_line_df_row = spp_lines[4, ], #ref_intercept_row = spp_lines$ref_intercept[5], 
	eqn_df = spp_sma_eqns[4, ], eqn_x = 700, eqn_y = 50, r2_x = 700, r2_y = 97, 
	n_x = 700, n_y = 170, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[1], 
	x_axis_text = TRUE, y_axis_text = TRUE, plot_title = "Cephalopholis argus") +
	geom_point(data = p_spp_dfs$CE.ARGU, aes(colour = dissected_by), size = 3)
ceurod <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$CE.UROD, 
	spp_line_df_row = spp_lines[5, ], #ref_intercept_row = spp_lines$ref_intercept[6], 
	eqn_df = spp_sma_eqns[5, ], eqn_x = 700, eqn_y = 50, r2_x = 700, r2_y = 97, 
	n_x = 700, n_y = 170, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[1], 
	x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Cephalopholis urodeta") +
	geom_point(data = p_spp_dfs$CE.UROD, aes(colour = dissected_by), size = 3)
valout <-
mk_multipanel_plots2(fg_point_df = p, spp_point_df = p_spp_dfs$VA.LOUT, 
	spp_line_df_row = spp_lines[6, ], #ref_intercept_row = spp_lines$ref_intercept[7],
	eqn_df = spp_sma_eqns[6, ], eqn_x = 700, eqn_y = 50, r2_x = 700, r2_y = 97, 
	n_x = 700, n_y = 170, x_axis_labels = FALSE, y_axis_labels = FALSE, 
	fg_line_intercept = fg_lines$ref_intercept[1], 
	x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Variola louti") +
	geom_point(data = p_spp_dfs$VA.LOUT, aes(colour = dissected_by), size = 3)



grid.arrange(paarca, mogran, painsu, ncol = 3,
	left = textGrob(expression(atop("", paste("log(gape area ", mm^2, ")", sep= ""))), 
		rot = 90), 
	sub = "log(standard length, mm) \n")

p_plot <- 
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = p, 
	spp_line_df_row = fg_lines[1, ], ref_intercept_row = fg_lines$ref_intercept[1])
b_plot <- 
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = b, 
	spp_line_df_row = fg_lines[2, ], ref_intercept_row = fg_lines$ref_intercept[2])
z_plot <-
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = z, 
	spp_line_df_row = fg_lines[3, ], ref_intercept_row = fg_lines$ref_intercept[3])
h_plot <- 
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = h, 
	spp_line_df_row = fg_lines[4, ], ref_intercept_row = fg_lines$ref_intercept[4])
c_plot <- 
mk_multipanel_plots2(fg_point_df = pento, spp_point_df = c, 
	spp_line_df_row = fg_lines[5, ], ref_intercept_row = fg_lines$ref_intercept[5])



ggplot(data = p_spp_dfs$CA.MELA, aes_string(x = "SL", y = "ga")) +
    geom_segment(data = spp_lines[1, ], aes_string(x = "from", xend = "to", 
     y = "yfrom", yend = "yto")) +
    geom_point( aes_string(colour = "SpeciesCode") ) +
    scale_y_log10() +
    scale_x_log10() +
    xlab("log(standard length, mm)") +     
    ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) + 
	geom_abline(intercept = spp_lines$ref_intercept[1], slope = 2, linetype = 2, 
		colour = "darkgrey") +
	theme_bw() +
	theme( plot.background = element_blank(), 
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.border = element_blank(), 
		panel.background = element_blank()
	) +
    theme(axis.line = element_line(color = 'black'))


m_x <- sqrt(fg_lines$from[1] * fg_lines$to[1])
m_y <- sqrt(fg_lines$yfrom[1] * fg_lines$yto[1])





















