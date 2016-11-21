library(lme4)
library(lmerTest)
library(multcomp)
library(beepr)
library(visreg)
library(dplyr)

# Fit random effects model with family as a random effect (phylogenetic non-independence)
#-------------------------------------------------------------------------------

boot_spp_summary_gh_mass <- 
  read.csv('gape_height_mass_species_bootstrapped_coefficients_10000.csv')

boot_spp_summary_gw_mass <- 
  read.csv('gape_width_mass_species_bootstrapped_coefficients_1000.csv')

spp_names <- c("Acanthurus nigricans", "Acanthurus olivaceus", "Aphareus furca", 
               "Lutjanus bohar", "Caesio teres", "Pterocaesio tile", 
               "Pseudanthias dispar", "Pseudanthias olivaceus", 
               "Cephalopholis argus", "Cephalopholis urodeta", "Variola louti", 
               "Chromis vanderbilti", "Monotaxis grandoculis", 
               "Chlorurus sordidus", "Scarus frenatus", "Scarus rubroviolaceus", 
               "Paracirrhites arcatus", "Caranx melampygus", 
               "Parupeneus insularis", "Centropyge flavissima")

spp_codes <- c('AC.NIGR', 'AC.OLIV', 'AP.FURC', 'LU.BOHA', 'CA.TERE', 
               'PT.TILE', 'PS.DISP', 'PS.OLIV', 'CE.ARGU', 'CE.UROD', 
               'VA.LOUT', 'CH.VAND', 'MO.GRAN', 'CH.SORD', 'SC.FREN', 
               'SC.RUBR', 'PA.ARCA', 'CA.MELA', 'PA.INSU', 'CE.FLAV')
# Create lookup dataframes
fg_lookup <- unique(data.frame('codes' = pento$Species, 'fg' = pento$j_fg))
spp_lookup_df <- data.frame('species' = spp_names, 'codes' = spp_codes)
fam_lookup <- ddply(fish, .(SpeciesCode), summarise, 'Family' = unique(Family), 
                    'Order' = unique(Order))

spp_lookup_df <- left_join(spp_lookup_df, fg_lookup) %>% 
                 left_join(., fam_lookup, by = c('codes' = 'SpeciesCode'))


# Gape height
all_spp_summ_df_gh <- 
  left_join(boot_spp_summary_gh_mass, spp_lookup_df, by = c('group' = 'codes')) %>% 
  as_data_frame() %>% 
  rename(codes = group) %>% 
  filter(!(codes %in% c('CH.ORNA', 'PS.BART')))

mean_spp_summ_gh <- ddply(all_spp_summ_df_gh, .(codes), summarise, 
                       'slope' = mean(Slope), 
                       'fg' = unique(fg), 
                       'species' = unique(species), 
                       'se' = sd(Slope),
                       'family' = unique(Family))

# Gape width
all_spp_summ_df_gw <- 
  left_join(boot_spp_summary_gw_mass, spp_lookup_df, by = c('group' = 'codes')) %>% 
  as_data_frame() %>% 
  rename(codes = group) %>% 
  filter(!(codes %in% c('CH.ORNA', 'PS.BART')))

mean_spp_summ_gw <- ddply(all_spp_summ_df_gw, .(codes), summarise, 
                       'slope' = mean(Slope), 
                       'fg' = unique(fg), 
                       'species' = unique(species), 
                       'se' = sd(Slope),
                       'family' = unique(Family))

test_gh <- lme4::lmer(slope ~ fg + 0 + (1 | family), data = mean_spp_summ_gh)
test_gh_ci <- confint(test_gh)

test_gh_summ <- data_frame(fg = factor(c('Pi', 'BI', 'ZP', 'He'), levels = c('Pi', 'BI', 'ZP', 'He')), 
                        estimate = summary(test_gh)$coefficients[, 1],
                        lwr = test_gh_ci[-(1:2), 1], 
                        upr = test_gh_ci[-(1:2), 2])

mean_spp_summ_gh2 <- 
  mean_spp_summ_gh %>% 
  mutate(id = as.numeric(codes)) %>% 
  mutate(fg = factor(fg, levels = levels(mean_spp_summ_gh$fg))) %>% 
  arrange(., fg)

spacing <- c(seq(0, 1, length.out = 6), seq(0, 1, length.out = 3), 
             seq(0, 1, length.out = 5), seq(0, 1, length.out = 6))
mean_spp_summ_gh2$spacing <- spacing

gh_plot <- 
ggplot(mean_spp_summ_gh) +
  geom_rect(data = test_gh_summ, aes(xmin = -0.05, xmax = 1.05, ymin = lwr, ymax = upr), fill = "gray85") +
  #geom_point(aes(x = -0.2, y = 1), alpha = 0) +
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
  geom_segment(data=distinct(test_gh_summ, fg, .keep_all = TRUE), aes(y = estimate, yend = estimate, x = -0.05, xend = 1.05), colour = '#099DFFFF', size = 1, lineend = 'round') +
  geom_hline(aes(yintercept = 1 / 3), colour = 'red', linetype = 2) +
  scale_x_continuous(breaks = 0.5) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), lim =c(0, 0.6)) +
  facet_grid(. ~ fg) 

test_gw <- lme4::lmer(slope ~ fg + 0 + (1 | family), data = mean_spp_summ_gw)
test_gw_ci <- confint(test_gw)

test_gw_summ <- data_frame(fg = factor(c('Pi', 'BI', 'ZP', 'He'), 
  levels = c('Pi', 'BI', 'ZP', 'He')), 
                        estimate = summary(test_gw)$coefficients[, 1],
                        lwr = test_gw_ci[-(1:2), 1], 
                        upr = test_gw_ci[-(1:2), 2])

mean_spp_summ_gw2 <- 
  mean_spp_summ_gw %>% 
  mutate(id = as.numeric(codes)) %>% 
  mutate(fg = factor(fg, levels = levels(mean_spp_summ_gw$fg))) %>% 
  arrange(., fg)

spacing <- c(seq(0, 1, length.out = 6), seq(0, 1, length.out = 3), 
             seq(0, 1, length.out = 5), seq(0, 1, length.out = 6))
mean_spp_summ_gw2$spacing <- spacing

gw_plot <- 
ggplot(mean_spp_summ_gw) +
  geom_rect(data = test_gw_summ, aes(xmin = -0.05, xmax = 1.05, ymin = lwr, ymax = upr), fill = "gray85") +
  #geom_point(aes(x = -0.2, y = 1), alpha = 0) +
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
  geom_segment(data=distinct(test_gw_summ, fg, .keep_all = TRUE), aes(y = estimate, yend = estimate, x = -0.05, xend = 1.05), colour = '#099DFFFF', size = 1, lineend = 'round') +
  geom_hline(aes(yintercept = 1 / 3), colour = 'red', linetype = 2) +
  scale_x_continuous(breaks = 0.5) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), lim =c(0, 0.7)) +
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
grid.text("Figure 1", vp = viewport(layout.pos.row = 4, layout.pos.col = 1:3),
    gp = gpar(fontsize = 9), vjust = 0, hjust = 1.3)

#dev.copy2eps(device = quartz, file = "panel_plots/DunicBaum_S2.eps")