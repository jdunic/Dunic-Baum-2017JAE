################################################################################
############            Gape Size ~ Body Size Analysis              ############
################################################################################

# Functional group level linear models:
fg_gh <- groupwise_lm_gh(pento, pento$j_fg)
fg_gw <- groupwise_lm_gw(pento, pento$j_fg)
fg_ga <- groupwise_lm_ga(pento, pento$j_fg)

fg_list <- rbind(fg_gh, fg_gw, fg_ga)

write.csv(fg_list, file = "fg_list.csv")

# Piscivore analysis:
# Linear mixed effects models for Functional groups:
library(lme4)
fm_p <- lmer(log(gh) ~ log(SLMM) + (1|SpeciesCode), p)
fm_p

fm_p2 <- lmer(log(gh) ~ log(SLMM) + (1|Family), p)
fm_p2

test <- lm(log(gh) ~ log(SLMM), data = p)
summary(test)

lm(log(gh) ~ log(SLMM), p)


x <- lme(log(gw) ~ log(SLMM), data = zp, random = ~ 1 | SpeciesCode)
y <- lme(log(gh) ~ log(SLMM), data = p, random = ~ 1 | Family)

summary(x)
intervals(x)
VarCorr(x)

# Level: Functional group - Piscivores ####

p1 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gh), colour=SpeciesCode)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
#  geom_text(data=data.frame(), aes(x=4.2, y=4.5, label=write_lme_gen(p, p$gh)), parse=TRUE) +
  ggtitle("Vertical Gape")

p2 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.2, y=4.5, label=write_lme_gen(p, p$gw)), parse=TRUE) +
  ggtitle("Horizontal Gape")

p3 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(ga))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.2, y=8.2, label=write_lme_gen(p, p$ga)), parse=TRUE) +
  ggtitle("Gape Area")

p_all <- grid.arrange(p1, p2, p3)

# Level: Family - Piscivores ####

pfam_gh <- groupwise_lm_gh(p, p$Family)
pfam_gw <- groupwise_lm_gw(p, p$Family)
pfam_ga <- groupwise_lm_ga(p, p$Family)

plist <- rbind(pfam_gh, pfam_gw, pfam_ga)
write.csv(plist, file = "plist.csv")

f1 <- write_lme_groups(pfam_gh, pfam_gh$variable)
f2 <- write_lme_groups(pfam_gw, pfam_gw$variable)
f3 <- write_lme_groups(pfam_ga, pfam_ga$variable)

fnames <- c("Family", "eqn")
colnames(f1) <- fnames
colnames(f2) <- fnames
colnames(f3) <- fnames

fam1 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f1, aes(x=5.0, y=5.1, label=eqn), size=4, parse=TRUE) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.05, vjust= -1)) +
  xlab("log(standard length)") +
  ylab("log(vertical gape)") 
facet_grid(. ~ cyl, labeller = label_both)

fam2 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f2, aes(x=5.0, y=5.2, label=eqn), size=4, parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.05, vjust= -1)) +
  xlab("log(standard length)") +
  ylab("log(horizontal gape)")
facet_grid(. ~ cyl, labeller = label_both)

fam3 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(ga))) + #, colour=SpeciesCode)) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f3, aes(x=4.8, y=9.9, label=eqn), size=4, parse=TRUE) +
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.05, vjust= -1)) +
  xlab("log(standard length)") +
  ylab("log(gape area)")
facet_grid(. ~ cyl, labeller = label_both)

grid.arrange(fam1, fam2)
fam3
p_fam <- grid.arrange(fam1, fam2, fam3)

# Level: Species - Piscivores ####
pspp_gh <- groupwise_lm_gh(p, p$SpeciesCode)
pspp_gw <- groupwise_lm_gw(p, p$SpeciesCode)
pspp_ga <- groupwise_lm_ga(p, p$SpeciesCode)

s1 <- write_lme_groups(pspp_gh, pspp_gh$variable)
s2 <- write_lme_groups(pspp_gw, pspp_gw$variable)
s3 <- write_lme_groups(pspp_ga, pspp_ga$variable)

snames <- c("SpeciesCode", "eqn")
colnames(s1) <- snames
colnames(s2) <- snames
colnames(s3) <- snames

spp1 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=4, drop=F) +
  geom_smooth(method=lm) +
  geom_text(data=s1, aes(x=4.85, y=4.8, label=eqn), size=3.5, parse=T) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.05, vjust= -1)) +
  xlab("log(standard length)") +
  ylab("log(vertical gape)")
  facet_grid(. ~ cyl, labeller = label_both)

spp2 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=4) +
  geom_smooth(method=lm) +
  ylim(c(1,5.3)) +
  geom_text(data=s2, aes(x=4.85, y=5.1, label=eqn), size=3.5, parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.05, vjust= -1)) +
  xlab("log(standard length)") +
  ylab("log(horizontal gape)")
facet_grid(. ~ cyl, labeller = label_both)

grid.arrange(spp1, spp2)

spp3 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(ga))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=4) +
  geom_smooth(method=lm) +
  ylim(c(1.9, 10.2)) +
  geom_text(data=s3, aes(x=5.0, y=10, label=eqn), size = 3.5, parse=TRUE) + 
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.05, vjust= -1)) +
  xlab("log(standard length)") +
  ylab("log(gape area)")
facet_grid(. ~ cyl, labeller = label_both)

spp3

p_spp <- grid.arrange(spp1, spp2, spp3)

pdf(file = "all_p_plots.pdf", height = 13, width = 9)
grid.arrange(p1, p2, p3)
grid.arrange(fam1, fam2, fam3)
grid.arrange(spp1, spp2, spp3)
dev.off()

# groupwise_lm_gh(df, variable)
# Functional Group level gape size ~ standard length relationships
fg_lm_summary_GH <- groupwise_lm_gh(df=pento, variable=pento$j_fg)
fg_lm_summary_GW <- groupwise_lm_gw(df=pento, variable=pento$j_fg)
fg_lm_summary_GA <- groupwise_lm_ga(df=pento, variable=pento$j_fg)

# Family level gape size ~ standard length relationships
fam_lm_summary_GH <- groupwise_lm_gh(df=pento, variable=pento$Family)
fam_lm_summary_GW <- groupwise_lm_gw(df=pento, variable=pento$Family)
fam_lm_summary_GA <- groupwise_lm_ga(df=pento, variable=pento$Family)

# Species level gape size ~ standard length relationships
sp_lm_summary_GH <- groupwise_lm_gh(df=pento, variable=pento$SpeciesCode)
sp_lm_summary_GW <- groupwise_lm_gw(df=pento, variable=pento$SpeciesCode)
sp_lm_summary_GA <- groupwise_lm_ga(df=pento, variable=pento$SpeciesCode)


# Benthic Invertivore analysis: ####
# Level: Functional group - Benthic Invertivore ####

b1 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.5, y=3.7, label=write_lme_gen(b, b$gh)), parse=TRUE) +
  ggtitle("Vertical Gape")

b2 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.5, y=3.7, label=write_lme_gen(b, b$gw)), parse=TRUE) +
  ggtitle("Horizontal Gape")

b3 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(ga))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.5, y=7.5, label=write_lme_gen(b, b$ga)), parse=TRUE) +
  ggtitle("Gape Area")

grid.arrange(b1, b2, b3)

# Level: Species - Benthic Invertivores ####
bspp_gh <- groupwise_lm_gh(b, b$SpeciesCode)
bspp_gw <- groupwise_lm_gw(b, b$SpeciesCode)
bspp_ga <- groupwise_lm_ga(b, b$SpeciesCode)

s1 <- write_lme_groups(bspp_gh, bspp_gh$variable)
s2 <- write_lme_groups(bspp_gw, bspp_gw$variable)
s3 <- write_lme_groups(bspp_ga, bspp_ga$variable)

snames <- c("SpeciesCode", "eqn")
colnames(s1) <- snames
colnames(s2) <- snames
colnames(s3) <- snames


spp1 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3, drop=F) +
  geom_smooth(method=lm) +
  geom_text(data=s1, aes(x=4.85, y=4.1, label=eqn), size=4, parse=T) +
  ggtitle("VerticalGape")
facet_grid(. ~ cyl, labeller = label_both)

spp2 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  #  ylim(c(1,5.3)) +
  geom_text(data=s2, aes(x=4.85, y=3.9, label=eqn), size=4, parse=TRUE) +
  ggtitle("HorizontalGape")
facet_grid(. ~ cyl, labeller = label_both)

spp3 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(ga))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  #  ylim(c(1.9, 10.2)) +
  geom_text(data=s3, aes(x=4.85, y=7.7, label=eqn), size=4, parse=TRUE) + 
  ggtitle("Gape Area")
facet_grid(. ~ cyl, labeller = label_both)

grid.arrange(spp1, spp2, spp3)

pdf(file = "all_b_plots.pdf", height = 13, width = 9)
grid.arrange(b1, b2, b3)
grid.arrange(spp1, spp2, spp3)
dev.off()


# Zooplanktivore analysis: ####
# Level: Functional group - Zooplanktivores ####

zp1 <-
  ggplot(data=zp, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4, y=3.2, label=write_lme_gen(zp, zp$gh)), parse=TRUE) +
  ggtitle("Vertical Gape")

zp2 <-
  ggplot(data=zp, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4, y=3.2, label=write_lme_gen(zp, zp$gw)), parse=TRUE) +
  ggtitle("Horizontal Gape")

zp3 <-
  ggplot(data=zp, aes(x=log(SLMM), y=log(ga))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4, y=6, label=write_lme_gen(zp, zp$ga)), parse=TRUE) +
  ggtitle("Gape Area")

grid.arrange(zp1, zp2, zp3)

# Level: Family - Zooplanktivores ####

zpfam_gh <- groupwise_lm_gh(zp, zp$Family)
zpfam_gw <- groupwise_lm_gw(zp, zp$Family)
zpfam_ga <- groupwise_lm_ga(zp, zp$Family)

zplist <- rbind(zpfam_gh, zpfam_gw, zpfam_ga)
write.csv(zplist, file = "zplist.csv")

f1 <- write_lme_groups(zpfam_gh, zpfam_gh$variable)
f2 <- write_lme_groups(zpfam_gw, zpfam_gw$variable)
f3 <- write_lme_groups(zpfam_ga, zpfam_ga$variable)

fnames <- c("Family", "eqn")
colnames(f1) <- fnames
colnames(f2) <- fnames
colnames(f3) <- fnames

fam1 <-
  ggplot(data=zp, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f1, aes(x=4.2, y=3.4, label=eqn), size=4, parse=TRUE) +
  ggtitle("VerticalGape")
facet_grid(. ~ cyl, labeller = label_both)

fam1

fam2 <-
  ggplot(data=zp, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f2, aes(x=4.2, y=3.4, label=eqn), size=4, parse=TRUE) +
  ggtitle("Horizontal Gape")
facet_grid(. ~ cyl, labeller = label_both)

fam2

fam3 <-
  ggplot(data=zp, aes(x=log(SLMM), y=log(ga))) + #, colour=SpeciesCode)) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f3, aes(x=4.2, y=6.4, label=eqn), size=4, parse=TRUE) +
  ggtitle("Gape Area")
facet_grid(. ~ cyl, labeller = label_both)

fam3

grid.arrange(fam1, fam2, fam3)

# Level: Species - Zooplanktivores ####

zpspp_gh <- groupwise_lm_gh(zp, zp$SpeciesCode)
zpspp_gw <- groupwise_lm_gw(zp, zp$SpeciesCode)
zpspp_ga <- groupwise_lm_ga(zp, zp$SpeciesCode)

s1 <- write_lme_groups(zpspp_gh, zpspp_gh$variable)
s2 <- write_lme_groups(zpspp_gw, zpspp_gw$variable)
s3 <- write_lme_groups(zpspp_ga, zpspp_ga$variable)

snames <- c("SpeciesCode", "eqn")
colnames(s1) <- snames
colnames(s2) <- snames
colnames(s3) <- snames

spp1 <-
  ggplot(data=zp, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3, drop=F) +
  geom_smooth(method=lm) +
  geom_text(data=s1, aes(x=4.2, y=3.3, label=eqn), size=4, parse=T) +
  ggtitle("VerticalGape")
facet_grid(. ~ cyl, labeller = label_both)

spp1

spp2 <-
  ggplot(data=zp, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  geom_text(data=s2, aes(x=4.1, y=3.25, label=eqn), size=4, parse=TRUE) +
  ggtitle("HorizontalGape")
facet_grid(. ~ cyl, labeller = label_both)

spp2

spp3 <-
  ggplot(data=zp, aes(x=log(SLMM), y=log(ga))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  geom_text(data=s3, aes(x=4.1, y=6.4, label=eqn), size=4, parse=TRUE) + 
  ggtitle("Gape Area")
facet_grid(. ~ cyl, labeller = label_both)

spp3

grid.arrange(spp1, spp2, spp3)

pdf(file = "all_zp_plots.pdf", height = 13, width = 9)
grid.arrange(zp1, zp2, zp3)
grid.arrange(fam1, fam2, fam3)
grid.arrange(spp1, spp2, spp3)
dev.off()


# Herbivore analysis: ####
# Level: Functional group - Herbivore ####

h1 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +  
  geom_text(data=data.frame(), aes(x=4.4, y=3.8, label=write_lme_gen(h, h$gh)), parse=TRUE) +
  ggtitle("Vertical Gape")

h2 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.4, y=3.4, label=write_lme_gen(h, h$gw)), parse=TRUE) +
  ggtitle("Horizontal Gape")

h3 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(ga))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.4, y=6.7, label=write_lme_gen(h, h$ga)), parse=TRUE) +
  ggtitle("Gape Area")

grid.arrange(h1, h2, h3)

# Level: Family - Herbivores ####
hfam_gh <- groupwise_lm_gh(h, h$Family)
hfam_gw <- groupwise_lm_gw(h, h$Family)
hfam_ga <- groupwise_lm_ga(h, h$Family)

hlist <- rbind(hfam_gh, hfam_gw, hfam_ga)
write.csv(hlist, file = "hlist.csv")

f1 <- write_lme_groups(hfam_gh, hfam_gh$variable)
f2 <- write_lme_groups(hfam_gw, hfam_gw$variable)
f3 <- write_lme_groups(hfam_ga, hfam_ga$variable)

fnames <- c("Family", "eqn")
colnames(f1) <- fnames
colnames(f2) <- fnames
colnames(f3) <- fnames


fam1 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f1, aes(x=5, y=4.2, label=eqn), size=4, parse=TRUE) +
  ggtitle("VerticalGape")
facet_grid(. ~ cyl, labeller = label_both)

fam2 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f2, aes(x=5, y=3.7, label=eqn), size=4, parse=TRUE) +
  ggtitle("Horizontal Gape")
facet_grid(. ~ cyl, labeller = label_both)

fam3 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(ga))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f3, aes(x=5, y=7.5, label=eqn), size=4, parse=TRUE) +
  ggtitle("Gape Area")
facet_grid(. ~ cyl, labeller = label_both)

grid.arrange(fam1, fam2, fam3)

# Level: Species - Herbivores ####
hspp_gh <- groupwise_lm_gh(h, h$SpeciesCode)
hspp_gw <- groupwise_lm_gw(h, h$SpeciesCode)
hspp_ga <- groupwise_lm_ga(h, h$SpeciesCode)

s1 <- write_lme_groups(hspp_gh, hspp_gh$variable)
s2 <- write_lme_groups(hspp_gw, hspp_gw$variable)
s3 <- write_lme_groups(hspp_ga, hspp_ga$variable)

snames <- c("SpeciesCode", "eqn")
colnames(s1) <- snames
colnames(s2) <- snames
colnames(s3) <- snames

spp1 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  geom_text(data=s1, aes(x=4.9, y=4.15, label=eqn), size=4, parse=TRUE) +
  ggtitle("VerticalGape")
facet_grid(. ~ cyl, labeller = label_both)

spp2 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  #  ylim(c(1,5.3)) +
  geom_text(data=s2, aes(x=4.8, y=3.65, label=eqn), size=4, parse=TRUE) +
  ggtitle("HorizontalGape")
facet_grid(. ~ cyl, labeller = label_both)

spp3 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(ga))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  #  ylim(c(1.9, 10.2)) +
  geom_text(data=s3, aes(x=4.8, y=7.45, label=eqn), size=4, parse=TRUE) + 
  ggtitle("Gape Area")
facet_grid(. ~ cyl, labeller = label_both)

grid.arrange(spp1, spp2, spp3)

pdf(file = "all_h_plots.pdf", height = 13, width = 9)
grid.arrange(h1, h2, h3)
grid.arrange(fam1, fam2, fam3)
grid.arrange(spp1, spp2, spp3)
dev.off()


# Corallivore analysis: ####
# Level: Functional Group (only 1 species) - Corallivore ####
c1 <-
  ggplot(data=c, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.6, y=2.8, label=write_lme_gen(c, c$gh)), parse=TRUE) +
  ggtitle("Vertical Gape")

c2 <-
  ggplot(data=c, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.6, y=2.8, label=write_lme_gen(c, c$gw)), parse=TRUE) +
  ggtitle("Horizontal Gape")

c3 <-
  ggplot(data=c, aes(x=log(SLMM), y=log(ga))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.6, y=5.35, label=write_lme_gen(c, c$ga)), parse=TRUE) +
  ggtitle("Gape Area")

grid.arrange(c1, c2, c3)

pdf(file = "all_c_plots.pdf", height = 13, width = 9)
grid.arrange(c1, c2, c3)
dev.off()



################################################################################
############            Gape Size ~ Body Size Isometry              ############
################################################################################

# Linear regression values for all species in the functional groups: 
# Pi, BI, ZP, He, C ####

write.csv(groupwise_lm_gh(pento, pento$SpeciesCode),
            file = "5_fg_gh_lm.csv")

write.csv(groupwise_lm_gw(pento, pento$SpeciesCode),
          file = "5_fg_gw_lm.csv")

write.csv(groupwise_lm_ga(pento, pento$SpeciesCode),
          file = "5_fg_ga_lm.csv")



x <- groupwise_lm_gh(zp, zp$Family)
x

y <- groupwise_lm_gw(zp, zp$Family)
y

z <- groupwise_lm_ga(zp, zp$Family)
z

# finding some averages for each functional group:
ddply(pento, ~ j_fg, summarise, mean = mean(gh_ratio), sd = sd(gh_ratio))
ddply(pento, ~ j_fg, summarise, mean = mean(gw_ratio), sd = sd(gw_ratio))

ddply(pento, ~ Family, summarise, mean = mean(gh_ratio), sd = sd(gh_ratio))
ddply(pento, ~ Family, summarise, mean = mean(gw_ratio), sd = sd(gh_ratio))

ddply(pento, ~ SpeciesCode, summarise, mean = mean(gh_ratio), sd = sd(gh_ratio))
ddply(pento, ~ SpeciesCode, summarise, mean = mean(gw_ratio), sd = sd(gw_ratio))


# Setting levels to functional groups: 
plot(gh_ratio~j_fg, pento)
plot(gh_ratio~Family, pento)


ggplot(pento, aes(x=SpeciesCode, y=gh_ratio, fill=Family)) +
  geom_boxplot() +
  xlab("Species") +
  ylab("Average Relative Gape Height") +
  theme(axis.text.x = element_text(angle=45, colour="black", vjust=0.5)) +
  theme(axis.text.y = element_text(colour="black")) 
  facet_wrap(~j_fg, scale="free_x", drop=TRUE)

ggplot(pento, aes(x=SpeciesCode, y=gh_ratio, fill=j_fg)) +
  geom_boxplot() +
  xlab("Species") +
  ylab("Average Relative Horizontal Gape") +
  theme(axis.text.x = element_text(angle=45, colour="black", vjust=0.5)) +
  theme(axis.text.y = element_text(colour="black")) +
  theme(axis.title.y = element_text(vjust = 0.3))

ggplot(pento, aes(x=SpeciesCode, y=gw_ratio, fill=j_fg)) +
  geom_boxplot() +
  xlab("Species") +
  ylab("Average Relative Vertical Gape") +
  theme(axis.text.x = element_text(angle=45, colour="black", vjust=0.5)) +
  theme(axis.text.y = element_text(colour="black")) +
  theme(axis.title.y = element_text(vjust = 0.3))

shapes <- c(1, 19, 17, 8, 18, 7, 2)

ggplot(p, aes(x=SLMM, y=gh_ratio, shape=SpeciesCode, colour=SpeciesCode)) +
  geom_point(size=3) +
  #scale_shape_manual(values=shapes) +
  geom_smooth(method=lm) 

ggplot(p, aes(x=SLMM, y=gh_ratio, shape=SpeciesCode, colour=SpeciesCode)) +
  geom_point(size=3) +
  scale_shape_manual(values=shapes) +
  geom_smooth(method=lm)


################################################################################
############            Predator ~ Prey Size Analysis               ############
################################################################################

# Actually doing the quantile regression:
# Have used the 10th and 90th percentiles because of sample size
# and following the general guidelines of n > (10/q) given in 
# Scharf et al, 1998 - Ecology

library('quantreg')
library('ggplot2')
library(gridExtra)
library(plyr)

# prey3 is a data.frame with just Pi, BI, ZP.
sl <- ggplot(prey2, aes(x = sl, y = psize)) +
  geom_point(aes(colour = ptype)) +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq")

gh <- ggplot(na.omit(prey2), aes(x = gh, y = psize)) +
  geom_point(aes(colour = ptype)) +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq")

grid.arrange(sl, gh)

df.n <- ddply(.data=prey3, .(fg), summarize, n=paste("n ==", length(fg)))

sl <- ggplot(prey3, aes(x = sl, y = psize)) +
  geom_point(aes(colour = ptype)) +
  geom_text(data = df.n, aes( x= 200, y = 230, label = n), parse = TRUE) +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq") +
  ylab("Prey total length (mm)") +
  xlab("Predator standard length (mm)") +
  labs(title = "a)") +
  theme(plot.title = element_text(hjust= -0.05, vjust= -1)) +
  #annotation_custom(
   # grob = textGrob(label = c("a)", "b)"), hjust = 7, vjust = -7,  gp = gpar(cex=1.5)),
    #ymin = -10,
    #ymax = 300,
    #xmin = -10, 
    #xmax = 600) +
  facet_wrap(~ fg) #, scales="free_x")

sl


df.n <- ddply(.data=na.omit(prey3), .(fg), summarize, n=paste("n ==", length(fg)))

gh <- ggplot(na.omit(prey3), aes(x = gh, y = psize)) +
  geom_point(aes(colour = ptype)) +
  geom_text(data = df.n, aes( x= 30, y = 230, label = n), parse = TRUE) +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq") +
  facet_wrap(~ fg,)
gh

gw <- ggplot(na.omit(prey3), aes(x = gw, y = psize)) +
  geom_point(aes(colour = ptype)) +
  geom_text(data = df.n, aes( x= 30, y = 230, label = n), parse = TRUE) +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq") +
  facet_wrap(~ fg)

grid.arrange(sl, gh) #, gw)

df.n <- ddply(.data=prey3, .(fg), summarize, n=paste("n ==", length(fg)))

ggplot(na.omit(prey3), aes(x = sl, y = psize)) +
  geom_point(aes(colour = species)) +
  geom_text(data = df.n, aes( x= 200, y = 200, label = n), parse = TRUE) +
  geom_smooth(method = "lm") +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.90), method = "rq") +
  facet_wrap(~ fg)

bi_qr <- rq(data = subset(prey3, prey3$fg == 'BI'), psize ~ sl, 
            tau = c(0.10, 0.50, 0.90))

p <- subset(prey3, prey3$fg == 'Pi')
b <- subset(prey3, prey3$fg == 'BI')

p <- subset(prey_gh, prey_gh$fg == 'Pi')
b <- subset(prey_gh, prey_gh$fg == 'BI')

summary(bi_qr, se = "rank")
bi_qr[[3]]

# Counting prey frequencies for Pi and BI species used in pred-prey analysis
# (in honours thesis) ####

colnames(myspp)
p <- subset(myspp, fg == 'Pi')
cnt_p <- length(p$SpeciesCode)

p_stuff <- subset(p, (invert != 0 | fish != 0))

p_empty <- subset(p, (invert == 0 & fish == 0))
cnt_pe <- length(p_empty$SpeciesCode)

p_invert <- subset(p, (invert == 1 & fish == 0))
cnt_pi <- length(p_invert$SpeciesCode)

p_fish <- subset(p, (invert == 0 & fish == 1))
cnt_pf <- length(p_fish$SpeciesCode)

p_both <- subset(p, (invert == 1 & fish ==1))
cnt_pb <- length(p_both$SpeciesCode)

p_mush <- subset(p, (invert == 2 | fish == 2))
cnt_mush <- length(p_mush$SpeciesCode)

ddply(p_empty, .(p_empty$SpeciesCode), function(x) {length(x$SpeciesCode)})

# total individual benthic invertivores:
cnt_p
# all prey combinations add up to total Pi individs :D
cnt_pe +
cnt_pi +
cnt_pf +
cnt_pb +
cnt_mush

# Benthic Invertivores:
b <- subset(myspp, fg == 'BI')
cnt_b <- length(b$SpeciesCode)

b_empty <- subset(b, (invert == 0 & fish == 0))
cnt_be <- length(b_empty$SpeciesCode)

b_invert <- subset(b, (invert == 1 & fish == 0))
cnt_bi <- length(b_invert$SpeciesCode)

b_fish <- subset(b, (invert == 0 & fish == 1))
cnt_bf <- length(b_fish$SpeciesCode)

b_both <- subset(b, (invert == 1 & fish ==1))
cnt_bb <- length(b_both$SpeciesCode)

b_mush <- subset(b, (invert == 2 | fish == 2))
cnt_mush <- length(b_mush$SpeciesCode)

# total individual benthic invertivores:
cnt_b 
# all prey combinations add up to total BI individs :D
cnt_be +
cnt_bi +
cnt_bf + 
cnt_bb +
cnt_mush

# how many fish had more than one prey item in their stomach?
p <- subset(prey3, fg == 'Pi')
length(which(duplicated(p$specimen)) == TRUE)
length(unique(p$specimen))

b <- subset(prey3, fg == 'BI')
length(which(duplicated(b$specimen)) == TRUE)
length(unique(b$specimen))

count_spp(p)
count_spp(b)





################################################################################
############        Overall Predator Counts for Thesis Table        ############
################################################################################


# Size ranges for each species in pento: ####
spp_list <- unique(pento$SpeciesCode)
m <- matrix(data = NA, nrow = 0, ncol = 4)

# populating matrix with min and max sizes for each pento spp
for (i in unique(pento$SpeciesCode)) {
  min <- min(fish[which(fish$SpeciesCode == i), 13])
  max <- max(fish[which(fish$SpeciesCode == i), 13])
  fg  <- unique(fish[which(fish$SpeciesCode == i), 4])
  m   <- rbind(m, c(as.character(i), as.character(fg), 
                    as.character(min), as.character(max)))
}

# getting this matrix setup for plotting
df <- as.data.frame(m)
colnames(df) <- c("SpeciesCode", "fg", "min", "max")
df$min <- as.numeric(as.character(df$min))
df$max <- as.numeric(as.character(df$max))
df$SpeciesCode <- factor(df$SpeciesCode, levels = pento_order)
df$SpeciesCode <- factor(df$SpeciesCode, levels = rev(levels(df$SpeciesCode)))
df$y <- as.numeric(seq(1:23)) # set y after order has been set
df <- merge(df, maxL)

spp_order <- df[with(df, order(SpeciesCode)),]
spp_order$y <- as.numeric(seq(1:23))


df$fg <- factor(df$fg, levels=c('Pi', 'BI', 'ZP', 'C', 'He'))

min_spp <- df[with(df, order(min)), ]

min_spp$y <- as.numeric(seq(1:23))

min <- ggplot(data=min_spp, aes(x=min, y=y, colour=fg)) +
  geom_segment(aes(xend=max, yend=y), lineend="round", size=0.8) +
  geom_segment(aes(x=max, xend=maxL, y=y, yend=y), linetype='dotted', 
               size=0.6, colour="black") +
  geom_point(aes(x=max, y=y), size=5, shape='|') +
  geom_point(aes(x=min, y=y), size=5, shape='|') +
  geom_point(aes(x=maxL, y=y), size=5, shape='|', colour="black") +
  xlim(-5, 1400) +
  xlab("Sampled body size range (mm)") +
  theme(axis.title.y=element_blank()) +
  theme(axis.ticks.y = element_blank(), axis.text.y=element_blank()) +
  geom_text(data=min_spp, aes(x=(maxL + 140), y=y, label=SpeciesCode)) 


spp <- ggplot(data = spp_order, aes(x=min, y=y, colour = fg)) +
  geom_segment(aes(xend=max, yend = y), lineend = "round", size=0.8) +
  geom_segment(aes(x=max, xend=maxL, y=y, yend=y), linetype='dotted', 
               size=0.6, colour="black") +
  geom_point(aes(x=max, y=y), size=5, shape='|') + 
  geom_point(aes(x=min, y=y), size=5, shape='|') +
  geom_point(aes(x=maxL, y=y), size=5, shape='|', colour="black") +
  xlim(-5, 1400) +
  xlab("Sampled body size range (mm)") +
  theme(axis.title.y=element_blank()) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_text(data=spp_order, aes(x=(maxL + 140), y=y, label=SpeciesCode))

g_legend <- function(a.gplot){
  tmp    <- ggplot_gtable(ggplot_build(a.gplot))
  leg    <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend<-g_legend(spp)

grid.arrange(spp + theme(legend.position = "none"), 
             min + theme(legend.position = "none"),
             legend, widths=c(5,5,0.8), nrow=1)


# Size ranges of all species included in qr for pred-prey size ####
# list of extra species
extra_spp <- c("AP.VIRE", "CA.ORTH", "EP.HEXA", "EP.MACU", "EP.SPIL", "EP.TAUV")

# populating matrix with min and max sizes for each of the 6 extra spp
e <- matrix(data = NA, nrow = 0, ncol = 4)
for (i in extra_spp) {
  min <- min(prey3[which(prey3$SpeciesCode == i), 5])
  max <- max(prey3[which(prey3$SpeciesCode == i), 5])
  fg  <- unique(prey3[which(prey3$SpeciesCode == i), 4])
  e   <- rbind(e, c(as.character(i), as.character(fg), 
                    as.character(min), as.character(max)))
}jjjj

# getting this matrix setup for plotting
e <- as.data.frame(e)
colnames(e) <- c("SpeciesCode", "fg", "min", "max")
e$min <- as.numeric(as.character(e$min))
e$max <- as.numeric(as.character(e$max))


pento_extra <- c("CA.MELA", "CA.ORTH", "AP.FURC", "AP.VIRE", "LU.BOHA", 
                 "LU.KASM", "CE.ARGU", "CE.UROD", "EP.HEXA", "EP.MACU",
                 "EP.SPIL", "EP.TAUV", "VA.LOUT", "PA.ARCA", "MO.GRAN",
                 "PA.INSU", "AC.NIGR", "AC.OLIV", "CE.FLAV", "CH.SORD",
                 "SC.FREN", "SC.RUBR", "CA.TERE", "PT.TILE", "CH.VAND", 
                 "PS.BART", "PS.DISP", "PS.OLIV", "CH.ORNA")


df$SpeciesCode <- factor(df$SpeciesCode, levels = pento_order)
df$y <- as.numeric(seq(1:23)) # set y after order has been set

all <- rbind(df, e)

spp_order <- df[with(df, order(min)),]
spp_order$y <- as.numeric(seq(1:23))

df$fg <- factor(df$fg, levels=c('Pi', 'BI', 'ZP', 'C', 'He'))

min_spp <- df[with(df, order(fg, min)), ]

min_spp$y <- as.numeric(seq(1:23))



ggplot(data = spp_order, aes(x=min, y=y, colour = fg)) +
  geom_segment(aes(xend=max, yend = y), lineend = "round", size=0.8) +
  geom_segement(data = maxL, aes(x=))
geom_point(aes(x=max, y=y), size=2, shape=19) + 
  geom_point(aes(x=min, y=y), size=2, shape=19) +
  xlim(-50, 720) +
  xlab("Size Range Sampled (mm)") +
  theme(axis.title.y=element_blank()) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_text(data=spp_order, aes(x=((min+max)/2), y=y+0.5, label=SpeciesCode))


ggplot(data = min_order, aes(x=min, y=y, colour = fg)) +
  geom_segment(aes(xend=max, yend = y), lineend = "round", size=0.8) +
  geom_segement(data = maxL, aes(x=))
geom_point(aes(x=max, y=y), size=2, shape=19) + 
  geom_point(aes(x=min, y=y), size=2, shape=19) +
  xlim(-50, 720) +
  xlab("Size Range Sampled (mm)") +
  theme(axis.title.y=element_blank()) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_text(data=spp_order, aes(x=((min+max)/2), y=y+0.5, label=SpeciesCode))


