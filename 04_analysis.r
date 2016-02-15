################################################################################
############            Gape Size ~ Body Size Analysis              ############
################################################################################

# Count number of fish per each species in df = fish
ddply(fish, .(SpeciesCode), summarize, n = length(SpeciesCode))

# Functional group level linear models:
fg_gh <- groupwise_lm_gh(pento, pento$j_fg)
fg_gw <- groupwise_lm_gw(pento, pento$j_fg)
fg_ga <- groupwise_lm_ga(pento, pento$j_fg)

fg_list <- rbind(fg_gh, fg_gw, fg_ga)

write.csv(fg_list, file = "fg_list.csv")

zp_gh <- groupwise_lm_gh(zp, zp$Family)
zp_gw <- groupwise_lm_gw(zp, zp$Family)
zp_ga <- groupwise_lm_ga(zp, zp$Family)

zp_list <- rbind(zp_gh, zp_gw, zp_ga)

write.csv(zp_list, file = "zp_fam_list.csv")

# Piscivore analysis:
# Linear mixed effects models for Functional groups:
library(lme4)
fm_p <- lmer(log(gh) ~ log(SL) + (1|SpeciesCode), p)
fm_p

fm_p2 <- lmer(log(gh) ~ log(SL) + (1|Family), p)
fm_p2

test <- lm(log(gh) ~ log(SL), data = p)
summary(test)

lm(log(gh) ~ log(SL), p)


x <- lme(log(gw) ~ log(SL), data = zp, random = ~ 1 | SpeciesCode)
y <- lme(log(gh) ~ log(SL), data = p, random = ~ 1 | Family)

summary(x)
intervals(x)
VarCorr(x)

# All Functional groups:
fg_gh <- groupwise_lm_gh(pento, pento$j_fg)
pspp_gw <- groupwise_lm_gw(p, p$SpeciesCode)
fg_ga <- groupwise_lm_ga(pento, pento$j_fg)

fg <- write_lme_groups(fg_gh, fg_gh$variable)
s2 <- write_lme_groups(pspp_gw, pspp_gw$variable)
fg3 <- write_lme_groups(fg_ga, fg_ga$variable)

fg_names <- c("j_fg", "eqn")
colnames(fg) <- fg_names
colnames(s2) <- snames
colnames(fg3) <- fg_names

pento2 <- pento
levels(pento2$j_fg)[levels(pento2$j_fg)=="Pi"] <- "Piscivore"
levels(pento2$j_fg)[levels(pento2$j_fg)=="BI"]   <- "Benthic Invertivore"
levels(pento2$j_fg)[levels(pento2$j_fg)=="ZP"] <- "Zooplanktivore"
levels(pento2$j_fg)[levels(pento2$j_fg)=="He"]   <- "Herbivore"
levels(pento2$j_fg)[levels(pento2$j_fg)=="C"]   <- "Corallivore"

fg2 <- fg3
levels(fg2$j_fg)[levels(fg2$j_fg)=="Pi"] <- "Piscivore"
levels(fg2$j_fg)[levels(fg2$j_fg)=="BI"]   <- "Benthic Invertivore"
levels(fg2$j_fg)[levels(fg2$j_fg)=="ZP"] <- "Zooplanktivore"
levels(fg2$j_fg)[levels(fg2$j_fg)=="He"]   <- "Herbivore"
levels(fg2$j_fg)[levels(fg2$j_fg)=="C"]   <- "Corallivore"

ggplot(data=pento2, aes(x=log(SL), y=log(ga))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  ylim(0, 11.5) +
  geom_text(data=fg2, aes(x=4.8, y=10.5, label=eqn), size=6, parse=T) +
  theme(strip.text.x = element_text(size=22)) +
  facet_wrap(~ j_fg, ncol = 3)


# Level: Functional group - Piscivores ####
  

p1 <-
  ggplot(data=p, aes(x=log(SL), y=log(gh))) +
  geom_point(shape=20) + #, aes(colour=SpeciesCode)) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.2, y=4.5, label=write_lme_gen(p, p$gh)), parse=TRUE) +
#  annotate("text", x = 4.2, y = 4.3, label = "-", size = 9)
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(vertical gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
  #ggtitle("Vertical Gape")

p2 <-
  ggplot(data=p, aes(x=log(SL), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.2, y=4.5, label=write_lme_gen(p, p$gw)), parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(horizontal gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
  #ggtitle("Horizontal Gape")

p3 <-
  ggplot(data=p, aes(x=log(SL), y=log(ga))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.2, y=8.2, label=write_lme_gen(p, p$ga)), parse=TRUE) +
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  theme(plot.margin = unit(c(0.5,1,1,1.3), "cm"))
  #ggtitle("Gape Area")

grid.arrange(p1, p2, p3)

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
  ggplot(data=p, aes(x=log(SL), y=log(gh))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f1, aes(x=4.8, y=5, label=eqn), size=4, parse=TRUE) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(vertical gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

fam2 <-
  ggplot(data=p, aes(x=log(SL), y=log(gw))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f2, aes(x=4.8, y=5.1, label=eqn), size=4, parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(horizontal gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

fam3 <-
  ggplot(data=p, aes(x=log(SL), y=log(ga))) + #, colour=SpeciesCode)) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f3, aes(x=4.8, y=9.9, label=eqn), size=4, parse=TRUE) +
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  theme(plot.margin = unit(c(0.5,1,1,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

grid.arrange(fam1, fam2)
fam3
grid.arrange(fam1, fam2, fam3)

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
  ggplot(data=p, aes(x=log(SL), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=4, drop=F) +
  geom_smooth(method=lm) +
  geom_text(data=s1, aes(x=4.85, y=4.8, label=eqn), size=3.5, parse=T) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(vertical gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
  facet_grid(. ~ cyl, labeller = label_both)

spp2 <-
  ggplot(data=p, aes(x=log(SL), y=log(gw))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=4) +
  geom_smooth(method=lm) +
  ylim(c(1,5.3)) +
  geom_text(data=s2, aes(x=4.85, y=5.1, label=eqn), size=3.5, parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(horizontal gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

spp3 <-
  ggplot(data=p, aes(x=log(SL), y=log(ga))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=4) +
  geom_smooth(method=lm) +
  ylim(c(1.9, 10.2)) +
  geom_text(data=s3, aes(x=4.8, y=9.5, label=eqn), size = 3.5, parse=TRUE) + 
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  theme(plot.margin = unit(c(0.5,1,1,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

spp3

grid.arrange(spp1, spp2, spp3)

pdf(file = "Gape_Body Plots/all_p_plots.pdf", height = 13, width = 9)
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
  ggplot(data=b, aes(x=log(SL), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.5, y=3.7, label=write_lme_gen(b, b$gh)), parse=TRUE) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(vertical gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

b2 <-
  ggplot(data=b, aes(x=log(SL), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.5, y=3.7, label=write_lme_gen(b, b$gw)), parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(horizontal gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

b3 <-
  ggplot(data=b, aes(x=log(SL), y=log(ga))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.5, y=7.5, label=write_lme_gen(b, b$ga)), parse=TRUE) +
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

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
  ggplot(data=b, aes(x=log(SL), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3, drop=F) +
  geom_smooth(method=lm) +
  geom_text(data=s1, aes(x=4.85, y=4.1, label=eqn), size=4, parse=T) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(vertical gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

spp2 <-
  ggplot(data=b, aes(x=log(SL), y=log(gw))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  #  ylim(c(1,5.3)) +
  geom_text(data=s2, aes(x=4.85, y=3.9, label=eqn), size=4, parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(horizontal gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

spp3 <-
  ggplot(data=b, aes(x=log(SL), y=log(ga))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  #  ylim(c(1.9, 10.2)) +
  geom_text(data=s3, aes(x=4.85, y=7.7, label=eqn), size=4, parse=TRUE) + 
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

grid.arrange(spp1, spp2, spp3)

pdf(file = "all_b_plots.pdf", height = 13, width = 9)
grid.arrange(b1, b2, b3)
grid.arrange(spp1, spp2, spp3)
dev.off()


# Zooplanktivore analysis: ####
# Level: Functional group - Zooplanktivores ####

zp1 <-
  ggplot(data=zp, aes(x=log(SL), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=3.5, y=3.2, label=write_lme_gen(zp, zp$gh)), parse=TRUE) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(vertical gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

zp2 <-
  ggplot(data=zp, aes(x=log(SL), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=3.5, y=3.2, label=write_lme_gen(zp, zp$gw)), parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(horizontal gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

zp3 <-
  ggplot(data=zp, aes(x=log(SL), y=log(ga))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=3.5, y=6, label=write_lme_gen(zp, zp$ga)), parse=TRUE) +
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

grid.arrange(zp1, zp2, zp3)

# Level: Family - Zooplanktivores ####

zpfam_gh <- groupwise_lm_gh(zp, zp$Family)
zpfam_gw <- groupwise_lm_gw(zp, zp$Family)
zpfam_ga <- groupwise_lm_ga(zp, zp$Family)

zplist <- rbind(zpfam_gh, zpfam_gw, zpfam_ga)
#write.csv(zplist, file = "zplist.csv")

f1 <- write_lme_groups(zpfam_gh, zpfam_gh$variable)
f2 <- write_lme_groups(zpfam_gw, zpfam_gw$variable)
f3 <- write_lme_groups(zpfam_ga, zpfam_ga$variable)

fnames <- c("Family", "eqn")
colnames(f1) <- fnames
colnames(f2) <- fnames
colnames(f3) <- fnames

fam1 <-
  ggplot(data=zp, aes(x=log(SL), y=log(gh))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f1, aes(x=4.2, y=3.4, label=eqn), size=4, parse=TRUE) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(vertical gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

fam1

fam2 <-
  ggplot(data=zp, aes(x=log(SL), y=log(gw))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f2, aes(x=4.2, y=3.4, label=eqn), size=4, parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(horizontal gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

fam2

fam3 <-
  ggplot(data=zp, aes(x=log(SL), y=log(ga))) + #, colour=SpeciesCode)) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f3, aes(x=4.2, y=6.4, label=eqn), size=4, parse=TRUE) +
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
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
  ggplot(data=zp, aes(x=log(SL), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3, drop=F) +
  geom_smooth(method=lm) +
  geom_text(data=s1, aes(x=4.2, y=3.3, label=eqn), size=4, parse=T) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(vertical gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

spp1

spp2 <-
  ggplot(data=zp, aes(x=log(SL), y=log(gw))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  geom_text(data=s2, aes(x=4.1, y=3.25, label=eqn), size=4, parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(horizontal gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

spp2

spp3 <-
  ggplot(data=zp, aes(x=log(SL), y=log(ga))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  geom_text(data=s3, aes(x=4.1, y=6.4, label=eqn), size=4, parse=TRUE) + 
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
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
  ggplot(data=h, aes(x=log(SL), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +  
  geom_text(data=data.frame(), aes(x=4.4, y=3.8, label=write_lme_gen(h, h$gh)), parse=TRUE) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(vertical gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

h2 <-
  ggplot(data=h, aes(x=log(SL), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.4, y=3.4, label=write_lme_gen(h, h$gw)), parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(horizontal gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

h3 <-
  ggplot(data=h, aes(x=log(SL), y=log(ga))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.4, y=6.7, label=write_lme_gen(h, h$ga)), parse=TRUE) +
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

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
  ggplot(data=h, aes(x=log(SL), y=log(gh))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f1, aes(x=5, y=4.2, label=eqn), size=4, parse=TRUE) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(vertical gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

fam2 <-
  ggplot(data=h, aes(x=log(SL), y=log(gw))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f2, aes(x=5, y=3.7, label=eqn), size=4, parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(horizontal gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

fam3 <-
  ggplot(data=h, aes(x=log(SL), y=log(ga))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=f3, aes(x=5, y=7.5, label=eqn), size=4, parse=TRUE) +
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
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
  ggplot(data=h, aes(x=log(SL), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  geom_text(data=s1, aes(x=4.9, y=4.15, label=eqn), size=4, parse=TRUE) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(vertical gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

spp2 <-
  ggplot(data=h, aes(x=log(SL), y=log(gw))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  #  ylim(c(1,5.3)) +
  geom_text(data=s2, aes(x=4.8, y=3.65, label=eqn), size=4, parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(horizontal gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
facet_grid(. ~ cyl, labeller = label_both)

spp3 <-
  ggplot(data=h, aes(x=log(SL), y=log(ga))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  #  ylim(c(1.9, 10.2)) +
  geom_text(data=s3, aes(x=4.8, y=7.45, label=eqn), size=4, parse=TRUE) + 
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
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
  ggplot(data=c, aes(x=log(SL), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  #geom_text(aes(position="jitter", label=SpecimenID, size=2)) +
  geom_text(data=data.frame(), aes(x=4.6, y=2.8, label=write_lme_gen(c, c$gh)), parse=TRUE) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(vertical gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

c2 <-
  ggplot(data=c, aes(x=log(SL), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.6, y=2.8, label=write_lme_gen(c, c$gw)), parse=TRUE) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  theme(axis.title.x = element_blank()) +
  ylab("log(horizontal gape, mm)") +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

c3 <-
  ggplot(data=c, aes(x=log(SL), y=log(ga))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.6, y=5.35, label=write_lme_gen(c, c$ga)), parse=TRUE) +
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.06, vjust= -1)) +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

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



x <- groupwise_lm_gh(p, p$Family)
x

y <- groupwise_lm_gw(p, p$Family)
y

z <- groupwise_lm_ga(p, p$Family)
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


gh <- ggplot(pento, aes(x=SpeciesCode, y=gh_ratio, fill=j_fg)) +
  geom_boxplot() +
  xlab("Species") +
  ylab("Relative vertical gape") +
  theme(axis.text.x = element_text(angle=45, colour="black", vjust=0.5)) +
  theme(axis.text.y = element_text(colour="black")) +
  theme(axis.title.y = element_text(vjust = 0.3)) +
  scale_fill_discrete(name = "  FG")

gw <- ggplot(pento, aes(x=SpeciesCode, y=gw_ratio, fill=j_fg)) +
  geom_boxplot() +
  xlab("Species") +
  ylab("Relative horizontal gape") +
  theme(axis.text.x = element_text(angle=45, colour="black", vjust=0.5)) +
  theme(axis.text.y = element_text(colour="black")) +
  theme(axis.title.y = element_text(vjust = 0.3)) +
  scale_fill_discrete(name = "  FG")

pdf(file = "pento_boxplots.pdf", width=9, height=6)
gh
gw
dev.off()

ggplot(pento, aes(x=SpeciesCode, y=gw_ratio, fill=j_fg)) +
  geom_boxplot() +
  xlab("Species") +
  ylab("Average Relative Horizontal Gape") +
  theme(axis.text.x = element_text(angle=45, colour="black", vjust=0.5)) +
  theme(axis.text.y = element_text(colour="black")) +
  theme(axis.title.y = element_text(vjust = 0.3))

shapes <- c(1, 4, 17, 8, 18, 7, 2)

gh <- ggplot(zp, aes(x=SL, y=gh_ratio, shape=SpeciesCode, colour=SpeciesCode)) +
  geom_point(size=2) +
#  scale_shape_manual(values=shapes) +
  ylab("Vertical gape (mm) / Standard length (mm)") +
  theme(axis.title.x = element_blank()) +
  xlim(0, 600) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.09, vjust= -1)) +
  theme(plot.margin = unit(c(0.5,1,0.1,1.3), "cm")) +
  geom_smooth(method=lm, fill=NA, size=1) 
gh

gw <- ggplot(zp, aes(x=SL, y=gw_ratio, shape=SpeciesCode, colour=SpeciesCode)) +
  geom_point(size=2) +
#  scale_shape_manual(values=shapes) +
  ylab("Horizontal gape (mm) / Standard length (mm)") +
  xlab("Standard length (mm)") +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.09, vjust= -1)) +
  theme(plot.margin = unit(c(0.1,1,0.5,1.3), "cm")) +
  geom_smooth(method=lm, fill=NA, size=1)
gw

pdf(file = "Ratio plots/zp_ratio_sl.pdf", width=8, height=9)
grid.arrange(gh, gw)
dev.off()



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

# Count number of fish per each species in df = preyX
ddply(prey5, .(SpeciesCode), summarize, n = length(SpeciesCode))
ddply(prey6, .(SpeciesCode), summarize, n = length(SpeciesCode))

# prey3 is a data.frame with just Pi, BI, ZP.
sl <- ggplot(prey2, aes(x = sl, y = psize)) +
  geom_point(aes(colour = ptype, shape = ptype)) +
  xlab("Predator standard length (mm)") +
  ylab("Prey minimum total length (mm)") +
  scale_colour_discrete(name = "Prey\nType") +
  scale_shape_discrete(name = "Prey\nType") +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.09, vjust= -1)) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm")) +
  #geom_smooth(method=lm, fill=NA, size=1)
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq")


gh <- ggplot(na.omit(prey2), aes(x = gh, y = psize)) +
  geom_point(aes(colour = ptype)) +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq")

grid.arrange(sl, gh)

df.n <- ddply(.data=prey3, .(fg), summarize, n=paste("n ==", length(fg)))

sl <- ggplot(prey3, aes(x = sl, y = psize)) +
  geom_point(aes(colour = ptype, shape = ptype)) +
  geom_text(data = df.n, aes( x= 200, y = 230, label = n), parse = TRUE) +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq") +
  ylab("Prey total length (mm)") +
  xlab("Predator standard length (mm)") +
  scale_colour_discrete(name = "Prey type") +
  scale_shape_discrete(name = "Prey type") +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.09, vjust= -1)) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm")) +
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
  geom_point(aes(colour = ptype, shape = ptype)) +
  geom_text(data = df.n, aes( x= 30, y = 230, label = n), parse = TRUE) +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq") +
  ylab("Prey total length (mm)") +
  xlab("Predator vertical gape (mm)") +
  scale_colour_discrete(name = "Prey type") +
  scale_shape_discrete(name = "Prey type") +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.09, vjust= -1)) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm")) +
  facet_wrap(~ fg)
gh

gw <- ggplot(na.omit(prey3), aes(x = gw, y = psize)) +
  geom_point(aes(colour = ptype, shape = ptype)) +
  geom_text(data = df.n, aes(x= 70, y = 220, label = n), parse = TRUE) +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq") +
  ylab("Prey total length (mm)") +
  xlab("Predator horizontal gape (mm)") +
  scale_colour_discrete(name = "Prey type") +
  scale_shape_discrete(name = "Prey type") +
  ggtitle("c)") +
  theme(plot.title = element_text(hjust= -0.09, vjust= -1)) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm")) +
  facet_wrap(~ fg, ncol=2)

p_prey <- subset(prey3, fg == "Pi")
df.n <- ddply(.data=na.omit(p_prey), .(fg), summarize, n=paste("n ==", length(fg)))
df.n <- ddply(.data=na.omit(prey3), .(fg), summarize, n=paste("n ==", length(fg)))

fgs <- list('Pi' = "Piscivore",
            'BI' = "Benthic Invertivore")

fg_labeller <- function(variable,value){
  return(fgs[value])
}

ga <- ggplot(na.omit(prey3), aes(x = pi*((gh/2)+ (gw/2)), y = psize)) +
  geom_point(aes(colour = ptype, shape = ptype), size = 4) +
  geom_text(data = df.n, aes(x= 125, y = 220, label = n), size = 10, parse = TRUE) +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq") +
  ylab("Prey total length (mm)") +
  xlab(expression(paste("Predator gape area ", "(", mm^2, ")", sep= ""))) +
  scale_colour_discrete(name = "Prey type") +
  scale_shape_discrete(name = "Prey type") +
  theme(axis.title.y = element_text(size = 30, vjust = 0.07)) +
  theme(axis.title.x = element_text(size = 30, vjust = 0.05)) +
  theme(legend.text = element_text(size = 28)) +
  theme(legend.key.height = unit(1.5, "line")) +
  theme(legend.title = element_text(size = 26)) +
  theme(strip.text.x = element_text(size = 34)) +
  theme(legend.position = c(0.9, 0.5)) +
  theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
  theme(axis.text.x = element_text(size = 28)) +
  theme(axis.text.y = element_text(size = 28)) +
#  ggtitle("b)") +
#  theme(plot.title = element_text(hjust= -0.09, vjust= -1)) +
  theme(plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))  +
  facet_grid(. ~ fg, labeller=fg_labeller)

pdf(file = "Pred_prey_plots/psize_ga_sl_fg.pdf", width = 8.5, height = 8.5)
grid.arrange(sl, ga) #, gw)
dev.off()

grid.arrange(gh, gw, ga)

df.n <- ddply(.data=prey5, .(fg), summarize, n=paste("n ==", length(fg)))

ggplot(na.omit(prey3), aes(x = sl, y = psize)) +
  geom_point(aes(colour = SpeciesCode)) +
  geom_text(data = df.n, aes( x= 200, y = 200, label = n), parse = TRUE) +
  #geom_smooth(method = "lm") +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.90), method = "rq") +
  facet_wrap(~ fg)

# Graphs for manuscript:
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
df.n <- ddply(.data=prey5, .(fg), summarize, n=paste("n ==", length(fg)))
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


master_layout <- 
grid.layout(nrow = 2, ncol = 3, 
      widths = unit(c(0.1, 1, 1), "null"),
      heights = unit(c(1, 0.1), "null"))
grid.newpage()
pushViewport(viewport(layout = master_layout))
print(pisc_prey, vp = set_vp(1, 2))
print(benth_prey, vp = set_vp(1, 3))
grid.text(
  expression( paste("Gape area (", mm^2, ")", sep = "") ), 
  vp = viewport(layout.pos.row = 2, layout.pos.col = 2:3), 
  gp = gpar(fontsize = 10), vjust = -0.25
  )
grid.text(
  "Prey total length (mm)",
  vp = viewport(layout.pos.row = 1, layout.pos.col = 1), 
  gp = gpar(fontsize = 10), rot = 90, vjust = 2
  )

grid.text(
  "a)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1), 
  gp = gpar(fontsize = 9), vjust = -7
  )
grid.text(
  "b)", vp = viewport(layout.pos.row = 4, layout.pos.col = 1), 
  gp = gpar(fontsize = 9), vjust = -7
  )

dev.copy2eps(device = quartz, file = "panel_plots/pred_prey_size.eps")




bi_qr <- rq(data = subset(prey3, prey3$fg == 'BI'), psize ~ sl, 
            tau = c(0.10, 0.50, 0.90))
bi_90 <- rq(data = subset(prey3, prey3$fg =='BI'), psize~gh, tau = 0.90)
bi_50 <- rq(data = subset(prey3, prey3$fg =='BI'), psize~gh, tau = 0.50)
bi_10 <- rq(data = subset(prey3, prey3$fg =='BI'), psize~gh, tau = 0.10)

p_90 <- rq(data = subset(prey3, prey3$fg =='Pi'), psize~sl, tau = 0.90)
p_50 <- rq(data = subset(prey3, prey3$fg =='Pi'), psize~sl, tau = 0.50)
p_10 <- rq(data = subset(prey3, prey3$fg =='Pi'), psize~sl, tau = 0.10)


p <- subset(prey3, prey3$fg == 'Pi')
b <- subset(prey3, prey3$fg == 'BI')

p <- subset(prey_gh, prey_gh$fg == 'Pi')
b <- subset(prey_gh, prey_gh$fg == 'BI')

summary(bi_qr, se = "rank")
bi_qr[[3]]

# QR for species with most pred-prey sizes to see if one species in particular 
# is responsible for eating certain prey items or prey sizes:
pento_extra <- c("CA.MELA", "CA.ORTH", "AP.FURC", "AP.VIRE", "LU.BOHA", 
                 "LU.KASM", "CE.ARGU", "CE.UROD", "EP.HEXA", "EP.MACU",
                 "EP.SPIL", "EP.TAUV", "VA.LOUT", "PA.ARCA", "MO.GRAN",
                 "PA.INSU", "AC.NIGR", "AC.OLIV", "CE.FLAV", "CH.SORD",
                 "SC.FREN", "SC.RUBR", "CA.TERE", "PT.TILE", "CH.VAND", 
                 "PS.BART", "PS.DISP", "PS.OLIV", "CH.ORNA")

prey3$SpeciesCode <- factor(prey3$SpeciesCode, levels=pento_extra)
six <- prey3[prey3$SpeciesCode %in% c('AP.FURC', 'CA.MELA', 'CE.ARGU', 
                                      'CE.UROD', 'LU.BOHA', 'VA.LOUT'), ]
six$SpeciesCode <- factor(six$SpeciesCode, levels = 'CA.MELA',
                                                    'AP.FURC',
                                                    'LU.BOHA',
                                                    'CE.ARGU',
                                                    'CE.UROD',
                                                    'VA.LOUT')

df.n <- ddply(.data=six, .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))
sl <- ggplot(six, aes(x = sl, y = psize)) +
  geom_point(aes(colour = ptype), name = "Prey Type") +
  geom_text(data = df.n, aes( x= 120, y = 230, label = n), parse = TRUE) +
  #stat_quantile(geom = "quantile", quantiles = c(0.5), method = "rq") +
  ylab("Prey total length (mm)") +
  theme(axis.title.y = element_text(vjust = 0.3)) +
  xlab("Predator standard length (mm)") +
  theme(axis.title.x = element_text(vjust = 0.1)) +
  labs(title = "a)") +
  theme(plot.title = element_text(hjust= -0.05, vjust= -1)) +
  labs(colour = "Prey Type") +
  facet_wrap(~ SpeciesCode) 
sl

df.n <- ddply(.data=na.omit(six), .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))

gh <- ggplot(na.omit(six), aes(x = gh, y = psize)) +
  geom_point(aes(colour = ptype)) +
  geom_text(data = df.n, aes( x= 30, y = 230, label = n), parse = TRUE) +
#  stat_quantile(geom = "quantile", quantiles = c(0.5), method = "rq") +
  ylab("Prey total length (mm)") +
  theme(axis.title.y = element_text(vjust = 0.3)) +
  xlab("Predator vertical gape (mm)") +
  theme(axis.title.x = element_text(vjust = 0.1)) +
  labs(title = "b)") +
  theme(plot.title = element_text(hjust= -0.05, vjust= -1)) +
  labs(colour = "Prey Type") +
  facet_wrap(~ SpeciesCode,)
gh

gw <- ggplot(na.omit(six), aes(x = gw, y = psize)) +
  geom_point(aes(colour = ptype)) +
  geom_text(data = df.n, aes( x= 30, y = 230, label = n), parse = TRUE) +
#  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq") +
  ylab("Prey total length (mm)") +
  theme(axis.title.y = element_text(vjust = 0.3)) +
  xlab("Predator horizontal gape (mm)") +
  theme(axis.title.x = element_text(vjust = 0.1)) +
  labs(title = "c)") +
  theme(plot.title = element_text(hjust= -0.05, vjust= -1)) +
  labs(colour = "Prey Type") +
  facet_wrap(~ SpeciesCode)
gw

ga <- ggplot(na.omit(six), aes(x = pi*(gh/2 + gw/2), y = psize)) +
  geom_point(aes(colour = ptype)) +
  geom_text(data = df.n, aes( x= 120, y = 230, label = n), parse = TRUE) +
  #  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq") +
  ylab("Prey total length (mm)") +
  theme(axis.title.y = element_text(vjust = 0.3)) +
  xlab(expression(paste("Predator gape area ", "(", mm^2, ")", sep= ""))) +
  theme(axis.title.x = element_text(vjust = 0.1)) +
  labs(title = "b)") +
  theme(plot.title = element_text(hjust= -0.05, vjust= -1)) +
  labs(colour = "Prey Type") +
  facet_wrap(~ SpeciesCode)
ga

pdf(file = "Pred_prey_plots/psize_six_ga_sl_fg.pdf", width = 8.5, height = 10)
grid.arrange(sl, ga) #, gw)
dev.off()

df.n <- ddply(.data=prey3, .(fg), summarize, n=paste("n ==", length(fg)))

ggplot(na.omit(prey3), aes(x = sl, y = psize)) +
  geom_point(aes(colour = species)) +
  geom_text(data = df.n, aes( x= 200, y = 200, label = n), parse = TRUE) +
  geom_smooth(method = "lm") +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.90), method = "rq") +
  facet_wrap(~ fg)


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

ddply(p_invert, .(p_invert$SpeciesCode), function(x) {length(x$SpeciesCode)})

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

# Displaying these counts as a graph: ####
# First, setting up dataframe so that presence/absence of prey types are a
# single set of factors:

myspp$stomach <- with(myspp,
                    ifelse((invert == 0 & fish == 0), 'E',
                          ifelse((invert == 0 & fish == 1), 'F',
                                 ifelse((invert == 1 & fish == 0), 'I',
                                        ifelse((invert == 1 & fish == 1), 'F + I', 'U')
                                        )
                                 )
                          )
                      )

myspp$stomach <- as.factor(myspp$stomach)
myspp$stomach <- factor(myspp$stomach, levels = c('E', 'F', 'I', 'F + I', 'U'))

p <- subset(myspp, fg == 'Pi')
b <- subset(myspp, fg == 'BI')

df.n <- ddply(.data=p, .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))

ggplot(data=p, aes(stomach, fill=SpeciesCode)) +
  geom_bar() +
  geom_text(data = df.n, aes(x = 3, y = 90, label = n), parse = TRUE) +
  facet_wrap(~ SpeciesCode)

# for plotting the most commmonly sampled species: CE.ARGU, CE.UROD, LU.BOHA: 
trio <- p[p$SpeciesCode %in% c("CE.ARGU", "CE.UROD", "LU.BOHA"), ]

df.n <- ddply(.data=trio, .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))

trio_plot <- ggplot(data=trio, aes(stomach, fill=SpeciesCode)) +
  geom_bar() +
  geom_text(data = df.n, aes(x = 4.5, y = 90, label = n), parse = TRUE) +
  xlab("Stomach contents") +
  ylab("Count") +
  theme(axis.title.x = element_text(vjust = 0.1)) +
#  theme(axis.text.x = element_text(angle = 30, vjust = 0.7)) +
  facet_wrap(~ SpeciesCode)

pdf(file = "Pred_prey_plots/trio_gut_contents_barchart.pdf", width = 9, height = 4)
trio_plot
dev.off()




################################################################################
############        Overall Predator Counts for Thesis Table        ############
################################################################################


# Size ranges for each species in pento: ####
pento$j_fg <- factor(pento$j_fg, levels=c('Pi', 'BI', 'He', 'ZP', 'C'))
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
#df$SpeciesCode <- factor(df$SpeciesCode, levels = pento_order)
#df$SpeciesCode <- factor(df$SpeciesCode, levels = rev(levels(df$SpeciesCode)))
#df$y <- as.numeric(seq(1:23)) # set y after order has been set
df <- merge(df, maxL)

# Size ranges of all species included in qr for pred-prey size ####
# list of EXTRA species
extra_spp <- c("AP.VIRE", "CA.ORTH", "EP.HEXA", "EP.MACU", "EP.SPIL", "EP.TAUV")
max_len_e <- c(750, 1120, 275, 605, 350, 750)
len_e <- rep("TL", 6)

maxLe <- data.frame(SpeciesCode=extra_spp, maxL=max_len_e, length=len_e)
# populating matrix with min and max sizes for each of the 6 extra spp
e <- matrix(data = NA, nrow = 0, ncol = 4)
for (i in extra_spp) {
  min <- min(prey3[which(prey3$SpeciesCode == i), 5])
  max <- max(prey3[which(prey3$SpeciesCode == i), 5])
  fg  <- unique(prey3[which(prey3$SpeciesCode == i), 4])
  e   <- rbind(e, c(as.character(i), as.character(fg), 
                    as.character(min), as.character(max)))
}

# getting this matrix setup for plotting
e <- as.data.frame(e)
colnames(e) <- c("SpeciesCode", "fg", "min", "max")
e$min <- as.numeric(as.character(e$min))
e$max <- as.numeric(as.character(e$max))
e$SpeciesCode <- as.character(e$SpeciesCode)
e <- merge(e, maxLe)

# Combining original species and extra (pred-prey) species for the 
# body size range plots. 

pento_extra <- c("CA.MELA", "CA.ORTH", "AP.FURC", "AP.VIRE", "LU.BOHA", 
                 "LU.KASM", "CE.ARGU", "CE.UROD", "EP.HEXA", "EP.MACU",
                 "EP.SPIL", "EP.TAUV", "VA.LOUT", "PA.ARCA", "MO.GRAN",
                 "PA.INSU", "AC.NIGR", "AC.OLIV", "CE.FLAV", "CH.SORD",
                 "SC.FREN", "SC.RUBR", "CA.TERE", "PT.TILE", "CH.VAND", 
                 "PS.BART", "PS.DISP", "PS.OLIV", "CH.ORNA")

# Combined original + extra species data frame: 
df <- rbind(df, e)
df$SpeciesCode <- factor(df$SpeciesCode, levels = pento_extra)
df$SpeciesCode <- factor(df$SpeciesCode, levels = rev(levels(df$SpeciesCode)))

spp_order <- df[with(df, order(SpeciesCode)),]
spp_order$y <- as.numeric(seq(1:29))
spp_order$fg <- factor(spp_order$fg, levels=c('Pi', 'BI', 'He', 'ZP', 'C'))

# max size of fish we measured are this % of max body size sampled:
ddply(spp_order, .(SpeciesCode), function(x) {(x$max/x$maxL)})

# custom colours with piscivores blue, herbivores green:
#"#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"
cc <- c("#00B0F6", "#F8766D", "#00BF7D", "#E76BF3", "#A3A500")

#     red        yellow     green     blue     purple
#[1] "#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"

min_spp <- df[with(df, order(min)), ]
min_spp$y <- as.numeric(seq(1:29))
min_spp$fg <- factor(min_spp$fg, levels=c('Pi', 'BI', 'He', 'ZP', 'C'))

max_spp <- df[with(df, order(max)), ]
max_spp$y <- as.numeric(seq(1:22))


min <- ggplot(data=min_spp, aes(x=min, y=y, colour=fg)) +
  geom_segment(aes(xend=max, yend=y), lineend="round", size=0.8) +
  geom_segment(aes(x=max, xend=maxL, y=y, yend=y), linetype='dotted', 
               size=0.6, colour="black") +
  geom_point(aes(x=max, y=y), size=5, shape='|') +
  geom_point(aes(x=min, y=y), size=5, shape='|') +
  geom_point(aes(x=maxL, y=y), size=5, shape='|', colour="black") +
  geom_point(data = data.frame(), aes(x=(255 + c(750, 1120, 750, 605, 330, 275)), 
                                      y=c(29.2, 28.2, 27.2, 26.2, 24.2, 21.2)), 
             shape = '*', size=6) +
  scale_colour_manual(values=cc) +
  xlim(-5, 1400) +
  xlab("Sampled body size range (mm)") +
  theme(axis.title.y=element_blank()) +
  theme(axis.ticks.y = element_blank(), axis.text.y=element_blank()) +
  geom_text(data=min_spp, aes(x=(maxL + 140), y=y, label=SpeciesCode)) +
  theme(plot.margin= unit(c(1, 0.3, 0.5, 0.4), "lines")) +
  labs(title = "b)") +
  theme(plot.title = element_text(hjust= -0.04, vjust= -1)) 
#min

spp <- ggplot(data = spp_order, aes(x=min, y=y, colour = fg)) +
  geom_segment(aes(xend=max, yend = y), lineend = "round", size=0.8) +
  geom_segment(aes(x=max, xend=maxL, y=y, yend=y), linetype='dotted', 
               size=0.6, colour="black") +
  geom_point(aes(x=max, y=y), size=5, shape='|') + 
  geom_point(aes(x=min, y=y), size=5, shape='|') +
  geom_point(aes(x=maxL, y=y), size=5, shape='|', colour="black") +
  geom_point(data = data.frame(), aes(x=(255 + c(1120, 750, 275, 605, 330, 750)), 
                                      y=c(28.2, 26.2, 21.2, 20.2, 19.2, 18.2)), 
             shape = '*', size=6) +
  scale_colour_manual(name="  FG",
                      values=cc) +
  xlim(-5, 1400) +
  xlab("Sampled body size range (mm)") +
  theme(axis.title.y=element_blank()) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_text(data=spp_order, aes(x=(maxL + 140), y=y, label=SpeciesCode)) +
  theme(plot.margin= unit(c(1, 0.4, 0.5, 0.5), "lines")) +
  labs(title = "a)") +
  theme(plot.title = element_text(hjust= -0.04, vjust= -1)) 
spp

legend<-extract_legend(spp)

grid.arrange(spp + theme(legend.position = "none"), 
             min + theme(legend.position = "none"),
             legend, widths=c(5,5,0.8), nrow=1)


# Size ranges of all species included in qr for pred-prey size ####
# list of extra species
extra_spp <- c("AP.VIRE", "CA.ORTH", "EP.HEXA", "EP.MACU", "EP.SPIL", "EP.TAUV")
max_len_e <- c(750, 1120, 275, 605, 350, 750)
len_e <- rep("TL", 6)

maxLe <- data.frame(SpeciesCode=extra_spp, maxL=max_len_e, length=len_e)
# populating matrix with min and max sizes for each of the 6 extra spp
e <- matrix(data = NA, nrow = 0, ncol = 4)
for (i in extra_spp) {
  min <- min(prey3[which(prey3$SpeciesCode == i), 5])
  max <- max(prey3[which(prey3$SpeciesCode == i), 5])
  fg  <- unique(prey3[which(prey3$SpeciesCode == i), 4])
  e   <- rbind(e, c(as.character(i), as.character(fg), 
                    as.character(min), as.character(max)))
}

# getting this matrix setup for plotting
e <- as.data.frame(e)
colnames(e) <- c("SpeciesCode", "fg", "min", "max")
e$min <- as.numeric(as.character(e$min))
e$max <- as.numeric(as.character(e$max))
e <- merge(e, maxLe)

# Combining original species and extra (pred-prey) species for the 
# body size range plots. 

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
  geom_segement(data = maxL, aes(x=)) +
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

fg_gh <- groupwise_lm_ga(pento, pento$j_fg)
fg_gh$num <- c(1,2,3,4,5)
# 
ggplot(data = fg_gh, aes(x=lw_conf_slp, y=num, colour = variable)) +
  geom_segment(aes(xend=up_conf_slp, yend=num),  size=0.8, colour="black") +
  geom_point(aes(x=slope, y=num), size=5, shape=19) +
  xlim(0.9,3.1) +
  geom_vline(xintercept=2) +
  xlab("Confidence Interval") +
  theme(axis.title.y=element_blank()) +
  theme(axis.ticks.y=element_blank()) +
  theme(axis.text.y=element_blank()) #+
  theme(panel.grid.minor=element_blank()) +
  theme(panel.grid.major=element_blank())

  
  geom_segment(aes(xend=max, yend = num), lineend = "round", size=0.8) +
  geom_segment(data = maxL, aes(x=)) +
  geom_point(aes(x=max, y=y), size=2, shape=19) + 
  geom_point(aes(x=min, y=y), size=2, shape=19) +
  xlim(0, 2) +
  xlab("Size Range Sampled (mm)") +
  theme(axis.title.y=element_blank()) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_text(data=spp_order, aes(x=((min+max)/2), y=y+0.5, label=SpeciesCode))
















