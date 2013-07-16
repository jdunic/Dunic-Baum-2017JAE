################################################################################
########                Individual functional group SMAs                ########
################################################################################
A2 <- par(mfrow = c(1, 1))
B2 <- par(mfrow = c(2, 1))
C2 <- par(mfrow = c(3, 2))

A3 <- par(mfrow = c(1, 1))
B3 <- par(mfrow = c(2, 1))
C3 <- par(mfrow = c(3, 2))

xh

# Piscivores
pGH <- sma(gh~SL, data=p, log="xy", method="SMA", robust=T, slope.test=1)
pGW <- sma(gw~SL, data=p, log="xy", method="SMA", robust=T, slope.test=1)
pGA <- sma(ga~SL, data=p, log="xy", method="SMA", robust=T, slope.test=2)

# verifying SMA assumptions
check_assump(pGH)
check_assump(pGW)
check_assump(pGA)

# Benthic Invertivores
bGH <- sma(gh~SL, data=b, log="xy", method="SMA", robust=T, slope.test=1)
bGW <- sma(gw~SL, data=b, log="xy", method="SMA", robust=T, slope.test=1)
bGA <- sma(ga~SL, data=b, log="xy", method="SMA", robust=T, slope.test=2)

# verifying SMA assumptions
check_assump(bGH)
check_assump(bGW)
check_assump(bGA)

# Zooplanktivores
zGH <- sma(gh~SL, data=zp, log="xy", method="SMA", robust=T, slope.test=1)
zGW <- sma(gw~SL, data=zp, log="xy", method="SMA", robust=T, slope.test=1)
zGA <- sma(ga~SL, data=zp, log="xy", method="SMA", robust=T, slope.test=2)
zGH

# verifying SMA assumptions
check_assump(zGH)
check_assump(zGW)
check_assump(zGA)

# Herbivores
hGH <- sma(gh~SL, data=h, log="xy", method="SMA", robust=T, slope.test=1)
hGW <- sma(gw~SL, data=h, log="xy", method="SMA", robust=T, slope.test=1)
hGA <- sma(ga~SL, data=h, log="xy", method="SMA", robust=T, slope.test=2)
hGH

# verifying SMA assumptions
check_assump(hGH)
check_assump(hGW)
check_assump(hGA)

# Corallivores
cGH <- sma(gh~SL, data=c, log="xy", method="SMA", robust=T, slope.test=1)
cGW <- sma(gw~SL, data=c, log="xy", method="SMA", robust=T, slope.test=1)
cGA <- sma(ga~SL, data=c, log="xy", method="SMA", robust=T, slope.test=2)
cGH

plot(cGH)

# verifying SMA assumptions
check_assump(cGH)
check_assump(cGW)
check_assump(cGA)


################################################################################
########                 Functional group SMA summaries                 ########
################################################################################
u <- mk_sma_summary(pGH, "pGH")
v <- mk_sma_summary(bGH, "bGH")
w <- mk_sma_summary(zGH, "zGH")
x <- mk_sma_summary(hGH, "hGH")
y <- mk_sma_summary(cGH, "cGH")
fg_height <- cbind(u, v, w, x, y)
fg_height

u <- mk_sma_summary(pGW, "pGW")
v <- mk_sma_summary(bGW, "bGW")
w <- mk_sma_summary(zGW, "zGW")
x <- mk_sma_summary(hGW, "hGW")
y <- mk_sma_summary(cGW, "cGW")
fg_width <- cbind(u, v, w, x, y)
fg_width

u <- mk_sma_summary(pGA, "pGA")
v <- mk_sma_summary(bGA, "bGA")
w <- mk_sma_summary(zGA, "zGA")
x <- mk_sma_summary(hGA, "hGA")
y <- mk_sma_summary(cGA, "cGA")
fg_area <- cbind(u, v, w, x, y)
fg_area

write.csv(fg_height, "fg_height.csv")
write.csv(fg_width, "fg_width.csv")
write.csv(fg_area, "fg_area.csv")

################################################################################
########                 Family level SMA analyses                      ########
################################################################################

# Vertical Gape
#===============================================================================
# Piscivores families:
carangidae_gh <- sma(gh~SL, data=p[which(p$Family=="Carangidae"), ], log="xy", 
                  method="SMA", robust=F, slope.test=1)
lutjanidae_gh <- sma(gh~SL, data=p[which(p$Family=="Lutjanidae"), ], log="xy",
                  method="SMA", robust=F, slope.test=1)
p_serranidae_gh <- sma(gh~SL, data=p[which(p$Family=="Serranidae"), ], log="xy",
                  method="SMA", robust=F, slope.test=1)

check_assump(carangidae_gh)
check_assump(lutjanidae_gh)
check_assump(p_serranidae_gh)

# Benthic invertivore families/species:
monotaxis_gh <- sma(gh~SL, data=b[which(b$SpeciesCode=="MO.GRAN"), ], log="xy",
                 method="SMA", robust=F, slope.test=1)
paracirrhites_gh <- sma(gh~SL, data=b[which(b$SpeciesCode=="PA.ARCA"), ], 
                     log="xy", method="SMA", robust=F, slope.test=1)
parupeneus_gh <- sma(gh~SL, data=b[which(b$SpeciesCode=="PA.INSU"), ], log="xy",
                  method="SMA", robust=F, slope.test=1)

check_assump(monotaxis_gh)
check_assump(paracirrhites_gh)
check_assump(parupeneus_gh)

# Zooplanktivore families:
caesionidae_gh <- sma(gh~SL, data=z[which(z$Family=="Caesionidae"), ], 
                   log="xy", method="SMA", robust=F, slope.test=1)
pomacentridae_gh <- sma(gh~SL, data=z[which(z$Family=="Pomacentridae"), ], 
                     log="xy", method="SMA", robust=F, slope.test=1)
z_serranidae_gh <- sma(gh~SL, data=z[which(z$Family=="Serranidae"), ],
                  log="xy", method="SMA", robust=F, slope.test=1)

check_assump(caesionidae_gh)
check_assump(pomacanthidae_gh)
check_assump(z_serranidae_gh)

# Herbivore families:
acanthuridae_gh <- sma(gh~SL, data=h[which(h$Family=="Acanthuridae"), ], 
                    log="xy", method="SMA", robust=F, slope.test=1)
pomacanthidae_gh <- sma(gh~SL, data=h[which(h$Family=="Pomacanthidae"), ],
                     log="xy", method="SMA", robust=F, slope.test=1)
scaridae_gh <- sma(gh~SL, data=h[which(h$Family=="Scaridae"), ], 
                log="xy", method="SMA", robust=F, slope.test=1)

check_assump(acanthuridae_gh)
check_assump(pomacanthidae_gh)
check_assump(scaridae_gh)

# Corallivore:
chaetodon_gh <- sma(gh~SL, data=c[which(c$SpeciesCode=="CH.ORNA"), ], 
                   log="xy", method="SMA", robust=F, slope.test=1)

u <- mk_sma_summary(carangidae_gh, "carangidae_gh")
v <- mk_sma_summary(lutjanidae_gh, "lutjanidae_gh")
w <- mk_sma_summary(p_serranidae_gh, "serranidae_gh")
p_height <- cbind(u, v, w)

u <- mk_sma_summary(monotaxis_gh, "MO.GRAN_gh")
v <- mk_sma_summary(paracirrhites_gh, "PA.ARCA_gh")
w <- mk_sma_summary(parupeneus_gh, "PA.INSU_gh")
b_height <- cbind(u, v, w)

u <- mk_sma_summary(caesionidae_gh, "caesionidae_gh")
v <- mk_sma_summary(pomacentridae_gh, "pomacentridae_gh")
w <- mk_sma_summary(z_serranidae_gh, "serranidae_gh")
z_height <- cbind(u, v, w)

u <- mk_sma_summary(acanthuridae_gh, "acanthuridae_gh")
v <- mk_sma_summary(pomacanthidae_gh, "pomacanthidae_gh")
w <- mk_sma_summary(scaridae_gh, "scaridae_gh")
h_height <- cbind(u, v, w)

c_height <- mk_sma_summary(chaetodon_gh, "chaetodon_gh")

write.csv(p_height, "p_family_gh.csv")
write.csv(b_height, "b_family_gh.csv")
write.csv(z_height, "z_family_gh.csv")
write.csv(h_height, "h_family_gh.csv")
write.csv(c_height, "c_family_gh.csv")


# Horizontal Gape
#===============================================================================

# Piscivores families:
carangidae_gw <- sma(gw~SL, data=p[which(p$Family=="Carangidae"), ], log="xy", 
                     method="SMA", robust=F, slope.test=1)
lutjanidae_gw <- sma(gw~SL, data=p[which(p$Family=="Lutjanidae"), ], log="xy",
                     method="SMA", robust=F, slope.test=1)
p_serranidae_gw <- sma(gw~SL, data=p[which(p$Family=="Serranidae"), ], log="xy",
                     method="SMA", robust=F, slope.test=1)

check_assump(carangidae_gw)
check_assump(lutjanidae_gw)
check_assump(p_serranidae_gw)

# Benthic invertivore families/species:
monotaxis_gw <- sma(gw~SL, data=b[which(b$SpeciesCode=="MO.GRAN"), ], log="xy",
                    method="SMA", robust=F, slope.test=1)
paracirrhites_gw <- sma(gw~SL, data=b[which(b$SpeciesCode=="PA.ARCA"), ], 
                        log="xy", method="SMA", robust=F, slope.test=1)
parupeneus_gw <- sma(gw~SL, data=b[which(b$SpeciesCode=="PA.INSU"), ], log="xy",
                     method="SMA", robust=F, slope.test=1)
check_assump(monotaxis_gw)
check_assump(paracirrhites_gw)
check_assump(parupeneus_gw)

# Zooplanktivore families:
caesionidae_gw <- sma(gw~SL, data=z[which(z$Family=="Caesionidae"), ], 
                      log="xy", method="SMA", robust=F, slope.test=1)
pomacentridae_gw <- sma(gw~SL, data=z[which(z$Family=="Pomacentridae"), ], 
                        log="xy", method="SMA", robust=F, slope.test=1)
z_serranidae_gw <- sma(gw~SL, data=z[which(z$Family=="Serranidae"), ],
                     log="xy", method="SMA", robust=F, slope.test=1)

check_assump(caesionidae_gw)
check_assump(pomacentridae_gw)
check_assump(z_serranidae_gw)

# Herbivore families:
acanthuridae_gw <- sma(gw~SL, data=h[which(h$Family=="Acanthuridae"), ], 
                       log="xy", method="SMA", robust=F, slope.test=1)
pomacanthidae_gw <- sma(gw~SL, data=h[which(h$Family=="Pomacanthidae"), ],
                        log="xy", method="SMA", robust=F, slope.test=1)
scaridae_gw <- sma(gw~SL, data=h[which(h$Family=="Scaridae"), ], 
                   log="xy", method="SMA", robust=F, slope.test=1)

check_assump(acanthuridae_gw)
check_assump(pomacanthidae_gw)
check_assump(scaridae_gw)

# Corallivore:
chaetodon_gw <- sma(gw~SL, data=c[which(c$SpeciesCode=="CH.ORNA"), ], 
                    log="xy", method="SMA", robust=F, slope.test=1)

u <- mk_sma_summary(carangidae_gw, "carangidae_gw")
v <- mk_sma_summary(lutjanidae_gw, "lutjanidae_gw")
w <- mk_sma_summary(p_serranidae_gw, "serranidae_gw")
p_width <- cbind(u, v, w)

u <- mk_sma_summary(monotaxis_gw, "MO.GRAN_gw")
v <- mk_sma_summary(paracirrhites_gw, "PA.ARCA_gw")
w <- mk_sma_summary(parupeneus_gw, "PA.INSU_gw")
b_width <- cbind(u, v, w)

u <- mk_sma_summary(caesionidae_gw, "caesionidae_gw")
v <- mk_sma_summary(pomacentridae_gw, "pomacentridae_gw")
w <- mk_sma_summary(z_serranidae_gw, "serranidae_gw")
z_width <- cbind(u, v, w)

u <- mk_sma_summary(acanthuridae_gw, "acanthuridae_gw")
v <- mk_sma_summary(pomacanthidae_gw, "pomacanthidae_gw")
w <- mk_sma_summary(scaridae_gw, "scaridae_gw")
h_width <- cbind(u, v, w)

c_width <- mk_sma_summary(chaetodon_gw, "chaetodon_gw")

write.csv(p_width, "p_family_gw.csv")
write.csv(b_width, "b_family_gw.csv")
write.csv(z_width, "z_family_gw.csv")
write.csv(h_width, "h_family_gw.csv")
write.csv(c_width, "c_family_gw.csv")


# Gape Area
#===============================================================================

# Piscivores families:
carangidae_ga <- sma(ga~SL, data=p[which(p$Family=="Carangidae"), ], log="xy", 
                     method="SMA", robust=F, slope.test=2)
lutjanidae_ga <- sma(ga~SL, data=p[which(p$Family=="Lutjanidae"), ], log="xy",
                     method="SMA", robust=F, slope.test=2)
p_serranidae_ga <- sma(ga~SL, data=p[which(p$Family=="Serranidae"), ], log="xy",
                     method="SMA", robust=F, slope.test=2)
check_assump(carangidae_ga)
check_assump(lutjanidae_ga)
check_assump(p_serranidae_ga)

# Benthic invertivore families/species:
monotaxis_ga <- sma(ga~SL, data=b[which(b$SpeciesCode=="MO.GRAN"), ], log="xy",
                    method="SMA", robust=F, slope.test=2)
paracirrhites_ga <- sma(ga~SL, data=b[which(b$SpeciesCode=="PA.ARCA"), ], 
                        log="xy", method="SMA", robust=F, slope.test=2)
parupeneus_ga <- sma(ga~SL, data=b[which(b$SpeciesCode=="PA.INSU"), ], log="xy",
                     method="SMA", robust=F, slope.test=2)

check_assump(monotaxis_ga)
check_assump(paracirrhites_ga)
check_assump(parupeneus_ga)

# Zooplanktivore families:
caesionidae_ga <- sma(ga~SL, data=z[which(z$Family=="Caesionidae"), ], 
                      log="xy", method="SMA", robust=F, slope.test=2)
pomacentridae_ga <- sma(ga~SL, data=z[which(z$Family=="Pomacentridae"), ], 
                        log="xy", method="SMA", robust=F, slope.test=2)
z_serranidae_ga <- sma(ga~SL, data=z[which(z$Family=="Serranidae"), ],
                     log="xy", method="SMA", robust=F, slope.test=2)

check_assump(caesionidae_ga)
check_assump(pomacentridae_ga)
check_assump(z_serranidae_ga)

# Herbivore families:
acanthuridae_ga <- sma(ga~SL, data=h[which(h$Family=="Acanthuridae"), ], 
                       log="xy", method="SMA", robust=F, slope.test=2)
pomacanthidae_ga <- sma(ga~SL, data=h[which(h$Family=="Pomacanthidae"), ],
                        log="xy", method="SMA", robust=F, slope.test=2)
scaridae_ga <- sma(ga~SL, data=h[which(h$Family=="Scaridae"), ], 
                   log="xy", method="SMA", robust=F, slope.test=2)

check_assump(acanthuridae_ga)
check_assump(pomacanthidae_ga)
check_assump(scaridae_ga)

# Corallivore:
chaetodon_ga <- sma(ga~SL, data=c[which(c$SpeciesCode=="CH.ORNA"), ], 
                    log="xy", method="SMA", robust=F, slope.test=2)

u <- mk_sma_summary(carangidae_ga, "carangidae_ga")
v <- mk_sma_summary(lutjanidae_ga, "lutjanidae_ga")
w <- mk_sma_summary(p_serranidae_ga, "serranidae_ga")
p_area <- cbind(u, v, w)

u <- mk_sma_summary(monotaxis_ga, "MO.GRAN_ga")
v <- mk_sma_summary(paracirrhites_ga, "PA.ARCA_ga")
w <- mk_sma_summary(parupeneus_ga, "PA.INSU_ga")
b_area <- cbind(u, v, w)

u <- mk_sma_summary(caesionidae_ga, "caesionidae_ga")
v <- mk_sma_summary(pomacentridae_ga, "pomacentridae_ga")
w <- mk_sma_summary(z_serranidae_ga, "serranidae_ga")
z_area <- cbind(u, v, w)

u <- mk_sma_summary(acanthuridae_ga, "acanthuridae_ga")
v <- mk_sma_summary(pomacanthidae_ga, "pomacanthidae_ga")
w <- mk_sma_summary(scaridae_ga, "scaridae_ga")
h_area <- cbind(u, v, w)

c_area <- mk_sma_summary(chaetodon_ga, "chaetodon_ga")

write.csv(p_area, "p_family_ga.csv")
write.csv(b_area, "b_family_ga.csv")
write.csv(z_area, "z_family_ga.csv")
write.csv(h_area, "h_family_ga.csv")
write.csv(c_area, "c_family_ga.csv")


# I think this is the same as what I have written above, but more simplified

gh_sma <- function(df) {
  sma(gh~SL, data=df, log="xy", method="SMA", robust=F, slope.test=1)
}

p_famGH <- dlply(p, .(Family), gh_sma)
z_famGH <- dlply(z, .(Family), gh_sma)
h_famGH <- dlply(h, .(Family), gh_sma)

p_famGH_summ <- ldply(p_famGH, .fun=mk_sma_summary(z, 
  group=as.character())))

p_famGH_assump <- l_ply(p_famGH, .fun=check_assump)
z_famGH_assump <- l_ply(z_famGH, .fun=check_assump)
h_famGH_assump <- l_ply(h_famGH, .fun=check_assump)

p_famGH_summ <- cbind(p_fam, mk_spp_summary(p_famGH, 3))
z_famGH_summ <- cbind(z_fam, mk_spp_summary(z_famGH, 3))
h_famGH_summ <- cbind(h_fam, mk_spp_summary(h_famGH, 3))


gw_sma <- function(df) {
  sma(gw~SL, data=df, log="xy", method="SMA", robust=F, slope.test=1)
}

p_famGW <- dlply(p, .(Family), gw_sma)
z_famGW <- dlply(z, .(Family), gw_sma)
h_famGW <- dlply(h, .(Family), gw_sma)

p_famGW_summ <- cbind(p_fam, mk_spp_summary(p_famGW, 3))
z_famGW_summ <- cbind(z_fam, mk_spp_summary(z_famGW, 3))
h_famGW_summ <- cbind(h_fam, mk_spp_summary(h_famGW, 3))


ga_sma <- function(df) {
  sma(ga~SL, data=df, log="xy", method="SMA", robust=F, slope.test=1)
}

p_famGA <- dlply(p, .(Family), ga_sma)
z_famGA <- dlply(z, .(Family), ga_sma)
h_famGA <- dlply(h, .(Family), ga_sma)

p_famGA_summ <- cbind(p_fam, mk_spp_summary(p_famGA, 3))
z_famGA_summ <- cbind(z_fam, mk_spp_summary(z_famGA, 3))
h_famGA_summ <- cbind(h_fam, mk_spp_summary(h_famGA, 3))


################################################################################
########                 Species level SMA analyses                     ########
################################################################################

# GH -----------------------------------------
gh_sma <- function(df) {
  sma(gh~SL, data=df, log="xy", method="SMA", robust=F, slope.test=1)
}

p_sppGH <- dlply(p, .(SpeciesCode), gh_sma)
b_sppGH <- dlply(b, .(SpeciesCode), gh_sma)
z_sppGH <- dlply(z, .(SpeciesCode), gh_sma)
h_sppGH <- dlply(h, .(SpeciesCode), gh_sma)
c_sppGH <- dlply(c, .(SpeciesCode), gh_sma)

p_sppGH_summ <- cbind(p_spp, mk_spp_summary(p_sppGH, 7))
b_sppGH_summ <- cbind(b_spp, mk_spp_summary(b_sppGH, 3))
z_sppGH_summ <- cbind(z_spp, mk_spp_summary(z_sppGH, 6))
h_sppGH_summ <- cbind(h_spp, mk_spp_summary(h_sppGH, 6))
c_sppGH_summ <- cbind("CH.ORNA", mk_spp_summary(c_sppGH, 1))

mk_smaSPP_graph_df(p_sppGH_summ, 7)

# Only for writing tables to csv:
p_sppGH_summ <- t(p_sppGH_summ)[-1, ]
colnames(p_sppGH_summ) <- p_spp
b_sppGH_summ <- t(b_sppGH_summ)[-1, ]
colnames(b_sppGH_summ) <- b_spp
z_sppGH_summ <- t(z_sppGH_summ)[-1, ]
colnames(z_sppGH_summ) <- z_spp
h_sppGH_summ <- t(h_sppGH_summ)[-1, ]
colnames(h_sppGH_summ) <- h_spp

write.csv(p_sppGH_summ, "p_spp_ghSMA.csv")
write.csv(z_sppGH_summ, "z_spp_ghSMA.csv")
write.csv(h_sppGH_summ, "h_spp_ghSMA.csv")


# GW -----------------------------------------

gw_sma <- function(df) {
  sma(gw~SL, data=df, log="xy", method="SMA", robust=F, slope.test=1)
}

p_sppGW <- dlply(p, .(SpeciesCode), gw_sma)
b_sppGW <- dlply(b, .(SpeciesCode), gw_sma)
z_sppGW <- dlply(z, .(SpeciesCode), gw_sma)
h_sppGW <- dlply(h, .(SpeciesCode), gw_sma)
c_sppGW <- dlply(c, .(SpeciesCode), gw_sma)

p_sppGW_summ <- cbind(p_spp, mk_spp_summary(p_sppGW, 7))
b_sppGW_summ <- cbind(b_spp, mk_spp_summary(b_sppGW, 3))
z_sppGW_summ <- cbind(z_spp, mk_spp_summary(z_sppGW, 6))
h_sppGW_summ <- cbind(h_spp, mk_spp_summary(h_sppGW, 6))
c_sppGW_summ <- cbind("CH.ORNA", mk_spp_summary(c_sppGW, 1))

# Only for writing tables to csv:
p_sppGW_summ <- t(p_sppGW_summ)[-1, ]
colnames(p_sppGW_summ) <- p_spp
b_sppGW_summ <- t(b_sppGW_summ)[-1, ]
colnames(b_sppGW_summ) <- b_spp
z_sppGW_summ <- t(z_sppGW_summ)[-1, ]
colnames(z_sppGW_summ) <- z_spp
h_sppGW_summ <- t(h_sppGW_summ)[-1, ]
colnames(h_sppGW_summ) <- h_spp

write.csv(p_sppGW_summ, "p_spp_gwSMA.csv")
write.csv(z_sppGW_summ, "z_spp_gwSMA.csv")
write.csv(h_sppGW_summ, "h_spp_gwSMA.csv")


# GA -----------------------------------------

ga_sma <- function(df) {
  sma(ga~SL, data=df, log="xy", method="SMA", robust=F, slope.test=2)
}

p_sppGA <- dlply(p, .(SpeciesCode), ga_sma)
b_sppGA <- dlply(b, .(SpeciesCode), ga_sma)
z_sppGA <- dlply(z, .(SpeciesCode), ga_sma)
h_sppGA <- dlply(h, .(SpeciesCode), ga_sma)
c_sppGA <- dlply(c, .(SpeciesCode), ga_sma)

p_sppGA_summ <- cbind(p_spp, mk_spp_summary(p_sppGA, 7))
b_sppGA_summ <- cbind(b_spp, mk_spp_summary(b_sppGA, 3))
z_sppGA_summ <- cbind(z_spp, mk_spp_summary(z_sppGA, 6))
h_sppGA_summ <- cbind(h_spp, mk_spp_summary(h_sppGA, 6))
c_sppGA_summ <- cbind("CH.ORNA", mk_spp_summary(c_sppGA, 1))

# Only for writing tables to csv:x`
p_sppGA_summ <- t(p_sppGA_summ)[-1, ]
colnames(p_sppGA_summ) <- p_spp
b_sppGA_summ <- t(b_sppGA_summ)[-1, ]
colnames(b_sppGA_summ) <- b_spp
z_sppGA_summ <- t(z_sppGA_summ)[-1, ]
colnames(z_sppGA_summ) <- z_spp
h_sppGA_summ <- t(h_sppGA_summ)[-1, ]
colnames(h_sppGA_summ) <- h_spp

write.csv(p_sppGA_summ, "p_spp_gaSMA.csv")
write.csv(z_sppGA_summ, "z_spp_gaSMA.csv")
write.csv(h_sppGA_summ, "h_spp_gaSMA.csv")


################################################################################
########                  Functional group SMA plots                    ########
################################################################################

# way to pass through SMA regression values to plot lines in ggplot2 - 
# Hadley Wickham answer from SO :)
library(ggplot2)
library(gridExtra)

height_df <- mk_sma_graph_df(fg_height) 
height_df$j_fg <- j_fg

width_df <- mk_sma_graph_df(fg_width)
width_df$j_fg <- j_fg

area_df <- mk_sma_graph_df(fg_area)
area_df$j_fg <- j_fg

fg_gh_plot <- mk_ghFG_SMAplot(pento, height_df)
fg_gh_facet <- mk_ghFG_SMAfacet(pento, height_df)

fg_gw_plot <-mk_gwFG_SMAplot(pento, width_df)
fg_gw_facet <-mk_gwFG_SMAfacet(pento, width_df)

fg_ga_plot <- mk_gaFG_SMAplot(pento, area_df)
fg_ga_facet <-mk_gaFG_SMAfacet(pento, area_df)


#
#
# Poster graph ####
ggplot(data = pento, aes(x = SL, y = ga, colour=j_fg)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  #scale_x_log10(limits=c(1, 1000)) +
  #xlab("log(standard length, mm)") +
  #ylab("log(vertical gape, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  scale_colour_discrete(name = "Functional \n Group") +
  #theme(axis.title.y = element_text(size = 30, vjust = 0.5)) +
  #theme(axis.title.x = element_text(size = 30, vjust = 0.1)) +
  #theme(axis.text.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 20)) +
  #theme(legend.text = element_text(size = 26)) +
  theme(legend.key.height = unit(1.5, "line")) +
  #theme(legend.title = element_text(size = 24)) +
  theme(legend.position = c(0.90, 0.35)) +
  theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
  geom_segment(data = area_df, aes(x = from, xend = to, y = yfrom, yend = yto)) 
#facet_wrap(~j_fg)



# Plots by Family 
# =============================================================================
# Piscivore fam plots
height_df <- mk_smaSPP_graph_df(p_famGH_summ, 3) 
names(height_df)[1] <- "Family"
p_famGH_plots <- mk_ghFam_SMAplot(p, height_df)
p_famGH_facets <- mk_ghFam_SMAfacets(p, height_df)

width_df <- mk_smaSPP_graph_df(p_famGW_summ, 3) 
names(width_df)[1] <- "Family"
p_famGW_plots <- mk_gwFam_SMAplot(p, width_df)
p_famGW_facets <- mk_gwFam_SMAfacets(p, width_df)

area_df <- mk_smaSPP_graph_df(p_famGA_summ, 3) 
names(area_df)[1] <- "Family"
p_famGA_plots <- mk_gaFam_SMAplot(p, area_df)
p_famGA_facets <- mk_gaFam_SMAfacets(p, area_df)

# Zooplanktivore fam plots
height_df <- mk_smaSPP_graph_df(z_famGH_summ, 3) 
names(height_df)[1] <- "Family"
z_famGH_plots <- mk_ghFam_SMAplot(z, height_df)
z_famGH_facets <- mk_ghFam_SMAfacets(z, height_df)

width_df <- mk_smaSPP_graph_df(z_famGW_summ, 3) 
names(width_df)[1] <- "Family"
z_famGW_plots <- mk_gwFam_SMAplot(z, width_df)
z_famGW_facets <- mk_gwFam_SMAfacets(z, width_df)

area_df <- mk_smaSPP_graph_df(z_famGA_summ, 3) 
names(area_df)[1] <- "Family"
z_famGA_plots <- mk_gaFam_SMAplot(z, area_df)
z_famGA_facets <- mk_gaFam_SMAfacets(z, area_df)

# Herbivore fam plots
height_df <- mk_smaSPP_graph_df(h_famGH_summ, 3) 
names(height_df)[1] <- "Family"
h_famGH_plots <- mk_ghFam_SMAplot(h, height_df)
h_famGH_facets <- mk_ghFam_SMAfacets(h, height_df)

width_df <- mk_smaSPP_graph_df(h_famGW_summ, 3) 
names(width_df)[1] <- "Family"
h_famGW_plots <- mk_gwFam_SMAplot(h, width_df)
h_famGW_facets <- mk_gwFam_SMAfacets(h, width_df)

area_df <- mk_smaSPP_graph_df(h_famGA_summ, 3) 
names(area_df)[1] <- "Family"
h_famGA_plots <- mk_gaFam_SMAplot(h, area_df)
h_famGA_facets <- mk_gaFam_SMAfacets(h, area_df)


# Plots by Species
# =============================================================================
# Piscivore spp plots
height_df <- mk_smaSPP_graph_df(p_sppGH_summ, 7) 
names(height_df)[1] <- "SpeciesCode"
p_sppGH_plots <- mk_gh_SMAplot(p, height_df)
p_sppGH_facets <- mk_gh_SMAplot_facets(p, height_df)

width_df <- mk_smaSPP_graph_df(p_sppGW_summ, 7) 
names(width_df)[1] <- "SpeciesCode"
width_df$SpeciesCode <- p_spp
p_sppGW_plots <- mk_gw_SMAplot(p, width_df)
p_sppGW_facets <- mk_gw_SMAplot_facets(p, width_df)

area_df <- mk_smaSPP_graph_df(p_sppGA_summ, 7) 
names(area_df)[1] <- "SpeciesCode"
mk_ga_SMAplot(p, area_df)
p_sppGA_plots <- mk_ga_SMAplot(p, area_df)
p_sppGA_facets <- mk_ga_SMAplot_facets(p, area_df)

# Benthic Invertivore spp plots
height_df <- mk_smaSPP_graph_df(b_sppGH_summ, 3) 
names(height_df)[1] <- "SpeciesCode"
mk_gh_SMAplot(b, height_df)
b_sppGH_plots <- mk_gh_SMAplot(b, height_df)
b_sppGH_facets <- mk_gh_SMAplot_facets(b, height_df)


width_df <- mk_smaSPP_graph_df(b_sppGW_summ, 3) 
names(width_df)[1] <- "SpeciesCode"
mk_gw_SMAplot(b, width_df)
b_sppGW_plots <- mk_gw_SMAplot(b, width_df)
b_sppGW_facets <- mk_gw_SMAplot_facets(b, width_df)

area_df <- mk_smaSPP_graph_df(b_sppGA_summ, 3) 
names(area_df)[1] <- "SpeciesCode"
mk_ga_SMAplot(b, area_df)
b_sppGA_plots <- mk_ga_SMAplot(b, area_df)
b_sppGA_facets <- mk_ga_SMAplot_facets(b, area_df)

# Zooplanktivore spp plots
height_df <- mk_smaSPP_graph_df(z_sppGH_summ, 6) 
names(height_df)[1] <- "SpeciesCode"
mk_gh_SMAplot(z, height_df)
z_sppGH_plots <- mk_gh_SMAplot(z, height_df)
z_sppGH_facets <- mk_gh_SMAplot_facets(z, height_df)

width_df <- mk_smaSPP_graph_df(z_sppGW_summ, 6) 
names(width_df)[1] <- "SpeciesCode"
mk_gw_SMAplot(z, width_df)
z_sppGW_plots <- mk_gw_SMAplot(z, width_df)
z_sppGW_facets <- mk_gw_SMAplot_facets(z, width_df)


area_df <- mk_smaSPP_graph_df(z_sppGA_summ, 6) 
names(area_df)[1] <- "SpeciesCode"
mk_ga_SMAplot(z, area_df)
z_sppGA_plots <- mk_ga_SMAplot(z, area_df)
z_sppGA_facets <- mk_ga_SMAplot_facets(z, area_df)


# Herbivore spp plots
height_df <- mk_smaSPP_graph_df(h_sppGH_summ, 6) 
names(height_df)[1] <- "SpeciesCode"
mk_gh_SMAplot(h, height_df)
h_sppGH_plots <- mk_gh_SMAplot(h, height_df)
h_sppGH_facets <- mk_gh_SMAplot_facets(h, height_df)

width_df <- mk_smaSPP_graph_df(h_sppGW_summ, 6) 
names(width_df)[1] <- "SpeciesCode"
mk_gw_SMAplot(h, width_df)
h_sppGW_plots <- mk_gw_SMAplot(h, width_df)
h_sppGW_facets <- mk_gw_SMAplot_facets(h, width_df)

area_df <- mk_smaSPP_graph_df(h_sppGA_summ, 6) 
names(area_df)[1] <- "SpeciesCode"
mk_ga_SMAplot(h, area_df)
h_sppGA_plots <- mk_ga_SMAplot(h, area_df)
h_sppGA_facets <- mk_ga_SMAplot_facets(h, area_df)

ggplot(data = p, aes(x = SL, y = gw, colour=SpeciesCode)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10(limits=c(1, 1000)) +
  xlab("log(standard length, mm)") +
  #ylab("log(vertical gape, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  scale_colour_discrete(name = "Functional \n Group") +
  theme(legend.key.height = unit(1.5, "line")) +
  theme(legend.position = c(0.90, 0.35)) +
  theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
  geom_segment(data = height_df, aes(x = from, xend = to, y = yfrom, yend = yto)) 

# Corallivore plots:
height_df <- mk_smaSPP_graph_df(c_sppGH_summ, 1)
names(height_df)[1] <- "SpeciesCode"
mk_gh_SMAplot(c, height_df)
c_sppGH_facets <- mk_gh_SMAplot_facets(c, height_df)

width_df <- mk_smaSPP_graph_df(c_sppGW_summ, 1)
names(width_df)[1] <- "SpeciesCode"
mk_gw_SMAplot(c, width_df)
c_sppGW_facets <- mk_gw_SMAplot_facets(c, width_df)

area_df <- mk_smaSPP_graph_df(c_sppGA_summ, 1)
names(area_df)[1] <- "SpeciesCode"
mk_ga_SMAplot(c, area_df)
c_sppGA_facets <- mk_ga_SMAplot_facets(c, area_df)


# Saving plots to PDFs
pdf(file = "sma_plots.pdf", height = 13, width = 9)

# Functional groups overall
grid.arrange(fg_gh_plot, fg_gw_plot, fg_ga_plot)
grid.arrange(fg_gh_facet, fg_gw_facet, fg_ga_facet)

# piscivore fam
grid.arrange(p_famGH_plots, p_famGW_plots, p_famGA_plots)
grid.arrange(p_famGH_facets, p_famGW_facets, p_famGA_facets)

# piscivore spp
grid.arrange(p_sppGH_plots, p_sppGW_plots, p_sppGA_plots)
grid.arrange(p_sppGH_facets, p_sppGW_facets, p_sppGA_facets)

# benthic invertivore spp
grid.arrange(b_sppGH_plots, b_sppGW_plots, b_sppGA_plots)
grid.arrange(b_sppGH_facets, b_sppGW_facets, b_sppGA_facets)

# zooplanktivore fam
grid.arrange(z_famGH_plots, z_famGW_plots, z_famGA_plots)
grid.arrange(z_famGH_facets, z_famGW_facets, z_famGA_facets)

# zooplanktivore spp
grid.arrange(z_sppGH_plots, z_sppGW_plots, z_sppGA_plots)
grid.arrange(z_sppGH_facets, z_sppGW_facets, z_sppGA_facets)

# herbivore fam
grid.arrange(h_famGH_plots, h_famGW_plots, h_famGA_plots)
grid.arrange(h_famGH_facets, h_famGW_facets, h_famGA_facets)

# herbivore spp
grid.arrange(h_sppGH_plots, h_sppGW_plots, h_sppGA_plots)
grid.arrange(h_sppGH_facets, h_sppGW_facets, h_sppGA_facets)

# corallivore spp
grid.arrange(c_sppGH_facets, c_sppGW_facets, c_sppGA_facets)

dev.off()


################################################################################
########                  Corallivore Random Effects                    ########
################################################################################

cGA <- sma(ga~SL*Region, data=c, log="xy", method="SMA", robust=T, slope.test=1)

cGA_summ <- mk_spp_summary(cGA, grouping=T)
cGA_reg_graphing <- mk_smaSPP_graph_df(cGA_summ, 5)
names(cGA_reg_graphing)[1] <- "Region"

A <- dev.cur()
B <- dev.cur()

A <- dev.set(which=2)
B <- dev.set(which=3)

plot.new()

cGA_SMAplot <- 
  ggplot(data = c, aes(x = SL, y = ga, colour=Region)) +
    geom_point() +
    geom_text(position=position_jitter(w=0.01, h=0.01), aes(label=dissected_by), 
      size=3) +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
    geom_segment(data = cGA_reg_graphing, aes(x = from, xend = to, y = yfrom, yend = yto))
cGA_SMAplot

mk_ghFG_SMAfacet(pento, height_df, gapeType="gh")

mk_ghReg_SMAplot(c, cGH_reg_graphing)


mk_ghReg_SMAplot <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = gh, colour = Region) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab("log(vertical gape, mm)") +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, 
      yend = yto))
}



ggplot(data = c, aes(x = SL, y = gh, colour=Region)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  xlab("log(standard length, mm)") +
  ylab("log(vertical gape, mm)") +
  #ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  scale_colour_discrete(name = "Functional \n Group") +
  theme(legend.key.height = unit(1.5, "line")) +
  theme(legend.position = c(0.90, 0.35)) +
  theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
  geom_segment(data = cGH_reg_graphing, aes(x = from, xend = to, y = yfrom, yend = yto)) 










































