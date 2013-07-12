################################################################################
############             Entering and Cleaning Data                 ############
################################################################################

# Remove PS.COOP and AC.TRIOS because low n, but more so because there is no 
# size gradient
# Which rows have PS.COOP and AC.TRIOS?
psc <- which(fish$SpeciesCode==('PS.COOP'))

#act <- which(fish$SpeciesCode==('AC.TRIOS')) 
# removing these rows from fish:
fish <- fish[-c(psc, act), ]

# Including area and relative gape sizes:
fish$ga <- with(fish, pi*(gh/2)*(gw/2))

fish$gh_ratio <- fish$gh/fish$SL
fish$gw_ratio <- fish$gw/fish$SL
fish$ga_ratio <- fish$ga/fish$SL



################################################################################
############             Organizing Piscivores                      ############
################################################################################

p <- fish[which(fish$j_fg=='Pi'), ]

# Setting factors for use in analyses
p$SpeciesCode <- factor(p$SpeciesCode)
p$Family <- factor(p$Family)

p_fam <- c("Carangidae", "Lutjanidae", "Serranidae")
p_fam <- factor(p_fam, levels=p_fam)

# Ordering by Family to make my life easier

p$SpeciesCode <- factor(p$SpeciesCode, levels = c("CA.MELA", 
                                                  "AP.FURC", 
                                                  "LU.BOHA",
                                                  "LU.KASM",
                                                  "CE.ARGU",
                                                  "CE.UROD",
                                                  "VA.LOUT")
)
p_spp <- c("CA.MELA", "AP.FURC", "LU.BOHA", "LU.KASM", "CE.ARGU", "CE.UROD", 
           "VA.LOUT")
p_spp <- factor(p_spp, levels=p_spp)

################################################################################
############             Organizing Benthic Invertivores            ############
################################################################################
b <- fish[which(fish$j_fg=='BI'), ]

b$SpeciesCode <- factor(b$SpeciesCode)
b$Family <- factor(b$Family)

# Ordering by Family to make my life easier 
b$SpeciesCode <- factor(b$SpeciesCode, levels = c("MO.GRAN",
                                                  "PA.ARCA",
                                                  "PA.INSU"
                                                  )
)

b_spp <- c("MO.GRAN", "PA.ARCA", "PA.INSU")
b_spp <- factor(b_spp, levels=b_spp)


################################################################################
############             Organizing Herbivores                      ############
################################################################################


h <- fish[which(fish$j_fg=='He'), ] 

h$SpeciesCode <- factor(h$SpeciesCode)
h$Family <- factor(h$Family)

h_fam <- c("Acanthuridae", "Pomacanthidae", "Scaridae")
h_fam <- factor(h_fam, levels=h_fam)

### ORDERING H by Family then Species --> this will make organizing facets 
# easier 

h$SpeciesCode <- factor(h$SpeciesCode, levels = c("AC.NIGR",
                                                  "AC.OLIV",
                                                  "CE.FLAV",
                                                  "CH.SORD",
                                                  "SC.FREN",
                                                  "SC.RUBR"
                                                  )
)

h_spp <- c("AC.NIGR", "AC.OLIV", "CE.FLAV", "CH.SORD", "SC.FREN", "SC.RUBR")
h_spp <- factor(h_spp, levels = h_spp)

################################################################################
############             Organizing Zooplanktivores                 ############
################################################################################

zp <- fish[which(fish$j_fg=='ZP'), ]

zp$SpeciesCode <- factor(zp$SpeciesCode)
zp$Family <- factor(zp$Family)

z_fam <- c("Caesionidae", "Pomacentridae", "Serranidae")
z_fam <- factor(z_fam, levels=z_fam)

# Ordering by Family to make my life easier
zp$SpeciesCode <- factor(zp$SpeciesCode, levels = c("CA.TERE",
                                                    "PT.TILE",
                                                    "CH.VAND",
                                                    "PS.BART",
                                                    "PS.DISP",
                                                    "PS.OLIV")
)

z_spp <- c("CA.TERE", "PT.TILE", "CH.VAND", "PS.BART", "PS.DISP", "PS.OLIV")
z_spp <- factor(z_spp, levels=z_spp)
z <- zp


################################################################################
############             Organizing Corallivores                 ############
################################################################################


c <- fish[fish$j_fg=='C', ]

c$SpeciesCode <- factor(c$SpeciesCode)
c$Family <- factor(c$Family)

# Ordering by Family to make my life easier 
c$SpeciesCode <- factor(c$SpeciesCode, levels = c("CH.ORNA")
)


################################################################################
############             Organizing 5 FGs SPP and Fam               ############
################################################################################

# picking just Pi, He, BI, and ZP and ordering them:
pento <- fish[fish$j_fg %in% c('Pi', 'He', 'BI', 'ZP', 'C'), ]

pento$j_fg <- factor(pento$j_fg, levels=c('Pi', 'BI', 'ZP', 'He', 'C'))

pento$Family <- factor(pento$Family, levels=c("Carangidae",
                                            "Lutjanidae",
                                            "Serranidae",
                                            "Cirrhitidae",
                                            "Lethrinidae",
                                            "Mullidae",
                                            "Acanthuridae",
                                            "Pomacanthidae",
                                            "Scaridae",
                                            "Caesionidae",
                                            "Pomacentridae",
                                            "Chaetodontidae")
)

pento$SpeciesCode <- factor(pento$SpeciesCode, levels = c("CA.MELA", 
                                                        "AP.FURC", 
                                                        "LU.BOHA",
                                                        "LU.KASM",
                                                        "CE.ARGU",
                                                        "CE.UROD",
                                                        "VA.LOUT",
                                                        "PA.ARCA",
                                                        "MO.GRAN",
                                                        "PA.INSU",
                                                        "AC.NIGR",
                                                        "AC.OLIV",
                                                        "CE.FLAV",
                                                        "CH.SORD",
                                                        "SC.FREN",
                                                        "SC.RUBR",
                                                        "CA.TERE",
                                                        "PT.TILE",
                                                        "CH.VAND",
                                                        "PS.BART",
                                                        "PS.DISP",
                                                        "PS.OLIV",
                                                        "CH.ORNA")
)

pento_order <- c("CA.MELA", "AP.FURC", "LU.BOHA", "LU.KASM", "CE.ARGU",
                 "CE.UROD", "VA.LOUT", "PA.ARCA", "MO.GRAN", "PA.INSU",
                 "AC.NIGR", "AC.OLIV", "CE.FLAV", "CH.SORD", "SC.FREN",
                 "SC.RUBR", "CA.TERE", "PT.TILE", "CH.VAND", "PS.BART",
                 "PS.DISP", "PS.OLIV", "CH.ORNA")

################################################################################
############             Cleaning Prey Size Data                    ############
################################################################################

to_remove <- which(prey$Prey!=('y'))
no_zps <- which(prey$j_fg==('ZP'))
prey$j_fg <- factor(prey$j_fg, levels=c('Pi', 'BI', 'ZP', 'C', 'He'))

# remove the EPITAU005A which I haven't figured out what its prey note should be:
# had a cephalopod beak, but presumably the original animal was bigger than the
# beak
epitau5 <- which(prey$SpecimenID==('EPITAU005A'))
prey_zps <- prey[-c(to_remove, epitau5, no_zps), ]

prey1 <- data.frame(specimen = prey$SpecimenID,
                            SpeciesCode = prey$SpeciesCode,
                            family = prey$Family,
                            fg = prey$j_fg,
                            sl = as.numeric(as.character(prey$SL)),
                            gh = as.numeric(as.character(prey$gh)),
                            gw = as.numeric(as.character(prey$gw)),
                            ptype = prey$pType,
                            psize = as.numeric(as.character(prey$pSize))
                            )

# without zoops <1mm
prey2 <- data.frame(specimen = prey_zps$SpecimenID,
                    SpeciesCode = prey_zps$SpeciesCode,
                    family = prey_zps$Family,
                    fg = prey_zps$j_fg,
                    sl = as.numeric(as.character(prey_zps$SL)),
                    gh = as.numeric(as.character(prey_zps$gh)),
                    gw = as.numeric(as.character(prey_zps$gw)),
                    ptype = prey_zps$pType,
                    psize = as.numeric(as.character(prey_zps$pSize))
)
prey3 <- subset(prey2, fg == 'Pi' | fg == 'BI')

prey_gh <- subset(prey3, is.na(prey3$gh)==FALSE)

prey4 <- subset(prey2, fg == 'Pi' | fg == 'BI' | fg == 'ZP' | fg == 'Om')


# Cleaning prey frequency data ####
# selecting just BI and PI from all data:
freq1 <- data.frame(specimen = freq$SpecimenID,
                    SpeciesCode = freq$SpeciesCode,
                    family = freq$Family,
                    fg = freq$j_fg,
                    sl = as.numeric(as.character(freq$SL)),
                    gh = as.numeric(as.character(freq$gh)),
                    gw = as.numeric(as.character(freq$gw)),
                    invert = freq$invert,
                    fish = freq$fish,
                    ptype = freq$pType,
                    psize = as.numeric(as.character(freq$pSize))
)

freq_bp <- subset(freq1, (fg == 'Pi' | fg == 'BI') & is.na(freq1$invert)==FALSE)

# including only species that were used in the quantile regression pred-prey 
# size analysis (for honours)
myspp <- freq_bp[freq_bp$SpeciesCode %in% c("CE.UROD", "CE.ARGU", "LU.BOHA", "LU.KASM",
                                "CA.MELA", "VA.LOUT", "MO.GRAN", "PA.INSU", 
                                "PA.ARCA", "EP.MACU", "EP.HEXA", "EP.TAUV", 
                                "EP.SPIL", "CA.ORTH", "AP.VIRE"), ]

######################
# Species specific dataframes:
ac.nigr <- fish[(which(fish$SpeciesCode == "AC.NIGR")), ]
ac.oliv <- fish[(which(fish$SpeciesCode == "AC.OLIV")), ]
ac.trios <- fish[(which(fish$SpeciesCode == "AC.TRIOS")), ]
ap.furc <- fish[(which(fish$SpeciesCode == "AP.FURC")), ]
ca.mela <- fish[(which(fish$SpeciesCode == "CA.MELA")), ]
ca.tere <- fish[(which(fish$SpeciesCode == "CA.TERE")), ]
ce.argu <- fish[(which(fish$SpeciesCode == "CE.ARGU")), ]
ce.flav <- fish[(which(fish$SpeciesCode == "CE.FLAV")), ]
ce.urod <- fish[(which(fish$SpeciesCode == "CE.UROD")), ]
ch.auri <- fish[(which(fish$SpeciesCode == "CH.AURI")), ]
ch.orna <- fish[(which(fish$SpeciesCode == "CH.ORNA")), ]
ch.sord <- fish[(which(fish$SpeciesCode == "CH.SORD")), ]
ch.vand <- fish[(which(fish$SpeciesCode == "CH.VAND")), ]
ct.marg <- fish[(which(fish$SpeciesCode == "CT.MARG")), ]
lu.boha <- fish[(which(fish$SpeciesCode == "LU.BOHA")), ]
lu.kasm <- fish[(which(fish$SpeciesCode == "LU.KASM")), ]
me.nige <- fish[(which(fish$SpeciesCode == "ME.NIGE")), ]
mo.gran <- fish[(which(fish$SpeciesCode == "MO.GRAN")), ]
pa.arca <- fish[(which(fish$SpeciesCode == "PA.ARCA")), ]
pa.insu <- fish[(which(fish$SpeciesCode == "PA.INSU")), ]
pl.dick <- fish[(which(fish$SpeciesCode == "PL.DICK")), ]
ps.bart <- fish[(which(fish$SpeciesCode == "PS.BART")), ]
ps.coop <- fish[(which(fish$SpeciesCode == "PS.COOP")), ]
ps.disp <- fish[(which(fish$SpeciesCode == "PS.DISP")), ]
ps.oliv <- fish[(which(fish$SpeciesCode == "PS.OLIV")), ]
pt.tile <- fish[(which(fish$SpeciesCode == "PT.TILE")), ]
sc.fren <- fish[(which(fish$SpeciesCode == "SC.FREN")), ]
sc.rubr <- fish[(which(fish$SpeciesCode == "SC.RUBR")), ]
va.lout <- fish[(which(fish$SpeciesCode == "VA.LOUT")), ]
























