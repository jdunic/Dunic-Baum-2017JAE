################################################################################
############             Entering and Cleaning Data                 ############
################################################################################

# Remove PS.COOP and AC.TRIOS because low n, but more so because there is no 
# size gradient
# Which rows have PS.COOP and AC.TRIOS?
ps_coop <- which(fish$SpeciesCode==('PS.COOP'))
ct_marg <- which(fish$SpeciesCode==('CT.MARG'))
ch_auri <- which(fish$SpeciesCode == 'CH.AURI')
me_nige <- which(fish$SpeciesCode == 'ME.NIGE')
ch_marg <- which(fish$SpeciesCode == 'CH.MARG')
ct_marg <- which(fish$SpeciesCode == 'CT.MARG')
ct_stri <- which(fish$SpeciesCode == 'CT.STRI')
ce_lori <- which(fish$SpeciesCode == 'CE.LORI')
pl_dick <- which(fish$SpeciesCode == 'PL.DICK')
#extra_sites <- which(fish$Region == 'EXTRA')
site_unknown <- which(fish$Site == 'Site Not Certain')
lu_kasm <- which(fish$SpeciesCode == 'LU.KASM')
kif12_130 <- which(fish$SpecimenID == 'KIF12_130')
kif11_319 <- which(fish$SpecimenID == 'KIF11_319')
kif12_050 <- which(fish$SpecimenID == "KIF12_050") # obvious outlier
ao <- which(fish$dissected_by == "AO" | fish$dissected_by == "angeleen")


#act <- which(fish$SpeciesCode==('AC.TRIOS')) 
# removing these rows from fish:
fish <- fish[-c(ps_coop, ch_marg, ch_auri, me_nige, ct_marg, ct_stri, ce_lori, 
                lu_kasm, site_unknown, kif12_130, kif11_319, kif12_050, ao, 
                pl_dick), ]

# Remove fish where site == "Site Not Certain"
#no_site <- which(fish$Site == "Site Not Certain")
#fish <- fish[-no_site, ]

# Including area and relative gape sizes:
fish$ga <- with(fish, pi*(gh/2)*(gw/2))

fish$gh_ratio <- fish$gh/fish$SL
fish$gw_ratio <- fish$gw/fish$SL
fish$ga_ratio <- fish$ga/fish$SL


fish$j_fg <- factor(fish$j_fg, levels = c("Pi", "BI", "ZP", "He", 
                                          "C", "Om", "De"))

# To allow observers to be referenced anonymously
fish$observer_id <- as.character(as.numeric(fish$dissected_by))

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
                                                  #"LU.KASM",
                                                  "CE.ARGU",
                                                  "CE.UROD",
                                                  "VA.LOUT")
)
#p_spp <- c("CA.MELA", "AP.FURC", "LU.BOHA", "LU.KASM", "CE.ARGU", "CE.UROD", 
#           "VA.LOUT")
p_spp <- c("CA.MELA", "AP.FURC", "LU.BOHA", "CE.ARGU", "CE.UROD", "VA.LOUT")
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
#ps_coop <- which(pento$SpeciesCode == "PS.COOP")
#ch_marg <- which(pento$SpeciesCode == "CH.MARG")
#ch_auri <- which(pento$SpeciesCode == "CH.AURI")
#me_nige <- which(pento$SpeciesCode == "ME.NIGE")
#ce_lori <- which(pento$SpeciesCode == "CE.LORI")

pento <- pento[-c(ps_coop, ch_marg, ch_auri, me_nige, ce_lori), ]

pento$j_fg <- factor(pento$j_fg, levels=c('Pi', 'BI', 'ZP', 'He', 'C'))

pento$Family <- factor(pento$Family, levels=c("Carangidae",
                                            "Lutjanidae",
                                            "Serranidae",
                                            "Cirrhitidae",
                                            "Lethrinidae",
                                            "Mullidae",
                                            "Caesionidae",
                                            "Pomacentridae",
                                            "Acanthuridae",
                                            "Pomacanthidae",
                                            "Scaridae",
                                            "Chaetodontidae")
)

pento$SpeciesCode <- factor(pento$SpeciesCode, levels = c("CA.MELA", 
                                                        "AP.FURC", 
                                                        "LU.BOHA",
                                                        #"LU.KASM",
                                                        "CE.ARGU",
                                                        "CE.UROD",
                                                        "VA.LOUT",
                                                        "PA.ARCA",
                                                        "MO.GRAN",
                                                        "PA.INSU",
                                                        "CA.TERE",
                                                        "PT.TILE",
                                                        "CH.VAND",
                                                        "PS.BART",
                                                        "PS.DISP",
                                                        "PS.OLIV",
                                                        "AC.NIGR",
                                                        "AC.OLIV",
                                                        "CE.FLAV",
                                                        "CH.SORD",
                                                        "SC.FREN",
                                                        "SC.RUBR",
                                                        "CH.ORNA")
)

#pento_order <- c("CA.MELA", "AP.FURC", "LU.BOHA", "LU.KASM", "CE.ARGU",
#                 "CE.UROD", "VA.LOUT", "PA.ARCA", "MO.GRAN", "PA.INSU",
#                 "CA.TERE", "PT.TILE", "CH.VAND", "PS.BART", "PS.DISP", 
#                 "PS.OLIV", "CH.ORNA", "AC.NIGR", "AC.OLIV", "CE.FLAV", 
#                 "CH.SORD", "SC.FREN", "SC.RUBR")

pento_order <- c("CA.MELA", "AP.FURC", "LU.BOHA", "CE.ARGU", "CE.UROD", 
                 "VA.LOUT", "PA.ARCA", "MO.GRAN", "PA.INSU", "CA.TERE", 
                 "PT.TILE", "CH.VAND", "PS.BART", "PS.DISP", "PS.OLIV", 
                 "CH.ORNA", "AC.NIGR", "AC.OLIV", "CE.FLAV", "CH.SORD", 
                 "SC.FREN", "SC.RUBR")

pento_order <- factor(pento_order, levels = pento_order)

p_spp_dfs <- split(p, p$SpeciesCode, drop=TRUE)
b_spp_dfs <- split(b, b$SpeciesCode, drop=TRUE)
z_spp_dfs <- split(z, z$SpeciesCode, drop=TRUE)
h_spp_dfs <- split(h, h$SpeciesCode, drop=TRUE)

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

prey5 <- subset(prey2, fg == 'Pi')

prey6 <- subset(prey2, fg == 'BI')


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
myspp <- freq_bp[freq_bp$SpeciesCode %in% c("CE.UROD", "CE.ARGU", "LU.BOHA", #"LU.KASM",
                                "CA.MELA", "VA.LOUT", "MO.GRAN", "PA.INSU", 
                                "PA.ARCA", "EP.MACU", "EP.HEXA", "EP.TAUV", 
                                "EP.SPIL", "CA.ORTH", "AP.VIRE"), ]

