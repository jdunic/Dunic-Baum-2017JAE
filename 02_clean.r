################################################################################
############             Entering and Cleaning Data                 ############
################################################################################

# loading libraries
library('plyr')

# Remove PS.COOP and AC.TRIOS because low n, but more so because there is no size gradient
# Which rows have PS.COOP and AC.TRIOS?
psc <- which(fish$SpeciesCode==('PS.COOP'))
act <- which(fish$SpeciesCode==('AC.TRIOS')) 
# removing these rows from fish:
fish <- fish[-c(psc, act),]

# Including area and relative gape sizes:
fish$ga <- with(fish, pi*(gh/2)*(gw/2))

fish$gh_ratio <- fish$gh/fish$SLMM
fish$gw_ratio <- fish$gw/fish$SLMM
fish$ga_ratio <- fish$ga/fish$SLMM



################################################################################
############             Organizing Piscivores                      ############
################################################################################

p <- fish[which(fish$j_fg=='Pi'),]

# Setting factors for use in analyses
p$SpeciesCode <- factor(p$SpeciesCode)
p$Family <- factor(p$Family)

# Ordering by Family to make my life easier
p$SpeciesCode <- factor(p$SpeciesCode, levels = c("CA.MELA", 
                                                  "AP.FURC", 
                                                  "LU.BOHA",
                                                  "LU.KASM",
                                                  "CE.ARGU",
                                                  "CE.UROD",
                                                  "VA.LOUT")
)

################################################################################
############             Organizing Benthic Invertivores            ############
################################################################################
b <- fish[which(fish$j_fg=='BI'),]

b$SpeciesCode <- factor(b$SpeciesCode)
b$Family <- factor(b$Family)

# Ordering by Family to make my life easier 
b$SpeciesCode <- factor(b$SpeciesCode, levels = c("PA.ARCA",
                                                  "MO.GRAN",
                                                  "PA.INSU"
                                                  )
)


################################################################################
############             Organizing Herbivores                      ############
################################################################################


h <- fish[which(fish$j_fg=='He'),] 

h$SpeciesCode <- factor(h$SpeciesCode)
h$Family <- factor(h$Family)

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

################################################################################
############             Organizing Zooplanktivores                 ############
################################################################################

zp <- fish[which(fish$j_fg=='ZP'),]

zp$SpeciesCode <- factor(zp$SpeciesCode)
zp$Family <- factor(zp$Family)

# Ordering by Family to make my life easier
zp$SpeciesCode <- factor(zp$SpeciesCode, levels = c("CA.TERE",
                                                    "PT.TILE",
                                                    "CH.VAND",
                                                    "PS.BART",
                                                    "PS.DISP",
                                                    "PS.OLIV")
)

################################################################################
############             Organizing Corallivores                 ############
################################################################################


c <- fish[fish$j_fg=='C',]

c$SpeciesCode <- factor(c$SpeciesCode)
c$Family <- factor(c$Family)

# Ordering by Family to make my life easier 
c$SpeciesCode <- factor(c$SpeciesCode, levels = c("CH.ORNA")
)


################################################################################
############             Organizing 5 FGs SPP and Fam               ############
################################################################################

# picking just Pi, He, BI, and ZP and ordering them:
pento <- fish[fish$j_fg %in% c('Pi', 'He', 'BI', 'ZP', 'C'),]

pento$j_fg <- factor(pento$j_fg, levels=c('Pi', 'BI', 'ZP', 'C', 'He'))

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

# remove the EPITAU005A which I haven't figured out what its prey note should be:
# had a cephalopod beak, but presumably the original animal was bigger than the
# beak
epitau5 <- which(prey$SpecimenID==('EPITAU005A'))
prey_zps <- prey[-c(to_remove, epitau5, no_zps),]

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

myspp <- freq_bp[freq_bp$SpeciesCode %in% c("CE.UROD", "CE.ARGU", "LU.BOHA", "LU.KASM",
                                "CA.MELA", "VA.LOUT", "MO.GRAN", "PA.INSU", 
                                "PA.ARCA", "EP.MACU", "EP.HEXA", "EP.TAUV", 
                                "EP.SPIL", "CA.ORTH", "AP.VIRE"), ]


























