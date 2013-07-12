##############################################################################################
################        Preliminary Gape Size - Body Size analysis           ################
##############################################################################################

# set working directory
setwd('/Users/jillian/R_projects/Allometry/')

########### Entering and Cleaning Data ##########

# Data for all specimens where a species has more than 3 samples with gape and SL data
fish <- read.csv('Allometry Data/fish_gape_cleaned_more_than_3_jan_24.csv', header = TRUE, sep = ",") 

# Enter Fish dimension data for gape - body size analysis
# imported as a dataframe

colnames(fish)

# Remove PS.COOP and AC.TRIOS because low n, but more so because there is no size gradient
  # Which rows have PS.COOP and AC.TRIOS?
psc <- which(fish$SpeciesCode==('PS.COOP'))
act <- which(fish$SpeciesCode==('AC.TRIOS')) 
  # removing these rows from fish:
fish <- fish[-c(psc, act),]

### Checking that there are no crazy or NA/NULL wt, gh, gw, or SLMM
max(fish$wt)
max(fish$gh)
max(fish$gw)
max(fish$SLMM)

##############################################################################################
################               Gape Size - Body Size Exploration              ################
##############################################################################################

# Gives me n for each guild
APCnt <- length(which(fish$Guild=='AP'))
BICnt <- length(which(fish$Guild=='BI'))
DeCnt <- length(which(fish$Guild=='De'))
HeCnt <- length(which(fish$Guild=='He'))
OmCnt <- length(which(fish$Guild=='Om'))
PiCnt <- length(which(fish$Guild=='Pi'))
ZPCnt <- length(which(fish$Guild=='ZP'))
CoCnt <- length(which(fish$Guild=='C'))


# Vector of n of each guild
c(APCnt, BICnt, DeCnt, HeCnt, OmCnt, PiCnt, ZPCnt, CoCnt)

# Dataframe of n specimens in each guild
cnt.fg <- c(APCnt, BICnt, DeCnt, HeCnt, OmCnt, PiCnt, ZPCnt, CoCnt)
fg.name <- c('AP', 'BI', 'De', 'He', 'Om', 'Pi', 'ZP', 'Co')
fg.cnt <- data.frame(fg=fg.name, cnt=cnt.fg)
fg.cnt

# Total number of each species - prints SpeciesCode and Count
unique(fish$SpeciesCode)

for (i in unique(fish$SpeciesCode)) {
    a <- length(which((fish$SpeciesCode)==i))
    print(c(i, a))
  }

# Sample size for each Species with n > 5 
# this function gives me the number of times each species appears :| that only took over an hour
# of searching for :(
spp.cnt <- as.data.frame(table(fish$SpeciesCode))



#### TIME YOUR CODE :D
start = proc.time()
## Code you wish to time goes here.
end = proc.time()
total = (end[3]-start[3])
cat(paste("Time Elapsed: ", total, "\n"))

############### Making some plots!! ####################

library("ggplot2")      # load ggplot
library('plyr')         # load plyr

################### Plots with all species ##################

# Creates a dataframe where there is a character expression n = freq(SpeciesCode) because of
# the 'summarize'
df.n <- ddply(.data=fish, .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))

# For each species, plot gape size~body size across all regions. Therefore use:
# <All specimens with gp and SL for exploration_cleaned.csv >
# include n = sample size

# Plotting Vertical Gape
ggplot(data=fish, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ SpeciesCode, ncol=5) +
  geom_text(data=df.n, aes(x=3.6, y=4.2, label=n), parse=TRUE) + 
  # adds n = on each facet
  ggtitle("Vertical Gape") 
  facet_grid(. ~ cyl, labeller = label_both)

# Plotting Horizontal Gape
ggplot(data=fish, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  facet_wrap(~ SpeciesCode, ncol = 5) +
  geom_text(data=df.n, aes(x=3.6, y = 4.2, label=n), parse=TRUE) +
  ggtitle("Horizontal Gape")
  facet_grid(. ~ cyl, labeller = label_both)

# Plotting Mouth Area 
# Ellipse equation: pi(a/2)(b/2)
fish$area <- with(fish, pi*(gh/2)*(gw/2)) ## adds fish$area to fish data frame

ggplot(data=fish, aes(x=log(SLMM), y=log(area))) +
  geom_point(shape=20) +
  facet_wrap(~ SpeciesCode, ncol = 5) +
  geom_text(data=df.n, aes(x=3.6, y =7.8, label=n), parse=TRUE) +
  ggtitle("Mouth Area")
  facet_grid(. ~ cyl, labeller = label_both)


################### Individual spp dataframes ######################

# Create a set of plots for a given species for each region
# Subsets of each species for regional analysis: 
ac.nigr <- fish[(which(fish$SpeciesCode == "AC.NIGR")),]
ac.oliv <- fish[(which(fish$SpeciesCode == "AC.OLIV")),]
ac.trios <- fish[(which(fish$SpeciesCode == "AC.TRIOS")),]
ap.furc <- fish[(which(fish$SpeciesCode == "AP.FURC")),]
ca.mela <- fish[(which(fish$SpeciesCode == "CA.MELA")),]
ca.tere <- fish[(which(fish$SpeciesCode == "CA.TERE")),]
ce.argu <- fish[(which(fish$SpeciesCode == "CE.ARGU")),]
ce.flav <- fish[(which(fish$SpeciesCode == "CE.FLAV")),]
ce.urod <- fish[(which(fish$SpeciesCode == "CE.UROD")),]
ch.auri <- fish[(which(fish$SpeciesCode == "CH.AURI")),]
ch.orna <- fish[(which(fish$SpeciesCode == "CH.ORNA")),]
ch.sord <- fish[(which(fish$SpeciesCode == "CH.SORD")),]
ch.vand <- fish[(which(fish$SpeciesCode == "CH.VAND")),]
ct.marg <- fish[(which(fish$SpeciesCode == "CT.MARG")),]
lu.boha <- fish[(which(fish$SpeciesCode == "LU.BOHA")),]
lu.kasm <- fish[(which(fish$SpeciesCode == "LU.KASM")),]
me.nige <- fish[(which(fish$SpeciesCode == "ME.NIGE")),]
mo.gran <- fish[(which(fish$SpeciesCode == "MO.GRAN")),]
pa.arca <- fish[(which(fish$SpeciesCode == "PA.ARCA")),]
pa.insu <- fish[(which(fish$SpeciesCode == "PA.INSU")),]
pl.dick <- fish[(which(fish$SpeciesCode == "PL.DICK")),]
ps.bart <- fish[(which(fish$SpeciesCode == "PS.BART")),]
ps.coop <- fish[(which(fish$SpeciesCode == "PS.COOP")),]
ps.disp <- fish[(which(fish$SpeciesCode == "PS.DISP")),]
ps.oliv <- fish[(which(fish$SpeciesCode == "PS.OLIV")),]
pt.tile <- fish[(which(fish$SpeciesCode == "PT.TILE")),]
sc.fren <- fish[(which(fish$SpeciesCode == "SC.FREN")),]
sc.rubr <- fish[(which(fish$SpeciesCode == "SC.RUBR")),]
va.lout <- fish[(which(fish$SpeciesCode == "VA.LOUT")),]

################### Plots with spp facetted by region ####################

# AC.NIGR by region

# n = per region
ac.nigr <- fish[(which(fish$SpeciesCode == "AC.NIGR")),] # AC.NIGR subset of df fish
df.ac.nigr.n <- ddply(.data=ac.nigr, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# AC.NIGR by region faceted plots
ggplot(data=ac.nigr, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("AC.NIGR") +
  geom_text(data=df.ac.nigr, aes(x=4.5, y=2.75, label=n), parse=TRUE)

# n = number of AC.NIGR in total
n.ac.nigr <- paste("n ==", length(ac.nigr$SpeciesCode))
print(n.ac.nigr)

# AC.NIGR with regions colour coded
ggplot(data=ac.nigr, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=20) +
  ggtitle("AC.NIGR") +
  annotate("text", x=4.5, y=2.75, label=n.ac.nigr, parse=TRUE)
            

# AC.OLIV by region
ac.oliv <- fish[(which(fish$SpeciesCode == "AC.OLIV")),] # AC.OLIV subset of fish
# n = per region
df.ac.oliv.n <- ddply(.data=ac.oliv, .(Fishing), summarize, n=paste("n ==", length(Fishing)))


# AC.OLIV by region faceted plots
ggplot(data=ac.oliv, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("AC.OLIV") +
  geom_text(data=df.ac.oliv.n, aes(x=5.1, y=2.75, label=n), parse=TRUE)

# AC.TRIOS by region
ac.trios <- fish[(which(fish$SpeciesCode == "AC.TRIOS")),] # AC.TRIOS subset of fish
# n = per region
df.ac.trios.n <- ddply(.data=ac.trios, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# AC.TRIOS by region faceted plots
ggplot(data=ac.trios, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("AC.TRIOS") +
  geom_text(data=df.ac.trios.n, aes(x=5.1, y=2.75, label=n), parse=TRUE)

# AP.FURC by region
ap.furc <- fish[(which(fish$SpeciesCode == "AP.FURC")),] # AP.FURC subset of df fish
# n = per region
df.ap.furc.n <- ddply(.data=ap.furc, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# AP.FURC by region faceted plots
ggplot(data=ap.furc, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("AP.FURC") +
  geom_text(data=df.ap.furc.n, aes(x=5.65, y=3.8, label=n), parse=TRUE)

# CA.MELA by region
# n = per region
ca.mela <- fish[(which(fish$SpeciesCode == "CA.MELA")),] # CA.MELA subset of df fish
df.ca.mela.n <- ddply(.data=ca.mela, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# CA.MELA by region faceted plots
ggplot(data=ca.mela, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("CA.MELA") +
  geom_text(data=df.ca.mela.n, aes(x=5.4, y=4, label=n), parse=TRUE)

# CA.TERE by region
# n = per region
ca.tere <- fish[(which(fish$SpeciesCode == "CA.TERE")),] # CA.TERE subset of df fish
df.ca.tere.n <- ddply(.data=ca.tere, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# CA.TERE by region faceted plots
ggplot(data=ca.tere, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("CA.TERE") +
  geom_text(data=df.ca.tere.n, aes(x=4.6, y=3.2, label=n), parse=TRUE)

# CE.ARGU by region
# n = per region
ce.argu <- fish[(which(fish$SpeciesCode == "CE.ARGU")),] # CE.ARGU subset of df fish
df.ce.argu.n <- ddply(.data=ce.argu, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# CE.ARGU by region faceted plots
ggplot(data=ce.argu, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("CE.ARGU") +
  geom_text(data=df.ce.argu.n, aes(x=5, y=4.5, label=n), parse=TRUE)

# CE.FLAV by region
# n = per region
ce.flav <- fish[(which(fish$SpeciesCode == "CE.FLAV")),] # CE.FLAV subset of df fish
df.ce.flav.n <- ddply(.data=ce.flav, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# CE.FLAV by region faceted plots
ggplot(data=ce.flav, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("CE.FLAV") +
  geom_text(data=df.ce.flav.n, aes(x=4.25, y=1.3, label=n), parse=TRUE)

# CE.UROD by region
# n = per region
ce.urod <- fish[(which(fish$SpeciesCode == "CE.UROD")),] # CE.UROD subset of df fish
df.ce.urod.n <- ddply(.data=ce.urod, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# CE.UROD by region faceted plots
ggplot(data=ce.urod, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("CE.UROD") +
  geom_text(data=df.ce.urod.n, aes(x=5, y=2.6, label=n), parse=TRUE)

# CH.AURI by region
# n = per region
ch.auri <- fish[(which(fish$SpeciesCode == "CH.AURI")),] # CH.AURI subset of df fish
df.ch.auri.n <- ddply(.data=ch.auri, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# CH.AURI by region faceted plots
ggplot(data=ch.auri, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("CH.AURI") +
  geom_text(data=df.ch.auri.n, aes(x=4.5, y=2.6, label=n), parse=TRUE)

# CH.ORNA by region
# n = per region
ch.orna <- fish[(which(fish$SpeciesCode == "CH.ORNA")),] # CH.ORNA subset of df fish
df.ch.orna.n <- ddply(.data=ch.orna, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# CH.ORNA by region faceted plots
ggplot(data=ch.orna, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("CH.ORNA") +
  geom_text(data=df.ch.orna.n, aes(x=4.57, y=2.8, label=n), parse=TRUE)


# CH.SORD by region
# n = per region
ch.sord <- fish[(which(fish$SpeciesCode == "CH.SORD")),] # CH.SORD subset of df fish
df.ch.sord.n <- ddply(.data=ch.sord, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# CH.SORD by region faceted plots
ggplot(data=ch.sord, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("CH.SORD") +
  geom_text(data=df.ch.sord.n, aes(x=4.5, y=3.8, label=n), parse=TRUE)


# CH.VAND by region
# n = per region
ch.vand <- fish[(which(fish$SpeciesCode == "CH.VAND")),] # CH.VAND subset of df fish
df.ch.vand.n <- ddply(.data=ch.vand, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# CH.VAND by region faceted plots
ggplot(data=ch.vand, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("CH.VAND") +
  geom_text(data=df.ch.vand.n, aes(x=3.6, y=0.6, label=n), parse=TRUE)


# CT.MARG by region
# n = per region
ct.marg <- fish[(which(fish$SpeciesCode == "CT.MARG")),] # CT.MARG subset of df fish
df.ct.marg.n <- ddply(.data=ct.marg, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# CT.MARG by region faceted plots
ggplot(data=ct.marg, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("CT.MARG") +
  geom_text(data=df.ct.marg.n, aes(x=5.25, y=1.5, label=n), parse=TRUE)


# LU.BOHA by region
# n = per region
lu.boha <- fish[(which(fish$SpeciesCode == "LU.BOHA")),] # LU.BOHA subset of df fish
df.lu.boha.n <- ddply(.data=lu.boha, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# LU.BOHA by region faceted plots
ggplot(data=lu.boha, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("LU.BOHA") +
  geom_text(data=df.lu.boha.n, aes(x=3.6, y=4.4, label=n), parse=TRUE)

# LU.KASM by region
lu.kasm <- fish[(which(fish$SpeciesCode == "LU.KASM")),] # LU.KASM subset of fish
# n = per region
df.lu.kasm.n <- ddply(.data=lu.kasm, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# LU.KASM by region faceted plots
ggplot(data=lu.kasm, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("LU.KASM") +
  geom_text(data=df.lu.kasm.n, aes(x=5.1, y=2.75, label=n), parse=TRUE)


# ME.NIGE by region
# n = per region
me.nige <- fish[(which(fish$SpeciesCode == "ME.NIGE")),] # ME.NIGE subset of df fish
df.me.nige.n <- ddply(.data=me.nige, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# ME.NIGE by region faceted plots
ggplot(data=me.nige, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("ME.NIGE") +
  geom_text(data=df.me.nige.n, aes(x=4.75, y=3, label=n), parse=TRUE)



# MO.GRAN by region
# n = per region
mo.gran <- fish[(which(fish$SpeciesCode == "MO.GRAN")),] # MO.GRAN subset of df fish
df.mo.gran.n <- ddply(.data=mo.gran, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# MO.GRAN by region faceted plots
ggplot(data=mo.gran, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("MO.GRAN") +
  geom_text(data=df.mo.gran.n, aes(x=5, y=4, label=n), parse=TRUE)


# PA.ARCA by region
# n = per region
pa.arca <- fish[(which(fish$SpeciesCode == "PA.ARCA")),] # PA.ARCA subset of df fish
df.pa.arca.n <- ddply(.data=pa.arca, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# PA.ARCA by region faceted plots
ggplot(data=pa.arca, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("PA.ARCA") +
  geom_text(data=df.pa.arca.n, aes(x=4.2, y=3, label=n), parse=TRUE)


# PA.INSU by region
# n = per region
pa.insu <- fish[(which(fish$SpeciesCode == "PA.INSU")),] # PA.INSU subset of df fish
df.pa.insu.n <- ddply(.data=pa.insu, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# PA.INSU by region faceted plots
ggplot(data=pa.insu, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("PA.INSU") +
  geom_text(data=df.pa.insu.n, aes(x=4.25, y=3.3, label=n), parse=TRUE)


# PL.DICK by region
# n = per region
pl.dick <- fish[(which(fish$SpeciesCode == "PL.DICK")),] # PL.DICK subset of df fish
df.pl.dick.n <- ddply(.data=pl.dick, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# PL.DICK by region faceted plots
ggplot(data=pl.dick, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("PL.DICK") +
  geom_text(data=df.pl.dick.n, aes(x=3.6, y=2, label=n), parse=TRUE)



# PS.BART by region
# n = per region
ps.bart <- fish[(which(fish$SpeciesCode == "PS.BART")),] # PS.BART subset of df fish
df.ps.bart.n <- ddply(.data=ps.bart, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# PS.BART by region faceted plots
ggplot(data=ps.bart, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("PS.BART") +
  geom_text(data=df.ps.bart.n, aes(x=3.3, y=2.3, label=n), parse=TRUE)



# PS.COOP by region
# n = per region
ps.coop <- fish[(which(fish$SpeciesCode == "PS.COOP")),] # PS.COOP subset of df fish
df.ps.coop.n <- ddply(.data=ps.coop, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# PS.COOP by region faceted plots
ggplot(data=ps.coop, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("PS.COOP") +
  geom_text(data=df.ps.coop.n, aes(x=3.925, y=2.225, label=n), parse=TRUE)


# PS.DISP by region
# n = per region
ps.disp <- fish[(which(fish$SpeciesCode == "PS.DISP")),] # PS.DISP subset of df fish
df.ps.disp.n <- ddply(.data=ps.disp, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# PS.DISP by region faceted plots
ggplot(data=ps.disp, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("PS.DISP") +
  geom_text(data=df.ps.disp.n, aes(x=3.3, y=2.25, label=n), parse=TRUE)


# PS.OLIV by region
# n = per region
ps.oliv <- fish[(which(fish$SpeciesCode == "PS.OLIV")),] # PS.OLIV subset of df fish
df.ps.oliv.n <- ddply(.data=ps.oliv, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# PS.OLIV by region faceted plots
ggplot(data=ps.oliv, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("PS.OLIV") +
  geom_text(data=df.ps.oliv.n, aes(x=3.3, y=2.3, label=n), parse=TRUE)


# PT.TILE by region
# n = per region
pt.tile <- fish[(which(fish$SpeciesCode == "PT.TILE")),] # PT.TILE subset of df fish
df.pt.tile.n <- ddply(.data=pt.tile, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# PT.TILE by region faceted plots
ggplot(data=pt.tile, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("PT.TILE") +
  geom_text(data=df.pt.tile.n, aes(x=3, y=2.9, label=n), parse=TRUE)


# SC.FREN by region
# n = per region
sc.fren <- fish[(which(fish$SpeciesCode == "SC.FREN")),] # SC.FREN subset of df fish
df.sc.fren.n <- ddply(.data=sc.fren, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# SC.FREN by region faceted plots
ggplot(data=sc.fren, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("SC.FREN") +
  geom_text(data=df.sc.fren.n, aes(x=4.85, y=3.5, label=n), parse=TRUE)


# SC.RUBR by region
# n = per region
sc.rubr <- fish[(which(fish$SpeciesCode == "SC.RUBR")),] # SC.RUBR subset of df fish
df.sc.rubr.n <- ddply(.data=sc.rubr, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# SC.RUBR by region faceted plots
ggplot(data=sc.rubr, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("SC.RUBR") +
  geom_text(data=df.sc.rubr.n, aes(x=4.6, y=3.8, label=n), parse=TRUE)



# VA.LOUT by region
# n = per region
va.lout <- fish[(which(fish$SpeciesCode == "VA.LOUT")),] # VA.LOUT subset of df fish
df.va.lout.n <- ddply(.data=va.lout, .(Fishing), summarize, n=paste("n ==", length(Fishing)))

# VA.LOUT by region faceted plots
ggplot(data=va.lout, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ Fishing) +
  ggtitle("VA.LOUT") +
  geom_text(data=df.va.lout.n, aes(x=5.3, y=4.6, label=n), parse=TRUE)






ggplot(dat, aes(x=date, y=value, color=location, group=location)) + 
  geom_line()+
  facet_grid(product ~ ., scale = "free_y")+
  geom_text(data=corr_plot, aes(x=date, y=value, label=label), 
            colour="black", inherit.aes=FALSE, parse=TRUE)


ggplot(fish, aes(x=log(Weight), y=log(GapeHeight))) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm, size=2) +
  opts(axis.title.x = theme_text(size=24, vjust=0, colour="black")) +
  opts(axis.title.y = theme_text(size=24, vjust=0.3, colour="black")) +
  opts(axis.text.x = theme_text(colour="black", size=16)) +
  opts(axis.text.y = theme_text(colour="black", size=16)) +
  opts(legend.text = theme_text(size=20, colour="black")) +
  opts(strip.text.x = theme_text(colour="black", size=16)) +
  opts(legend.title = theme_blank()) +  
  guides(colour = guide_legend(override.aes=list(size=3))) + 
  facet_wrap(~ Guild, nrow=2)


################### Plots of individ. spp coloured by region #####################

# AC.NIGR by region
# n = number of AC.NIGR
n.ac.nigr <- paste("n ==", length(ac.nigr$SpeciesCode))

# AC.NIGR with regions colour coded
ggplot(data=ac.nigr, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("AC.NIGR") +
  annotate("text", x=4.5, y=2.75, label=n.ac.nigr, parse=TRUE)


# AC.OLIV by region
# n = number of AC.OLIV
n.ac.oliv <- paste("n ==", length(ac.oliv$SpeciesCode))

# AC.OLIV with regions colour coded
ggplot(data=ac.oliv, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("AC.OLIV") +
  annotate("text", x=5.1, y=2.85, label=n.ac.oliv, parse=TRUE)

# AC.TRIOS by region
# n = number of AC.TRIOS
n.ac.trios <- paste("n ==", length(ac.trios$SpeciesCode))

# AC.TRIOS with regions colour coded
ggplot(data=ac.trios, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("AC.TRIOS") +
  annotate("text", x=4.8, y=2.45, label=n.ac.trios, parse=TRUE)


# AP.FURCA by region
# n = number of AP.FURCA
n.ap.furca <- paste("n ==", length(ap.furca$SpeciesCode))

# AP.FURCA with regions colour coded
ggplot(data=ap.furca, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("AP.FURCA") +
  annotate("text", x=4.8, y=2.45, label=n.ap.furca, parse=TRUE)



# CA.MELA by region
# n = number of CA.MELA
n.ca.mela <- paste("n ==", length(ca.mela$SpeciesCode))

# CA.MELA with regions colour coded
ggplot(data=ca.mela, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("CA.MELA") +
  annotate("text", x=5.5, y=4, label=n.ca.mela, parse=TRUE)



# CA.TERE by region
# n = number of CA.TERE
n.ca.tere <- paste("n ==", length(ca.tere$SpeciesCode))

# CA.TERE with regions colour coded
ggplot(data=ca.tere, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("CA.TERE") +
  annotate("text", x=4.6, y=3.3, label=n.ca.tere, parse=TRUE)


# CE.ARGU by region
# n = number of CE.ARGU
n.ce.argu <- paste("n ==", length(ce.argu$SpeciesCode))

# CE.ARGU with regions colour coded
ggplot(data=ce.argu, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("CE.ARGU") +
  annotate("text", x=5, y=4.4, label=n.ce.flav, parse=TRUE)


# CE.FLAV by region
# n = number of CE.FLAV
n.ce.flav <- paste("n ==", length(ce.flav$SpeciesCode))

# CE.FLAV with regions colour coded
ggplot(data=ce.flav, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("CE.FLAV") +
  annotate("text", x=3.95, y=1.8, label=n.ce.flav, parse=TRUE)


# CE.UROD by region
# n = number of CE.UROD
n.ce.urod <- paste("n ==", length(ce.urod$SpeciesCode))

# CE.UROD with regions colour coded
ggplot(data=ce.urod, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("CE.UROD") +
  annotate("text", x=4.58, y=3.6, label=n.ce.urod, parse=TRUE)


# CH.AURI by region
# n = number of CH.AURI
n.ch.auri <- paste("n ==", length(ch.auri$SpeciesCode))

# CH.AURI with regions colour coded
ggplot(data=ch.auri, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("CH.AURI") +
  annotate("text", x=4.5, y=2.45, label=n.ch.auri, parse=TRUE)


# CH.ORNA by region
# n = number of CH.ORNA
n.ch.orna <- paste("n ==", length(ch.orna$SpeciesCode))

# CH.ORNA with regions colour coded
ggplot(data=ch.orna, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("CH.ORNA") +
  annotate("text", x=4.6, y=2.75, label=n.ch.orna, parse=TRUE)


# CH.SORD by region
# n = number of CH.SORD
n.ch.sord <- paste("n ==", length(ch.sord$SpeciesCode))

# CH.SORD with regions colour coded
ggplot(data=ch.sord, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("CH.SORD") +
  annotate("text", x=4.5, y=3.5, label=n.ch.sord, parse=TRUE)


# CH.VAND by region
# n = number of CH.VAND
n.ch.vand <- paste("n ==", length(ch.vand$SpeciesCode))

# CH.VAND with regions colour coded
ggplot(data=ch.vand, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("CH.VAND") +
  annotate("text", x=3.1, y=1.4, label=n.ch.vand, parse=TRUE)


# CT.MARG by region
# n = number of CT.MARG
n.ct.marg <- paste("n ==", length(ct.marg$SpeciesCode))

# CT.MARG with regions colour coded
ggplot(data=ct.marg, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("CT.MARG") +
  annotate("text", x=3.9, y=3, label=n.ct.marg, parse=TRUE)


# LU.BOHA by region
# n = number of LU.BOHA
n.lu.boha <- paste("n ==", length(lu.boha$SpeciesCode))

# LU.BOHA with regions colour coded
ggplot(data=lu.boha, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("LU.BOHA") +
  annotate("text", x=4.8, y=4.3, label=n.lu.boha, parse=TRUE)


# LU.KASM by region
# n = number of LU.KASM
n.lu.kasm <- paste("n ==", length(lu.kasm$SpeciesCode))

# LU.KASM with regions colour coded
ggplot(data=lu.kasm, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("LU.KASM") +
  annotate("text", x=3.7, y=3, label=n.lu.kasm, parse=TRUE)


# ME.NIGE by region
# n = number of ME.NIGE
n.me.nige <- paste("n ==", length(me.nige$SpeciesCode))

# ME.NIGE with regions colour coded
ggplot(data=me.nige, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("ME.NIGE") +
  annotate("text", x=4.8, y=2.8, label=n.me.nige, parse=TRUE)


# MO.GRAN by region
# n = number of MO.GRAN
n.mo.gran <- paste("n ==", length(mo.gran$SpeciesCode))

# MO.GRAN with regions colour coded
ggplot(data=mo.gran, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("MO.GRAN") +
  annotate("text", x=5, y=3.8, label=n.mo.gran, parse=TRUE)



# PA.ARCA by region
# n = number of PA.ARCA
n.pa.arca <- paste("n ==", length(pa.arca$SpeciesCode))

# PA.ARCA with regions colour coded
ggplot(data=pa.arca, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("PA.ARCA") +
  annotate("text", x=4.2, y=3, label=n.pa.arca, parse=TRUE)



# PA.INSU by region
# n = number of PA.INSU
n.pa.insu <- paste("n ==", length(pa.insu$SpeciesCode))

# PA.INSU with regions colour coded
ggplot(data=pa.insu, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("PA.INSU") +
  annotate("text", x=4.4, y=3.2, label=n.pa.insu, parse=TRUE)


# PL.DICK by region
# n = number of PL.DICK
n.pl.dick <- paste("n ==", length(pl.dick$SpeciesCode))

# PL.DICK with regions colour coded
ggplot(data=pl.dick, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("PL.DICK") +
  annotate("text", x=3.6, y=2, label=n.pl.dick, parse=TRUE)



# PS.BART by region
# n = number of PS.BART
n.ps.bart <- paste("n ==", length(ps.bart$SpeciesCode))

# PS.BART with regions colour coded
ggplot(data=ps.bart, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("PS.BART") +
  annotate("text", x=3.5, y=2.25, label=n.ps.bart, parse=TRUE)



# PS.COOP by region
# n = number of PS.COOP
n.ps.coop <- paste("n ==", length(ps.coop$SpeciesCode))

# PS.COOP with regions colour coded
ggplot(data=ps.coop, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("PS.COOP") +
  annotate("text", x=3.9, y=2.25, label=n.ps.coop, parse=TRUE)



# PS.DISP by region
# n = number of PS.DISP
n.ps.disp <- paste("n ==", length(ps.disp$SpeciesCode))

# PS.DISP with regions colour coded
ggplot(data=ps.disp, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("PS.DISP") +
  annotate("text", x=3.4, y=2.2, label=n.ps.disp, parse=TRUE)



# PS.OLIV by region
# n = number of PS.OLIV
n.ps.oliv <- paste("n ==", length(ps.oliv$SpeciesCode))

# PS.OLIV with regions colour coded
ggplot(data=ps.oliv, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("PS.OLIV") +
  annotate("text", x=3.3, y=2.4, label=n.ps.oliv, parse=TRUE)



# PT.TILE by region
# n = number of PT.TILE
n.pt.tile <- paste("n ==", length(pt.tile$SpeciesCode))

# PT.TILE with regions colour coded
ggplot(data=pt.tile, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("PT.TILE") +
  annotate("text", x=4.2, y=2.8, label=n.pt.tile, parse=TRUE)



# SC.FREN by region
# n = number of SC.FREN
n.sc.fren <- paste("n ==", length(sc.fren$SpeciesCode))

# SC.FREN with regions colour coded
ggplot(data=sc.fren, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("SC.FREN") +
  annotate("text", x=4.8, y=3.5, label=n.sc.fren, parse=TRUE)


# SC.RUBR by region
# n = number of SC.RUBR
n.sc.rubr <- paste("n ==", length(sc.rubr$SpeciesCode))

# SC.RUBR with regions colour coded
ggplot(data=sc.rubr, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("SC.RUBR") +
  annotate("text", x=4.8, y=3.5, label=n.sc.rubr, parse=TRUE)



# VA.LOUT by region
# n = number of VA.LOUT
n.va.lout <- paste("n ==", length(va.lout$SpeciesCode))

# VA.LOUT with regions colour coded
ggplot(data=va.lout, aes(x=log(SLMM), y=log(gh), colour=Fishing)) +
  geom_point(shape=19) +
  ggtitle("VA.LOUT") +
  annotate("text", x=5.3, y=4.5, label=n.va.lout, parse=TRUE)

################### Modelling each species ####################
# General info: 
  # dlply: takes a data frame, SPLITS up in the same way as ddply, APPLIES
  # function to each piece and COMBINES the results into a list
  
  # ld.ply: takes a list, SPLITS up into elements, APPLIES function to each
  # piece then COMBINES the results into a dataframe


# dlply works like ddply but returns a LIST
# the function(z) is applied to each group of SpeciesCode in the dataframe

#### spp_lms_xxx runs the lm function for each Spp in the dataframe ####
spp_lms_gh <- dlply(fish, .(SpeciesCode), function(z) {
  lm(log(gh)~log(SLMM), data=z)
})

spp_lms_gw <- dlply(fish, .(SpeciesCode), function(z) {
  lm(log(gw)~log(SLMM), data=z)
})

fish$area <- with(fish, pi*(gh/2)*(gw/2)) # just making sure I have another area
# calculation with the relevant data frame
spp_lms_ga <- dlply(fish, .(SpeciesCode), function(z) {
  lm(log(area)~log(SLMM), data=z)
})

##### spp_lms_xxx coefficients dataframes #####
coefs_gh <- ldply(spp_lms_gh, coef)
coefs_gw <- ldply(spp_lms_gw, coef)
coefs_ga <- ldply(spp_lms_ga, coef)


#### Master regression dataframes!!!! #####
####For FUTURE REF: 
  # str() tells you what the internal structure of an R object... 
  # i.e., how the hell to reference R object columns lik:
  # $r.square

##### spp_lms_gh regression dataframes with r^2, SE, p-val ####

#### these pull the rsqd, SE, and p-value for each spp regression (and rename the
  # columns as such)
rsq <- function(spp_lms_gh) c("rsqd"=summary(spp_lms_gh)$r.squared)
se <- function(spp_lms_gh) c("SE"=summary(spp_lms_gh)$coefficients[2,2])
p_val <- function(spp_lms_gh) c("p_val"=summary(spp_lms_gh)$coefficients[2,4])

# This dumps the rsqd, SE, and p-values intoa list of data frames ready to be 
  # join_all'd
dfs <- list(
  rsqs <- ldply(spp_lms_gh, rsq),
  ses <- ldply(spp_lms_gh, se),
  p_vals <- ldply(spp_lms_gh, p_val)
)

# This produces my super data frame!! Containing each species, the rsqd, SE, and
  # p-value of the regression!
gh_summ <- join_all(dfs, by="SpeciesCode", match="first")

#### spp_lms_gw regression dataframes with r^2, SE, p-val####
rsq <- function(spp_lms_gw) c("rsqd"=summary(spp_lms_gw)$r.squared)
se <- function(spp_lms_gw) c("SE"=summary(spp_lms_gw)$coefficients[2,2])
p_val <- function(spp_lms_gw) c("p_val"=summary(spp_lms_gw)$coefficients[2,4])

dfs <- list(
  rsqs <- ldply(spp_lms_gw, rsq),
  ses <- ldply(spp_lms_gw, se),
  p_vals <- ldply(spp_lms_gw, p_val)
)

gw_summ <- join_all(dfs, by="SpeciesCode", match="first")

#### spp_lms_ga regression dataframes with r^2, SE, p-val####
rsq <- function(spp_lms_ga) c("rsqd"=summary(spp_lms_ga)$r.squared)
se <- function(spp_lms_ga) c("SE"=summary(spp_lms_ga)$coefficients[2,2])
p_val <- function(spp_lms_ga) c("p_val"=summary(spp_lms_ga)$coefficients[2,4])

dfs <- list(
  rsqs <- ldply(spp_lms_ga, rsq),
  ses <- ldply(spp_lms_ga, se),
  p_vals <- ldply(spp_lms_ga, p_val)
)

ga_summ <- join_all(dfs, by="SpeciesCode", match="first")

#### spp_lms_gh regression MASTER dataframe (coefs and summ) #####

dfs <- list(
  coefs_gh <- ldply(spp_lms_gh, coef),
  rsqs <- ldply(spp_lms_gh, rsq),
  ses <- ldply(spp_lms_gh, se),
  p_vals <- ldply(spp_lms_gh, p_val)
)

gh_master <- join_all(dfs, by="SpeciesCode", match="first")

#### spp_lms_gw regression MASTER dataframe (coefs and summ) #####

dfs <- list(
  coefs_gw <- ldply(spp_lms_gw, coef),
  rsqs <- ldply(spp_lms_gw, rsq),
  ses <- ldply(spp_lms_gw, se),
  p_vals <- ldply(spp_lms_gw, p_val)
)

gw_master <- join_all(dfs, by="SpeciesCode", match="first")

#### spp_lms_ga regression MASTER dataframe (coefs and summ) #####

dfs <- list(
  coefs_ga <- ldply(spp_lms_ga, coef),
  rsqs <- ldply(spp_lms_ga, rsq),
  ses <- ldply(spp_lms_ga, se),
  p_vals <- ldply(spp_lms_ga, p_val)
)

ga_master <- join_all(dfs, by="SpeciesCode", match="first")


#### Random groupwise modelling stuff from H. Wickam's lecture ####
length(spp_lms_gh)
length(spp_lms_gw)
length(spp_lms_ga)

summary(spp_lms_gh)

head(coef(spp_lms_gh[[1:27]]))

# Just looking at some of the summary outputs from the dlply crap
summary(spp_lms_gh[[1]])
anova(spp_lms_gh[[1]])
coef(spp_lms_gh[[1]])
resid(spp_lms_gh[[1]])

# ldply is the inverse of dlply - it takes a list and combines the ouput into a
  # dataframe
coefs1 <- ldply(spp_lms_gh, coef)
coefs2 <- ldply(spp_lms_gw, coef)
coefs3 <- ldply(spp_lms_ga, coef)

resids <- ldply(spp_lms_gh, resid)

################### Plots for linear regressions for each spp ######################

### Dataframe of unique species
distinct.fish <- as.data.frame(as.factor((unique(fish$SpeciesCode))))

### Plots for each spp. residuals of gh ~ SLMM
pdf() ## fuck yes. I have just printed 29 graphs into a single pdf


names <- t(c("spp", "fg", "int", "slope", "r^2"))
x <- matrix(names, ncol=5, nrow=1)

for (i in distinct.fish[1:29,1]) { # gives list of SpeciesCodes
  a <- fish[(which(fish$SpeciesCode==i)),]
  lm_sp <- with(a, lm(log(gh) ~ log(SLMM)))
  fg_spp <- as.character(unique(fish[(which(fish$SpeciesCode==i)),3]))
  ord <- as.character(unique(fish[(which(fish$SpeciesCode==i)),]))
  x <- rbind(x, t(c(i, fg_spp, coef(lm_sp)[1], coef(lm_sp)[2], summary(lm_sp)$r.squared))) 
  #  code for plotting residual plots of each spp:
  #  resid <- with(a, plot(SLMM, residuals(lm_sp), xlab="SLMM", ylab="Residuals", main=i))
  #  abline(h=0, col='red')
}

x.df <- as.data.frame(x)
ggplot(data=x.df, aes(x=[2:30,4], y=x.df[2:30,2])) +
  geom_point(shape=20)


ggplot(data=x, aes(x=slope, y=fg, colour))

ggplot(data=fish, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~ SpeciesCode, ncol=5) +
  geom_text(data=df.n, aes(x=3.6, y=4.2, label=n), parse=TRUE) + 
  # adds n = on each facet
  ggtitle("Vertical Gape") 
facet_grid(. ~ cyl, labeller = label_both)

dev.off() # don't forget to turn off pdf() so that I can see plots in RStudio again

### Attempt to put regression values into matrix....
regress=matrix(NA, ncol=3, nrow=29)

for (r in distinct.fish[1:29,1]) {
    regress[r,1] = r                 # SpeciesCode
#    regress[r,2] = coef(lm_sp)[1]     # Intercept
#    regress[r,3] = coef(lm_sp)[2]     # Slope
  }
}

write.table(x,file="regress.csv",sep=",")   #### writes x into csv


lm_sp <- with(ac.nigr, lm(gh ~ SLMM))


  m = lm(log(gh) ~)


lm_eqn3 = function(x){
  m = lm(log(GapeWidth) ~ log(SL), x)
  eqn <- substitute(italic(y) == a%.% italic(x) + b*","~~italic(r)^2~"="~r2,
                    list(a = format(coef(m)[2], digits = 3),
                         b = format(coef(m)[1], digits = 3),
                         r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eqn))
}



########## Learning how to use ggplot2 ############
library("ggplot2")


lm_eqn = function(df){
  m = lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

p1 = p + geom_text(aes(x = 25, y = 300, label = lm_eqn(df)), parse = TRUE)


as

eq <- italic(y) == a %.% italic(x)^b
as.expression(y == ax^b)

# Creates hypothesis graph.
ggplot(data.frame(x=c(0,5), y=c(0,5)), aes(x, size=2)) +
  opts(axis.text.x = theme_text(size=16, colour="black")) +
  opts(axis.text.y = theme_text(size=16, colour="black")) +
  xlab("log(Standard Length)") +
  ylab("log(Gape Height)") +
#  annotate("text", x=1, y = 2.5,label="y == ax^b", parse=T, size=18) +
  opts(axis.title.x = theme_text(size=28, vjust = 0.2, hjust=0.2)) +
  opts(axis.title.y = theme_text(size=28, vjust = 0.3, hjust=0.2)) +
  #stat_function(fun=function(x) {0.7*x + 0.3 }, geom="line", aes(colour="red")) +
  stat_function(fun=function(x) { x }, geom="line", aes(colour="black")) +
#  stat_function(fun=function(x) {0.9*x }, geom="line", aes(colour="green")) +
  xlim(0,5) +
  ylim(0,5) +
  theme(legend.position="none") +
  scale_x_discrete(breaks=NULL) +
  scale_y_discrete(breaks=NULL) +
  theme(plot.margin = unit(c(3,1,1,2), "cm"))

scale_x_discrete(breaks=NULL)

  scale_colour_manual("Function", value=c("blue","red"), breaks=c("square","exp"))

ggplot(data=fishSL, aes(x=log(SL), y=log(GapeHeight), colour=Guild)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  scale_colour_hue(l=70, c=175) +
  annotate("text", x=3.5, y=4.5, label=lm_eqn2(fishSL), parse=TRUE)

ggplot(data=fishSL, aes(x=SL, y=GapeHeight, colour=Guild)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  scale_colour_hue(l=70, c=175) #+
#  annotate("text", x=3.5, y=4.5, label=lm_eqn2(fishSL), parse=TRUE)


ggplot(data=fish, aes(x=log(Weight), y=log(GapeHeight), colour=Guild)) +
  geom_point(shape=19) +
  opts(axis.title.x = theme_text(size=24, vjust=0, colour="black")) +
  opts(axis.title.y = theme_text(size=24, vjust=0.3, colour="black")) +
  opts(axis.text.x = theme_text(size=16, colour="black")) +
  opts(axis.text.y = theme_text(size=16, colour="black")) +
  opts(legend.text = theme_text(size=20, colour="black")) +
  opts(legend.title = theme_blank()) +
  guides(colour = guide_legend(override.aes=list(shape=20, size=7)))

ggplot(data=fish, aes(x=Weight, y=GapeHeight, colour=Guild)) +
  geom_point(shape=19) +
  opts(axis.title.x = theme_text(size=24, vjust=0, colour="black")) +
  opts(axis.title.y = theme_text(size=24, vjust=0.3, colour="black")) +
  opts(axis.text.x = theme_text(size=16, colour="black")) +
  opts(axis.text.y = theme_text(size=16, colour="black")) +
  opts(legend.text = theme_text(size=20, colour="black")) +
  opts(legend.title = theme_blank()) +
  guides(colour = guide_legend(override.aes=list(shape=20, size=7)))


ggplot(data=fishSL, aes(x=log(SL), y=log(GapeHeight), colour=Guild)) +
  geom_point(shape=19) +
  opts(axis.title.x = theme_text(size=24, vjust=0, colour="black")) +
  opts(axis.title.y = theme_text(size=24, vjust=0.3, colour="black")) +
  opts(axis.text.x = theme_text(size=16, colour="black")) +
  opts(axis.text.y = theme_text(size=16, colour="black")) +
  opts(legend.text = theme_text(size=20, colour="black")) +
  opts(legend.title = theme_blank()) +
  guides(colour = guide_legend(override.aes=list(shape=20, size=7)))


opts(legend.text = theme_text(colour = 'red', angle = 45, size = 10, hjust = 3, vjust = 3, face = 'bold')
  
  #geom_smooth(method=lm) +
  scale_colour_hue(l=70, c=175) #+
#  annotate("text", x=1.5, y=4.3, label=lm_eqn(fish), parse=TRUE)

ggplot(data=fish, aes(x=log(Weight), y=log(GapeWidth), colour=Guild)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  scale_colour_hue(l=70, c=175) +
  annotate("text", x=1.5, y=4.3, label=lm_eqn3(fishGW), parse=TRUE)

ggplot(data=fish, aes(x=log(Weight), y=log(GapeHeight), colour=Guild)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  scale_colour_hue(l=70, c=175) +
  geom_abline(intercept = 0, slope = 1, lty=2)


ggplot(fish, aes(x=log(Weight), y=log(GapeHeight))) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method="lm") 
  
ggplot(fish, aes(x=log(Weight), y=log(GapeHeight))) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm, size=2) +
  opts(axis.title.x = theme_text(size=24, vjust=0, colour="black")) +
  opts(axis.title.y = theme_text(size=24, vjust=0.3, colour="black")) +
  opts(axis.text.x = theme_text(colour="black", size=16)) +
  opts(axis.text.y = theme_text(colour="black", size=16)) +
  opts(legend.text = theme_text(size=20, colour="black")) +
  opts(strip.text.x = theme_text(colour="black", size=16)) +
  opts(legend.title = theme_blank()) +  
  guides(colour = guide_legend(override.aes=list(size=3))) + 
  facet_wrap(~ Guild, nrow=2)
     
  opts(strip.text.x = theme_text(colour = 'red', angle = 45, size = 10, hjust = 0.5, vjust = 0.5, face = 'bold'))
  
  facet_wrap(~ Guild, nrow=2)

lm(log(GapeHeight) ~ log(Weight), fish)


####### Linear equations for labelling on plots above ########
# The following function is from: http://stackoverflow.com/questions/7549694/ggplot2-adding-regression-line-equation-and-r2-on-graph
# It generates the regression equation that can be used to annotate a ggplot.
lm_eqn = function(x){
  m = lm(log(gh) ~ log(SLMM), x)
  eqn <- substitute(italic(y) == a%.% italic(x) + b*","~~italic(r)^2~"="~r2,
                    list(a = format(coef(m)[2], digits = 3),
                         b = format(coef(m)[1], digits = 3),
                         r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eqn))
}

print(lm_eqn(fish))

     
which(fish$SL == 0)
fishSL <- fish[c(-321, -344),]

lm_eqn2 = function(x){
  m = lm(log(GapeHeight) ~ log(SL), x)
  eqn <- substitute(italic(y) == a%.% italic(x) + b*","~~italic(r)^2~"="~r2,
                    list(a = format(coef(m)[2], digits = 3),
                         b = format(coef(m)[1], digits = 3),
                         r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eqn))
}

which(is.na(fish$GapeWidth)==TRUE)
fishGW <- fish[-884,]

lm_eqn3 = function(x){
  m = lm(log(GapeWidth) ~ log(SL), x)
  eqn <- substitute(italic(y) == a%.% italic(x) + b*","~~italic(r)^2~"="~r2,
                    list(a = format(coef(m)[2], digits = 3),
                         b = format(coef(m)[1], digits = 3),
                         r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eqn))
}


qplot(fish$Guild, data=fish, geom='bar', colour=I("black"), fill=SpeciesCode)
qplot(fish$Guild, data=fish, geom='bar', colour=I("black"), fill=SizeClass) # is a histogram of fish frequency by guild, with colours representing the number of samples we have in given size classes

qplot(log(Weight), log(GapeHeight), data=fish, colour=SizeClass)

qplot(log(Weight), log(GapeHeight/Weight), data=fish, colour=Guild)

qplot(log(Weight), log(GapeHeight), data=fish, colour=Guild)
qplot(log(Weight), log(GapeHeight), data=fish, colour=SpeciesCode)

qplot(aes(x=log(Weight), y=log(GapeHeight) #+
  geom_point(shape=1)# Use hollow circles
  geom_smooth(method=lm)


qplot(log(Weight), log(GapeHeight), data=fish, colour=SizeClass)

qplot(log(Weight), log(GapeHeight), data=fish, colour=SpeciesCode)
