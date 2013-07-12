# set working directory
setwd('/Users/jillian/R_projects/Allometry/')

# Data for all specimens where a species has more than 3 samples with gape and SL data
fish <- read.csv('Allometry Data/fish_gape_cleaned_more_than_3_jan_24.csv', 
                 header = TRUE, sep = ",") 

library('ggplot2')
library('plyr')
library('gridExtra')

########### Entering and Cleaning Data ##########

# Remove PS.COOP and AC.TRIOS because low n, but more so because there is no size gradient
# Which rows have PS.COOP and AC.TRIOS?
psc <- which(fish$SpeciesCode==('PS.COOP'))
act <- which(fish$SpeciesCode==('AC.TRIOS')) 
# removing these rows from fish:
fish <- fish[-c(psc, act),]

# Setting up dataframe:
fish$area <- with(fish, pi*(gh/2)*(gw/2))

fish$gh_ratio <- fish$gh/fish$SLMM
fish$gw_ratio <- fish$gw/fish$SLMM
fish$ga_ratio <- fish$area/fish$SLMM
##### Functional group level ANOVA ####

# picking just Pi, He, BI, and ZP and ordering them:
pento <- fish[fish$j_fg %in% c('Pi', 'He', 'BI', 'ZP') | 
               fish$SpeciesCode %in% c("CH.ORNA", "CH.AURI"),]


# changing j_fg to C for CH.AURI
pento <- within(fish, j_fg[SpeciesCode=='CH.AURI'] <- 'C')

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
                                                        "CH.AURI",
                                                        "CH.ORNA")
)

# finding some averages for each functional group:

ddply(pento,~j_fg,summarise, mean=mean(gh_ratio), sd=sd(gh_ratio))
ddply(pento,~j_fg,summarise, mean=mean(gw_ratio), sd=sd(gw_ratio))

ddply(pento,~SpeciesCode,summarise,mean=mean(gh_ratio), sd=sd(gh_ratio))
ddply(pento,~SpeciesCode,summarise,mean=mean(gw_ratio), sd=sd(gw_ratio))


# Setting levels to functional groups: 

plot(gh_ratio~j_fg, pento)
plot(gh_ratio~Family, pento)

TukeyHSD(aov(gh_ratio~Family, pento))

# One-way ANOVAs
#for gh_ratio and functional group
fitgh <- aov(gh_ratio~j_fg, pento)
summary(fitgh)

# for gw_ratio and functional group
fitgw <- aov(gw_ratio~j_fg, pento)
summary(fitgw)

# for ga_ratio and functional group
fitga <- aov(ga_ratio~j_fg, pento)
summary(fitga)

# One-way ANOVAs
#for gh_ratio and family
fitgh <- aov(gh_ratio~Family, trio)
summary(fitgh)

# for gw_ratio and family
fitgw <- aov(gw_ratio~Family, trio)
summary(fitgw)

# for ga_ratio and family
fitga <- aov(ga_ratio~Family, trio)
summary(fitga)

# One-way ANOVAs
#for gh_ratio and species
fitgh <- aov(gh_ratio~SpeciesCode, trio)
summary(fitgh)

# for gw_ratio and species
fitgw <- aov(gw_ratio~SpeciesCode, trio)
summary(fitgw)

# for ga_ratio and species
fitga <- aov(ga_ratio~SpeciesCode, trio)
summary(fitga)


ggplot(trio, aes(x=j_fg, y=gh_ratio, fill=j_fg)) +
  geom_boxplot()

ggplot(trio, aes(x=Family, y=gh_ratio, fill=j_fg)) +
  geom_boxplot() +
  scale_fill_manual(values=cbPalette) +
  facet_wrap(~j_fg, scale="free_x")

ggplot(trio, aes(x=Family, y=gw_ratio, fill=j_fg)) +
  geom_boxplot() +
  scale_fill_manual(values=cbPalette) +
  facet_wrap(~j_fg, scale="free_x")

ggplot(trio, aes(x=Family, y=gw_ratio, fill=j_fg)) +
  geom_boxplot() +
  scale_fill_manual(values=cbPalette) +
  facet_wrap(~j_fg, scale="free_x")



cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbPalette <- c("#D55E00", "#0072B2", "#009E73", "#999999", "#E69F00", "#56B4E9", "#F0E442", "#CC79A7")


ggplot(trio, aes(x=SpeciesCode, y=gh_ratio, fill=j_fg)) +
  geom_boxplot() +
  scale_fill_manual(values=cbPalette) +
  facet_grid(~Family, scale="free_x", space="free")



ggplot(trio, aes(x=SpeciesCode, y=gw_ratio, fill=j_fg)) +
  geom_boxplot() +
  scale_fill_manual(values=cbPalette) +
  facet_grid(~Family, scale="free_x", space="free")


ggplot(trio, aes(x=SpeciesCode, y=gh_ratio, fill=j_fg)) +
  geom_boxplot() +
  scale_fill_manual(values=cbPalette) +
  xlab("Species") +
  ylab("Average Relative Gape Height") +
  theme(axis.text.x = element_text(angle=45, size=24, colour="black", vjust=0.5)) +
  theme(axis.text.y = element_text(size=24, colour="black")) +
  theme(axis.title.x = element_text(size=33, vjust = -2, hjust=0.53)) +
  theme(axis.title.y = element_text(size=33, vjust = 0)) +
  facet_grid(~Family, scale="free_x", space="free", labeller=family_labeller) +
  theme(strip.text.x = element_text(angle=50, size=25)) +
  theme(legend.key.size = unit(2, "cm")) +
  scale_fill_manual(values=c("#D55E00", "#0072B2", "#009E73"),
                    #name="Functional\nGroup",
                    labels=c("Piscivore", "Benthic\nInvertivore", "Herbivore")
  ) +
  labs(fill="") +
  theme(legend.text = theme_text(size=24)) +
  theme(plot.margin = unit(c(3,1,1,2), "cm")) 



guide_legend(title = NULL)

scale_fill_discrete(name="Experimental\nCondition")

plot + opts(legend.title = theme_text(colour = 'red', angle = 45, size = 10, hjust = 3, vjust = 7, face = 'italic'))


plot + opts(legend.key.size = unit(2, "cm"))

opts(strip.text.x = theme_text(colour = 'red', angle = 45, size = 10, hjust = 0.5, vjust = 0.5, face = 'bold'))

axis.text.x  = element_text(angle=90, vjust=0.5, size=16)
stat_boxplot(geom="boxplot", position = "dodge", width = 0.60, na.rm = TRUE) +  facet_grid(.~ind




# Creates hypothesis graph.
ggplot(data.frame(x=c(0,5), y=c(0,5)), aes(x, size=2)) +                                                                     
 theme(axis.text.x = theme_text(size=16, colour="black")) +                                                                                 
 theme(axis.text.y = theme_text(size=16, colour="black")) +
 xlab("log(Weight)") +
 ylab("log(GapeHeight)") +
 theme(axis.title.x = theme_text(size=28, vjust = 0)) +
 theme(axis.title.y = theme_text(size=28, vjust = 0.3)) +
 #stat_function(fun=function(x) {0.7*x + 0.3 }, geom="line", aes(colour="red")) +
 stat_function(fun=function(x) {0.9*x }, geom="line", aes(colour="blue")) +
 stat_function(fun=function(x) {0.2*x }, geom="line", aes(colour="green")) +
 xlim(0,5) +
 ylim(0,5) +
 theme(legend.position="none")


family_names <- list(
 'Carangidae'="Jacks",
 'Lutjanidae'="Snappers",
 'Serranidae'="Groupers",
 'Cirrhitidae'="Hawkfishes",
 'Lethrinidae'="Emperors",
 'Mullidae'="Goatfishes",
 'Acanthuridae'="Surgeonfishes",
 'Pomacentridae'="Angelfishes",
 'Scaridae'="Parrotfishes"
)

family_labeller <- function(variable,value) {
 return(family_names[value])
}  

hospital_names <- list(
 'Hospital#1'="Some Hospital",
 'Hospital#2'="Another Hospital"
)


hospital_labeller <- function(variable,value){
 return(hospital_names[value])
}

ggplot(survey,aes(x=age)) + stat_bin(aes(n=nrow(h3),y=..count../n), binwidth=10) + 
 facet_grid(hospital ~ ., labeller=hospital_labeller)

+ facet_grid(hospital ~ ., labeller=family_labeller)

