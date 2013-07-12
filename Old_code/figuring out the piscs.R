### Figuring out what's going on with the piscivores

# set working directory
setwd('/Users/jillian/R_projects/Allometry/')

# Data for all specimens where a species has more than 3 samples with gape and SL data
fish <- read.csv('Allometry Data/fish_gape_cleaned_more_than_3_jan_24.csv', 
                 header = TRUE, sep = ",") 

########### Entering and Cleaning Data ##########

# Remove PS.COOP and AC.TRIOS because low n, but more so because there is no size gradient
# Which rows have PS.COOP and AC.TRIOS?
psc <- which(fish$SpeciesCode==('PS.COOP'))
act <- which(fish$SpeciesCode==('AC.TRIOS')) 
# removing these rows from fish:
fish <- fish[-c(psc, act),]


##### Selecting only piscivores #######

p <- fish[which(fish$j_fg=='Pi'),]

p$SpeciesCode <- factor(p$SpeciesCode)
p$Family <- factor(p$Family)

##### Ordering by Family to make my life easier ####
p$SpeciesCode <- factor(p$SpeciesCode, levels = c("CA.MELA", 
                                                  "AP.FURC", 
                                                  "LU.BOHA",
                                                  "LU.KASM",
                                                  "CE.ARGU",
                                                  "CE.UROD",
                                                  "VA.LOUT")
)


#### Making some plots!!! ######

library("ggplot2")
library('plyr')
library('gridExtra') # allows me to make multiplots

#### calculating mouth area (and adding it to pisc dataframe): ####
p$area <- with(p, pi*(gh/2)*(gw/2))

#### gh, gw, ga for all piscivore species combined (WITH EQNS!!!!!!) ####

caranx <- fish[which(fish$SpeciesCode=='CA.MELA'),]

df.n <- ddply(.data=p, .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))

lm_eqn = function(x, y){
  m = lm(log(y) ~ log(SLMM), x);
  l <- list(a = format(coef(m)[1], digits = 2), 
            b = format(coef(m)[2], digits = 2), 
            r2 = format(summary(m)$r.squared, digits = 3)
  )
  if (l$a >= 0) {
  eq <- substitute(italic(y) == b %.% italic(x) + a*","~~italic(r)^2~"="~r2, l) 
  }
  
  else {
    l <- list(a = format(abs(coef(m)[1]), digits = 2), 
              b = format(coef(m)[2], digits = 2), 
              r2 = format(summary(m)$r.squared, digits = 3)
    )
    eq <- substitute(italic(y) == b %.% italic(x) - a*","~~italic(r)^2~"="~r2, l) 
  }

  as.character(as.expression(eq));                 
}

# LOG TRANSFORMED
caranx_l <-
  ggplot(data=caranx, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=3.3, y=4, label=lm_eqn(caranx, caranx$gh)), parse=TRUE) +
  ggtitle("CA.MELA") +
  ylim(0, 4.5) +
  xlim(0, 7)

# NOT TRANSFORMED
lm_eqn_notlog = function(x, y){
  m = lm(y ~ SLMM, x);
  l <- list(a = format(coef(m)[1], digits = 2), 
            b = format(coef(m)[2], digits = 2), 
            r2 = format(summary(m)$r.squared, digits = 3)
  )
  if (l$a >= 0) {
    eq <- substitute(italic(y) == b %.% italic(x) + a*","~~italic(r)^2~"="~r2, l) 
  }
  
  else {
    l <- list(a = format(abs(coef(m)[1]), digits = 2), 
              b = format(coef(m)[2], digits = 2), 
              r2 = format(summary(m)$r.squared, digits = 3)
    )
    eq <- substitute(italic(y) == b %.% italic(x) - a*","~~italic(r)^2~"="~r2, l) 
  }
  
  as.character(as.expression(eq));                 
}

caranx_nl <-
  ggplot(data=caranx, aes(x=SLMM, y=gh)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=250, y=65, label=lm_eqn_notlog(caranx, caranx$gh)), parse=TRUE) +
  ggtitle("CA.MELA") +
  ylim(0,100) +
  xlim(0,620)

grid.arrange(caranx_l, caranx_nl, ncol=2)


## AP.FURC
aphareus <- fish[which(fish$SpeciesCode=='AP.FURC'),]

# LOG TRANSFORMED
aphar_l <-
  ggplot(data=aphareus, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=2.5, y=4.3, label=lm_eqn(aphareus, aphareus$gh)), parse=TRUE) +
  ggtitle("AP.FURC") +
  ylim(0,4.4) +
  xlim(0,5.8)

# NOT TRANSFORMED
aphar_nl <-
  ggplot(data=aphareus, aes(x=SLMM, y=gh)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=250, y=85, label=lm_eqn_notlog(aphareus, aphareus$gh)), parse=TRUE) +
  ggtitle("AP.FURC") +
  ylim(0,100) +
  xlim(0,620)

grid.arrange(aphar_l, aphar_nl, ncol=2)


## LU.BOHA
bohar <- fish[which(fish$SpeciesCode=='LU.BOHA'),]

# LOG TRANSFORMED
bohar_l <-
  ggplot(data=bohar, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=2.5, y=4.3, label=lm_eqn(bohar, bohar$gh)), parse=TRUE) +
  ggtitle("LU.BOHA") +
  ylim(0,4.6) +
  xlim(0,6.5)

# NOT TRANSFORMED
bohar_nl <-
  ggplot(data=bohar, aes(x=SLMM, y=gh)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=250, y=95, label=lm_eqn_notlog(bohar, bohar$gh)), parse=TRUE) +
  ggtitle("LU.BOHA") +
  ylim(0,100) +
  xlim(0,550)

grid.arrange(bohar_l, bohar_nl, ncol=2)


## LU.KASM
kasmira <- fish[which(fish$SpeciesCode=='LU.KASM'),]

# LOG TRANSFORMED
kasmira_l <-
  ggplot(data=kasmira, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=2.5, y=3.5, label=lm_eqn(kasmira, kasmira$gh)), parse=TRUE) +
  ggtitle("LU.KASM") +
  ylim(0,4) +
  xlim(0,5.5)

# NOT TRANSFORMED
kasmira_nl <-
  ggplot(data=kasmira, aes(x=SLMM, y=gh)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=75, y=38, label=lm_eqn_notlog(kasmira, kasmira$gh)), parse=TRUE) +
  ggtitle("LU.KASM") +
  ylim(0,40) +
  xlim(0,175)

grid.arrange(kasmira_l, kasmira_nl, ncol=2)



## CE.ARGU
argus <- fish[which(fish$SpeciesCode=='CE.ARGU'),]

# LOG TRANSFORMED
argus_l <-
  ggplot(data=argus, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=2.5, y=4.5, label=lm_eqn(argus, argus$gh)), parse=TRUE) +
  ggtitle("CE.ARGU") +
  ylim(0,5) +
  xlim(0,6)

# NOT TRANSFORMED
argus_nl <-
  ggplot(data=argus, aes(x=SLMM, y=gh)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=150, y=85, label=lm_eqn_notlog(argus, argus$gh)), parse=TRUE) +
  ggtitle("CE.ARGU") +
  ylim(0,100) +
  xlim(0,400)


grid.arrange(argus_l, argus_nl, ncol=2)


## CE.UROD
urod <- fish[which(fish$SpeciesCode=='CE.UROD'),]

# LOG TRANSFORMED
urod_l <-
  ggplot(data=urod, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=2.5, y=3.5, label=lm_eqn(urod, urod$gh)), parse=TRUE) +
  ggtitle("CE.UROD") +
  ylim(0,4.1) +
  xlim(0,5.2)

# NOT TRANSFORMED
urod_nl <-
  ggplot(data=urod, aes(x=SLMM, y=gh)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=60, y=45, label=lm_eqn_notlog(urod, urod$gh)), parse=TRUE) +
  ggtitle("CE.UROD") +
  ylim(0,52) +
  xlim(0,200)

grid.arrange(urod_l, urod_nl, ncol=2)

## VA.LOUT
variola <- fish[which(fish$SpeciesCode=='VA.LOUT'),]

# LOG TRANSFORMED
variola_l <-
  ggplot(data=urod, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=2.5, y=3.5, label=lm_eqn(variola, variola$gh)), parse=TRUE) +
  ggtitle("VA.LOUT") +
  ylim(0,4.1) +
  xlim(0,5.2)

# NOT TRANSFORMED
variola_nl <-
  ggplot(data=variola, aes(x=SLMM, y=gh)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=200, y=115, label=lm_eqn_notlog(variola, variola$gh)), parse=TRUE) +
  ggtitle("VA.LOUT") +
  ylim(0,135) +
  xlim(0,620)

grid.arrange(variola_l, variola_nl, ncol=2)

pdf(file="log_non_log_pisc.pdf", width=10, height=6)
grid.arrange(caranx_l, caranx_nl, ncol=2)
grid.arrange(aphar_l, aphar_nl, ncol=2)
grid.arrange(bohar_l, bohar_nl, ncol=2)
grid.arrange(kasmira_l, kasmira_nl, ncol=2)
grid.arrange(argus_l, argus_nl, ncol=2)
grid.arrange(urod_l, urod_nl, ncol=2)
grid.arrange(variola_l, variola_nl, ncol=2)
dev.off()


#### Making relative gape height, width, area? plots for all
#### of the piscivores ####

# add relative gape height, width, and area to the dataframe p

p$gh_ratio <- p$gh/p$SLMM
p$gw_ratio <- p$gw/p$SLMM
p$ga_ratio <- p$area/p$SLMM

rel_gh <- dlply(p, .(SpeciesCode), function(z) {
  lm(gh_ratio/SLMM~SLMM, data=p)
})

ggplot(data=p, aes(x=SLMM, y=gh_ratio), colour=Family) + 
  geom_point(shape=20) +
  geom_smooth(method=lm)
  


spp1 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3, drop=F) +
  geom_smooth(method=lm) +
  geom_text(data=gh_spp, aes(x=4.85, y=4.1, label=lm_eq), parse=T) +
  ggtitle("VerticalGape")
facet_grid(. ~ cyl, labeller = label_both)














