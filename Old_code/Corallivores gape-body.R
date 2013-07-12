###### Corallivore gape-body size analysis ######

################        Corallivore analysis             ##############

# Multipanel plot with log(gw), log(gh), log(MA) against log(SL)

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


##### Selecting only corallivores #######
c <- fish[fish$SpeciesCode %in% c("CH.ORNA", "CH.AURI"),]


c$SpeciesCode <- factor(c$SpeciesCode)
c$Family <- factor(c$Family)

##### Ordering by Family to make my life easier ####
c$SpeciesCode <- factor(c$SpeciesCode, levels = c("CH.AURI",
                                                  "CH.ORNA")
)

#### Making some plots!!! ######

library("ggplot2")
library('plyr')
library('gridExtra') # allows me to make multiplots

#### calculating mouth area (and adding it to corallivore dataframe): ####
c$area <- with(c, pi*(gh/2)*(gw/2))

#### gh, gw, ga for all corallivore species combined (WITH EQNS!!!!!!) ####

df.n <- ddply(.data=c, .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))

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

c1 <-
  ggplot(data=c, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.6, y=2.8, label=lm_eqn(c, c$gh)), parse=TRUE) +
  ggtitle("Vertical Gape")

c2 <-
  ggplot(data=c, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.6, y=2.8, label=lm_eqn(c, c$gw)), parse=TRUE) +
  ggtitle("Horizontal Gape")

c3 <-
  ggplot(data=c, aes(x=log(SLMM), y=log(area))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.6, y=5.35, label=lm_eqn(c, c$area)), parse=TRUE) +
  ggtitle("Gape Area")

grid.arrange(c1, c2, c3)

#### gh, gw, ga for all corallivore species colour-coded ####
# (has regression lines in colour for each species)
cc1 <-
  ggplot(data=c, aes(x=log(SLMM), y=log(gh), colour=SpeciesCode)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  ggtitle("Vertical Gape")

cc2 <-
  ggplot(data=c, aes(x=log(SLMM), y=log(gw), colour=SpeciesCode)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  ggtitle("Horizontal Gape")

cc3 <-
  ggplot(data=c, aes(x=log(SLMM), y=log(area), colour=SpeciesCode)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  ggtitle("Gape Area")

grid.arrange(cc1, cc2, cc3)

#### gh, gw, ga for all corallivores for EACH SPP (facets) #####
# (SET OF 3 SECTIONS)
#### (1) gh plots with regression lines and equations ####
## gh lm_eqns for plotting:
spp_lms_gh <- dlply(c, .(SpeciesCode), function(z) {
  lm(log(gh)~log(SLMM), data=z)
})

spp.n <- ddply(.data=c, .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))

## pieces for gh_summ:
coefs_gh <- ldply(spp_lms_gh, coef)
rsq <- function(spp_lms_gh) c("rsqd"=summary(spp_lms_gh)$r.squared)
se <- function(spp_lms_gh) c("SE"=summary(spp_lms_gh)$coefficients[2,2])
p_val <- function(spp_lms_gh) c("p_val"=summary(spp_lms_gh)$coefficients[2,4])
spp.n <- ddply(.data=c, .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))

# master df for gh lms 
dfs <- list(
  coefs_gh,
  rsqs <- ldply(spp_lms_gh, rsq),
  ses <- ldply(spp_lms_gh, se),
  p_vals <- ldply(spp_lms_gh, p_val),
  spp.n
)

gh_summ <- join_all(dfs, by="SpeciesCode")

## generating lm_eqns to put onto graphs
m <- gh_summ
gh_spp = matrix(data=NA, nrow=0, ncol=2)
for (i in (1:2)) {
  l <- list(a = format(m[i, 3], digits=2),
            b = format(m[i, 2], digits=2), 
            r2 = format(m[i, 4], digits=2)
  )
  if (l$b >= 0) {
    eqn <- substitute(italic(y) ==a%.% italic(x) + b*~~italic(r)^2~"="~r2, l)
  }
  
  else {
    l <- list(a = format(m[i, 3], digits=2),
              b = format(abs(m[i, 2]), digits=2), 
              r2 = format(m[i, 4], digits=2)
    )
    eqn <- substitute(italic(y) ==a%.% italic(x) - b*~~italic(r)^2~"="~r2, l)
  }
  lm_eq = as.character(as.expression(eqn))
  gh_spp <- rbind(gh_spp, c(as.character(m[i,1]), lm_eq))
}



# gh_spp is the data frame containing the spp name (for ref) and the lm_eqn for
# each regression
gh_spp <- as.data.frame(gh_spp)
colnames(gh_spp) <- c("SpeciesCode", "eqn")

spp1 <-
  ggplot(data=c, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3, drop=F) +
  geom_smooth(method=lm) +
  geom_text(data=gh_spp, aes(x=4.6, y=2.8, label=eqn), parse=T) +
  ggtitle("VerticalGape")
facet_grid(. ~ cyl, labeller = label_both)

spp1
#### (2) gw plots with regression lines and equations ####
## gw lm_eqns for plotting:
spp_lms_gw <- dlply(c, .(SpeciesCode), function(z) {
  lm(log(gw)~log(SLMM), data=z)
})
## pieces for gw_summ:
coefs_gw <- ldply(spp_lms_gw, coef)
rsq <- function(spp_lms_gw) c("rsqd"=summary(spp_lms_gw)$r.squared)
se <- function(spp_lms_gw) c("SE"=summary(spp_lms_gw)$coefficients[2,2])
p_val <- function(spp_lms_gw) c("p_val"=summary(spp_lms_gw)$coefficients[2,4])

# master df for gw lms 
dfs <- list(
  coefs_gw,
  rsqs <- ldply(spp_lms_gw, rsq),
  ses <- ldply(spp_lms_gw, se),
  p_vals <- ldply(spp_lms_gw, p_val)
)

gw_summ <- join_all(dfs, by="SpeciesCode", match="first")

## generating lm_eqns to put onto graphs
m <- gw_summ
gw_spp = matrix(data=NA, nrow=0, ncol=2)
for (i in (1:2)) {
  l <- list(a = format(m[i, 3], digits=2),
            b = format(m[i, 2], digits=2), 
            r2 = format(m[i, 4], digits=2)
  )
  if (l$b >= 0) {
    eqn <- substitute(italic(y) ==a%.% italic(x) + b*","~~italic(r^2)~"="~r2,l)
  }
  
  else {
    l <- list(a = format(m[i, 3], digits=2),
              b = format(abs(m[i, 2]), digits=2), 
              r2 = format(m[i, 4], digits=2)
    )
    eqn <- substitute(italic(y) ==a%.% italic(x) - b*","~~italic(r^2)~"="~r2,l)
  }
  lm_eq = as.character(as.expression(eqn))
  gw_spp <- rbind(gw_spp, c(as.character(m[i,1]), lm_eq))
}
# gw_spp is the data frame containing the spp name (for ref) and the lm_eqn for
# each regression
gw_spp <- as.data.frame(gw_spp)
colnames(gw_spp) <- c("SpeciesCode", "eqn")

spp2 <-
  ggplot(data=c, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  geom_text(data=gw_spp, aes(x=4.6, y=2.8, label=eqn), parse=TRUE) +
  ggtitle("HorizontalGape")
facet_grid(. ~ cyl, labeller = label_both)

spp2
#### (3) ga plots with regression lines and equations ####
## ga lm_eqns for plotting:
spp_lms_ga <- dlply(c, .(SpeciesCode), function(z) {
  lm(log(area)~log(SLMM), data=z)
})
## pieces for ga_summ:
coefs_ga <- ldply(spp_lms_ga, coef)
rsq <- function(spp_lms_ga) c("rsqd"=summary(spp_lms_ga)$r.squared)
se <- function(spp_lms_ga) c("SE"=summary(spp_lms_ga)$coefficients[2,2])
p_val <- function(spp_lms_ga) c("p_val"=summary(spp_lms_ga)$coefficients[2,4])

# master df for ga lms 
dfs <- list(
  coefs_ga,
  rsqs <- ldply(spp_lms_ga, rsq),
  ses <- ldply(spp_lms_ga, se),
  p_vals <- ldply(spp_lms_ga, p_val)
)

ga_summ <- join_all(dfs, by="SpeciesCode", match="first")

## generating lm_eqns to put onto graphs
m <- ga_summ
ga_spp = matrix(data=NA, nrow=0, ncol=2)
for (i in (1:2)) {
  l <- list(a = format(m[i, 3], digits=2),
            b = format(m[i, 2], digits=2), 
            r2 = format(m[i, 4], digits=2)
  )
  if (l$b >= 0) {
    eqn <- substitute(italic(y) ==a%.% italic(x) + b*","~~italic(r^2)~"="~r2,l)
  }
  
  else {
    l <- list(a = format(m[i, 3], digits=2),
              b = format(abs(m[i, 2]), digits=2), 
              r2 = format(m[i, 4], digits=2)
    )
    eqn <- substitute(italic(y) ==a%.% italic(x) - b*","~~italic(r^2)~"="~r2,l)
  }
  lm_eq = as.character(as.expression(eqn))
  ga_spp <- rbind(ga_spp, c(as.character(m[i,1]), lm_eq))
}
# ga_spp is the data frame containing the spp name (for ref) and the lm_eqn for
# each regression
ga_spp <- as.data.frame(ga_spp)
colnames(ga_spp) <- c("SpeciesCode", "eqn")

spp3 <-
  ggplot(data=c, aes(x=log(SLMM), y=log(area))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  geom_text(data=ga_spp, aes(x=4.6, y=5.4, label=eqn), parse=TRUE) + 
  ggtitle("Gape Area")
facet_grid(. ~ cyl, labeller = label_both)

spp3



































