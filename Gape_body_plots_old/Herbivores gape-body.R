###### Herbivores gape-body size analysis ######

################        Herbivores analysis             ##############

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


##### Selecting only herbivores #######

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


#### Making some plots!!! ######

library("ggplot2")
library('plyr')
library('gridExtra') # allows me to make multiplots

#### calculating mouth area (and adding it to herb dataframe): ####
h$area <- with(h, pi*(gh/2)*(gw/2))

#### gh, gw, ga for all herbivore species combined ####
#### (1) gh plots with regression lines and equations COMBINED SPP ####
## gh lm_eqns for plotting:
all_lms_gh <- lm(log(gh)~log(SLMM), data=h)
               
## pieces for gh_summ:
dfs <- list(int <- coef(all_lms_gh)[[1]],
            slope <- coef(all_lms_gh)[[2]],
            rsq <- summary(all_lms_gh)$r.squared,
            se <- summary(all_lms_gh)$coefficients[2,2],
            p_val <- summary(all_lms_gh)$coefficients[2,4]
)
# master df for gh lms 
df <- data.frame(t(rep(NA, 5)))
names(df) <- c("int", "slope", "rsq", "se", "p_val")
df <- df[-1,]
df[1,] <- v_gh

v_gh <- c(int, slope, rsq, se, p_val)

gh_summ <- join_all(df)

all_lms <- lm(log(gh)~log(SLMM), data=h)


## FUNCTION to generate lm_eqns to put onto total spp graphs
lm_all_eqn <- function(df, gape) {
  all_lms <- lm(log(gape)~log(SLMM), data=df)
  
  ## pieces for gh_summ:
  int <- coef(all_lms)[[1]]
  slope <- coef(all_lms)[[2]]
  rsq <- summary(all_lms)$r.squared
  se <- summary(all_lms)$coefficients[2,2]
  p_val <- summary(all_lms)$coefficients[2,4]
  v_gh <- c(int, slope, rsq, se, p_val)
  
  l <- list(a = format(slope, digits=2),
            b = format(int, digits=2),
            r2 = format(rsq, digits=2)
  )
  
  if (int >= 0) {
    eqn <- substitute(italic(y) ==a%.% italic(x) + b*","~~italic(r^2)~"="~r2, l)
    as.character(as.expression(eqn))
  }
  
  else {
    l <- list(a = format(slope, digits=2),
              b = format(abs(int), digits=2),
              r2 = format(rsq, digits=2)
    )
    eqn <- substitute(italic(y) ==a%.% italic(x) - b*","~~italic(r^2)~"="~r2,l)
    as.character(as.expression(eqn))
  }
}

# gh_spp is the data frame containing the spp name (for ref) and the lm_eqn for
# each regression
gh_spp <- as.data.frame(gh_spp)
colnames(gh_spp) <- c("SpeciesCode", "eqn")

c1 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +  
  geom_text(data=data.frame(), aes(x=4.4, y=3.8, label=lm_all_eqn(h, h$gh)), parse=TRUE) +
  ggtitle("Vertical Gape")

c2 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.4, y=3.4, label=lm_all_eqn(h, h$gw)), parse=TRUE) +
  ggtitle("Horizontal Gape")

c3 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(area))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.4, y=6.7, label=lm_all_eqn(h, h$area)), parse=TRUE) +
  ggtitle("Gape Area")

grid.arrange(c1, c2, c3)


#### gh, gw, ga for all piscivore species colour-coded ####
# (has regression lines in colour for each species)
cc1 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gh), colour=SpeciesCode)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  ggtitle("Vertical Gape")

cc2 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gw), colour=SpeciesCode)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  ggtitle("Horizontal Gape")

cc3 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(area), colour=SpeciesCode)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  ggtitle("Gape Area")

grid.arrange(cc1, cc2, cc3)

#### gh, gw, ga for all herbivores for EACH SPP (facets) #####
# (SET OF 3 SECTIONS)
#### (1) gh plots with regression lines and equations ####
## gh lm_eqns for plotting:
spp_lms_gh <- dlply(h, .(SpeciesCode), function(z) {
  lm(log(gh)~log(SLMM), data=z)
})

## pieces for gh_summ:
coefs_gh <- ldply(spp_lms_gh, coef)
rsq <- function(spp_lms_gh) c("rsqd"=summary(spp_lms_gh)$r.squared)
se <- function(spp_lms_gh) c("SE"=summary(spp_lms_gh)$coefficients[2,2])
p_val <- function(spp_lms_gh) c("p_val"=summary(spp_lms_gh)$coefficients[2,4])

# master df for gh lms 
dfs <- list(
  coefs_gh,
  rsqs <- ldply(spp_lms_gh, rsq),
  ses <- ldply(spp_lms_gh, se),
  p_vals <- ldply(spp_lms_gh, p_val)
)

gh_summ <- join_all(dfs, by="SpeciesCode")

## generating lm_eqns to put onto graphs

m <- gh_summ
gh_spp = matrix(data=NA, nrow=0, ncol=2)
for (i in (1:6)) {
  l <- list(a = format(m[i,3], digits=2),
            b = format(m[i,2], digits=2),
            r2 = format(m[i,4], digits=2)
            )
  if (l$b >= 0) {
    eqn <- substitute(italic(y) ==a%.% italic(x) + b*","~~italic(r^2)~"="~r2, l)
  }
  
  else {
    l <- list(a = format(m[i,3], digits=2),
              b = format(abs(m[i,2]), digits=2),
              r2 = format(m[i,4], digits=2)
    )
    eqn <- substitute(italic(y) ==a%.% italic(x) - b*","~~italic(r^2)~"="~r2, l)
  }

  lm_eq = as.character(as.expression(eqn))
  gh_spp <- rbind(gh_spp, c(as.character(m[i,1]), lm_eq))
}
# gh_spp is the data frame containing the spp name (for ref) and the lm_eqn for
# each regression
gh_spp <- as.data.frame(gh_spp)
colnames(gh_spp) <- c("SpeciesCode", "eqn")

spp1 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
  geom_text(data=gh_spp, aes(x=4.9, y=4.2, label=eqn), parse=TRUE) +
  ggtitle("VerticalGape")
  facet_grid(. ~ cyl, labeller = label_both)

#### (2) gw plots with regression lines and equations ####
## gw lm_eqns for plotting:
spp_lms_gw <- dlply(h, .(SpeciesCode), function(z) {
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
for (i in (1:6)) {
  l <- list(a = format(m[i,3], digits=2),
            b = format(m[i,2], digits=2),
            r2 = format(m[i,4], digits=2)
  )
  if (l$b >= 0) {
    eqn <- substitute(italic(y) ==a%.% italic(x) + b*","~~italic(r^2)~"="~r2, l)
  }
  
  else {
    l <- list(a = format(m[i,3], digits=2),
              b = format(abs(m[i,2]), digits=2),
              r2 = format(m[i,4], digits=2)
    )
    eqn <- substitute(italic(y) ==a%.% italic(x) - b*","~~italic(r^2)~"="~r2, l)
  }
  lm_eq = as.character(as.expression(eqn))
  gw_spp <- rbind(gw_spp, c(as.character(m[i,1]), lm_eq))
}
# gw_spp is the data frame containing the spp name (for ref) and the lm_eqn for
# each regression
gw_spp <- as.data.frame(gw_spp)
colnames(gw_spp) <- c("SpeciesCode", "eqn")

spp2 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
#  ylim(c(1,5.3)) +
  geom_text(data=gw_spp, aes(x=4.8, y=3.7, label=eqn), parse=TRUE) +
  ggtitle("HorizontalGape")
  facet_grid(. ~ cyl, labeller = label_both)

#### (3) ga plots with regression lines and equations ####
## ga lm_eqns for plotting:
spp_lms_ga <- dlply(h, .(SpeciesCode), function(z) {
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
  
for (i in (1:6)) {
  l <- list(a = format(m[i,3], digits=2),
            b = format(m[i,2], digits=2),
            r2 = format(m[i,4], digits=2)
  )
  if (l$b >= 0) {
    eqn <- substitute(italic(y) ==a%.% italic(x) + b*","~~italic(r^2)~"="~r2, l)
  }
  
  else {
    l <- list(a = format(m[i,3], digits=2),
              b = format(abs(m[i,2]), digits=2),
              r2 = format(m[i,4], digits=2)
    )
    eqn <- substitute(italic(y) ==a%.% italic(x) - b*","~~italic(r^2)~"="~r2, l)
  }
  lm_eq = as.character(as.expression(eqn))
  ga_spp <- rbind(ga_spp, c(as.character(m[i,1]), lm_eq))
}
# ga_spp is the data frame containing the spp name (for ref) and the lm_eqn for
# each regression
ga_spp <- as.data.frame(ga_spp)
colnames(ga_spp) <- c("SpeciesCode", "eqn")

spp3 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(area))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
#  ylim(c(1.9, 10.2)) +
  geom_text(data=ga_spp, aes(x=4.7, y=7.5, label=eqn), parse=TRUE) + 
  ggtitle("Gape Area")
  facet_grid(. ~ cyl, labeller = label_both)

spp1
spp2
spp3

#### gh, gw, ma for all herbivores for EACH FAMILY (facets) #####
# (SET OF 3 SECTIONS)
#### (1) gh plots with regression lines and equations ####
## gh lm_eqns for plotting:
fam_lms_gh <- dlply(h, .(Family), function(z) {
  lm(log(gh)~log(SLMM), data=z)
})
## pieces for gh_summ:
coefs_gh <- ldply(fam_lms_gh, coef)
rsq <- function(fam_lms_gh) c("rsqd"=summary(fam_lms_gh)$r.squared)
se <- function(fam_lms_gh) c("SE"=summary(fam_lms_gh)$coefficients[2,2])
p_val <- function(fam_lms_gh) c("p_val"=summary(fam_lms_gh)$coefficients[2,4])

# master df for gh lms 
dfs <- list(
  coefs_gh,
  rsqs <- ldply(fam_lms_gh, rsq),
  ses <- ldply(fam_lms_gh, se),
  p_vals <- ldply(fam_lms_gh, p_val)
)

gh_summ <- join_all(dfs, by="Family", match="first")

## generating lm_eqns to put onto graphs
m <- gh_summ
gh_fam = matrix(data=NA, nrow=0, ncol=2)
for (i in (1:3)) {
  l <- list(a = format(m[i,3], digits=2),
            b = format(m[i,2], digits=2),
            r2 = format(m[i,4], digits=2)
  )
  if (l$b >= 0) {
    eqn <- substitute(italic(y) ==a%.% italic(x) + b*","~~italic(r^2)~"="~r2, l)
  }
  
  else {
    l <- list(a = format(m[i,3], digits=2),
              b = format(abs(m[i,2]), digits=2),
              r2 = format(m[i,4], digits=2)
    )
    eqn <- substitute(italic(y) ==a%.% italic(x) - b*","~~italic(r^2)~"="~r2, l)
  }
  
  lm_eq = as.character(as.expression(eqn))
  gh_fam <- rbind(gh_fam, c(as.character(m[i,1]), lm_eq))
}
# gh_fam is the data frame containing the fam name (for ref) and the lm_eqn for
# each regression
gh_fam <- as.data.frame(gh_fam)
colnames(gh_fam) <- c("Family", "eqn")

fam1 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=gh_fam, aes(x=4.8, y=4.2, label=eqn), parse=TRUE) +
  ggtitle("VerticalGape")
facet_grid(. ~ cyl, labeller = label_both)


#### (2) gw plots with regression lines and equations ####
## gw lm_eqns for plotting:
fam_lms_gw <- dlply(h, .(Family), function(z) {
  lm(log(gw)~log(SLMM), data=z)
})
## pieces for gw_summ:
coefs_gw <- ldply(fam_lms_gw, coef)
rsq <- function(fam_lms_gw) c("rsqd"=summary(fam_lms_gw)$r.squared)
se <- function(fam_lms_gw) c("SE"=summary(fam_lms_gw)$coefficients[2,2])
p_val <- function(fam_lms_gw) c("p_val"=summary(fam_lms_gw)$coefficients[2,4])

# master df for gw lms 
dfs <- list(
  coefs_gw,
  rsqs <- ldply(fam_lms_gw, rsq),
  ses <- ldply(fam_lms_gw, se),
  p_vals <- ldply(fam_lms_gw, p_val)
)

gw_summ <- join_all(dfs, by="Family", match="first")

## generating lm_eqns to put onto graphs
m <- gw_summ
gw_fam = matrix(data=NA, nrow=0, ncol=2)
for (i in (1:3)) {
  l <- list(a = format(m[i,3], digits=2),
            b = format(m[i,2], digits=2),
            r2 = format(m[i,4], digits=2)
  )
  if (l$b >= 0) {
    eqn <- substitute(italic(y) ==a%.% italic(x) + b*","~~italic(r^2)~"="~r2, l)
  }
  
  else {
    l <- list(a = format(m[i,3], digits=2),
              b = format(abs(m[i,2]), digits=2),
              r2 = format(m[i,4], digits=2)
    )
    eqn <- substitute(italic(y) ==a%.% italic(x) - b*","~~italic(r^2)~"="~r2, l)
  }
  
  lm_eq = as.character(as.expression(eqn))
  gw_fam <- rbind(gw_fam, c(as.character(m[i,1]), lm_eq))
}
# gw_fam is the data frame containing the fam name (for ref) and the lm_eqn for
# each regression
gw_fam <- as.data.frame(gw_fam)
colnames(gw_fam) <- c("Family", "eqn")

fam2 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=gw_fam, aes(x=4.9, y=3.7, label=eqn), parse=TRUE) +
  ggtitle("Horizontal Gape")
facet_grid(. ~ cyl, labeller = label_both)

#### (3) ga plots with regression lines and equations ####
## ga lm_eqns for plotting:
fam_lms_ga <- dlply(h, .(Family), function(z) {
  lm(log(area)~log(SLMM), data=z)
})
## pieces for ga_summ:
coefs_ga <- ldply(fam_lms_ga, coef)
rsq <- function(fam_lms_ga) c("rsqd"=summary(fam_lms_ga)$r.squared)
se <- function(fam_lms_ga) c("SE"=summary(fam_lms_ga)$coefficients[2,2])
p_val <- function(fam_lms_ga) c("p_val"=summary(fam_lms_ga)$coefficients[2,4])

# master df for ga lms 
dfs <- list(
  coefs_ga,
  rsqs <- ldply(fam_lms_ga, rsq),
  ses <- ldply(fam_lms_ga, se),
  p_vals <- ldply(fam_lms_ga, p_val)
)

ga_summ <- join_all(dfs, by="Family", match="first")

## generating lm_eqns to put onto graphs
m <- ga_summ
ga_fam = matrix(data=NA, nrow=0, ncol=2)
for (i in (1:3)) {
  l <- list(a = format(m[i,3], digits=2),
            b = format(m[i,2], digits=2),
            r2 = format(m[i,4], digits=2)
  )
  if (l$b >= 0) {
    eqn <- substitute(italic(y) ==a%.% italic(x) + b*","~~italic(r^2)~"="~r2, l)
  }
  
  else {
    l <- list(a = format(m[i,3], digits=2),
              b = format(abs(m[i,2]), digits=2),
              r2 = format(m[i,4], digits=2)
    )
    eqn <- substitute(italic(y) ==a%.% italic(x) - b*","~~italic(r^2)~"="~r2, l)
  }
  
  lm_eq = as.character(as.expression(eqn))
  ga_fam <- rbind(ga_fam, c(as.character(m[i,1]), lm_eq))
}
# ga_fam is the data frame containing the fam name (for ref) and the lm_eqn for
# each regression
ga_fam <- as.data.frame(ga_fam)
colnames(ga_fam) <- c("Family", "eqn")

fam3 <-
  ggplot(data=h, aes(x=log(SLMM), y=log(area))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=ga_fam, aes(x=4.8, y=7.5, label=eqn), parse=TRUE) +
  ggtitle("Gape Area")
facet_grid(. ~ cyl, labeller = label_both)

grid.arrange(fam1, fam2, fam3)