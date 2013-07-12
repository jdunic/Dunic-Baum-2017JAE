###### Piscivores gape-body size analysis ######

################        Piscvores analysis             ##############

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

df.n <- ddply(.data=p, .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))

lm_eqn = function(x, y){
  m = lm(log(y) ~ log(SLMM), x);
  l <- list(a = format(coef(m)[1], digits = 2), 
            b = format(coef(m)[2], digits = 2), 
            r2 = format(summary(m)$r.squared, digits = 3)
            )
  
  if (l$a >= 0) {
    eq <- substitute(log~italic(y) == b %.% log~italic(x) + a*","~~italic(r)^2~"="~r2, l) 
  }
  
  else {
    l <- list(a = format(abs(coef(m)[1]), digits = 2), 
              b = format(coef(m)[2], digits = 2), 
              r2 = format(summary(m)$r.squared, digits = 3)
    )
    eq <- substitute(log~italic(y) == b %.% log~italic(x) - a*","~~italic(r)^2~"="~r2, l)
  }
  as.character(as.expression(eq));                 
}

c1 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.2, y=4.5, label=lm_eqn(p, p$gh)), parse=TRUE) +
  ggtitle("Vertical Gape")

c2 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.2, y=4.5, label=lm_eqn(p, p$gw)), parse=TRUE) +
  ggtitle("Horizontal Gape")

c3 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(area))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.2, y=8.2, label=lm_eqn(p, p$area)), parse=TRUE) +
  ggtitle("Gape Area")

grid.arrange(c1, c2, c3)

#### gh, gw, ga for all piscivore species colour-coded ####
  # (has regression lines in colour for each species)
cc1 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gh), colour=SpeciesCode)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  ggtitle("Vertical Gape")

cc2 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gw), colour=SpeciesCode)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  ggtitle("Horizontal Gape")

cc3 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(area), colour=SpeciesCode)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  ggtitle("Gape Area")

grid.arrange(cc1, cc2, cc3)

#### gh, gw, ga for all piscivores for EACH SPP (facets) #####
    # (SET OF 3 SECTIONS)
#### (1) gh plots with regression lines and equations ####
## gh lm_eqns for plotting:
spp_lms_gh <- dlply(p, .(SpeciesCode), function(z) {
  lm(log(gh)~log(SLMM), data=z)
})

## pieces for gh_summ:
int <- function(spp_lms_gh) c("int"=coefficients(spp_lms_gh)[1])
slope <- function(spp_lms_gh) c("slope"=coefficients(spp_lms_gh)[2])
rsq <- function(spp_lms_gh) c("rsqd"=summary(spp_lms_gh)$r.squared)
se <- function(spp_lms_gh) c("SE"=summary(spp_lms_gh)$coefficients[2,2])
p_val <- function(spp_lms_gh) c("p_val"=summary(spp_lms_gh)$coefficients[2,4])
spp.n <- ddply(.data=p, .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))

# master df for gh lms 
dfs <- list(
  ints <- ldply(spp_lms_gh, int),
  slopes <- ldply(spp_lms_gh, slope),
  rsqs <- ldply(spp_lms_gh, rsq),
  ses <- ldply(spp_lms_gh, se),
  p_vals <- ldply(spp_lms_gh, p_val),
  spp.n
)

gh_summ <- join_all(dfs, by="SpeciesCode")


## generating lm_eqns to put onto graphs
m <- gh_summ
gh_spp = matrix(data=NA, nrow=0, ncol=2)
for (i in (1:7)) {
  l <- list(a = format(m[i, 3], digits=2),
              b = format(m[i, 2], digits=2), 
              r2 = format(m[i, 4], digits=2)
    )
  if (l$b >= 0) {
    eqn <- substitute(log~italic(y) ==a%.% log~italic(x) + b*~~italic(r)^2~"="~r2, l)
  }
  
  else {
    l <- list(a = format(m[i, 3], digits=2),
              b = format(abs(m[i, 2]), digits=2), 
              r2 = format(m[i, 4], digits=2)
    )
    eqn <- substitute(log~italic(y) ==a%.% log~italic(x) - b*~~italic(r)^2~"="~r2, l)
  }
  lm_eq = as.character(as.expression(eqn))
  gh_spp <- rbind(gh_spp, c(as.character(m[i,1]), lm_eq))
}



# gh_spp is the data frame containing the spp name (for ref) and the lm_eqn for
  # each regression
gh_spp <- as.data.frame(gh_spp)
colnames(gh_spp) <- c("SpeciesCode", "eqn")

spp1 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=4, drop=F) +
  geom_smooth(method=lm) +
  geom_text(data=gh_spp, aes(x=4.7, y=4.7, label=eqn), parse=T) +
  ggtitle("VerticalGape")
  facet_grid(. ~ cyl, labeller = label_both)



spp1
#### (2) gw plots with regression lines and equations ####
## gw lm_eqns for plotting:
spp_lms_gw <- dlply(p, .(SpeciesCode), function(z) {
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
for (i in (1:7)) {
  l <- list(a = format(m[i, 3], digits=2),
       b = format(m[i, 2], digits=2), 
       r2 = format(m[i, 4], digits=2)
  )
  if (l$b >= 0) {
    eqn <- substitute(log(italic(y)) ==a%.% italic(x) + b*","~~italic(r^2)~"="~r2,l)
  }
  
  else {
    l <- list(a = format(m[i, 3], digits=2),
              b = format(abs(m[i, 2]), digits=2), 
              r2 = format(m[i, 4], digits=2)
    )
    eqn <- substitute(log(italic(y)) ==a%.% italic(x) - b*","~~italic(r^2)~"="~r2,l)
  }
  lm_eq = as.character(as.expression(eqn))
  gw_spp <- rbind(gw_spp, c(as.character(m[i,1]), lm_eq))
}
# gw_spp is the data frame containing the spp name (for ref) and the lm_eqn for
# each regression
gw_spp <- as.data.frame(gw_spp)
colnames(gw_spp) <- c("SpeciesCode", "eqn")

spp2 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=4) +
  geom_smooth(method=lm) +
  ylim(c(1,5.3)) +
  geom_text(data=gw_spp, aes(x=4.8, y=5.1, label=eqn), parse=TRUE) +
  ggtitle("HorizontalGape")
  facet_grid(. ~ cyl, labeller = label_both)

spp2

grid.arrange(spp1, spp2)
#### (3) ga plots with regression lines and equations ####
## ga lm_eqns for plotting:
spp_lms_ga <- dlply(p, .(SpeciesCode), function(z) {
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
for (i in (1:7)) {
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
  ggplot(data=p, aes(x=log(SLMM), y=log(area))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=4) +
  geom_smooth(method=lm) +
  ylim(c(1.9, 10.2)) +
  geom_text(data=ga_spp, aes(x=5.0, y=10, label=eqn), parse=TRUE) + 
  ggtitle("Gape Area")
  facet_grid(. ~ cyl, labeller = label_both)

spp3

grid.arrange(spp1, spp2, spp3)

#### gh, gw, ma for all piscivores for EACH FAMILY (facets) #####
# (SET OF 3 SECTIONS)
#### (1) gh plots with regression lines and equations ####
## gh lm_eqns for plotting:
fam_lms_gh <- dlply(p, .(Family), function(z) {
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
  gh_fam <- rbind(gh_fam, c(as.character(m[i,1]), lm_eq))
}
# gh_fam is the data frame containing the fam name (for ref) and the lm_eqn for
# each regression
gh_fam <- as.data.frame(gh_fam)
colnames(gh_fam) <- c("Family", "eqn")

fam1 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=gh_fam, aes(x=5.0, y=5.1, label=eqn), parse=TRUE) +
  ggtitle("VerticalGape")
facet_grid(. ~ cyl, labeller = label_both)






grid.arrange(fam1, fam2, fam3)

#### (2) gw plots with regression lines and equations ####
## gw lm_eqns for plotting:
fam_lms_gw <- dlply(p, .(Family), function(z) {
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
  gw_fam <- rbind(gw_fam, c(as.character(m[i,1]), lm_eq))
}
# gw_fam is the data frame containing the fam name (for ref) and the lm_eqn for
# each regression
gw_fam <- as.data.frame(gw_fam)
colnames(gw_fam) <- c("Family", "eqn")

fam2 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=gw_fam, aes(x=5.0, y=5.2, label=eqn), parse=TRUE) +
  ggtitle("Horizontal Gape")
facet_grid(. ~ cyl, labeller = label_both)

#### (3) ga plots with regression lines and equations ####
## ga lm_eqns for plotting:
fam_lms_ga <- dlply(p, .(Family), function(z) {
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
  ga_fam <- rbind(ga_fam, c(as.character(m[i,1]), lm_eq))
}
# ga_fam is the data frame containing the fam name (for ref) and the lm_eqn for
# each regression
ga_fam <- as.data.frame(ga_fam)
colnames(ga_fam) <- c("Family", "eqn")

fam3 <-
  ggplot(data=p, aes(x=log(SLMM), y=log(area))) + #, colour=SpeciesCode)) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=ga_fam, aes(x=4.8, y=9.9, label=eqn), parse=TRUE) +
  ggtitle("Gape Area")
  facet_grid(. ~ cyl, labeller = label_both)

grid.arrange(fam1, fam2, fam3)

#### gh, gw, ma regression for all piscivore species combined ####
reg.all.p <- lm(formula=log(gh)~log(SLMM), data=p)

ggplot(data=reg.all.p, aes())

#### gh, gw, ga lm for each spp of piscivore in dataframe ####

spp_lms_gh <- dlply(p, .(SpeciesCode), function(z) {
  lm(log(gh)~log(SLMM), data=z)
})

spp_lms_gw <- dlply(p, .(SpeciesCode), function(z) {
  lm(log(gw)~log(SLMM), data=z)
})

spp_lms_ga <- dlply(p, .(SpeciesCode), function(z) {
  lm(log(area)~log(SLMM), data=z)
})

#### lm coefficients (in df) for each piscivore spp (gh, gw, ga) ####
coefs_gh <- ldply(spp_lms_gh, coef)
coefs_gw <- ldply(spp_lms_gw, coef)
coefs_ga <- ldply(spp_lms_ga, coef)

# summary stats of gh lm 
rsq <- function(spp_lms_gh) c("rsqd"=summary(spp_lms_gh)$r.squared)
se <- function(spp_lms_gh) c("SE"=summary(spp_lms_gh)$coefficients[2,2])
p_val <- function(spp_lms_gh) c("p_val"=summary(spp_lms_gh)$coefficients[2,4])

# master df for gh lms 
dfs <- list(
  coefs_gh <- ldply(spp_lms_gh, coef),
  rsqs <- ldply(spp_lms_gh, rsq),
  ses <- ldply(spp_lms_gh, se),
  p_vals <- ldply(spp_lms_gh, p_val)
)

gh_summ <- join_all(dfs, by="SpeciesCode", match="first")

# Function to produce lm expression
# Where: x = data.frame, y = formula

##### lm_eqn function number 1 #####
lm_eqn = function(x, m, rsq) {
  eqn <- substitute(italic(y) == a%.% italic(x) + b*","~~italic(r)^2~"="~r2)
    list(a = format(coef(m)[2], digits = 3),
         b = format(coef(m)[1], digits = 3),
         r2 = format(rsq, digits = 3))
  as.character(as.expression(eqn))
}

spp_gh_eqns <- lm_eqn(p, spp_lms_gh, rsq) 

all.p.lm <- lm_eqn(p, log(gh)~log(SLMM))
print(all.p.lm)

# Fix this lm_eqn2 :(, the insides work, but as a function it doesn't????
#### lm_eqn function number 2... using my MASTER dfs!!!! #####
lm_eqn2 = function(m) { #### FUNCTION NOT WORKING :( :( :( 

  gh.spp = matrix(data=NA, nrow=0, ncol=2)
  for (i in (1:7)) {
    eqn <- substitute(italic(y) ==a%.% italic(x) + b*","~~italic(r^2)~"="~r2, 
                      list(a = format(m[i, 3], digits=3),
                           b = format(m[i, 2], digits=3), 
                          r2 = format(m[i, 4], digits=3)
                      )
            )
    lm_eq = as.character(as.expression(eqn))
    gh.spp <- rbind(gh.spp, c(as.character(m[i,1]), lm_eq))
  }
}

y <- gh_summ

lm_eqn2(gh_summ, length(gh_summ))

cat(lm_eqn2(gh_summ, length(gh_summ)))

lm_eqn3 = function(x){
  m = lm(log(GapeWidth) ~ log(SL), x)
  eqn <- substitute(italic(y) == a%.% italic(x) + b*","~~italic(r)^2~"="~r2,
                    list(a = format(coef(m)[2], digits = 3),
                         b = format(coef(m)[1], digits = 3),
                         r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eqn))
}





































