###### Benthic Invertivores gape-body size analysis ######

################        Benthic Invertivores analysis             ##############

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


##### Selecting only benthic invertivores #######

b <- fish[which(fish$j_fg=='BI'),]

b$SpeciesCode <- factor(b$SpeciesCode)
b$Family <- factor(b$Family)

##### Ordering by Family to make my life easier ####
b$SpeciesCode <- factor(b$SpeciesCode, levels = c("PA.ARCA",
                                                  "MO.GRAN",
                                                  "PA.INSU"
                                                  )
)

#### Making some plots!!! ######

library("ggplot2")
library('plyr')
library('gridExtra') # allows me to make multiplots

#### calculating mouth area (and adding it to BI dataframe): ####
b$area <- with(b, pi*(gh/2)*(gw/2))

#### gh, gw, ga for all BI species combined (WITH EQNS!!!!!!) ####

df.n <- ddply(.data=b, .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))

lm_eqn = function(x, y){
  m = lm(log(y) ~ log(SLMM), x);
  l <- list(b = format(coef(m)[1], digits = 2), 
            a = format(coef(m)[2], digits = 2), 
            r2 = format(summary(m)$r.squared, digits = 3)
  )
  
  if (l$b >= 0) {
    eq <- substitute(italic(y) ==  a %.% italic(x) + b*","~~italic(r)^2~"="~r2, l) 
  }
  
  else {
    l <- list(b = format(abs(coef(m)[1]), digits = 2), 
              a = format(coef(m)[2], digits = 2), 
              r2 = format(summary(m)$r.squared, digits = 3)
    )
    eq <- substitute(italic(y) == a %.% italic(x) - b*","~~italic(r)^2~"="~r2, l)
  }
  as.character(as.expression(eq));                 
}

c1 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.5, y=3.7, label=lm_eqn(b, b$gh)), parse=TRUE) +
  ggtitle("Vertical Gape")

c2 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.5, y=3.7, label=lm_eqn(b, b$gw)), parse=TRUE) +
  ggtitle("Horizontal Gape")

c3 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(area))) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=4.5, y=7.5, label=lm_eqn(b, b$area)), parse=TRUE) +
  ggtitle("Gape Area")

grid.arrange(c1, c2, c3)

#### gh, gw, ga for all piscivore species colour-coded ####
# (has regression lines in colour for each species)
cc1 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(gh), colour=SpeciesCode)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  ggtitle("Vertical Gape")

cc2 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(gw), colour=SpeciesCode)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  ggtitle("Horizontal Gape")

cc3 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(area), colour=SpeciesCode)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  ggtitle("Gape Area")

grid.arrange(cc1, cc2, cc3)

#### gh, gw, ga for all BIs for EACH SPP (facets) #####
# (SET OF 3 SECTIONS)
#### (1) gh plots with regression lines and equations ####
## gh lm_eqns for plotting:
spp_lms_gh <- dlply(b, .(SpeciesCode), function(z) {
  lm(log(gh)~log(SLMM), data=z)
})

spp.n <- ddply(.data=b, .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))

## pieces for gh_summ:
coefs_gh <- ldply(spp_lms_gh, coef)
rsq <- function(spp_lms_gh) c("rsqd"=summary(spp_lms_gh)$r.squared)
se <- function(spp_lms_gh) c("SE"=summary(spp_lms_gh)$coefficients[2,2])
p_val <- function(spp_lms_gh) c("p_val"=summary(spp_lms_gh)$coefficients[2,4])
spp.n <- ddply(.data=b, .(SpeciesCode), summarize, n=paste("n ==", length(SpeciesCode)))

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
for (i in (1:3)) {
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

#### screwing around trying to get the regression eqn on one line and r^2 and n on another ####
m <- gh_summ
gh_spp = data.frame(data=NA, nrow=0, ncol=4)

gh_spp <- NULL

for (i in (1:3)) {
  l <- list(a = format(m[i, 3], digits=2),
            b = format(abs(m[i, 2]), digits=2), 
            r2 = format(m[i, 4], digits=2),
            n = m[i,7]
  )
  eqn <- as.character(as.expression(substitute(italic(y) ==a%.% italic(x) + b, l)))
  line2 <- as.character(as.expression(substitute(italic(r)^2~"="~r2 ~~ n, l)))
  exp <- capture.output(expression(atop(paste(eqn,line2))))
  gh_spp <- rbind(gh_spp, data.frame("SpeciesCode"=as.character(m[i,1]), "eqn"=eqn, "line2"=line2, "exp"=exp))
}

# gh_spp is the data frame containing the spp name (for ref) and the lm_eqn for

# each regression ####
gh_spp <- as.data.frame(gh_spp)
colnames(gh_spp) <- c("SpeciesCode", "lm_eq")

spp1 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(gh))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3, drop=F) +
  geom_smooth(method=lm) +
  geom_text(data=gh_spp, aes(x=4.85, y=4.1, label=lm_eq), parse=T) +
  ggtitle("VerticalGape")
facet_grid(. ~ cyl, labeller = label_both)


spp1
#### (2) gw plots with regression lines and equations ####
## gw lm_eqns for plotting:
spp_lms_gw <- dlply(b, .(SpeciesCode), function(z) {
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
  gw_spp <- rbind(gw_spp, c(as.character(m[i,1]), lm_eq))
}
# gw_spp is the data frame containing the spp name (for ref) and the lm_eqn for
# each regression
gw_spp <- as.data.frame(gw_spp)
colnames(gw_spp) <- c("SpeciesCode", "eqn")

spp2 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(gw))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
#  ylim(c(1,5.3)) +
  geom_text(data=gw_spp, aes(x=4.85, y=3.9, label=eqn), parse=TRUE) +
  ggtitle("HorizontalGape")
facet_grid(. ~ cyl, labeller = label_both)

spp2
#### (3) ga plots with regression lines and equations ####
## ga lm_eqns for plotting:
spp_lms_ga <- dlply(b, .(SpeciesCode), function(z) {
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
  ga_spp <- rbind(ga_spp, c(as.character(m[i,1]), lm_eq))
}
# ga_spp is the data frame containing the spp name (for ref) and the lm_eqn for
# each regression
ga_spp <- as.data.frame(ga_spp)
colnames(ga_spp) <- c("SpeciesCode", "eqn")

spp3 <-
  ggplot(data=b, aes(x=log(SLMM), y=log(area))) +
  geom_point(shape=20) +
  facet_wrap(~SpeciesCode, ncol=3) +
  geom_smooth(method=lm) +
#  ylim(c(1.9, 10.2)) +
  geom_text(data=ga_spp, aes(x=4.85, y=7.7, label=eqn), parse=TRUE) + 
  ggtitle("Gape Area")
facet_grid(. ~ cyl, labeller = label_both)

spp3

grid.arrange(spp1, spp2, spp3)

#### Family level analysis not necessary --> one spp per family ####

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




































