###### Benthic Invertivores gape-body size analysis ISOMETRY ######

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

b$gh_ratio <- b$gh/b$SLMM
b$gw_ratio <- b$gw/b$SLMM
b$ga_ratio <- b$area/b$SLMM

#### Making some plots!!! ######

library("ggplot2")
library('plyr')
library('gridExtra') # allows me to make multiplots

#### Plots of isometry ###

lm_all_eqn <- function(df, gape) {
  all_lms <- lm(gape~SLMM, data=df)
  
  ## pieces for gh_summ:
  int <- coef(all_lms)[[1]]
  slope <- coef(all_lms)[[2]]
  rsq <- summary(all_lms)$r.squared
  se <- summary(all_lms)$coefficients[2,2]
  p_val <- summary(all_lms)$coefficients[2,4]
  v_gh <- c(int, slope, rsq, se, p_val)
  
  l <- list(a = format(slope, digits=2),
            b = format(abs(int), digits=2),
            r2 = format(rsq, digits=2)
  )
  
  if (int >= 0) {
    eqn <- substitute(italic(y) ==a%.% italic(x) + b*","~~italic(r^2)~"="~r2, l)
    as.character(as.expression(eqn))
  }
  
  else {
    eqn <- substitute(italic(y) ==a%.% italic(x) - b*","~~italic(r^2)~"="~r2,l)
    as.character(as.expression(eqn))
  }
}

# gh_spp is the data frame containing the spp name (for ref) and the lm_eqn for
# each regression
gh_spp <- as.data.frame(gh_spp)
colnames(gh_spp) <- c("SpeciesCode", "eqn")

#### non-transformed ALL SPP plots ####
c1 <-
  ggplot(data=b, aes(x=SLMM, y=gh)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +  
  geom_text(data=data.frame(), aes(x=100, y=50, label=lm_all_eqn(b, b$gh)), parse=TRUE) +
  ggtitle("Vertical Gape")

c2 <-
  ggplot(data=b, aes(x=SLMM, y=gw)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=100, y=31, label=lm_all_eqn(b, b$gw)), parse=TRUE) +
  ggtitle("Horizontal Gape")

c3 <-
  ggplot(data=b, aes(x=SLMM, y=area)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  geom_text(data=data.frame(), aes(x=100, y=1200, label=lm_all_eqn(b, b$area)), parse=TRUE) +
  ggtitle("Gape Area")

grid.arrange(c1, c2, c3)

#### Gives me an awesome df with lots of important lm values ####
awesome = function(lm) {
  ldply(lm, function(model) {
    c(
      "coefs" = coef(model),
      "confint" = confint(model, level=0.95),
      "rsq" = summary(model)$r.squared,
      "se" = summary(model)$coefficients[2,2],
      "p_val" = summary(model)$coefficients[2,4]
    )
  })
}

#### NON-transformed regressions (h$family) ####
fam_lms_gh <- dlply(b, .(Family), function(z) {
  lm(gh~SLMM, data=z)
})

fam_lms_gw <- dlply(b, .(Family), function(z) {
  lm(gw~SLMM, data=z)
})

fam_lms_ga <- dlply(b, .(Family), function(z) {
  lm(area~SLMM, data=z)
})

## generating lm_eqns to put onto graphs
m <- awesome(fam_lms_gh)

gh_fam = matrix(data=NA, nrow=0, ncol=2)

# builds eqn and stores in gh_fam
for (i in (1:3)) {
  l <- list(a = format(m[i, 3], digits=2),
            b = format(abs(m[i, 2]), digits=2), 
            r2 = format(m[i, 8], digits=2)
  )
  if (l$b >= 0) {
    eqn <- substitute(italic(y) ==a%.% italic(x) + b*","~~italic(r^2)~"="~r2,l)
  }
  
  else {
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
  ggplot(data=b, aes(x=SLMM, y=gh)) +
  geom_point(shape=20, aes(colour=SpeciesCode)) +
  facet_wrap(~Family) +
  geom_smooth(method=lm) +
  geom_text(data=gh_fam, aes(x=280, y=120, label=eqn), parse=TRUE) +
  ggtitle("VerticalGape")
facet_grid(. ~ cyl, labeller = label_both)

# Checking standardized residuals for normality: 
l_ply(fam_lms_gh, function(model) {plot(rstandard(model)); abline(0,0)})

# Normal probability plots of standardized residuals for fam lms:
l_ply(fam_lms_gh, function(model) {qqnorm(rstandard(model)); qqline(rstandard(model))})

# Checking standardized residuals for normality: 
l_ply(fam_lms_gw, function(model) {plot(rstandard(model)); abline(0,0)})

# Normal probability plots of standardized residuals for fam lms:
l_ply(fam_lms_gw, function(model) {qqnorm(rstandard(model)); qqline(rstandard(model))})

# Checking standardized residuals for normality: 
l_ply(fam_lms_ga, function(model) {plot(rstandard(model)); abline(0,0)})

# Normal probability plots of standardized residuals for fam lms:
l_ply(fam_lms_ga, function(model) {qqnorm(rstandard(model)); qqline(rstandard(model))})

##### Mean ratio of gh, gw, ga to SLMM ####
# Standard Error of a mean function:
se_mean <- function(x) {
  sqrt(var(x)/length(x))
}

# All species combined:
gh_all <- mean(b$gh_ratio)
gh_all_se <- se_mean(b$gh_ratio)

gw_all <- mean(b$gw_ratio)
gw_all_se <- se_mean(b$gw_ratio)

ga_all <- mean(b$ga_ratio)
ga_all_se <- se_mean(b$ga_ratio)

# By family:
fam_gh_rats <- ddply(b, .(Family), summarise, "mean"=mean(gh_ratio), "se"=se_mean(gh_ratio))
fam_gw_rats <- ddply(b, .(Family), summarise, "mean"=mean(gw_ratio), "se"=se_mean(gw_ratio))
fam_ga_rats <- ddply(b, .(Family), summarise, "mean"=mean(ga_ratio), "se"=se_mean(ga_ratio))

fam_gh_rats
fam_gw_rats
fam_ga_rats

# By Species:
spp_gh_rats <- ddply(b, .(SpeciesCode), summarise, "mean"=mean(gh_ratio), "se"=se_mean(gh_ratio))
spp_gw_rats <- ddply(b, .(SpeciesCode), summarise, "mean"=mean(gw_ratio), "se"=se_mean(gw_ratio))
spp_ga_rats <- ddply(b, .(SpeciesCode), summarise, "mean"=mean(ga_ratio), "se"=se_mean(ga_ratio))

spp_gh_rats
spp_gw_rats
spp_ga_rats
