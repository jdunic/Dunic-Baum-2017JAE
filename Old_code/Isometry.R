###### Piscivores testing Isometry vs Allometry ######

################        Piscvores analysis             ##############

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

#### Adding mouth area to p ####
p$area <- with(p, pi*(gh/2)*(gw/2))

#### Linear models for ALL SPP -- Checking Isometry ####

lm_all = function(y,df) {
  lm(log(y)~log(SLMM), df)
}

lm_all_gh <- lm_all(p$gh, p)
summary(lm_all(p$gh, p))
lm_all_gw <- lm_all(p$gw, p)
summary(lm_all(p$gw, p))
lm_all_ga <- lm_all(p$area, p)
summary(lm_all(p$area, p))

# Confidence interval for ALL SPP:
confint(lm_all_gh, level=0.95)
confint(lm_all_gw, level=0.95)
confint(lm_all_ga, level=0.95)

# Awesome dataframe function. Pulls the confints and coefficients I want from 
# groupwise analyses 
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

#confint1 = lwr intercept
#confint2 = lwr slope
#confint3 = upr intercept
#confint4 = upr slope

# Family level awesome lm dataframe
fam_lms_gh <- dlply(p, .(Family), function(z) {
  lm(log(gh)~log(SLMM), data=z)
})

fam_lms_gw <- dlply(p, .(Family), function(z) {
  lm(log(gw)~log(SLMM), data=z)
})

fam_lms_ga <- dlply(p, .(Family), function(z) {
  lm(log(area)~log(SLMM), data=z)
})

awesome(fam_lms_gh)
awesome(fam_lms_gw)
awesome(fam_lms_ga)


# Species level awesome lm dataframes:
spp_lms_gh <- dlply(p, .(SpeciesCode), function(z) {
  lm(log(gh)~log(SLMM), data=z)
})

spp_lms_gw <- dlply(p, .(SpeciesCode), function(z) {
  lm(log(gw)~log(SLMM), data=z)
})

spp_lms_ga <- dlply(p, .(SpeciesCode), function(z) {
  lm(log(area)~log(SLMM), data=z)
})

awesome(spp_lms_gh)
awesome(spp_lms_gw)
awesome(spp_lms_ga)



###############       Herbivores Analysis        ##########################

# Data for all specimens where a species has more than 3 samples with gape and SL data
fish <- read.csv('Allometry Data/fish_gape_cleaned_more_than_3_jan_24.csv', 
                 header = TRUE, sep = ",") 

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

##### Ordering by Family to make my life easier ####
h$SpeciesCode <- factor(h$SpeciesCode, levels = c("AC.NIGR", 
                                                  "AC.OLIV", 
                                                  "CE.FLAV",
                                                  "CH.SORD",
                                                  "SC.FREN",
                                                  "SC.RUBR"
)
)
#### Adding mouth area to h ####
h$area <- with(h, pi*(gh/2)*(gw/2))

#### Linear models for ALL SPP -- Checking Isometry ####

lm_all = function(y,df) {
  lm(log(y)~log(SLMM), df)
}

lm_all_gh <- lm_all(h$gh, h)
summary(lm_all(h$gh, h))
lm_all_gw <- lm_all(h$gw, h)
lm_all_ga <- lm_all(h$area, h)

# Confidence interval for ALL SPP:
confint(lm_all_gh, level=0.95)
confint(lm_all_gw, level=0.95)
confint(lm_all_ga, level=0.95)

# Awesome dataframe function. Pulls the confints and coefficients I want from 
# groupwise analyses 
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

#confint1 = lwr intercept
#confint2 = lwr slope
#confint3 = upr intercept
#confint4 = upr slope

# Family level awesome lm dataframe
fam_lms_gh <- dlply(h, .(Family), function(z) {
  lm(log(gh)~log(SLMM), data=z)
})

fam_lms_gw <- dlply(h, .(Family), function(z) {
  lm(log(gw)~log(SLMM), data=z)
})

fam_lms_ga <- dlply(h, .(Family), function(z) {
  lm(log(area)~log(SLMM), data=z)
})

awesome(fam_lms_gh)
awesome(fam_lms_gw)
awesome(fam_lms_ga)


# Species level awesome lm dataframes:
spp_lms_gh <- dlply(h, .(SpeciesCode), function(z) {
  lm(log(gh)~log(SLMM), data=z)
})

spp_lms_gw <- dlply(h, .(SpeciesCode), function(z) {
  lm(log(gw)~log(SLMM), data=z)
})

spp_lms_ga <- dlply(h, .(SpeciesCode), function(z) {
  lm(log(area)~log(SLMM), data=z)
})

awesome(spp_lms_gh)
awesome(spp_lms_gw)
awesome(spp_lms_ga)


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
#### Adding mouth area to h ####
b$area <- with(b, pi*(gh/2)*(gw/2))

##### Ordering by Family to make my life easier ####
b$SpeciesCode <- factor(b$SpeciesCode, levels = c("PA.ARCA",
                                                  "MO.GRAN",
                                                  "PA.INSU"
                                                  )
)

#### Linear models for ALL SPP -- Checking Isometry ####

lm_all = function(y,df) {
  lm(log(y)~log(SLMM), df)
}

lm_all_gh <- lm_all(b$gh, b)
summary(lm_all(b$gh, b))
lm_all_gw <- lm_all(b$gw, b)
lm_all_ga <- lm_all(b$area, b)

# Confidence interval for ALL SPP:
confint(lm_all_gh, level=0.95)
confint(lm_all_gw, level=0.95)
confint(lm_all_ga, level=0.95)

# Awesome dataframe function. Pulls the confints and coefficients I want from 
# groupwise analyses 
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

#confint1 = lwr intercept
#confint2 = lwr slope
#confint3 = upr intercept
#confint4 = upr slope

# Family level awesome lm dataframe ####
fam_lms_gh <- dlply(b, .(Family), function(z) {
  lm(log(gh)~log(SLMM), data=z)
})

fam_lms_gw <- dlply(b, .(Family), function(z) {
  lm(log(gw)~log(SLMM), data=z)
})

fam_lms_ga <- dlply(b, .(Family), function(z) {
  lm(log(area)~log(SLMM), data=z)
})

awesome(fam_lms_gh)
awesome(fam_lms_gw)
awesome(fam_lms_ga)


# Species level awesome lm dataframes: ####
spp_lms_gh <- dlply(b, .(SpeciesCode), function(z) {
  lm(log(gh)~log(SLMM), data=z)
})

spp_lms_gw <- dlply(b, .(SpeciesCode), function(z) {
  lm(log(gw)~log(SLMM), data=z)
})

spp_lms_ga <- dlply(b, .(SpeciesCode), function(z) {
  lm(log(area)~log(SLMM), data=z)
})

awesome(spp_lms_gh)
awesome(spp_lms_gw)
awesome(spp_lms_ga)


























