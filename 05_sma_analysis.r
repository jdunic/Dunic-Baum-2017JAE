# all functional groups together:
pento_gh <- sma(gh~SL, data=pento, log="xy", method="SMA", robust=T, slope.test=1)
plot(pento_gh)
pento_gw <- sma(gw~SL, data=pento, log="xy", method="SMA", robust=T, slope.test=1)
plot(pento_gw)
pento_ga <- sma(ga~SL, data=pento, log="xy", method="SMA", robust=T, slope.test=2)
plot(pento_ga)

# sma of piscivores as a whole:
pisc <- sma(gh~SL, data=p, log="xy", method="SMA", robust=T)
plot(pisc)
plot(pisc, which="residual")
abline(h=0)
hist(pisc$resid)
str(pisc)

# Testing for common slope across the different functional groups:
pento_gh <- sma(gh~SL*j_fg, data=pento, log="xy", robust=T, method="SMA", slope.test=1)
pento_gw <- sma(gw~SL*j_fg, data=pento, log="xy", robust=T, method="SMA", slope.test=1)
pento_ga <- sma(ga~SL*j_fg, data=pento, log="xy", robust=T, method="SMA", slope.test=2)

plot(pento_gh)
plot(pento_gw)
plot(pento_ga)
legend("topleft", legend=unique(pento$j_fg))

# Get results for individual functional groups:
pento_gh$groupsummary
pento_gw$groupsummary
pento_ga$groupsummary


# Checking if slopes are common across family within functional groups:
pfam_gh  <- sma(gh~SL*Family, data=p, log="xy", method="SMA", robust=T)
bfam_gh  <- sma(gh~SL*Family, data=b, log="xy", method="SMA", robust=T)
zpfam_gh <- sma(gh~SL*Family, data=zp, log="xy", method="SMA", robust=T)
hfam_gh  <- sma(gh~SL*Family, data=h, log="xy", method="SMA", robust=T)

pfam_gw  <- sma(gw~SL*Family, data=p, log="xy", method="SMA", robust=T)
bfam_gw  <- sma(gw~SL*Family, data=b, log="xy", method="SMA", robust=T)
zpfam_gw <- sma(gw~SL*Family, data=zp, log="xy", method="SMA", robust=T)
hfam_gw  <- sma(gw~SL*Family, data=h, log="xy", method="SMA", robust=T)

pfam_ga  <- sma(ga~SL*Family, data=p, log="xy", method="SMA", robust=T) #, slope.test=1)
bfam_ga  <- sma(ga~SL*Family, data=b, log="xy", method="SMA", robust=T)
zpfam_ga <- sma(ga~SL*Family, data=zp, log="xy", method="SMA", robust=T)
hfam_ga  <- sma(ga~SL*Family, data=h, log="xy", method="SMA", robust=T)


pfam_gh_slp <- with(p, slope.com(gh, SL, Family, robust=T, slope.test=1))


# Checking if slopes are common across SPECIES within functional groups:
pspp_gh <- sma(gh~SL*SpeciesCode, data=p, log="xy", method="SMA", robust=T)
bspp_gh <- sma(gh~SL*SpeciesCode, data=b, log="xy", method="SMA", robust=T)
zpspp_gh <- sma(gh~SL*SpeciesCode, data=zp, log="xy", method="SMA", robust=T)
hspp_gh <- sma(gh~SL*SpeciesCode, data=h, log="xy", method="SMA", robust=T,
               slope.test=1)


pspp_gw <- sma(gw~SL*SpeciesCode, data=p, log="xy", method="SMA", robust=T)
bspp_gw <- sma(gw~SL*SpeciesCode, data=b, log="xy", method="SMA", robust=T)
zpspp_gw <- sma(gw~SL*SpeciesCode, data=zp, log="xy", method="SMA", robust=T)
hspp_gw <- sma(gw~SL*SpeciesCode, data=h, log="xy", method="SMA", robust=T)


pspp_ga <- sma(ga~SL*SpeciesCode, data=p, log="xy", method="SMA", robust=T,
               slope.test=2)
bspp_ga <- sma(ga~SL*SpeciesCode, data=b, log="xy", method="SMA", robust=T)
zpspp_ga <- sma(ga~SL*SpeciesCode, data=zp, log="xy", method="SMA", robust=T)
hspp_ga <- sma(ga~SL*SpeciesCode, data=h, log="xy", method="SMA", robust=T)




# checking SMA assumptions for each functional group: 
plot(p_sma, which="residual")
abline(h=0)
plot(p_sma, which="qq")
plot(b_sma, which="residual")
abline(h=0)
plot(b_sma, which="qq")
plot(zp_sma, which="residual")
abline(h=0)
plot(zp_sma, which="qq")
plot(h_sma, which="residual")
abline(h=0)
plot(h_sma, which="qq")

# dataframe for functional group SMAs:
j_fg <- c('Pi', 'BI', 'ZP', 'H', 'C')
p_sma_df <- mk_sma_df(p_sma)
b_sma_df <- mk_sma_df(b_sma)
zp_sma_df <- mk_sma_df(zp_sma)
h_sma_df <- mk_sma_df(h_sma)
c_sma_df <- mk_sma_df(c_sma)
master <- rbind(p_sma_df, b_sma_df, zp_sma_df, h_sma_df, c_sma_df)
master <- cbind(j_fg, master)

# checking SMA assumptions for each functional group, with Species groups:
p_sma_spp <- sma(gh~SL*SpeciesCode, data=p, log="xy", method="SMA", slope.test=1)
b_sma_spp <- sma(gh~SL*SpeciesCode, data=b, log="xy", method="SMA", slope.test=1)
zp_sma_spp <- sma(gh~SL*SpeciesCode, data=zp, log="xy", method="SMA", slope.test=1)
h_sma_spp <- sma(gh~SL*SpeciesCode, data=h, log="xy", method="SMA", slope.test=1)

plot(p_sma_spp, which="residual")
abline(h=0)
plot(p_sma_spp, which="qq")
plot(b_sma_spp, which="residual")
abline(h=0)
plot(b_sma_spp, which="qq")
plot(zp_sma_spp, which="residual")
abline(h=0)
plot(zp_sma_spp, which="qq")
plot(h_sma_spp, which="residual")
abline(h=0)
plot(h_sma_spp, which="qq")


# Checking if all species have a common slope:
p_gh <- sma(gh~SL*SpeciesCode, data=p, log="xy", method="SMA", robust=T, slope.test=1)
b_gh <- sma(gh~SL*SpeciesCode, data=b, log="xy", method="SMA", robust=T, slope.test=1)
zp_gh <- sma(gh~SL*SpeciesCode, data=zp, log="xy", method="SMA", robust=T, slope.test=1)
h_gh <- sma(gh~SL*SpeciesCode, data=h, log="xy", method="SMA", robust=T, slope.test=1)

p_gw <- sma(gw~SL*SpeciesCode, data=p, log="xy", method="SMA", robust=T, slope.test=1)
b_gw <- sma(gw~SL*SpeciesCode, data=b, log="xy", method="SMA", robust=T, slope.test=1)
zp_gw <- sma(gw~SL*SpeciesCode, data=zp, log="xy", method="SMA", robust=T, slope.test=1)
h_gw <- sma(gw~SL*SpeciesCode, data=h, log="xy", method="SMA", robust=T, slope.test=1)

p_ga <- sma(ga~SL*SpeciesCode, data=p, log="xy", method="SMA", robust=T, slope.test=2)
b_ga <- sma(ga~SL*SpeciesCode, data=b, log="xy", method="SMA", robust=T, slope.test=2)
zp_ga <- sma(ga~SL*SpeciesCode, data=zp, log="xy", method="SMA", robust=T, slope.test=2)
h_ga <- sma(ga~SL*SpeciesCode, data=h, log="xy", method="SMA", robust=T, slope.test=2)

plot(p_gh)
plot(b_gh)



p_sma$r2
p_sma$n
p_sma$slopetest

all_sma_gh <- groupwise_sma_gh(pento, pento$SpeciesCode)

no_caranx <- subset(p, SpeciesCode != 'CA.MELA')


p_spp_sma <- sma(gh~SL*SpeciesCode, data=p, log="xy", method="SMA", slope.test=1, robust=T)

p_spp_sma <- sma(gh~SL+SpeciesCode, data=b, log='xy')
plot(p_spp_sma, which="residual")
plot(p_spp_sma)
legend("topleft", legend=unique(p$SpeciesCode))

p_sma <- sma(gh~SL, data=p, log="xy", method="SMA", alpha=0.05, robust=T,
             slope.test=1, group=)

caranx <- subset(p, SpeciesCode=="CA.MELA")

caranx_sma <- sma(gh~SL, data=caranx, log="xy", robust=T, slope.test=1)

plot(caranx_sma, which="residual")
plot(caranx_sma, which='qq')

y <- rnorm(12, mean=3.89, sd=0.35)
x <- rnorm(12, mean=6, sd=0.37)

b_sma <- sma(gh~SL, data=b, log="xy", method="SMA", alpha=0.05, robust=T,
             slope.test=1)
zp_sma <- sma(gh~SL, data=zp, log="xy", method="SMA", alpha=0.05, robust=T,
              slope.test=1)
h_sma <- sma(gh~SL, data=h, log="xy", method="SMA", alpha=0.05, robust=T,
             slope.test=1)
c_sma <- sma(gh~SL, data=c, log="xy", method="SMA", alpha=0.05, robust=T,
             slope.test=1)




# way to pass through SMA regression values to plot lines in ggplot2 - HW answer
# from stackoverflow :)
library(ggplot2)
library(smatr)

dat <- data.frame(a=log10(rnorm(50, 30, 10)), b=log10(rnorm(50, 20, 2)))
mod <- lmodel2(a ~ b, data=dat,"interval", "interval", 99)

smatr

reg <- mod$regression.results
names(reg) <- c("method", "intercept", "slope", "angle", "p-value")



sma <- pento_ga$groupsummary
for (i in 1:5) {
  
  sma$from <- pento_ga$from[[i]]
  sma$to[i]   <- pento_ga$to[[i]]
}
sma$yfrom <- 10^(sma$Slope*log10(sma$from) + sma$Int)
sma$yto   <- 10^(sma$Slope*log10(sma$to) + sma$Int)
names(sma)[1] <- 'j_fg'

ggplot(data = pento, aes(x = SL, y = ga, colour=j_fg)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  scale_colour_discrete(name = "Functional \n Group") +
  theme(axis.title.y = element_text(size = 30, vjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30, vjust = 0.1)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(legend.text = element_text(size = 26)) +
  theme(legend.key.height = unit(1.5, "line")) +
  theme(legend.title = element_text(size = 24)) +
  theme(legend.position = c(0.85, 0.24)) +
  theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
  geom_segment(data = sma, aes(x = from, xend = to, y = yfrom, yend = yto)) # +
  #facet_wrap(~j_fg)


sma <- pspp_ga$groupsummary
for (i in 1:7) {
  sma$from[i] <- pspp_ga$from[[i]]
  sma$to[i]   <- pspp_ga$to[[i]]
}
sma$yfrom <- 10^(sma$Slope*log10(sma$from) + sma$Int)
sma$yto   <- 10^(sma$Slope*log10(sma$to) + sma$Int)
names(sma)[1] <- 'SpeciesCode'

ggplot(data = p, aes(x = SL, y = ga, colour=SpeciesCode) +
  geom_point(size = 3 ) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("log(standard length, mm)") +
  ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
  scale_colour_discrete(name = "Species") +
  theme(axis.title.y = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.key.height = unit(1.5, "line")) +
  theme(legend.title = element_text(size = 20)) +
  theme(legend.position = c(0.9, 0.5)) +
  geom_segment(data = sma, aes(x = from, xend = to, y = yfrom, yend = yto), size = 1)
  #facet_wrap(~SpeciesCode)






















