# Prey Size analysis

colnames(prey1)

data.frame()

df.n <- ddply(.data=prey2, .(species), summarize, n=paste("n ==", length(species)))

# all species
ggplot(prey2, aes(x=sl, y=psize, colour=ptype)) +
  geom_point(shape = 20)

ggplot(prey2, aes(x=sl, y=psize)) +
  geom_point(aes(colour = ptype)) + 
  geom_text(data=df.n, aes(x=200, y=200, label=n), parse=TRUE) +
  facet_wrap(~ species)

df.n <- ddply(.data=prey2, .(fg), summarize, n=paste("n ==", length(fg)))

ggplot(prey2, aes(x=sl, y=psize)) +
  geom_point(aes(colour=ptype)) +
  geom_text(data=df.n, aes(x=200, y=200, label=n), parse=TRUE) +
  facet_wrap(~ fg)  

# Counting unique specimens for each group: p, bi, zp to look at 
# how that might affect quantile selection ####
goatfish <- subset(prey2, species == 'PA.INSU')
length(unique(goatfish$specimen))

b <- subset(prey2, fg == 'BI')
View(b)

p <- subset(prey2, fg == 'Pi')
View(p)

length(unique(subset(p, species == 'CE.UROD')$specimen))

for (i in unique(p$species)) {
  a <- length(unique(subset(p, species == i)$specimen))
  print(c(i, a))
}

zp <- subset(prey2, fg == 'ZP')
View(zp)

for (i in unique(zp$species)) {
  a <- length(unique(subset(zp, species == i)$specimen))
  print(c(i, a))
}

# Actually doing the quantile regression:
# Have used the 10th and 90th percentiles because of sample size
# and following the general guidelines of n > (10/q) given in 
# Scharf et al, 1998 - Ecology

library('quantreg')
library('ggplot2')

# prey3 is a data.frame with just Pi, BI, ZP.
df.n <- ddply(.data=prey3, .(fg), summarize, n=paste("n ==", length(fg)))

ggplot(prey3, aes(x = sl, y = psize)) +
  geom_point(aes(colour = species)) +
  geom_text(data = df.n, aes( x= 200, y = 200, label = n), parse = TRUE) +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.5, 0.90), method = "rq") +
  facet_wrap(~ fg)

df.n <- ddply(.data=prey3, .(fg), summarize, n=paste("n ==", length(fg)))

ggplot(na.omit(prey3), aes(x = sl, y = psize)) +
  geom_point(aes(colour = species)) +
  geom_text(data = df.n, aes( x= 200, y = 200, label = n), parse = TRUE) +
  geom_smooth(method = "lm") +
  stat_quantile(geom = "quantile", quantiles = c(0.10, 0.90), method = "rq") +
  facet_wrap(~ fg)

bi_qr <- rq(data = subset(prey3, prey3$fg == 'BI'), psize ~ sl, 
            tau = c(0.10, 0.50, 0.90))

p <- subset(prey3, prey3$fg == 'Pi')
b <- subset(prey3, prey3$fg == 'BI')

is.na(p$sl)

lm.p <- lm(data = p, formula = psize ~ sl, na.action=na.exclude)
p.res <- resid(lm.p)
resid <- with(p, plot(sl, residuals(lm.p)))
abline(0,0, col="red")

lm.b <- lm(data = b, formula = psize ~ sl, na.action=na.exclude)
b.res <- resid(lm.b)
resid <- with(b, plot(sl, b.res))
abline(0,0, col="red")

summary(bi_qr)
bi_qr[[3]]

groupwise_rq <- function(df, variable) {
  #browser()
  rq  <- ddply(df, .(variable), function(z) {
    #browser()
    r <- qr_func(z) 
  })
}
groupwise_rq(prey3, prey3$fg)

rq  <- ddply(prey3, .(prey3$fg), function(z) {
  r <- rq(psize ~ sl, tau = c(0.10, 0.90), data = prey3) 
  })


