##############################################################################################
################        Preliminary Gape Size - Body Size analysis           ################
##############################################################################################

# set working directory
setwd('/Users/jillian/R_projects/Allometry/')

########### Entering and Cleaning Data ###########

fish <- read.csv('Fish Dimensions_cleaned.csv', header = TRUE, sep = ",") 
# Enter Fish dimension data for gape - body size analysis
# imported as a dataframe

names(fish)[9] <- "Weight" # changes name of column 9 heading

grep("GapeHeight", colnames(fish)) # finds the index of a column with a given name

sapply(fish, mode) # checks the type of each column e.g.,  numeric, factor

fish[is.na(fish$GapeWidth),] #Shows me the row where GapeWidth has an NA
# Logic breakdown: is.na(fish$GapeWidth) specifies the row where this case is TRUE
# So you are selecting this row and all of the columns where this is TRUE [row case, all col]

fish[(fish[13]==445),] # This finds me the row where the max GapeHeight is
# This is an anomalously high number that I have yet to verify from the original data so
# I would like to remove it from the dataframe

fish <- fish[!(fish[13]==445,)] # makes a dataframe excluding the row where GapeHeight = 445 
is.numeric(fish[13]) # tells me if fish[13] is numeric
colnames(fish) # gives me the column names of the dataframe 'fish'

# Finding the anomalously large GH values
# METHOD 1:
max(fish[13]) # gives me the FIRST largest value of GapeHeight
# Below gives me the SECOND largest value of GapeHeight
n <- length(unique(fish$GapeHeight)) # makes n the length of the list of unique GapeHeights
which(fish$GapeHeight == sort(unique(fish$GapeHeight), partial=n-2)[n-2])
# unique selects all of the distinct values of GapeHeight
# sort(unique(fish$GapeHeight)) gives ascending order of GapeHeights
# 'partial' is an option for sort() that specifies the subset of x to be sorted
# sort(unique(fish$GapeHeight), partial=n-1) then indexed by [n-1] 
#### same as x[y]
# Then fish$GapeHeight == sort(unique(fish$GapeHeight), partial=n-1)[n-1]) is a comparison of 
# all of the values in fish$GapeHeight with the the n-1 value of the partial sort array
# which() gives the TRUE indices of a logical object

# METHOD 2:
which.max(fish[,13]) # index of first maximum value
which.max((fish[-(which.max(fish[,13])),])[,13]) # index of second maximum value


(which(fish[13]==max(fish[13]))) # gives me the index of the max fish GapeHeight
big_val <- c(64, 1066) # first and second largest values which have been identified as anomalous

fish <- fish[-big_val,] # new dataframe 'fish' with the 2 anomalously large GH values excluded

which(fish[9]==0)
  fish <- fish[-853,]
#needed to remove line with a zero because cannot do 'method=lm' in ggplot() when there is an INF value (log transformed zero)

which(is.na(fish["Guild"]))

fish[1064, 13] <- 29.1


##### how many of x do I have? #####
a <- table(fish$SpeciesCode)

ch.orna <- which(fish[,4]=="CH.ORNA")
fish(


##############################################################################################
######### Exploratory Analysis #########

wtgape.lm <- lm(log(GapeHeight) ~ log(SL), data=fishSL)
wtgape.stdres <- rstandard(wtgape.lm)

plot(log(fishSL$GapeHeight), wtgape.stdres, ylab="Standardized Residuals", xlab="Weight", main="confused")

> eruption.lm = lm(eruptions ~ waiting, data=faithful) 
> eruption.stdres = rstandard(eruption.lm)
We now plot the standardized residual against the observed values of the variable waiting.

> plot(faithful$waiting, eruption.stdres, 
       +     ylab="Standardized Residuals", 
       +     xlab="Waiting Time", 
       +     main="Old Faithful Eruptions") 
> abline(0, 0)                  # the horizon


#### Testing for Normality #####
with(fish, qqnorm(log(Weight)))
with(fish, qqline(log(Weight)))
with(fish, shapiro.test(log(Weight)))

with(fish, qqnorm(log(GapeHeight)))
with(fish, qqline(log(GapeHeight)))
with(fish, shapiro.test(log(GapeHeight)))

with(fish, qqnorm(GapeHeight))
with(fish, qqline(GapeHeight))
with(fish, shapiro.test(GapeHeight))

with(fish, qqnorm(SL))
with(fish, qqline(SL))
with(fish, shapiro.test(SL))

with(fish, plot(Weight, GapeHeight, xlab = "Weight (g)", ylab = "Gape Height (mm)", pch = 20))
# standard plot() scatterplot of GapeHeight and Weight with small black dots (pch=20)
# This originally showed me two anomalously large GH values which I removed from the current 'fish' above.

# Gives me n for each guild
APCnt <- length(which(fish$Guild=='AP'))
BICnt <- length(which(fish$Guild=='BI'))
DeCnt <- length(which(fish$Guild=='De'))
HeCnt <- length(which(fish$Guild=='He'))
OmCnt <- length(which(fish$Guild=='Om'))
PiCnt <- length(which(fish$Guild=='Pi'))
ZPCnt <- length(which(fish$Guild=='ZP'))
CoCnt <- length(which(fish$Guild=='C'))

# Vector of n of each guild
c(APCnt, BICnt, DeCnt, HeCnt, OmCnt, PiCnt, ZPCnt, CoCnt)


fish[which(fish)


unique(fish$SpeciesCode)
for (i in unique(fish$SpeciesCode)) {
  a <- length(which((fish$SpeciesCode)==i))
  print(c(i, a))
}


########## Learning how to use ggplot2 ############
library("ggplot2")

ggplot(data.frame(x=c(0,5), y=c(0,5)), aes(x, size=2)) +
  opts(axis.text.x = theme_text(size=16, colour="black")) +
  opts(axis.text.y = theme_text(size=16, colour="black")) +
  xlab("log(Weight)") +
  ylab("log(GapeHeight)") +
  opts(axis.title.x = theme_text(size=28, vjust = 0)) +
  opts(axis.title.y = theme_text(size=28, vjust = 0.3)) +
  #stat_function(fun=function(x) {0.7*x + 0.3 }, geom="line", aes(colour="red")) +
  stat_function(fun=function(x) {0.9*x }, geom="line", aes(colour="blue")) +
  stat_function(fun=function(x) {0.2*x }, geom="line", aes(colour="green")) +
  xlim(0,5) +
  ylim(0,5) +
  theme(legend.position="none")



  scale_colour_manual("Function", value=c("blue","red"), breaks=c("square","exp"))

ggplot(data=fishSL, aes(x=log(SL), y=log(GapeHeight), colour=Guild)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  scale_colour_hue(l=70, c=175) +
  annotate("text", x=3.5, y=4.5, label=lm_eqn2(fishSL), parse=TRUE)

ggplot(data=fishSL, aes(x=SL, y=GapeHeight, colour=Guild)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  scale_colour_hue(l=70, c=175) #+
#  annotate("text", x=3.5, y=4.5, label=lm_eqn2(fishSL), parse=TRUE)


ggplot(data=fish, aes(x=log(Weight), y=log(GapeHeight), colour=Guild)) +
  geom_point(shape=19) +
  opts(axis.title.x = theme_text(size=24, vjust=0, colour="black")) +
  opts(axis.title.y = theme_text(size=24, vjust=0.3, colour="black")) +
  opts(axis.text.x = theme_text(size=16, colour="black")) +
  opts(axis.text.y = theme_text(size=16, colour="black")) +
  opts(legend.text = theme_text(size=20, colour="black")) +
  opts(legend.title = theme_blank()) +
  guides(colour = guide_legend(override.aes=list(shape=20, size=7)))

ggplot(data=fish, aes(x=Weight, y=GapeHeight, colour=Guild)) +
  geom_point(shape=19) +
  opts(axis.title.x = theme_text(size=24, vjust=0, colour="black")) +
  opts(axis.title.y = theme_text(size=24, vjust=0.3, colour="black")) +
  opts(axis.text.x = theme_text(size=16, colour="black")) +
  opts(axis.text.y = theme_text(size=16, colour="black")) +
  opts(legend.text = theme_text(size=20, colour="black")) +
  opts(legend.title = theme_blank()) +
  guides(colour = guide_legend(override.aes=list(shape=20, size=7)))


ggplot(data=fishSL, aes(x=log(SL), y=log(GapeHeight), colour=Guild)) +
  geom_point(shape=19) +
  opts(axis.title.x = theme_text(size=24, vjust=0, colour="black")) +
  opts(axis.title.y = theme_text(size=24, vjust=0.3, colour="black")) +
  opts(axis.text.x = theme_text(size=16, colour="black")) +
  opts(axis.text.y = theme_text(size=16, colour="black")) +
  opts(legend.text = theme_text(size=20, colour="black")) +
  opts(legend.title = theme_blank()) +
  guides(colour = guide_legend(override.aes=list(shape=20, size=7)))


opts(legend.text = theme_text(colour = 'red', angle = 45, size = 10, hjust = 3, vjust = 3, face = 'bold')
  
  #geom_smooth(method=lm) +
  scale_colour_hue(l=70, c=175) #+
#  annotate("text", x=1.5, y=4.3, label=lm_eqn(fish), parse=TRUE)

ggplot(data=fish, aes(x=log(Weight), y=log(GapeWidth), colour=Guild)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  scale_colour_hue(l=70, c=175) +
  annotate("text", x=1.5, y=4.3, label=lm_eqn3(fishGW), parse=TRUE)

ggplot(data=fish, aes(x=log(Weight), y=log(GapeHeight), colour=Guild)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  scale_colour_hue(l=70, c=175) +
  geom_abline(intercept = 0, slope = 1, lty=2)


ggplot(fish, aes(x=log(Weight), y=log(GapeHeight))) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method="lm") 
  
ggplot(fish, aes(x=log(Weight), y=log(GapeHeight))) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm, size=2) +
  opts(axis.title.x = theme_text(size=24, vjust=0, colour="black")) +
  opts(axis.title.y = theme_text(size=24, vjust=0.3, colour="black")) +
  opts(axis.text.x = theme_text(colour="black", size=16)) +
  opts(axis.text.y = theme_text(colour="black", size=16)) +
  opts(legend.text = theme_text(size=20, colour="black")) +
  opts(strip.text.x = theme_text(colour="black", size=16)) +
  opts(legend.title = theme_blank()) +  
  guides(colour = guide_legend(override.aes=list(size=3))) + 
  facet_wrap(~ Guild, nrow=2)
     
  opts(strip.text.x = theme_text(colour = 'red', angle = 45, size = 10, hjust = 0.5, vjust = 0.5, face = 'bold'))
  
  facet_wrap(~ Guild, nrow=2)

lm(log(GapeHeight) ~ log(Weight), fish)




# The following function is from: http://stackoverflow.com/questions/7549694/ggplot2-adding-regression-line-equation-and-r2-on-graph
# It generates the regression equation that can be used to annotate a ggplot.
lm_eqn = function(x){
  m = lm(log(GapeHeight) ~ log(Weight), x)
  eqn <- substitute(italic(y) == a%.% italic(x) + b*","~~italic(r)^2~"="~r2,
                    list(a = format(coef(m)[2], digits = 3),
                         b = format(coef(m)[1], digits = 3),
                         r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eqn))
}

which(fish$SL == 0)
fishSL <- fish[c(-321, -344),]

lm_eqn2 = function(x){
  m = lm(log(GapeHeight) ~ log(SL), x)
  eqn <- substitute(italic(y) == a%.% italic(x) + b*","~~italic(r)^2~"="~r2,
                    list(a = format(coef(m)[2], digits = 3),
                         b = format(coef(m)[1], digits = 3),
                         r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eqn))
}

which(is.na(fish$GapeWidth)==TRUE)
fishGW <- fish[-884,]

lm_eqn3 = function(x){
  m = lm(log(GapeWidth) ~ log(SL), x)
  eqn <- substitute(italic(y) == a%.% italic(x) + b*","~~italic(r)^2~"="~r2,
                    list(a = format(coef(m)[2], digits = 3),
                         b = format(coef(m)[1], digits = 3),
                         r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eqn))
}


qplot(fish$Guild, data=fish, geom='bar', colour=I("black"), fill=SpeciesCode)
qplot(fish$Guild, data=fish, geom='bar', colour=I("black"), fill=SizeClass) # is a histogram of fish frequency by guild, with colours representing the number of samples we have in given size classes

qplot(log(Weight), log(GapeHeight), data=fish, colour=SizeClass)

qplot(log(Weight), log(GapeHeight/Weight), data=fish, colour=Guild)

qplot(log(Weight), log(GapeHeight), data=fish, colour=Guild)
qplot(log(Weight), log(GapeHeight), data=fish, colour=SpeciesCode)

qplot(aes(x=log(Weight), y=log(GapeHeight) #+
  geom_point(shape=1)# Use hollow circles
  geom_smooth(method=lm)


qplot(log(Weight), log(GapeHeight), data=fish, colour=SizeClass)

qplot(log(Weight), log(GapeHeight), data=fish, colour=SpeciesCode)
