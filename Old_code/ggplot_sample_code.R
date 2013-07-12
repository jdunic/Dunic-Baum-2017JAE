########## Learning how to use ggplot2 ############
library("ggplot2")

# gives a scatter plot of logWeight and logGape. I used opts() here to make the figures
# presentable for the big screen (e.g., increasing axes titles, label, and legend text size)

ggplot(data=fish, aes(x=log(Weight), y=log(GapeHeight), colour=Guild)) +
  geom_point(shape=19) +
  opts(axis.title.x = theme_text(size=24, vjust=0, colour="black")) +
  opts(axis.title.y = theme_text(size=24, vjust=0.3, colour="black")) +
  opts(axis.text.x = theme_text(size=16, colour="black")) +
  opts(axis.text.y = theme_text(size=16, colour="black")) +
  opts(legend.text = theme_text(size=20, colour="black")) +
  opts(legend.title = theme_blank()) +
  guides(colour = guide_legend(override.aes=list(shape=20, size=7)))



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


# Example of scatterplot with regression eqn written on the plot
ggplot(data=fish, aes(x=log(Weight), y=log(GapeHeight), colour=Guild)) +
  geom_point(shape=20) +
  geom_smooth(method=lm) +
  scale_colour_hue(l=70, c=175) +
  annotate("text", x=3.5, y=4.5, label=lm_eqn(fish), parse=TRUE)  ## here is where i use the lm_eqn label


# An example of where I used facets to show log(Weight)~log(Gape) relationships for each of the different
# functional groups.   
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
  facet_wrap(~ Guild, nrow=2)  # gives me two rows of






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
