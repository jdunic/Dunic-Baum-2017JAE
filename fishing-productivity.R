library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)

pento <- 
separate(data = pento, col = Region, into = c('productivity', 'fishing')) %>% 
  mutate(fishing = factor(fishing, levels = c('HF', 'MF', 'LF')), 
         productivity = as.factor(productivity), 
         log_gh = log10(gh), 
         log_gw = log10(gw), 
         log_SL = log10(SL), 
         log_m = log10(wt))

pento %>% 
  filter(SpeciesCode != 'CH.ORNA') %>% 
  group_by(j_fg, fishing, SpeciesCode) %>% 
  dplyr::summarise(count = n()) %>%
ggplot(data = ., aes(x = SpeciesCode, y = count, fill = j_fg)) + 
  geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 45)) + 
  facet_grid(fishing ~ .)

pento %>% 
  filter(SpeciesCode != 'CH.ORNA') %>%  
  group_by(j_fg, productivity, SpeciesCode) %>% 
  dplyr::summarise(count = n()) %>% 
ggplot(data = ., aes(x = SpeciesCode, y = count, fill = j_fg)) + 
  geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 45)) + 
  facet_grid(productivity ~ .)

test <- lmer(log_gh ~ log_m + log_m:fishing + 0 + (1 | SpeciesCode), 
                   data = filter(pento, j_fg == 'Pi'))
summary(test)

test <- lmer(log_gh ~ log_m + log_m:fishing + 0 + (1 | SpeciesCode), 
                   data = filter(pento, j_fg == 'BI'))
summary(test)

test <- lmer(log_gh ~ log_m + log_m:fishing + 0 + (1 | SpeciesCode), 
                   data = filter(pento, j_fg == 'ZP'))
summary(test)

test <- lmer(log_gh ~ log_m + log_m:fishing + 0 + (1 | SpeciesCode), 
                   data = filter(pento, j_fg == 'He'))
summary(test)

test <- lmer(log_gh ~ log_m + log_m:productivity + 0 + (1 | SpeciesCode), 
                   data = filter(pento, j_fg == 'Pi'))
summary(test)

test <- lmer(log_gh ~ log_m + log_m:productivity + 0 + (1 | SpeciesCode), 
                   data = filter(pento, j_fg == 'BI'))
summary(test)

test <- lmer(log_gh ~ log_m + log_m:productivity + 0 + (1 | SpeciesCode), 
                   data = filter(pento, j_fg == 'ZP'))
summary(test)

test <- lmer(log_gh ~ log_m + log_m:productivity + 0 + (1 | SpeciesCode), 
                   data = filter(pento, j_fg == 'He'))
summary(test)

test <- 
  filter(pento, j_fg == 'Pi') %>%
  as.data.frame(.) %>% 
sma(gh~SL * fishing, data = ., log = "xy", method = "SMA", robust = T,
      slope.test = 2, multcomp = T, multcompmethod = "adjusted") %>% 

p_site_summ <- mk_spp_summary(test, length(test$groups), grouping = TRUE)
p_site_graph_df <- mk_smaSPP_graph_df(p_site_summ, length(test$groups), "fishing")


filter(pento, j_fg == 'BI') %>% 
  as.data.frame(.) %>% 
sma(gh~SL * fishing, data = ., log = "xy", method = "SMA", robust = T,
      slope.test = 2, multcomp = T, multcompmethod = "adjusted") %>% 
mk_SMAplot(df_points = ., df_lines = p_site_graph_df, facets = FALSE, 
         x = "SL", gapeType = "gh", grouping = "Region", labels = "None", 
         axis_labels = FALSE) +
  theme(axis.title = element_text(size = 10)) + 
  theme(axis.title.x = element_text(vjust = -0.5)) +
  theme(axis.title.y = element_text(vjust = 0.4))


filter(pento, j_fg == 'ZP') %>% 
  as.data.frame(.) %>% 
sma(gh~SL * fishing, data = ., log = "xy", method = "SMA", robust = T,
      slope.test = 2, multcomp = T, multcompmethod = "adjusted") %>% 
  summary(.)

separate(data = pento, col = Region, into = c('productivity', 'fishing')) %>% 
  filter(j_fg == 'He') %>% 
  mutate(fishing = as.factor(fishing)) %>% 
  as.data.frame(.) %>% 
sma(gh~SL * fishing, data = ., log = "xy", method = "SMA", robust = T,
      slope.test = 2, multcomp = T, multcompmethod = "adjusted") %>% 
  summary(.)


separate(data = pento, col = Region, into = c('productivity', 'fishing')) %>% 
  filter(j_fg == 'Pi') %>% 
  mutate(productivity = as.factor(productivity)) %>% 
  as.data.frame(.) %>% 
sma(gh~SL * productivity, data = ., log = "xy", method = "SMA", robust = T,
      slope.test = 2, multcomp = T, multcompmethod = "adjusted")

separate(data = pento, col = Region, into = c('productivity', 'fishing')) %>% 
  filter(j_fg == 'BI') %>% 
  mutate(productivity = as.factor(productivity)) %>% 
  as.data.frame(.) %>% 
sma(gh~SL * productivity, data = ., log = "xy", method = "SMA", robust = T,
      slope.test = 2, multcomp = T, multcompmethod = "adjusted") %>% 
  summary(.)

separate(data = pento, col = Region, into = c('productivity', 'fishing')) %>% 
  filter(j_fg == 'ZP') %>% 
  mutate(productivity = as.factor(productivity)) %>% 
  as.data.frame(.) %>% 
sma(gh~SL * productivity, data = ., log = "xy", method = "SMA", robust = T,
      slope.test = 2, multcomp = T, multcompmethod = "adjusted") %>% 
  #summary(.)
  plot(., col = c('red', 'blue'))
  legend('topleft', c('HP', 'LP'), col = c('red', 'blue'), pch=rep(1,2))


separate(data = pento, col = Region, into = c('productivity', 'fishing')) %>% 
  filter(j_fg == 'ZP') %>% 
  mutate(fishing = as.factor(fishing)) %>% 
  as.data.frame(.) %>% 
sma(gh~SL * fishing, data = ., log = "xy", method = "SMA", robust = T,
      slope.test = 2, multcomp = T, multcompmethod = "adjusted") %>% 
  #summary(.)
  plot(., col = c('red', 'blue', 'orange'))
  legend('topleft', c('HF', 'LF', 'MF'), col = c('red', 'blue', 'orange'), pch=rep(1, 2, 3))

separate(data = pento, col = Region, into = c('productivity', 'fishing')) %>% 
  filter(j_fg == 'ZP') %>% 
  mutate(fishing = as.factor(fishing), 
         productivity = as.factor(productivity), 
         log_gh = log(gh), 
         log_SL = log(SL), 
         log_m = log(wt)) %>% 
  as.data.frame(.) %>% 
lm(log_gh ~ log_SL + log_SL:fishing + log_SL:productivity +  0, data = .) %>% 
  summary(.)


separate(data = pento, col = Region, into = c('productivity', 'fishing')) %>% 
  filter(j_fg == 'He') %>% 
  mutate(productivity = as.factor(productivity)) %>% 
  as.data.frame(.) %>% 
sma(gh~SL * productivity, data = ., log = "xy", method = "SMA", robust = T,
      slope.test = 2, multcomp = T, multcompmethod = "adjusted") %>% 
  summary(.)


separate(data = pento, col = Region, into = c('productivity', 'fishing')) %>% 
  filter(j_fg == 'ZP') %>% 
  distinct(fishing) %>% 
  select(fishing)

separate(data = pento, col = Region, into = c('productivity', 'fishing')) %>% 
  filter(j_fg == 'BI') %>% 
  distinct(fishing) %>% 
  select(fishing)


# Looking at fishing and productivity as moderators
filter(pento, j_fg == 'Pi') %>% 
  as.data.frame(.) %>% 
lm(log_gh ~ log_SL * (fishing + productivity) +  0, data = .) %>% 
  car::Anova(.)

filter(pento, j_fg == 'BI') %>% 
  as.data.frame(.) %>% 
lm(log_gh ~ log_SL * (fishing + productivity) +  0, data = .) %>% 
  car::Anova(.)

filter(pento, j_fg == 'ZP') %>% 
  as.data.frame(.) %>% 
lm(log_gh ~ log_SL * (fishing + productivity) +  0, data = .) %>% 
  car::Anova(.)

filter(pento, j_fg == 'He') %>% 
  as.data.frame(.) %>% 
lm(log_gh ~ log_SL * (fishing + productivity) +  0, data = .) %>% 
  car::Anova(.)


# gape width ~ SL
#-------------------------------------------------------------------------------

filter(pento, j_fg == 'Pi') %>% 
  as.data.frame(.) %>% 
lm(log_gw ~ log_SL * (fishing + productivity) +  0, data = .) %>% 
  car::Anova(.)

filter(pento, j_fg == 'BI') %>% 
  as.data.frame(.) %>% 
lm(log_gw ~ log_SL * (fishing + productivity) +  0, data = .) %>% 
  car::Anova(.)

filter(pento, j_fg == 'ZP') %>% 
  as.data.frame(.) %>% 
lm(log_gw ~ log_SL + log_SL:fishing + log_SL:productivity +  0, data = .) %>% 
  car::Anova(.)

filter(pento, j_fg == 'He') %>% 
  as.data.frame(.) %>% 
lm(log_gw ~ log_SL + log_SL:fishing + log_SL:productivity +  0, data = .) %>% 
  car::Anova(.)


# gape height ~ weight
#-------------------------------------------------------------------------------
fit <- 
lme(log_gh ~ log_m * (fishing + productivity), data = filter(pento, j_fg == 'Pi'))
summary(fit)
car::Anova(fit)

fit <- 
lm(log_gh ~ log_m * (fishing + productivity), data = filter(pento, j_fg == 'BI'))
summary(fit)
car::Anova(fit)

fit <- 
lm(log_gh ~ log_m * (fishing + productivity), data = filter(pento, j_fg == 'ZP'))
summary(fit)
car::Anova(fit)


filter(pento, j_fg == 'BI') %>% 
  as.data.frame(.) %>% 
lm(log_gh ~ log_m * (fishing + productivity) +  0, data = .) %>% 
  car::Anova(.)

filter(pento, j_fg == 'ZP') %>% 
  as.data.frame(.) %>% 
lm(log_gh ~ log_m * (fishing + productivity) +  0, data = .) %>% 
  car::Anova(.)

filter(pento, j_fg == 'He') %>% 
  as.data.frame(.) %>% 
lm(log_gh ~ log_m * (fishing + productivity) +  0, data = .) %>% 
  car::Anova(.)

# gape width ~ mass
#-------------------------------------------------------------------------------
filter(pento, j_fg == 'Pi') %>% 
  as.data.frame(.) %>% 
lm(log_gw ~ log_m : (fishing + productivity) +  0, data = .) %>% 
  car::Anova(.)

filter(pento, j_fg == 'BI') %>% 
  as.data.frame(.) %>% 
lm(log_gw ~ log_m * (fishing + productivity) +  0, data = .) %>% 
  car::Anova(.)

filter(pento, j_fg == 'ZP') %>% 
  as.data.frame(.) %>% 
lm(log_gw ~ log_m * (fishing + productivity) +  0, data = .) %>% 
  car::Anova(.)

filter(pento, j_fg == 'He') %>% 
  as.data.frame(.) %>% 
lm(log_gw ~ log_m * (fishing + productivity) +  0, data = .) %>% 
  car::Anova(.)


test <- sma(log_gw ~ log_m * fishing, data = filter(pento, j_fg == 'Pi'))
plot(test, col = c('blue', 'red', 'green'))
legend(2, 4, as.character(unique(filter(pento, j_fg == 'Pi')$fishing)), col=colours, pch=rep(1, 1, 1))

test <- sma(log_gw ~ log_m * fishing, data = filter(pento, j_fg == 'BI'))
plot(test, col = c('blue', 'red', 'green'))
legend(2, 4, as.character(unique(filter(pento, j_fg == 'BI')$fishing)), col=colours, pch=rep(1, 1, 1))


test <- sma(log_gw ~ log_m * fishing, data = filter(pento, j_fg == 'ZP'))
plot(test, col = c('blue', 'red', 'green'))
legend(2, 4, as.character(unique(filter(pento, j_fg == 'ZP')$fishing)), col=colours, pch=rep(1, 1, 1))

test <- sma(log_gw ~ log_m * fishing, data = filter(pento, j_fg == 'He'))
plot(test, col = c('blue', 'red', 'green'))
legend(2, 4, as.character(unique(filter(pento, j_fg == 'He')$fishing)), col=colours, pch=rep(1, 1, 1))

