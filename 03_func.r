library("ggplot2")
library('plyr')
library('gridExtra') # allows me to make multiplots
library(nlme)

write_lme_gen <- function(df, y){
  m = lm(log(y) ~ log(SLMM), df);
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
  as.character(as.expression(eq))
}


count_spp <- function(df) {
  ddply(.data = df, .(SpeciesCode), summarize, 
        len = length(SpeciesCode),
        n = paste("n ==", len)
  )
}


fg_lm <- function(df, gape) {
    lm <- lm(log(gape)~log(SLMM), data=df)
    data.frame(int    = coefficients(lm)[1],
               slope  = coefficients(lm)[2],
               rsq    = summary(lm)$r.squared,
               se     = summary(lm)$coefficients[2,2],
               p_val  = summary(lm)$coef[2,4],
               lw_conf_slp = confint(lm, level = 0.95)[2,1],
               up_conf_slp = confint(lm, level = 0.95)[2,2],
               lw_conf_int = confint(lm, level = 0.95)[1,1],
               up_conf_int = confint(lm, level = 0.95)[1,2]
               )
}

groupwise_lm_gh <- function(df, variable) {
  lm  <- with(data=df, ddply(df, .(variable), function(z) {
    t <- lm(log(gh)~log(SLMM), data=z)
    data.frame(int    = coefficients(t)[1],
               slope  = coefficients(t)[2],
               rsq    = summary(t)$r.squared,
               se     = summary(t)$coefficients[2,2],
               p_val  = summary(t)$coef[2,4],
               lw_conf_slp = confint(t, level = 0.95)[2,1],
               up_conf_slp = confint(t, level = 0.95)[2,2],
               lw_conf_int = confint(t, level = 0.95)[1,1],
               up_conf_int = confint(t, level = 0.95)[1,2])
  }))
}

groupwise_lm_gw <- function(df, variable) {
  lm  <- with(data=df, ddply(df, .(variable), function(z) {
    t <- lm(log(gw)~log(SLMM), data=z)
    data.frame(int    = coefficients(t)[1],
               slope  = coefficients(t)[2],
               rsq    = summary(t)$r.squared,
               se     = summary(t)$coefficients[2,2],
               p_val  = summary(t)$coef[2,4],
               lw_conf_slp = confint(t, level = 0.95)[2,1],
               up_conf_slp = confint(t, level = 0.95)[2,2],
               lw_conf_int = confint(t, level = 0.95)[1,1],
               up_conf_int = confint(t, level = 0.95)[1,2])
  }))
}

groupwise_lm_ga <- function(df, variable) {
  lm  <- with(data=df, ddply(df, .(variable), function(z) {
    t <- lm(log(ga)~log(SLMM), data=z)
    data.frame(int    = coefficients(t)[1],
               slope  = coefficients(t)[2],
               rsq    = summary(t)$r.squared,
               se     = summary(t)$coefficients[2,2],
               p_val  = summary(t)$coef[2,4],
               lw_conf_slp = confint(t, level = 0.95)[2,1],
               up_conf_slp = confint(t, level = 0.95)[2,2],
               lw_conf_int = confint(t, level = 0.95)[1,1],
               up_conf_int = confint(t, level = 0.95)[1,2])
  }))
}


awesome <- function(lm) {
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

# groupwise_lm_gh <- function(df, variable) {
#   lm  <- with(data=df, ddply(df, .(variable), function(z) {
#     t <- lm(log(gh)~log(SLMM), data=z)
#     data.frame(int   = coefficients(t)[1],
#                slope = coefficients(t)[2],
#                rsq   = summary(t)$r.squared,
#                se    = summary(t)$coefficients[2,2],
#                p_val = summary(t)$coef[2,4])
#   }))
# }
# 
# groupwise_lm_gw <- function(df, variable) {
#   lm  <- with(data=df, ddply(df, .(variable), function(z) {
#     t <- lm(log(gw)~log(SLMM), data=z)
#     data.frame(int   = coefficients(t)[1],
#                slope = coefficients(t)[2],
#                rsq   = summary(t)$r.squared,
#                se    = summary(t)$coefficients[2,2],
#                p_val = summary(t)$coef[2,4])
#   }))
# }
# 
# groupwise_lm_ga <- function(df, variable) {
#   lm  <- with(data=df, ddply(df, .(variable), function(z) {
#     t <- lm(log(ga)~log(SLMM), data=z)
#     data.frame(int   = coefficients(t)[1],
#                slope = coefficients(t)[2],
#                rsq   = summary(t)$r.squared,
#                se    = summary(t)$coefficients[2,2],
#                p_val = summary(t)$coef[2,4])
#   }))
# }
# 
# call_master_df <- function(lm) {
#     coefs <- ldply(lm, coef)
#     rsq   <- function(lm) c("rsqd"  = summary(lm)$r.squared)
#     se    <- function(lm) c("SE"    = summary(lm)$coefficients[2,2])
#     p_val <- function(lm) c("p_val" = summary(lm)$coefficients[2,4])
#     spp.n <- ddply(.data = p, .(SpeciesCode), summarize, 
#                    n     = paste("n ==", length(SpeciesCode))
#                    )
#     dfs <- list(coefs,
#                 rsqs   <- ldply(lm, rsq),
#                 ses    <- ldply(lm, se),
#                 p_vals <- ldply(lm, p_val),
#                 spp.n)
#     all <- join_all(dfs, by="SpeciesCode")
# }

# Function that generates data

# Example dataframe produced by using above groupwise function
# test <- groupwise_lm_gh(fish, fish$SpeciesCode)


write_lme_groups <- function(summ_df, variable) {
  df <- summ_df
  m = matrix(data=NA, nrow=0, ncol=2)
  len <- length(variable)
  for (i in (1:len)) {
    l <- list(slp = format(df[[3]][i], digits=2),
              int = format(df[[2]][i], digits=2), 
              r2 = format(df[[4]][i], digits=2)
    )
    if (l$int >= 0) {
      eqn <- substitute(log~italic(y) ==
                          slp%.% log~italic(x) + int*~~italic(r)^2~"="~r2, l)
    }
    
    else {
      l <- list(slp = format(df[[3]][i], digits=2),
                int = format(abs(df[[2]][i]), digits=2), 
                r2 = format(df[[4]][i], digits=2)
      )
      eqn <- substitute(log~italic(y) ==
                          slp%.% log~italic(x) - int*~~italic(r)^2~"="~r2, l)
    }
    #browser()
    lm_eq <- as.character(as.expression(eqn)) 
    m <- rbind(m, c(as.character(df[[1]][i]), lm_eq))
    #m <- rbind(m, c(as.character(df[i,1]), lm_eq))
  }
  m <- as.data.frame(m)
}


################################################################################
#############             Predator - Prey Size Functions            ############
################################################################################
# groupwise_rq <- function(df, variable) {
#   rq  <- ddply(df, .(variable), function(z) {
#     r <- rq(psize ~ sl, tau = c(0.10, 0.90), data = z) 
#     })
# }
# 
# rq(formula, tau=.5, data, subset, weights, na.action,
#    method="br", model = TRUE, contrasts, ...) 
# 
# groupwise_lm_gw <- function(df, variable) {
#   lm  <- with(data=df, ddply(df, .(variable), function(z) {
#     t <- lm(log(gw)~log(SLMM), data=z)
#     data.frame(int   = coefficients(t)[1],
#                slope = coefficients(t)[2],
#                rsq   = summary(t)$r.squared,
#                se    = summary(t)$coefficients[2,2],
#                p_val = summary(t)$coef[2,4])
#   }))
# }























