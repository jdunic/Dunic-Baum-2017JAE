# ggplot default colour list for n colours
gg_colour_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

################################################################################
#############             Graph equations and formatting            ############
################################################################################
# extract ggplot legend to stick beside gridExtra plots as desired:
extract_legend <- function(a.gplot){
  tmp    <- ggplot_gtable(ggplot_build(a.gplot))
  leg    <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

write_lme_gen <- function(df, y){
  m = lm(log(y) ~ log(SL), df);
  l <- list(a = format(coef(m)[1], digits = 2), 
            b = format(coef(m)[2], digits = 2), 
            r2 = format(summary(m)$r.squared, digits = 3)
  )
  
  if (l$a >= 0) {
    eq <- substitute(italic(y) == b %.% italic(x) + a*","~~italic(r)^2~"="~r2, l) 
  } else {
      l <- list(a = format(abs(coef(m)[1]), digits = 2), 
                b = format(coef(m)[2], digits = 2), 
                r2 = format(summary(m)$r.squared, digits = 3)
    )
    eq <- substitute(italic(y) == b %.% italic(x) - a*","~~italic(r)^2~"="~r2, l)
  }
  as.character(as.expression(eq))
}


count_spp <- function(df) {
  ddply(.data = df, .(SpeciesCode), summarize, 
        len = length(SpeciesCode),
        n = paste("n ==", len)
  )
}


################################################################################
#############                    Linear Models                      ############
################################################################################

fg_lm <- function(df, gape) {
    lm <- lm(log(gape)~log(SL), data=df)
    data.frame(int   = coefficients(lm)[1],
               slope = coefficients(lm)[2],
               rsq   = summary(lm)$r.squared,
               se    = summary(lm)$coefficients[2,2],
               p_val = summary(lm)$coef[2,4],
               lw_conf_slp = confint(lm, level = 0.95)[2,1],
               up_conf_slp = confint(lm, level = 0.95)[2,2],
               lw_conf_int = confint(lm, level = 0.95)[1,1],
               up_conf_int = confint(lm, level = 0.95)[1,2]
               )
}

groupwise_lm_gh <- function(df, variable) {
  lm  <- with(data=df, ddply(df, .(variable), function(z) {
    t <- lm(log(gh)~log(SL), data=z)
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
    t <- lm(log(gw)~log(SL), data=z)
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
    t <- lm(log(ga)~log(SL), data=z)
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

################################################################################
#############                    SMA Functions                      ############
################################################################################
mk_sma_df <- function(t) {
  data.frame(elevation  = t$coef[[1]][1,1],
             lw_ci_elev = t$coef[[1]][1,2],
             up_ci_elev = t$coef[[1]][1,3],
             slope      = t$coef[[1]][2,1],
             lw_ci_slp  = t$coef[[1]][2,2],
             up_ci_slp  = t$coef[[1]][2,3],
             r2         = t$r2[[1]],
             n          = t$n[[1]],
             pval       = t$pval[[1]]
  )
}

mk_sma_graph_df <- function(sma_summary_df) {
  sma_graph_df <- data.frame(group=numeric(), slp=numeric(), int=numeric(), 
                             from=numeric(), to=numeric(), yfrom=numeric(), 
                             yto=numeric()
  )
  for (i in 1:5) {
    from  <- sma_summary_df[10, i]
    to    <- sma_summary_df[11, i]
    slp   <- sma_summary_df[3, i]
    int   <- sma_summary_df[1, i]
    yfrom <- 10^(slp*log10(from) + int)
    yto   <- 10^(slp*log10(to) + int)
    group <- colnames(sma_summary_df)[i]
    
    row <- t(c(group=group, slp=slp, int=int, from=from, to=to, yfrom=yfrom,
               yto=yto))
    sma_graph_df <- rbind(sma_graph_df, row)
  }
  sma_graph_df[, 2] <- as.numeric(as.character(sma_graph_df[, 2]))
  sma_graph_df[, 3] <- as.numeric(as.character(sma_graph_df[, 3]))
  sma_graph_df[, 4] <- as.numeric(as.character(sma_graph_df[, 4]))
  sma_graph_df[, 5] <- as.numeric(as.character(sma_graph_df[, 5]))
  sma_graph_df[, 6] <- as.numeric(as.character(sma_graph_df[, 6]))
  sma_graph_df[, 7] <- as.numeric(as.character(sma_graph_df[, 7]))
  return(sma_graph_df)
}


mk_sma_summary <- function(sma_object, group="column_name", grouping=F) {
  rows <- c('elev', 'slp_test', 'slope', 'lower_ci', 'upper_ci',
            'slp_p_value', 'xy_r^2', 'xy_corr_p_value', 'n', 'from', 'to')
  if (grouping==F) {
    elev = coef(sma_object)[[1]]
    slp_test = sma_object$slopetest[[1]][[4]]
    slope  = sma_object$slopetest[[1]][[5]]
    lower  = sma_object$slopetest[[1]][6][[1]][[1]]
    upper  = sma_object$slopetest[[1]][6][[1]][[2]]
    slp_p_val  = sma_object$slopetest[[1]][[3]]
    xy_r2  = sma_object$r2[[1]]
    xy_cor = sma_object$pval[[1]]
    n = sma_object$n[[1]]
    from = sma_object$from[[1]]
    to   = sma_object$to[[1]]
  } else {
    for (i in 1:length(sma_object$groups)) {
      elev = sma_object$coef[[i]][[1]][1]
      slp_test = sma_object$slopetest[[i]]$test.value
      slope = sma_object$slopetest[[i]]$b
      lower = sma_object$slopetest[[i]]$ci[1, 1]
      upper = sma_object$slopetest[[i]]$ci[1, 2]
      slp_p_val = sma_object$slopetest[[i]]$p
      xy_r2 = sma_object$r
      xy_cor = sma_object$pval
      n = sma_object$n
      from = sma_object$from
      to = sma_object$to
    }
  }
  
  columns <- c(elev, slp_test, slope, lower, upper, slp_p_val, xy_r2, xy_cor, n,
               from, to)
  
  sma_df <- data.frame(columns, row.names=rows)
  #names(sma_df) <- sma_object$groups
  #sma_df <- format(sma_df[1], digits=3, sci=F)
  return(sma_df)
}

mk_spp_summary <- function(sma_object, num_spp) {

  sma_df <- data.frame(group=character(), elev=numeric(), slp_test=numeric(), 
                      slope=numeric(), upper=numeric(), lower=numeric(), 
                      slp_p_val=numeric(), xy_r2=numeric(), xy_cor=numeric(), 
                      n=numeric(), from=numeric(), to=numeric()
  )
  
  for (i in 1:num_spp) {
    #spp = as.factor(names(sma_object[[1]]))
    #spp = as.character(attr(sma_object[i], which="split_labels")[[i]][[1]])
    elev = coef(sma_object[[i]])[[1]]
    slp_test = sma_object[[i]]$slopetest[[1]][[4]]
    slope  = sma_object[[i]]$slopetest[[1]][[5]]
    lower  = sma_object[[i]]$slopetest[[1]][6][[1]][[1]]
    upper  = sma_object[[i]]$slopetest[[1]][6][[1]][[2]]
    slp_p_val  = sma_object[[i]]$slopetest[[1]][[3]]
    xy_r2  = sma_object[[i]]$r2[[1]]
    xy_cor = sma_object[[i]]$pval[[1]]
    n = sma_object[[i]]$n[[1]]
    from = sma_object[[i]]$from[[1]]
    to   = sma_object[[i]]$to[[1]]
    group = sma_object$groups[i]  
    
    row <- c(group, elev, slp_test, slope, lower, upper, slp_p_val, xy_r2, xy_cor, 
             n, from, to)
    sma_df <- rbind(sma_df, row)
  }

  columns <- c("elev", "slp_test", "slope", "lower", "upper", "slp_p_val", 
               "xy_r2", "xy_cor", "n", "from", "to")
  
  #sma_df <- data.frame(columns, row.names=rows)
  names(sma_df) <- columns
  #sma_df <- format(sma_df[1], digits=3, sci=F)
  return(sma_df)
}

mk_smaSPP_graph_df <- function(sma_summary_df, num_spp) {
  sma_graph_df <- data.frame(group=character(), slp=numeric(), int=numeric(), 
                             from=numeric(), to=numeric(), yfrom=numeric(), 
                             yto=numeric()
  )
  for (i in 1:num_spp) {
    from  <- sma_summary_df[i, 11]
    to    <- sma_summary_df[i, 12]
    slp   <- sma_summary_df[i, 4]
    int   <- sma_summary_df[i, 2]
    yfrom <- 10^(slp*log10(from) + int)
    yto   <- 10^(slp*log10(to) + int)
    group <- as.character(sma_summary_df[i, 1])
    
    row <- t(c(group=group, slp=slp, int=int, from=from, to=to, yfrom=yfrom,
               yto=yto)
             )
    sma_graph_df <- rbind(sma_graph_df, row)
  }
  sma_graph_df[, 2] <- as.numeric(as.character(sma_graph_df[, 2]))
  sma_graph_df[, 3] <- as.numeric(as.character(sma_graph_df[, 3]))
  sma_graph_df[, 4] <- as.numeric(as.character(sma_graph_df[, 4]))
  sma_graph_df[, 5] <- as.numeric(as.character(sma_graph_df[, 5]))
  sma_graph_df[, 6] <- as.numeric(as.character(sma_graph_df[, 6]))
  sma_graph_df[, 7] <- as.numeric(as.character(sma_graph_df[, 7]))
  return(sma_graph_df)
}

# Makes SMA plots for Families all on one graph
# ==============================================================================
mk_ghFG_SMAplot <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = gh, colour=j_fg)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab("log(vertical gape, mm)") +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) 
}

mk_gwFG_SMAplot <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = gw, colour=j_fg)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab("log(horizontal gape, mm)") +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) 
}

mk_gaFG_SMAplot <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = ga, colour=j_fg)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) 
}

# Makes SMA plots for FGs facetted
# ==============================================================================

mk_ghFG_SMAfacet <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = gh)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab("log(vertical gape, mm)") +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) +
    facet_wrap(~ j_fg)
}

mk_gwFG_SMAfacet <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = gw)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab("log(horizontal gape, mm)") +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) +
    facet_wrap(~ j_fg)
}

mk_gaFG_SMAfacet <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = ga)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) +
    facet_wrap(~ j_fg)
}

# Makes SMA plots for Families all on one graph
# ==============================================================================
mk_ghFam_SMAplot <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = gh, colour=Family)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab("log(vertical gape, mm)") +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) 
}

mk_gwFam_SMAplot <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = gw, colour=Family)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab("log(horizontal gape, mm)") +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) 
}

mk_gaFam_SMAplot <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = ga, colour=Family)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) 
}

# Makes SMA plots for families all on separate facets (colour = black)
# ==============================================================================
mk_ghFam_SMAfacets <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = gh)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab("log(vertical gape, mm)") +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, 
                 aes(x = from, xend = to, y = yfrom, yend = yto)
    ) +
    facet_wrap(~ Family)
}
mk_gwFam_SMAfacets <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = gw)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab("log(horizontal gape, mm)") +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) +
    facet_wrap(~ Family)
}

mk_gaFam_SMAfacets <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = ga)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
    scale_colour_discrete(name = "Functional \n Group") +
    theme(legend.key.height = unit(1.5, "line")) +
    theme(legend.position = c(0.90, 0.35)) +
    theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) +
    facet_wrap(~ Family)
}


# Makes SMA plots for species all on one graph
# ==============================================================================
mk_gh_SMAplot <- function(df_points, df_lines, colour_by) {
  ggplot(data = df_points, aes(x = SL, y = gh, colour=SpeciesCode)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab("log(vertical gape, mm)") +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) 
}

mk_gw_SMAplot <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = gw, colour=SpeciesCode)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab("log(horizontal gape, mm)") +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) 
}

mk_ga_SMAplot <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = ga, colour=SpeciesCode)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) 
}

# Makes SMA plots for species all on separate facets (colour = black)
# ==============================================================================
mk_gh_SMAplot_facets <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = gh)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab("log(vertical gape, mm)") +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, 
                 aes(x = from, xend = to, y = yfrom, yend = yto)
                 ) +
    facet_wrap(~ SpeciesCode)
}
 mk_gw_SMAplot_facets <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = gw)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab("log(horizontal gape, mm)") +
    scale_colour_discrete(name = "Functional \n Group") +
    #theme(legend.key.height = unit(1.5, "line")) +
    #theme(legend.position = c(0.90, 0.35)) +
    #theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) +
    facet_wrap(~ SpeciesCode)
}

mk_ga_SMAplot_facets <- function(df_points, df_lines) {
  ggplot(data = df_points, aes(x = SL, y = ga)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    #scale_x_log10(limits=c(1, 1000)) +
    xlab("log(standard length, mm)") +
    ylab(expression(paste("log(gape area ", mm^2, ")", sep= ""))) +
    scale_colour_discrete(name = "Functional \n Group") +
    theme(legend.key.height = unit(1.5, "line")) +
    theme(legend.position = c(0.90, 0.35)) +
    theme(legend.background = element_rect(fill = "#FFFFFFaa", colour = 'NA')) +
    geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, yend = yto)) +
    facet_wrap(~ SpeciesCode)
}

# groupwise_lm_gh <- function(df, variable) {
#   lm  <- with(data=df, ddply(df, .(variable), function(z) {
#     t <- lm(log(gh)~log(SL), data=z)
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
#     t <- lm(log(gw)~log(SL), data=z)
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
#     t <- lm(log(ga)~log(SL), data=z)
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
      eqn <- substitute(italic(y) ==
                          slp%.%italic(x) + int*~~italic(r)^2~"="~r2, l)
    } else {
        l <- list(slp = format(df[[3]][i], digits=2),
                  int = format(abs(df[[2]][i]), digits=2), 
                  r2 = format(df[[4]][i], digits=2)
      )
      eqn <- substitute(italic(y) ==
                          slp%.% italic(x) - int*~~italic(r)^2~"="~r2, l)
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
#     t <- lm(log(gw)~log(SL), data=z)
#     data.frame(int   = coefficients(t)[1],
#                slope = coefficients(t)[2],
#                rsq   = summary(t)$r.squared,
#                se    = summary(t)$coefficients[2,2],
#                p_val = summary(t)$coef[2,4])
#   }))
# }























