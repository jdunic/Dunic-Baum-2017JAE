# Change active plot window (quartz())
A <- function(PlotWindow) {
  dev.set(which=PlotWindow)
}

# ggplot default colour list for n colours
gg_colour_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

################################################################################
#############             Graph equations and formatting            ############
################################################################################
# Get geometric mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# extract ggplot legend to stick beside gridExtra plots as desired:
extract_legend <- function(a.gplot){
  tmp    <- ggplot_gtable(ggplot_build(a.gplot))
  leg    <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

write_sma_eqn <- function(df, y){
  m = lm(log(y) ~ log(SL), df);
  l <- list(a = format(coef(m)[1], digits = 2), 
            b = format(coef(m)[2], digits = 2), 
            r2 = format(summary(m)$r.squared, digits = 2, nsmall = 2)
  )  
  if (l$a >= 0) {
    eq <- substitute(italic(y) == b %.% italic(x) + a*","~~italic(r)^2~"="~r2, l) 
  } else {
      l <- list(a = format(abs(coef(m)[1]), digits = 2), 
                b = format(coef(m)[2], digits = 2), 
                r2 = format(summary(m)$r.squared, digits = 2, nsmall = 2)
    )
    eq <- substitute(italic(y) == b %.% italic(x) - a*","~~italic(r)^2~"="~r2, l)
  }
  as.character(as.expression(eq))
}

write_group_sma_eqn <- function(sma_summary_df, group_column) {
  df <- sma_summary_df
  m = matrix(data=NA, nrow=0, ncol=5)
  count <- length(group_column)
  for (i in (1:count)) {
    l <- list(slp = format(sma_summary_df$slope[i], digits=2),
              int = format(sma_summary_df$elev[i], digits=2), 
              r2  = format(sma_summary_df$xy_r2[i], digits=2, nsmall = 2),
              count = sma_summary_df$n[i], 
              sig = sma_summary_df$sig[i]
    )
    if (l$int >= 0) {
      eqn_r2 <- substitute(atop(italic(r)^2~"="~r2, italic(y) ==
                          slp%.%italic(x) + int), l)
      eqn <- substitute(italic(y) == slp*italic(x) + int*sig, l)
      r2  <- substitute(italic(r)^2~"="~r2, l)
      n   <- substitute(italic(n) ~ "=" ~ count, l)
    } else {
        l <- list(slp = format(sma_summary_df$slope[i], digits=2),
                  int = format(abs(sma_summary_df$elev[i]), digits=2),
                  r2  = format(sma_summary_df$xy_r2[i], digits=2, nsmall = 2),
                  count   = sma_summary_df$n[i], 
                  sig = sma_summary_df$sig[i]
      )
      eqn_r2 <- substitute(atop(italic(r)^2~"="~r2, italic(y) ==
                          slp%.% italic(x) - int), l)
      eqn <- substitute(italic(y) == slp*italic(x) - int*sig, l)
      r2  <- substitute(italic(r)^2~"="~r2, l)
      n   <- substitute(italic(n) ~ "=" ~ count, l)
    }    
    #browser()
    eqn_r2 <- as.character(as.expression(eqn_r2)) 
    eqn    <- as.character(as.expression(eqn))
    r2     <- as.character(as.expression(r2))
    n      <- as.character(as.expression(n))
    m <- rbind(m, c(as.character(df[[1]][i]), eqn_r2, eqn, r2, n))
    #m <- rbind(m, c(as.character(df[i,1]), lm_eq))
  }
  m <- as.data.frame(m)
}

count_spp <- function(df) {
  ddply(.data = df, .(SpeciesCode), summarize, 
        len = length(SpeciesCode),
        n = paste("n ==", len)
  )
}

# test whether a slope is allometric 
get_allometry <- function(slope, p_val, iso_val = 1) {
    allometry <- 'I'
    if (p_val <= 0.05 & slope < iso_val) allometry <- 'N'
    if (p_val <= 0.05 & slope > iso_val) allometry <- 'P'
    return(allometry)
}

get_sig <- function(p_val) {
    sig <- ''
    if (p_val < 0.05) sig <- "*"
    return(sig)
}

################################################################################
#############                    SMA Functions                      ############
################################################################################
run_sma <- function(df, gapeType=c("gh", "gw", "ga"), robust=TRUE) {
  if (robust == TRUE) {
    switch(gapeType,
    "gh" = { sma(gh ~ SL, data = df, log = "xy", method = "SMA", robust = TRUE, 
      slope.test = 1) },
    "gw" = { sma(gw ~ SL, data = df, log = "xy", method = "SMA", robust = TRUE, 
      slope.test = 1) },
    "ga" = { sma(ga ~ SL, data = df, log = "xy", method = "SMA", robust = TRUE, 
      slope.test = 2) }
    )
  } else if (robust == FALSE) {
    switch(gapeType, 
    "gh" = { sma(gh ~ SL, data = df, log = "xy", method = "SMA", robust = FALSE, 
      slope.test = 1) },
    "gw" = { sma(gw ~ SL, data = df, log = "xy", method = "SMA", robust = FALSE, 
      slope.test = 1) },
    "ga" = { sma(ga ~ SL, data = df, log = "xy", method = "SMA", robust = FALSE, 
      slope.test = 2) }
    )
  }
}

check_assump <- function(sma_object, plotTitle) {
  plot(sma_object, which = "qq")
  plot(sma_object, which = "residual")
  title(main = plotTitle)
  abline(h=0, col="red")
}

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


mk_sma_graph_df <- function(sma_summary_df, num_groups, group_name) {
  sma_graph_df <- data.frame(group=character(), slp=numeric(), int=numeric(), 
                             from=numeric(), to=numeric(), yfrom=numeric(), 
                             yto=numeric(), stringsAsFactors=FALSE
  )
  for (i in 1:num_groups) {
    from  <- sma_summary_df[10, i]
    to    <- sma_summary_df[11, i]
    slp   <- sma_summary_df[3, i]
    int   <- sma_summary_df[1, i]
    yfrom <- 10^(slp*log10(from) + int)
    yto   <- 10^(slp*log10(to) + int)
    group <- colnames(sma_summary_df)[i]
    midpoint_y <- sqrt(yfrom * yto)
    midpoint_x <- sqrt(from * to)
    ref_intercept <- log10(midpoint_y/(midpoint_x^2))
    
    row <- t(c(group=group, slp=slp, int=int, from=from, to=to, yfrom=yfrom,
               yto=yto, midpoint_x=midpoint_x, midpoint_y=midpoint_y, 
               ref_intercept=ref_intercept))
    sma_graph_df <- rbind(sma_graph_df, row)
  }
  sma_graph_df[, 2] <- as.numeric(as.character(sma_graph_df[, 2]))
  sma_graph_df[, 3] <- as.numeric(as.character(sma_graph_df[, 3]))
  sma_graph_df[, 4] <- as.numeric(as.character(sma_graph_df[, 4]))
  sma_graph_df[, 5] <- as.numeric(as.character(sma_graph_df[, 5]))
  sma_graph_df[, 6] <- as.numeric(as.character(sma_graph_df[, 6]))
  sma_graph_df[, 7] <- as.numeric(as.character(sma_graph_df[, 7]))
  sma_graph_df[, 8] <- as.numeric(as.character(sma_graph_df[, 8]))
  sma_graph_df[, 9] <- as.numeric(as.character(sma_graph_df[, 9]))
  sma_graph_df[, 10] <- as.numeric(as.character(sma_graph_df[, 10]))
  names(sma_graph_df)[1] <- group_name
  return(sma_graph_df)
}

# Used for single group
mk_sma_summary <- function(sma_object, group="column_name") {
  rows <- c('elev', 'slp_test', 'slope', 'lower_ci', 'upper_ci',
            'slp_p_value', 'xy_r^2', 'xy_corr_p_value', 'n', 'from', 'to')
  #if (grouping==F) {
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

  columns <- c(elev, slp_test, slope, lower, upper, slp_p_val, xy_r2, xy_cor, n,
               from, to)
  
  sma_df <- data.frame(columns, row.names=rows)
  names(sma_df) <- sma_object$groups
  #sma_df <- format(sma_df[1], digits=3, sci=F)
  return(sma_df)
}

mk_spp_summary <- function(sma_object, num_spp=NA, grouping=F, group_name) {
# Use (grouping == F) when multiple sma_objects are generated using dlply
# Use (grouping == T) when the sma_object was generated using x~y*group
  if (grouping==F) {
    sma_df <- data.frame(elev=numeric(), slp_test=numeric(), slope=numeric(), 
                      upper=numeric(), lower=numeric(), slp_p_val=numeric(), 
                      xy_r2=numeric(), xy_cor=numeric(), n=numeric(), 
                      from=numeric(), to=numeric()
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
      
      row <- c(elev, slp_test, slope, lower, upper, slp_p_val, xy_r2, xy_cor, 
               n, from, to)
      sma_df <- rbind(sma_df, row)
      columns <- c("elev", "slp_test", "slope", "lower", "upper", "slp_p_val", 
               "xy_r2", "xy_cor", "n", "from", "to")
  
      #sma_df <- data.frame(columns, row.names=rows)
      names(sma_df) <- columns
      #sma_df <- format(sma_df[1], digits=3, sci=F)
    }
    return(sma_df)
  } else if (grouping==T) {

    sma_df <- data.frame(group=character(), elev=numeric(), slp_test=numeric(), 
                      slope=numeric(), upper=numeric(), lower=numeric(), 
                      slp_p_val=numeric(), xy_r2=numeric(), xy_cor=numeric(), 
                      n=numeric(), from=numeric(), to=numeric(),
                      stringsAsFactors=FALSE
    )
    for (i in 1:length(sma_object$groups)) {

      elev = sma_object$coef[[i]][[1]][1]
      slp_test = sma_object$slopetest[[i]]$test.value
      slope = sma_object$slopetest[[i]]$b
      lower = sma_object$slopetest[[i]]$ci[1, 1]
      upper = sma_object$slopetest[[i]]$ci[1, 2]
      slp_p_val = sma_object$slopetest[[i]]$p
      xy_r2 = sma_object$r[i][[1]]
      xy_cor = sma_object$pval[i][[1]]
      n = sma_object$n[i][[1]]
      from = sma_object$from[i][[1]]
      to = sma_object$to[i][[1]]
      group = sma_object$groups[i][[1]]

      row <- c("group"=as.character(group), "elev"=as.numeric(elev), 
               "slp_test"=as.numeric(slp_test), "slope"=as.numeric(slope),
               "lower"=as.numeric(lower), "upper"=as.numeric(upper), 
               "slp_p_val"=as.numeric(slp_p_val), "xy_r2"=as.numeric(xy_r2), 
               "xy_cor"=as.numeric(xy_cor), "n"=as.numeric(n), 
               "from"=as.numeric(from), "to"=as.numeric(to))
      sma_df[i, ] <- row
    }
  columns <- c("group", "elev", "slp_test", "slope", "lower", "upper", 
         "slp_p_val", "xy_r2", "xy_cor", "n", "from", "to")
  #sma_df <- data.frame(columns, row.names=rows)
  names(sma_df) <- columns
  #sma_df <- format(sma_df[1], digits=3, sci=F)
  
  for (x in 2:12) {
    sma_df[, x] <- as.numeric(sma_df[, x])
  }
  return(sma_df)
  }
}



mk_smaSPP_graph_df <- function(sma_summary_df, num_spp, group_name, iso_slope = 1) {
  sma_graph_df <- data.frame(group=character(), slp=numeric(), int=numeric(), 
                             from=numeric(), to=numeric(), yfrom=numeric(), 
                             yto=numeric(),
                             stringsAsFactors=FALSE
  )
  for (i in 1:num_spp) {
    from  <- sma_summary_df[i, 11]
    to    <- sma_summary_df[i, 12]
    slp   <- sma_summary_df[i, 4]
    int   <- sma_summary_df[i, 2]
    yfrom <- 10^(slp*log10(from) + int)
    yto   <- 10^(slp*log10(to) + int)
    group <- as.character(sma_summary_df[i, 1])
    midpoint_y <- sqrt(yfrom * yto)
    midpoint_x <- sqrt(from * to)
    ref_intercept_iso <- log10(midpoint_y / midpoint_x ^ (iso_slope))
    slope_test <- sma_summary_df[i, 3]
#    
    row <- t(c(group=group, slp=slp, int=int, from=from, to=to, yfrom=yfrom,
               yto=yto, midpoint_x=midpoint_x, midpoint_y=midpoint_y, 
               ref_intercept_iso=ref_intercept_iso, slope_test = slope_test)
             )
    sma_graph_df <- rbind(sma_graph_df, row)
  }
  sma_graph_df[, 2]  <- as.numeric(as.character(sma_graph_df[, 2]))
  sma_graph_df[, 3]  <- as.numeric(as.character(sma_graph_df[, 3]))
  sma_graph_df[, 4]  <- as.numeric(as.character(sma_graph_df[, 4]))
  sma_graph_df[, 5]  <- as.numeric(as.character(sma_graph_df[, 5]))
  sma_graph_df[, 6]  <- as.numeric(as.character(sma_graph_df[, 6]))
  sma_graph_df[, 7]  <- as.numeric(as.character(sma_graph_df[, 7]))
  sma_graph_df[, 8]  <- as.numeric(as.character(sma_graph_df[, 8]))
  sma_graph_df[, 9]  <- as.numeric(as.character(sma_graph_df[, 9]))
  sma_graph_df[, 10] <- as.numeric(as.character(sma_graph_df[, 10]))
  sma_graph_df[, 11] <- as.numeric(as.character(sma_graph_df[, 11]))
  names(sma_graph_df)[1] <- group_name
  return(sma_graph_df)
}

# Makes SMA plots for Families all on one graph
# ==============================================================================
mk_SMAfacets <- function( df_points, df_lines, gapeType = c("gh", "gw", "ga"), 
  point_colour = c("j_fg", "Family", "SpeciesCode", "Region", "dissected_by", 
                   "observer_id"),
  labels = c("dissected_by", "Region", "SpecimenID", "None"),
  facetting = c("j_fg", "Family", "SpeciesCode", "Region", "dissected_by", 
                "observer_id"), 
  facet_columns ) { 
  plot_base <- ggplot(data = df_points, aes_string(x = "SL", y = gapeType)) +
       geom_point( aes_string(colour = point_colour)) +
       geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, 
        yend = yto)) +
       scale_y_log10() +
       scale_x_log10() +
       xlab("log(standard length, mm)") +
       theme_bw()
  switch(gapeType,
    "gh" = { plot_base <- plot_base + ylab("log(vertical gape, mm)") },
      "gw" = { plot_base <- plot_base + ylab("log(horizontal gape, mm)") },
      "ga" = { plot_base <- plot_base + ylab(expression(
        paste("log(gape area ", mm^2, ")", sep= ""))) }
  )
  if (labels == "None") {
    plot1 <- plot_base
  } else {
    plot1 <- plot_base + geom_text(position = position_jitter(w = 0.02, 
      h = 0.02), aes_string(label = labels), size = 2)
  }
  plot1 + facet_wrap( as.formula(sprintf('~ %s', facetting)), ncol = facet_columns )
}

theme_L_border <- function(colour = "black", size = 1, linetype = 1) {
  structure(
    function(x = 0, y = 0, width = 1, height = 1, ...) {
      polylineGrob(
        x=c(x+width, x, x), y=c(y,y,y+height), ..., default.units = "npc",
        gp=gpar(lwd=size, col=colour, lty=linetype),
      )
    },
    class = "theme",
    type = "box",
    call = match.call()
  )
}

mk_SMAfacets2 <- function( df_points, df_lines, gapeType = c("gh", "gw", "ga"), 
    #point_colour = c("j_fg", "Family", "SpeciesCode", "Region", "dissected_by"),
    labels = c("dissected_by", "Region", "SpecimenID", "None"),
    facetting = c("j_fg", "Family", "SpeciesCode", "Region", "dissected_by"), 
    facet_columns, eqn_df
    ) {
  plot_base <- ggplot(data = df_points, aes_string(x = "SL", y = gapeType)) +
       geom_point(shape = 1, colour = "grey") +
       geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, 
        yend = yto)) +
       scale_y_log10() +
       scale_x_log10() +
       xlab("log(standard length, mm)") +
       theme_classic() +
       theme(strip.background = element_blank(),
         plot.background = element_blank(), 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), 
         panel.border = element_blank(), 
         panel.background = element_blank(),
         axis.line = element_line(colour = "black"),
         axis.line.x = element_line(colour = "black")
        ) +
        geom_point(aes(x = 10, y = 1), alpha = 0) +
        geom_point(aes(x = 650, y = 12000), alpha = 0) +
        geom_text(data = eqn_df, aes(x=280, y=3.5, 
          label=eqn), parse=TRUE, size = 3.5) +
        geom_abline(data = df_lines, aes_string(intercept = "ref_intercept"), 
          slope = 2, linetype = 2, colour = "grey50") 
  switch(gapeType,
    "gh" = { plot_base <- plot_base + ylab("log(vertical gape, mm)") },
      "gw" = { plot_base <- plot_base + ylab("log(horizontal gape, mm)") },
      "ga" = { plot_base <- plot_base + ylab(expression(
        paste("log(gape area ", mm^2, ")", sep= ""))) }
  )
  if (labels == "None") {
    plot1 <- plot_base
  } else {
    plot1 <- plot_base + geom_text(position = position_jitter(w = 0.02, 
      h = 0.02), aes_string(label = labels), size = 2)
  }
  plot1 + facet_wrap( as.formula(sprintf('~ %s', facetting)), ncol = facet_columns, 
  scales = "free")
}

mk_multipanel_plots2 <- function(fg_point_df, spp_point_df, spp_line_df_row, 
  #ref_intercept_row, 
  eqn_df, eqn_x, eqn_y, r2_x, r2_y, n_x, n_y, x_axis_labels=TRUE, 
  y_axis_labels=TRUE, fg_line_intercept, y_axis_text = TRUE, x_axis_text = TRUE,
  plot_title = "", y_value, gape_dim = 'gh') 
  {
  plotTitle <- substitute(italic(plot_title), list(plot_title = plot_title))
  plot_base <- 
      ggplot(data = fg_point_df, aes_string(x = "SL", y = gape_dim)) +
        geom_point(shape = 1, colour = "grey") +
        geom_segment(data = spp_line_df_row, aes_string(x = "from", xend = "to", 
         y = "yfrom", yend = "yto")) +
        geom_point(data = spp_point_df, colour = "black", shape = 1) +
        scale_y_log10() +
        scale_x_log10() +
      geom_abline(intercept = fg_line_intercept, slope = 1, linetype = 2, 
        colour = "darkgrey") +
      theme_bw() +
      theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
      theme(axis.line.x = element_line(color = 'black'), 
            axis.line.y = element_line(color = 'black')) +
      geom_text(data = eqn_df, aes_string(x = eqn_x, y = eqn_y, 
        label = "eqn"), parse = TRUE, size = 3, hjust = 1) +
      geom_text(data = eqn_df, aes_string(x = r2_x, y = r2_y, 
        label = "r2"), parse = TRUE, size = 3, hjust = 1) +
      geom_text(data = eqn_df, aes_string(x = n_x, y = n_y, 
        label = "n"), parse = TRUE, size = 3, hjust = 1) +
      labs(title = bquote(plain(.(plotTitle)))) +
      theme(plot.title = element_text(size = 9), 
            axis.text = element_text(size = 8),
            axis.ticks.length = unit(-0.1, "cm"),
            axis.text.y = element_text(margin = margin(0, 5, 0, 0)), 
            axis.text.x = element_text(margin = margin(5, 0, 0, 0), vjust = 1))
      #plot <- plot_base +
  if (x_axis_labels == TRUE) {
    plot1 <- plot_base + xlab("standard length, mm")
  } else if (x_axis_labels == FALSE) {
    plot1 <- plot_base + theme(axis.title.x = element_blank())
  }
  if (y_axis_labels == TRUE) {
    plot2 <- plot1 + ylab(expression(paste("gape height, ", mm, "", sep= "")))
  } else if (y_axis_labels == FALSE) {
    plot2 <- plot1 + theme(axis.title.y = element_blank())
  } 
  if (y_axis_text == TRUE) {
    plot3 <- plot2
  } else if (y_axis_text == FALSE) {
    plot3 <- plot2 + theme(axis.text.y = element_blank())
  }
  if (x_axis_text == TRUE) {
    plot4 <- plot3
  } else if (x_axis_text == FALSE) {
    plot4 <- plot3 + theme(axis.text.x = element_blank())
  }
  plot4  
}

mk_multipanel_plots_mass <- function(fg_point_df, spp_point_df, spp_line_df_row, 
  ref_intercept_row, 
  eqn_df, eqn_x, eqn_y, r2_x, r2_y, n_x, n_y, x_axis_labels=TRUE, 
  y_axis_labels=TRUE, fg_line_intercept, y_axis_text = TRUE, x_axis_text = TRUE,
  plot_title = "", y_value, gape_dim = 'gh') 
  {
  plotTitle <- substitute(italic(plot_title), list(plot_title = plot_title))
  plot_base <- 
      ggplot(data = fg_point_df, aes_string(x = "wt", y = gape_dim)) +
        geom_point(shape = 1, colour = "grey") +
        geom_segment(data = spp_line_df_row, aes_string(x = "from", xend = "to", 
         y = "yfrom", yend = "yto")) +
        geom_point(data = spp_point_df, colour = "black", shape = 1) +
        scale_y_log10() +
        scale_x_log10() +
      geom_abline(intercept = fg_line_intercept, slope = 1 / 3, linetype = 2, 
        colour = "darkgrey") +
      theme_bw() +
      theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
      theme(axis.line = element_line(color = 'black')) +
      geom_text(data = eqn_df, aes_string(x = eqn_x, y = eqn_y, 
        label = "eqn"), parse = TRUE, size = 3, hjust = 1) +
      geom_text(data = eqn_df, aes_string(x = r2_x, y = r2_y, 
        label = "r2"), parse = TRUE, size = 3, hjust = 1) +
      geom_text(data = eqn_df, aes_string(x = n_x, y = n_y, 
        label = "n"), parse = TRUE, size = 3, hjust = 1) +
      labs(title = bquote(plain(.(plotTitle)))) +
      theme(plot.title = element_text(size = 9), 
            axis.text = element_text(size = 8),
            axis.ticks.length = unit(-0.1, "cm"),
            axis.text.y = element_text(margin = margin(0, 5, 0, 0)), 
            axis.text.x = element_text(margin = margin(5, 0, 0, 0), vjust = 1))
      #plot <- plot_base +
  if (x_axis_labels == TRUE) {
    plot1 <- plot_base + xlab("standard length, mm")
  } else if (x_axis_labels == FALSE) {
    plot1 <- plot_base + theme(axis.title.x = element_blank())
  }
  if (y_axis_labels == TRUE) {
    plot2 <- plot1 + ylab(expression(paste("gape height, ", mm, "", sep= "")))
  } else if (y_axis_labels == FALSE) {
    plot2 <- plot1 + theme(axis.title.y = element_blank())
  } 
  if (y_axis_text == TRUE) {
    plot3 <- plot2
  } else if (y_axis_text == FALSE) {
    plot3 <- plot2 + theme(axis.text.y = element_blank())
  }
  if (x_axis_text == TRUE) {
    plot4 <- plot3
  } else if (x_axis_text == FALSE) {
    plot4 <- plot3 + theme(axis.text.x = element_blank())
  }
  plot4  
}

mk_corallivore_plot <- function(fg_point_df, spp_point_df, 
  #ref_intercept_row, 
  eqn_df, eqn_x, eqn_y, r2_x, r2_y, n_x, n_y, x_axis_labels=TRUE, 
  y_axis_labels=TRUE, fg_line_intercept, y_axis_text = TRUE, x_axis_text = TRUE,
  plot_title = "", y_value, gape_dim = 'gh') 
  {
  plotTitle <- substitute(italic(plot_title), list(plot_title = plot_title))
  plot_base <- 
      ggplot(data = fg_point_df, aes_string(x = "SL", y = gape_dim)) +
        geom_point(shape = 1, colour = "grey") +
        geom_point(data = spp_point_df, colour = "black", shape = 1) +
        scale_y_log10() +
        scale_x_log10() +
      geom_abline(intercept = fg_line_intercept, slope = 1, linetype = 2, 
        colour = "darkgrey") +
      theme_bw() +
      theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
      theme(axis.line = element_line(color = 'black')) +
      geom_text(data = eqn_df, aes_string(x = eqn_x, y = eqn_y, 
        label = "eqn"), parse = TRUE, size = 3, hjust = 1) +
      geom_text(data = eqn_df, aes_string(x = r2_x, y = r2_y, 
        label = "r2"), parse = TRUE, size = 3, hjust = 1) +
      geom_text(data = eqn_df, aes_string(x = n_x, y = n_y, 
        label = "n"), parse = TRUE, size = 3, hjust = 1) +
      labs(title = bquote(plain(.(plotTitle)))) +
      theme(plot.title = element_text(size = 9), 
            axis.text = element_text(size = 8),
            axis.ticks.length = unit(-0.1, "cm"),
            axis.text.y = element_text(margin = margin(0, 5, 0, 0)), 
            axis.text.x = element_text(margin = margin(5, 0, 0, 0), vjust = 1))
      #plot <- plot_base +
  if (x_axis_labels == TRUE) {
    plot1 <- plot_base + xlab("standard length, mm")
  } else if (x_axis_labels == FALSE) {
    plot1 <- plot_base + theme(axis.title.x = element_blank())
  }
  if (y_axis_labels == TRUE) {
    plot2 <- plot1 + ylab(expression(paste("gape height, ", mm, "", sep= "")))
  } else if (y_axis_labels == FALSE) {
    plot2 <- plot1 + theme(axis.title.y = element_blank())
  } 
  if (y_axis_text == TRUE) {
    plot3 <- plot2
  } else if (y_axis_text == FALSE) {
    plot3 <- plot2 + theme(axis.text.y = element_blank())
  }
  if (x_axis_text == TRUE) {
    plot4 <- plot3
  } else if (x_axis_text == FALSE) {
    plot4 <- plot3 + theme(axis.text.x = element_blank())
  }
  plot4  
}


mk_SMAplot <- function(df_points, df_lines, facets = TRUE, x = "SL", gapeType = 
  c("gh", "gw", "ga"), grouping = c("j_fg", "Family", "SpeciesCode", "Region", 
  "dissected_by"), labels = c("Region", "Region_colour", "dissected_by",  
  "dissected_colour", "SpecimenID", "None"), axis_labels) {
  plot_base <- ggplot(data = df_points, aes_string(x = x, y = gapeType)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("standard length (mm)") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  #labs(title = bquote(plain(.(plotTitle)))) +
      #labs(title = bquote(italic(.(plotTitle)))) +
  theme(axis.ticks.length = unit(-0.2, "cm")) +
  theme(axis.text = element_text(size = 8),
        axis.ticks.length = unit(-0.1, "cm"),
        axis.text.y = element_text(margin = margin(0, 5, 0, 0)), 
        axis.text.x = element_text(margin = margin(5, 0, 0, 0), vjust = 1)) + 
  scale_x_log10(breaks=scales::pretty_breaks(3)) +
  scale_y_log10(breaks=scales::pretty_breaks(3))
  if (facets == FALSE) {
    plot1 <- plot_base + geom_point( aes_string(colour = grouping), size = 1.5 ) +
             geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, 
                yend = yto)) + aes_string(colour = grouping)
    switch(labels,
    "dissected_by" = { plot1 <- plot1 + 
      geom_text(position=position_jitter(w=0.03, h=0.03), 
        aes(label=dissected_by), size=3) },
    "Region" = { plot1 <- plot1 + 
      geom_text(position=position_jitter(w=0.03, h=0.03), 
        aes(label=dissected_by), size=3) },
    "None" = { plot1 <- plot1 }
      )
    switch(gapeType,
      "gh" = { plot1 + ylab("gape height (mm)") },
      "gw" = { plot1 + ylab("gape width, (mm)") },
      "ga" = { plot1 + ylab(expression(paste("gape area (", mm^2, ")",  
        sep= ""))) }
      )
  } else if (facets == TRUE) {
    plot2 <- plot_base + geom_point() +
             geom_segment(data = df_lines, aes(x = from, xend = to, y = yfrom, 
                yend = yto))
    switch(labels,
    "dissected_by" = { plot2 <- plot2 + 
      geom_point( aes_string(colour=labels) ) +
      geom_text(position=position_jitter(w=0.03, h=0.03), 
        aes(label=dissected_by), size=3) },
    "dissected_colour" = { plot2 <- plot2 +
      geom_point( aes(colour=dissected_by)) },
    "Region" = { plot2 <- plot2 + 
      geom_point( aes_string(colour=labels) ) +
      geom_text(position=position_jitter(w=0.03, h=0.03), 
        aes(label=Region), size=3) },
    "Region_colour" = { plot2 <- plot2 +
      geom_point( aes(colour=Region)) },
    "None" = { plot2 <- plot2 }
      )
    switch(gapeType,
      "gh" = { plot3 <- plot2 + ylab("gape height (mm)") },
      "gw" = { plot3 <- plot2 + ylab("gape width (mm)") },
      "ga" = { plot3 <- plot2 + ylab(expression(paste("gape area (", mm^2, ")", 
        sep= ""))) }
      )
    switch(grouping,
      "j_fg" = { plot3 + facet_wrap( ~ j_fg) },
      "Family" = { plot3 + facet_wrap( ~ Family) },
      "SpeciesCode" = { plot3 + facet_wrap( ~ SpeciesCode) },
      "Region" = { plot3 + facet_wrap( ~ Region) },
      "dissected_by" = { plot3 + facet_wrap( ~ dissected_by) }
      )
  }
}

# Plotting using grid, with a set viewport
set_vp <- function(row, column) {
  viewport(layout.pos.row = row, layout.pos.col = column)
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

# Mapping function taken from SO answer by Joris Meys from:
# http://stackoverflow.com/questions/5353184/fixing-maps-library-data-for-pacific-centred-0-360-longitude-display
# Used to adjust polygons so that they are not left 'open' on the cut when the 
# ends (for "world" and "worldHiRes") when the map is pacific ocean-centric.
# xlimits have been added to the final call for map() at the end of the function 
# because they were causing islands in the pacific to disappear in first part of 
# the function where the polygons are moved around. 
plot.map <- function(database,center, xlimits, ...){
    Obj <- map(database,...,plot=F)
    coord <- cbind(Obj[[1]],Obj[[2]])

    # split up the coordinates
    id <- rle(!is.na(coord[,1]))
    id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
    polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})

    # split up polygons that differ too much
    polygons <- lapply(polygons,function(x){
        x[,1] <- x[,1] + center
        x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
        if(sum(diff(x[,1])>300,na.rm=T) >0){
          id <- x[,1] < 0
          x <- rbind(x[id,],c(NA,NA),x[!id,])
       }
       x
    })
    # reconstruct the object
    polygons <- do.call(rbind,polygons)
    Obj[[1]] <- polygons[,1]
    Obj[[2]] <- polygons[,2]

    map(Obj,..., xlim=xlimits)
}

