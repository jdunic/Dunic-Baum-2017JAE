library(plyr)
library(ggplot2)
library(doMC)

load('01_load.r')
load('02_clean.r')
load('03_func.r')

# bootstrapping sma()
ddply(data = p, .(SpeciesCode), fun(x) sample(x, size = 12, replace = TRUE, prob = NULL))

x <- ddply(.data = pento, .(j_fg, SpeciesCode), summarise, len = length(SpeciesCode))
ddply(x, .(j_fg), summarise, min(len))

#-------------------------------------------------------------------------------
# Bootstrapping functions
# Hadley Wickham answer on SO!!! 
# How to sub sample data by group efficiently using plyr
# http://stackoverflow.com/questions/16912186/how-do-i-sub-sample-data-by-group-efficiently
sample_by_spp <- function(df, dfgroup, size, replace) {
  ddply(df, dfgroup, function(x) {
    x[sample(nrow(x), size = size, replace = replace), , drop = FALSE]
  })
}

mk_all_spp_samp_df <- function(p_df, b_df, z_df, h_df, c_df){
  p_samp <- sample_by_spp(df = p_df, dfgroup = 'SpeciesCode', size = 12, replace = TRUE)
  b_samp <- sample_by_spp(df = b_df, dfgroup = 'SpeciesCode', size = 24, replace = TRUE)
  z_samp <- sample_by_spp(df = z_df, dfgroup = 'SpeciesCode', size = 8, replace = TRUE)
  h_samp <- sample_by_spp(df = h_df, dfgroup = 'SpeciesCode', size = 7, replace = TRUE)
  all_spp_samp <- rbind.fill(p_samp, b_samp, z_samp, h_samp, c_df)
  return(all_spp_samp)
}

slp_check <- function(x) {
  check_df <- data.frame(j_fg = c("Pi", "BI", "ZP", "He", "C"), iso = 0, 
               neg = 0, pos = 0)
  for (i in 1:5) {
    if (x$groupsummary$Slope_test_p[i] > 0.05) {
      check_df$iso[i] <- 1
    } else {
        if (x$groupsummary$Slope_lowCI[i] > 2) {
          check_df$pos[i] <- 1
        } else {
            if (x$groupsummary$Slope_highCI[i] < 2) {
              check_df$neg[i] <- 1
          }
        }
      }
   }
   return(check_df)
}

# For each bootSMA and the multiple pairwise comparison between functional groups
# count if there was a significant difference
multcomp_check <- function(x) {
  multcomp_check_df <- data.frame(j_fg_1 = x$multcompresult[, 1], 
                  j_fg_2 = x$multcompresult[, 2],
                  sig_diff = 0
             )
  for (i in 1:length(x$multcompresult[, 1])) {
    if (x$multcompresult$Pval[i] <= 0.05) {
      multcomp_check_df$sig_diff[i] <- 1
    }
  }
  return(multcomp_check_df)
}

# Make a summary df for each bootSMA object
mk_boot_summ_df <- function(x) {
  #browser()
  #print(x)
  p_from <- x$from$Pi
  b_from <- x$from$BI
  z_from <- x$from$ZP
  h_from <- x$from$He
  c_from <- x$from$C
#
  p_to <- x$to$Pi
  b_to <- x$to$BI
  z_to <- x$to$ZP
  h_to <- x$to$He
  c_to <- x$to$C
#
  p_slp <- x$groupsummary$Slope[1]
  b_slp <- x$groupsummary$Slope[2]
  z_slp <- x$groupsummary$Slope[3]
  h_slp <- x$groupsummary$Slope[4]
  c_slp <- x$groupsummary$Slope[5]
#
  p_int <- x$groupsummary$Int[1]
  b_int <- x$groupsummary$Int[2]
  z_int <- x$groupsummary$Int[3]
  h_int <- x$groupsummary$Int[4]
  c_int <- x$groupsummary$Int[5]
#
  p_r2  <- x$groupsummary$r2[1]
  b_r2  <- x$groupsummary$r2[2]
  z_r2  <- x$groupsummary$r2[3]
  h_r2  <- x$groupsummary$r2[4]
  c_r2  <- x$groupsummary$r2[5]
#
  group <- c("Pi", "BI", "ZP", "He", "C")
#
  summ_df <- data.frame(group = group, 
              from  = c(p_from, b_from, z_from, h_from, c_from),
              to    = c(p_to, b_to, z_to, h_to, c_to), 
              slope = c(p_slp, b_slp, z_slp, h_slp, c_slp), 
              int   = c(p_int, b_int, z_int, h_int, c_int),
              r2    = c(p_r2, b_r2, z_r2, h_r2, c_r2)
         )
  return(summ_df)
  summ_df
}


get_mean_slope_int_r2 <- function(x) {
  mean_slope <- mean(x$slope)
  mean_int <- mean(x$int)
  mean_r2 <- mean(x$r2)
  df <- data.frame(slope = mean_slope, int = mean_int, r2 = mean_r2)
  return(df)
}

get_y_FromTo <- function(x) {
  yfrom = 10^(x$slope * log10(x$from) + x$int)
  yto   = 10^(x$slope * log10(x$to) + x$int)
  df = data.frame(yfrom = yfrom, yto = yto)
  return(df)
}

write_bootSMA_eqn <- function(boot_df) {
    #browser()
    l <- list(slope = format(boot_df$slope, digits=2),
              int = format(boot_df$int, digits=2), 
              r2  = format(boot_df$r2, digits=2),
              count   = boot_df$n
    )
    print(l)
    if (l$int >= 0) {
      eqn_r2 <- substitute(atop(italic(r)^2~"="~r2, italic(y) ==
                          slope%.%italic(x) + int), l)
      eqn <- substitute(italic(y) == slope*italic(x) + int, l)
      r2  <- substitute(italic(r)^2~"="~r2, l)
      n   <- substitute(italic(n) ~ "=" ~ count, l)
    } else {
        l <- list(slope = format(boot_df$slope, digits=2),
                  int = format(abs(boot_df$int), digits=2),
                  r2  = format(boot_df$r2, digits=2),
                  count   = boot_df$n
      )
      eqn_r2 <- substitute(atop(italic(r)^2~"="~r2, italic(y) ==
                          slope%.% italic(x) - int), l)
      eqn <- substitute(italic(y) == slope*italic(x) - int, l)
      r2  <- substitute(italic(r)^2~"="~r2, l)
      n   <- substitute(italic(n) ~ "=" ~ count, l)
    }    
    #browser()
    eqn_r2 <- as.character(as.expression(eqn_r2)) 
    eqn    <- as.character(as.expression(eqn))
    r2     <- as.character(as.expression(r2))
    n      <- as.character(as.expression(n))
    df <- data.frame(eqn_r2 = eqn_r2, eqn = eqn, r2 = r2, n = n)
    return(df)
  }




#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# The bootstrapping
# The actual bootstrapping code has been commented out so that if this file is
# called with source, it will not spend 10 hours rerunning the code, it will 
# just load the SMA bootstrapped data from bootSMA_10000.RData located in the
# ~/R_projects/Allometry directory.

#m <- 10000  ## bootstrap replicates to use for offical MS values
#bootX <- rlply(m, mk_all_spp_samp_df(p, b, z, h, c), .progress = "text")

#registerDoMC()
#start = proc.time()
#bootSMA <- llply(bootX, function(x) {
#	sma(ga ~ SL * j_fg, data = x, log = "xy", method = "SMA", robust = T, 
#		slope.test = 2, multcomp = T, multcompmethod = "adjusted")
#	}, .parallel = TRUE)
#end = proc.time()
#total = (end[3]-start[3])
#cat(paste("Time Elapsed: ", total, "\n"))

#save(bootX, bootSMA, file="Allometry/bootSMA_10000.RData")
load("bootSMA_10000.RData")
# Note to even more future self, this object is called bootSMA.

#-------------------------------------------------------------------------------
# Allometry Check
# for each bootSMA and each functiona group in that bootSMA determine whether the 
# sma regression was negatively allometric, isometric, or positively allometric
slp_check_df <- data.frame(j_fg = c("Pi", "BI", "ZP", "He", "C"), iso = 0, 
						   neg = 0, pos = 0)

# Create dataframe to determine whether sma regression was neg, iso, or pos
# Checking allometry
check_df <- ldply(bootSMA, function(x) slp_check(x))
allo_check_df <- ddply(check_df, .(j_fg), summarise, 
             iso = sum(iso), neg = sum(neg), pos = sum(pos)
             )
allo_check_df$j_fg <- factor(allo_check_df$j_fg, levels = c("Pi", "BI", "ZP", "He", "C"))
allo_check_df <- arrange(allo_check_df, j_fg)

#-------------------------------------------------------------------------------
# Multiple Comparison Check
# For each bootSMA and the multiple pairwise comparison between functional groups
# count if there was a significant difference
multcomp_all_df <- ldply(bootSMA, function(x) multcomp_check(x))
multcomp_df <- ddply(multcomp_all_df, .(j_fg_1, j_fg_2), summarise, sig_diff = sum(sig_diff))
multcomp_df$prop <- multcomp_df$sig_diff / 10000
multcomp_df

#-------------------------------------------------------------------------------
# Average functional group slopes, intercepts, and r2
boot_summ  <- ldply(bootSMA, mk_boot_summ_df)
boot_means <- ddply(boot_summ, .(j_fg), get_mean_slope_int_r2)
boot_means
#  j_fg    slope        int         r2
#1   Pi 1.622224 -0.5110697 0.66575368
#2   BI 1.414475 -0.3602625 0.75836879
#3   ZP 1.935321 -1.8798773 0.83952185
#4   He 2.539374 -3.4064609 0.90062850
#5    C 4.003380 -6.3401771 0.07056937

# Use the overall functional group to find the midpoint of the cluster of points
allGA <- sma(ga~SL, data = pento, log = "xy", method = "SMA", robust = T, slope.test = 2)
#check_assump(allGA, "Pento Gape Area All")
allGA_summ <- mk_sma_summary(allGA, 1)
allGA_graph_df <- mk_sma_graph_df(allGA_summ, 1, "j_fg")
names(allGA_graph_df)

# SMA regression for each functional group
all_fg_GA <- sma(ga~SL * j_fg, data = pento, log = "xy", method = "SMA", robust = T, 
  slope.test = 2, multcomp = T, multcompmethod = "adjusted")
all_fg_GA_summ <- mk_spp_summary(all_fg_GA, 5, grouping=TRUE)
all_fg_GA_graph_df <- mk_smaSPP_graph_df(all_fg_GA_summ, 5, "j_fg")

# Get x and y to and from values to plot geom_segment for each functional group

# Get x to and from values
to_from_n_df <- ddply(pento, .(j_fg), summarise, from = min(SL), to = max(SL), 
  n = length(j_fg))
boot_graph_df <- merge(x = to_from_n_df, y = boot_means, by = "j_fg")
boot_graph_df <- arrange(boot_graph_df, j_fg)

# Get y to and from values
boot_y_FromTo <- ddply(boot_graph_df, .(j_fg), get_y_FromTo)
boot_graph_df <- merge(x = boot_y_FromTo, y = boot_graph_df, by = "j_fg")
boot_graph_df <- arrange(boot_graph_df, j_fg)
boot_graph_df

# Get bootSMA equations for the graph
bootSMA_eqns <- ddply(boot_graph_df, .(j_fg), function(x) write_bootSMA_eqn(x))
bootSMA_eqns

#save(boot_graph_df, file="~/R_projects/Allometry/boot_graph_df_10000.RData")

mk_bootSMA_plots <- function(fg_point_df, spp_point_df, spp_line_df_row, 
  x_axis_labels=TRUE, y_axis_labels=TRUE, fg_line_intercept, y_axis_text = TRUE, 
  x_axis_text = TRUE, plot_title = "") 
  {
  plotTitle <- substitute(italic(plot_title), list(plot_title = plot_title))
  plot_base <- 
      ggplot(data = fg_point_df, aes_string(x = "SL", y = "ga")) +
        geom_point(shape = 1, colour = "grey") +
        geom_segment(data = spp_line_df_row, aes_string(x = "from", xend = "to", 
         y = "yfrom", yend = "yto")) +
        geom_point(data = spp_point_df, colour = "black", shape = 1) +
        scale_y_log10() +
        scale_x_log10() +
      geom_abline(intercept = fg_line_intercept, slope = 2, linetype = 2, 
        colour = "darkgrey") +
      theme_bw() +
      theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
      theme(axis.line = element_line(color = 'black')) +
      labs(title = bquote(plain(.(plotTitle)))) +
      #labs(title = bquote(italic(.(plotTitle)))) +
      theme(plot.title = element_text(size = 9)) +
      theme(axis.text = element_text(size = 8)) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(axis.ticks.length = unit(-0.2, "cm")) +
      theme(axis.ticks.margin = unit(0.3, "cm"))
      #plot <- plot_base +
  if (x_axis_labels == TRUE) {
    plot1 <- plot_base + xlab("standard length, mm")
  } else if (x_axis_labels == FALSE) {
    plot1 <- plot_base + theme(axis.title.x = element_blank())
  }
  if (y_axis_labels == TRUE) {
    plot2 <- plot1 + ylab(expression(paste("gape area, ", mm^2, sep= "")))
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



allGA <- sma(ga~SL, data = pento, log = "xy", method = "SMA", robust = T, slope.test = 2)
#check_assump(allGA, "Pento Gape Area All")
allGA_summ <- mk_sma_summary(allGA, 1)
pento_line <- mk_sma_graph_df(allGA_summ, 1, "j_fg")

pisc_plot <- 
	mk_bootSMA_plots(fg_point_df = pento, spp_point_df = p, 
		spp_line_df_row = bootSMA_graph_df[1, ], x_axis_labels = FALSE, 
		y_axis_labels = FALSE, fg_line_intercept = pento_line$ref_intercept, 
		x_axis_text = TRUE, y_axis_text = TRUE, plot_title = "Piscivores") +
	  geom_text(data = bootSMA_eqns[1, ], aes(x = 570, y = 1.4, label = eqn), 
	  	parse = TRUE, size = 3, hjust = 1)
benth_plot <- 
	mk_bootSMA_plots(fg_point_df = pento, spp_point_df = b, 
		spp_line_df_row = bootSMA_graph_df[2, ], x_axis_labels = FALSE, 
		y_axis_labels = FALSE, fg_line_intercept = pento_line$ref_intercept, 
		x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Benthic Invertivores") +
	  geom_text(data = bootSMA_eqns[2, ], aes(x = 570, y = 1.4, label = eqn), 
	  	parse = TRUE, size = 3, hjust = 1)
zoop_plot <- 
	mk_bootSMA_plots(fg_point_df = pento, spp_point_df = z, 
		spp_line_df_row = bootSMA_graph_df[3, ], x_axis_labels = FALSE, 
		y_axis_labels = FALSE, fg_line_intercept = pento_line$ref_intercept, 
		x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Zooplanktivores") +
	  geom_text(data = bootSMA_eqns[3, ], aes(x = 570, y = 1.4, label = eqn), 
	  	parse = TRUE, size = 3, hjust = 1)
herb_plot <- 
	mk_bootSMA_plots(fg_point_df = pento, spp_point_df = h, 
		spp_line_df_row = bootSMA_graph_df[4, ], x_axis_labels = FALSE, 
		y_axis_labels = FALSE, fg_line_intercept = pento_line$ref_intercept, 
		x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Herbivores") +
	  geom_text(data = bootSMA_eqns[4, ], aes(x = 570, y = 1.4, label = eqn), 
	  	parse = TRUE, size = 3, hjust = 1)
coral_plot <- 
	mk_bootSMA_plots(fg_point_df = pento, spp_point_df = c, 
		spp_line_df_row = bootSMA_graph_df[5, ], x_axis_labels = FALSE, 
		y_axis_labels = FALSE, fg_line_intercept = pento_line$ref_intercept, 
		x_axis_text = TRUE, y_axis_text = FALSE, plot_title = "Corallivore") +
	  geom_text(data = bootSMA_eqns[5, ], aes(x = 570, y = 1.4, label = eqn), 
	  	parse = TRUE, size = 3, hjust = 1)
dev.new(height = 2.5, width = 9)
master_layout <- 
grid.layout(nrow = 2, ncol = 6, 
			widths = unit(c(0.1, 1, 0.9, 0.9, 0.9, 0.9), "null"), 
			heights = unit(c(1, 0.1), "null")
			)
grid.newpage()
pushViewport(viewport(layout = master_layout))
print(pisc_plot, vp = set_vp(1, 2))
print(benth_plot, vp = set_vp(1, 3))
print(zoop_plot, vp = set_vp(1, 4))
print(herb_plot, vp = set_vp(1, 5))
print(coral_plot, vp = set_vp(1, 6))
grid.text(
	expression( paste("log(gape area, ", mm^2, ")", sep = "") ), 
	vp = viewport(layout.pos.row = 1, layout.pos.col = 1),
	rot = 90, gp = gpar(fontsize = 9), 
	vjust = 1
	)
grid.text(
	"log(standard length, mm)",
	vp = viewport(layout.pos.row = 2, layout.pos.col = 4),
	vjust = -0.3, gp = gpar(fontsize = 9)
	)

#dev.copy2eps(device = quartz, file = "panel_plots/boot_fg_plot.eps")


