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

#-------------------------------------------------------------------------------
# Bootstrap results for piscivores without CA.MELA:
p_no_caranx <- p[!p$SpeciesCode == "CA.MELA", ]
m <- 1000
boot_noCA <- rlply(m, mk_all_spp_samp_df(p_no_caranx, b, z, h, c), .progress = "text")

# Do the SMA calculations (parallelised)
registerDoMC()
start = proc.time()
bootSMA_noCA <- llply(boot_noCA, function(x) {
	sma(ga ~ SL * j_fg, data = x, log = "xy", method = "SMA", robust = T, 
		slope.test = 2, multcomp = T, multcompmethod = "adjusted")
	}, .progress = "text", .parallel = TRUE)
end = proc.time()
total = (end[3]-start[3])
cat(paste("Time Elapsed: ", total, "\n"))

# Checking allometry
check_df <- ldply(bootSMA_noCA, function(x) slp_check(x))
allo_check_df <- ddply(check_df, .(j_fg), summarise, 
					   iso = sum(iso), neg = sum(neg), pos = sum(pos)
					   )
allo_check_df$j_fg <- factor(allo_check_df$j_fg, levels = c("Pi", "BI", "ZP", "He", "C"))
allo_check_df <- arrange(allo_check_df, j_fg)

# Calculating functional group slopes
boot_summ_noCA  <- ldply(bootSMA_noCA, mk_boot_summ_df)
boot_means_noCA <- ddply(boot_summ_noCA, .(j_fg), get_mean_slope_int_r2)

boot_means_noCA
  j_fg    slope        int         r2
1   Pi 1.761515 -0.7805252 0.83365686
2   BI 1.424071 -0.3908142 0.83492085
3   ZP 2.047717 -2.0882218 0.86691138
4   He 2.348532 -2.9849965 0.85529891
5    C 4.003380 -6.3401771 0.07056937

#-------------------------------------------------------------------------------
# Bootstrap results for piscivores without LU.BOHA:
p_no_bohar  <- p[!p$SpeciesCode == "LU.BOHA", ]
m <- 1000
boot_noBO <- rlply(m, mk_all_spp_samp_df(p_no_bohar, b, z, h, c), .progress = "text")

# Do the SMA calculations (parallelised)
registerDoMC()
start = proc.time()
bootSMA_noBO <- llply(boot_noBO, function(x) {
	sma(ga ~ SL * j_fg, data = x, log = "xy", method = "SMA", robust = T, 
		slope.test = 2, multcomp = T, multcompmethod = "adjusted")
	}, .progress = "text", .parallel = TRUE)
end = proc.time()
total = (end[3]-start[3])
cat(paste("Time Elapsed: ", total, "\n"))

# Checking allometry
check_df <- ldply(bootSMA_noBO, function(x) slp_check(x))
allo_check_df <- ddply(check_df, .(j_fg), summarise, 
					   iso = sum(iso), neg = sum(neg), pos = sum(pos)
					   )
allo_check_df$j_fg <- factor(allo_check_df$j_fg, levels = c("Pi", "BI", "ZP", "He", "C"))
allo_check_df <- arrange(allo_check_df, j_fg)

# Calculating functional group slopes
boot_summ_noBO  <- ldply(bootSMA_noBO, mk_boot_summ_df)
boot_means_noBO <- ddply(boot_summ_noBO, .(j_fg), get_mean_slope_int_r2)
boot_means_noBO
  j_fg    slope        int         r2
1   Pi 1.510585 -0.2548429 0.69587067
2   BI 1.431133 -0.4320715 0.80220320
3   ZP 1.988395 -1.9722458 0.85620520
4   He 2.440499 -3.1469848 0.87890711
5    C 4.003380 -6.3401771 0.07056937

#-------------------------------------------------------------------------------
# Bootstrap results for piscivores without LU.BOHA:
p_no_boca  <- p[!p$SpeciesCode %in% c("LU.BOHA", "CA.MELA"), ]
m <- 1000
boot_noBOCA <- rlply(m, mk_all_spp_samp_df(p_no_boca, b, z, h, c), .progress = "text")

# Do the SMA calculations (parallelised)
registerDoMC()
start = proc.time()
bootSMA_noBOCA <- llply(boot_noBOCA, function(x) {
	sma(ga ~ SL * j_fg, data = x, log = "xy", method = "SMA", robust = T, 
		slope.test = 2, multcomp = T, multcompmethod = "adjusted")
	}, .progress = "text", .parallel = TRUE)
end = proc.time()
total = (end[3]-start[3])
cat(paste("Time Elapsed: ", total, "\n"))

# Checking allometry
check_df <- ldply(bootSMA_noBOCA, function(x) slp_check(x))
allo_check_df <- ddply(check_df, .(j_fg), summarise, 
					   iso = sum(iso), neg = sum(neg), pos = sum(pos)
					   )
allo_check_df$j_fg <- factor(allo_check_df$j_fg, levels = c("Pi", "BI", "ZP", "He", "C"))
allo_check_df <- arrange(allo_check_df, j_fg)

# Calculating functional group slopes
boot_summ_noBOCA  <- ldply(bootSMA_noBOCA, mk_boot_summ_df)
boot_means_noBOCA <- ddply(boot_summ_noBOCA, .(j_fg), get_mean_slope_int_r2)
boot_means_noBOCA
  j_fg    slope        int         r2
1   Pi 1.686903 -0.6008857 0.83836460
2   BI 1.485623 -0.5418033 0.76484747
3   ZP 2.040010 -2.0714786 0.86727724
4   He 2.733723 -3.8600253 0.94047144
5    C 4.003380 -6.3401771 0.07056937

# Timing plyr function with .parallel = TRUE to see how much faster it is.
registerDoMC()
start = proc.time()
parallelSMA <- llply(parallel_test, function(x) {
	sma(ga ~ SL * j_fg, data = x, log = "xy", method = "SMA", robust = T, 
	  slope.test = 2, multcomp = T, multcompmethod = "adjusted") 
	}, .progress = "text", .parallel = TRUE)
end = proc.time()
total = (end[3]-start[3])
cat(paste("Time Elapsed: ", total, "\n"))

# Hadley Wickham answer on SO!!! 
# How to sub sample data by group efficiently using plyr
# http://stackoverflow.com/questions/16912186/how-do-i-sub-sample-data-by-group-efficiently
sample_by_spp <- function(df, dfgroup, size, replace) {
  ddply(df, dfgroup, function(x) {
    x[sample(nrow(x), size = size, replace = replace), , drop = FALSE]
  })
}

mk_all_spp_samp_df <- function(){
  p_samp <- sample_by_spp(df = p, dfgroup = 'SpeciesCode', size = 12, replace = TRUE)
  b_samp <- sample_by_spp(df = b, dfgroup = 'SpeciesCode', size = 24, replace = TRUE)
  z_samp <- sample_by_spp(df = z, dfgroup = 'SpeciesCode', size = 8, replace = TRUE)
  h_samp <- sample_by_spp(df = h, dfgroup = 'SpeciesCode', size = 7, replace = TRUE)
  all_spp_samp <- rbind.fill(p_samp, b_samp, z_samp, h_samp, c)
  return(all_spp_samp)
}

m <- 10 
m <- 1000
m <- 10000  ## bootstrap replicates to use for offical MS values
bootX2 <- rlply(m, mk_all_spp_samp_df(), .progress = "text")
ddply(x, .(j_fg), summarise, min(SL))

bootX3 <- rlply(m, mk_all_spp_samp_df(), .progress = "text")
bootSMA_test <- llply(bootX3, function(x) {
	sma(ga ~ SL * j_fg, data = x, log = "xy", method = "SMA", robust = T, 
		slope.test = 2, multcomp = T, multcompmethod = "adjusted")
	}, .progress = "text")

#-------------------------------------------------------------------------------
# The bootstrapping

m <- 10000  ## bootstrap replicates to use for offical MS values
bootX <- rlply(m, mk_all_spp_samp_df(), .progress = "text")

bootSMA <- llply(bootX, function(x) {
	sma(ga ~ SL * j_fg, data = x, log = "xy", method = "SMA", robust = T, 
		slope.test = 2, multcomp = T, multcompmethod = "adjusted")
	}, .progress = "text")

m <- 1000
bootX2 <- rlply(m, mk_all_spp_samp_df(), .progress = "text")
bootSMA2 <- llply(bootX2, function(x) {
	sma(ga ~ SL * j_fg, data = x, log = "xy", method = "SMA", robust = T, 
		slope.test = 2, multcomp = T, multcompmethod = "adjusted")
	}, .progress = "text")



# for each bootSMA and each functiona group in that bootSMA determine whether the 
# sma regression was negatively allometric, isometric, or positively allometric
slp_check <- function(x, slp_check_df) {
  for (i in 1:5) {
    if (x$groupsummary$Slope_test_p[i] > 0.05) {
      slp_check_df$iso[i] <<- slp_check_df$iso[i] + 1
    } else {
        if (x$groupsummary$Slope_lowCI[i] > 2) {
          slp_check_df$pos[i] <<- slp_check_df$pos[i] + 1
        } else {
            if (x$groupsummary$Slope_highCI[i] < 2) {
              slp_check_df$neg[i] <<- slp_check_df$neg[i] + 1
          }
        }
      }
   }
}

slp_check_df <- data.frame(j_fg = c("Pi", "BI", "ZP", "He", "C"), iso = 0, 
						   neg = 0, pos = 0)

ldply(bootSMA_noCA, function(x) slp_check(x, slp_check_df))

# For each bootSMA and the multiple pairwise comparison between functional groups
# count if there was a significant difference
multcomp_check <- function(x, multcomp_df) {
	for (i in 1:length(x$multcompresult[, 1])) {
		if (x$multcompresult$Pval[i] <= 0.05) {
			multcomp_df$sig_diff[i] <<- multcomp_df$sig_diff[i] + 1
		}
	}
}
multcomp_df <- bootSMA[[1]]$multcompresult[, 1:2]
multcomp_df$sig_diff <- 0
lapply(bootSMA, function(x) multcomp_check(x, multcomp_df))
multcomp_df$prop <- multcomp_df$sig_diff / 1000


#-------------------------------------------------------------------------------
# Making graphing dataframe for averaged SMA bootstrap

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

boot_summ_df <- ldply(bootSMA, mk_boot_summ_df)

boot_graph_df <- ddply(boot_summ_df, .(group), mk_bootSMA_eqn) 

# Mean slope and intercept values
boot_slp_int_r2_df <- ddply(boot_summ_df, .(group), get_mean_slope_int_r2) 
names(boot_slp_int_r2_df)[1] <- "j_fg" 
boot_slp_int_r2_df$j_fg <- factor(boot_slp_int_r2_df$j_fg, levels = c("Pi", "BI", "ZP", "He", "C"))
boot_slp_int_r2_df <- arrange(boot_slp_int_r2_df, j_fg)

boot_summ_df_test <- ldply(bootSMA_test, mk_boot_summ_df)
boot_graph_df_test <- ddply(boot_summ_df_test, .(group), mk_bootSMA_eqn) 

# Get the mean slope for each functional group from the bootSMA list
get_mean_slope_int_r2 <- function(x) {
	mean_slope <- mean(x$slope)
	mean_int <- mean(x$int)
	mean_r2 <- mean(x$r2)
	df <- data.frame(slope = mean_slope, int = mean_int, r2 = mean_r2)
	return(df)
}
# actually getting the mean slope and int df
boot_summ_df <- ldply(bootSMA, mk_boot_summ_df)
noCA_boot_summ_df <- ldply(bootSMA_noCA, mk_boot_summ_df)

# Use the overal functional group to find the midpoint of the cluster of points
allGA <- sma(ga~SL, data = pento, log = "xy", method = "SMA", robust = T, slope.test = 2)
#check_assump(allGA, "Pento Gape Area All")
allGA_summ <- mk_sma_summary(allGA, 1)
allGA_graph_df <- mk_sma_graph_df(allGA_summ, 1, "j_fg")
names(allGA_graph_df)

# SMA regression for each functional group
all_fg_GA <- sma(ga~SL * j_fg, data = pento, log = "xy", method = "SMA", robust = T, 
	slope.test = 2, multcomp = T, multcompmethod = "adjusted")
#check_assump(allGA, "Pento Gape Area All")
all_fg_GA_summ <- mk_spp_summary(all_fg_GA, 5, grouping=TRUE)
all_fg_GA_graph_df <- mk_smaSPP_graph_df(all_fg_GA_summ, 5, "j_fg")

bootSMA_graph_df <- all_fg_GA_graph_df
bootSMA_graph_df$slp <- boot_slp_int_df$slope
bootSMA_graph_df$int <- boot_slp_int_df$int
bootSMA_graph_df$r2  <- boot_slp_int_df$r2

ddply()

# Get the yfrom and yto values
to_from_df <- ddply(pento, .(j_fg), summarise, from = min(SL), to = max(SL))

line_data <- merge(x = to_from_df, y = boot_slp_int_df, by.x = "j_fg", by.y = "group")
line_data <- arrange(line_data, j_fg)

bootSMA_NEWgraph_df <- get_y_and_midpoint(line_data)


noCA_boot_graph_df <- get_y_FromTo(noCA_boot_summ_df, 100)
names(noCA_boot_graph_df)[1] <- "j_fg"
noCA_boot_graph_df$j_fg <- factor(noCA_boot_graph_df$j_fg, levels = c("Pi", "BI", "ZP", "He", "C"))



ggplot(data = pento, aes(x = SL, y = ga)) +
  geom_segment(data = noCA_boot_graph_df, aes(x = from, xend = to, y = yfrom, yend = yto), alpha = 0.2) +
  geom_abline(data = all_fg_GA_graph_df, aes(intercept = ref_intercept, slope = 2), colour = "grey") +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(~ j_fg, scales = "free")

get_y_FromTo <- function(x, num_of_rows) {
	for (i in 1:num_of_rows) {
	yfrom = 10^(x$slope * log10(x$from) + x$int)
	yto   = 10^(x$slope * log10(x$to) + x$int)
	midpoint_y = sqrt(yfrom * yto)
	midpoint_x = sqrt(x$from * x$to)
	ref_intercept = log10(midpoint_y / midpoint_x^2)
	x$yfrom = yfrom
	x$yto = yto
	x$ref_intercept = ref_intercept
	}
	return(x)
}


get_line_data <- function(x) {
	from <- min(SL)
	to   <- max(SL)
	yfrom = 10^(mean_slope * log10(mean_from) + mean_int)
	yto   = 10^(mean_slope * log10(mean_to) + mean_int)
	midpoint_y = sqrt(yfrom * yto)
	midpoint_x = sqrt(mean_from * mean_to)
	ref_intercept = log10(midpoint_y / midpoint_x^2)
	df <- data.frame(from = mean_from, to = mean_to, slope = mean_slope, 
					 int = mean_int, yfrom = yfrom, yto = yto, 
					 midpoint_x = midpoint_x, midpoint_y = midpoint_y, 
					 ref_intercept = ref_intercept)
	return(df)
}


mk_bootSMA_eqn <- function(x) {
	mean_from <- mean(x$from)
	mean_to <- mean(x$to)
	mean_slope <- mean(x$slope)
	mean_int <- mean(x$int)
	yfrom = 10^(mean_slope * log10(mean_from) + mean_int)
	yto   = 10^(mean_slope * log10(mean_to) + mean_int)
	midpoint_y = sqrt(yfrom * yto)
	midpoint_x = sqrt(mean_from * mean_to)
	ref_intercept = log10(midpoint_y / midpoint_x^2)
	df <- data.frame(from = mean_from, to = mean_to, slope = mean_slope, 
					 int = mean_int, yfrom = yfrom, yto = yto, 
					 midpoint_x = midpoint_x, midpoint_y = midpoint_y, 
					 ref_intercept = ref_intercept)
	return(df)
}

mk_bootSMA_eqn <- function(x) {
	mean_from <- mean(x$from)
	mean_to <- mean(x$to)
	mean_slope <- mean(x$slope)
	mean_int <- mean(x$int)
	yfrom = 10^(mean_slope * log10(mean_from) + mean_int)
	yto   = 10^(mean_slope * log10(mean_to) + mean_int)
	midpoint_y = sqrt(yfrom * yto)
	midpoint_x = sqrt(mean_from * mean_to)
	ref_intercept = log10(midpoint_y / midpoint_x^2)
	df <- data.frame(from = mean_from, to = mean_to, slope = mean_slope, 
					 int = mean_int, yfrom = yfrom, yto = yto, 
					 midpoint_x = midpoint_x, midpoint_y = midpoint_y, 
					 ref_intercept = ref_intercept)
	return(df)
}


write_bootSMA_eqns <- function(sma_summary_df, group_column) {
  dfm <- data.frame("j_fg" = rep(NA, 5), "eqn" = rep(NA, 5))
  m = matrix(data=NA, nrow=0, ncol=5)
  count <- length(group_column)
  for (i in (1:count)) {
    l <- list(slp = format(sma_summary_df$slp[i], digits = 2, nsmall = 2),
              int = format(sma_summary_df$int[i], digits = 2, nsmall = 2), 
              r2  = "NA",
              count   = "NA"
    )
    if (l$int >= 0) {
      eqn_r2 <- substitute(atop(italic(r)^2~"="~r2, italic(y) ==
                          slp%.%italic(x) + int), l)
      eqn <- substitute(italic(y) == slp*italic(x) + int, l)
      r2  <- substitute(italic(r)^2~"="~r2, l)
      n   <- substitute(italic(n) ~ "=" ~ count, l)
    } else {
        l <- list(slp = format(sma_summary_df$slp[i], digits = 2, nsmall = 2),
                  int = format(abs(sma_summary_df$int[i]), digits = 2, nsmall = 2),
                  r2  = "NA",
                  count   = "NA"
      )
      eqn_r2 <- substitute(atop(italic(r)^2~"="~r2, italic(y) ==
                          slp%.% italic(x) - int), l)
      eqn <- substitute(italic(y) == slp*italic(x) - int, l)
      r2  <- substitute(italic(r)^2~"="~r2, l)
      n   <- substitute(italic(n) ~ "=" ~ count, l)
    }    
    #browser()
    eqn_r2 <- as.character(as.expression(eqn_r2)) 
    eqn    <- as.character(as.expression(eqn))
    r2     <- as.character(as.expression(r2))
    n      <- as.character(as.expression(n))
    dfm$j_fg[i] <- as.character(sma_summary_df[i, 1])
    dfm$eqn[i]  <- eqn
    #m <- rbind(m, c(as.character(df[i,1]), lm_eq))
  }
  return(dfm)
}

bootSMA_eqns <- write_bootSMA_eqns(sma_summary_df = bootSMA_graph_df, 
								   group_column = j_fg)


bootSMA_eqns$j_fg <- factor(bootSMA_eqns$j_fg, levels = c("Pi", "BI", "ZP", "He", "C"))

y <- arrange(bootSMA_eqns, j_fg)

y <- arrange(boot_summ_df, group, from)



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

dev.copy2eps(device = quartz, file = "panel_plots/boot_fg_plot.eps")


