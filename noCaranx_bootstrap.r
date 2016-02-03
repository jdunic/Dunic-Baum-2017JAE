
# Piscivore slope estimate without Caranx
m <- 10000
p_no_caranx <- p[!p$SpeciesCode == "CA.MELA", ]
boot_noCA <- rlply(m, sample_by_spp(df = p_no_caranx, dfgroup = 'SpeciesCode', 
									size = 12, replace = TRUE))
registerDoMC()
start = proc.time()
bootSMA_noCA <- llply(boot_noCA, function(x) {
	sma(ga ~ SL, data = x, log = "xy", method = "SMA", robust = T, 
		slope.test = 2)
	}, .parallel = TRUE)
end = proc.time()
total = (end[3]-start[3])
cat(paste("Time Elapsed: ", total, "\n"))

save(boot_noCA, bootSMA_noCA, file = "~/R_projects/Allometry/bootSMA_noCA_10000.RData")

boot_p_noCA_slopes <- function(x) {
	slope <- x$groupsummary$Slope
	int <- x$groupsummary$Int[1]
	r2  <- x$groupsummary$r2[1]
	df <- data.frame(j_fg = "Pi", slope = slope, int = int, r2 = r2)
	return(df)
}
boot_noCA_slopes <- ldply(bootSMA_noCA, function(x) boot_p_noCA_slopes(x))

# Checking for symmetric distribution of slope estimates
ggplot(data = boot_noCA_slopes, aes(x = slope)) +
  geom_histogram()

boot_summ_noCA <- get_mean_slope_int_r2(boot_noCA_slopes)
boot_summ_noCA$lowerCI <- boot_summ_noCA$slope - (1.96 * boot_summ_noCA$slope_se)
boot_summ_noCA$upperCI <- boot_summ_noCA$slope + (1.96 * boot_summ_noCA$slope_se)
boot_summ_noCA


  j_fg    slope        int         r2   slope_se  lowerCI  upperCI
1   Pi 1.531990 -0.3067777 0.67674821 0.09725793 1.341365 1.722616
2   BI 1.521871 -0.6038637 0.82719876 0.08753882 1.350295 1.693447
3   ZP 2.042337 -2.0824801 0.85913588 0.08884415 1.868203 2.216472
4   He 2.449285 -3.2007591 0.89168539 0.10539008 2.242721 2.655850
5    C 4.003380 -6.3401771 0.07056937 0.00000000 3.040136 5.271820

> boot_summ_noCA
     slope        int        r2   slope_se  lowerCI  upperCI
1 1.752957 -0.7699998 0.8460138 0.08614923 1.584104 1.921809