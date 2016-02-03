library(ape)
library(nlme)
library(phytools)

# Load data
fishTree <- read.tree('Jillian tree.txt')
mean_spp_summ <- read.csv('mean_spp_summ.csv')[, -1]
row.names(mean_spp_summ) <- mean_spp_summ$species

# Prepare tree
fishTree_no_out <- drop.tip(fishTree, "OUTGROUP_Myxine_glutinosa")

# Make tree ultrametric
chronos_tree1 <- chronos(fishTree_no_out, lambda = 1)

# Load function
lk <- function (sig2, y, X, C, v=NULL, opt=TRUE) {
    n <- nrow(C)
    
#    browser()

    if (is.null(v)) { v <- rep(0, n) }
    
    V <- sig2 * C + diag(v)
    beta <- solve(t(X) %*% solve(V) %*% X) %*% (t(X) %*% 
        solve(V) %*% y)
    logL <- -(1/2) * t(y - X %*% beta) %*% solve(V) %*% 
        (y - X %*% beta) - (1/2) * determinant(V, 
        logarithm = TRUE)$modulus - (n / 2) * log(2 * pi)
    
    if (opt == TRUE) { 
        return(-logL)
        } else { 
            return(list(beta = beta, sig2e = sig2, logL = logL))
        }
}

#------------------------------------------------------------------------------
# Comparison of pgls and custom function using simulated tree
#------------------------------------------------------------------------------
# Simulate tree
tree <- pbtree(n = 21, tip.label = as.character(mean_spp_summ$species), scale = 1)

# Create design matrix
design_mat <- model.matrix(~ fg - 1, data = mean_spp_summ)
slopes <- mean_spp_summ$slope

# New method
fit_lik <- optimize(lk, interval = c(0, 1000), y = slopes, 
                    X = design_mat, C = vcv(tree), v = NULL, 
                    opt = TRUE)
fitted_test <- lk(sig2 = fit_lik$minimum, y = slopes, 
             X = design_mat, C = vcv(tree), v = NULL, opt = FALSE)

# pgls
fit_gls <- gls(slope ~ fg - 1, data = mean_spp_summ, correlation=corBrownian(1, tree), method="ML")

# Compare these estimates to the gls with only the phylogenetic correlations
t(fitted_test$beta)
coef(fit_gls)
# These results are the same. 


#------------------------------------------------------------------------------
# Comparison of pgls and custom function using actual tree
#------------------------------------------------------------------------------
# The tree
chronos_tree1

# Create design matrix
design_mat <- model.matrix(~ fg - 1, data = mean_spp_summ)
slopes <- mean_spp_summ$slope

# New method
fit_lik <- optimize(lk, interval = c(0, 1000), y = slopes, 
                    X = design_mat, C = vcv(chronos_tree1), v = NULL, 
                    opt = TRUE)
fitted_test <- lk(sig2 = fit_lik$minimum, y = slopes, 
             X = design_mat, C = vcv(chronos_tree1), v = NULL, opt = FALSE)

# pgls
fit_gls <- gls(slope ~ fg - 1, data = mean_spp_summ, correlation=corBrownian(1, chronos_tree1), method="ML")

# Compare these estimates to the gls with only the phylogenetic correlations
t(fitted_test$beta)
coef(fit_gls)
# These results are different...