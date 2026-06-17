# 
# n <- 1e5
# 
# thetas_baseline <- gtools::rdirichlet(n, c(95, 3, 1, 1));
# 
# 
# omega_thetas <- abs(rnorm(n, 0, 0.05))
# 
# alphas = thetas_baseline / omega_thetas + c(1,1,1,1)
# 
# 
# thetas_idio <- lapply(1:nrow(alphas), function(i){gtools::rdirichlet(1, alphas[i,])})
# thetas_idio <- do.call(rbind, thetas_idio)
# 
# par(mfrow=c(2,2))
# hist(thetas_baseline[,1], xlim = c(0,1), breaks = seq(0,1,0.01))
# hist(thetas_baseline[,2], xlim = c(0,1), breaks = seq(0,1,0.01))
# hist(thetas_baseline[,3], xlim = c(0,1), breaks = seq(0,1,0.01))
# hist(thetas_baseline[,4], xlim = c(0,1), breaks = seq(0,1,0.01))
# 
# par(mfrow=c(2,2))
# hist(thetas_idio[,1], xlim = c(0,1), ylim = c(0, 15000), breaks = seq(0,1,0.01))
# text(x = 0.2, y = 13000, labels = paste(round(quantile(thetas_idio[,1], c(0.05, 0.95)), 1), collapse =' - '))
# hist(thetas_idio[,2], xlim = c(0,1), ylim = c(0, 15000), breaks = seq(0,1,0.01))
# text(x = 0.2, y = 13000, labels = paste(round(quantile(thetas_idio[,2], c(0.05, 0.95)), 1), collapse =' - '))
# hist(thetas_idio[,3], xlim = c(0,1), ylim = c(0, 15000), breaks = seq(0,1,0.01))
# text(x = 0.2, y = 13000, labels = paste(round(quantile(thetas_idio[,3], c(0.05, 0.95)), 1), collapse =' - '))
# hist(thetas_idio[,4], xlim = c(0,1), ylim = c(0, 15000), breaks = seq(0,1,0.01))
# text(x = 0.2, y = 13000, labels = paste(round(quantile(thetas_idio[,4], c(0.05, 0.95)), 1), collapse =' - '))
# 
# 
# 

n <- 1e5
thetas_baseline <- c(2, 100, 10)
omega_thetas <- 6
alphas = thetas_baseline / omega_thetas + c(1, 1, 1)
thetas_conc <- gtools::rdirichlet(n, alphas)

par(mfrow=c(2,2))
hist(thetas_conc[,1], xlim = c(0,1), ylim = c(0, 15000), breaks = seq(0,1,0.01))
text(x = 0.2, y = 13000, labels = paste(round(quantile(thetas_conc[,1], c(0.05, 0.95)), 2), collapse =' - '))
hist(thetas_conc[,2], xlim = c(0,1), ylim = c(0, 15000), breaks = seq(0,1,0.01))
text(x = 0.2, y = 13000, labels = paste(round(quantile(thetas_conc[,2], c(0.05, 0.95)), 2), collapse =' - '))
hist(thetas_conc[,3], xlim = c(0,1), ylim = c(0, 15000), breaks = seq(0,1,0.01))
text(x = 0.2, y = 13000, labels = paste(round(quantile(thetas_conc[,3], c(0.05, 0.95)), 2), collapse =' - '))

# hist(thetas_idio[,4], xlim = c(0,1), ylim = c(0, 15000), breaks = seq(0,1,0.01))
# text(x = 0.2, y = 13000, labels = paste(round(quantile(thetas_idio[,4], c(0.05, 0.95)), 2), collapse =' - '))

