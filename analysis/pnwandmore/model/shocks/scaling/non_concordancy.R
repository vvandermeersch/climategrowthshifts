

par(mfrow=c(3,1), cex.lab = 1.2)
util$plot_expectand_pushforward(t(matrix(old_param_samples[1:1000,1:4,paste0('phi_sck[',s,']')], nrow = 1000, ncol = 4)), 50,
                                flim = c(0,1), ylim = c(0,20),
                                display_name= expression(phi[shock]),
                                col = '#27278f')
util$plot_expectand_pushforward(t(matrix(param_samples[1:1000,1:4,paste0('phi_sck[',s,']')], nrow = 1000, ncol = 4)), 50,
                                flim = c(0,1), ylim = c(0,20),
                                display_name=expression(phi[shock]),
                                col = '#278f27', add = TRUE)

util$plot_expectand_pushforward(t(matrix(old_param_samples[1:1000,1:4,paste0('omega_conc_sck[',s,']')], nrow = 1000, ncol = 4)), 50,
                                flim = c(0,1), ylim = c(0,20),
                                display_name= expression(omega[shock]^concordant),
                                col = '#27278f')
util$plot_expectand_pushforward(t(matrix(param_samples[1:1000,1:4,paste0('omega_conc_sck[',s,']')], nrow = 1000, ncol = 4)), 50,
                                flim = c(0,1), ylim = c(0,20),
                                display_name=expression(omega[shock]^concordant),
                                col = '#278f27', add = TRUE)


util$plot_expectand_pushforward(t(matrix(param_samples[1:1000,1:4,paste0('omega_nonconc_sck[',s,']')], nrow = 1000, ncol = 4)), 50,
                                flim = c(0,1), ylim = c(0,20),
                                display_name=expression(omega[shock]^{non-concordant}),
                                col = '#278f27')

test <- list(t(matrix(param_samples[1:1000,1:4,paste0('omega_conc_sck[12]')], nrow = 1000, ncol = 4)))
names(test) <- paste0('omega_conc_sck[12]')
util$check_all_expectand_diagnostics(test)

s <- 12
util$plot_pairs_by_chain(t(matrix(param_samples[1:1000,1:4,paste0('omega_conc_sck[',s,']')], nrow = 1000, ncol = 4)), expression(omega[sck]^concordant),
                         t(matrix(param_samples[1:1000,1:4,paste0('omega_nonconc_sck[',s,']')], nrow = 1000, ncol = 4)), expression(omega[sck]^{non-concordant}))


samples <- lapply(1:dim(param_samples)[3], function(k){t(matrix(param_samples[1:dim(param_samples)[1],1:dim(param_samples)[2],k], 
                                                                nrow = dim(param_samples)[1], ncol = dim(param_samples)[2]))})
names(samples) <- dimnames(param_samples)$variable
util$check_all_expectand_diagnostics(samples)


par(mfrow=c(1,1), cex.lab = 1.2)
util$plot_expectand_pushforward(t(matrix(param_samples[1:1000,1:4,paste0('omega_nonconc_sck[',1,']')], nrow = 1000, ncol = 4)), 50,
                                flim = c(0,0.2), ylim = c(0,200),
                                display_name=expression(omega[shock]^{non-concordant}),
                                col = '#278f27')
for(s in 1:data$N_stand_species){
  if(s == 12 | s == 1){next}
  util$plot_expectand_pushforward(t(matrix(param_samples[1:1000,1:4,paste0('omega_nonconc_sck[',s,']')], nrow = 1000, ncol = 4)), 50,
                                  flim = c(0,0.2), ylim = c(0,200),
                                  display_name=expression(omega[shock]^{non-concordant}),
                                  col = '#278f27', add = T)
}

