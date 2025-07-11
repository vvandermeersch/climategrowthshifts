

ggplot(phy.plants.here) + 
  ggtree::geom_tree() +
  ggtree::geom_tiplab(as_ylab=TRUE, color='firebrick')


phy.plants.here$tip.label




td <- data.frame(node = phy.plants.here$node.label,
                 trait = NA)
nd <- data.frame(node = phy.plants.here$tip.label, trait = medrhos$rho)

d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(phy.plants.here, d, by = 'node')

color <- c(medrhos$rho, rep(NA, 29))

color <- 1:59


phy.plants.here

phy.plants.here$rho <-  medrhos$rho

kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF",
                  "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")


estimates <- 
  data.frame(
    rho = sapply(1:data$N_species, function(s) median(samples[[paste0('rho_sp[',s,']')]])),
    kappa = sapply(1:data$N_species, function(s) median(samples[[paste0('kappa_sh[',s,']')]])),
    shortname = uniq_species_ids
  )
estimates <- merge(estimates, sppfull[,c('shortname', 'phylo.name')])
names(estimates)[names(estimates)=='phylo.name'] <- 'label'
tree <- ggtree(phy.plants.here, size = 0.5) 
datatree <- tree$data
datatree <- merge(datatree, estimates, all.x = TRUE)

ggtree(datatree) +
  geom_tiplab(aes(color=rho)) +
  ggplot2::scale_color_gradientn(colors = kippenberger, limits = c(4,32)) +
  xlim(c(0,200))

ggplot(data = datatree) +
  geom_point(aes(x = rho, y = y))
  
  # ggplot2::scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  geom_tiplab(aes(color=color)) +
  xlim(c(0,200)) +
  ggplot2::scale_color_gradientn(colors = kippenberger, limits = c(4,32))
