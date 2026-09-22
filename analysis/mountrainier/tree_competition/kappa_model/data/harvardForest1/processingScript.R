# started on Aug 23, 2026 by Devina
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tidyverse)
library(dbscan)

setwd("D:/ubc_study/udergrad_research/temporal ecology lab/climategrowthshifts/analysis/mountrainier/tree_competition/kappa_model/data/harvardForest1")

# processing species code, clean up individuals with unknown genus
# unify the unknown species to "spp."
species_code <- read_csv("hf253-02-species-codes.csv") %>%
  filter(genus!="Unidentified") %>%
  mutate(species = ifelse(species %in% c("unknown","unk_hardwood","unk_conifer","species"),
                          "spp.",species)) %>%
  mutate(latin = paste0(genus,"_",species))%>%
  select(sp,latin,genus,species) 

stemdataName <- list("hf253-05-stems-2014.csv","hf253-06-stems-2019.csv")
stem_list <- map(stemdataName,function(f){
  dat <- read_csv(f,show_col_types = FALSE)
  dat%>%
    filter(status=='A') %>%
    select(tree.id,stem.id,sp,gx,gy,dbh)
})

#treeid, stemid, sp, gx, gy, dbh
# alive status
# calculate tree growth from 2014-2019 as the radial growth
growthchange <- stem_list[[1]]%>%
  inner_join(stem_list[[2]],
             by = c("tree.id","stem.id"),
             suffix = c("_2014","_2019"))%>%
  filter(sp_2014==sp_2019) %>%
  mutate(growth = 10*(dbh_2019-dbh_2014)/2, # tree growth (mm)
         locationchange = sqrt((gx_2014-gx_2019)^2+(gy_2014-gy_2019)^2))%>%
  filter(locationchange==0)%>%
  select(tree.id,stem.id,sp_2014,gx_2014,gy_2014,growth,dbh_2019,dbh_2014)%>%
  rename(sp=sp_2014,gx=gx_2014,gy=gy_2014)%>%
  inner_join(species_code,by="sp")

# calculate the neighborhood in 2014 and 2019
neighborhood_list <- purrr::map(stem_list,function(d) {
  dat <- d %>% filter(!is.na(gx),!is.na(gy))
  coords <- as.matrix(dat[,c("gx","gy")])
  neighbors <- frNN(coords,eps=10,sort=TRUE)
  
  neighborhood <- map_dfr(seq_len(nrow(dat)),function(i){
    neighbor_idx <- neighbors$id[[i]]
    neighbor_idx <- neighbor_idx[neighbor_idx!=i]
    
    if (length(neighbor_idx)==0) {
      return (NULL)
    }
    
    tibble(
      focal_treeid = dat$tree.id[i],
      focal_stemid = dat$stem.id[i],
      
      neighbor_treeid = dat$tree.id[neighbor_idx],
      neighbor_stemid = dat$stem.id[neighbor_idx],
      neighbor_sp = dat$sp[neighbor_idx],
      neighbor_dbh = dat$dbh[neighbor_idx],
      
      dist = sqrt((dat$gx[neighbor_idx]-dat$gx[i])^2+
                    (dat$gy[neighbor_idx]-dat$gy[i])^2)
    )
    }) %>%
    filter(focal_treeid!=neighbor_treeid)
  
  return(neighborhood)
})

# find the latin name of neighboring species
neighbor_2014 <- neighborhood_list[[1]] %>%
  inner_join(species_code,by = c("neighbor_sp"="sp"))%>%
  rename(neighbor_latin = latin,
         neighbor_genus = genus,
         neighbor_species = species)
neighbor_2019 <- neighborhood_list[[2]] %>%
  inner_join(species_code,by = c("neighbor_sp"="sp")) %>%
  rename(neighbor_latin = latin,
         neighbor_genus = genus,
         neighbor_species = species)

# save processed files
write.csv(neighbor_2014,"processed/neighborhood_2014.csv")
write.csv(neighbor_2019,"processed/neighborhood_2019.csv")
write.csv(growthchange,"processed/tree_growth.csv")

species_2014 <- growthchange %>%
  left_join(neighbor_2014,by= c("tree.id"="focal_treeid",
                                "stem.id"="focal_stemid"))%>%
  mutate(ba = pi*(neighbor_dbh/2)^2) %>%
  group_by(tree.id,stem.id,neighbor_latin) %>%
  summarize(ba_sum = sum(ba)/1e3) %>% # basal area (*1000cm^2)
  left_join(growthchange) %>%
  mutate(radius_2014 = dbh_2014/20,
         radius_2019 = dbh_2019/20) %>% #radius unit (*10 cm)
  select(tree.id,stem.id,ba_sum,neighbor_latin,radius_2014,latin,genus,species,growth) %>%
  filter(growth > 0)%>%
  rename("speciesName"=species,
         "genusName"=genus)

species_2019 <- growthchange %>%
  left_join(neighbor_2019,by= c("tree.id"="focal_treeid",
                                "stem.id"="focal_stemid"))%>%
  mutate(ba = pi*(neighbor_dbh/2)^2) %>%
  group_by(tree.id,stem.id,neighbor_latin) %>%
  summarize(ba_sum = sum(ba)/1e3) %>% # basal area (*1000cm^2)
  left_join(growthchange) %>%
  mutate(radius_2014 = dbh_2014/20,
         radius_2019 = dbh_2019/20) %>% #radius unit (*10 cm)
  select(tree.id,stem.id,ba_sum,neighbor_latin,radius_2019,latin,genus,species,growth) %>%
  filter(growth > 0)%>%
  rename("speciesName"=species,
         "genusName"=genus)

write.csv(species_2014,"processed/species_2014.csv")
write.csv(species_2019,"processed/species_2019.csv")
