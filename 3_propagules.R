##############################
##### EVANS ET AL. 2022. #####
##############################

#################################
#####Propagule analysis #########
#################################

##load propagule shapefile

library(dplyr)
sp_prop<-read.csv("sp_list_prop.csv")
prop_dat<-read.csv("prop-dat.csv")

uni_prop<-unique(sp_prop$Spp)


n_prop <- unique(sp_prop$Spp) %>%
  length()

tester<-NULL;
for (i in 1:385208){
  newdat<-prop_dat[sample(1:nrow(prop_dat), 1), ]
  tester<-rbind(tester, newdat)
  #sp_groups<-sp_prop%>%group_by(Spp)
  # rand_df<-prop_dat[sample(nrow(prop_dat), size = prop_tmp)]
  # tester<-rbind(tester, rand_dat)
}

tester1<-tester[-c(1:20),]

prop_combo<-cbind(sp_prop, tester1)

write.csv(prop_combo, "Propgule_random_sp.csv")


library(rgdal)
summary(species_loc)
props<-readOGR("propagules.shp")
props$location.l->props$y
props$location_1->props$x
sp_dat<-as.data.frame(props)
summary(sp_dat)
sp_dat<-sp_dat[-c(1:2, 5:6)]
write.csv(sp_dat, "props.csv")
plot(prop_combo, add = T)
values<-extract(world, tester1)
prop_bind<-cbind.data.frame(prop_combo, values)
prop_bind
write.csv(prop_bind, "Propagule_bind.csv")


names(prop_bind)[names(prop_bind)=="bio1"] <- "tmp"
names(prop_bind)[names(prop_bind)=="bio5"] <- "tmax"
names(prop_bind)[names(prop_bind)=="bio6"] <- "tmin"
names(prop_bind)[names(prop_bind)=="bio12"] <- "p"

prop_bind<-prop_bind[-c(8:15)]
write.csv(prop_bind, "prop_clean.csv")

prop<-read.csv("prop_clean.csv")

names(prop)
summary(prop)

prop_select<- prop %>%
  group_by(spp_clean)%>%
  dplyr::summarize(sum=n())%>%
  filter(sum>=5)

prop<-prop%>%
  filter(spp_clean%in%prop_select$spp_clean)

nat<-read.csv("natbind.csv")

names(nat)[names(nat)=="bio1"] <- "tmp"
names(nat)[names(nat)=="bio5"] <- "tmax"
names(nat)[names(nat)=="bio6"] <- "tmin"
names(nat)[names(nat)=="bio12"] <- "p"

summary(nat)
nat<-nat[-c(8:14)]


nat<-na.omit(nat)
prop<-na.omit(prop)
nat_select<-nat%>%
  group_by(spp_clean)%>%
  dplyr::summarize(sum=n())%>%
  filter(sum>=5)

nat<-nat%>%
  filter(spp_clean%in%nat_select$spp_clean)
write.csv(nat, "nat.csv")
n_sp <- unique(nat$spp_clean) %>%
  length()


uni_nat<-unique(nat$spp_clean)
uni_prop<-unique(prop$spp_clean)
i
summary(nat)
summary(prop)
library(ecospat)
sp_overlap<-NULL;
pca_mapping<-NULL;

for (i in 1:n_sp) {
  prop_tmp<-prop %>% 
    filter(spp_clean==uni_prop[i])
  nat_tmp<-nat %>%
    filter(spp_clean==uni_nat[i])
  
  pca.env<-dudi.pca(rbind(nat_tmp, prop_tmp)[,4:7],scannf=F,nf=2)
  #pca_mapping<-rbind(pca.env, prop_tmp$x, prop_tmp$y)
  #pca.env<-na.omit(pca.env)
  ecospat.plot.contrib(contrib = pca.env$co, eigen=pca.env$eig)
  globclim<-pca.env$li
  #globclim<-globclim[-2,]
  scores.sp.nat<-suprow(pca.env,nat_tmp[which(nat_tmp[,8]==1),4:7])$li
  scores.clim.nat <- suprow(pca.env,nat_tmp[,4:7])$li
  scores.sp.inv<-suprow(pca.env,prop_tmp[which(prop_tmp[,8]==1),4:7])$li
  scores.clim.inv <- suprow(pca.env,prop_tmp[,4:7])$li
  #globclim<-globclim[complete.cases(globclim),]
  print(i)
  grid.clim.nat <- ecospat.grid.clim.dyn(globclim,
                                         scores.clim.nat,
                                         scores.sp.nat, R=100,
                                         th.sp=0)
  
  grid.clim.inv <- ecospat.grid.clim.dyn(globclim,
                                         scores.clim.inv,
                                         scores.sp.inv, R=100,
                                         th.sp=0)
  D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor = TRUE)$D
  print(D.overlap)
  tmp_D<-D.overlap
  sp_overlap<-rbind(sp_overlap, tmp_D)
}


propagule_sp_overlap_future_dat<-as.data.frame(sp_overlap)
write.csv(propagule_sp_overlap_future_dat, "Propagule_D_overlap_all_sp.csv")

prop_d<-read.csv("Propagule_D_overlap_all_sp.csv")

plot(prop_d$x, prop_d$V1)


##### Propgaules - future climate conditions####





