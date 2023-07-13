##############################
##### EVANS ET AL. 2022. #####
##############################

install.packages("pacman")

pacman::p_load(rgbif, 
               raster, 
               sf, 
               stringr, 
               dplyr, 
               tidyverse, 
               spocc, 
               read.table, 
               readr,
               doParallel, 
               rgdal, 
               ecospat, 
               ggplot2,
               taxizedb, 
               terra, 
               taxise)

###rgbif location extraction####

install.packages("rgbif")
library(raster)
library(sf)
library(stringr)
spec_tax<-read.csv("Specific_taxonomy.csv")
sp_lst<-list(df_species_unique)
library(rgbif)
library(dplyr)
??rgbif
df_subset
devtools::install_github("ropensci/rgbif")

keys <- sapply(sp_list, fun ction(x) name_suggest(x)$key[1], USE.NAMES=T)

mydata<-occ_data(taxonKey=keys, hasCoordinate=TRUE)


occ_search(taxonKey = keys, limit = 5) 


specific<-read.csv("Specific_taxonomy.csv")

specific_species<-unique(remove_error$species)
write.csv(specific_species, "spec_sp.csv")
remove_error<-specific[!grepl('sp.', specific$species),]

library(tidyverse)
library(spocc)
species<-sp_list
# Species to serche for
species<-c(df1)
# empty data frame
df<- data.frame(x=NA, longitude=NA, latitude=NA)
df2<-as.list(df1)
df2$?..species
# species to loop over 
sp_list<- unique(species$species)

df1<-read.csv("spec_sp.csv")

for (i in df2$x){
  out<-spocc::occ(query=i , from="gbif", has_coords = T, limit=500)
  print(i)
  out_df<-spocc::occ2df(out)
  write.csv(out_df, paste0(i,'.csv'))
  #out_df <- dplyr::select(out_df, out$gbif$data$name, out$gbif$data$longitude, out$gbif$data$latitude)
  #print(out_df)
  #complete<-dplyr::bind_rows(df, out_df) %>% na.omit()
}

library(dplyr)
library(data.table)
library(readr)
csvs<-list.files(path=".", full.names = T)



install.packages("pacman")
library(pacman)
p_load(doParallel,data.table,stringr)

# get the file name
dir() %>% str_subset("\\.csv$") -> fn

# use parallel setting
(cl = detectCores() %>% 
    makeCluster()) %>% 
  registerDoParallel()

# read and bind
system.time({
  big_df = foreach(i = fn,
                   .packages = "data.table") %dopar% {
                     fread(i,colClasses = "chracter")
                   } %>% 
    rbindlist(fill = T)
})

# end of parallel work
stopImplicitCluster(cl)

write.csv(big_df, "Species_ranges_combined.csv")


#### add bioclim variables #####
library(raster)
library(sf)
dir() %>% str_subset("\\.tif$") -> ras

stack<-stack(ras)

loc<-st_as_sf(big_df, coords = c("latitude", "longitude"), crs = WGS84)  
  
world<-getData("worldclim",var="bio", res = 5)
plot(world[[1]])
world_rast<-subset(world, c(1,5:6,12))
summary(world_rast)
writeRaster(world_rast[[1]], "temp.tif")
writeRaster(world_rast[[2]], "tmax.tif")
writeRaster(world_rast[[3]], "tmin.tif")
writeRaster(world_rast[[4]], "precip.tif")


library(rgdal)
summary(species_loc)
species_loc<-readOGR("Species_range_clean_new.shp")
sp_dat<-as.data.frame(species_loc)
values<-extract(world, species_loc)
bind<-cbind.data.frame(species_loc, values)
bind
write.csv(bind, "species_loc_bioclim_clean.csv")
bind<-read.csv("species_loc_bioclim_clean.csv")
##### ADD biolcim to random points
rand<-readOGR("Random_point_clean.shp")

randdat<-as.data.frame(rand)
randvalues<-extract(world, rand)
bindrand2<-cbind.data.frame(randdat, rand, randvalues)
bindrand2
write.csv(bindrand2, "random_clean.csv")
bindrand<-read.csv("random_clean.csv")
tail(bindrand)
###### tidy data for ease of niche overlap analysis####

summary(bind)
bind<-write.csv(bind, "bind.csv")
bind<-read.csv("bind.csv")
#drop columns 
bind<-bind[-c(8:14)]
bind<-bind[-c(3, 5:7, 10:14, 16:22)]
bind<-bind[-c(4:9)]
bind<-bind[-c(7:13)]
#rename columns
names(bind)[names(bind)=="coords.x1"] <- "x"
names(bind)[names(bind)=="coords.x2"] <- "y"

names(bind)[names(bind)=="bio1"] <- "tmp"
names(bind)[names(bind)=="bio5"] <- "tmax"
names(bind)[names(bind)=="bio6"] <- "tmin"
names(bind)[names(bind)=="bio12"] <- "p"

summary(bindrand)
#rename columns random
names(bindrand)[names(bindrand)=="coords.x1"] <- "x"
names(bindrand)[names(bindrand)=="coords.x2"] <- "y"

names(bindrand)[names(bindrand)=="bio1"] <- "tmp"
names(bindrand)[names(bindrand)=="bio5"] <- "tmax"
names(bindrand)[names(bindrand)=="bio6"] <- "tmin"
names(bindrand)[names(bindrand)=="bio12"] <- "p"

bindrand<-bindrand[-c(8:14)]

bindrand<-bindrand[-c(5:7, 10:14, 16:22)]
bindrand<-bindrand[-c(6:10, 12:18)]

#### Make spatial

wgs84CRS<-st_crs(bind)
wgs84CRS


###### Niche overlap construction######CURRENT#####

library(ecospat)

bindrand<-read.csv("bindrand.csv") ###In working_runs folder
bind<-read.csv("spp_number.csv")
inv<-bindrand
nat<-bind

nat<-na.omit(nat)
inv<-na.omit(inv)

###summarize and filter by number of gbif occurrences >50
inv_select<-inv%>%
  group_by(spp_clean)%>%
  summarize(sum=n())%>%
  filter(sum>=50)

inv<-inv%>%
  filter(spp_clean%in%inv_select$spp_clean)

nat_select<-nat%>%
  group_by(spp_clean)%>%
  summarize(sum=n())%>%
  filter(sum>=50)

nat<-nat%>%
  filter(spp_clean%in%nat_select$spp_clean)
###
n_sp <- unique(nat$spp_clean) %>%
  length()

uni_nat<-unique(nat$spp_clean)
uni_inv<-unique(inv$spp_clean)
i
summary(nat)

sp_overlap<-NULL;

for (i in 1:n_sp) {
  inv_tmp<-inv %>% 
    filter(spp_clean==uni_inv[i])
  nat_tmp<-nat %>%
    filter(spp_clean==uni_nat[i])
  
  pca.env<-dudi.pca(rbind(nat_tmp, inv_tmp)[,4:7],scannf=F,nf=2)
  #pca.env<-na.omit(pca.env)
  ecospat.plot.contrib(contrib = pca.env$co, eigen=pca.env$eig)
  globclim<-pca.env$li
  #globclim<-globclim[-2,]
  scores.sp.nat<-suprow(pca.env,nat_tmp[which(nat_tmp[,8]==1),4:7])$li
  scores.clim.nat <- suprow(pca.env,nat_tmp[,4:7])$li
  scores.sp.inv<-suprow(pca.env,inv_tmp[which(inv_tmp[,8]==1),4:7])$li
  scores.clim.inv <- suprow(pca.env,inv_tmp[,4:7])$li
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


sp_overlap_dat<-as.data.frame(sp_overlap)
write.csv(sp_overlap_dat, "D_overlap_all_sp.csv")

nat_list<-as.data.frame(uni_nat)

write.csv(nat_list, "species_list.csv")

overlap_proc<-read.csv("D_overlap_all_sp.csv")

plot(overlap_proc$V1, overlap_proc$order)
plot(overlap_proc$order, overlap_proc$V1)
plot(overlap_proc$V1, overlap_proc$order)

library(ggplot2)

# Basic line plot
ggplot(data=overlap_proc, aes(x=order, y=V1))+
  geom_line()

#########rerun with climate projections######

futurestack<-getData("CMIP5",var="bio", res = 5, rcp=85, model ="AC", year=50)

rcp85<-subset(futurestack, c(1,5:6,12))
summary(rcp85)
writeRaster(rcp85[[1]], "temp.tif")
writeRaster(rcp85[[2]], "tmax.tif")
writeRaster(rcp85[[3]], "tmin.tif")
writeRaster(rcp85[[4]], "precip.tif")

dir() %>% stringr::str_subset("\\.tif$") -> future


library(epwshiftr)
 
rastlist<-list.files(path = ".", pattern = '.tif$', 
                     all.files = T, full.names = F)

futurestack<-lapply(rastlist, raster)
futurestack<-stack(futurestack)

library(rgdal)
summary(species_loc)
future_inv<-readOGR("Future_inv.shp")
sp_dat<-as.data.frame(future_inv)

values<-extract(futurestack, future_inv)
futurerandbind<-cbind.data.frame(future_inv, values)
futurerandbind
write.csv(futurerandbind, "future_inv_bind.csv")
bind<-read.csv("species_loc_bioclim_clean.csv")
##### ADD biolcim to random points
futurerand<-readOGR("Future_inv.shp")
futurerand<-read.csv("future_inv.csv")
randdat<-as.data.frame(futurerand)
randvalues<-extract(futurestack, futurerand)
bindrand2<-cbind.data.frame(futurerand, randvalues)
bindrand2
write.csv(bindrand2, "future_rand.csv")
bindrand<-read.csv("random_clean.csv")
tail(bindrand2)

###### Future niche overlap construction######
library(dplyr)
library(ecospat)

futurerand<-read.csv("future_inv_bind.csv")
futurebind<-read.csv("spp_number.csv")
finv<-futurerand
fnat<-futurebind
summary(fnat)
summary(finv)
nat<-na.omit(fnat)
inv<-na.omit(finv)
#write.csv(nat,"nat.csv")
# nat<-nat[complete.cases(nat),]
# inv<-inv[complete.cases(inv),]
# nat<-nat[!apply(is.na(nat) | nat == " ", 1, all),]
# inv<-inv[!apply(is.na(inv) | inv == " ", 1, all),]
# summary(nat)
#inv<-inv[-c(300174:301873)]
#write.csv(inv, "inv.csv")
# nat$species<-as.factor(nat$species)
# inv$species<-as.factor(inv$species)


finv_select<-finv%>%
  group_by(spp_clean)%>%
  summarize(sum=n())%>%
  filter(sum>=50)

finv<-finv%>%
  filter(spp_clean%in%finv_select$spp_clean)

fnat_select<-fnat%>%
  group_by(spp_clean)%>%
  summarize(sum=n())%>%
  filter(sum>=50)

fnat<-fnat%>%
  filter(spp_clean%in%fnat_select$spp_clean)

n_sp <- unique(fnat$spp_clean) %>%
  length()

funi_nat<-unique(fnat$spp_clean)
funi_inv<-unique(finv$spp_clean)
i
summary(fnat)

sp_overlap<-NULL;
pca_mapping<-NULL;

for (i in 1:n_sp) {
  finv_tmp<-finv %>% 
    filter(spp_clean==funi_inv[i])
  fnat_tmp<-fnat %>%
    filter(spp_clean==funi_nat[i])
  
  pca.env<-dudi.pca(rbind(fnat_tmp, finv_tmp)[,4:7],scannf=F,nf=2)
  pca_mapping<-rbind(pca.env, finv_tmp$x, finv_tmp$y)
  #pca.env<-na.omit(pca.env)
  ecospat.plot.contrib(contrib = pca.env$co, eigen=pca.env$eig)
  globclim<-pca.env$li
  #globclim<-globclim[-2,]
  scores.sp.nat<-suprow(pca.env,fnat_tmp[which(fnat_tmp[,8]==1),4:7])$li
  scores.clim.nat <- suprow(pca.env,fnat_tmp[,4:7])$li
  scores.sp.inv<-suprow(pca.env,finv_tmp[which(finv_tmp[,8]==1),4:7])$li
  scores.clim.inv <- suprow(pca.env,finv_tmp[,4:7])$li
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

grid.clim.nat
sp_overlap_future_dat<-as.data.frame(sp_overlap)
write.csv(sp_overlap_future_dat, "Future_D_overlap_all_sp.csv")
write.csv(sp_overlap_dat, "D_overlap_all_sp.csv")
pca_mapping_dat<-as.data.frame(pca_mapping)
pca_mapping_dat
future_d<-read.csv("Future_D_overlap_all_sp.csv")

"Current climate" <- overlap_proc
"Future climate" <- future_d

# Basic line plot

d_overlap_combo<-read.csv("d_overlap_combo_clean.csv")


p<-ggplot(data = d_overlap_combo, aes(x = order)) +
  geom_line(aes(y = current_d, colour = "Current Climate"), size =1.25) +
  geom_line(aes(y = future_d, colour = "Climate 2050"), size = 1.25) +
  scale_colour_manual("", 
                      breaks = c("Current Climate", "Climate 2050"),
                      values = c("Current Climate"="red", "Climate 2050"="blue")) +
  xlab("Species") +
  ylab("Niche overlap")+
  geom_ribbon(aes(ymin= current_d,
                  ymax= future_d),
                      fill = "grey")+
  theme(aspect.ratio = 1)
p

p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))

###Generate class taxonomy 
library(taxizedb)

ids<-name2taxid(c(sp2class), out_type = "summary")
??name2taxid
ids
test=classification(ids$id)
library(taxize)

test=classification(ids, db = 'itis')

test1=as.data.frame(t(sapply(names(test), function (x) test[[x]] [,1])[c(20,30),]))

test$`52216`


write.csv(test, "sp_taxonomy.csv")

sp2class<-d_overlap_combo$uni_nat
sp2class

class_delin<-tax_name(sci = c(sp2class), get = "class")
####spatial niche mapping 

#######Current conditions

library(dplyr)
library(ecospat)

range_summary<-read.csv("nat.csv")

range_summary_select<-range_summary%>%
  group_by(spp_clean)%>%
  summarize(sum=n())%>%
  filter(sum>=50)

range_summary<-range_summary%>%
  filter(spp_clean%in%range_summary_select$spp_clean)

n_sp <- unique(range_summary$spp_clean) %>%
  length()

uni_range_summary<-unique(range_summary$spp_clean)
uni_inv<-unique(inv$spp_clean)
i
summary(range_summary)




results <- data.frame(mn_tmp=NA, lwr_tmp=NA, upr_tmp=NA, mn_tmax=NA, lwr_tmax=NA, upr_tmax=NA,
                      mn_tmin=NA, lwr_tmin=NA, upr_tmin=NA, mn_p=NA, lwr_p=NA, upr_p=NA)

for (i in 1:n_sp) {
  # inv_tmp<-inv %>% 
  #   filter(spp_clean==uni_inv[i])
  nat_tmp<-range_summary %>%
    filter(spp_clean==uni_range_summary[i])
  ###tmp
  mn_tmp<-mean(nat_tmp$tmp)
  sder_tmp <- sd(nat_tmp$tmp)/sqrt(60)
  lwr_tmp <- mn_tmp - 1.96 * sder_tmp
  upr_tmp <- mn_tmp + 1.96 * sder_tmp
  ###tmax
  mn_tmax<-mean(nat_tmp$tmax)
  sder_tmax <- sd(nat_tmp$tmax)/sqrt(60)
  lwr_tmax <- mn_tmax - 1.96 * sder_tmax
  upr_tmax <- mn_tmax + 1.96 * sder_tmax
  ###tmax
  mn_tmin<-mean(nat_tmp$tmin)
  sder_tmin <- sd(nat_tmp$tmin)/sqrt(60)
  lwr_tmin <- mn_tmin - 1.96 * sder_tmin
  upr_tmin <- mn_tmin + 1.96 * sder_tmin
  ###precip
  mn_p<-mean(nat_tmp$p)
  sder_p <- sd(nat_tmp$p)/sqrt(60)
  lwr_p <- mn_p - 1.96 * sder_p
  upr_p <- mn_p + 1.96 * sder_p
  
  results[i, ] <- c(mn_tmp, lwr_tmp, upr_tmp, mn_tmax, lwr_tmax, upr_tmax,
                    mn_tmin, lwr_tmin, upr_tmin, mn_p, lwr_p, upr_p) 
  min = lwr_p
  max = upr_p
  
  stack2<-stack[[1]]
  stack2[stack2 <=min] <-NA
  stack2[stack2 >=max] <-NA
  
  stack2[stack2 >0] <- 1
  
  min = lwr_tmp
  max = upr_tmp
  stack5<-stack[[2]]
  stack5[stack5<=min]<-NA
  stack5[stack5>=max]<-NA
  
  stack5[stack5>0] <-1
  
  
  stack_merge<-mosaic(stack2, stack5, fun = sum)
  
  stack_merge <- stack_merge-1
  stack_merge[stack_merge < 1]<-NA
  plot(stack_merge)
  print(i)
  f <- paste0(i, '.tif')
  writeRaster(stack_merge, filename = f)
  
}

write.csv(results, "Climate_stdevs.csv")

####merge raster layers

dir() %>% str_subset("\\.tif$") -> stack_sp

stack_lt<-list(stack_sp)
stack_species<-stack(stack_sp)

sp_sum<-calc(stack, sum)

brick_sp<-brick(stack_species)

list.get<-lapply(stack_species, get)
mosaic<-do.call(merge, list.get)


writeRaster(stack_species, "Stack_species.tif")

sp_mosaic_global<-mosaic(stack_lt, fun = sum)

####Spatial RCP 8.5 predictions

library(dplyr)
library(ecospat)

range_summary<-read.csv("future_nat_bind.csv")

range_summary_select<-range_summary%>%
  group_by(spp_clean)%>%
  summarize(sum=n())%>%
  filter(sum>=50)

range_summary<-na.omit(range_summary)

range_summary<-range_summary%>%
  filter(spp_clean%in%range_summary_select$spp_clean)

n_sp <- unique(range_summary$spp_clean) %>%
  length()

uni_range_summary<-unique(range_summary$spp_clean)
uni_inv<-unique(inv$spp_clean)
i
summary(range_summary)


#RCP8.5 rasters

dir() %>% str_subset("\\.tif$") -> future_stack_sp

future_stack<-stack(future_stack_sp)

stack<-future_stack

summary(future_stack)
results_future <- data.frame(mn_tmp=NA, lwr_tmp=NA, upr_tmp=NA, mn_tmax=NA, lwr_tmax=NA, upr_tmax=NA,
                      mn_tmin=NA, lwr_tmin=NA, upr_tmin=NA, mn_p=NA, lwr_p=NA, upr_p=NA)

for (i in 1:n_sp) {
  # inv_tmp<-inv %>% 
  #   filter(spp_clean==uni_inv[i])
  nat_tmp<-range_summary %>%
    filter(spp_clean==uni_range_summary[i])
  ###tmp
  mn_tmp<-mean(nat_tmp$tmp)
  sder_tmp <- sd(nat_tmp$tmp)/sqrt(60)
  lwr_tmp <- mn_tmp - 1.96 * sder_tmp
  upr_tmp <- mn_tmp + 1.96 * sder_tmp
  ###tmax
  mn_tmax<-mean(nat_tmp$tmax)
  sder_tmax <- sd(nat_tmp$tmax)/sqrt(60)
  lwr_tmax <- mn_tmax - 1.96 * sder_tmax
  upr_tmax <- mn_tmax + 1.96 * sder_tmax
  ###tmax
  mn_tmin<-mean(nat_tmp$tmin)
  sder_tmin <- sd(nat_tmp$tmin)/sqrt(60)
  lwr_tmin <- mn_tmin - 1.96 * sder_tmin
  upr_tmin <- mn_tmin + 1.96 * sder_tmin
  ###precip
  mn_p<-mean(nat_tmp$p)
  sder_p <- sd(nat_tmp$p)/sqrt(60)
  lwr_p <- mn_p - 1.96 * sder_p
  upr_p <- mn_p + 1.96 * sder_p
  
  results_future[i, ] <- c(mn_tmp, lwr_tmp, upr_tmp, mn_tmax, lwr_tmax, upr_tmax,
                    mn_tmin, lwr_tmin, upr_tmin, mn_p, lwr_p, upr_p) 

  min = lwr_p
  max = upr_p
  
  stack2<-stack[[1]]
  stack2[stack2 <=min] <-NA
  stack2[stack2 >=max] <-NA
  
  stack2[stack2 >0] <- 1
  
  
  
  min = lwr_tmp
  max = upr_tmp
  stack5<-stack[[2]]
  stack5[stack5<=min]<-NA
  stack5[stack5>=max]<-NA
  
  stack5[stack5>0] <-1
  
  stack_merge<-mosaic(stack2, stack5, fun = sum)
  
  stack_merge <- stack_merge-1
  stack_merge[stack_merge < 1]<-NA
plot(stack_merge)
print(i)
f <- paste0(i, '.tif')
writeRaster(stack_merge, filename = f)

}

write.csv(results, "Climate_stdevs.csv")

####merge raster layers

dir() %>% str_subset("\\.tif$") -> stack_sp

stack_lt<-list(stack_sp)
stack_species<-stack(stack_sp)

sp_sum<-calc(stack, sum)

brick_sp<-brick(stack_species)

list.get<-lapply(stack_species, get)
mosaic<-do.call(merge, list.get)


writeRaster(stack_species, "Stack_species.tif")

sp_mosaic_global<-mosaic(stack_lt, fun = sum)

##### Create distribution map of potential invasives
##all pet trade sp

library(rgdal)
library(raster)
library(terra)

shp_files <- list.files(".", pattern = "\\.shp$")

total_nat<-readOGR("Total_native.shp")

r<-raster(ncol = 5000, nrow=5000)
s<-raster(ncol = 5000, nrow=5000)

r<-setValues(r, 1)
s<-setValues(s, 1)

extent(r)<-extent(total_inv)
extent(s)<-extent(total_inv)


for (i in 1:length(total_inv)) {
  x <- crop(r, total_nat[i,])
  x <- mask(x, total_nat[i,])
  print(i)
  s<-mosaic(s, x, fun = 'sum')
  plot(s)
  #f <- paste0(i, '.tif')
  #writeRaster(x, filename = f)

}

plot(s)
writeRaster(s, "IUCN_ranges_native.tif")

###just pot invasives

total_inv<-readOGR("Total_inv_non_native.shp")

r<-raster(ncol = 5000, nrow=5000)
w<-raster(ncol = 5000, nrow=5000)

r<-setValues(r, 1)
w<-setValues(s, 1)

extent(r)<-extent(total_inv)
extent(w)<-extent(total_inv)


for (i in 1:length(total_inv)) {
  x <- crop(r, total_inv[i,])
  x <- mask(x, total_inv[i,])
  print(i)
  w<-mosaic(w, x, fun = 'sum')
  plot(w)
  #f <- paste0(i, '.tif')
  #writeRaster(x, filename = f)
  
}

plot(w)
writeRaster(w, "IUCN_ranges_nn.tif")

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

library(dplyr)
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

##### Future propagule niche overlaps ####
library(stringr)
##### ADD biolcim to random points
futurerand<-readOGR("Future_inv.shp")
futurerand<-read.csv("Future_prop_xy.csv")
randdat<-as.data.frame(futurerand)

#RCP8.5 rasters

dir() %>% str_subset("\\.tif$") -> future_stack_sp

futurestack<-stack(future_stack_sp)
futurestack<-futurestack/10
futurestack
randvalues<-extract(futurestack, futurerand)
bindrand2<-cbind.data.frame(futurerand, randvalues)
summary(bindrand2)
write.csv(bindrand2, "Prop_future_rand1.csv")
bindrand<-read.csv("Prop_future_rand1.csv")
tail(bindrand)

###### Future niche overlap construction
library(dplyr)
library(ecospat)

#futurerand<-read.csv("inv.csv")
futurebind<-read.csv("Future_nat.csv")
inv<-bindrand
nat<-futurebind
summary(nat)
summary(inv)

inv$tmp<-as.numeric(inv$tmp)
inv$tmax<-as.numeric(inv$tmax)
inv$tmin<-as.numeric(inv$tmin)

nat<-na.omit(nat)
inv<-na.omit(inv)
#write.csv(nat,"nat.csv")
# nat<-nat[complete.cases(nat),]
# inv<-inv[complete.cases(inv),]
# nat<-nat[!apply(is.na(nat) | nat == " ", 1, all),]
# inv<-inv[!apply(is.na(inv) | inv == " ", 1, all),]
# summary(nat)
#inv<-inv[-c(300174:301873)]
#write.csv(inv, "inv.csv")
# nat$species<-as.factor(nat$species)
# inv$species<-as.factor(inv$species)


inv_select<-inv%>%
  group_by(spp_clean)%>%
  summarize(sum=n())%>%
  filter(sum>=50)

inv<-inv%>%
  filter(spp_clean%in%inv_select$spp_clean)

nat_select<-nat%>%
  group_by(spp_clean)%>%
  summarize(sum=n())%>%
  filter(sum>=50)

nat<-nat%>%
  filter(spp_clean%in%nat_select$spp_clean)

n_sp <- unique(nat$spp_clean) %>%
  length()

uni_nat<-unique(nat$spp_clean)
uni_inv<-unique(inv$spp_clean)
i
summary(nat)
summary(inv)
sp_overlap<-NULL;
pca_mapping<-NULL;

for (i in 1:n_sp) {
  inv_tmp<-inv %>% 
    filter(spp_clean==uni_inv[i])
  nat_tmp<-nat %>%
    filter(spp_clean==uni_nat[i])
  
  pca.env<-dudi.pca(rbind(nat_tmp, inv_tmp)[,4:7],scannf=F,nf=2)
  pca_mapping<-rbind(pca.env, inv_tmp$x, inv_tmp$y)
  #pca.env<-na.omit(pca.env)
  ecospat.plot.contrib(contrib = pca.env$co, eigen=pca.env$eig)
  globclim<-pca.env$li
  #globclim<-globclim[-2,]
  scores.sp.nat<-suprow(pca.env,nat_tmp[which(nat_tmp[,8]==1),4:7])$li
  scores.clim.nat <- suprow(pca.env,nat_tmp[,4:7])$li
  scores.sp.inv<-suprow(pca.env,inv_tmp[which(inv_tmp[,8]==1),4:7])$li
  scores.clim.inv <- suprow(pca.env,inv_tmp[,4:7])$li
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

grid.clim.nat
sp_overlap_future_prop<-as.data.frame(sp_overlap)
write.csv(sp_overlap_future_prop, "Future_D_overlap_all_sp_prop.csv")

pca_mapping_dat<-as.data.frame(pca_mapping)
pca_mapping_dat
overlap_proc_prop<-read.csv("Propagule_D_overlap_current.csv")
future_d_prop<-read.csv("Future_D_overlap_all_sp_prop.csv")
"Current climate" <- overlap_proc_prop
"Future climate" <- future_d_prop



# Basic line plot
library(ggplot2)
prop_overlap_combo<-read.csv("prop_d_combo.csv")###Propagule_run folder



pp<- ggplot(data = prop_overlap_combo, aes(x = order)) +
  geom_line(aes(y = current_d, colour = "Current Climate"), size =1.25) +
  geom_line(aes(y = future_d, colour = "Climate 2050"), size = 1.25) +
  scale_colour_manual("", 
                      breaks = c("Current Climate", "Climate 2050"),
                      values = c("Current Climate"="red", "Climate 2050"="blue")) +
  geom_line(data = d_overlap_combo, aes(x = order, y = current_d, colour = "Current Climate"), size =1.25) +
  geom_line(data = d_overlap_combo, aes(x = order, y = future_d, colour = "Climate 2050"), size = 1.25) +
  xlab("Species") +
  ylab("Niche overlap")+
  geom_ribbon(aes(ymin= current_d,
                  ymax= future_d),
              fill = "grey30")+
  geom_ribbon(data = d_overlap_combo, aes(ymin= current_d,
                  ymax= future_d),
              fill = "grey50")+
  theme(aspect.ratio = 1)
pp

pp<-pp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

