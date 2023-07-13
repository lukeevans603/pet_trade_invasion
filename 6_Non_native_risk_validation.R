##############################
##### EVANS ET AL. 2022. #####
##############################

#### Non-native risk validation#####


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
spec_tax<-read.csv("Specific_taxonomy_non_native.csv")
sp_list<-list(spec_tax)
library(rgbif)
library(dplyr)
??rgbif
df_subset
devtools::install_github("ropensci/rgbif")

keys <- sapply(sp_list, function(x) name_suggest(x)$key[1], USE.NAMES=T)

mydata<-occ_data(taxonKey=keys, hasCoordinate=TRUE, country = "US")


occ_search(taxonKey = keys, limit = 5) 


specific<-read.csv("Specific_taxonomy_non_native.csv")

specific_species<-unique(spec_tax$species)
write.csv(specific_species, "spec_sp.csv")
remove_error<-specific[!grepl('sp.', specific$species),]

library(tidyverse)
library(spocc)
species<-sp_list
# Species to seearch for
species<-c(df1)
# empty data frame
df<- data.frame(x=NA, longitude=NA, latitude=NA)
df2<-as.list(sp_list)

# species to loop over 
sp_list<- unique(df1$species)
df2
#df1<-read.csv("spec_sp.csv")
df1<-remove_error
for (i in sp_list){
  out<-spocc::occ(query=i , from="gbif", has_coords = T, limit=500, gbifopts = list(country = "US"))
  print(i)
  out_df<-spocc::occ2df(out)
  write.csv(out_df, paste0(i,'.csv'))
  #out_df <- dplyr::select(out_df, out$gbif$data$name, out$gbif$data$longitude, out$gbif$data$latitude)
  #print(out_df)
  #complete<-dplyr::bind_rows(df, out_df) %>% na.omit()
}
??spocc::occ
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


###plot

centroid<-read.csv("centroid_current_climate_suitability_2.csv")


p<-ggplot(data = centroid, aes(x = order)) +
  geom_line(aes(y = Currentclim), size=1)+#, colour = "Current Climate"), size =1.25) #+
  # geom_line(aes(y = future_d, colour = "Climate 2050"), size = 1.25) +
  #scale_colour_manual("", 
  # breaks = c("Current Climate", "Climate 2050"),
  #values = c("Current Climate"="red", "Climate 2050"="blue")) +
  xlab("Species") +
  ylab("Current Climate Suitability")+
  #geom_hline(yintercept = 0)+
  #geom_ribbon(aes(ymin= current_d,
  #   ymax= future_d),
  # fill = "grey")+
  theme(aspect.ratio = 1)
p

q<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"))

q
