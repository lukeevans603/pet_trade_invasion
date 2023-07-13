##############################
##### EVANS ET AL. 2022. #####
##############################

####spatial niche mapping 
library(dplyr)
library(ecospat)
library(raster)

range_summary<-read.csv("nat.csv")

range_summary_select<-range_summary%>%
  group_by(spp_clean)%>%
  summarize(sum=n())%>%
  filter(sum>=5)

range_summary<-range_summary%>%
  filter(spp_clean%in%range_summary_select$spp_clean)

n_sp <- unique(range_summary$spp_clean) %>%
  length()

uni_range_summary<-unique(range_summary$spp_clean)
uni_inv<-unique(inv$spp_clean)
i
summary(range_summary)



results <- data.frame(mean=NA, lwr=NA, upr=NA)   # Not "" which makes the variables character strings
set.seed(1234)
for (i in 1:samples) {
  x <- rgamma(n, shape = s, rate = r)
  mn <- mean(x)
  sder <- sd(x)/sqrt(n)
  lwr <- mn - 1.96 * sder
  upr <- mn + 1.96 * sder
  results[i, ] <- c(mn, lwr, upr) 
}

??rgamma

sp_overlap<-NA;
pca_mapping<-NA;

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
  min = lwr_tmax
  max = upr_tmax
  
  stack2<-stack[[1]]
  stack2[stack2 <=min] <-NA
  stack2[stack2 >=max] <-NA
  
  plot(stack2)
  
  min = lwr_tmp
  max = upr_tmp
  
  stack5<-stack[[4]]
  stack5[stack5<=min]<-NA
  stack5[stack5>=max]<-NA
  
  plot(stack5)
  print(i)
  stack_merge<-mosaic(stack2, stack5, fun = mean)
  f <- paste0(i, '.tif')
  writeRaster(stack_merge, filename = f)
  
}

write.csv(results, "Climate_stdevs.csv")
sp_plot_list<-as.data.frame(uni_range_summary)
write.csv(sp_plot_list, "sp_plot_list.csv")
###raster_rescaling


library(raster)
library(sf)
dir() %>% str_subset("\\.tif$") -> processing

stack<-stack(processing)

stack_merge_proc<-mosaic(stack, fun = mean)

plot(stack)

library(terra)
img <- list.files(getwd(), "tif$", full.names=TRUE)  
ic <- src(lapply(img, rast))
r <- mosaic(ic)

plot(r)
writeRaster(r, "merge_test.tif")

rescale <- function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
  if(is.null(x.min)) x.min = min(x)
  if(is.null(x.max)) x.max = max(x)
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

rescale(200, x.min = 0, x.max = 200, new.min = 0, new.max = 1)
rescale(100, x.min = 0, x.max = 200, new.min = 0, new.max = 1)
rescale(0, x.min = 0, x.max = 200, new.min = 0, new.max = 1)






#####raster stacks

library(raster)
library(sf)
library(stringr)

dir() %>% str_subset("\\.tif$") -> ras

stack<-stack(ras)
stack
plot(stack[[2]])

stack[[1]]<-stack[[1]]/10

stack[[2]]<-stack[[2]]/10
stack[[3]]<-stack[[3]]/10
stack[[4]]<-stack[[4]]/10



min = lwr_p
max = upr_p

stack2<-stack[[1]]
stack2[stack2 <=min] <-NA
stack2[stack2 >=max] <-NA

stack2[stack2 >0] <- 1

plot(stack2)

min = lwr_tmp
max = upr_tmp
plot(stack[[2]])
stack5<-stack[[2]]
stack5[stack5<=min]<-NA
stack5[stack5>=max]<-NA

stack5[stack5>0] <-1

plot(stack5)
stack_merge<-mosaic(stack2, stack5, fun = sum)

stack_merge <- stack_merge-1
stack_merge[stack_merge < 1]<-NA


plot(stack_merge)

writeRaster(stack_merge, "species_range_trial.tif", overwrite =T)
