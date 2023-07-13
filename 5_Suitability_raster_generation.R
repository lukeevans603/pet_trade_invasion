##############################
##### EVANS ET AL. 2022. #####
##############################

####bivariate mapping

# Install:

install.packages("classInt")
install.packages("raster")
install.packages("rgdal")
install.packages("dismo")
install.packages("XML")
install.packages("maps")
install.packages("sp")
install.packages("basemaps")
install.packages("RgoogleMaps")
# Load:

library(classInt)
library(raster)
library(terra)
library(rgdal)
library(dismo)
library(XML)
library(maps)
library(sp)
library(dplyr)
library(stringr)
library(RgoogleMaps)
library(ggmap)
# lat<-c(22,50)
# lon<-c(-65,-127)
# center = c(mean(lat), mean(lon))
# zoom<-3
# terrmap<-GetMap(center=center, zoom = zoom, maptype = "terrain", destfile = "base.png")

us_map<-get_map(location = "USA", zoom = "auto", maptype = "terrain-background",
                color = "bw")
  
  

# The function that produces the colour matrix
colmat <- function(nbreaks = 3, breakstyle = "quantile",
                   upperleft = "#0096EB", upperright = "#820050", 
                   bottomleft = "#BEBEBE", bottomright = "#FFE60F",
                   xlab = "x label", ylab = "y label", plotLeg = TRUE,
                   saveLeg = FALSE) {
  # TODO - replace any tidyr, dplyr etc. functions with data.table #
  library(tidyverse)
  require(ggplot2)
  require(classInt)
  if (breakstyle == "sd") {
    warning("SD breaks style cannot be used.\nWill not always return the correct number of breaks.\nSee classInt::classIntervals() for details.\nResetting to quantile",
            call. = FALSE, immediate. = FALSE)
    breakstyle <- "quantile"}

  my.data <- seq(0, 1, .01)
  my.class <- classInt::classIntervals(my.data,
                                       n = nbreaks,
                                       style = breakstyle,
  )
  my.pal.1 <- classInt::findColours(my.class, c(upperleft, bottomleft))
  my.pal.2 <- classInt::findColours(my.class, c(upperright, bottomright))
  col.matrix <- matrix(nrow = 101, ncol = 101, NA)
  for (i in 1:101) {
    my.col <- c(paste(my.pal.1[i]), paste(my.pal.2[i]))
    col.matrix[102 - i, ] <- classInt::findColours(my.class, my.col)
  }
  col.matrix.plot <- col.matrix %>%
    as.data.frame(.) %>% 
    mutate("Y" = row_number()) %>%
    mutate_at(.tbl = ., .vars = vars(starts_with("V")), .funs = list(as.character)) %>% 
    pivot_longer(data = ., cols = -Y, names_to = "X", values_to = "HEXCode") %>% 
    mutate("X" = as.integer(sub("V", "", .$X))) %>%
    distinct(as.factor(HEXCode), .keep_all = TRUE) %>%
    mutate(Y = rev(.$Y)) %>% 
    dplyr::select(-c(4)) %>%
    mutate("Y" = rep(seq(from = 1, to = nbreaks, by = 1), each = nbreaks),
           "X" = rep(seq(from = 1, to = nbreaks, by = 1), times = nbreaks)) %>%
    mutate("UID" = row_number())
  # Use plotLeg if you want a preview of the legend
  if (plotLeg) {
    p <- ggplot(col.matrix.plot, aes(X, Y, fill = HEXCode)) +
      geom_tile() +
      scale_fill_identity() +
      coord_equal(expand = FALSE) +
      theme_void() +
      theme(aspect.ratio = 1,
            axis.title = element_text(size = 12, colour = "black",hjust = 0.5, 
                                      vjust = 1),
            axis.title.y = element_text(angle = 90, hjust = 0.5)) +
      xlab(bquote(.(xlab) ~  symbol("\256"))) +
      ylab(bquote(.(ylab) ~  symbol("\256")))
    print(p)
    assign(
      x = "BivLegend",
      value = p,
      pos = .GlobalEnv
    )
  }
  # Use saveLeg if you want to save a copy of the legend
  if (saveLeg) {
    ggsave(filename = "bivLegend.pdf", plot = p, device = "pdf",
           path = "./", width = 4, height = 4, units = "in",
           dpi = 300)
  }
  seqs <- seq(0, 100, (100 / nbreaks))
  seqs[1] <- 1
  col.matrix <- col.matrix[c(seqs), c(seqs)]
  attr(col.matrix, "breakstyle") <- breakstyle
  attr(col.matrix, "nbreaks") <- nbreaks
  return(col.matrix)
}

# Function to assign colour-codes to raster data
# As before, by default assign tercile breaks
bivariate.map <- function(rasterx, rastery, colourmatrix = col.matrix,
                          export.colour.matrix = TRUE,
                          outname = paste0("colMatrix_rasValues", names(rasterx))) {
  # TO DO - replace raster with terra #
  require(raster)
  require(classInt)
  # export.colour.matrix will export a data.frame of rastervalues and RGB codes 
  # to the global environment outname defines the name of the data.frame
  quanx <- getValues(rasterx)
  tempx <- data.frame(quanx, quantile = rep(NA, length(quanx)))
  brks <- with(tempx, classIntervals(quanx,
                                     n = attr(colourmatrix, "nbreaks"),
                                     style = attr(colourmatrix, "breakstyle"))$brks)
  ## Add (very) small amount of noise to all but the first break
  ## https://stackoverflow.com/a/19846365/1710632
  brks[-1] <- brks[-1] + seq_along(brks[-1]) * .Machine$double.eps
  r1 <- within(tempx, quantile <- cut(quanx,
                                      breaks = brks,
                                      labels = 2:length(brks),
                                      include.lowest = TRUE))
  quantr <- data.frame(r1[, 2])
  quany <- getValues(rastery)
  tempy <- data.frame(quany, quantile = rep(NA, length(quany)))
  brksy <- with(tempy, classIntervals(quany,
                                      n = attr(colourmatrix, "nbreaks"),
                                      style = attr(colourmatrix, "breakstyle"))$brks)
  brksy[-1] <- brksy[-1] + seq_along(brksy[-1]) * .Machine$double.eps
  r2 <- within(tempy, quantile <- cut(quany,
                                      breaks = brksy,
                                      labels = 2:length(brksy),
                                      include.lowest = TRUE
  ))
  quantr2 <- data.frame(r2[, 2])
  as.numeric.factor <- function(x) {
    as.numeric(levels(x))[x]
  }
  col.matrix2 <- colourmatrix
  cn <- unique(colourmatrix)
  for (i in 1:length(col.matrix2)) {
    ifelse(is.na(col.matrix2[i]),
           col.matrix2[i] <- 1, col.matrix2[i] <- which(
             col.matrix2[i] == cn
           )[1]
    )
  }
  # Export the colour.matrix to data.frame() in the global env
  # Can then save with write.table() and use in ArcMap/QGIS
  # Need to save the output raster as integer data-type
  if (export.colour.matrix) {
    # create a dataframe of colours corresponding to raster values
    exportCols <- as.data.frame(cbind(
      as.vector(col.matrix2), as.vector(colourmatrix),
      t(col2rgb(as.vector(colourmatrix)))
    ))
    # rename columns of data.frame()
    colnames(exportCols)[1:2] <- c("rasValue", "HEX")
    # Export to the global environment
    assign(
      x = outname,
      value = exportCols,
      pos = .GlobalEnv
    )
  }
  cols <- numeric(length(quantr[, 1]))
  for (i in 1:length(quantr[, 1])) {
    a <- as.numeric.factor(quantr[i, 1])
    b <- as.numeric.factor(quantr2[i, 1])
    cols[i] <- as.numeric(col.matrix2[b, a])
  }
  r <- rasterx
  r[1:length(r)] <- cols
  return(r)
}

dir() %>% str_subset("\\.tif$") -> ras

stack<-stack(ras)
stack

bivmapQ <- bivariate.map(rasterx = stack[["current_climate_rescaled"]], rastery = r[["Propagule_pressure"]],
                         export.colour.matrix = FALSE,
                         colourmatrix = col.matrixQ)


#### TESTS ####
example <- FALSE
if (example) {
  library(data.table)
  library(tidyverse)
  library(raster)
  library(classInt)
  library(patchwork)
  
  # Retrieve some raster data from worldclim
  # and clip to SE Asia and Northern Australia
  clipExt <- extent(stack)
  # r <- getData("worldclim", var = "bio", res = 10)
  # r <- crop(r[[c(1, 12)]], clipExt)
  # r[[1]] <- r[[1]]/10
  
  names(stack) <- c("Climate","Propagules")
  
  # Define the number of breaks
  nBreaks <- 5
  
  # Create the colour matrix
  col.matrixQ <- colmat(nbreaks = nBreaks, breakstyle = "quantile",
                        xlab = "Climate Suitability", ylab = "Propagule Pressure", 
                        bottomright = "#F7900A", upperright = "#993A65",
                        bottomleft = "#44B360", upperleft = "#3A88B5",
                        saveLeg = FALSE, plotLeg = TRUE)
  
  # create the bivariate raster
  bivmapQ <- bivariate.map(rasterx = stack[["Climate"]], rastery = stack[["Propagules"]],
                           export.colour.matrix = FALSE,
                           colourmatrix = col.matrixQ)
  
  # Convert to dataframe for plotting with ggplot
  bivMapDFQ <- setDT(as.data.frame(bivmapQ, xy = TRUE))
  colnames(bivMapDFQ)[3] <- "BivValue"
  bivMapDFQ <- melt(bivMapDFQ, id.vars = c("x", "y"),
                    measure.vars = "BivValue",
                    value.name = "bivVal",
                    variable.name = "Variable")
  
  # Make the map using ggplot
  map_q <- ggplot(bivMapDFQ, aes(x = x, y = y)) +
    geom_raster(aes(fill = bivVal)) +
    scale_y_continuous(breaks = seq(-20, 60, by = 10), 
                       labels = paste0(seq(-20, 60, 10), "?")) +
    scale_x_continuous(breaks = seq(50,175,25), 
                       labels = paste0(seq(50,175,25), "?")) +
    scale_fill_gradientn(colours = col.matrixQ, na.value = "white") + 
    theme_bw() +
    theme(text = element_text(size = 10, colour = "black")) +
    #borders(colour = "black", size = 0.5) +
    coord_quickmap(expand = FALSE, xlim = clipExt[1:2], ylim = clipExt[3:4])# +
    #theme(legend.position = "none",
     #     plot.background = element_blank(),
      #    strip.text = element_text(size = 12, colour = "black"),
       #   axis.text.y = element_text(angle = 90, hjust = 0.5),
        #  axis.text = element_text(size = 12, colour = "black"),
         # axis.title = element_text(size = 12, colour = "black")) +
    #labs(x = "Longitude", y = "Latitude")
  map_q
  # Using different breaks algorithm
  # Create the colour matrix
  col.matrixF <- colmat(nbreaks = nBreaks, breakstyle = "fisher",
                        xlab = "Climate Suitability", ylab = "Propagule Likelihood", 
                        bottomright = "#F7900A", upperright = "#993A65",
                        bottomleft = "#44B360", upperleft = "#3A88B5",
                        saveLeg = FALSE, plotLeg = TRUE)
  
  # create the bivariate raster
  bivmapF <- bivariate.map(rasterx = stack[["Climate"]], rastery = stack[["Propagules"]],
                           export.colour.matrix = FALSE,
                           colourmatrix = col.matrixF)
  bivMapDFF <- setDT(as.data.frame(bivmapF, xy = TRUE))
  colnames(bivMapDFF)[3] <- "BivValue"
  bivMapDFF <- melt(bivMapDFF, id.vars = c("x", "y"),
                    measure.vars = "BivValue",
                    value.name = "bivVal",
                    variable.name = "Variable")
  
  # Make the map using ggplot
  map_F <- ggplot(bivMapDFF, aes(x = x, y = y)) +
    geom_raster(aes(fill = bivVal)) +
    scale_y_continuous(breaks = seq(-20, 60, by = 10), 
                       labels = paste0(seq(-20, 60, 10), "?")) +
    scale_x_continuous(breaks = seq(50,175,25), 
                       labels = paste0(seq(50,175,25), "?")) +
    scale_fill_gradientn(colours = col.matrixF, na.value = "transparent") + 
    theme_bw() +
    $theme(text = element_text(size = 10, colour = "black")) +    
    #borders(colour = "black", size = 0.5) +
    #coord_quickmap(expand = FALSE, xlim = clipExt[1:2], ylim = clipExt[3:4]) +
    # theme(legend.position = "none",
    #       plot.background = element_blank(),
    #       strip.text = element_text(size = 12, colour = "black"),
    #       axis.text.y = element_text(angle = 90, hjust = 0.5),
    #       axis.text = element_text(size = 12, colour = "black"),
    #       axis.title = element_text(size = 12, colour = "black")) +
    # labs(x = "Longitude", y = "Latitude")
    # 
  fig <- {{map_q + ggtitle("Quantile breaks")} / {map_F + ggtitle("Fisher breaks")} } + 
    inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                   colour = NA)), 
                  left = 0.5, bottom = 0.5, right = 1.9, top = 1.20,
                  align_to = "full") +
    plot_annotation(caption = "Both plots are made on the same data, but have breaks defined differently")
  
  map_F
  
  # Save
  ggsave(plot = fig,
         filename = "BivariatePlot_ggsave.pdf",
         device = "pdf", path = "./",
         width = 6, height = 7, units = "in",
         dpi = 320)
}

map_q
  fig <- {{map_q}} #/ {map_F + ggtitle("Fisher breaks")} } + 
    inset_element(BivLegend + theme(plot.background = element_rect(fill = "white",
                                                                   colour = NA)), 
                  left = 0.5, bottom = 0.5, right = 1.9, top = 1.20,
                  align_to = "full") +
      #plot_annotation(caption = "Both plots are made on the same data, but have breaks defined differently")
    fig
    
    
pacman::p_load("tidyverse",
                   "ggthemes",
                   "scales",
                   "raster",
                   "rgeos",
                   "viridis",
                   "tigris")
    
states <- tigris::states(class = "sp") %>%
      subset(!NAME %in% c("Puerto Rico",
                          "American Samoa",
                          "Guam",
                          "Puerto Rico",
                          "United States Virgin Islands",
                          "Commonwealth of the Northern Mariana Islands",
                          "Hawaii",
                          "Alaska"))
    
gdem <- raster("demera5.tif") %>%
      raster::crop(states) %>%
      raster::mask(states)
    

base<-ggplot() + 
    geom_raster(data = as.data.frame(gdem, xy = T), aes(x = x, y = y, fill = demera5)) +
      # geom_polygon(data = states, aes(x = long, y = lat, group = group), fill = NA, col = "grey40", size = 1) +
    coord_cartesian(xlim = c(-125, -60), ylim = c(20, 50)) + 
    scale_fill_gradientn("Elevation (m)", colors = rev(grey.colors(256)), na.value = "white") + 
    theme_map() + theme(legend.position = "bottom",
                          legend.key=element_blank())
version
