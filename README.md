
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Spatially Explicit Model of Landscape Exposure (SEMLE)

## Table of contents
1. [Overview](#Overview)
2. [Citing information](#Citing information) 
3. [Installation](#Installation) 
4. [Vignette](#Vingette)
5. [Performance](#Performance)

## **Overview** <a name="Overview"></a>

This repository includes the R code of the model used for predicting the exposure of bees to pesticides at landscape scales. The associated paper has been been published as an article in *Science of the Total Environment*, and is available through the publisher [here](XXX). In this repo you can find code that was used to produce the findings in the manuscript, as well as a vignette using the model in a simulated landscape.  



## **Citing information** <a name="Citing information"></a>

Lonsdorf, E. V., Nicholson, C. C., Rundlöf, M., & Williams, N. M. (2023). A spatially explicit model of landscape pesticide exposure to bees: development, exploration, and evaluation. *Science of the Total Environment* doi: TBD .

The model code has been released on Zenodo and has a DOI:  
[![DOI](https://zenodo.org/badge/184192167.svg)](https://zenodo.org/badge/latestdoi/184192167)  

## **Installation** <a name="Installation"></a>

```
install.packages("devtools")
library(devtools)
install_github("CCNicholson/SEMLE")
```

Or just copy and paste the code  

## **Vignette** <a name="Vignette"></a>

The model requires the following inputs:  

1. A raster of classified land use, where cell values correspond to land use types  
2. A raster of pesticide load, where cell values correspond to an observed or epxected quantity of pesticide (e.g., *g*/$m^2$)  
3. A look up table that translates landuse codes to habitat resource values (e.g., floral resources for a bee)  
4. The column name in this lookup table that contains the land use codes  
5. The column names in this lookup table that contains the habitat resource values  
6. A numerical value for $\gamma$ - the foraging range of the organism. This is used to specify the dimensions of internal moving window functions  
7. A description of organism foraging behavior as either  
* "simple" - where the organism visits patches based on distance only  
* "complex" - where the organism visits patches based on their distance and quality relative to all other patches. 


To demonstrate the use of the model we first simulate classified landscapes using the `landscapeR` and `raster`packages


```{r message=FALSE, warning=FALSE}
require(landscapeR)
require(raster)

### Landscape Params
set.seed(2024)
mapsize = 1000
patch_fract=8 # patch number param
patch_size=5000 # patch size param

### Create an empty raster
m = matrix(0, mapsize, mapsize) 
r = raster(m, xmn=0, xmx=1000, ymn=0, ymx=1000)

### Create landscape with three classified land uses
LU.ras <- makeClass(r, npatch = patch_fract, 
                    size = patch_size, 
                    val=1)
LU.ras <- makeClass(LU.ras, npatch = patch_fract, 
                    size = patch_size, 
                    val=2)
LU.ras <- makeClass(LU.ras, npatch = patch_fract, 
                    size = patch_size, 
                    val=3)
plot(LU.ras)

```

For the raster of pesticide load we will assume that one of the land use classes (e.g., class 1) receives an arbitrary value of pesticide application (0.5), while the other patches are untreated
```{r message=FALSE, warning=FALSE}

# create load look up table
load.df <- data.frame(LU_Code = c(0, 1, 2, 3), 
                      LU_Pval = c(0.0, 0.5, 0.0, 0.0))

# reclassify base landscape based on table
load.ras <- reclassify(LU.ras, load.df)

plot(load.ras, col = heat.colors(2, rev = T))

```

Internally, the model reclassifies the land use raster into a map of resources based on a user provided lookup table. It then uses this map to create a raster of distance weighted resources. Having this process internal to the model is self-contained but inefficient (discussed below). For the model to generate these resource maps, we create a data frame with a column of the 4 land use values (0, 1, 2, 3, 4) and a column of corresponding resource values (0.0, 0.9, 0.4, 0.2). Note, one could have multiple resource value columns (e.g., corresponding to different seasons) that one could iterate through.  

```{r message=FALSE, warning=FALSE}

# create resource look up table
RvalLookup.df <- data.frame(LU_Code = c(0, 1, 2, 3), 
                            LU_Rval = c(0.0, 0.9, 0.4, 0.2))  


```  

Now we are ready to run the model. We provide the two raster inputs (land use: `LU.ras` & pesticide load: `load.ras`). We provide the data frame (`RvalLookup.df`) with land use and resource columns specified (`LU_Code` & `LU_Rval`). We expect this organism to forage optimally based on patch distance *and quality*, so we select `"complex"` behavior. Lastly, we specify a foraging range. Note, this range value is numerical and should make sense given the map units of the coordinate reference system (e.g., EPSG:3857 & meters vs. EPSG:4326 & degrees). The model has a few control structures to make sure the user is providing appropriate spatial inputs, but it is advised to make sure the CRS and overlay of maps is correct. 

```{r message=FALSE, warning=FALSE, echo= F}
SEMLE <- function(LU_Map,
                  Load_Map,
                  LU_lookup, 
                  LU_ID_col, 
                  Rsrc_val_col, 
                  gamma_param, 
                  behavior = c("simple", "complex")){
  
  #### DEPENDENCIES ####
  require(smoothie)
  require(raster)
  
  #### CONTROL STRUCTURES ####

  if(!compareCRS(LU_Map, Load_Map)){
    stop("Raster coordinate reference systems do not match")
  }
  
  if(!all(res(LU_Map) == res(Load_Map))){
    stop("Raster resolutions do not match")
  }
  
  e <- intersect(extent(LU_Map), extent(Load_Map))
  if (is.null(e)) {
    stop("Raster extents do not overlap")
  }

  if (!is.numeric(LU_lookup[,LU_ID_col])) {
    stop("Land use codes are not numerical values")
  }
  
  if (!is.numeric(LU_lookup[,Rsrc_val_col])) {
    stop("Land use codes are not numerical values")
  }
  
  if(!(behavior %in% c("simple", "complex"))){
    return("behavior must be either 'simple' or 'complex'")
  }
  ### Make moving window object based on foraging range
  {maxforage.dist = 2*gamma_param
    c.size = res(LU_Map)[1]
    radius <- round(maxforage.dist/c.size) # matrix boundary for moving window
    
    # create matrices for weight arguments in focal() and setting moving window distance matrix
    weight.m = dist.m = matrix(data=1, 
                               nrow= radius * 2 + 1, 
                               ncol= radius * 2 + 1) 
    focal.c = median(1:nrow(dist.m)) # focal.c is defined as row number of focal cell in the distance matrix
    
    # calculating distance all cells from the focal cell in the nested matrix
    for(i in 1:nrow(dist.m)) {
      for (j in 1:ncol(dist.m)) {
        dist.m[i,j] = sqrt(((i-0.5)*c.size - (focal.c - 0.5)*c.size)^2 +
                             ((j-0.5)*c.size - (focal.c - 0.5)*c.size)^2 )
      }
    }
    
    dist.v = as.vector(dist.m)  # change distance matrix to vector
    effdist.m = exp(-dist.m / gamma_param) # calculate effective distance matrix (for FFT kernal smoothing)
    
    # set 0, when effdist > (2xforaging distance)
    effdist.m[which(dist.m > (2*gamma_param) | dist.m > maxforage.dist)] <- 0 # (for FFT kernal smoothing)
  }
  
  #### RESOURCES ####
  ### Reclassify landuse map with resource values based on lookup table
  {
    floral <- raster::reclassify(LU_Map, 
                         LU_lookup[,c(LU_ID_col, 
                                      Rsrc_val_col)])
    #plot(floral, main="Resources")
  }
  
  #### FORAGE ####
  ## Create rasters of distance weighted forage suitability based on floral resource maps
  {
    FFT_matrix <- smoothie::kernel2dsmooth(x=raster::as.matrix(floral),
                                 K=effdist.m, setup=T)
    forage_dw <- smoothie::kernel2dsmooth(x=raster::as.matrix(floral), 
                                          W=FFT_matrix)
    forage <- raster::raster(forage_dw, 
                     template=floral)
    
    #convert raster > 0 into all ones (will use to clip output raster and calc moving window normalization values)
    mask_land <- raster::crop(reclassify(floral, 
                                         cbind(0, 
                                               max(raster::values(floral), na.rm=T), 
                                               1),
                                         include.lowest =T), 
                      extent(floral))
    
    #calculate moving window weight sums across raster (will be the same value for most of the raster except along edges)
    mask_dw <- smoothie::kernel2dsmooth(x=as.matrix(mask_land), 
                              W=FFT_matrix )
    window_sum <- raster::raster(mask_dw, 
                         template=mask_land)
    
    #divide the distance weighted forage raster by the moving window sum & 
    #clip distance weighted raster to boundary of land use raster
    forage <- (forage/window_sum)* mask_land
    
    
    # plot(forage, main="Distance weighted resources")
    
  }
  
  #### EXPOSURE ####
  
  ### STEP 1: PREPARE PESTICIDE LOAD ###
  
  ### STEP 1a: PESTICIDE LOAD ###
  if (behavior == "simple") { load <- Load_Map }
  
  ### STEP 1b: FLORAL WEIGHTED PESTICIDE LOAD ###
  if (behavior == "complex") { load <- floral * Load_Map} 
  
  
  ### STEP 2: CALCULATE DISTANCE- (& FLORAL-) WEIGHTED PESTICIDE LOAD ###
  FFT_matrix <- smoothie::kernel2dsmooth(x=raster::as.matrix(load),
                               K=effdist.m,
                               setup=T)
  load_dw <- smoothie::kernel2dsmooth(x=raster::as.matrix(load),
                                      W=FFT_matrix)
  pestiLoad <- raster::raster(load_dw, template=load)
  #plot(pestiLoad)
  
  #convert land use raster > 0 into all ones (will use to clip output raster and calc moving window normalization values)
  mask_land <- raster::crop(raster::reclassify(forage, 
                               cbind(0, 
                                     max(raster::values(forage), na.rm=T), 
                                     1)),  
                    extent(forage))
  
  #calculate moving window weight sums across raster (will be the same value for most of the raster except along edges)
  mask_dw <- smoothie::kernel2dsmooth(x=as.matrix(mask_land),
                            W=FFT_matrix )
  window_sum <- raster::raster(mask_dw, template=mask_land)
  
  #divide the distance weighted insecticide raster by the moving window sum & 
  #clip distance weighted raster to boundary of land use raster
  pestiLoad <- (pestiLoad/window_sum)* mask_land
  
  ### STEP 3: CALCULATE EXPOSURE ###
  
  ### STEP 1a: DISTANCE WEIGHTED EXPOSURE ###
  if (behavior == "simple") {exposure <- pestiLoad}
  
  ### STEP 1b: DISTANCE- (& FLORAL-) WEIGHTED EXPOSURE ###
  if (behavior == "complex") {exposure <- pestiLoad/forage}
  
  # plot(exposure, col = heat.colors(n = 100, rev = T),main="Pesticide exposure")
  
  return(exposure)
}
```


```{r message=FALSE, warning=FALSE}

# Use the SEMLE function to create a exposure raster object
exposure.ras <- SEMLE(
  LU_Map = LU.ras, 
  Load_Map = load.ras, 
  LU_lookup = RvalLookup.df, 
  LU_ID_col = "LU_Code",
  Rsrc_val_col = "LU_Rval",
  behavior = "complex",
  gamma_param = 300)  

```  
  
The model outputs a raster of expected exposure. One could extract values from any point in this map, for example the location of a bee colony (depicted as a triangle).

```{r message=FALSE, warning=FALSE, echo= F}
require(pals)
require(ggplot2)
require(ggthemes)
require(tidyverse)
exposure.df <- as.data.frame(as(exposure.ras, "SpatialPixelsDataFrame")) %>% 
  rename("Exposure" = "layer" )
ggplot()+
  geom_tile(data=exposure.df, aes(x=x, y=y, fill = Exposure), alpha=0.8)+
  geom_point(aes(x = 500, 
                 y = 500),
             color = "white",
             fill = "black",
             alpha = 0.75,
             shape = 24,
              size = 3)+
  scale_fill_gradientn(colours=parula(100), 
                       guide = "colourbar")+
  scale_color_gradientn(colours=parula(100), 
                        guide = "colourbar")+
  coord_equal() +
  theme_map() +
  theme(legend.position="right")
```  
  
## **Performance** <a name="Performance"></a>  
  
The model is slow (this sample landscape ran in c. 3.5 seconds on my cpu), and run time will increase with map scale and foraging range. The model is also slow because it contains a few processes that could be run externally if one is iterating with the model (for example different seasons in the same place). The moving window object could be specified outside of and removed from the function. Similarly, in some instances it might make sense to create the resource maps and distance weighted forage maps outside of the function. These adjustments would speed things up. Finally, the model deploys moving window functions with a Fast Fourier Transform (thanks to [Melanie Kamerrer](https://github.com/melaniekamm) for the FFT code), which speeds things up quite a bit, but there may be even more efficient methods.
