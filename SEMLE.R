
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
