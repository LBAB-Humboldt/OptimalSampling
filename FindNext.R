FindNext<-function(occ.table, env.vars, add.site=NULL){
  #Finds the next most complementary or disimilar site to sampled sites
  #within a set of target sites.
  
  #Args:
  #   occ.table:  data frame with species occurrence information. 
  #               Must contain lon,lat,species and site columns.
  #   env.vars:   RasterStack or RasterBrick of environmental 
  #               variables that define gradients of interest
  #   add.site:   data frame with x and y coordinates of sites
  #               that should be added to the sampled site set
  #               for iterative selection of survey sites. See
  #               example.
  #
  #Returns:
  #   A list object containing:
  #     dist.raster:  Raster of ED complementarity values
  #     dist.table:   data frame of x,y and ED complementarity
  #                   values for all target sites.
  #Example:
  #   Pending development.
  
  #Compute dissimilarity model from occurrence data
  gdm.table <- formatsitepair(occ.table, bioFormat=2, XColumn="lon", YColumn="lat", 
                              sppColumn="species", siteColumn="site", 
                              predData=env.vars)
  gdm.table <- na.omit(gdm.table)
  gdm.rast <- gdm(gdm.table, geo=T)
  
  #Create a cell id raster
  t.raster <- env.vars[[1]]
  t.raster[1:ncell(t.raster)]<-1:ncell(t.raster)
  t.raster[is.na(prod(env.vars))]<-NA
  
  #Prepare table for gdm prediction in unsampled sites
  sampled <- na.omit(unique(cellFromXY(t.raster, occ.table[,c("lon","lat")])))
  if(!is.null(add.site)){
    cell2add <- cellFromXY(t.raster,add.site)
    sampled <- c(sampled, cell2add)
  }
  t.raster[sampled] <- NA
  usmpled <- na.omit(unique(getValues(t.raster)))
  site.grid <- expand.grid(sampled, usmpled)
  gdm.pred.table <- cbind(distance=0, weights=1, xyFromCell(t.raster, site.grid$Var1), 
                          xyFromCell(t.raster, site.grid$Var2), env.vars[site.grid$Var1],
                          env.vars[site.grid$Var2])
  colnames(gdm.pred.table)<-colnames(gdm.table)
  
  #Compute GDM to new table
  gdm.pred <- predict(gdm.rast, gdm.pred.table)
  res <- data.frame(c.sampled = site.grid$Var1, c.unsampl = site.grid$Var2, dis=gdm.pred)
  res2 <- dcast(res, c.unsampl~c.sampled, value.var="dis")
  min.dis.vals <- apply(res2[, 2:ncol(res2)],1,min)  #For all candidate sites, 
  #find the minimum distance
  #to a surveyed site
  
  #Plot distances from sampled to unsampled sites
  dist2sampled <- t.raster
  dist2sampled[t.raster] <- 0
  dist2sampled[res2$c.unsampl] <- min.dis.vals
  plot(dist2sampled,col=rev(heat.colors(10)))
  
  #Order sites in decreasing order of predicted biological distance
  #First site on list is the most disimilar site to surveyed sites
  dists <- data.frame(xyFromCell(t.raster,res2$c.unsampl), dist=min.dis.vals)
  dists <- dists[order(-dists$dist),]
  points(dists[1,1:2],col="blue",pch=18)
  points(occ.table[,c("lon","lat")],pch=18,cex=0.6)
  return(list(dist.raster=dist2sampled, dist.table=dists))
}