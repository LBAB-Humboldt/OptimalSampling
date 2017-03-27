FindSimilar<-function(site.coords, threshold, occ.table, env.vars){
  require(gdm)
  require(reshape2)
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
  sampled <- na.omit(unique(cellFromXY(t.raster, site.coords[,c("lon","lat")])))
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
  if(ncol(res2)>2){
    min.dis.vals <- apply(res2[, 2:ncol(res2)],1,min)  #For all candidate sites, 
  } else {
    min.dis.vals <- res2[, 2]
  }
  #find the minimum distance
  #to a surveyed site
  
  #Plot distances from sampled to unsampled sites
  dist2sampled <- t.raster
  dist2sampled[res2$c.unsampl] <- (1 - min.dis.vals)
  plot(dist2sampled)
  threshold<-quantile(na.omit(getValues(dist2sampled)),threshold)
  plot(dist2sampled > threshold, 
       col=c(rgb(255, 255, 255, 0, maxColorValue=255),"cyan"), add=TRUE)
  points(site.coords,col="blue",pch=18)
  return(list(dist.raster=dist2sampled, dist.table=dists))
}