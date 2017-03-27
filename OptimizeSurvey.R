#Iteratively find best locations
OptimizeSurvey<-function(occ.table, env.vars){
  #Interactive tool to find the optimum locations to fill gaps in sampling.
  
  #Args:
  #   occ.table:  data frame with species occurrence information. 
  #               Must contain lon,lat,species and site columns.
  #   env.vars:   RasterStack or RasterBrick of environmental 
  #               variables that define gradients of interest
  #
  #Returns:
  #   A list object containing:
  #     ed.rasters:  ED complementarity surfaces for each iteration
  #     sites:       SpatialPointsData frame with coordinates, ED total
  #                  complementarity and remarks for each site selected.
  #Example:
  #   result <- OptimizeSurvey(occ.table, env.boy)
  
  #Initialize parameters
  response <- "Y"
  new.sites<-data.frame()
  par(mfrow=c(2,1))
  ed.rasters<-stack()
  iter <- 1
  
  #Start interactive loop
  while(response=="Y"){
    if(nrow(new.sites)!=0){
      iter<-FindNext(occ.table, env.vars, add.site=new.sites[,1:2])
    } else {
      iter<-FindNext(occ.table, env.vars, add.site=NULL)
    }
    ed.total<- round(cellStats(iter$dist.raster,"sum"))
    response <- readline(prompt="Choose option: (I)gnore location, (A)dd this location, Add (C)ustom location\n")

    #Case 1: ignore site
    if(response=="I"){
      ignore.cell<-cellFromXY(iter$dist.raster, iter$dist.table[1,1:2])
      env.vars[ignore.cell]<-NA
    }
    
    #Case 2: add site to survey locations
    if(response=="A"){
      add.row<-data.frame(iter$dist.table[1, 1:2], ed.total,"Choosen by software")
      colnames(add.row)<-c("x","y","EDTotal","Remarks")
      new.sites <- rbind(new.sites, add.row)
      ed.rasters<-addLayer(ed.rasters, iter$dist.raster)
      plot(iter$dist.raster, col=rev(heat.colors(10)))
      points(iter$dist.table[1, 1:2],col="blue",pch=18)
      points(occ.table[,c("lon","lat")],pch=18,cex=0.6)
      plot(1:nrow(new.sites),new.sites[,3],type="l",ylab="EDtotal",xlab="sitios")
      points(1:nrow(new.sites),new.sites[,3])
      if(nrow(new.sites) > 1){
        print(new.sites)
        sites<-nrow(new.sites)
        change <- round((new.sites$EDTotal[sites]-new.sites$EDTotal[sites-1])*100/new.sites$EDTotal[sites], 5)
        legend("topright",legend=change, bty="n",pch=NA)
      }
    }
    
    #Case 3: choose custom survey location
    if(response=="C"){
      coords <- readline(prompt="Enter comma separated coordinates(lon, lat): ")
      coords <- as.numeric(strsplit(coords,",")[[1]])
      add.row <- data.frame(x=coords[1], y=coords[2], EDTotal=ed.total, Remarks="Choosen by user")
      colnames(add.row)<-c("x","y","EDTotal","Remarks")
      new.sites <- rbind(new.sites, add.row)
      ed.rasters<-addLayer(ed.rasters, iter$dist.raster)
      plot(iter$dist.raster, col=rev(heat.colors(10)))
      points(coords[1],coords[2],col="blue",pch=18)
      points(occ.table[,c("lon","lat")],pch=18,cex=0.6)
      plot(1:nrow(new.sites),new.sites[,3],type="l",ylab="EDtotal",xlab="sitios")
      points(1:nrow(new.sites),new.sites[,3])
      if(nrow(new.sites) > 1){
        print(new.sites)
        sites<-nrow(new.sites)
        change <- round((new.sites$EDTotal[sites]-new.sites$EDTotal[sites-1])*100/new.sites$EDTotal[sites], 5)
        legend("topright",legend=change, bty="n",pch=NA)
      }
    }
    
    #End loop?
    response<-readline(prompt="Find next optimum site: (Y)es, (N)o: ")
  }
  coordinates(new.sites)<- ~x+y
  return(list(sites=new.sites, ed.rasters=ed.rasters))
}
