#Toy implementation with Boyacá data
library(raster)
library(rgdal)
library(gdm)
sp.data<-read.csv("C:/Modelos/BirdRichness/Modelos/summary/ebdcol_v2.csv",as.is=T)
env.data<-stack(paste0("Z:/Capas/NewExtent/","bio_",c(1,2,3,4,12,15,18),".asc"))
col<-readOGR("D:/Datos/BaseLayers","COL_adm1")

##Perform analysis on Boyacá data only

#Filter to Boyaca
boy<-col[col@data$NAME_1=="BoyacÃ¡",]
env.boy<-crop(env.data,boy)
env.boy<-mask(env.boy,boy)
sp.boy<-sp.data[which(!is.na(extract(env.boy,sp.data[,c(3,2)]))),]
template<-env.boy[[1]]
template[1:ncell(template)]<-1:ncell(template)
template[is.na(env.boy[[1]])]<-NA
sp.boy$inCol<-extract(template,sp.boy[,c(3,2)])
sp.boy<-na.omit(sp.boy)
sp.env<-extract(env.boy,sp.boy[,c(3,2)])
sp.env<-cbind(sp.boy[,2:4],sp.env)

#occ.table: lon, lat, species and site data frame
occ.table<-sp.boy
colnames(occ.table)<-c("species","lat","lon","site")
rm(list=setdiff(ls(), c("occ.table","env.boy")))

#Iteratively find best locations
while(response!="S"){
  if(exists(new.sites)){
    iter<-FindNext(occ.table, env.vars, add.site=new.sites)
  } else {
    new.sites<-data.frame()
    iter<-FindNext(occ.table, env.vars, add.site=NULL)
  }
  title(paste("ED = ",round(cellStats(iter$dist.raster,"sum"))))
  response <- readline(prompt="Choose option: (I)gnore location, (A)dd this location, Add (C)ustom location, (S)top ")

  if(response=="I"){
    #Code to remove site from target areas
  }
  if(response=="A"){
    new.sites <- rbind(new.sites, iter$dist.table[1, 1:2])
  }
  if(response=="C"){
    coords <- readline(prompt="Enter comma separated coordinates(lon, lat): ")
    coords <- as.numeric(strsplit(coords,",")[[1]])
    new.sites <- rbind(new.sites, coords) 
  }
  #If not any of the above set response to S
}
