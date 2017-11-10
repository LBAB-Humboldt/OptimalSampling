#Estimate richness from ebird data
load("M:/BirdRichness/rawEBD/ebd.all.RData")
library(lubridate)
library(raster)
library(reshape2)
library(vegan)
ebd.all$OBSERVATION.DATE<-ymd(ebd.all$OBSERVATION.DATE)
ebd.all$year<-year(ebd.all$OBSERVATION.DATE)
ebd.all$LATITUDE<-as.numeric(ebd.all$LATITUDE)
ebd.all$LONGITUDE<-as.numeric(ebd.all$LONGITUDE)
ebd.sub<-ebd.all[which(ebd.all$CATEGORY=="species"),]
cellid<-raster("C:/Modelos/BirdRichness/Modelos/summary/template2.5m.tif")
ebd.sub$cell<-extract(cellid, ebd.sub[,c("LONGITUDE","LATITUDE")])
ebd.sub<-ebd.sub[-which(is.na(ebd.sub$cell)),]
rm(ebd.all)

#Estimation from transect data
#Id cells with sufficient sampling (more than 4hrs of sampling)
ind.recs<-which(ebd.sub$DURATION.MINUTES>=180 & 
                  ebd.sub$PROTOCOL.TYPE=="eBird - Traveling Count" & 
                  ebd.sub$EFFORT.DISTANCE.KM < 5 &
                  ebd.sub$year>=2012)

ebd.sub1<-ebd.sub[ind.recs,]

#Reformat matrix 
M<-dcast(ebd.sub1,OBSERVATION.DATE+cell~SCIENTIFIC.NAME,fun.aggregate=length)
Mi<-M
Mi[,3:ncol(Mi)] <- (Mi[,3:ncol(Mi)]>0)*1 #Incidence-based matrix

#Find cells sufficiently sampled
singNdoub<-function(cell,siteBySp){
  sub.m<-siteBySp[siteBySp$cell==cell, ]
  nrecs<-apply(sub.m[,3:ncol(sub.m)], 2, sum)
  nobs<-sum(nrecs>0) #Species observed
  smps<-nrow(sub.m) #Samples
  sing<-sum(nrecs==1) #Singletons
  doub<-sum(nrecs==2) #Doubletons
  suf<-(sing/nobs<0.5)*1 #Criteria to define a place as sufficiently sampled according to Anne Chao in EstimateS User Guide
  return(list(nobs=nobs, smps=smps, sing=sing, doub=doub, suf=suf))
}

results<-sapply(unique(Mi$cell), singNdoub, siteBySp=Mi)
results<-t(results)
results<-as.data.frame(results)
results$cell<-unique(Mi$cell)
results[,"nobs"]<-unlist(results[,"nobs"])
results[,"smps"]<-unlist(results[,"smps"])
results[,"sing"]<-unlist(results[,"sing"])
results[,"doub"]<-unlist(results[,"doub"])
results[,"suf"]<-unlist(results[,"suf"])
results[,"cell"]<-unlist(results[,"cell"])

selected.cells<-results$cell[results$suf==1]
Msub<-Mi[Mi$cell%in%selected.cells, ]

#Estimate richness
Rtransect<-specpool(Msub[,3:ncol(M)], Msub[,2])
Rtransect$cell<-as.numeric(row.names(Rtransect))
xyR<-xyFromCell(cellid,Rtransect$cell)
Rtransect<-cbind(Rtransect,xyR)

#Compute total richness
spp.obs<-rasterize(ebd.sub[,c("LONGITUDE","LATITUDE")], cellid, 
          field=ebd.sub$SCIENTIFIC.NAME,
          fun=function(x, ...){length(unique(x))})
writeRaster(spp.obs,"M:/BirdRichness/observedRichness.tif")

Rtransect$all.spp<-extract(spp.obs,Rtransect[,c("x","y")])
write.csv(Rtransect,"M:/BirdRichness/estimatedRichness.csv",row.names=F)
#Processed some columns in excel
Rtransect2<-read.csv("M:/BirdRichness/estimatedRichness_v2.csv",as.is=T)

#Interpolate richness
env.data<-stack(paste0("Z:/Capas/NewExtent/","bio_",c(1,2,3,4,12,15,18),".asc"))
clouds<-stack(paste0("C:/Capas/aoi/cloud/",c("MODCF_intraannualSD","MODCF_meanannual",
                     "MODCF_spatialSD_1deg"),".tif"))
tree<-raster("C:/Capas/aoi/forestHeight.tif")
env.data<-addLayer(env.data,clouds,tree)
env.data5km<-aggregate(env.data,fact=5,fun=mean)
veg<-stack(list.files("Z:/Capas/landModis","*.tif$",full.names=T))
env.data5km<-addLayer(env.data5km,veg)

covs<-extract(env.data5km, Rtransect2[,c("x","y")])

df2<-cbind(Rtransect2, covs)
df2$logjack2<-log(df2$jack2)
df2$roundjack2<-round(df2$jack2)
ind<-which(Rtransect2$jack2/Rtransect2$all.spp>0.95) #select cells that have more than 
#95% congruence between estimates from jack2 and species recorded
df3<-df2[ind,]

#points(Rtransect2[which(Rtransect2$jack2/Rtransect2$all.spp>0.95),c("x","y")],pch=18)

#Partition data into spatial folds
chk2.pts <- get.checkerboard2(df3[,c("x","y")],env.data5km, df3[,c("x","y")], 4)
df3$group<-chk2.pts$occ.grp

#Normal brt all variables
df3$pred_n1<-NA
for(i in 1:4){
  n1<-gbm(logjack2~bio_1+bio_2+bio_3+bio_4+bio_12+bio_15+bio_18+
            MODCF_intraannualSD+MODCF_meanannual+landClass_0+landClass_1+
            landClass_2+landClass_3+landClass_4+landClass_5+landClass_6+
            landClass_7+landClass_8+landClass_9+landClass_10+landClass_11+
            landClass_12+landClass_13+landClass_14+landClass_15+landClass_16+
            landShannoInd,
          data=df3[df3$group!=i, ],
          distribution="gaussian",
          n.trees=5000,
          shrinkage=0.01,
          interaction.depth=3,
          cv.folds = 5,
          keep.data=TRUE,
          verbose=FALSE,
          n.cores=5)
  best.iter <- gbm.perf(n1, method="cv")
  p1<-exp(predict(env.data5km,n1,n.trees=best.iter,type="response"))
  if(exists("n1.stack")){
    n1.stack<-addLayer(n1.stack,p1)
  } else {
    n1.stack<-stack(p1)
  }
  df3[df3$group==i, "pred_n1"]<-extract(p1, df3[df3$group==i,c("x","y")])
}

mean(abs(df3$pred_n1-df3$jack2)) #Mean absolute error 68.27366

#Normal brt without precp variables
df3$pred_n2<-NA
for(i in 1:4){
  n2<-gbm(logjack2~bio_1+bio_2+bio_3+bio_4+
            MODCF_intraannualSD+MODCF_meanannual+landClass_0+landClass_1+
            landClass_2+landClass_3+landClass_4+landClass_5+landClass_6+
            landClass_7+landClass_8+landClass_9+landClass_10+landClass_11+
            landClass_12+landClass_13+landClass_14+landClass_15+landClass_16+
            landShannoInd,
          data=df3[df3$group!=i, ],
          distribution="gaussian",
          n.trees=5000,
          shrinkage=0.01,
          interaction.depth=3,
          cv.folds = 5,
          keep.data=TRUE,
          verbose=FALSE,
          n.cores=5)
  best.iter <- gbm.perf(n2, method="cv")
  p2<-exp(predict(env.data5km,n2,n.trees=best.iter,type="response"))
  if(exists("n2.stack")){
    n2.stack<-addLayer(n2.stack,p2)
  } else {
    n2.stack<-stack(p2)
  }
  df3[df3$group==i, "pred_n2"]<-extract(p2, df3[df3$group==i,c("x","y")])
}

mean(abs(df3$pred_n2-df3$jack2)) #Mean absolute error jack2 68.30411
mean(abs(df3$pred_n2-df3$all.spp)) #Mean absolute error all.sp 82.5029

#Poisson BRT all variables
df3$pred_n3<-NA
for(i in 1:4){
  n3<-gbm(roundjack2~bio_1+bio_2+bio_3+bio_4+bio_12+bio_15+bio_18+
            MODCF_intraannualSD+MODCF_meanannual+forestHeight,
          data=df3[df3$group!=i, ],
          distribution="poisson",
          n.trees=5000,
          shrinkage=0.01,
          interaction.depth=3,
          cv.folds = 5,
          keep.data=TRUE,
          verbose=FALSE,
          n.cores=5)
  best.iter <- gbm.perf(n3, method="cv")
  p3<-predict(env.data5km,n3,n.trees=best.iter,type="response")
  if(exists("n3.stack")){
    n3.stack<-addLayer(n3.stack,p3)
  } else {
    n3.stack<-stack(p3)
  }
  df3[df3$group==i, "pred_n3"]<-extract(p3, df3[df3$group==i,c("x","y")])
}

mean(abs(df3$pred_n3-df3$jack2)) #Mean absolute error jack2 67.85034

#Poisson BRT all variables
df3$pred_n4<-NA
for(i in 1:4){
  n4<-gbm(roundjack2~bio_1+bio_2+bio_3+bio_4+bio_12+bio_15+bio_18+
            MODCF_intraannualSD+MODCF_meanannual+forestHeight,
          data=df3[df3$group!=i, ],
          distribution="poisson",
          n.trees=5000,
          shrinkage=0.01,
          interaction.depth=3,
          cv.folds = 5,
          keep.data=TRUE,
          verbose=FALSE,
          n.cores=5)
  best.iter <- gbm.perf(n4, method="cv")
  p4<-predict(env.data5km,n4,n.trees=best.iter,type="response")
  if(exists("n4.stack")){
    n4.stack<-addLayer(n4.stack,p4)
  } else {
    n4.stack<-stack(p4)
  }
  df3[df3$group==i, "pred_n4"]<-extract(p4, df3[df3$group==i,c("x","y")])
}

mean(abs(df3$pred_n4-df3$jack2)) #Mean absolute error jack2 67.93849

write.csv(df3, "M:/BirdRichness/interpolationResults.csv",row.names=F)





#Poisson brt

df2$logjack1<-log(df2$jack1)
r2 <- gbm.step(data=df2, gbm.x = 12:18, gbm.y = 19,
  family = "gaussian", tree.complexity = 5,
  learning.rate = 0.001, bag.fraction = 0.5)

df2$jack1.round<-round(df2$jack1)
df2$logjack1<-log(df2$jack1)

r3 <- gbm.step(data=df2[ind,], gbm.x = c(12:15,19:22), gbm.y = 25,
               family = "gaussian", tree.complexity = 3,
               learning.rate = 0.01, bag.fraction = 0.5)

#Climate+Clouds+Land (all)

r4<-gbm(jack1.round~bio_1+bio_2+bio_3+bio_4+MODCF_intraannualSD+
          MODCF_meanannual+forestHeight,data=df2,
        distribution="poisson",
        n.trees=5000,
        shrinkage=0.01,
        interaction.depth=3,
        cv.folds = 10,
        keep.data=TRUE,
        verbose=FALSE,
        n.cores=5)
best.iter <- gbm.perf(r4,method="cv")
p3<-predict(env.data,r4,n.trees=best.iter,type="response")
plot(exp(p2))
points(df[,c("x","y")])

#Report square MSE
sqrt(sum((df2$jack1.round-r4$cv.fitted)^2)/nrow(df2))
#194.5233
#Report COR
cor(df2$jack1.round,r4$cv.fitted)
#0.1760886

#Climate+Land+Clouds(C>80%)
r5<-gbm(jack1.round~bio_1+bio_2+bio_3+bio_4+MODCF_intraannualSD+
          MODCF_meanannual+forestHeight,data=df2[ind,],
        distribution="poisson",
        n.trees=5000,
        shrinkage=0.01,
        interaction.depth=3,
        cv.folds = 10,
        keep.data=TRUE,
        verbose=FALSE,
        n.cores=5)
best.iter <- gbm.perf(r5,method="cv")


#Report square MSE
sqrt(sum((df2$jack1.round[ind]-r5$cv.fitted)^2)/nrow(df2[ind,]))
#254.0505
#Report COR
cor(df2$jack1.round[ind],r5$cv.fitted)
#0.002899747

#Climate+Cloud(all)
r6<-gbm(jack1.round~bio_1+bio_2+bio_3+bio_4+MODCF_intraannualSD+
          MODCF_meanannual,data=df2,
        distribution="poisson",
        n.trees=5000,
        shrinkage=0.01,
        interaction.depth=3,
        cv.folds = 10,
        keep.data=TRUE,
        verbose=FALSE,
        n.cores=5)
best.iter <- gbm.perf(r6,method="cv")


#Report square MSE
sqrt(sum((df2$jack1.round-r6$cv.fitted)^2)/nrow(df2))
#254.0505
#Report COR
cor(df2$jack1.round[ind],r5$cv.fitted)

#****IDEAS*****


# Use only cells with complete surveys according to three estimators 
# (chao2, jack1,jack2)
# Find out whats the length of survey (time, distance) after which
# a day-survey is nearly complete, perhaps per habitat type.
# Use internal validation with spatial cross-validation
# Use external validation at different resolutions using sites with
# complete inventories or guesstimates from experts
# Use rs data for interpolation
# produce uncertainty maps






df3$response<-log(df3$response)
r2 <- gbm.step(data=df, gbm.x = 12:18, gbm.y = 19,
                   family = "gaussian", tree.complexity = 5,
                   learning.rate = 0.001, bag.fraction = 0.5)
p2<-predict(env.data,r2,n.trees=r2$gbm.call$best.trees,type="response")
writeRaster(exp(p2),"c:/workspace/potential_richness_jack1a.tif")
col<-readOGR("D:/Datos/BaseLayers","COL_adm1")

###


ind<-which(ebd.all$ALL.SPECIES.REPORTED==1&
             (ebd.all$PROTOCOL.TYPE=="eBird - Traveling Count"|
                ebd.all$PROTOCOL.TYPE=="eBird - Stationary Count")&
             ebd.all$CATEGORY=="species")

ebd.sub<-ebd.all[ind, c("LONGITUDE","LATITUDE","SCIENTIFIC.NAME",
                     "SAMPLING.EVENT.IDENTIFIER","OBSERVATION.DATE")]


#Filter to Colombia
env.col<-crop(env.data,col)
env.col<-mask(env.col,col)
sp.col<-ebd.sub[which(!is.na(extract(env.col, ebd.sub[,1:2]))),]
template<-env.col[[1]]
template[1:ncell(template)]<-1:ncell(template)
template[is.na(env.col[[1]])]<-NA
sp.col$inCol<-extract(template, sp.col[,1:2])
sp.col<-na.omit(sp.col)
sp.col<-unique(sp.col)

#Structure sample by species matrix for all sites
M<-dcast(sp.col,OBSERVATION.DATE+inCol~SCIENTIFIC.NAME,fun.aggregate=length)

R<-specpool(M[,3:1663], M[,2])
xyR<-xyFromCell(template,as.numeric(row.names(R)))
R<-cbind(R,xyR)
R$cv<-sqrt(R$n)*R$chao.se/R$chao
write.csv(R,"C:/Workspace/RichnessEst2.csv",row.names = FALSE)

#Include abundance estimations
#Include other census types (e.g. incidental observations?)
#Or on the contrary, only censuses beyond a certain temporal length
#Or maybe a minimum number of censuses?


SamplesSpeciesM<-function(x){
  dcast(x,OBSERVATION.DATE~SCIENTIFIC.NAME,fun.aggregate=length)
  return()
}


dcast(sp.boy,inCol~SCIENTIFIC.NAME,fun.aggregate = function(x){
  return((length(x)>0)*1)
})

#Select only transect or point count with more than 4 h (vs. everything?)
#Category all
#make sure date is after 1999

incidence<-dcast(sp.boy,inCol~SCIENTIFIC.NAME,fun.aggregate = function(x){
  return((length(x)>0)*1)
})

singles<-dcast(sp.boy,inCol~SCIENTIFIC.NAME,fun.aggregate = function(x){
  return((length(x)==1)*1)
  })

doubles<-dcast(sp.boy,inCol~SCIENTIFIC.NAME,fun.aggregate = function(x){
  return((length(x)==2)*1)
})

samples<-dcast(sp.boy,inCol~.,value.var="OBSERVATION.DATE",
               fun.aggregate = function(x){
                                  return(nsamples=length(unique(x)))
})

ibycell<-apply(incidence[,2:ncol(incidence)], 1, sum)
sbycell<-apply(singles[,2:ncol(singles)], 1, sum)
dbycell<-apply(doubles[,2:ncol(doubles)], 1, sum)
samples<-samples[,2]
est.richness<-ibycell+((samples-1)/samples)*(sbycell*(sbycell-1))/(2*(dbycell+1))
completeness<-ibycell*100/est.richness

#Compute lower and upper confidence bounds for completeness and est.richness
result<-cbind(xyFromCell(template,incidence$inCol),Nsamples=samples, Singles=sbycell,
              Doubles=dbycell,Sobs=ibycell,Sest=est.richness,C=completeness)
write.csv(result,"c:/workspace/result3.csv",row.names=FALSE)
write.csv(sp.boy,"c:/workspace/sp.boy.csv",row.names=FALSE)
