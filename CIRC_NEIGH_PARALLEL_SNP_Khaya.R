# Circular neighborhood///////////////////////////////////////////////////
library(raster)
library(parallel)

#########################################
#########################################
#########################################
# only adjust this data based on your situation (also see lower)
inputdir <- getwd()
outputdir <- getwd()

#########################
# create a folder called 'maps'in outputdir (with windows explorer), and withing that folder tow additional folders called "asc" "kml"
#########################
setwd(inputdir)
inputfile<-"Khaya2.csv"
ssr_data<-read.csv(inputfile, stringsAsFactors = FALSE)
nutpoints<-ssr_data[, 1:4]
# nutpoints<-nutpoints[, 1:16]
a<-5
while (a<ncol(ssr_data)){
  nutpoints<-cbind(nutpoints, paste(ssr_data[, a], "/", ssr_data[, a+1], sep="")) 
  names(nutpoints)[ncol(nutpoints)]<-names(ssr_data)[a]
  nutpoints[,ncol(nutpoints)]<-as.character(nutpoints[,ncol(nutpoints)])
  a=a+2
}


Nloci<-ncol(nutpoints)-4
RESOLUTION_minutes<-10
CIRCULARNEIGHBORHOOD_diam_minutes<-60
sp_name<-"Khaya"
#########################################
#########################################
#########################################




names(nutpoints)[3:4]<-c("X","Y")

RESOLUTION<-RESOLUTION_minutes/60
CIRCULARNEIGHBORHOOD<-CIRCULARNEIGHBORHOOD_diam_minutes/60/2

r<-raster(xmn=min(nutpoints$X)-2, xmx=max(nutpoints$X)+2, ymn=min(nutpoints$Y)-2, ymx=max(nutpoints$Y)+2)
res(r)<-RESOLUTION
projection(r)<-"+proj=longlat +datum=WGS84"
r



circ.neigh<-function(nutpoints, x, r, CIRCULARNEIGHBORHOOD){
  cn<-rep(0, ncol(nutpoints)+3)
  po_co<-subset(nutpoints[x,], select=c(X,Y))
  exte<-extent(po_co$X-CIRCULARNEIGHBORHOOD,po_co$X+CIRCULARNEIGHBORHOOD,
               po_co$Y-CIRCULARNEIGHBORHOOD,po_co$Y+CIRCULARNEIGHBORHOOD)
  window<-crop(r, exte)
  cellcenters<-xyFromCell(window, cell=1:ncell(window), spatial=F)
  distpointcells<-pointDistance(po_co, cellcenters, longlat=F)
  names(distpointcells)<-cellFromXY(r,cellcenters)
  temp1<- distpointcells[distpointcells<=CIRCULARNEIGHBORHOOD]
  if (length(temp1)>=1) {
    for (l in 1:length(temp1)){
      n<-as.numeric(names(temp1[l]))
      temp2<-cbind(nutpoints[x,],n,xyFromCell(r,cell=n, spatial=F))
      cn<-rbind(cn,temp2)
    }
  }
  return(as.data.frame(cn[-1,]))
  return(print(paste(x , 'out of', nrow(nutpoints), sep="")))}



library(parallel)
cl <- makeCluster(2)
clusterExport(cl, c("nutpoints", "r", "CIRCULARNEIGHBORHOOD", "circ.neigh", "extent", 
                    "cellFromXY", "xyFromCell", "crop", "ncell", "pointDistance"))
intermRes<-parSapply(cl, 1:nrow(nutpoints), simplify = F, FUN=function(x) circ.neigh(nutpoints, x, r, CIRCULARNEIGHBORHOOD))
stopCluster(cl)    

cn_t<-intermRes[[1]]
for(n in 2:length(intermRes)){
  cn_t<-rbind(cn_t, intermRes[[n]])
}


cn<-cn_t
cn_trees<-data.frame(cbind(1:nrow(cn), cn[,(Nloci+5):(Nloci+7)], cn[,5:(Nloci+4)]))
colnames(cn_trees)<-c("ID", "cell","Xnew", "Ynew", colnames(cn_trees[5:ncol(cn_trees)]))
plot(cn_trees[,3:4])

setwd(outputdir)
write.csv(cn_trees, file=paste(sp_name, "_cn_ssr_data",RESOLUTION_minutes,"min_cn", CIRCULARNEIGHBORHOOD_diam_minutes,"min.csv"),row.names = F)
cn_trees<-read.csv(paste(sp_name, "_cn_ssr_data",RESOLUTION_minutes,"min_cn", CIRCULARNEIGHBORHOOD_diam_minutes,"min.csv"))
# Bootstrap/////////////////////////////////////////////////////////////
# setwd("C:/Users/calcazar/Dropbox/SAP - Spatial analyses/output/genetic diversity/Phaseolus")
# cn_trees<-read.csv("cn_ssr_data 10 min_cn 60 min.csv", header=T)



library(adegenet)
library(plyr)
library(stringr)
plants_cn<-cn_trees
names(plants_cn)[5:ncol(plants_cn)]<-paste("L", 1:(ncol(plants_cn)-4), sep="")
for (i in 4:69){
  plants_cn[,i]<-as.character(plants_cn[,i])
}
# setwd("C:/Users/calcazar/Dropbox/SAP - Spatial analyses/script")
source("df2genind2.R")


system.time(plants_cn1<-df2genind(plants_cn[,5:ncol(plants_cn)], sep="/",ind.names=plants_cn[,1], loc.names=NULL, pop=plants_cn[,2], NA.char=NA, ploidy=2, type="codom"))

plants_cn2<-genind2genpop(plants_cn1)
A<-plants_cn1
B<-plants_cn2


# Observed & Expected heterozygocity, Fixation index, Locally Common Alleles
pops<-as.numeric(as.matrix(pop(A)))
t1<-cbind(pops,A$tab)
t1<-as.data.frame(t1)
t2 <-B$tab

#########################################
#########################################
#########################################

subsample<-3
# subsample<-floor(median(table(cn_trees[2])))
iterations<-1:1000

#########################################
#########################################
######## read in functions###############

prop.of.pops <- function(x) length(x[x>0])/length(x) 
prop.pops <- apply(t2, 2, prop.of.pops)


bootstr <- function(t1, x, subsample, iterations){
  temp <- t1[t1[,1] == x,2:ncol(t1)]
  
 
  n_He<- rep(NA, length(iterations))
  sd_n_He<-rep(NA, length(iterations))
  n_all<- rep(NA, length(iterations))
  sd_n_all<- rep(NA, length(iterations))
  n_lca<-rep(NA, length(iterations))
  sd_n_lca<-rep(NA, length(iterations))
  n_lca10<-rep(NA, length(iterations))
  sd_n_lca10<-rep(NA, length(iterations))
  n_Het<-rep(NA, length(iterations))
  sd_n_Het<-rep(NA, length(iterations))
  n_F<-rep(NA, length(iterations))
  sd_n_F<-rep(NA, length(iterations))
  Sh<-rep(NA, length(iterations))
  sd_Sh<-rep(NA, length(iterations))
  
  if (nrow(temp) >= subsample){
    for (r1 in iterations){
      temporary<- temp[sample(1:nrow(temp), subsample ,replace=F),]
      allelespop<-apply(temporary,2, sum, na.rm=TRUE)
      NAl<-allelespop>0
      NAl<-sapply(unique(str_sub(names(NAl), 1,3)), function(loci) sum(NAl[grep(loci, names(NAl), fixed=T)],na.rm=TRUE))
      n_all[r1] <- mean(NAl[NAl!=0],na.rm=T)
      sd_n_all[r1]<-sd(NAl, na.rm=T)
      tmp<-apply(temporary,2, function(x) sum(x==1, na.rm=TRUE))
      Ns<-apply(temporary,2, function(x) length(x[!is.na(x)]))
      tmp<-sapply(str_sub(names(tmp), 1,3), function(loci) sum(tmp[grep(loci, names(allelespop), fixed=T)], na.rm=T))/2/Ns
      tmp<-tmp[unique(names(tmp)) ]
      n_Het[r1]<-mean(tmp, na.rm=T)
      sd_n_Het[r1]<-sd(tmp, na.rm=T)
      afpop <-allelespop / sapply(str_sub(names(allelespop), 1,3), function(loci) sum(allelespop[ grep(loci, names(allelespop), fixed=T)],na.rm=T)) 
      lc <- prop.pops <= 0.25 & afpop >= 0.05
      lc25<-sapply(unique(str_sub(names(lc), 1,3)), function(loci) sum(lc[grep(loci, names(lc),fixed=T)],na.rm=TRUE))
      n_lca[r1]<- mean(lc25,na.rm=T)
      sd_n_lca[r1]<- sd(lc25,na.rm=T)
      lc10 <- prop.pops <= 0.1 & afpop >= 0.05
      lc10<-sapply(unique(str_sub(names(lc10), 1,3)), function(loci) sum(lc10[grep(loci, names(lc10), fixed=T)],na.rm=TRUE))
      n_lca10[r1]<-mean(lc10,na.rm=T)
      sd_n_lca10[r1]<-sd(lc10,na.rm=T)
      HZ <- afpop^2
      HZ <- sapply(unique(str_sub(names(HZ), 1,3)), function(loci) sum(HZ[grep(loci, names(HZ),fixed=T)],na.rm=TRUE))
      eHZ <- 1-HZ
      n_He[r1]<-mean(eHZ, na.rm=TRUE)
      sd_n_He[r1]<-sd(eHZ, na.rm=TRUE)
      n_F[r1]<-mean((eHZ-tmp)/eHZ,na.rm=T)
      sd_n_F[r1]<-sd((eHZ-tmp)/eHZ,na.rm=T)
      S <- afpop*log(afpop)
      S <- sapply(unique(str_sub(names(S), 1,3)), function(loci) -sum(S[grep(loci, names(S), fixed=T)],na.rm=TRUE))
      Sh[r1]<-mean(S,na.rm=TRUE)
      sd_Sh[r1]<-sd(S,na.rm=TRUE)
    } 
    return(c(mean(n_Het, na.rm=TRUE), mean(sd_n_Het, na.rm=TRUE), mean(n_all, na.rm=TRUE), mean(sd_n_all, na.rm=TRUE), 
             mean(n_He, na.rm=TRUE), mean(sd_n_He, na.rm=TRUE), mean(n_lca, na.rm=TRUE), mean(sd_n_lca, na.rm=TRUE), 
             mean(n_lca10, na.rm=TRUE), mean(sd_n_lca10, na.rm=TRUE), mean(n_F, na.rm=TRUE),  mean(sd_n_F, na.rm=TRUE),
             mean(Sh,na.rm=TRUE), mean(sd_Sh,na.rm=TRUE)))
  } else return(rep(NA, 7))
  print(paste(which(unique(t1[,1])==x), 'out of', length(unique(t1[,1])), sep=""))
}




# uncorrected parameters
He_u<-rep(NA, length(unique(t1[,1])))
lca_u<-rep(NA, length(unique(t1[,1])))
lca_u10<-rep(NA, length(unique(t1[,1])))
Ho_u<-rep(NA, length(unique(t1[,1])))
Fix_u<-rep(NA, length(unique(t1[,1])))
richness_u<-rep(NA, length(unique(t1[,1])))
I_u<-rep(NA, length(unique(t1[,1])))
sd_He_u<-rep(NA, length(unique(t1[,1])))
sd_lca_u<-rep(NA, length(unique(t1[,1])))
sd_lca_u10<-rep(NA, length(unique(t1[,1])))
sd_Ho_u<-rep(NA, length(unique(t1[,1])))
sd_Fix_u<-rep(NA, length(unique(t1[,1])))
sd_richness_u<-rep(NA, length(unique(t1[,1])))
sd_I_u<-rep(NA, length(unique(t1[,1])))

m<-1

for (i in unique(t1[,1])) {
  temp <- t1[t1[,1] == i,2:ncol(t1)]
  
      
  # ///////////////////uncorrected values/////////////////// 
  allelespop<-apply(temp,2, sum, na.rm=TRUE)
  NAl<-allelespop>0
  NAl<-sapply(unique(str_sub(names(NAl), 1,3)), function(loci) sum(NAl[grep(loci, names(NAl), fixed=T)],na.rm=TRUE))
  richness_u[m]<- mean(NAl[NAl!=0], na.rm=T)
  sd_richness_u[m]<-sd(NAl, na.rm=T)
  tmp<-apply(temp,2, function(x) sum(x==1, na.rm=TRUE))
  Ns<-apply(temp,2, function(x) length(x[!is.na(x)]))
  tmp<-sapply(str_sub(names(tmp), 1,3), function(loci) sum(tmp[grep(loci, names(allelespop), fixed=T)], na.rm=T))/2/Ns
  tmp<-tmp[unique(names(tmp)) ]
  Ho_u[m]<-mean(tmp, na.rm=T)
  sd_Ho_u[m]<-sd(tmp, na.rm=T)
  afpop <-allelespop / sapply(str_sub(names(allelespop), 1,3), function(loci) sum(allelespop[ grep(loci, names(allelespop), fixed=T)],na.rm=T)) 
  lc <- prop.pops <= 0.25 & afpop >= 0.05
  lc25<-sapply(unique(str_sub(names(lc), 1,3)), function(loci) sum(lc[grep(loci, names(lc), fixed=T)],na.rm=TRUE))
  lca_u[m]<- mean(lc25, na.rm=T)
  sd_lca_u[m]<- sd(lc25, na.rm=T)
  lc10 <- prop.pops <= 0.1 & afpop >= 0.05
  lc10<-sapply(unique(str_sub(names(lc10), 1,3)), function(loci) sum(lc10[grep(loci, names(lc10), fixed=T)],na.rm=TRUE))
  lca_u10[m]<-mean(lc10,na.rm=T)
  sd_lca_u10[m]<-sd(lc10,na.rm=T)
  HZ <- afpop^2
  HZ <- sapply(unique(str_sub(names(HZ), 1,3)), function(loci) sum(HZ[grep(loci, names(HZ), fixed=T)],na.rm=TRUE))
  eHZ <- 1-HZ
  He_u[m]<-mean(eHZ, na.rm=TRUE)
  sd_He_u[m]<-sd(eHZ, na.rm=TRUE)
  Fix_u[m]<-mean((eHZ-tmp)/eHZ, na.rm=T)
  sd_Fix_u[m]<-sd((eHZ-tmp)/eHZ, na.rm=T)
  S <- afpop*log(afpop)
  S <- sapply(unique(str_sub(names(S), 1,3)), function(loci) -sum(S[grep(loci, names(S), fixed=T)],na.rm=TRUE))
  I_u[m]<-mean(S,na.rm=TRUE)
  sd_I_u[m]<-sd(S,na.rm=TRUE)
  
  print (m)
  m<-m+1
}
  #///////////////////////bootstrap///////////////  

# Make cluster
library(parallel)
cl <- makeCluster(2)
clusterExport(cl, c("t1", "subsample", "Nloci", "str_sub", "prop.pops", "bootstr", "iterations"))
intermRes<-parSapply(cl, unique(t1[,1]), FUN=function(x) bootstr(t1, x, subsample, iterations))
stopCluster(cl)     

bootres<-t(as.data.frame(intermRes))
colnames(bootres)<-c(paste("Ho_",subsample,"trees", sep=""), paste("sd_Ho_",subsample,"trees", sep=""),
                     paste("rich",subsample,"trees", sep=""),  paste("sd_rich",subsample,"trees", sep=""),
                     paste("He_",subsample,"trees", sep=""), paste("sd_He_",subsample,"trees", sep=""), 
                     paste("lca25_",subsample,"trees", sep=""), paste("sd_lca25_",subsample,"trees", sep=""),
                     paste("lca10_",subsample,"trees", sep=""), paste("sd_lca10_",subsample,"trees", sep=""),
                     paste("Fix_",subsample,"trees", sep=""),paste("sd_Fix_",subsample,"trees", sep=""),
                     paste("Shannon_",subsample,"trees", sep=""), paste("sd_Shannon_",subsample,"trees", sep=""))
  
setwd(outputdir)
# write.csv(bootres, 'bootstrapAstronium3trees_cn10min.csv')


dupli<-bootres
dumi<-cbind(richness_u, sd_richness_u,lca_u, sd_lca_u, lca_u10, sd_lca_u10, Ho_u, sd_Ho_u, He_u, sd_He_u, I_u, sd_I_u, Fix_u, sd_Fix_u)
colnames(dumi)<-c('rich_u', 'sd_rich_u', 'lca_u', 'sd_lca_u', 'lca_u10', 'sd_lca_u10', 'Ho_u', 'sd_Ho_u', 'He_u', 'sd_He_u',
                  'I_u', 'sd_I_u', 'Fix_u', 'sd_Fix_u')


coords<-cn_trees[,2:4][!duplicated(cn_trees$cell), ]
results<-data.frame(coords, bootres, dumi)
# colnames(results)<-c("cellnr","x", "Y", "rich_no_correction", paste("rich",subsample,"trees", sep=""), "Shannon_no_correction", paste("Shannon_",subsample,"trees", sep=""), "lca25_no_correction", "lca10_no_correction", paste("lca25_",subsample,"trees", sep=""), paste("lca10_",subsample,"trees", sep=""), "Ho_no_correction", paste("Ho_",subsample,"trees", sep=""), "He_no_correction", paste("He_",subsample,"trees", sep=""), "Fix_no_correction", paste("Fix_",subsample,"trees", sep="")  )
write.csv(results, file=paste("boot", sp_name, "_", subsample, "_plants_cn", CIRCULARNEIGHBORHOOD_diam_minutes, "min_res", RESOLUTION_minutes, "complete.csv", sep=""),row.names = F)
nas<-which(is.na(results[,5]))
results_noNA<-results[-nas,]
write.csv(results_noNA, file=paste("boot", sp_name, "_", subsample, "_plants_cn", CIRCULARNEIGHBORHOOD_diam_minutes, "min_res", RESOLUTION_minutes, "noNA_complete.csv", sep=""),row.names = F)


# ramp<-colorRamp(c(rgb(255,255,128, max=255),rgb(242,167,46, max=255),rgb(107,0,0, max=255)))
# ramp<-colorRamp(c("yellow", "brown"))
for (a in 4:31){
  
  val<- results[,a]
  coor<-results[,2:3]
  coor<-coor[!is.na(val),]
  val<-val[!is.na(val)]
  rast<-rasterize(coor,r,field=val, ext=extent(r))
  rast<-trim(rast)
  writeRaster(rast, filename=paste("maps/", colnames(results)[a], "_", sp_name, "_", subsample, "_plants_cn", CIRCULARNEIGHBORHOOD_diam_minutes, "min_res", RESOLUTION_minutes, ".asc", sep=""),overwrite=T)
#   KML(rast, col=rgb(ramp(seq(0, 1, length = 40)), max = 255), blur=15, file=paste("maps/kml/", colnames(results)[a], "_", sp_name, "_", subsample, "_plants_cn", CIRCULARNEIGHBORHOOD_diam_minutes, "min_res", RESOLUTION_minutes, ".kml", sep=""))
  
}
library(maptools)
coordinates(plantpoints)<- ~X+Y
kmlPoints(plantpoints,kmlfile="Vigourouxpoints.kml",icon="http://google.com/mapfiles/kml/paddle/red-square-lv.png", name="")