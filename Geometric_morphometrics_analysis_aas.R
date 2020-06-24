

####### CULTURAL CONTINUITIES AND DISCONTINUITIES AT THE NEOLITHIC TRANSITION IN EASTERN
####### IBERIA: AN ANLYSIS OF THE MORPHOMETRY OF THE GEOMETRIC MICROLITHS 

#### Archaeological and Anthropological Sciences

#### Article authors: A. Cortell-Nicolau, O. García-Puchol, S. Shennan

#### Code author: A. Cortell-Nicolau

#### London, 2018-2019 

#### Reviewed: València, 2020

#### License: Permission is granted to use and adapt this code. Please do not distribute
#### without permission of the author

#### Warranty: No warranty is expressed or implied. Use at your own risk

#### About: This performs an outline geometric morphometrics analysis. Outlines are 
#### provided as online resource 2. 

#### Obtained results may be slightly different from the publication. See code line 49.

#### GEOMETRIC MORPHOMETRICS ANALYSIS

## library
library(Momocs)

## Data
geo_df <- read.csv('/Users/acortell/Desktop/GeoMorph/geo_df.csv', header = FALSE) ## Your route
geo_df <- geo_df[,c(1,6)]
colnames(geo_df)<-c("Fase","ID")

geo_df<-geo_df[order(geo_df$ID),]

## Group C phases
geo_df$Fase<-sub("C1","C",geo_df$Fase)
geo_df$Fase<-sub("C2","C",geo_df$Fase)

geo_df$Fase<-as.factor(geo_df$Fase) ## Phases as factors

## Extract outlines
setwd('/Users/acortell/Desktop/GeoMorph/outlines') ## Folder where the outlines are stored
geo.list<-list.files(pattern="*", recursive = TRUE) ## Order outlines
#geo.list

geo2<-Out(import_jpg(geo.list,auto.notcentered= T), fac=geo_df) ## Because auto.notcenterd = T, results may not be exactly the same. See package documentation

## overlay
# stack(coo_center(coo_scale(geo2)))

## Elliptical fourier analysis at n=20 harmonics
geo2F <- efourier(coo_center(coo_scale(geo2)), nb.h = 20, norm=F, start=F)

## PCA
geo2.pca<-PCA(geo2F)
#plot(geo2.pca, "Fase")
plot(geo2.pca, 1, ellipses=TRUE, ellipsesax=FALSE, pch=c(4,5))
#plot(geo2.pca, 1, chull=TRUE, pos.shp="full_axes", abbreviate.labelsgroups=TRUE, points= FALSE, labelspoints=TRUE)
#plot(geo2.pca, 1, pos.shp="circle",stars=TRUE, palette=col_qual)

#scree_plot(geo2.pca)

## LDA
geo2.lda<-LDA(geo2F, fac=geo_df$Fase)
#plot_LDA(geo2.lda,points=TRUE, chull = FALSE, zoom = 0.7, labelgroups = TRUE)
plot_LDA(geo2.lda, chullfilled = TRUE, zoom = 1, box = TRUE, points = TRUE)

## Iterations (n=25) for cross validation table
cvn40<-LDA(geo2F, fac=geo_df$Fase)
cvn40<-cvn40$CV.tab

for (i in 1:25){
  provcv40<-LDA(geo2F, fac=geo_df$Fase)
  provcv40<-provcv40$CV.tab
  cvn40<-cvn40+provcv40
}


plot_CV(cvn40, freq = FALSE)
plot_CV2(geo2.lda)

## Pariwise MANOVA
geo2.man<-MANOVA_PW(geo2.pca, fac=geo_df$Fase)
geo2tab<-geo2.man$summary

## K-Means
KMEANS(geo2.pca, centers=5)

