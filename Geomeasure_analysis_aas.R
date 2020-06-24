

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

#### About: This script uses the information obtained using the beta version of the R
#### package Geomeasure (in development, expected for summer 2020) for performing
#### different ecological and statistical analysis. Article data is provided in 
#### online material 1

#### GEOMEASURE ANALYSIS

## 1. LOAD PACKAGES AND DATA
## packages
library(ggplot2)
library(plyr) ##Changes names of factors
library(devtools) ## Necessary for instaling ggbiplot
library(ggbiplot)
library(MASS) ## for huber M estimator
library(pheatmap) ## jaccard heatmaps
library(vegan) ## jaccard
library(RColorBrewer) ## nice colors
library(reshape2)## melt function

## Data can be found in online resource 1
G <- read.csv("~/gs_data.csv", header = TRUE) ## Your route
G <- G[,-1]

## 2. COMPUTE INTERVARIANCE FOR L-MEASURES. 
## Extract max intervariance for N-lines, S-lines, E-lines, W-lines, ED-lines, EP-lines, WD-lines and WP-lines

## Create list with groups by phases
Interv_groups <- list("A1_int_group" = G[G$Phase %in% "A1",], "A2_int_group" = G[G$Phase %in% "A2",],
                      "B1_int_group" = G[G$Phase %in% "B1",], "B2_int_group" = G[G$Phase %in% "B2",],
                      "B3_int_group" = G[G$Phase %in% "B3",], "C_int_group" = G[G$Phase %in% "C",],
                      "Or_Ia_int_group" = G[G$Phase %in% "Or Ia",], "Or_Ib_int_group" = G[G$Phase %in% "Or Ib",])

## Organise values
Ig_full <- lapply(Interv_groups,select_full_lengths) ## Group only full length L-lines
Ig_half <- lapply(Interv_groups,select_half_lengths) ## Group only half length L-lines

## Set nas equal to column mean 
## NAs have been given for incomplete lines (microlith is not complete on that L-measure) for efficiency
## Use Huber M-estimator for means (see function)
Ig_full <- na_to_mean(Ig_full)
Ig_half <- na_to_mean(Ig_half)

## Select only columns that have 5 observations or more
Ig_full <- five_or_more(Ig_full)
Ig_half <- five_or_more(Ig_half)

## Intervariance for full length L-measures
## Select only L-measures that are present in all phases
Ig_full[[1]] <- Ig_full[[1]][,-c(11,24,25,29:32)]
Ig_full[[2]] <- Ig_full[[2]][,-c(11,22,25,26,30:32)]
Ig_full[[3]] <- Ig_full[[3]][,-c(11,24,28:30)]
Ig_full[[4]] <- Ig_full[[4]][,-c(11,12,23,24,27,31:33)]
Ig_full[[5]] <- Ig_full[[5]][,-c(26)]
Ig_full[[6]] <- Ig_full[[6]][,-c(11,12,23)]
Ig_full[[7]] <- Ig_full[[7]][,-c(11,22:24,27:29,33:36)]
Ig_full[[8]] <- Ig_full[[8]][,-c(11,12,23:26,29:31,35:39)]

## Compute Intervariance
Interv_Ig_full <- intervariance(Ig_full)

## Intervariance for half length L-measures
## Select only L-measures that are present in all phases
Ig_half[[1]] <- Ig_half[[1]][,-c(4,5,9,10,14:17,21:23)]
Ig_half[[2]] <- Ig_half[[2]][,-c(4,5,9,10,14:16,20:22)]
Ig_half[[3]] <- Ig_half[[3]][,-c(4,8,12:14,18:20)]
Ig_half[[4]] <- Ig_half[[4]][,-c(4,5,9,13:15,19:21)]
Ig_half[[5]] <- Ig_half[[5]][,-c(10,14)]
Ig_half[[7]] <- Ig_half[[7]][,-c(4:6,10:12,16:19,23:26)]
Ig_half[[8]] <- Ig_half[[8]][,-c(4:6,10:12,16:19,23:27)]

## Compute Intervariance
Interv_Ig_half <- intervariance(Ig_half)
length(Interv_Ig_full)
## 3. PCAs
Phase<-G$Phase #For visualisation

## 4-dimension groups for general measures
Gmea<-data.frame(G$MaxThick,G$Area,G$Long,G$Width)
colnames(Gmea)<-c("Max Thick","Area","Length","Width")

Gmea2<-data.frame(G$A1,G$A2,G$M_Wi,G$L_Le)
colnames(Gmea2)<-c("Angle 1","Angle 2","L-Width","L-Length")

## 4-dimension groups for intervariance L-measures
## NAs have been given for incomplete lines (microlith is not complete on that L-measure) for efficiency
L_full <- data.frame(G$N4,G$S4,G$E2,G$W3)
colnames(L_full)<-c("N4","S4","E2","W3")

L_half <- data.frame(G$ED2,G$EP3,G$WD3,G$WP2)
colnames(L_half)<-c("ED2","EP3","WD3","WP2")

##General measures
Gmea<-scale(Gmea)
Gmea.pca <- prcomp(Gmea,
                   center = TRUE,
                   scale = TRUE) 

plot(Gmea.pca, type="l")

Gmeap <- ggbiplot(Gmea.pca, obs.scale = 1, var.scale = 1, 
                  groups = Phase, ellipse = TRUE, 
                  circle = TRUE)
Gmeap <- Gmeap + scale_color_discrete(name = '')
Gmeap <- Gmeap + theme(legend.direction = 'horizontal', 
                       legend.position = 'top')
print(Gmeap)

##General measures 2
mg2<-mean(Gmea2[,4], na.rm = TRUE)## Non robust mean. Gmea2 variables present less outliers
Gmea2[is.na(Gmea2)]<-mg2
Gmea2<-scale(Gmea2)
Gmea2.pca <- prcomp(Gmea2,
                    center = TRUE,
                    scale = TRUE) 

plot(Gmea2.pca, type="l")

Gmea2p <- ggbiplot(Gmea2.pca, obs.scale = 1, var.scale = 1, 
                   groups = Phase, ellipse = TRUE, 
                   circle = TRUE)
Gmea2p <- Gmea2p + scale_color_discrete(name = '')
Gmea2p <- Gmea2p + theme(legend.direction = 'horizontal', 
                         legend.position = 'top')
print(Gmea2p)

## Intervariance full measures
L_full<-apply(L_full,2,na_to_mean_vec)
L_full<-scale(L_full)
L_full.pca <- prcomp(L_full,
                     center = TRUE,
                     scale = TRUE) 

plot(L_full.pca, type="l")

L_fullp <- ggbiplot(L_full.pca, obs.scale = 1, var.scale = 1, 
                    groups = Phase, ellipse = TRUE, 
                    circle = TRUE)
L_fullp <- L_fullp + scale_color_discrete(name = '')
L_fullp <- L_fullp + theme(legend.direction = 'horizontal', 
                           legend.position = 'top')
print(L_fullp)

## Intervariance half measures
L_half<-apply(L_half,2,na_to_mean_vec)
L_half<-scale(L_half)
L_half.pca <- prcomp(L_half,
                     center = TRUE,
                     scale = TRUE) 

plot(L_half.pca, type="l")

L_halfp <- ggbiplot(L_half.pca, obs.scale = 1, var.scale = 1, 
                    groups = Phase, ellipse = TRUE, 
                    circle = TRUE)
L_halfp <- L_halfp + scale_color_discrete(name = '')
L_halfp <- L_halfp + theme(legend.direction = 'horizontal', 
                           legend.position = 'top')
print(L_halfp)



### 4. JACCARD ANALYSIS
## Prepare data
N4 <- c()
for (i in 1:length(Ig_full)){
  current <- Ig_full[[i]]$N4
  N4 <- append(N4,current)
}

S4 <- c()
for (i in 1:length(Ig_full)){
  current <- Ig_full[[i]]$S5
  S4 <- append(S4,current)
}

E2<- c()
for (i in 1:length(Ig_full)){
  current <- Ig_full[[i]]$E2
  E2 <- append(E2,current)
}

W3 <- c()
for (i in 1:length(Ig_full)){
  current <- Ig_full[[i]]$W3
  W3 <- append(W3,current)
}

GJac<-data.frame(G$Phase,G$R1,G$MaxThick,G$Area,G$A1,G$A2,N4,S4,E2,W3)
GJac<-as.data.frame(GJac)
Phase<-as.factor(G$Phase)
colnames(GJac)<-c("Phase","R1","MaxThick","Area","A1","A2","N4","S4","E2","W3")

## Use mean values for easier visualisation
## For A1
HA1<-subset(GJac,Phase == "A1")
HA1$Phase <- NULL
MeA1<-c()
for (i in 1:9){
  me<-mean(HA1[,i])
  MeA1<-append(MeA1,me)
}

## For A2
HA2<-subset(GJac,Phase == "A2")
HA2$Phase <- NULL
MeA2<-c()
for (i in 1:9){
  me<-mean(HA2[,i])
  MeA2<-append(MeA2,me)
}

## For B1
HB1<-subset(GJac,Phase == "B1")
HB1$Phase <- NULL
MeB1<-c()
for (i in 1:9){
  me<-mean(HB1[,i])
  MeB1<-append(MeB1,me)
}

## For B2
HB2<-subset(GJac,Phase == "B2")
HB2$Phase <- NULL
MeB2<-c()
for (i in 1:9){
  me<-mean(HB2[,i])
  MeB2<-append(MeB2,me)
}

## For B3
HB3<-subset(GJac,Phase == "B3")
HB3$Phase <- NULL
MeB3<-c()
for (i in 1:9){
  me<-mean(HB3[,i])
  MeB3<-append(MeB3,me)
}

## For C
HC<-subset(GJac,Phase == "C")
HC$Phase <- NULL
MeC<-c()
for (i in 1:9){
  me<-mean(HC[,i])
  MeC<-append(MeC,me)
}

## For OrIa
HOrIa<-subset(GJac,Phase == "Or Ia")
HOrIa$Phase <- NULL
MeOrIa<-c()
for (i in 1:9){
  me<-mean(HOrIa[,i])
  MeOrIa<-append(MeOrIa,me)
}

## For OrIb
HOrIb<-subset(GJac,Phase == "Or Ib")
HOrIb$Phase <- NULL
MeOrIb<-c()
for (i in 1:9){
  me<-mean(HOrIb[,i])
  MeOrIb<-append(MeOrIb,me)
}


Medph<-data.frame(MeA1,MeA2,MeB1,MeB2,MeB3,MeC,MeOrIa,MeOrIb)
Medph<-t(Medph)
colnames(Medph)<-c("R1","MaxThick","Area","A1","A2","N4","S5","E2","W3")
row.names(Medph)<-c("A1","A2","B1","B2","B3","C","OrIa","OrIb")
GJac

## Jaccard distance and plots
jac<-vegdist(Medph,method = "jaccard")
jach<-as.matrix(jac)

# Plot
djac<-as.dendrogram(hclust(jac))
plot(djac,main = "Sites clustered by Jaccard similarity",horiz = TRUE,ylab = "",
     edgePar = list(col="Dark blue",lwd=2:2))

## Jaccard heatmap
pheatmap(jach,color=brewer.pal(9,"Blues"))

### 5. DIACHRONIC VISUALISATION OF MEANS OF L-LINES
## N-lines (N1-N10)
G_A1nl <-Ig_full$A1_int_group[,c(1:10)]
G_A1nlm <- apply(G_A1nl,2,mean,na.rm = TRUE)
G_A2nl <-Ig_full$A2_int_group[,c(1:10)]
G_A2nlm <- apply(G_A2nl,2,mean,na.rm = TRUE)
G_B1nl <-Ig_full$B1_int_group[,c(1:10)]
G_B1nlm <- apply(G_B1nl,2,mean,na.rm = TRUE)
G_B2nl <-Ig_full$B2_int_group[,c(1:10)]
G_B2nlm <- apply(G_B2nl,2,mean,na.rm = TRUE)
G_B3nl <-Ig_full$B3_int_group[,c(1:10)]
G_B3nlm <- apply(G_B3nl,2,mean,na.rm = TRUE)
G_Cnl<-Ig_full$C_int_group[,c(1:10)]
G_Cnlm <- apply(G_Cnl,2,mean,na.rm = TRUE)
G_Or_Ianl <-Ig_full$Or_Ia_int_group[,c(1:10)]
G_Or_Ianlm <- apply(G_Or_Ianl,2,mean,na.rm = TRUE)
G_Or_Ibnl <-Ig_full$Or_Ib_int_group[,c(1:10)]
G_Or_Ibnlm <- apply(G_Or_Ibnl,2,mean,na.rm = TRUE)

## Data frame for plotting
N_LinesO<-as.data.frame(rbind(G_A1nlm,G_A2nlm,G_B1nlm,G_B2nlm,G_B3nlm,G_Or_Ianlm,G_Or_Ibnlm))
Names <- c("A1","A2","B1","B2","B3","Or_Ia","Or_Ib") ## Without phase C
N_LinesO$Names<-Names

N_LinesO_long<-melt(N_LinesO,id="Names")
NLO<-ggplot(data=N_LinesO_long,aes(x=Names,y=value, color = variable,
                                   group = variable)) +
  ylab("Mean") + xlab("Phases") + labs(color="N_Lines") +
  geom_line(size=2, linejoin = "round") 
NLO

## S-lines (S1-S10)
G_A1sl <-Ig_full$A1_int_group[,c(11:20)]
G_A1slm <- apply(G_A1sl,2,mean,na.rm = TRUE)
G_A2sl <-Ig_full$A2_int_group[,c(11:20)]
G_A2slm <- apply(G_A2sl,2,mean,na.rm = TRUE)
G_B1sl <-Ig_full$B1_int_group[,c(11:20)]
G_B1slm <- apply(G_B1sl,2,mean,na.rm = TRUE)
G_B2sl <-Ig_full$B2_int_group[,c(11:20)]
G_B2slm <- apply(G_B2sl,2,mean,na.rm = TRUE)
G_B3sl <-Ig_full$B3_int_group[,c(11:20)]
G_B3slm <- apply(G_B3sl,2,mean,na.rm = TRUE)
G_Csl<-Ig_full$C_int_group[,c(11:20)]
G_Cslm <- apply(G_Csl,2,mean,na.rm = TRUE)
G_Or_Iasl <-Ig_full$Or_Ia_int_group[,c(11:20)]
G_Or_Iaslm <- apply(G_Or_Iasl,2,mean,na.rm = TRUE)
G_Or_Ibsl <-Ig_full$Or_Ib_int_group[,c(11:20)]
G_Or_Ibslm <- apply(G_Or_Ibsl,2,mean,na.rm = TRUE)

## Data frame for plotting
S_LinesO<-as.data.frame(rbind(G_A1slm,G_A2slm,G_B1slm,G_B2slm,G_B3slm,G_Or_Iaslm,G_Or_Ibslm))
Names <- c("A1","A2","B1","B2","B3","Or_Ia","Or_Ib") ## Without phase C
S_LinesO$Names<-Names

S_LinesO_long<-melt(S_LinesO,id="Names")
SlO<-ggplot(data=S_LinesO_long,aes(x=Names,y=value, color = variable,
                                   group = variable)) +
  ylab("Mean") + xlab("Phases") + labs(color="S_Lines") +
  geom_line(size=2, linejoin = "round") 
SlO

## E-lines (E1-E2)
G_A1el <-Ig_full$A1_int_group[,c(21:22)]
G_A1elm <- apply(G_A1el,2,mean,na.rm = TRUE)
G_A2el <-Ig_full$A2_int_group[,c(21:22)]
G_A2elm <- apply(G_A2el,2,mean,na.rm = TRUE)
G_B1el <-Ig_full$B1_int_group[,c(21:22)]
G_B1elm <- apply(G_B1el,2,mean,na.rm = TRUE)
G_B2el <-Ig_full$B2_int_group[,c(21:22)]
G_B2elm <- apply(G_B2el,2,mean,na.rm = TRUE)
G_B3el <-Ig_full$B3_int_group[,c(21:22)]
G_B3elm <- apply(G_B3el,2,mean,na.rm = TRUE)
G_Cel<-Ig_full$C_int_group[,c(21:22)]
G_Celm <- apply(G_Cel,2,mean,na.rm = TRUE)
G_Or_Iael <-Ig_full$Or_Ia_int_group[,c(21:22)]
G_Or_Iaelm <- apply(G_Or_Iael,2,mean,na.rm = TRUE)
G_Or_Ibel <-Ig_full$Or_Ib_int_group[,c(21:22)]
G_Or_Ibelm <- apply(G_Or_Ibel,2,mean,na.rm = TRUE)

## Data frame for plotting
E_LinesO<-as.data.frame(rbind(G_A1elm,G_A2elm,G_B1elm,G_B2elm,G_B3elm,G_Or_Iaelm,G_Or_Ibelm))
Names <- c("A1","A2","B1","B2","B3","Or_Ia","Or_Ib") ## Without phase C
E_LinesO$Names<-Names

E_LinesO_long<-melt(E_LinesO,id="Names")
ElO<-ggplot(data=E_LinesO_long,aes(x=Names,y=value, color = variable,
                                   group = variable)) +
  ylab("Mean") + xlab("Phases") + labs(color="E_Lines") +
  geom_line(size=2, linejoin = "round") 
ElO


## W-lines (W1-W3)
G_A1wl <-Ig_full$A1_int_group[,c(23:25)]
G_A1wlm <- apply(G_A1wl,2,mean,na.rm = TRUE)
G_A2wl <-Ig_full$A2_int_group[,c(23:25)]
G_A2wlm <- apply(G_A2wl,2,mean,na.rm = TRUE)
G_B1wl <-Ig_full$B1_int_group[,c(23:25)]
G_B1wlm <- apply(G_B1wl,2,mean,na.rm = TRUE)
G_B2wl <-Ig_full$B2_int_group[,c(23:25)]
G_B2wlm <- apply(G_B2wl,2,mean,na.rm = TRUE)
G_B3wl <-Ig_full$B3_int_group[,c(23:25)]
G_B3wlm <- apply(G_B3wl,2,mean,na.rm = TRUE)
G_Cwl<-Ig_full$C_int_group[,c(23:25)]
G_Cwlm <- apply(G_Cwl,2,mean,na.rm = TRUE)
G_Or_Iawl <-Ig_full$Or_Ia_int_group[,c(23:25)]
G_Or_Iawlm <- apply(G_Or_Iawl,2,mean,na.rm = TRUE)
G_Or_Ibwl <-Ig_full$Or_Ib_int_group[,c(23:25)]
G_Or_Ibwlm <- apply(G_Or_Ibwl,2,mean,na.rm = TRUE)

## Data frame for plotting
W_LinesO<-as.data.frame(rbind(G_A1wlm,G_A2wlm,G_B1wlm,G_B2wlm,G_B3wlm,G_Or_Iawlm,G_Or_Ibwlm))
Names <- c("A1","A2","B1","B2","B3","Or_Ia","Or_Ib") ## Witout phase C
W_LinesO$Names<-Names

W_LinesO_long<-melt(W_LinesO,id="Names")
WlO<-ggplot(data=W_LinesO_long,aes(x=Names,y=value, color = variable,
                                   group = variable)) +
  ylab("Mean") + xlab("Phases") + labs(color="W_Lines") +
  geom_line(size=2, linejoin = "round") 
WlO







