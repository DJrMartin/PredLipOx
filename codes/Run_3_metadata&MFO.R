# rm(list=ls())

rep = "raw_data/TXT_V3/"
dir(rep)
calcul <- function(x, duree, width, applied_function="M"){
  D_lip=D_glu=NULL
  Time_corrected <- x$TEMPS-min(x$TEMPS)
  if(applied_function == "M"){
    for (i in seq(0,duree, by=width)){
      fenetre <- which(Time_corrected>=i&Time_corrected<i+width)
      D_lip=c(D_lip, mean((1.695*(x$VO2[fenetre]/1000)- 1.701*(x$VCO2[fenetre]/1000))))
      D_glu=c(D_glu, mean((4.210*(x$VCO2[fenetre]/1000)- 2.962*(x$VO2[fenetre]/1000))))
    }
  }
  if(applied_function == "SD"){
    for (i in seq(0,duree, by=width)){
      fenetre <- which(Time_corrected>=i&Time_corrected<i+width)
      D_lip=c(D_lip, sd((1.695*(x$VO2[fenetre]/1000)- 1.701*(x$VCO2[fenetre]/1000))))
      D_glu=c(D_glu, sd((4.210*(x$VCO2[fenetre]/1000)- 2.962*(x$VO2[fenetre]/1000))))
    }
  }
  return(list(D_lip, D_glu))
}

CHO_mean =CHO_sd =NULL
Lipids_mean = Lipids_sd = Puissance = NULL
for (i in 1:52){
  META <- read.table(paste0(rep, dir(rep)[i]), dec=",", sep="\t")
  
  time=as.character(META$V1)
  META=data.frame(t(META))
  parameters=as.character(META$X1)
  
  exercice = which(time=='Départ Exercice')
  time=time[(exercice+1):dim(META)[2]]
  
  time=as.numeric(gsub(":", ".", gsub("\\.", "", time)))
  NA_omit=which(is.na(time))
  time=time[-NA_omit]
  
  Meta=META[,(exercice+1):dim(META)[2]]
  Meta=data.frame(t(Meta[,-NA_omit]))
  colnames(Meta)=parameters
  Meta=apply(Meta, 2, function(x) as.numeric(gsub(",", ".", gsub("\\.", "", x))))
  Meta=as.data.frame(Meta)
  Meta$TEMPS=time
  
  Puissance <- cbind(Puissance, Meta$PUIS.)
  Lipids_mean <- cbind(Lipids_mean, calcul(Meta, 19, 2)[[1]])
  Lipids_sd <- cbind(Lipids_sd, calcul(Meta, 19, 2, "SD")[[1]])
  CHO_mean <- cbind(CHO_mean, calcul(Meta, 19, 2)[[2]])
  CHO_sd <- cbind(CHO_sd, calcul(Meta, 19, 2, "SD")[[2]])
  
}

colnames(Lipids_mean) = colnames(Lipids_sd) = 
  colnames(CHO_mean) = colnames(CHO_sd) = as.character(substr(dir(rep), 1,4))

Lipids_mean[which(Lipids_mean<0)]=0
CHO_mean[which(CHO_mean<0)]=0

Lipids_mean[which(is.na(Lipids_mean))]=0
CHO_mean[which(is.na(CHO_mean))]=0

Lipids_mean <- as.matrix(Lipids_mean)
ID <- data.frame(CALIBRATION$ID, CALIBRATION$Groupe)
rownames(ID) = CALIBRATION$ID

df_spectro <- (merge(ID,t(Lipids_mean),  by="row.names", all=T))
df_spectro_sd <- (merge(ID,t(Lipids_sd),  by="row.names", all=T))
col.group = c("#999933","#AA4466","#88CCEE")

mean_bike <- apply(df_spectro[which(df_spectro$CALIBRATION.Groupe=="BIKE"),-c(1:3)], 2, function(x) mean(x[is.na(x)==FALSE]))
sd_bike <- apply(df_spectro_sd[which(df_spectro_sd$CALIBRATION.Groupe=="BIKE"),-c(1:3)], 2, function(x) mean(x[is.na(x)==FALSE]))

## PLOT FOR BIKE GROUP.
svg("figures/figure1.svg", width = 7, height = 4)
par(mar = c(5, 5, 2, 2))
plot(mean_bike, type='b', ylim=c(0,1.5), xlab='Duration of the submaximal test (minutes)', 
     ylab = expression("Maximal Fat Oxidation (g." * min^{-1} * ")"), col=col.group[3], pch=15, axes=F)
axis(1, at=seq(0,10, by=1), labels = seq(0,20, by=2))
axis(2)
segments(x0 = 1:10, y0=mean_bike-sd_bike, y1=mean_bike+sd_bike, col=col.group[3])

## POINTS FOR FOOT GROUP.
mean_foot <- apply(df_spectro[which(df_spectro$CALIBRATION.Groupe=="FOOT"),-c(1:3)], 2, function(x) mean(x[is.na(x)==FALSE]))
sd_foot <- apply(df_spectro_sd[which(df_spectro_sd$CALIBRATION.Groupe=="FOOT"),-c(1:3)], 2, function(x) mean(x[is.na(x)==FALSE]))
points(mean_foot,type='b', ylim=c(0,1.5), col=col.group[2], pch=16)
segments(x0 = 1:10, y0=mean_foot-sd_foot, y1=mean_foot+sd_foot, col=col.group[2])

## POINTS FOR CTL GROUP.
mean_ctl <- apply(df_spectro[which(df_spectro$CALIBRATION.Groupe=="CTL"),-c(1:3)], 2, function(x) mean(x[is.na(x)==FALSE]))
sd_ctl <- apply(df_spectro_sd[which(df_spectro_sd$CALIBRATION.Groupe=="CTL"),-c(1:3)], 2, function(x) mean(x[is.na(x)==FALSE]))
points(mean_ctl,type='b', ylim=c(0,1.5), col=col.group[1], pch=17)
segments(x0 = 1:10, y0=mean_ctl-sd_ctl, y1=mean_ctl+sd_ctl, col=col.group[1])

## EXPERIMENTAL DESIGN
abline(v=5, lty=2)
text(x=1.2, y=1.4, "VT1", cex=0.8)
text(x=6.2, y=1.4, "90% VT2", cex=0.8)
legend("topright", col= col.group[c(3,2,1)], pch=c(15,16,17), 
       legend = c("Cyclists", "Soccer players", "Non Athletes"), bty='n', cex=0.8)
abline(v=c(3.5, 4.5), lty=2, col="#DD0000")
text(x=4, y=1.4, "Predicting \n MFO", cex=0.8, srt=90, xpd=NA, col="#DD0000")
dev.off()

###### TABLE #############
df.clin <- read.table("data/meta_data.csv", sep=";", row.names = 1, dec=".", header=T)
rownames(df.clin) <- substr(rownames(df.clin), 4, 7)
df.clinical <- df.clin[df_spectro$Row.names[-which(is.na(df_spectro$CALIBRATION.ID))],-c(1,7,9,10,12,14,20,21,25,26)]
df.clinical$Lipides_SV1 = df.clinical$Lipides_SV1/9
df.clinical$Lipides_SV2 = df.clinical$Lipides_SV2/9
df.clinical$Glucides_SV1 = df.clinical$Glucides_SV1/4
df.clinical$Glucides_SV2 = df.clinical$Glucides_SV2/4

clinical.data=NULL
for (i in c("CTL", "FOOT", "BIKE")){
  w.g <- which(df.clin[df_spectro$Row.names[-which(is.na(df_spectro$CALIBRATION.ID))],]$Groupe==i)
  clinical.data<- cbind(clinical.data,paste0(round(colMeans(df.clinical[w.g,]),2), 
                                             " (+/- ",round(apply(df.clinical[w.g,], 2, sd),2),")"))
  
}
rownames(clinical.data) = colnames(df.clinical)
colnames(clinical.data) = c("CTL", "FOOT", "BIKE")

