rm(list=ls())
library(prospectr)
library(fda)

setwd(dir="you_path")

MIR.1 = read.csv("sp_brutes.csv", sep=';')
meta_data <- read.csv("meta_data.csv",sep=';', header=T)

# Data preparation
meta=data.frame(ID=substr(meta_data$ID,4,7),Y_lipides=as.numeric(meta_data$Lipides_SV1),
                Y_glucides=as.numeric(meta_data$Glucides_SV1),Y_ratio=as.numeric(meta_data$R_SV1), 
                Groupe=meta_data$Groupe, Time=meta_data$Date_sampling_IR)
meta$Y_lipides[which(meta$Y_lipides<0)]=meta$Y_ratio[which(meta$Y_ratio<0)]=0
meta$Time=as.numeric(substr(meta$Time, 4,5))*30+as.numeric(substr(meta$Time, 1,2))
meta$Groupe=factor(meta$Groupe,levels=c("CTL","FOOT","BIKE"))

# Extraction of the intrested wanenumbers
freq=as.numeric(substr(colnames(MIR.1[,-1]),2,9))

fingerprint=c(which(freq<1800&freq>800))

matplot(t(MIR.1[,fingerprint]), typ="l",lty=1,ylab="Absorbances",xlab="freq", xaxt="n")
axis(1,at=seq(0,dim(MIR.1[,fingerprint])[2], by=dim(MIR.1[,fingerprint])[2]/10 ), labels=as.character(seq(1800,800, by=-100 )))

spectra_1=data.frame(apply(MIR.1[,fingerprint], 2, as.numeric))
rownames(spectra_1)=MIR.1$ID

# First Outliers
res.PCA=FactoMineR::PCA(MIR.1[,c(which(freq<3000&freq>2900), which(freq<1600&freq>1500))], graph = T)
outliers=c(which(as.numeric(res.PCA$ind$coord[,1])>20), which(as.numeric(res.PCA$ind$coord[,1])<(-20)))

corrected=cbind(spectra_1[-outliers,])

matplot(t(corrected), typ="l",lty=1,ylab="Absorbances",xlab="freq", xaxt='n')
axis(1,at=seq(0,dim(MIR.1[,fingerprint])[2], by=dim(MIR.1[,fingerprint])[2]/14 ), labels=c(as.character(seq(3200,800, by=-170))))

AV=corrected[grep("AV",rownames(corrected)),]
AP=corrected[grep("AP",rownames(corrected)),]

rownames(AV)=rownames(corrected)[grep("AV",rownames(corrected))]
rownames(AP)=rownames(corrected)[grep("AP",rownames(corrected))]

# Normalisation
AV_FOR_AN=AP_FOR_AN=NULL
for (i in unique(substr(rownames(corrected), 3, 6))){
  new_data=apply(AV[which(substr(rownames(AV), 3, 6)==i), ], 2, mean)
  AV_FOR_AN=rbind(AV_FOR_AN,new_data)
  new_data=apply(AP[which(substr(rownames(AP), 3, 6)==i), ], 2, mean)
  AP_FOR_AN=rbind(AP_FOR_AN,new_data)
}
matplot(t(AP_FOR_AN), typ="l",lty=1,ylab="Absorbances",xlab="freq", xaxt='n')
axis(1,at=seq(0,dim(MIR.1[,fingerprint])[2], by=dim(MIR.1[,fingerprint])[2]/14 ), labels=c(as.character(seq(3200,800, by=-170))))

corrected_3=EMSC::EMSC(rbind(AV_FOR_AN,AP_FOR_AN))
w=which(FactoMineR::PCA(corrected_3$corrected)$ind$coord[,1]>50)
matplot(t(corrected_3$corrected[-w,]), typ="l",lty=1,ylab="Absorbances",xlab="freq", xaxt='n')
axis(1,at=seq(0,dim(MIR.1[,fingerprint])[2], by=dim(MIR.1[,fingerprint])[2]/14 ), labels=c(as.character(seq(3200,800, by=-170))))

corrected_4=corrected_3$corrected[-w,]

projRecomp = function(xdata, nBasis, t = 1:dim(xdata)[2], basis = "Splines"){
  t = sort((t-min(t))/(max(t)-min(t)))
  if (basis == "Fourier") {basisobj = create.fourier.basis(nbasis = nBasis)}
  if (basis == "Splines") {basisobj = create.bspline.basis(norder = 3, breaks = seq(head(t,1), tail(t,1), length = nBasis-1))}
  BFunction = getbasismatrix(t, basisobj) 
  Fdata = t(sapply(1:nrow(xdata), 
                   FUN = function(i) t(solve(t(BFunction)
                                             %*% BFunction)%*% t(BFunction) %*% xdata[i,])))
  FdataRec = t(sapply(1:nrow(xdata), 
                      FUN = function(i) t(BFunction %*% solve(t(BFunction)%*% BFunction)%*% 
                                            t(BFunction) %*% xdata[i,])))
  return(list(coeffProj = Fdata, foncRecon = FdataRec, BFunction = BFunction, basisobj = basisobj))
}

# Spline projections
b=c(7,5,3)
p = dim(corrected_4)[2]
VI = matrix(0,length(b),p)
x_all = NULL
cnt = 1
L = 0
x_recom=NULL
for (i in b){
  nBasis = round(p/i)
  L = L + nBasis
  x = projRecomp(as.matrix(corrected_4), nBasis = nBasis)
  print(nBasis)
  if (cnt==1){
    x_all = x$coeffProj
    x_recom= x$foncRecon
  } else {
    x_all = cbind(x_all,x$coeffProj)
    x_recom= cbind(x_recom, x$foncRecon)
  }
  cnt = cnt+1
}

matplot(t(x_recom), typ="l",lty=1,ylab="Absorbances",xlab="freq", xaxt='n')
x_recom = data.frame(x_recom)
splines_ID = c(paste0(unique(substr(rownames(corrected), 3, 6)), "_AV"), paste0(unique(substr(rownames(corrected), 3, 6)), "_AP"))[-w]
x_recom$ID = substr(splines_ID, 1, 4)

AV_FOR_ANALYSIS=merge(meta,x_recom[grep("AV", splines_ID),],  by='ID', all=F)
AP_FOR_ANALYSIS=merge(meta,x_recom[grep("AP", splines_ID),],  by='ID', all=F)

