library(caret)
library(pROC)
library(randomForest)

x_all_1 = data.frame(apply(x_all, 2, as.numeric))
x_all_1$ID = substr(splines_ID, 1, 4)

diff <- merge(x_all_1[grep("AP", splines_ID),],x_all_1[grep("AV", splines_ID),],  by='ID', all=F)
Delta <- diff[,grep("x",colnames(diff))]/diff[,grep("y",colnames(diff))]
Delta$ID = diff[,1]

Y=data.frame(ID=substr(meta_data$ID,4,7), Y=meta_data$Lipides_SV1, Groupe=meta_data$Groupe, Time=meta_data$Date_sampling_IR)
CALIBRATION=data.frame(merge(meta,x_all_1[grep("AV", splines_ID),],  by='ID', all=F))

CALIBRATION$Y_lipides=CALIBRATION$Y_lipides/9

#========== Manual Search ====================
# control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
# tunegrid <- expand.grid(.mtry=c(10,50, 100, 200))
# modellist <- list()
# metric='Accuracy'
# for (maxnodes in c(2,5,15)) {
#   set.seed(122)
#   fit <- train(Y_lipides~.,data = CALIBRATION, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control, ntree=2000, maxnodes=maxnodes)
#   key <- toString(maxnodes)
#   modellist[[key]] <- fit
# }
# # compare results
# results <- resamples(modellist)
# summary(results)
# dotplot(results)

Y_TEST=PRED=IMP=T_test=NULL
set.seed(122)
#Classification test using RANDOM FOREST ALGORITHM
for (rep in (1:50)){
  #PARTITION
  inTraining <- caret::createDataPartition(CALIBRATION$Y_lipides, p=0.70, list = FALSE)
  training <- CALIBRATION[ inTraining,-c(1,3,4, 5, 6)]
  testing  <- CALIBRATION[-inTraining,-c(1,3,4, 5, 6)]
  #MODEL
  rf=randomForest::randomForest(Y_lipides~.,data = training, mtry=200,
                                ntree=2000, maxnodes=7)
  #PREDICTION ON VALIDATION SETS
  PRED = c(PRED,predict(rf, newdata = testing, type='response'))
  Y_TEST=c(Y_TEST,testing$Y_lipides)
  T_test = c(T_test, CALIBRATION$Time[-inTraining])
  
  IMP=cbind(IMP, rank(caret::varImp(rf)[,1]))
}


layout(matrix(c(1,2), nrow=1))
## Prediction en fonction de l'oxidation lipidique observée.
id.pred <- as.numeric(unique(rownames(as.matrix(PRED))))
col.group = as.character(factor(CALIBRATION$Groupe[id.pred], c('CTL',"FOOT","BIKE"),c("#99993390","#AA446690","#88CCEE90")))

boxplot(PRED~Y_TEST, col=col.group[order(CALIBRATION$Y_lipides[id.pred])], boxwex=0.2,
        axes=F, ylab='Predict lipid oxidation (g/min)',outlines=FALSE,xlim=c(-.2, 1.55),
        xlab="Observed lipid oxidation (g/min)", at=sort(unique(CALIBRATION$Y_lipides)))
abline(v=0.40, col="#55AE99")
abline(h=0.40, col="#55AE99")
axis(2)
axis(1)
legend("bottomright", fill= unique(col.group),
       legend = c("Cyclists", "Soccer", "Non Athl."), bty='n', cex=0.8)

ROC = pROC::roc(Y_TEST>0.40~PRED)
plot(ROC$specificities, ROC$sensitivities, lwd=2,
     type='l', col="#55AE99", xlab="Specificities", ylab="1-Sensitivities")
ROC = pROC::roc(Y_TEST>0.18~PRED)
points(ROC$specificities, ROC$sensitivities, lwd=2,
     type='l', col="tomato")
abline(1,-1)
legend("bottomleft", fill=c("#55AE99","tomato"), 
       legend=c(paste("AUROC of > 0.40 g/min =", round(pROC::roc(Y_TEST>0.40~PRED)$auc,3)),
                paste("AUROC of > 0.18 g/min =", round(pROC::roc(Y_TEST>0.18~PRED)$auc,3))), 
       bty="n", cex=0.7)


seuil <- which(Y_TEST<1)
sqrt(mean((PRED[seuil]-Y_TEST[seuil])^2))
caret::confusionMatrix(data=as.factor(PRED>0.4), 
                       reference = as.factor(Y_TEST>0.4))


sIMP = apply(IMP, 1, sum)
qimp = quantile(sIMP,.97)

freq[fingerprint]
L = length(fingerprint)
L7 = length(seq(4,L,by=7))
L5 = length(seq(3,L,by=5))
L3 = length(seq(2,L,by=3))
L5+L7+L3 == length(sIMP)

w7 = which(sIMP[1:L7]>qimp)
w5 = which(sIMP[L7+1:L5]>qimp)
w3 = which(sIMP[L7+L5+1:L3]>qimp)

freq_imp <- c(freq[fingerprint][seq(4,L,by=7)][w7],
  freq[fingerprint][seq(3,L,by=5)][w5],
  freq[fingerprint][seq(2,L,by=3)][w3])

layout(matrix(c(1,1,2,3,
                1,1,2,4), nrow = 2, byrow=T))
par(mar=c(5,4,5,2))
plot(apply(t(x_recom[grep("AV", splines_ID),][,1:L]),1,mean), 
     typ="l",lty=1,ylab="Absorbances",xlab=expression(paste("Wavenumbers (",cm^-1,")")), xaxt='n', axes=F, ylim=c(0.05, 0.2))
axis(1,at=seq(0,dim(MIR.1[,fingerprint])[2], by=dim(MIR.1[,fingerprint])[2]/7 ), 
     labels=c(as.character(round(seq(1800,800, by=-1000/7)))))
axis(2)

arrows(x0=(1:length(corrected_3$corrected[9,]))[seq(4,L,by=7)][w7],length=0.05, col="#55AE99",
       y0=rep(max(corrected_3$corrected),length(w7)), y1=rep(max(corrected_3$corrected)-0.015,length(w7)))
arrows(x0=(1:length(corrected_3$corrected[9,]))[seq(3,L,by=5)][w5],length=0.05, col="tomato",
       y0=rep(max(corrected_3$corrected),length(w7)), y1=rep(max(corrected_3$corrected)-0.01,length(w5)))
arrows(x0=(1:length(corrected_3$corrected[9,]))[seq(2,L,by=3)][w3],length=0.05, col="#AA4499",
       y0=rep(max(corrected_3$corrected),length(w3)), y1=rep(max(corrected_3$corrected)-0.005,length(w3)))
text(x=c((1:length(corrected_3$corrected[9,]))[seq(4,L,by=7)][w7], (1:length(corrected_3$corrected[9,]))[seq(3,L,by=5)][w5]),
     y=rep(0.2,3),srt=90,cex=0.9,pos=1,
  c(round(freq[fingerprint][seq(4,L,by=7)][w7]), round(freq[fingerprint][seq(3,L,by=5)][w5])))
legend("bottomright", legend = c('Spline of Basis 7','Spline of Basis 5', 'Spline of Basis 3'),
       fill = c("#55AE99","tomato", "#AA4499"), bty="n", cex=0.9, title = "Var. Imp.")

col.splines <- data.frame(IMP = rf$importance[c(w7,w5,w3)], 
                          colors = c(rep("#55AE99", length(w7)), rep("tomato",length(w5)) ,rep("#AA4499",length(w3))),
                          FREQ = c(freq[fingerprint][seq(4,L,by=7)][w7],freq[fingerprint][seq(3,L,by=5)][w5],freq[fingerprint][seq(2,L,by=3)][w3]))
order.imp <- order(col.splines$IMP)
barplot(col.splines$IMP[order.imp], horiz = T, xlab='Mean Decrease Gini', 
        col=col.splines$colors[order.imp] )
text(y=seq(0.6,12.5, by=10.9/9), x=rep(0.05, length(c(w3,w5,w7))), 
     labels = round(col.splines$FREQ)[order.imp])
par(mar=c(5,4,2,2))
plot(CALIBRATION$Y_lipides~CALIBRATION[,-c(1,3,4, 5, 6)][,which(sIMP>qimp)[1]], 
     pch=as.numeric(as.character(factor(CALIBRATION$Groupe,c("CTL", 'FOOT', 'BIKE'),c(17,16,15)))),
     xlab=paste(round(freq_imp[1]),"cm-1"),ylab="Lipid oxidation (g/min)",
     col=as.character(factor(CALIBRATION$Groupe,c("CTL", 'FOOT', 'BIKE'),c("#999933","#AA4466","#88CCEE"))))
legend("topleft", col= c("#88CCEE","#AA4466","#999933"), pch=c(15,16,17), 
       legend = c("Cyclists", "Soccer", "Non Athl."), bty='n', cex=0.8)
plot(CALIBRATION$Y_lipides~CALIBRATION[,-c(1,3,4, 5, 6)][,which(sIMP>qimp)[2]], 
     pch=as.numeric(as.character(factor(CALIBRATION$Groupe,c("CTL", 'FOOT', 'BIKE'),c(17,16,15)))),
     xlab=paste(round(freq_imp[2]),"cm-1"), ylab="Lipid oxidation (g/min)",
     col=as.character(factor(CALIBRATION$Groupe,c("CTL", 'FOOT', 'BIKE'),c("#999933","#AA4466","#88CCEE"))))
legend("topright", col= c("#88CCEE","#AA4466","#999933"), pch=c(15,16,17), 
       legend = c("Cyclists", "Soccer", "Non Athl."), bty='n', cex=0.8)

