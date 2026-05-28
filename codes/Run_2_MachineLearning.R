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

####################### ITERATION PROCEDURE #########################
Y_TEST = PRED = IMP = T_test = unique_PRED_DF = NULL
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
  PRED = c(PRED, predict(rf, newdata = testing, type='response'))
  Y_TEST = c(Y_TEST, testing$Y_lipides)
  T_test = c(T_test, CALIBRATION$Time[-inTraining])
  
  PRED_DF <- rep(NA, 48)
  PRED_DF[-inTraining] = predict(rf, newdata = testing, type='response')
  unique_PRED_DF = cbind(unique_PRED_DF, PRED_DF)
  
  IMP=cbind(IMP, rank(caret::varImp(rf)[,1]))
}

####################### LOO #########################
pred = NULL
for (rep in 1:48){
  training <- CALIBRATION[-rep,-c(1,3,4, 5, 6)]
  testing  <- CALIBRATION[rep,-c(1,3,4, 5, 6)]
  #MODEL
  rf = randomForest::randomForest(Y_lipides~., data = training, mtry = 200,
                                ntree = 2000, maxnodes = 7)
  #PREDICTION ON VALIDATION SETS
  pred = c(pred, predict(rf, newdata = testing, type='response'))
}
#####################################################
seuil <- which(CALIBRATION$Y_lipides<1)
sqrt(mean((pred[seuil]-CALIBRATION$Y_lipides[seuil])^2))
caret::confusionMatrix(data=as.factor(pred>0.35), 
                       reference = as.factor(CALIBRATION$Y_lipides>0.35))
pROC::roc(CALIBRATION$Y_lipides>0.35~pred)

svg("figures/figure2.svg", width = 12, height = 6)
layout(matrix(c(1,2), nrow=1))
## Prediction en fonction de l'oxidation lipidique observée.
id.pred <- as.numeric(unique(rownames(as.matrix(PRED))))
col.group = as.character(factor(CALIBRATION$Groupe[id.pred], c('CTL',"FOOT","BIKE"),c("#99993390","#AA446690","#88CCEE90")))
par(mar = c(5,5, 2, 2))
boxplot(PRED~Y_TEST, col = col.group[order(CALIBRATION$Y_lipides[id.pred])], boxwex = 0.2,
        axes = F, ylab = expression("Predicted MFO (g." * min^{-1} * ")"), outlines = FALSE,xlim = c(-.2, 1.55),
        xlab = expression("Observed MFO (g." * min^{-1} * ")"), at = sort(unique(CALIBRATION$Y_lipides)),
        ylim = c(-.2, 1.55))

legend("topright", "TP", bty = "n", cex = 0.8)
legend("bottomright", "FN", bty = "n", cex = 0.8)
legend("bottomleft", "TN", bty = "n", cex = 0.8)
legend("topleft", "FP", bty = "n", cex = 0.8)
abline(v=0.35, col="#55AE99")
abline(h=0.35, col="#55AE99")
axis(2)
axis(1)
abline(0, 1, lty="dashed")
legend("top", fill= unique(col.group),
       legend = c("Cyclists", "Soccer players", "Non Athletes"), bty='n', cex=0.8)

mean_pred <- apply(unique_PRED_DF, 1, function(x) mean(x[which(x!="NA")]))

ROC_T1 = pROC::roc(CALIBRATION$Y_lipides>0.35~mean_pred)
plot(ROC_T1$specificities, ROC_T1$sensitivities, lwd=2,
     type='l', col="#55AE99", xlab="Specificities", ylab="1-Sensitivities")

ROC_T2 = pROC::roc(CALIBRATION$Y_lipides>0.7~mean_pred)
points(ROC_T2$specificities, ROC_T2$sensitivities, lwd=2,
     type='l', col="tomato")

ROC_T3 = pROC::roc(CALIBRATION$Y_lipides>0.49~mean_pred)
points(ROC_T3$specificities, ROC_T3$sensitivities, lwd=2,
       type='l', col="purple")

abline(1,-1)
legend("bottomleft", fill=c("#55AE99","tomato", "purple"), 
       legend=c(paste("AUROC of < 0.35 g/min =", round(ROC_T1$auc,2)),
                paste("AUROC of > 0.7 g/min =", round(ROC_T2$auc,2)),
                paste("AUROC of > 0.5 g/min =", round(ROC_T3$auc,2))), 
       bty="n", cex=0.8)
dev.off()

seuil <- which(CALIBRATION$Y_lipides<1)
sqrt(mean((mean_pred[seuil]-CALIBRATION$Y_lipides[seuil])^2))
caret::confusionMatrix(data=as.factor(mean_pred>0.35), 
                       reference = as.factor(CALIBRATION$Y_lipides>0.35))

sIMP = apply(IMP, 1, sum)
qimp = quantile(sIMP,.97)

freq[fingerprint]
L = length(fingerprint)

L7 = length(seq(4,L,by=7))
L5 = length(seq(3,L,by=5))
L3 = length(seq(2,L,by=3))
L5+L7+L3 == length(sIMP)
CS = cumsum(c(L7, L5, L3))

w7 = which(sIMP[1:L7]>qimp)
w5 = which(sIMP[(L7+1):CS[2]]>qimp)
w3 = which(sIMP[(L7+L5+1):CS[3]]>qimp)

freq_imp <- c(freq[fingerprint][seq(4,L,by=7)][w7],
  freq[fingerprint][seq(3,L,by=5)][w5],
  freq[fingerprint][seq(2,L,by=3)][w3])

svg("figures/figure3.svg", width = 12, height = 6)
layout(matrix(c(1,1,2,2,
                1,1,3,3), nrow = 2, byrow=T))
par(mar=c(5,4,5,2))
plot(apply(t(x_recom[grep("AV", splines_ID),][,1:L]),1,mean), 
     typ="l",lty=1,ylab="Absorbances",xlab = expression(paste("Wavenumbers (",cm^-1,")")), xaxt='n', axes=F, ylim=c(0.05, 0.2))
axis(1,at=seq(0,dim(MIR.1[,fingerprint])[2], by=dim(MIR.1[,fingerprint])[2]/7 ), 
     labels=c(as.character(round(seq(1800,800, by=-1000/7)))))
axis(2)

arrows(x0=(1:length(corrected_3$corrected[9,]))[seq(4,L,by=7)][w7],length=0.05, col="#55AE99",
       y0=rep(max(corrected_3$corrected),length(w7)), y1=rep(max(corrected_3$corrected)-0.015,length(w7)))
arrows(x0=(1:length(corrected_3$corrected[9,]))[seq(3,L,by=5)][w5],length=0.05, col="tomato",
       y0=rep(max(corrected_3$corrected),length(w5)), y1=rep(max(corrected_3$corrected)-0.01,length(w5)))
arrows(x0=(1:length(corrected_3$corrected[9,]))[seq(2,L,by=3)][w3],length=0.05, col="#AA4499",
       y0=rep(max(corrected_3$corrected),length(w3)), y1=rep(max(corrected_3$corrected)-0.005,length(w3)))
text(x=c((1:length(corrected_3$corrected[9,]))[seq(4,L,by=7)][w7], (1:length(corrected_3$corrected[9,]))[seq(3,L,by=5)][w5]),
     y=rep(0.2,3),srt=90,cex=1,pos=1,
  c(round(freq[fingerprint][seq(4,L,by=7)][w7]), round(freq[fingerprint][seq(3,L,by=5)][w5])))

legend("bottomright", legend = c('Spline of Basis 7','Spline of Basis 5', 'Spline of Basis 3'),
       fill = c("#55AE99","tomato", "#AA4499"), bty="n", cex=0.9, title = "Var. Imp.")

data.frame(
  "wavenumbers" = c(freq[fingerprint][seq(4,L,by=7)], freq[fingerprint][seq(3,L,by=5)], freq[fingerprint][seq(2,L,by=3)])[which(sIMP>qimp)],
  "imp" = rf$importance[which(sIMP>qimp)])

par(mar=c(5,5,2,2))
plot(CALIBRATION$Y_lipides ~ CALIBRATION[,-c(1,3,4, 5, 6)][,which(sIMP>qimp)[6]], 
     pch=as.numeric(as.character(factor(CALIBRATION$Groupe,c("CTL", 'FOOT', 'BIKE'),c(17,16,15)))),
     xlab = expression(paste("1518 ",cm^-1,)), 
     ylab = expression("Fat oxidation (g." * cm^{-1} * ")"),
     col = as.character(factor(CALIBRATION$Groupe,c("CTL", 'FOOT', 'BIKE'),c("#999933","#AA4466","#88CCEE"))))

legend("topleft", col= c("#88CCEE","#AA4466","#999933"), pch=c(15,16,17), 
       legend = c("Cyclists", "Soccer players", "Non Athletes"), bty='n', cex=1)

summary(model <- lm(scale(CALIBRATION$Y_lipides)~scale(CALIBRATION[,-c(1,3,4, 5, 6)][,which(sIMP>qimp)[6]])))
legend("left", legend="Pearson's r = 0.55, p<0.001", bty="n")

plot(CALIBRATION$Y_lipides~CALIBRATION[,-c(1,3,4, 5, 6)][,which(sIMP>qimp)[3]], 
     pch=as.numeric(as.character(factor(CALIBRATION$Groupe,c("CTL", 'FOOT', 'BIKE'),c(17,16,15)))),
     xlab=expression(paste("1166 ",cm^-1,)), ylab=expression("Fat oxidation (g." * cm^{-1} * ")"),
     col=as.character(factor(CALIBRATION$Groupe,c("CTL", 'FOOT', 'BIKE'),c("#999933","#AA4466","#88CCEE"))))

legend("topright", col= c("#88CCEE","#AA4466","#999933"), pch=c(15,16,17), 
       legend = c("Cyclists", "Soccer players", "Non Athletes"), bty='n', cex=1)

model <- lm(scale(CALIBRATION$Y_lipides)~scale(CALIBRATION[,-c(1,3,4, 5, 6)][,which(sIMP>qimp)[3]]))
summary(model)
legend("right", legend="Pearson's r = -0.41, p=0.003", bty="n")
dev.off()

