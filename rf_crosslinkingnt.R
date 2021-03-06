source("script/functions.R")
data.file <- "raw/featuretable_v1.txt"
cv.folds <- 10
repeat.times <- 1
mtry.seq <- seq(1,20)
ntree.seq <- c(10,100,200,500,1000,2000,5000)
plot <- TRUE
write <- TRUE
seed <- 30
permute.time <- 2000

feature.fill.color <- c("#8dd3c7", "#b3de69","#ffffb3", "#ae017e", "#dd3497", "#f768a1","#fa9fb5", "#fcc5c0","#fde0dd", "#fb8072", "#d7191c","#fdae61",  '#abd9e9',  "#2c7bb6")
aa.color.code <- c("#991122", "#ffddee", "#aa88bb", "#774499", "#3377bb",  "#ddeeff",  "#77aacc", "#225577", "#449944", "#bbddbb", "#224422", "#77bb77", "#999933",  "#dddd33", "#eeee88", "#ffeedd", "#bb6622", "#dd2222", "#ee9999", "#992222")

summary.file <- "result/sample_summary.txt"
aa.freq.file <- "result/aa_interaction_frequency.txt"
prediction.file <- "result/RF_prediction_result.txt"
aaprediction.file <-  "result/RF_AA_prediction_result.txt"
featurerankres.file <- "result/feature_importance_summary.txt"
aafeaturerankres.file <- "result/aa_feature_importance_summary.txt"
aafeatureimp.file <- "result/aa_feature_localimp.txt"
rf.tune.roc.mat.file <- "result/mtry_ntree_roc_mat.txt"



modelcomapareplot.file <- "plot/RF_onevstwolayer_roc_plot.pptx"
featurerankplot.file <- "plot/feature_importance_scatterplot.pptx"
aafeaturerankplot.file <- "plot/aa_feature_importance_scatterplot.pptx"
featurelocalimpplot.file <- "plot/feature_localimp_heatmap.pptx"
featurerankplot.stackwitharomatic.file <- "plot/feature_importance_scatterplot_stackwitharomatic.pptx"

### preprocess prediction and feature table ###
process.obj <- preProcessData(data.file, varcutoff = 0)
rawdata <- process.obj$rawData
rawdata.full <- extendrawdata(rawdata)
cleandata <- process.obj$processedData
interact.nt.idx <- process.obj$interact.nt.idx 
aa.code <- unlist(lapply(colnames(cleandata)[grep("interact_", colnames(cleandata))], function(x) return(unlist(strsplit(as.character(x), split="_"))[2])))
interaction.code <- c("interact", "hbond_po4.sidechain", "hbond_po4.backbone","hbond_sugar.sidechain",  "hbond_sugar.backbone", "hbond_base.sidechain", "hbond_base.backbone", "base_aa_pair", "base_aa_stack")
conform.code <- colnames(cleandata)[2:13]
xl.idx <- which(cleandata$group == "xl")

interact.nt.idx.out <- rep(0, dim(rawdata)[1])
interact.nt.idx.out[interact.nt.idx] <- 1
write.table(interact.nt.idx.out, file="result/final.keep.idx.txt", sep="\t", quote = F, col.names = F, row.names = F)

## check domain nt composition
domain.xl.nt.freq <- matrix(0, nrow = length(levels(as.factor(rawdata.full$domain))), ncol = 4)
rownames(domain.xl.nt.freq) <- levels(as.factor(rawdata.full$domain))
colnames(domain.xl.nt.freq) <- c("A", "C", "G", "U")
for (i in 1:length(levels(as.factor(rawdata.full$domain)))){
  print(i)
  domain.xl.nt.freq[i,] <- table(rawdata.full$nt[interact.nt.idx][xl.idx][which(rawdata.full$domain[interact.nt.idx][xl.idx] == levels(as.factor(rawdata.full$domain))[i])])
  
}
write.table(domain.xl.nt.freq, "result/doman_xl_nt_frequency.txt" , quote = F, sep="\t")


### summarize features ###

for (i in interaction.code){
  print(i)
  aa.freq <- aafreqCalculation(cleandata, write = T,aa.freq.file = paste("result/", i, ".freq.txt", sep=""), key = i, aa.code = aa.code)
}

conform.freq.mat <- NULL
for (i in 2:13){
  print(i)
  conform.freq.mat <- cbind(conform.freq.mat, conformfreqCalculation(cleandata, i))
}
conform.freq.mat <- rbind(conform.freq.mat, percentage = apply(conform.freq.mat, 2, sum)/rep(c(length(xl.idx), (dim(cleandata)[1] - length(xl.idx))), 12))

write.table(conform.freq.mat, "result/conformation_frequency.txt" , quote = F, sep="\t")


## check the upstream context and downstream context ## 

upnt.xl.count <-  table(rawdata.full$upstreamNt[process.obj$interact.nt.idx][xl.idx])
dnnt.xl.count <-  table(rawdata.full$dnstreamNt[process.obj$interact.nt.idx][xl.idx])

upnt.nonxl.count <-  table(rawdata.full$upstreamNt[process.obj$interact.nt.idx][-xl.idx])
dnnt.nonxl.count <-  table(rawdata.full$dnstreamNt[process.obj$interact.nt.idx][-xl.idx])

write.table(cbind(upnt.xl.count,dnnt.xl.count, upnt.nonxl.count, dnnt.nonxl.count), file = "result/updnstream.nt.count.txt", sep = "\t", quote = F, col.names = T, row.names = T)



### predict with randome forest ###


## add up/dn stream code
cleandata.bak <- cleandata
cleandata <- cbind(cleandata,rawdata.full[interact.nt.idx, c(259,262, 266) ] )
cleandata$upstreamNt <- as.factor(cleandata$upstreamNt)
cleandata$dnstreamNt <- as.factor(cleandata$dnstreamNt)
cleandata$nt <- as.factor(cleandata$nt)


numericl.feature.mat <- cbind(model.matrix(~ nt -1 , cleandata),
                              model.matrix(~ upstreamNt-1 , cleandata),
                              model.matrix(~ dnstreamNt-1 , cleandata))
#cleandata.bak <- cleandata
cleandata <- cbind(cleandata[,1:180],as.data.frame(numericl.feature.mat) )

set.seed(seed)
rf.obj <- rfPredict(cleandata, mtry.seq, ntree.seq, cv.folds, repeat.times, smote = T)

rf.obj.nosmote.stackwitharomatic <- rfPredict(cleandata[which(cleandata$base_aa_stack_aromatic!=0),], mtry.seq, ntree.seq, cv.folds, repeat.times, smote = F)

rf.par <-  rf.obj$fit$bestTune
rf.stackwitharomatic.par <- rf.obj.nosmote.stackwitharomatic$fit$bestTune
### permute features to get feature rank pvalue ###
set.seed(seed)
rf.permute <- rf.permutation(cleandata, SMOTE = T, permute.time = permute.time, rf.par = rf.par)
rf.stackwitharomatic.permute <- rf.permutation(cleandata[which(cleandata$base_aa_stack_aromatic!=0),], SMOTE = F, permute.time = permute.time, rf.par = rf.stackwitharomatic.par)

### calculate GLM weight ###  
set.seed(seed)
glmGrid <- expand.grid(alpha=seq(0,1,0.1), lambda = seq(0,1, length =50))
glm.fit.obj <- calculateGLMweight(cbind(group=cleandata[,1], as.data.frame(scale(cleandata[, -1]))), cv.folds,repeat.times = repeat.times, metric = "AUC")
glm.fit.obj.stackwitharomatic <- calculateGLMweight(cbind(group=cleandata[,1], as.data.frame(scale(cleandata[, -1])))[which(cleandata$base_aa_stack_aromatic!=0),], cv.folds, repeat.times = repeat.times, SMOTE = F, metric = "AUC", tuneGrid= glmGrid)

## remove stacking features
set.seed(seed)
rf.obj.rm.stacking <- rfPredict(cleandata[,-c(167:180)], mtry.seq, ntree.seq, cv.folds, repeat.times, smote = T)
rf.rm.stacking.par <- rf.obj.rm.stacking$fit$bestTune
rf.rm.stacking.permute <- rf.permutation(cleandata[,-c(167:180)], SMOTE = T, permute.time = permute.time, rf.par = rf.par)

glmGrid <- expand.grid(alpha=seq(0,1,0.1), lambda = seq(0,1, length =50))
glm.fit.obj.rm.stacking <- calculateGLMweight(cbind(group=cleandata[,1], as.data.frame(scale(cleandata[, -c(1,167:180)]))), cv.folds,repeat.times = repeat.times, metric = "AUC")
feature.data.rm.stacking <- summarizeFeature(glm.fit.obj.rm.stacking$glm.feature.coef, rf.obj.rm.stacking$fit, rf.rm.stacking.permute)
feature.data.rm.stacking$group <- as.character(feature.data.rm.stacking$group)
feature.data.rm.stacking$group[166:209] <- c(rep("nt", 4), rep("upstreamNt",20 ), rep("dnstreamNt",20 ))

## remove stacking features
set.seed(seed)
rf.obj.rm.stacking <- rfPredict(cleandata[,-c(167:224)], mtry.seq, ntree.seq, cv.folds, repeat.times, smote = T)
rf.rm.stacking.par <- rf.obj.rm.stacking$fit$bestTune
rf.rm.stacking.permute <- rf.permutation(cleandata[,-c(167:224)], SMOTE = T, permute.time = permute.time, rf.par = rf.par)

glmGrid <- expand.grid(alpha=seq(0,1,0.1), lambda = seq(0,1, length =50))
glm.fit.obj.rm.stacking <- calculateGLMweight(cbind(group=cleandata[,1], as.data.frame(scale(cleandata[, -c(1,167:224)]))), cv.folds,repeat.times = repeat.times, metric = "AUC")
feature.data.rm.stacking <- summarizeFeature(glm.fit.obj.rm.stacking$glm.feature.coef, rf.obj.rm.stacking$fit, rf.rm.stacking.permute)


### summary prediction and feature
feature.data <- summarizeFeature(glm.fit.obj$glm.feature.coef, rf.obj$fit, rf.permute)
feature.data.stackwitharomatic <- summarizeFeature(glm.fit.obj.stackwitharomatic$glm.feature.coef, rf.obj.nosmote.stackwitharomatic$fit, rf.stackwitharomatic.permute)

feature.data$group <- as.character(feature.data$group)
feature.data$group[180:223] <- c(rep("nt", 4), rep("upstreamNt",20 ), rep("dnstreamNt",20 ))
feature.data$group <- as.factor(feature.data$group)


feature.data.stackwitharomatic$group <- feature.data$group


### output summary results ###
if(write){
  summarizeSample(cleandata, summary.file = summary.file)
  # aafreqCalculation(cleandata, aa.freq.file = aa.freq.file)
  write.table(feature.data, file = featurerankres.file, row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(feature.data.stackwitharomatic, file = "result/feature_stackwitharomatic_importance_summary.txt", row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(feature.data.rm.stacking, file = "result/feature_rm.stacking_rm.nt_importance_summary.txt", row.names = T, col.names = T, quote = F, sep = "\t")
  
  write.table(matrix(rf.obj$fit$results$ROC, ncol=length(mtry.seq), nrow=length(ntree.seq),dimnames = list(ntree.seq, mtry.seq)), file = rf.tune.roc.mat.file,row.names = T, col.names = T, quote = F, sep = "\t")
}

### plot figures ###
if(plot){
  width <- 3.42*1.5
  height <- 3.42*1.5
  # ## plot model tune
  # plot(rf.obj$fit, ylim=c(0,1))
  # graph2office(file= modeltuneplot.file, type = "ppt")
  # plotRF(rf.obj,tune.plot = T, modeltuneplot.file = modeltuneplot.file, rocplot.file = rocplot.file, width, height, seed)
  set.seed(seed)
  rocplot.file <- paste("plot/rf_rocwithCI_mtry",rf.obj$fit$bestTune$mtry, "_ntree", rf.obj$fit$bestTune$ntree, ".pptx", sep="")
  plotROCwithCI(rf.obj$fit$pred$obs, rf.obj$fit$pred$xl, plot.file = rocplot.file, width, height)
  
  rocplot.file.stackwitharomatic <- paste("plot/rf_nosmote_stackwitharomatic_rocwithCI_mtry",rf.obj.nosmote.stackwitharomatic$fit$bestTune$mtry, "_ntree", rf.obj.nosmote.stackwitharomatic$fit$bestTune$ntree, ".pptx", sep="")
  plotROCwithCI(rf.obj.nosmote.stackwitharomatic$fit$pred$obs, rf.obj.nosmote.stackwitharomatic$fit$pred$xl, plot.file =  rocplot.file.stackwitharomatic, width, height)
  
  rocplot.file.rm.stacking <- paste("plot/rf_rm.stacking_rocwithCI_mtry",rf.obj.rm.stacking$fit$bestTune$mtry, "_ntree", rf.obj.rm.stacking$fit$bestTune$ntree, ".pptx", sep="")
  plotROCwithCI(rf.obj.rm.stacking$fit$pred$obs, rf.obj.rm.stacking$fit$pred$xl, plot.file =  rocplot.file.rm.stacking, width, height)
  
  glm.rocplot.file <- paste("plot/glm_rocwithCI_alpha",glm.fit.obj$glm.fit$bestTune$alpha, "_lambda", glm.fit.obj$glm.fit$bestTune$lambda, ".pptx", sep="")
  plotROCwithCI(glm.fit.obj$glm.fit$pred$obs, glm.fit.obj$glm.fit$pred$xl, plot.file = glm.rocplot.file, width, height)

  plotFeature.volcano(feature.data=feature.data, featurerankplot.file=featurerankplot.file, width=8, height=4)
}