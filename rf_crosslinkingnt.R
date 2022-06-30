source("script/functions.R")
source("script/calculate_twolevelmodel_localImp_new.R")
data.file <- "raw/featuretable_v1.txt"
cv.folds <- 10
repeat.times <- 1
mtry.seq <- seq(1,20)
ntree.seq <- c(10,100,200,500,1000,2000,5000)
plot <- TRUE
write <- TRUE
seed <- 30
permute.time <- 2000
# feature.fill.color1 <- c("#f9c7a8", "#f39a9b" ,  "#ec6e8d" , "#e072a2", "#a93790", "#7b2397", "#5c1886", "#461266", "#300c46")
# feature.fill.color <- c("#FFC75F","#7DA136", "#845EC2", "#f9c7a8", "#f39a9b" ,  "#ec6e8d" , "#e072a2", "#a93790", "#7b2397","#2C73D2", "#0089BA", "#008F7A", "#BA3B64")
# feature.fill.color <- c("#c51b7d", "#de77ae", "#ffffbf" ,"#f1b6da","#fde0ef", "#f7f7f7", "#e6f5d0","#b8e186", "#7fbc41", "#4d9221","#d7191c", "#fdae61",  '#abd9e9',  "#2c7bb6")
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
# rf.obj.scale <- rfPredict(cbind(group=cleandata[,1], as.data.frame(scale(cleandata[, -1]))), mtry.seq, ntree.seq, cv.folds, repeat.times, smote = T)
#rf.obj.nosmote <- rfPredict(cleandata, mtry.seq, ntree.seq, cv.folds, repeat.times, smote = F)
## only classify the samples has base_aa_stack_aromatic larger than 0
#rf.obj.stackwitharomatic <- rfPredict(cleandata[which(cleandata$base_aa_stack_aromatic!=0),], mtry.seq, ntree.seq, cv.folds, repeat.times)
rf.obj.nosmote.stackwitharomatic <- rfPredict(cleandata[which(cleandata$base_aa_stack_aromatic!=0),], mtry.seq, ntree.seq, cv.folds, repeat.times, smote = F)


# set.seed(10)
### try without caret package
# rf.obj.nocaret <- randomforest_cv_tune(cleandata, mtry.seq = mtry.seq, ntree.seq = ntree.seq,  SMOTE = T, seed = seed)
# rf.obj.nocaret.nosmote <- randomforest_cv_tune(cleandata, mtry.seq = mtry.seq, ntree.seq = ntree.seq,  SMOTE = F , seed = seed)
# rf.obj.nocaret.stackwitharomatic <- randomforest_cv_tune(cleandata[which(cleandata$base_aa_stack_aromatic!=0),], mtry.seq = mtry.seq, ntree.seq = ntree.seq,  SMOTE = T, seed = seed)


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
# glm.fit.obj <- calculateGLMweight(cleandata, cv.folds,repeat.times = repeat.times, metric = "AUC",preProcess = c("scale","center"))
#tmp.glm.fit.obj <- calculateGLMweight(cbind(group=cleandata[,1], as.data.frame(scale(tmp[, -c(1:2)]))), cv.folds,repeat.times = repeat.times, metric = "AUC")

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
#feature.data.rm.stacking$group <- as.character(feature.data.rm.stacking$group)
#feature.data.rm.stacking$group[166:209] <- c(rep("nt", 4), rep("upstreamNt",20 ), rep("dnstreamNt",20 ))

### summary prediction and feature
feature.data <- summarizeFeature(glm.fit.obj$glm.feature.coef, rf.obj$fit, rf.permute)
feature.data.stackwitharomatic <- summarizeFeature(glm.fit.obj.stackwitharomatic$glm.feature.coef, rf.obj.nosmote.stackwitharomatic$fit, rf.stackwitharomatic.permute)

feature.data$group <- as.character(feature.data$group)
feature.data$group[180:223] <- c(rep("nt", 4), rep("upstreamNt",20 ), rep("dnstreamNt",20 ))
feature.data$group <- as.factor(feature.data$group)


feature.data.stackwitharomatic$group <- feature.data$group


###################not run below
set.seed(seed)
rf.local <- randomForest(group ~., cleandata, localImp = T, mtry= rf.obj$fit$bestTune$mtry, ntree = rf.obj$fit$bestTune$ntree)
localimp <- t(rf.local$localImportance)

#=====================================================
###one layer with local importance ###
## check the distribution of RRM, other and nonxl 
xl.idx <- which(cleandata$group == "xl")
xl.rrm.idx <- intersect(xl.idx, grep("RRM", rownames(cleandata)))
xl.other.idx <- xl.idx[-match(xl.rrm.idx, xl.idx)]
sample.domain <- rawdata.full$domain[interact.nt.idx]


### boxplot
ttest.pvalue <- list()
wilcox.pvalue <- list()
for (i in 1:dim(localimp)[2]){
  # tmp <- wilcox.test(localimp[xl.rrm.idx, i], localimp[xl.other.idx, i])
  # if(!is.na(tmp$p.value) & tmp$p.value<0.05){
  #   print(colnames(localimp)[i])
  # }

  melt.data <- melt(list(RRM=localimp[xl.rrm.idx, i], Other=localimp[xl.other.idx, i], Nonxl=localimp[-xl.idx, i]), value.name = "localimp")
  melt.data$domain <- as.factor(c(sample.domain[c(xl.rrm.idx, xl.other.idx)], sample.domain[-xl.idx]))
  ttest.pvalue[[colnames(localimp)[i]]] <- pairwise.t.test(melt.data$localimp, melt.data$L1,p.adjust.method = "BH")$p.value
  wilcox.pvalue[[colnames(localimp)[i]]] <- pairwise.wilcox.test(melt.data$localimp, melt.data$L1,p.adjust.method = "BH")$p.value
  
  ggplot(melt.data, aes(x=L1, y=localimp)) +
    geom_boxplot(outlier.shape = NA)+
    # geom_dotplot(binaxis='y', stackdir='center')+
    geom_jitter(aes(shape=domain, color=L1),position=position_jitter(0.2), size=2)+
    # geom_sina()+
    # stat_summary(fun.data="mean_sdl", width=0.5, geom = "pointrange")+
    scale_color_manual(values = c( "#4d9221","#de77ae","#c51b7d"))+
    scale_shape_manual(values = c(17, 18, 16, 3, 4,  15))+
    scale_x_discrete(limits=c("RRM", "Other", "Nonxl"))+
    # geom_boxplot()+
    theme_classic()+
    labs(x="", y = "Local importance")+
    ylim(-0.1,0.1)
    # theme(legend.position="none")
    ggsave(paste("plot/tmp/",colnames(localimp)[i], ".png", sep=""))
  # graph2office(file=paste("plot/",colnames(localimp)[i], "_v1.pptx", sep=""), width = 5, height=3)
  dev.off()
}
localimp.pvalue <- wilcox.test.localimp(localimp, xl.idx)
localimp.boxplot.obj <- boxplot.localimp(localimp[,which(localimp.pvalue < 0.05)], xl.idx)
show(localimp.boxplot.obj)

ttest.anytrue.feature <- which(unlist(lapply(ttest.pvalue, function(x) return(any(x< 0.05)))))
wilcox.anyture.feature <- which(unlist(lapply(wilcox.pvalue, function(x) return(any(x< 0.05)))))
## output the test pvalue for all features
feature.allp <- cbind(localimp.pvalue, t(as.data.frame(lapply(ttest.pvalue, function(x) {return(as.numeric(x))})))[, -3], t(as.data.frame(lapply(wilcox.pvalue, function(x) {return(as.numeric(x))})))[, -3])
colnames(feature.allp) <- c("xlvsNonxl", "othervsNonxl.ttest", "RRMvsNonxl.ttest", "RRMvsOther.ttest", "othervsNonxl.wilcoxtest", "RRMvsNonxl.wilcoxtest", "RRMvsOther.wilcoxtest")
write.table(feature.allp, file = "result/localimp.feaute.test.pvalue.txt", sep="\t", quote = F)

### get the first five features for each prediction
n <- dim(localimp)[2]
first.n.feautre <- matrix(0, nrow = dim(localimp)[1], ncol = n)
rownames(first.n.feautre) <- rownames(localimp)
for (i in 1:dim(localimp)[1]){
  first.n.feautre[i, ] <- colnames(localimp)[order(localimp[i,], decreasing = T)[1:n]]
}

write.table(first.n.feautre[,1:5], file = "result/firstfive.localimp.feaute.txt", sep="\t", quote = F)

#=====================================================
### two layer classfier combining with local importance ###



#**********test************
set.seed(seed)
rf.twolevel.cv.tune.obj  <- rf_twolevel_cv_tune(cleandata, mtry.twolevel.seq = 1:10, ntree.twolevel.seq = c(100,500,1000,2000,5000), cv.folds = 10, aa.code= aa.code, seed = seed, mtry.seq = 1:5, ntree.seq=c(100,500,1000), perc.over = 1000, perc.under = 200, metric = "AUC")
save.image("rf.twolevel.cv.tune.0904.RData")
## test the interaction code ## 
rf.twolevel.cv.tune.obj.interaction  <- rf_twolevel_cv_tune(cleandata, mtry.twolevel.seq = 1:5, ntree.twolevel.seq =c(100,500,1000,2000,5000), cv.folds = 10, aa.code= interaction.code, seed = seed, mtry.seq = c(1,3,5,7), ntree.seq=c(100,500,1000), perc.over = 1000, perc.under = 200, metric = "AUC")

rf.twolevel.cv.tune.obj.best <- rf.twolevel.cv.tune.obj[[which(rf.twolevel.cv.tune.obj$acc$mtry==rf.twolevel.cv.tune.obj$besttune$mtry & rf.twolevel.cv.tune.obj$acc$ntree == rf.twolevel.cv.tune.obj$besttune$ntree)]]

rf.twolevel.cv.tune.obj.interaction.best <- rf.twolevel.cv.tune.obj.interaction[[which(rf.twolevel.cv.tune.obj.interaction$acc$mtry==rf.twolevel.cv.tune.obj.interaction$besttune$mtry & rf.twolevel.cv.tune.obj.interaction$acc$ntree == rf.twolevel.cv.tune.obj.interaction$besttune$ntree)]]
# rf.twolevel.cv.obj <- rf_twolevel_cv(cleandata, mtry = rf.twolevel.cv.tune.obj$besttune$mtry, ntree = rf.twolevel.cv.tune.obj$besttune$ntree, aa.code = aa.code ,cv.folds= 10,  seed = 10, mtry.seq = 1:5, ntree.seq=c(100,500,1000),perc.over = 1000, perc.under = 100)


rf.twolevel.obj <- rf_twolevel(cleandata, aa.code, mtry = rf.twolevel.cv.tune.obj$besttune$mtry, ntree = rf.twolevel.cv.tune.obj$besttune$ntree,cv.folds= cv.folds,  seed = seed, mtry.seq = 1:5, ntree.seq=c(100,500,1000),perc.over = 1000, perc.under = 200)

rf.twolevel.interaction.obj <- rf_twolevel(cleandata, interaction.code, mtry = rf.twolevel.cv.tune.obj.interaction$besttune$mtry, ntree = rf.twolevel.cv.tune.obj.interaction$besttune$ntree,cv.folds= cv.folds,  seed = seed, mtry.seq = c(1,3,5,7), ntree.seq=c(100,500,1000),perc.over = 1000, perc.under = 200)
# idx <- sort(match(rownames(cleandata), rownames(tmp$localimp))[!is.na(match(rownames(cleandata), rownames(tmp$localimp)))])
# localimp <- tmp$localimp[idx, ]
# tmp.xl.idx <- match(rownames(cleandata)[xl.idx], rownames(localimp))
# 
# heatmap.2(rbind(localimp[-xl.idx, ], localimp[xl.idx,]), RowSideColors = c(rep("black", dim(localimp)[1] - length(xl.idx)), rep("red", length(xl.idx))), labRow = NA , trace = "none", col = bluered, scale = "row")
# 
# new.rf.permute <- rf.permutation(tmp$new.feature.mat, SMOTE = T, permute.time = permute.time, rf.par = tmp1$besttune)
# new.glm.fit.obj <- calculateGLMweight(tmp.new.feature.mat, cv.folds)

#**********test************
  
###sth is wrong here 
#=========================================================================
# rf.twolevel.fit <- randomforest_twolevel_cv(cleandata, aa.code = aa.code)
# set.seed(seed)
# rf.local.twolevel.obj <- calculate2LRFLocalImp(cleandata, aa.code, seed, ntree = c(100, 200, 500, 1000, 2000, 5000), mtry.seq = 1:10, ntree.seq =  c(100, 200, 500, 1000, 2000, 5000), cv.folds = cv.folds, repeat.times = repeat.times)
sample.idx <- match(rownames(cleandata), rownames(rf.twolevel.obj$localimp))[!is.na(match(rownames(cleandata), rownames(rf.twolevel.obj$localimp)))]
localimp <- rf.twolevel.obj$localimp[sample.idx, ]
localimp.xl.idx <- match(rownames(cleandata)[which(cleandata$group == "xl")], rownames(localimp))
localimp.boxplot.obj <- boxplot.localimp(localimp = localimp, xl.idx = localimp.xl.idx, aa.color.code = c(rep("black", 12), aa.color.code, rep("black", 6)))
show(localimp.boxplot.obj)
localimp.pvalue <- wilcox.test.localimp(localimp, xl.idx)

sample.idx <- match(rownames(cleandata), rownames(rf.twolevel.interaction.obj$localimp))[!is.na(match(rownames(cleandata), rownames(rf.twolevel.interaction.obj$localimp)))]
localimp <- rf.twolevel.interaction.obj$localimp[sample.idx, ]
localimp.xl.idx <- match(rownames(cleandata)[which(cleandata$group == "xl")], rownames(localimp))
localimp.interaction.boxplot.obj <- boxplot.localimp(localimp =localimp[c(intersect(grep("RRM", rownames(localimp)), localimp.xl.idx), intersect(grep("KH", rownames(localimp)), localimp.xl.idx)), ], xl.idx = 1:length(intersect(grep("RRM", rownames(localimp)),localimp.xl.idx)), aa.color.code = c(rep("black", 12), aa.color.code[1:9]))
show(localimp.interaction.boxplot.obj)
localimp.pvalue <- wilcox.test.localimp(localimp[c(intersect(grep("RRM", rownames(localimp)), localimp.xl.idx), intersect(grep("KH", rownames(localimp)), localimp.xl.idx)), ], xl.idx = 1:length(intersect(grep("RRM", rownames(localimp)),localimp.xl.idx)))
# new.feature.mat <- rf.local.twolevel.obj$new.feature.data
# new.rf.permute <- rf.permutation(new.feature.mat, SMOTE = T, permute.time = permute.time, rf.par = rf.local.twolevel.obj$newdata.rf.obj$rf.tune.obj$fit$bestTune)
# new.glm.fit.obj <- calculateGLMweight(new.feature.mat, cv.folds)
# 
# 
# ##summary predication and feature on new feature set
# new.feature.data <- summarizeFeature(new.glm.fit.obj$glm.feature.coef, rf.local.twolevel.obj$newdata.rf.obj$rf.tune.obj$fit, new.rf.permute)
new.pred.dat <- summarizePrediction(rf.local.twolevel.obj$newdata.rf.obj$rf.fit, rawdata, interact.nt.idx)

#=========================================================================######
################not run done

### output summary results ###
if(write){
  summarizeSample(cleandata, summary.file = summary.file)
  # aafreqCalculation(cleandata, aa.freq.file = aa.freq.file)
  write.table(feature.data, file = featurerankres.file, row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(feature.data.stackwitharomatic, file = "result/feature_stackwitharomatic_importance_summary.txt", row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(feature.data.rm.stacking, file = "result/feature_rm.stacking_rm.nt_importance_summary.txt", row.names = T, col.names = T, quote = F, sep = "\t")
  # write.table(new.feature.data, file = aafeaturerankres.file, row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(cbind(pred.dat, localimp), file = prediction.file, row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(matrix(rf.obj$fit$results$ROC, ncol=length(mtry.seq), nrow=length(ntree.seq),dimnames = list(ntree.seq, mtry.seq)), file = rf.tune.roc.mat.file,row.names = T, col.names = T, quote = F, sep = "\t")
  # write.table(new.pred.dat, file = aaprediction.file, row.names = T, col.names = T, quote = F, sep = "\t")
  # write.table(localimp, file = aafeatureimp.file, row.names = T, col.names = T, quote = F, sep = "\t")
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
  
  #plotFeature(feature.data, feature.fill.color = feature.fill.color, featurerankplot.file, 8, 4)
  plotFeature.volcano(feature.data=feature.data, featurerankplot.file=featurerankplot.file, width=8, height=4)
  #plotFeature(feature.data.stackwitharomatic, feature.fill.color = feature.fill.color, featurerankplot.stackwitharomatic.file, 8, 4)
  
  # plotFeature(new.feature.data[which(sapply(as.character(new.feature.data$group), nchar)==0), ], featurerankplot.file =  aafeaturerankplot.file, width = 9, height = height)
  
 ## oompare two layer and full model
  
  roc(rf.twolevel.cv.tune.obj.best$rf.pred.prob$obs,rf.twolevel.cv.tune.obj.best$rf.pred.prob$xl, legacy.axes = T, asp = 0, plot = T, col="#ca0020")
  plot(roc(rf.obj$fit$pred$obs, rf.obj$fit$pred$xl),add = T, col="black")
  plot(roc(rf.twolevel.cv.tune.obj.interaction.best$rf.pred.prob$obs, rf.twolevel.cv.tune.obj.interaction.best$rf.pred.prob$xl),add = T, col="#0571b0")
  legend(x=0.5, y=0.2, legend = c(paste("Full model : AUC = ", round(roc(rf.obj$fit$pred$obs, rf.obj$fit$pred$xl)$auc, 3), sep=""), paste("Two-layer AA model: AUC = ",round(roc(rf.twolevel.cv.tune.obj.best$rf.pred.prob$obs,rf.twolevel.cv.tune.obj.best$rf.pred.prob$xl)$auc, 3), sep =""), paste("Two-layer interaction model: AUC = ", round(roc(rf.twolevel.cv.tune.obj.interaction.best$rf.pred.prob$obs, rf.twolevel.cv.tune.obj.interaction.best$rf.pred.prob$xl)$auc, 3), sep = "")), col =c("black", "#ca0020", "#0571b0"), lty=1)
  
  graph2ppt(file = modelcomapareplot.file, width=width, height=height)

  heatmap.2(rbind(localimp[-xl.idx,  ], localimp[xl.idx,]), RowSideColors = c(rep("black", dim(cleandata)[1] - length(xl.idx)), rep("red", length(xl.idx))), trace = "none", col = bluered, scale = "row", Rowv = F, dendrogram = F)
  graph2ppt(file = featurelocalimpplot.file)
  local.plot.obj <- boxplot.localimp(localimp = rf.twolevel.obj$localimp, xl.idx = xl.idx,aa.color.code = aa.color.code[1:20])
  

}