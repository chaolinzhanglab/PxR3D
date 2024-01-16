source("functions.R")
data.file <- "example/featuretable.aa.xlsx"
raw.data <- read.xlsx(data.file,sheet=1)
pdb.id.info <- read.xlsx(data.file, sheet=2)
aa.code.name <- read.xlsx(data.file, sheet=3)
cv.folds <- 10
repeat.times <- 1
mtry.seq <- seq(1,10)
ntree.seq <- c(500,1000,2000,5000)
plot <- TRUE
write <- TRUE
seed <- 30
permute.time <- 2000

feature.fill.color <- c("#8dd3c7", "#b3de69","#ffffb3", "#ae017e", "#dd3497", "#f768a1","#fa9fb5", "#fcc5c0","#fde0dd", "#fb8072", "#d7191c","#fdae61",  '#abd9e9',  "#2c7bb6")
aa.color.code <- c("#991122", "#ffddee", "#aa88bb", "#774499", "#3377bb",  "#ddeeff",  "#77aacc", "#225577", "#449944", "#bbddbb", "#224422", "#77bb77", "#999933",  "#dddd33", "#eeee88", "#ffeedd", "#bb6622", "#dd2222", "#ee9999", "#992222")

summary.file <- "result/bae_sample_summary.txt"
prediction.file <- "result/bae_RF_prediction_result.txt"
featurerankres.file <- "result/bae_feature_importance_summary.txt"
 rf.tune.roc.mat.file <- "result/bae_mtry_ntree_roc_mat.txt"

 ## add aa id to get closest nt
raw.data$aa_id <- unlist(lapply(raw.data$idstr, function(x) {unlist(strsplit(x,split = ".", fixed = T))[2]}))

pdb.list <- unique(raw.data$pdb)

raw.data$aa.code.simple <- raw.data$aacode
for (i in 1:dim(aa.code.name)[1]){
  raw.data$aa.code.simple[which(!is.na(match(raw.data$aacode, aa.code.name$Triple[i])))] <- aa.code.name$Single[i]
}
## preprocess data
## add aa identity and up/dn aa as feature
raw.data$up_dipeptide[1] <- "A"
raw.data$dn_dipeptide[1] <- "AS"
raw.data$up_dipeptide[dim(raw.data)[1]] <- "QK"
raw.data$dn_dipeptide[dim(raw.data)[1]] <- "K"
for (i in 2:(dim(raw.data)[1]-1)){
  print(i)
  if(raw.data$pdb[i]==raw.data$pdb[i-1]&raw.data$chain[i]==raw.data$chain[i-1]){
    raw.data$up_dipeptide[i] <- paste(raw.data$aa.code.simple[i-1], raw.data$aa.code.simple[i], sep="")
  }else{
    raw.data$up_dipeptide[i] <- raw.data$aa.code.simple[i]
  }
  if(raw.data$pdb[i]==raw.data$pdb[i+1]&raw.data$chain[i]==raw.data$chain[i+1]){
    raw.data$dn_dipeptide[i] <- paste(raw.data$aa.code.simple[i], raw.data$aa.code.simple[i+1], sep="")
  }else{
    raw.data$dn_dipeptide[i] <- raw.data$aa.code.simple[i]
  }
}



raw.data.withnt <- NULL
for (i in 1:length(pdb.list)){
  print(i)
  raw.data.withnt <- rbind(raw.data.withnt, get_closet_nt_rbp(pdb.list[i], raw.data[raw.data$pdb==pdb.list[i],]))
}

raw.data$RBD <- pdb.id.info$RBD[match(raw.data$pdb, pdb.id.info$struct_id)]
## check the feature value distriubution
feature.idx <- 11:46
for (i in feature.idx){
  hist(raw.data[,i],main=colnames(raw.data)[i], xlab=NULL)
}

## count the xl scores and number of chains for each pdb id
filter.data <- function(raw.data){
  pdb.info <- raw.data %>% 
    group_by(pdb) %>%
    summarise(gene.symbol =unique(gene.symbol),no_xl =sum(xl_score), no_chain=length(unique(chain))) 
  
  ## count how many non-zero value for each pdb id
  pdb.info1 <- raw.data %>% 
    group_by(pdb) %>%
    summarise_at(vars(interact_A:base_aa_stack_U), function(x) length(which(x!=0)))
  pdb.info1 <- as.data.frame(pdb.info1)
  pdb.info1$non_zero <- rowSums(pdb.info1[, -1])
  pdb.info <- cbind(pdb.info, non_zero=pdb.info1$non_zero)
  ## calculare the mean information content for each pdb id by adding xl scores and non zero features divided by number of chains
  pdb.info$mean_info <- (pdb.info$no_xl+pdb.info$non_zero)/pdb.info$no_chain
  pdb.info
}

filter.data.chain <- function(raw.data){
  pdb.info <- raw.data %>% 
    group_by(chain) %>%
    summarise(pdb=unique(pdb),gene.symbol =unique(pdb),no_xl =sum(xl_score), no_chain=length(unique(chain))) 
  
  ## count how many non-zero value for each pdb id
  pdb.info1 <- raw.data %>% 
    group_by(chain) %>%
    summarise_at(vars(interact_A:base_aa_stack_U), function(x) length(which(x!=0)))
  pdb.info1 <- as.data.frame(pdb.info1)
  pdb.info1$non_zero <- rowSums(pdb.info1[, -1])
  pdb.info <- cbind(pdb.info, non_zero=pdb.info1$non_zero)
  ## calculare the mean information content for each pdb id by adding xl scores and non zero features divided by number of chains
  pdb.info$mean_info <- (pdb.info$no_xl+pdb.info$non_zero)/pdb.info$no_chain
  pdb.info
}

## keep the pdbid with maxiunm mean_info
pdb.info <- filter.data(raw.data)
pdb.keep <- pdb.info %>% group_by(gene.symbol) %>% top_n(mean_info, n=1) %>% arrange(gene.symbol)
pdb.keep.id <- pdb.id.info$struct_id[pdb.id.info$keep.idx==1]
info[match(pdb.keep.id, pdb.info$pdb), ]
pdb.keep$RBD <- pdb.id.info$RBD[match(pdb.keep$pdb, pdb.id.info$struct_id)]
# keep only one chain for each pdb id
pdb.chain.info <- NULL
for (i in 1:length(pdb.keep$pdb)){
  pdb.chain.info <- rbind(pdb.chain.info, filter.data.chain(raw.data[raw.data$pdb==pdb.keep$pdb[i], ]))
}
chain.keep <- pdb.chain.info %>% group_by(pdb) %>% top_n(non_zero, n=1) 
chain.keep <- as.data.frame(chain.keep)[c(1:11,13:61), ]

keep.idx <- which(paste(raw.data$chain, raw.data$pdb,sep="_") %in% paste(chain.keep$chain, chain.keep$pdb,sep="_"))
raw.data.filtered <- raw.data[keep.idx, ]
raw.data.withnt.filtered <- raw.data.withnt[keep.idx, ]


## remove the amino acids not contacting any RNA
nonzerofeature.idx <- which(!apply(raw.data.withnt.filtered[,11:46], 1, function(x){all(x==0)}))
final.keep.idx <- keep.idx[nonzerofeature.idx]

keep.idx.out <- rep(0,dim(raw.data)[1])
keep.idx.out[keep.idx] <- 1
final.keep.idx.out <- rep(0,dim(raw.data)[1])
final.keep.idx.out[final.keep.idx] <- 1

raw.data.filtered <- raw.data.filtered[nonzerofeature.idx, ]
raw.data.withnt.filtered <- raw.data.withnt.filtered[nonzerofeature.idx, ]

## reorginze the data for input in RF prediction
cleandata <- cbind(group=raw.data.filtered$xl, raw.data.filtered[,c(11:46, 9, 49:50)])
rownames(cleandata) <- paste(raw.data.filtered$gene.symbol, raw.data.filtered$pdb, raw.data.filtered$idstr, sep="_")
cleandata$group[raw.data.filtered$xl==0] <- "nonxl"
cleandata$group[raw.data.filtered$xl==1] <- "xl"
cleandata$group <- factor(cleandata$group, levels = c("xl", "nonxl"))
cleandata$aacode <- factor(cleandata$aacode)
cleandata$up_dipeptide <- factor(cleandata$up_dipeptide)
cleandata$dn_dipeptide <- factor(cleandata$dn_dipeptide)

numericl.feature.mat <- cbind(model.matrix(~ aacode -1 , cleandata),
model.matrix(~ up_dipeptide -1 , cleandata),
model.matrix(~ dn_dipeptide -1 , cleandata))

cleandata.bak <- cleandata
cleandata <- cbind(cleandata[,1:37],as.data.frame(numericl.feature.mat) )
### predict with randome forest ###
set.seed(seed)
rf.obj <- rfPredict(cleandata, mtry.seq, ntree.seq, cv.folds, repeat.times, smote = T)
### permute features to get feature rank pvalue ###
rf.par <-  rf.obj$fit$bestTune
set.seed(seed)
rf.permute <- rf.permutation(cleandata, SMOTE = T, permute.time = permute.time, rf.par = rf.par)
### calculate GLM weight ###
set.seed(seed)
glmGrid <- expand.grid(alpha=seq(0,1,0.1), lambda = seq(0,10, length =20))
glm.fit.obj <- calculateGLMweight(cbind(group=cleandata[,1], as.data.frame(scale(cbind(cleandata[, 2:37], numericl.feature.mat)))), cv.folds,repeat.times = repeat.times, metric = "AUC",tuneGrid = glmGrid)
### summary prediction and feature
feature.data <- summarizeFeature(glm.fit.obj$glm.feature.coef, rf.obj$fit, rf.permute)
feature.data$p.adjust <- p.adjust(feature.data$pvalue, method="fdr")

pred.dat <- cbind(raw.data.withnt.filtered, rf.obj$fit$pred[order(rf.obj$fit$pred$rowIndex),])


pdb.info$keep.idx <- 0
pdb.info$keep.idx[match(pdb.keep$pdb, pdb.info$pdb)] <- 1

if(write){
  write.table(pdb.info, file = "result/bae.pdb.info.summary.txt", row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(feature.data, file = featurerankres.file, row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(pred.dat, file = prediction.file, row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(matrix(rf.obj$fit$results$ROC, ncol=length(mtry.seq), nrow=length(ntree.seq),dimnames = list(ntree.seq, mtry.seq)), file = rf.tune.roc.mat.file,row.names = T, col.names = T, quote = F, sep = "\t")
}


if(plot){
  width <- 3.42*1.5
  height <- 3.42*1.5
  set.seed(seed)
  rocplot.file <- paste("plot/bae_rf_rocwithCI_mtry",rf.obj$fit$bestTune$mtry, "_ntree", rf.obj$fit$bestTune$ntree, ".pptx", sep="")
  plotROCwithCI(rf.obj$fit$pred$obs, rf.obj$fit$pred$xl, plot.file = rocplot.file, width, height)
  
  glm.rocplot.file <- paste("plot/bae_glm_rocwithCI_alpha",glm.fit.obj$glm.fit$bestTune$alpha, "_lambda", glm.fit.obj$glm.fit$bestTune$lambda, ".pptx", sep="")
  plotROCwithCI(glm.fit.obj$glm.fit$pred$obs, glm.fit.obj$glm.fit$pred$xl, plot.file = glm.rocplot.file, width, height)
  
  plotFeature(feature.data, feature.fill.color = feature.fill.color[c(4:9, 1:3)],featurerankplot.file="plot/bae_feature_importance_scatterplot.pptx", width=8, height=4)
  
}
