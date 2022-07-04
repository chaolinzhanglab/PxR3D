library(openxlsx)
library(dplyr)
require(resha)
source("script/functions.R")

get_closet_nt_rbp <- function(pdbid, pdb.mat){
  file <- paste("raw/protein2rna/pdb", tolower(pdbid), ".interface.aa.txt", sep="")
  mat <- read.table(file, header=T, sep="\t", check.names = F, stringsAsFactors = F, comment.char = "", quote = "")
  mat$id_edited <- gsub("-", "",paste(mat$chain_id, mat$aa_id, sep="."))
  #print(dim(pdb.mat))
  #print(dim(mat))
  idx <- match( pdb.mat$idstr,mat$id_edited)
  if(!identical(idx, 1:dim(pdb.mat)[1])){
    print("Some aa does not have a match.")
  }
  mat.match <- as.data.frame(matrix(0, ncol = 9, nrow=dim(pdb.mat)[1]))
  colnames(mat.match) <- colnames(mat)[5:12]
  mat.match[which(!is.na(idx)),] <- mat[idx, 5:12]
  return(cbind(pdb.mat,mat.match))
}


data.file <- "raw/PDB.Bae.XL.3dna.summary.xlsx"
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
# aa.freq.file <- "result/aa_interaction_frequency.txt"
prediction.file <- "result/bae_RF_prediction_result.txt"
# aaprediction.file <-  "result/RF_AA_prediction_result.txt"
featurerankres.file <- "result/bae_feature_importance_summary.txt"
# aafeaturerankres.file <- "result/aa_feature_importance_summary.txt"
# aafeatureimp.file <- "result/aa_feature_localimp.txt"
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
## no need for standilization

##########remove reduanct strutures and chain#############

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
## pum1 and pum2 has two id with exact same mean_info, only keep one randomly
pdb.keep <- pdb.info[match(pdb.keep.id, pdb.info$pdb), ]
pdb.keep$RBD <- pdb.id.info$RBD[match(pdb.keep$pdb, pdb.id.info$struct_id)]
# keep only one chain for each pdb id
pdb.chain.info <- NULL
for (i in 1:length(pdb.keep$pdb)){
  pdb.chain.info <- rbind(pdb.chain.info, filter.data.chain(raw.data[raw.data$pdb==pdb.keep$pdb[i], ]))
}
chain.keep <- pdb.chain.info %>% group_by(pdb) %>% top_n(non_zero, n=1) 
chain.keep <- as.data.frame(chain.keep)[c(1:11,13:61), ]

#paste(chain.keep$chain, chain.keep$pdb,sep="_")
keep.idx <- which(paste(raw.data$chain, raw.data$pdb,sep="_") %in% paste(chain.keep$chain, chain.keep$pdb,sep="_"))
raw.data.filtered <- raw.data[keep.idx, ]
raw.data.withnt.filtered <- raw.data.withnt[keep.idx, ]

## check feature variance
feature.var <- apply(raw.data.filtered[,11:46], 2, var)
## all larger than 0, no need to filter 


##########check the basic charateristics of the data#############
raw.data.filtered$xl <- 0
raw.data.filtered$xl[raw.data.filtered$xl_score>0] <- 1
raw.data.withnt.filtered$nt <- unlist(lapply(raw.data.withnt.filtered$closest_nuc, function(x) return((strsplit(x, "-")[[1]][1]))))

raw.data.filtered.bak <- raw.data.filtered
raw.data.withnt.filtered.bak <- raw.data.withnt.filtered 


## remove the amino acids not contacting any RNA
nonzerofeature.idx <- which(!apply(raw.data.withnt.filtered[,11:46], 1, function(x){all(x==0)}))
final.keep.idx <- keep.idx[nonzerofeature.idx]

keep.idx.out <- rep(0,dim(raw.data)[1])
keep.idx.out[keep.idx] <- 1
final.keep.idx.out <- rep(0,dim(raw.data)[1])
final.keep.idx.out[final.keep.idx] <- 1

write.table(keep.idx.out, quote=F, sep="\t", file="result/bae.remove.reduandancy.txt")
write.table(final.keep.idx.out, quote=F, sep="\t", file="result/bae.remove.reduandancy.remove.uncontacted.txt")

raw.data.filtered <- raw.data.filtered[nonzerofeature.idx, ]
raw.data.withnt.filtered <- raw.data.withnt.filtered[nonzerofeature.idx, ]

#number of positive and negative samples
table(raw.data.filtered$xl)
## closet nt in the xl and nonxl group
nt.comp.cnt <- cbind(nonxl=table(raw.data.withnt.filtered$nt[raw.data.withnt.filtered$xl_score==0]),xl=
table(raw.data.withnt.filtered$nt[raw.data.withnt.filtered$xl_score>0]))

## aa comp in the xl and nonxl group
aa.names <- names(table(raw.data.withnt.filtered$aacode[raw.data.withnt.filtered$xl_score==0]))
raw.data.withnt.filtered$aacode <- factor(raw.data.withnt.filtered$aacode, levels = aa.names)

aa.comp.cnt <- cbind(nonxl=table(raw.data.withnt.filtered.bak$aacode[raw.data.withnt.filtered.bak$xl_score==0]),xl=
                       table(raw.data.withnt.filtered.bak$aacode[raw.data.withnt.filtered.bak$xl_score>0]))

write.table(nt.comp.cnt, quote=F, sep="\t", file="result/bae.nt.comp.cnt.txt")
write.table(aa.comp.cnt, quote=F, sep="\t", file="result/bae.aa.comp.cnt.txt")
## check close nt associated with each aa

xl.aa.nt.cnt <- raw.data.withnt.filtered[raw.data.withnt.filtered$xl_score>0,-58] %>% 
   group_by(aacode) %>%
   summarise(A=length(which(nt=="A")), C=length(which(nt=="C")), G=length(which(nt=="G")), U=length(which(nt=="U")) )

nonxl.aa.nt.cnt <- raw.data.withnt.filtered[raw.data.withnt.filtered$xl_score==0,-58] %>% 
    group_by(aacode) %>%
    summarise(A=length(which(nt=="A")), C=length(which(nt=="C")), G=length(which(nt=="G")), U=length(which(nt=="U")) )

write.table(xl.aa.nt.cnt , quote=F, sep="\t", file="result/bae.xl.aa.nt.cnt.txt")
write.table(nonxl.aa.nt.cnt , quote=F, sep="\t", file="result/bae.nonxl.aa.nt.cnt.txt")


xl.aa.basestacking.cnt <- raw.data.withnt.filtered[raw.data.withnt.filtered$xl_score>0,c(9,43:46,60)] %>% 
  group_by(aacode) %>%
  summarise(A=length(which(base_aa_stack_A>0)), C=length(which(base_aa_stack_C>0)),G=length(which(base_aa_stack_G>0)),U=length(which(base_aa_stack_U>0)) )

nonxl.aa.basestacking.cnt <- raw.data.withnt.filtered[raw.data.withnt.filtered$xl_score==0, c(9,43:46,60)] %>% 
  group_by(aacode) %>%
  summarise(A=length(which(base_aa_stack_A>0)), C=length(which(base_aa_stack_C>0)),G=length(which(base_aa_stack_G>0)),U=length(which(base_aa_stack_U>0)) )

write.table(xl.aa.basestacking.cnt , quote=F, sep="\t", file="result/bae.xl.aa.basestacking.cnt.txt")
write.table(nonxl.aa.basestacking.cnt , quote=F, sep="\t", file="result/bae.nonxl.aa.basestacking.cnt.txt")
##check the features by group. no significant difference was found
feature.idx <- 11:46
for (i in feature.idx){
  boxplot(list(xl=raw.data.filtered[raw.data.filtered$xl==1,i], nonxl=raw.data.filtered[raw.data.filtered$xl==0,i]),main=colnames(raw.data.filtered)[i], xlab=NULL, outline = F)
}
feature.idx <- 11:46
for (i in feature.idx){
  plot(density(raw.data.filtered[raw.data.filtered$xl==1,i]), col="red",main=colnames(raw.data.filtered)[i], xlab=NULL)
  lines(density(raw.data.filtered[raw.data.filtered$xl==0,i]), add=T)
}

group.feature.mean <- raw.data.filtered %>% 
   group_by(xl) %>%
   summarise_at(vars(interact_A:base_aa_stack_U), function(x) mean(x[x>0]))
as.data.frame(group.feature.mean)
write.table(as.data.frame(group.feature.mean), quote=F, sep="\t", file="result/bae.group.feature.mean.txt")

group.feature.cnt <- raw.data.filtered %>%
  group_by(xl) %>%
  summarise_at(vars(interact_A:base_aa_stack_U), function(x) length(which(x>0)))
# write.table(as.data.frame(group.feature.cnt), quote=F, sep="\t", file="result/bae.group.feature.cnt.txt")

#View(group.feature.cnt)

feature.group.cnt <- cbind(apply(raw.data.filtered[,11:14], 1, function(x) {as.numeric(any(x>0))}),
apply(raw.data.filtered[,15:18], 1, function(x) {as.numeric(any(x>0))}),
apply(raw.data.filtered[,19:22], 1, function(x) {as.numeric(any(x>0))}),
apply(raw.data.filtered[,23:26], 1, function(x) {as.numeric(any(x>0))}),
apply(raw.data.filtered[,27:30], 1, function(x) {as.numeric(any(x>0))}),
apply(raw.data.filtered[,31:34], 1, function(x) {as.numeric(any(x>0))}),
apply(raw.data.filtered[,35:38], 1, function(x) {as.numeric(any(x>0))}),
apply(raw.data.filtered[,39:42], 1, function(x) {as.numeric(any(x>0))}),
apply(raw.data.filtered[,43:46], 1, function(x) {as.numeric(any(x>0))}))

colnames(feature.group.cnt) <- c("interact", "hbond_po4:sidechain",
                                 "hbond_po4:backbone",
                                 "hbond_sugar:sidechain",
                                 "hbond_sugar:backbone",
                                 "hbond_base:sidechain",
                                 "hbond_base:backbone",
                                "base_aa_pair",
                                "base_aa_stack")

feature.group.cnt <- cbind(feature.group.cnt, raw.data.filtered$xl)
colnames(feature.group.cnt)[10] <- "xl"



feature.group.cnt <- as.data.frame(feature.group.cnt)

group.feature.cnt <- feature.group.cnt %>% 
  group_by(xl) %>%
  summarise_at(vars(interact:base_aa_stack), sum)

write.table(as.data.frame(group.feature.cnt), quote=F, sep="\t", file="result/bae.group.feature.cnt.txt")

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
  #summarizeSample(cleandata, summary.file = summary.file)
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

raw.data.out <- cbind(raw.data, final.keep.idx=final.keep.idx.out)
write.table(raw.data.out, file = "result/bae.sample.feature.matrix.txt", row.names = T, col.names = T, quote = F, sep = "\t")
save.image("bae_analysis_v2.Rdata")


write.table(matrix(rf.obj$fit$results$ROC, ncol=length(mtry.seq), nrow=length(ntree.seq),dimnames = list(ntree.seq, mtry.seq)), file = "result/Bae.rf.model.tune.txt",row.names = T, col.names = T, quote = F, sep = "\t")
