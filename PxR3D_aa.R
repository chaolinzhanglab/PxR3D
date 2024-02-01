source("functions.R")

options(warn=-1)
suppressMessages(library(getopt, quietly=T))
#argument mask: 0 - no argument, 1 - required argument, 2 - optional argument

optionSpec = matrix(c(
  'data.file',                 'i',     1, "character",
  'cv.folds',						'c',		 2, "character",
  'permute.time',           'p',     2, "character",
  'model.file',           'o',    1,  "character",
  'plot.file',             'f',     2, "character",
  'verbose',              'v',     2, "integer",
  'help'   ,              'h',     0, "logical"
), byrow=TRUE, ncol=4)
opt = getopt(optionSpec)

## default parameters
data.file <- opt$data.file
cv.folds <- 10
model.file <- opt$model.file
permute.time <- 2000
verbose <- 0
seed <- 30

if (!is.null(opt$data.file)) {data.file <- opt$data.file}
if (!is.null(opt$cv.folds)) {cv.folds <- opt$cv.folds}
if (!is.null(opt$plot.file)) {plot.file <- opt$plot.file}
if (!is.null(opt$model.file)) {model.file <- opt$model.file}
if (!is.null(opt$permute.time)) {permute.time <- opt$permute.time}
if (!is.null(opt$verbose)) {verbose <- opt$verbose}

if ( !is.null(opt$help) | is.null(opt$data.file) | is.null(opt$model.file))
{
  cat (
    'Predict crosslinking amino acids by PxR3D\n',
    'Usage: Rscript ', get_Rscript_filename(),"[options] -i <feature mat file> -o <model file>\n",
    ' -i, feature mat file     feature matrix file with crosslinking information\n',
    ' -o, model file(Rds format)           final model file\n',
    '[options]\n',
    ' -c, CV folds             Integer(default: 10)\n',
    ' -p, permutation time     Integer(default: 2000)\n',
    ' -f, plot file            ROC plot file\n',
    ' -v, verbose              Verbose mode\n',
    ' -h, help                 Print usage\n'
  );
  q(status=1);
}

### read feature table ###
if(verbose){
  cat("Preproces data\n")
}
cleandata <- read.table(data.file, stringsAsFactors = F, check.names = F, header = T, sep = "\t")

## set mtry and ntree
mtry.seq <- seq(1,10)
ntree.seq <- c(50,1000,2000,5000)
set.seed(seed)
repeat.times <- 1

if(verbose){
  cat("Model training\n")
}
## Predict crosslinking nt using random forest with SMOTE for imbalanced datasets
rf.obj <- rfPredict(cleandata, mtry.seq, ntree.seq, cv.folds, repeat.times, smote = T)
rf.par <-  rf.obj$fit$bestTune

### permute features to get feature rank pvalue ###
set.seed(seed)
rf.permute <- rf.permutation(cleandata, SMOTE = T, permute.time = permute.time, rf.par = rf.par)

### calculate GLM weight ###  
set.seed(seed)
glmGrid <- expand.grid(alpha=seq(0,1,0.1), lambda = seq(0,1, length =20))
glm.fit.obj <- calculateGLMweight(cbind(group=cleandata[,1], as.data.frame(scale(cleandata[,-1]))), cv.folds,repeat.times = repeat.times, metric = "AUC")

if(verbose){
  cat("Saving model\n")
}
if(!is.null(model.file)){
  saveRDS(rf.obj$fit$finalModel, model.file)
}

if(verbose){
  cat("Plotting ROC curve\n")
}
### plot figures ###
if(!is.null(plot.file)){
  set.seed(seed)
  plotROCwithCI(rf.obj$fit$pred$obs, rf.obj$fit$pred$xl, plot.file = rocplot.file)
}