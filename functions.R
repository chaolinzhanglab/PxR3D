# library(openxlsx)
library(caret)
# library(plotROC)
library(randomForest)
library(pROC)
# library(corrplot)
library(doParallel)
library(devtools)
library(magrittr)
library(glmnet)
library(DMwR)
# devtools::install_github("tomwenseleers/export")
library(export)
library(gplots)
library(MLmetrics)
library(reshape2)
library(ggplot2)
library(ggforce)
library(parallel) 
#============================================================================================
#Cross validataion using Gradient boosting machine
#============================================================================================
gbm_cv <- function(gbm.dat, n.trees=1000, cv.folds=5, interaction.depth=1, 
                   shrinkage= 0.001, bag.fraction=1, n.cores=1, verbose=TRUE, weights=NULL, SMOTE = T)
{
  # split samples
  sample.idx.group <- split(1:nrow(gbm.dat), gbm.dat$group)
  cv.idx.group <- list()
  for (i in 1:length(sample.idx.group)){
    cv.idx.group[[i]] <- sample(1:cv.folds, length(sample.idx.group[[i]]), replace=TRUE)
  }
  
  cv.idx <- as.integer(rep(0, nrow(gbm.dat)))
  cv.idx[unlist(sample.idx.group)] <- unlist(cv.idx.group)
  #gbm.dat.cv <- cbind(data.frame(cv.idx=cv.idx), gbm.dat)
  
  # cross validation
  rl.gbm.cv <- list()
  gbm.pred <- character(nrow(gbm.dat))
  for (i in 1:cv.folds){
    gbm.dat.train <- gbm.dat[cv.idx != i, ]
    # gbm.dat.test <- gbm.dat[cv.idx == i, ]
    if(SMOTE){
      gbm.dat.train.SMOTE <- SMOTE(group ~ ., gbm.dat.train, perc.over = 1000, perc.under =100)
      gbm.dat.train <- gbm.dat.train.SMOTE
    }
    
    tmp.dat <- filterbyvar(gbm.dat.train)
    gbm.dat.train <- tmp.dat$dat
    var.idx <- tmp.dat$idx
    gbm.dat.test <- gbm.dat[cv.idx == i, var.idx]
    if(is.null(weights)){
      cur.train <- gbm.fit(as.matrix(gbm.dat.train[,-1]), gbm.dat.train[,1], 
                           distribution = 'multinomial', n.trees = n.trees,
                           interaction.depth = interaction.depth, shrinkage = shrinkage,
                           bag.fraction = bag.fraction, verbose = verbose, n.minobsinnode=3)
    }else{
      cur.train <- gbm.fit(as.matrix(gbm.dat.train[,-1]), gbm.dat.train[,1], 
                           distribution = 'multinomial', n.trees = n.trees,
                           interaction.depth = interaction.depth, shrinkage = shrinkage,
                           bag.fraction = bag.fraction, verbose = verbose, n.minobsinnode=3, w=weights)
    }
    
    cur.test <- predict(cur.train, gbm.dat.test, n.trees=n.trees, type="response")
    cur.pred <- colnames(cur.test)[apply(cur.test, 1, which.max)]
    gbm.pred[cv.idx == i] <- cur.pred
    rl.gbm.cv[[i]] <- cur.train
  }
  acc <- sum(as.character(gbm.dat$group) == gbm.pred) / nrow(gbm.dat)
  acc.group <- split(as.character(gbm.dat$group) == gbm.pred, as.character(gbm.dat$group))
  acc.group <- sapply(acc.group, function(x) sum(x)/length(x))
  
  list(rl.gbm.cv= rl.gbm.cv, gbm.pred = gbm.pred, cv.idx = cv.idx, acc=acc, acc.group=acc.group)	
}

#============================================================================================
#Cross validataion using supporting vector machine
#============================================================================================
svm_cv <- function(gbm.dat, cv.folds=5, cost, gamma, n.cores=1, verbose=TRUE, SMOTE=FALSE, ...)
{
  # split samples
  sample.idx.group <- split(1:nrow(gbm.dat), gbm.dat$group)
  cv.idx.group <- list()
  for (i in 1:length(sample.idx.group)){
    cv.idx.group[[i]] <- sample(1:cv.folds, length(sample.idx.group[[i]]), replace=TRUE)
  }
  
  cv.idx <- as.integer(rep(0, nrow(gbm.dat)))
  cv.idx[unlist(sample.idx.group)] <- unlist(cv.idx.group)
  #gbm.dat.cv <- cbind(data.frame(cv.idx=cv.idx), gbm.dat)
  
  # cross validation
  rl.gbm.cv <- list()
  gbm.pred <- character(nrow(gbm.dat))
  for (i in 1:cv.folds){
    cat('fold:',i,'\r')
    gbm.dat.train <- gbm.dat[cv.idx != i, ]
    # gbm.dat.test <- gbm.dat[cv.idx == i, ]
    ##select features have var larger than var.cutoff
    if(SMOTE){
      gbm.dat.train.SMOTE <- SMOTE(group ~ ., gbm.dat.train, perc.over = 1000, perc.under =100)
      gbm.dat.train <- gbm.dat.train.SMOTE
    }
    tmp.dat <- filterbyvar(gbm.dat.train)
    gbm.dat.train <- tmp.dat$dat
    var.idx <- tmp.dat$idx
    gbm.dat.test <- gbm.dat[cv.idx == i, var.idx]
    
    cur.train <- svm(as.matrix(gbm.dat.train[,-1]), gbm.dat.train[,1], gamma = gamma, cost = cost, ...)
    cur.pred <- predict(cur.train, as.matrix(gbm.dat.test[,-1]))
    gbm.pred[cv.idx == i] <- as.character(cur.pred)
    rl.gbm.cv[[i]] <- cur.train
  }	
  cat('fold:',i,'\n')
  acc <- sum(as.character(gbm.dat$group) == gbm.pred) / nrow(gbm.dat)
  acc.group <- split(as.character(gbm.dat$group) == gbm.pred, as.character(gbm.dat$group))
  acc.group <- sapply(acc.group, function(x) sum(x)/length(x))
  
  list(rl.svm.cv= rl.gbm.cv, svm.pred = gbm.pred, cv.idx = cv.idx, acc=acc, acc.group=acc.group)	
}


#============================================================================================
#Cross validataion using generalized linear regression model
#============================================================================================
glm_cv <- function(gbm.dat, cv.folds=5, verbose=TRUE, SMOTE=F, ...)
{
  # split samples
  sample.idx.group <- split(1:nrow(gbm.dat), gbm.dat$group)
  cv.idx.group <- list()
  for (i in 1:length(sample.idx.group)){
    cv.idx.group[[i]] <- sample(1:cv.folds, length(sample.idx.group[[i]]), replace=TRUE)
  }
  
  cv.idx <- as.integer(rep(0, nrow(gbm.dat)))
  cv.idx[unlist(sample.idx.group)] <- unlist(cv.idx.group)
  #gbm.dat.cv <- cbind(data.frame(cv.idx=cv.idx), gbm.dat)
  
  # cross validation
  rl.gbm.cv <- list()
  gbm.pred <- character(nrow(gbm.dat))
  for (i in 1:cv.folds){
    cat('fold:',i,'\r')
    gbm.dat.train <- gbm.dat[cv.idx != i, ]
    # gbm.dat.test <- gbm.dat[cv.idx == i, ]
    ## select features that has var larger than 0.05
    
    if(SMOTE){
      gbm.dat.train.SMOTE <- SMOTE(group ~ ., gbm.dat.train, perc.over = 1000, perc.under =100)
      gbm.dat.train <- gbm.dat.train.SMOTE
    }
    tmp.dat <- filterbyvar(gbm.dat.train)
    gbm.dat.train <- tmp.dat$dat
    var.idx <- tmp.dat$idx
    gbm.dat.test <- gbm.dat[cv.idx == i, var.idx]
    cur.train <- glmnet(as.matrix(gbm.dat.train[,-1]), as.numeric(as.character(gbm.dat.train[,1])), ...)
    cur.pred <- predict(cur.train, as.matrix(gbm.dat.test[,-1]), type = "class", s=0.09)
    gbm.pred[cv.idx == i] <- as.character(cur.pred)
    rl.gbm.cv[[i]] <- cur.train
  }	
  cat('fold:',i,'\n')
  acc <- sum(as.character(gbm.dat$group) == gbm.pred) / nrow(gbm.dat)
  acc.group <- split(as.character(gbm.dat$group) == gbm.pred, as.character(gbm.dat$group))
  acc.group <- sapply(acc.group, function(x) sum(x)/length(x))
  
  list(rl.glm.cv= rl.gbm.cv, glm.pred = gbm.pred, cv.idx = cv.idx, acc=acc, acc.group=acc.group)	
}

#============================================================================================
#Cross validataion using decision tree
#============================================================================================
rtree_cv <- function(gbm.dat, cv.folds=5, verbose=TRUE, SMOTE=F, ...)
{
  # split samples
  sample.idx.group <- split(1:nrow(gbm.dat), gbm.dat$group)
  cv.idx.group <- list()
  for (i in 1:length(sample.idx.group)){
    cv.idx.group[[i]] <- sample(1:cv.folds, length(sample.idx.group[[i]]), replace=TRUE)
  }
  
  cv.idx <- as.integer(rep(0, nrow(gbm.dat)))
  cv.idx[unlist(sample.idx.group)] <- unlist(cv.idx.group)
  #gbm.dat.cv <- cbind(data.frame(cv.idx=cv.idx), gbm.dat)
  
  # cross validation
  rl.gbm.cv <- list()
  gbm.pred <- character(nrow(gbm.dat))
  for (i in 1:cv.folds){
    cat('fold:',i,'\r')
    gbm.dat.train <- gbm.dat[cv.idx != i, ]
    # gbm.dat.test <- gbm.dat[cv.idx == i, ]
    ## select features that has var larger than 0.05
    
    if(SMOTE){
      gbm.dat.train.SMOTE <- SMOTE(group ~ ., gbm.dat.train, perc.over = 1000, perc.under =100)
      gbm.dat.train <- gbm.dat.train.SMOTE
    }
    tmp.dat <- filterbyvar(gbm.dat.train)
    gbm.dat.train <- tmp.dat$dat
    var.idx <- tmp.dat$idx
    gbm.dat.test <- gbm.dat[cv.idx == i, var.idx]
    cur.train <- rpart(group~.,as.data.frame(gbm.dat.train), ...)
    cur.pred <- predict(cur.train, gbm.dat.test[,-1], type="class")
    gbm.pred[cv.idx == i] <- as.character(cur.pred)
    rl.gbm.cv[[i]] <- cur.train
  }	
  cat('fold:',i,'\n')
  acc <- sum(as.character(gbm.dat$group) == gbm.pred) / nrow(gbm.dat)
  acc.group <- split(as.character(gbm.dat$group) == gbm.pred, as.character(gbm.dat$group))
  acc.group <- sapply(acc.group, function(x) sum(x)/length(x))
  
  list(rl.rtree.cv= rl.gbm.cv, rtree.pred = gbm.pred, cv.idx = cv.idx, acc=acc, acc.group=acc.group)	
}

#============================================================================================
# filter features by variance
#============================================================================================
filterbyvar <- function(dat, var.cutoff=0.05){
  ##the first column is group
  res.list <- list()
  feature.var <- apply(dat[,-1], 2, var)
  var.idx <- which(feature.var>var.cutoff)+1
  res.list$dat <- dat[,c(1, var.idx)]
  res.list$idx <- c(1, var.idx)
  return(res.list)
}

#============================================================================================
#using caret for randome forest prediction with no cross validation
#============================================================================================
caret.wrap.nocv <- function(dat,method,SMOTE=TRUE, ...){
  if(!SMOTE){
    fitControl <- trainControl(method = "none", savePredictions = "final", summaryFunction = twoClassSummary, classProbs = T )
  }else{
    fitControl <- trainControl(method = "none", savePredictions = "final", sampling = "smote", summaryFunction = twoClassSummary, classProbs = T )
  }
  
  fit.obj <- train(group ~ ., data = dat, method=method, trControl = fitControl, verbose = FALSE, ...)
  bestobj.res <- fit.obj$pred
  acc <- sum(as.character(bestobj.res$obs) == bestobj.res$pred) / nrow(bestobj.res)
  acc.group <- split(as.character(bestobj.res$obs) == bestobj.res$pred, as.character(bestobj.res$obs))
  acc.group <- sapply(acc.group, function(x) sum(x)/length(x))
  list(fit=fit.obj, acc=acc, acc.group=acc.group)
}

#============================================================================================
#Leave-one-out Cross validataion using randome forest from Caret
#============================================================================================
caret.wrap.leaveoneout <- function(dat, idx, cv.folds, repeat.times, method, ...){
    fitControl <- trainControl(method = "repeatedcv", number = cv.folds, repeats = repeat.times, savePredictions = "final", summaryFunction = twoClassSummary, classProbs = T )
  
  y <- dat[, idx]
  y [y==1] <- "xl"
  y [y==0] <- "nonxl"
  y <- as.factor(y)
  # levels(y) <- c("xl", "nonxl")
  
  fit.obj <- train(dat[,-idx], y, method=method, trControl = fitControl, verbose = FALSE, ...)
  bestobj.res <- fit.obj$pred
  acc <- sum(as.character(bestobj.res$obs) == bestobj.res$pred) / nrow(bestobj.res)
  acc.group <- split(as.character(bestobj.res$obs) == bestobj.res$pred, as.character(bestobj.res$obs))
  acc.group <- sapply(acc.group, function(x) sum(x)/length(x))
  list(fit=fit.obj, acc=acc, acc.group=acc.group)
}

#============================================================================================
#Cross validataion using randome forest from Caret
#============================================================================================

caret.wrap <- function(dat,method, cv.folds, repeat.times,SMOTE=TRUE, ...){
  # #library(parallel) 
  # # Calculate the number of cores
  # no_cores <- detectCores() - 1
  # 
  # library(doParallel)
  # # create the cluster for caret to use
  # cl <- makePSOCKcluster(no_cores)
  # registerDoParallel(cl)
  
  # do your regular caret train calculation enabling
  # allowParallel = TRUE for the functions that do
  # use it as part of their implementation. This is
  # determined by the caret package.
  
  if(!SMOTE){
    fitControl <- trainControl(method = "repeatedcv", number = cv.folds,repeats = repeat.times, savePredictions = "final", summaryFunction = twoClassSummary, classProbs = T )
  }else{
    fitControl <- trainControl(method = "repeatedcv", number = cv.folds,repeats = repeat.times, savePredictions = "final", sampling = "smote", summaryFunction = twoClassSummary, classProbs = T )
  }
  
  fit.obj <- train(group ~ ., data = dat, method=method, trControl = fitControl, verbose = FALSE,allowParallel=T, ...)
  
  # stopCluster(cl)
  # registerDoSEQ()
  bestobj.res <- fit.obj$pred
  acc <- sum(as.character(bestobj.res$obs) == bestobj.res$pred) / nrow(bestobj.res)
  acc.group <- split(as.character(bestobj.res$obs) == bestobj.res$pred, as.character(bestobj.res$obs))
  acc.group <- sapply(acc.group, function(x) sum(x)/length(x))
  list(fit=fit.obj, acc=acc, acc.group=acc.group)
}


caret.wrap.oob <- function(dat,method, SMOTE=TRUE, ...){
  
  if(!SMOTE){
    fitControl <- trainControl(method = "oob", savePredictions = "final", summaryFunction = twoClassSummary, classProbs = T,allowParallel=T )
  }else{
    fitControl <- trainControl(method = "oob", savePredictions = "final", sampling = "smote", summaryFunction = twoClassSummary, classProbs = T, allowParallel=T)
  }
  
  fit.obj <- train(group ~ ., data = dat, method=method, trControl = fitControl, verbose = FALSE,allowParallel=T, ...)
  
  # stopCluster(cl)
  # registerDoSEQ()
  bestobj.res <- fit.obj$pred
  acc <- sum(as.character(bestobj.res$obs) == bestobj.res$pred) / nrow(bestobj.res)
  acc.group <- split(as.character(bestobj.res$obs) == bestobj.res$pred, as.character(bestobj.res$obs))
  acc.group <- sapply(acc.group, function(x) sum(x)/length(x))
  list(fit=fit.obj, acc=acc, acc.group=acc.group)
}

#============================================================================================
#Customize randome forest function allowing tuning mtry and ntree
#============================================================================================
customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

#============================================================================================
#permute sample to get feautre importance p value 
#============================================================================================
rf.permutation <- function(dat, SMOTE, permute.time, rf.par, ... ){
  # set.seed(10)
  feature.imp.mat <- matrix(0, nrow = dim(dat)[2]-1, ncol = permute.time)
  rownames(feature.imp.mat) <- colnames(dat)[-1]
  feature.imp.scale.mat <- feature.imp.mat
  pb <- txtProgressBar(min = 0, max = permute.time-1, style = 3)
  for (i in 1:permute.time){
    # progress(i-1, progress.bar = TRUE)
    #Sys.sleep(0.01)
    permute.fit <- rf.permutation.sinlge(dat, rf.par, permute = T, SMOTE = SMOTE )
    feature.imp.mat[,i] <- permute.fit$imp
    feature.imp.scale.mat[,i] <- permute.fit$imp.scale
    setTxtProgressBar(pb, i-1)
    if (i == permute.time) cat("Done!\n")
  }
  close(pb)
  feature.imp.true <- rf.permutation.sinlge(dat, rf.par, permute = F, SMOTE = SMOTE)
  feature.imp.mat <- cbind(feature.imp.true$imp, feature.imp.mat)
  feature.imp.scale.mat <- cbind(feature.imp.true$imp.scale, feature.imp.scale.mat)
  
  feature.imp.rank <- apply(feature.imp.mat, 2, function(x) return(rank(-x)))
  rownames(feature.imp.rank) <- colnames(dat)[-1]
  
  feature.pvalue.rank <- apply(feature.imp.rank, 1, function(x) { length(which(x[-1]<=x[1]))/permute.time})
  feature.pvalue.imp <- apply(feature.imp.mat, 1, function(x) { length(which(x[-1]>=x[1]))/permute.time})
  feature.pvalue.imp.scale <- apply(feature.imp.scale.mat, 1, function(x) { length(which(x[-1]>=x[1]))/permute.time}) 
  
  list(importance.matrix = feature.imp.mat, rank.pvalue = feature.pvalue.rank, importance.pvalue = feature.pvalue.imp, scaled.importance.pvalue = feature.pvalue.imp.scale)
}

#============================================================================================
rf.permutation.sinlge <- function(dat, permute, SMOTE = TRUE, rf.par, ...){
  if(permute){
    idx.rd <- sample(1:nrow(dat), table(dat$group)[1])
    dat$group[idx.rd] <- levels(dat$group)[1]
    dat$group[-idx.rd] <- levels(dat$group)[2]
  }
  fitCtrl <- trainControl(method = "none", classProbs = TRUE)
  if(SMOTE){
    fitCtrl <- trainControl(method = "none", classProbs = TRUE, sampling = "smote")
  }
  rd.fit <- train(group ~., data= dat, method = customRF, trControl = fitCtrl, tuneGrid = rf.par, ...)
  feature.imp <- varImp(rd.fit, scale = F)$importance
  feature.imp.scale <- varImp(rd.fit, scale = T)$importance
  list(imp = feature.imp[,1], imp.scale = feature.imp.scale [, 1])
}

#============================================================================================
#customized glm function allowing tuning lamda and alpha
#============================================================================================
my_glmnet <- getModelInfo("glmnet") %>% magrittr::extract2("glmnet")
my_glmnet$grid <- function (x, y, len = NULL, search = "grid") {
  if (search == "grid") {
    numLev <- if (is.character(y) | is.factor(y)) 
      length(levels(y))
    else NA
    if (!is.na(numLev)) {
      fam <- ifelse(numLev > 2, "multinomial", "binomial")
    }
    else fam <- "gaussian"
    init <- glmnet(as.matrix(x), y, family = fam, nlambda = 52, alpha = 0.5)
    lambda <- unique(init$lambda)
    lambda <- lambda[-c(1, length(lambda))]
    l_seq <- seq(1, length(lambda), length = len) %>% round %>% unique
    lambda <- lambda[l_seq]
    out <- expand.grid(alpha = seq(0.1, 1, length = len), 
                       lambda = lambda)
  }
  else {
    out <- data.frame(alpha = runif(len, min = 0, 1), lambda = 2^runif(len, 
                                                                       min = -10, 3))
  }
  out
}
my_glmnet$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
my_glmnet$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "response")
my_glmnet$sort <- function(x) x[order(x[,1]),]
my_glmnet$levels <- function(x) x$classes

#============================================================================================
randomforest_cv_pred.imp <- function(gbm.dat, cv.folds=5, verbose=TRUE, aa.code, ...)
{
  # set.seed(30)
  # split samples
  sample.idx.group <- split(1:nrow(gbm.dat), gbm.dat$group)
  cv.idx.group <- list()
  for (i in 1:length(sample.idx.group)){
    cv.idx.group[[i]] <- sample(1:cv.folds, length(sample.idx.group[[i]]), replace=TRUE)
  }
  
  cv.idx <- as.integer(rep(0, nrow(gbm.dat)))
  cv.idx[unlist(sample.idx.group)] <- unlist(cv.idx.group)
  #gbm.dat.cv <- cbind(data.frame(cv.idx=cv.idx), gbm.dat)
  
  # cross validation
  rfGrid <- expand.grid(.mtry= seq(1,5), .ntree=c(10,100,1000))
  res <- list()
  for (x in 0:length(aa.code)){
    if(x ==0 ){
      # res[["fullmodel"]]$fit <- NULL
      res[["fullmodel"]] <- data.frame(pred=character(nrow(gbm.dat)),obs = character(nrow(gbm.dat)),xl=numeric(nrow(gbm.dat)), nonxl=numeric(nrow(gbm.dat)), stringsAsFactors = F)
    }else{
      # res[[aa.code[x]]]$fit <- NULL
      res[[aa.code[x]]] <- data.frame(pred=character(nrow(gbm.dat)),obs = character(nrow(gbm.dat)),xl=numeric(nrow(gbm.dat)), nonxl=numeric(nrow(gbm.dat)), stringsAsFactors = F)
    }
    
  }
  
  for (i in 1:cv.folds){
    cat('fold:',i,'\r')
    gbm.dat.train <- gbm.dat[cv.idx != i, ]
    gbm.dat.test <- gbm.dat[cv.idx == i, ]
   
    for (j in 0:length(aa.code)){
      if(j == 0){
        cur.train <- caret.wrap(gbm.dat.train, method=customRF, cv.folds = 5, repeat.times = 1, SMOTE = T, tuneGrid=rfGrid, metric = "ROC")
        cur.pred <- predict( cur.train$fit$finalModel, gbm.dat.test[, -1])
        pred.prob <- predict(cur.train$fit$finalModel, gbm.dat.test[, -1], type = "prob", n.trees = cur.pred.sub$fit$bestTune$n.trees)
      res[["fullmodel"]]$xl[cv.idx == i ] <- pred.prob[, which(colnames(pred.prob)=="xl")]
      res[["fullmodel"]]$nonxl[cv.idx == i ] <- pred.prob[, which(colnames(pred.prob)=="nonxl")]
      
      res[["fullmodel"]]$pred[cv.idx == i] <- as.character(cur.pred)
      res[["fullmodel"]]$obs[cv.idx == i] <- gbm.dat.test$group
      
      }else{
        rm.idx <- grep(aa.code[j], colnames(gbm.dat.train))
        cur.train <- caret.wrap(gbm.dat.train[, -rm.idx], method=customRF, cv.folds = 5, repeat.times = 1, SMOTE = T, tuneGrid=rfGrid, metric = "ROC")
        cur.pred <- predict( cur.train$fit$finalModel, gbm.dat.test[, -c(1, rm.idx)])
        pred.prob <- predict(cur.train$fit$finalModel, gbm.dat.test[, -c(1, rm.idx)], type = "prob", n.trees = cur.pred.sub$fit$bestTune$n.trees)
        res[[aa.code[j]]]$xl[cv.idx == i ] <- pred.prob[, which(colnames(pred.prob)=="xl")]
        res[[aa.code[j]]]$nonxl[cv.idx == i ] <- pred.prob[, which(colnames(pred.prob)=="nonxl")]
        
        res[[aa.code[j]]]$pred[cv.idx == i] <- as.character(cur.pred)
        res[[aa.code[j]]]$obs[cv.idx == i] <- gbm.dat.test$group
        
      }
    }
    
  }	
  
  cat('fold:',i,'\n')
 res
}
### rowIndex is wrong

#============================================================================================
#Plot roc curve with confidence interval using roc function
#============================================================================================
plotROCwithCI <- function(obs, pred, plot.file, width, height, ...){
  rf.fit.rocobj <- roc(obs, pred, ci=T, legacy.axes = T, asp = 0,identity.lty=6, plot = T)
  ciobj <- ci.se(rf.fit.rocobj, specificities=seq(0, 1, 0.02)) 
  plot(ciobj, type="shape", border = NA)
  ### needs to re-add roc curve because of the overlay of color
  plot(rf.fit.rocobj, col="steelblue", add = T)
  text(x=c(0.15,0.15), y=c(0.1, 0.05), labels = c(paste("AUC = ", round(rf.fit.rocobj$auc, 3)), paste( "95% CI : ", round(rf.fit.rocobj$ci[1],3), "-", round(rf.fit.rocobj$ci[3],3), sep="")))
  graph2ppt(file= plot.file, width = width, height = height, ...)
}

#============================================================================================
#Wrapper for random forest prediction
#============================================================================================
rfPredict <- function(pdbrna.data.clean, mtry.seq, ntree.seq, cv.folds, repeat.tim, smote = T,...){
  
  rfGrid <- expand.grid(.mtry=mtry.seq, .ntree=ntree.seq)
  rf.tune.obj <- caret.wrap(pdbrna.data.clean, method=customRF, cv.folds = cv.folds , repeat.times = repeat.times , SMOTE = smote, tuneGrid=rfGrid, metric = "ROC",...)
  
  ### decide final model ###
  # rf.par <- rf.tune.obj$fit$bestTune
  # fitCtrl <- trainControl(method = "cv", summaryFunction=twoClassSummary, savePredictions = T, preProcOptions = c("center", "scale"),  classProbs = T, sampling = "smote",verboseIter = TRUE)
  # ## weights or not??
  # # model_weights <- ifelse(pdbrna.data.clean$group == levels(pdbrna.data.clean$group)[1],
  # # (1/table(pdbrna.data.clean$group)[1]) * 0.5,
  # # (1/table(pdbrna.data.clean$group)[2]) * 0.5)
  # rf.fit <- train(group ~., data= pdbrna.data.clean, method = customRF, trControl = fitCtrl, tuneGrid = rf.par)
  # list(rf.tune.obj = rf.tune.obj, rf.fit = rf.fit) 
  rf.tune.obj
}


rfPredict.oob <- function(pdbrna.data.clean, mtry.seq, ntree.seq, smote = T,...){
  
  #rfGrid <- expand.grid(.mtry=mtry.seq, .ntree=ntree.seq)
  rfGrid <- expand.grid(.mtry=mtry.seq)
  
  rf.tune.obj <- caret.wrap.oob(pdbrna.data.clean, method="rf", SMOTE = smote, tuneGrid=rfGrid, metric = "ROC",...)
  rf.tune.obj
}

#============================================================================================
#Plot tunning plot and model roc for random forest
#============================================================================================
plotRF <- function(rf.tune.obj, tune.plot , modeltuneplot.file, rocplot.file, width, height, seed){
  if(tune.plot){
    # pdf(modeltuneplot.file)
    plot(rf.tune.obj$fit, ylim=c(0, 1))
    graph2ppt(modeltuneplot.file, width, height)
  }
  set.seed(seed)
  plotROCwithCI(rf.tune.obj$fit$pred$obs, rf.tune.obj$fit$pred$xl, plot.file = rocplot.file, width, height)
}

#============================================================================================
#Preprocess raw data to get rid of uncontacted nucleotides
#============================================================================================
preProcessData <- function(data.file, varcutoff = NULL){
  
  pdbrna.data <- read.table(data.file, sep="\t", check.names = T, header = T, comment.char = "")
  ### preprocess prediction and feature table ###
  rownames(pdbrna.data) <- paste(pdbrna.data$protein,pdbrna.data$PDB,pdbrna.data$idstr, sep = "_")
  pdbrna.data.clean <- subset(pdbrna.data, subset = is.na(pdbrna.data$excluded), select = c(grep("XL", colnames(pdbrna.data)), grep("_", colnames(pdbrna.data))))
  # rownames(pdbrna.data.clean) <- paste(pdbrna.data$protein,pdbrna.data$PDB,pdbrna.data$idstr, sep = "_")
  pdbrna.data.clean$XL[is.na(pdbrna.data.clean$XL)] <- 0
  colnames(pdbrna.data.clean)[1] <- "group"
  pdbrna.data.clean$group <- as.factor(pdbrna.data.clean$group)
  xl.idx <- which(pdbrna.data.clean$group==1)
  notxl.idx <- which(pdbrna.data.clean$group==0)
  levels(pdbrna.data.clean$group) <- c("xl", "nonxl")
  pdbrna.data.clean$group[xl.idx] <- "xl"
  pdbrna.data.clean$group[notxl.idx] <- "nonxl"
  
  ### keep samples interacting with AA ###
  aa.interaction.idx <- grep("interact", colnames(pdbrna.data.clean))
  interact.nt.idx <- which(!apply(pdbrna.data.clean[,aa.interaction.idx], 1, function(x){return(all(x==0))}))
  pdbrna.data.clean <- pdbrna.data.clean[interact.nt.idx,]
  if(!is.null(varcutoff)){
    pdbrna.data.clean <- filterbyvar(pdbrna.data.clean, var.cutoff = varcutoff)$dat
  }
  list(rawData = pdbrna.data, processedData = pdbrna.data.clean, interact.nt.idx = which(is.na(pdbrna.data$excluded))[interact.nt.idx])
  
}

#============================================================================================
#Calculate GLM weight for samples 
#============================================================================================

calculateGLMweight <- function(pdbrna.data.clean, cv.folds,repeat.times, SMOTE= TRUE, ...){
  if(!SMOTE){
    fitControl <- trainControl(method = "repeatedcv", number = cv.folds,repeats = repeat.times, savePredictions = "final", summaryFunction = prSummary, classProbs = T)
  }else{
    fitControl <- trainControl(method = "repeatedcv", number = cv.folds,repeats = repeat.times, savePredictions = "final", sampling = "smote", summaryFunction = prSummary, classProbs = T )
  }
  # fitCtrl <- trainControl(method = "repeatedcv" ,number = cv.folds, savePredictions = T,  classProbs = T,sampling = "smote",verboseIter = F, preProcOptions = c("center", "scale"), summaryFunction = prSummary)
  glm.fit <- train(group ~., data= pdbrna.data.clean, method = "glmnet", trControl = fitControl,...)
  glm.feature.coef <- as.matrix(coef(glm.fit$finalModel, s= glm.fit$bestTune$lambda))[-1,]
  list(glm.fit = glm.fit, glm.feature.coef = glm.feature.coef)
}

#============================================================================================
#scatter plot to combine feature importance, GLM weight and permutation p value
#============================================================================================
plotFeature <- function(feature.data, feature.fill.color = NULL, featurerankplot.file, width, height){
  
  if(!is.null(feature.fill.color)){
    g <- ggplot(data = feature.data, aes(x = MeanDecreaseGini, y = coef, color=group, size=1/pvalue ))
  }else{
    g <- ggplot(data = feature.data, aes(x = MeanDecreaseGini, y = coef, size=1/pvalue ))
  }
  g  <- g + geom_point() +
    ylab("GLM Coefficient")+
    scale_size_continuous(trans="log10", range=c(1,8)) +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 1, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 1, linetype = "solid", colour = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  if(!is.null(feature.fill.color)){
    g <- g + scale_color_manual(values = feature.fill.color)
  }else{
    g <- g
  }
  # geom_text(aes(label = rownames(feature.data)), size=3 ,show_guide=F, angle=30)
  graph2ppt(g, file = featurerankplot.file, width=width, height=height)
}

plotFeature.volcano <- function(feature.data, feature.fill.color = NULL, featurerankplot.file, width, height){
  
  if(!is.null(feature.fill.color)){
    g <- ggplot(data = feature.data, aes(x = sign(coef) * MeanDecreaseGini, y = -log(pvalue) , color=group))
  }else{
    g <- ggplot(data = feature.data, aes(x = sign(coef) * MeanDecreaseGini, y = -log(pvalue), color=group))
  }
  g  <- g + geom_point() +
    ylab("Feature robustness")+
    xlab("Gini decrease with direction")+
    #scale_size_continuous(trans="log10", range=c(1,8)) +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 1, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 1, linetype = "solid", colour = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
    #+geom_text(aes(label = rownames(feature.data)))
  
  if(!is.null(feature.fill.color)){
    g <- g + scale_color_manual(values = feature.fill.color)
  }else{
    g <- g
  }

  graph2ppt(g, file = featurerankplot.file, width=width, height=height)
}

#============================================================================================
#calculate AA frequency for features
#============================================================================================
aafreqCalculation <- function(pdbrna.data.clean, write = T, key, aa.freq.file, aa.code){
  nt <- sapply(rownames(pdbrna.data.clean), function(x){
    return(substr(strsplit(x, split = "\\.")[[1]][2], 1,1))
  })
  pdbrna.data.clean$nt <- nt
  aa.freq.group <- matrix(0, nrow = length(aa.code), ncol = 2 * length(names(table(nt))))
  idx <- match(sapply(colnames(pdbrna.data.clean)[grep(key, colnames(pdbrna.data.clean))], function(x){return(rev(strsplit(x, split = "_")[[1]])[1])}), aa.code)
  
  for (i in 1:length(names(table(nt)))){
    print(i)
    xl.interact <- pdbrna.data.clean[which(pdbrna.data.clean$group=="xl" & pdbrna.data.clean$nt == names(table(nt))[i]),grep(key, colnames(pdbrna.data.clean))]
    notxl.interact <- pdbrna.data.clean[which(pdbrna.data.clean$group=="nonxl" & pdbrna.data.clean$nt == names(table(nt))[i]),grep(key, colnames(pdbrna.data.clean))]
    ifinteract <- function(x){
      return(length(which(x>=1)))
    }
    tmp <- cbind(apply(xl.interact, 2, ifinteract), apply(notxl.interact, 2,ifinteract))
    # colnames(aa.freq.group)[c(2*i-1, 2*i)] <- c(paste("xl", i, sep="_"), paste("nonxl", i, sep="_"))
    aa.freq.group[idx,c(2*i-1, 2*i)] <- tmp
  }
  
  # aa.freq.group <- aa.freq.group[, c(1,3,5,7,2,4,6,8)]

  aa.freq.group <- aa.freq.group[, c(1,3,5,7,2,4,6,8)]
  colnames(aa.freq.group) <- c(paste("xl", names(table(nt)), sep = "_"), paste("nonxl", names(table(nt)), sep = "_"))
  aa.freq.group <- cbind(aa.freq.group, xl=apply(aa.freq.group[,1:4], 1, sum), nonxl=apply(aa.freq.group[,5:8], 1, sum))
  rownames(aa.freq.group) <- paste(key, aa.code, sep="_")

  if(write){
    write.table(aa.freq.group, aa.freq.file , quote = F, sep="\t")
    
  }
  aa.freq.group
}
conformfreqCalculation <- function(cleandata, idx){
  nt <- sapply(rownames(cleandata), function(x){
    return(substr(strsplit(x, split = "\\.")[[1]][2], 1,1))
  })
  cleandata$nt <- nt
  ifinteract <- function(x){
    return(length(which(x>0)))
  }
  res <- matrix(0, nrow = length(names(table(nt))), ncol= 2)
  rownames(res) <- names(table(nt))
  colnames(res) <- paste(colnames(cleandata)[idx], c("xl", "nonxl"), sep = "_")
  
  for (i in 1:length(names(table(nt)))){
    xl.interact <- cleandata[which(cleandata$group=="xl" & cleandata$nt == names(table(nt))[i]),idx]
    notxl.interact <- cleandata[which(cleandata$group=="nonxl" & cleandata$nt == names(table(nt))[i]),idx]
    # print(xl.interact)
    ifinteract <- function(x){
      return(length(which(x>=1)))
    }
   res[i, 1] <- ifinteract(xl.interact)
   res[i, 2] <- ifinteract(notxl.interact)
  }
  return(res)
}

#============================================================================================
summarizeFeature <- function(glm.feature.coef, rf.fit, rf.permute){
  feature.group <- unlist(lapply(names(glm.feature.coef), function(x) {return(paste(strsplit(x,"_")[[1]][1:length(unlist(strsplit(x,"_")))-1],collapse = "_"))}))
  feature.group <- gsub("^puckering_C[2,3]", "puckering", feature.group)
  # feature.group[grep("hbond",feature.group)] <- "hbond"
  feature.data <- data.frame(imp=importance(rf.fit$finalModel), coef=-glm.feature.coef, pvalue= rf.permute$rank.pvalue, group=feature.group)
  feature.data
}

#============================================================================================
summarizePrediction <- function(rf.fit, pdbrna.data, interact.nt.idx){
  rf.best.prediction <- rf.fit$pred
  return(cbind(pdbrna.data[interact.nt.idx, 2:8], rf.best.prediction[order(rf.best.prediction$rowIndex), 1:6]))
}

#============================================================================================
summarizeSample <- function(pdbrna.data.clean, write = T, summary.file){
  rbp <- unique(t(as.data.frame(strsplit(rownames(pdbrna.data.clean), split = "_")))[ ,1])
  rbd <- unique(t(as.data.frame(strsplit(rbp, split = " ")))[ ,2])
  xl.idx <- which(pdbrna.data.clean$group == "xl")
  nt <- as.character(gsub('[0-9]+', '', t(as.data.frame(strsplit(rownames(pdbrna.data.clean), split='\\.')))[,2]))
  nt[nt == "GMP"] <- "G"
  nt.list <- c("A", "C", "U", "G")
  
  if(write){
    con <- file(summary.file, "w")
    cat(paste("XL # :", table(pdbrna.data.clean$group)[1], sep = " "), "\n", file=con)
    cat(paste("nonXL # :", table(pdbrna.data.clean$group)[2], sep = " "), "\n",file=con, append = T)
    cat(paste("RBP # :", length(unique(t(as.data.frame(strsplit(rbp, split = " ")))[ ,1])), sep = " "), "\n",file=con, append = T)
    cat(paste("PDBid # :", length(unique(t(as.data.frame(strsplit(rownames(pdbrna.data.clean), split = "_")))[ ,2])), sep = " "), "\n",file=con, append = T)
    cat(paste("RRM # :", length(grep("RRM", rbd)), sep = " "), "\n",file=con, append = T)
    cat(paste("KH # :", length(grep("KH", rbd)), sep = " "), "\n",file=con, append = T)
    cat("XL nt composition", "\n",file=con, append = T)
    for (i in 1:length(nt.list)){
      cat(paste(nt.list[i], ":", length(which(nt[xl.idx] == nt.list[i])), sep =" "), "\n", file = con, append = T)
    }
    cat("nonXL nt composition", "\n",file=con, append = T)
    for (i in 1:length(nt.list)){
      cat(paste(nt.list[i], ":", length(which(nt[-xl.idx] == nt.list[i])), sep =" "), "\n", file = con, append = T)
    }
    close(con)
  }
  
}
#============================================================================================


rf_twolevel <- function(data, aa.code, mtry, ntree, mtry.seq, smote = T, perc.over, perc.under,...){
  if(smote){
    # set.seed(10)
    data.SMOTE <- SMOTE(group ~ ., data, perc.over = perc.over, perc.under =perc.under)
    data <- data.SMOTE
  }
  rf.firstlevel.obj <- list()
  new.data <- NULL
  for (j in 1:length(aa.code)){
    print(aa.code[j])
    aa.idx <- grep(aa.code[j], colnames(data))
    
    rf.firstlevel.obj[[aa.code[j]]] <- randomForest_feature_construction(data, aa.idx = aa.idx, mtry.seq = mtry.seq, ... )
    new.data <- cbind(new.data, rf.firstlevel.obj[[aa.code[j]]]$new.train.data)
  }
  new.data[new.data == "xl"] <- 1
  new.data[new.data == "nonxl"] <- 0
  colnames(new.data) <- aa.code
  new.data <- cbind(data[, 1:13], new.data)
  rf.twolevel.obj <- randomForest(as.matrix(new.data[,-1]), new.data[,1], localImp = T, mtry= mtry, ntree = ntree)
  localimp <- t(rf.twolevel.obj$localImportance)
  list(rf.firstlevel.obj = rf.firstlevel.obj, rf.twolevel.obj = rf.twolevel.obj, localimp = localimp, new.feature.mat = new.data)
}
#============================================================================================

rf_twolevel_cv_tune <- function(data, mtry.twolevel.seq, ntree.twolevel.seq,metric, ...){
  res.list <- list()
  acc.res <- matrix(0, ncol=24, nrow=length(mtry.twolevel.seq) * length(ntree.twolevel.seq))
  colnames(acc.res) <- c("mtry", "ntree", "acc", "nonxl", "xl", c("Accuracy", "Kappa", "AccuracyLower","AccuracyUpper","AccuracyNull","AccuracyPValue","McnemarPValue","Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision","Recall","F1","Prevalence","Detection Rate","Detection Prevalence","Balanced Accuracy"), "AUC")
  acc.res <- as.data.frame(acc.res)
  for (i in 1:length(mtry.twolevel.seq)){
    for (j in 1:length(ntree.twolevel.seq)){
      # set.seed(seed)
      print(mtry.twolevel.seq[i])
      print(ntree.twolevel.seq[j])
      tmp <- rf_twolevel_cv(data, mtry=mtry.twolevel.seq[i], ntree=ntree.twolevel.seq[j], ...)
      res.list[[paste("mtry=",mtry.twolevel.seq[i]," - ntree=",ntree.twolevel.seq[j],sep="" )]] <- tmp
      acc.res[(i-1)*length(ntree.twolevel.seq) + j, ] <- c(mtry.twolevel.seq[i], ntree.twolevel.seq[j], tmp$acc, tmp$acc.group, tmp$evaluation, tmp$auc)
      
    }
  }
  
  res.list$acc <- acc.res
  res.list$besttune <- acc.res[which.max(acc.res[,match(metric, colnames(acc.res))]), 1:2]
  # res.list$besttune <- acc.res[which.max(acc.res$xl), 1:2]
  res.list
}

#============================================================================================

rf_twolevel_cv <- function(data, smote = T, aa.code, cv.folds, mtry, ntree, mtry.seq,perc.over, perc.under, ...){
  
  sample.idx.group <- split(1:nrow(data), data$group)
  cv.idx.group <- list()
  # set.seed(10)
  for (i in 1:length(sample.idx.group)){
    cv.idx.group[[i]] <- sample(1:cv.folds, length(sample.idx.group[[i]]), replace=TRUE)
  }
  cv.idx <- as.integer(rep(0, nrow(data)))
  cv.idx[unlist(sample.idx.group)] <- unlist(cv.idx.group)
  
  # cross validation
  rf.cv <- list()
  rf.pred <- character(nrow(data))
  names(rf.pred) <- rownames(data)
  pred.prob <- matrix(0, nrow = nrow(data), ncol = 3)
  colnames(pred.prob) <- c("obs", levels(data$group))# rowIndex <- NULL
  rownames(pred.prob) <- rownames(data)
  pred.prob <- as.data.frame(pred.prob)
  pred.prob$obs <- data$group
  
  for (i in 1:cv.folds){
    print(i)
    dat.train <- data[cv.idx != i, ]
    dat.test <- data[cv.idx == i, ]
    if(smote){
      # set.seed(10)
      dat.train.SMOTE <- SMOTE(group ~ ., dat.train, perc.over = perc.over, perc.under = perc.under)
      dat.train <- dat.train.SMOTE
    }
    
    cur.train.sub <- list()
    new.train.data <- NULL
    new.test.data <- NULL
    
    ## in each fold, determine the best parameters using CV in randomforest model
    for (j in 1:length(aa.code)){
      print(aa.code[j])
      aa.idx <- grep(aa.code[j], colnames(dat.train))
      cur.train.sub[[aa.code[j]]] <- randomForest_feature_construction(dat.train, dat.test, aa.idx,mtry.seq = mtry.seq, ... )
      new.train.data <- cbind(new.train.data, cur.train.sub[[aa.code[j]]]$new.train.data)
      new.test.data <- cbind(new.test.data, cur.train.sub[[aa.code[j]]]$new.test.data)
      
    }
    
    ## get new feature data 
    new.train.data[new.train.data == "xl"] <- 1
    new.train.data[new.train.data == "nonxl"] <- 0
    new.test.data[new.test.data == "xl"] <- 1
    new.test.data[new.test.data == "nonxl"] <- 0
    colnames(new.train.data) <- aa.code
    new.train.data <- cbind(dat.train[, 1:13], new.train.data)
    colnames(new.test.data) <- aa.code
    new.test.data <- cbind(dat.test[,1:13], new.test.data)
    
    ## train new feature data on second level
    cur.train.leveltwo <- randomForest(as.matrix(new.train.data[,-1]), new.train.data[,1], mtry= mtry, ntree=ntree)
    cur.pred <- predict(cur.train.leveltwo, as.matrix(new.test.data[,-1]))
    cur.pred.prob <- predict(cur.train.leveltwo, as.matrix(new.test.data[,-1]), type="prob")
    rf.pred[cv.idx == i] <- as.character(cur.pred)
    pred.prob[cv.idx == i, -1] <- cur.pred.prob
    rf.cv[[i]] <- cur.train.leveltwo
  }	
  
  
  acc <- sum(as.character(data$group) == rf.pred) / nrow(data)
  acc.group <- split(as.character(data$group) == rf.pred, as.character(data$group))
  acc.group <- sapply(acc.group, function(x) sum(x)/length(x))
  
  auc <- roc(pred.prob$obs, pred.prob$xl)$auc[1]
  ## calculate precesion and recall
  cm <- confusionMatrix(factor(rf.pred, levels=levels(data$group)), data$group)
  evaluation <- c(cm$overall, cm$byClass)
  list(rf.cv= rf.cv, rf.pred = rf.pred, rf.pred.prob = pred.prob, cv.idx = cv.idx, acc=acc, acc.group=acc.group, auc = auc, evaluation = evaluation)
}


#============================================================================================

randomForest_feature_construction <- function(dat.train, dat.test = NULL, aa.idx, mtry.seq, ...){
  
  res <- list()
  train.data <- dat.train[, c(1, aa.idx )]
  if(! is.null(dat.test)){
    test.data <- dat.test[, aa.idx]
  }
  if(max(mtry.seq) > length(aa.idx)){
    mtry.seq <- 1:length(aa.idx)
  }else{
    mtry.seq <- mtry.seq
  }
  cur.pred.sub <- randomforest_cv_tune(train.data,mtry.seq = mtry.seq,...)
  ## construct new feature set using the best tune model 
  new.train.data <- as.character(predict(cur.pred.sub$finalmodel, train.data[, -1]))
  if(! is.null(dat.test)){
    new.test.data <- as.character(predict(cur.pred.sub$finalmodel, test.data))
  }
  res$model <- cur.pred.sub
  res$new.train.data <- new.train.data
  if(! is.null(dat.test)){
    res$new.test.data <- new.test.data
  }
  res
}
#============================================================================================


randomforest_cv_tune <- function(data, mtry.seq, ntree.seq, seed, metric = "acc",...){
  res.list <- list()
  if( metric != "acc"){
    acc.res <- matrix(0, ncol=24, nrow=length(mtry.seq) * length(ntree.seq))
    colnames(acc.res) <- c("mtry", "ntree", "acc", "nonxl", "xl", c("Accuracy", "Kappa", "AccuracyLower","AccuracyUpper","AccuracyNull","AccuracyPValue","McnemarPValue","Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision","Recall","F1","Prevalence","Detection Rate","Detection Prevalence","Balanced Accuracy"), "AUC")
  }else{
    acc.res <- matrix(0, ncol=5, nrow=length(mtry.seq) * length(ntree.seq))
    colnames(acc.res) <- c("mtry", "ntree", "acc", "nonxl", "xl")
  }
  
  if(is.na(match(metric,colnames(acc.res)))){
    stop("Metric is not existed")
  }
  acc.res <- as.data.frame(acc.res)
  for (i in 1:length(mtry.seq)){
    for (j in 1:length(ntree.seq)){
      set.seed(seed)
      tmp <- randomforest_cv(data, mtry=mtry.seq[i], ntree=ntree.seq[j], ...)
      res.list[[paste("mtry=",mtry.seq[i]," - ntree=",ntree.seq[j],sep="" )]] <- tmp
      if( metric != "acc"){
        acc.res[(i-1)*length(ntree.seq) + j, ] <- c(mtry.seq[i], ntree.seq[j], tmp$acc, tmp$acc.group, tmp$evaluation, tmp$auc)
      }else{
        acc.res[(i-1)*length(ntree.seq) + j, ] <- c(mtry.seq[i], ntree.seq[j], tmp$acc, tmp$acc.group)
      }
    }
  }
  res.list$acc <- acc.res
  res.list$besttune <- acc.res[which.max(acc.res[,match(metric, colnames(acc.res))]), 1:2]
  res.list$finalmodel <- randomForest(as.matrix(data[,-1]), data[,1], mtry= res.list$besttune$mtry, ntree=res.list$besttune$ntree)
  res.list
}

#============================================================================================

randomforest_cv <- function(gbm.dat, cv.folds=5, verbose=TRUE, SMOTE=F, calAUC = F,...)
{
  res <- list()
  # split samples
  sample.idx.group <- split(1:nrow(gbm.dat), gbm.dat$group)
  cv.idx.group <- list()
  # set.seed(10)
  for (i in 1:length(sample.idx.group)){
    cv.idx.group[[i]] <- sample(1:cv.folds, length(sample.idx.group[[i]]), replace=TRUE)
  }
  
  cv.idx <- as.integer(rep(0, nrow(gbm.dat)))
  cv.idx[unlist(sample.idx.group)] <- unlist(cv.idx.group)
  #gbm.dat.cv <- cbind(data.frame(cv.idx=cv.idx), gbm.dat)
  
  # cross validation
  rl.gbm.cv <- list()
  gbm.pred <- character(nrow(gbm.dat))
  pred.prob <- matrix(0, nrow = nrow(gbm.dat), ncol = 3)
  colnames(pred.prob) <- c("obs", levels(gbm.dat$group))# rowIndex <- NULL
  rownames(pred.prob) <- rownames(gbm.dat)
  pred.prob <- as.data.frame(pred.prob)
  pred.prob$obs <- gbm.dat$group
  for (i in 1:cv.folds){
    # cat('fold:',i,'\r')
    gbm.dat.train <- gbm.dat[cv.idx != i, ]
    # gbm.dat.test <- gbm.dat[cv.idx == i, ]
    ## select features that has var larger than 0.05
    
    if(SMOTE){
      print("enter smote")
      gbm.dat.train.SMOTE <- SMOTE(group ~ ., gbm.dat.train, perc.over = 1000, perc.under =100)
      gbm.dat.train <- gbm.dat.train.SMOTE
    }
    gbm.dat.test <- gbm.dat[cv.idx == i, ]
    cur.train <- randomForest(as.matrix(gbm.dat.train[,-1]), gbm.dat.train[,1], ...)
    cur.pred <- predict(cur.train, as.matrix(gbm.dat.test[,-1]))
    gbm.pred[cv.idx == i] <- as.character(cur.pred)
    rl.gbm.cv[[i]] <- cur.train
    cur.pred.prob <- predict(cur.train, as.matrix(gbm.dat.test[,-1]), type="prob")
    pred.prob[cv.idx == i, -1] <- cur.pred.prob
  }	
  # cat('fold:',i,'\n')
  acc <- sum(as.character(gbm.dat$group) == gbm.pred) / nrow(gbm.dat)
  
  acc.group <- split(as.character(gbm.dat$group) == gbm.pred, as.character(gbm.dat$group))
  acc.group <- sapply(acc.group, function(x) sum(x)/length(x))
  
  if(calAUC){
    auc <- roc(pred.prob$obs, pred.prob$xl)$auc[1]
    ## calculate precesion and recall
    cm <- confusionMatrix(factor(gbm.pred, levels=levels(gbm.dat$group)), gbm.dat$group)
    evaluation <- c(cm$overall, cm$byClass)
    res <- list(rl.rf.cv= rl.gbm.cv, rf.pred = gbm.pred, rf.pred.prob = pred.prob, cv.idx = cv.idx, acc=acc, acc.group=acc.group, evaluation=evaluation, auc = auc)	
  }else{
    res <- list(rl.rf.cv= rl.gbm.cv, rf.pred = gbm.pred, rf.pred.prob = pred.prob, cv.idx = cv.idx, acc=acc, acc.group=acc.group)
  }
  
  res
}
#============================================================================================


gbm_twolevel_cv_tune <- function(data, mtry.twolevel.seq, ntree.twolevel.seq, ...){
  res.list <- list()
  acc.res <- matrix(0, ncol=5, nrow=length(mtry.twolevel.seq) * length(ntree.twolevel.seq))
  colnames(acc.res) <- c("mtry", "ntree", "acc", "nonxl", "xl")
  acc.res <- as.data.frame(acc.res)
  for (i in 1:length(mtry.twolevel.seq)){
    for (j in 1:length(ntree.twolevel.seq)){
      # set.seed(seed)
      print(mtry.twolevel.seq[i])
      print(ntree.twolevel.seq[j])
      tmp <- gbm_twolevel_cv(data, mtry=mtry.twolevel.seq[i], ntree=ntree.twolevel.seq[j], ...)
      res.list[[paste("mtry=",mtry.twolevel.seq[i]," - ntree=",ntree.twolevel.seq[j],sep="" )]] <- tmp
      acc.res[(i-1)*length(ntree.twolevel.seq) + j, ] <- c(mtry.twolevel.seq[i], ntree.twolevel.seq[j], tmp$acc, tmp$acc.group)
      
    }
  }
  res.list$acc <- acc.res
  res.list$besttune <- acc.res[which.max(acc.res$xl), 1:2]
  res.list
}
#============================================================================================


gbm_twolevel_cv <- function(gbm.dat, smote= T, cv.folds=5, verbose=TRUE, aa.code, mtry, ntree,perc.over, perc.under, ...)
{
  # set.seed(30)
  # split samples
  sample.idx.group <- split(1:nrow(gbm.dat), gbm.dat$group)
  cv.idx.group <- list()
  for (i in 1:length(sample.idx.group)){
    cv.idx.group[[i]] <- sample(1:cv.folds, length(sample.idx.group[[i]]), replace=TRUE)
  }
  
  cv.idx <- as.integer(rep(0, nrow(gbm.dat)))
  cv.idx[unlist(sample.idx.group)] <- unlist(cv.idx.group)
  #gbm.dat.cv <- cbind(data.frame(cv.idx=cv.idx), gbm.dat)
  
  # cross validation
  rf.cv <- list()
  rf.pred <- character(nrow(gbm.dat))
  names(rf.pred) <- rownames(gbm.dat)
  pred.prob <- matrix(0, nrow = nrow(gbm.dat), ncol = 3)
  colnames(pred.prob) <- c("obs", levels(gbm.dat$group))# rowIndex <- NULL
  rownames(pred.prob) <- rownames(gbm.dat)
  pred.prob <- as.data.frame(pred.prob)
  pred.prob$obs <- gbm.dat$group
  
  
  for (i in 1:cv.folds){
    print(i)
    gbm.dat.train <- gbm.dat[cv.idx != i, ]
    gbm.dat.test <- gbm.dat[cv.idx == i, ]
    new.train.data <- NULL
    new.test.data <- NULL
    if(smote){
      # set.seed(10)
      gbm.dat.train.SMOTE <- SMOTE(group ~ ., gbm.dat.train, perc.over = perc.over, perc.under = perc.under)
      gbm.dat.train <- gbm.dat.train.SMOTE
    }
    
    ### train first level
    for (j in 1:length(aa.code)){
      
      print(aa.code[j])
      aa.idx <- grep(aa.code[j], colnames(gbm.dat.train))
      train.data <- gbm.dat.train[, c(1, aa.idx )]
      test.data <- gbm.dat.test[, aa.idx]
      
      # first level using gbm
      gbmGrid <-  expand.grid(interaction.depth = seq(1,5,2), 
                              n.trees = 10^seq(1,3), 
                              shrinkage = 0.1,
                              n.minobsinnode = 1:3)
      ## ntree 100, minbosinnode 5 interaction.depth 1
      cur.pred.sub <- caret.wrap(train.data, method="gbm", cv.folds = 5, repeat.times = 1, SMOTE = F, tuneGrid=gbmGrid, metric = "ROC")
      
      #print("done")
      
      new.train.data <- cbind(new.train.data,   cur.pred.sub$fit$pred$xl[match(1:dim(train.data)[1], cur.pred.sub$fit$pred$rowIndex)])
      new.test.data <- cbind(new.test.data, predict(cur.pred.sub$fit$finalModel, as.matrix(test.data), n.trees = cur.pred.sub$fit$bestTune$n.trees, type = "response"))
      
    }
    
    new.train.data[new.train.data >= 0.5] <- 1
    new.train.data[new.train.data < 0.5] <- 0
    colnames(new.train.data) <- aa.code
    new.train.data <- cbind(gbm.dat.train[, 1:13], new.train.data)
    
    ## glmnet give the smaller prob as group xl
    new.test.data[new.test.data < 0.5] <- 0
    new.test.data[new.test.data >= 0.5] <- 1
    colnames(new.test.data) <- aa.code
    new.test.data <- cbind(gbm.dat.test[,1:13], new.test.data)
    # 
    # rfGrid <- expand.grid(.mtry= seq(1,5), .ntree=c(10,100,1000))
    # cur.train <- caret.wrap(new.train.data, method=customRF, cv.folds = 5, repeat.times = 1, SMOTE = T, tuneGrid=rfGrid, metric = "ROC")
    # 
    # cur.pred <- predict( cur.train$fit$finalModel, new.test.data[, -1])
    # 
    cur.train.leveltwo <- randomForest(as.matrix(new.train.data[,-1]), new.train.data[,1], mtry= mtry, ntree=ntree,...)
    cur.pred <- predict(cur.train.leveltwo, as.matrix(new.test.data[,-1]))
    cur.pred.prob <- predict(cur.train.leveltwo, as.matrix(new.test.data[,-1]), type="prob")
    rf.pred[cv.idx == i] <- as.character(cur.pred)
    pred.prob[cv.idx == i, -1] <- cur.pred.prob
    rf.cv[[i]] <- cur.train.leveltwo
    
  }
  
  cat('fold:',i,'\n')
  acc <- sum(as.character(gbm.dat$group) == rf.pred) / nrow(gbm.dat)
  auc <- auc(roc(gbm.dat$group, pred.prob$xl))
  acc.group <- split(as.character(gbm.dat$group) == rf.pred, as.character(gbm.dat$group))
  acc.group <- sapply(acc.group, function(x) sum(x)/length(x))
  
  list(rf.cv= rf.cv, rf.pred = rf.pred, rf.pred.prob = pred.prob, cv.idx = cv.idx, acc=acc, acc.group=acc.group, auc=auc)	
}
#============================================================================================

boxplot.localimp <- function(localimp, xl.idx, aa.color.code= NULL, outfile = NULL,...){
  plot.data <- data.frame()
  
  for (i in 1:dim(localimp)[2]){
    
    plot.data <- rbind(plot.data, data.frame(group=paste(colnames(localimp)[i],"nonxl",  sep="."), Prob=localimp[-xl.idx, i], xl="nonxl"))
    plot.data <- rbind(plot.data, data.frame(group=paste(colnames(localimp)[i],"xl", sep="."), Prob=localimp[xl.idx, i], xl="xl"))
  }
  
    ### plot the probality median distribution to show the most informative amino acid
    g <- ggplot(plot.data, aes(x = group, y = Prob, fill = group)) +
      # geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
      geom_boxplot(size = 0.1, outlier.shape = NA)+
      # scale_fill_manual(values = as.vector(sapply(aa.color.code, function(x){return(rep(x,2))}))) +
      # scale_x_discrete(limits = conf.order.new) +
      # geom_boxplot(width=0.1, fill="grey")+
      # scale_color_manual(values = c("grey", "black"))+
      theme_classic() +
      theme(legend.position = "none") +
      # theme(
      #   axis.title = element_blank(),
      #   axis.text = element_blank(),
      #   axis.ticks = element_blank(),
      #   axis.line = element_blank()
      # ) +
      theme(axis.text.x = element_text(angle=45, size=8))+
      # scale_x_discrete(breaks = labels=colnames(new.feature.data)[14:33])+
      stat_summary(fun.y = median,
                   geom = "point",
                   size = 1)
    if(!is.null(outfile)){
      graph2ppt(g,
                file = outfile,
                width = width,
                height = height,...)
      
    }
    
 g
}
wilcox.test.localimp <- function(p.aa.fulldata, xl.idx){
  pvalue <- NULL
  for (i in 1:dim(p.aa.fulldata)[2]){
    print(colnames(p.aa.fulldata)[i])
    pvalue <- c(pvalue, wilcox.test(p.aa.fulldata[xl.idx, i], p.aa.fulldata[-xl.idx, i])$p.value)
  }
  names(pvalue) <- colnames(p.aa.fulldata)
  pvalue
}


extendrawdata <- function(rawdata){
  ### get the upstream and downstream nt
  synoranti <- function(x){
    return(strsplit(colnames(x)[which.max(x)], split = "_")[[1]][3])
  }
  c2orc3 <- function(x){
    return(strsplit(colnames(x)[which.max(x)], split = "_")[[1]][2])
  }
  for (i in unique(rawdata$PDB)){
    idx <- which(rawdata$PDB == i)
    print(length(idx))
    for (j in idx){
      if(j == min(idx)){
        rawdata$upstreamNt[j] <- toupper(as.character(rawdata$ntcode[j]))
        rawdata$upstreamNtConf[j] <- NA
        rawdata$upstreamNtPuckering[j] <- NA
        rawdata$dnstreamNt[j] <- paste(toupper(as.character(rawdata$ntcode[j])),toupper(as.character(rawdata$ntcode)[j+1]),sep="")
        rawdata$dnstreamNtConf[j] <- synoranti(rawdata[j+1, 12:14])
        rawdata$dnstreamNtPuckering[j] <- c2orc3(rawdata[j+1, 15:17])
        
      }else if(j == max(idx)){
        rawdata$upstreamNt[j] <- paste(toupper(as.character(rawdata$ntcode[j-1])),toupper(as.character(rawdata$ntcode)[j]),sep="")
        rawdata$upstreamNtConf[j] <- synoranti(rawdata[j-1, 12:14])
        rawdata$upstreamNtPuckering[j] <- c2orc3(rawdata[j-1, 15:17])
        rawdata$dnstreamNt[j] <- toupper(as.character(rawdata$ntcode[j]))
        rawdata$dnstreamNtConf[j] <- NA
        rawdata$dnstreamNtPuckering[j] <- NA
      }else{
        rawdata$upstreamNt[j] <-  paste(toupper(as.character(rawdata$ntcode[j-1])),toupper(as.character(rawdata$ntcode)[j]),sep="")
        rawdata$upstreamNtConf[j] <- synoranti(rawdata[j-1, 12:14])
        rawdata$upstreamNtPuckering[j] <- c2orc3(rawdata[j-1, 15:17])
        rawdata$dnstreamNt[j] <- paste(toupper(as.character(rawdata$ntcode[j])),toupper(as.character(rawdata$ntcode)[j+1]),sep="")
        rawdata$dnstreamNtConf[j] <- synoranti(rawdata[j+1, 12:14])
        rawdata$dnstreamNtPuckering[j] <- c2orc3(rawdata[j+1, 15:17])
        
      }
      
    }
    ### get the RBP domain  
    sample.domain <- unlist(lapply(rownames(rawdata), function(x) {return(unlist(strsplit(x, "_"))[1])}))
    
    sample.domain[grep("RRM", sample.domain)] <- "RRM"
    sample.domain[grep("Zn", sample.domain)] <- "ZnF"
    sample.domain[grep("KH", sample.domain)] <- "KH"
    sample.domain[grep("let-7g", sample.domain)] <- "CSD/ZKD"
    sample.domain[grep("STAR", sample.domain)] <- "STAR"
    sample.domain[grep("QKI", sample.domain)] <- "STAR"
    rawdata$domain <- sample.domain
    rawdata$nt <- as.factor(toupper(rawdata$ntcode))
  }
  rawdata
}