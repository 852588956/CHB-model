
library(pROC)
library(glmnet)
library(MASS)
library(e1071)
library(xgboost) #用于拟合xgboost模型
library(randomForest)
library(tidyverse)
library(randomForestSRC)
library(gbm)
library(caret)


AUC <- function(data,RS){
  dfroc1 <- roc(data[, 1], RS)
  return(dfroc1$auc[1])
}


####LASSO
set.seed(seed = 123)
my_lasso <- function(training, val, name){
  y <- as.matrix(training[, 1])
  x <- as.matrix(training[, -1])
  
  cv.fit = cv.glmnet(x, y, family = "gaussian", alpha = 1, nfolds = 10)
  fit = glmnet(x = x,y = y, family = "gaussian", alpha = 1, lambda = cv.fit$lambda.min)
  coef.min = coef(cv.fit, s = "lambda.min")  
  active.min = which(coef.min@i != 0) ## 找出那些回归系数没有被惩罚为0的
  lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] ## 提取基因名称
  fit$subFeature <- lasso_geneids
  if(fit$subFeature[1] == "(Intercept)"){
    fit$subFeature <- fit$subFeature[2:length(fit$subFeature)]
  }
  RS <- as.numeric(predict(fit, as.matrix(training[, -1])))
  auc <- AUC(training,RS)
  fit$AUC <- auc
  
  auc <- c(auc)
  for (i in 1:length(val)) {
    data <- val[[i]]
    RS <- as.numeric(predict(fit, as.matrix(data[, -1])))
    au <- AUC(data,RS)
    auc <- c(auc, au)
  }
  fit$name <- name
  fit$auc <- auc
  return(list(fit))
}




#stepAIC
my_stepAIC <- function(training, val, name){
  y <- as.matrix(training[, 1])
  x <- as.matrix(training[, -1])
  direction = c("both", "backward", "forward")
  models <- list()
  for (i in direction) {
    full.model <- lm(Y~., data = training)
    fit <- stepAIC(full.model, direction = i, trace = 0)
    fit$subFeature = names(coef(fit))
    if(fit$subFeature[1] == "(Intercept)"){
      fit$subFeature <- fit$subFeature[2:length(fit$subFeature)]
    }
    
    RS <- as.numeric(predict(fit, training))
    auc <- AUC(training,RS)
    fit$AUC <- auc
    
    
    name2 <- paste0(name, "[", i, "]")
    auc <- c(auc)
    for (j in 1:length(val)) {
      data <- val[[j]]
      y <- as.matrix(data[, 1])
      x <- as.matrix(data[, -1])
      RS <- as.numeric(predict(fit, data))
      au <- AUC(data,RS)
      auc <- c(auc, au)
    }
    fit$name <- name2
    fit$auc <- auc
    
    models[[i]] <- fit
  }
  
  return(models)
}

#SVM
my_svm <- function(training, val, name){
  y <- as.matrix(training[, 1])
  x <- as.matrix(training[, -1])
  fit <- svm(Y~.,data = training)
  fit$subFeature = colnames(x)
  RS <- as.numeric(predict(fit, training[, -1]))
  auc <- AUC(training,RS)
  fit$AUC <- auc
  
  
  auc <- c(auc)
  for (i in 1:length(val)) {
    data <- val[[i]]
    
    RS <- as.numeric(predict(fit, data[, -1]))
    au <- AUC(data, RS)
    auc <- c(auc, au)
  }
  fit$name <- name
  fit$auc <- auc
  return(list(fit))
}



#importance <- xgb.importance(colnames(training[, -1]), model = bst)
#head(importance)
#xgb.ggplot.importance(importance)


#xgBoost
my_xgBoost <- function(training, val, name){
  fit <- xgboost(data = as.matrix(training[, -1]), 
                 label = training[, 1], 
                 max.depth = 6, 
                 eta = 0.5, 
                 nthread = 2, 
                 nrounds = 50,verbose = 2)
  
  fit$subFeature <- colnames(training[, -1])
  RS <- as.numeric(predict(fit, as.matrix(training[, -1])))
  auc <- AUC(training,RS)
  fit$AUC <- auc
  
  auc <- c(auc)
  for (i in 1:length(val)) {
    data <- val[[i]]
    
    RS <- as.numeric(predict(fit, as.matrix(data[, -1])))
    au <- AUC(data,RS)
    auc <- c(auc, au)
  }
  fit$name <- name
  fit$auc <- auc
  return(list(fit))
}





#Rideg
my_Rideg <- function(training, val, name){
  y <- data.matrix(training[, 1])
  x <- training[, -1]
  
  cv.fit = cv.glmnet(x = as.matrix(x), y = y, family = "gaussian", alpha = 0, nfolds = 10)
  fit = glmnet(x = x,y = y, family = "gaussian", alpha = 0, lambda = cv.fit$lambda.min)
  coef.min = coef(cv.fit, s = "lambda.min")  
  active.min = which(coef.min@i != 0) ## 找出那些回归系数没有被惩罚为0的
  lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] ## 提取基因名称
  fit$subFeature <- lasso_geneids
  if(fit$subFeature[1] == "(Intercept)"){
    fit$subFeature <- fit$subFeature[2:length(fit$subFeature)]
  }
  RS <- as.numeric(predict(fit, as.matrix(training[, -1])))
  auc <- AUC(training,RS)
  fit$AUC <- auc
  
  auc <- c(auc)
  for (i in 1:length(val)) {
    data <- val[[i]]
    RS <- as.numeric(predict(fit, as.matrix(data[, -1])))
    au <- AUC(data,RS)
    auc <- c(auc, au)
  }
  fit$name <- name
  fit$auc <- auc
  return(list(fit))
}






#RF
my_RF <- function(training, val, name){
  fit <- randomForest(Y~., data=training, proximity=TRUE,
                      importance= T,
                      forest = T,
                      mtry=1,
                      nodesize=10,
                      ntree=200
                      )
  
  fit$subFeature = rownames(fit$importance)
  
  
  RS <- predict(fit, as.data.frame(training)) 
  auc <- AUC(training,RS)
  fit$AUC <- auc
  
  auc <- c(auc)
  for (i in 1:length(val)) {
    data <- val[[i]]
    RS <- predict(fit, as.data.frame(data)) 
    au <- AUC(data,RS)
    auc <- c(auc, au)
  }
  fit$name <- name
  fit$auc <- auc
  return(list(fit))
}


#Enet
my_Enet <- function(training, val, name){
  y <- data.matrix(training[, 1])
  x <- training[, -1]
  models <- list()
  for (i in (1:9)/10) {
    cv.fit = cv.glmnet(x = as.matrix(x), y = y, family = "gaussian", alpha = i, nfolds = 10)
    fit = glmnet(x = as.matrix(x), y = y, family = "gaussian", alpha = i, lambda = cv.fit$lambda.min)
    coef.min = coef(cv.fit, s = "lambda.min")  
    active.min = which(coef.min@i != 0) ## 找出那些回归系数没有被惩罚为0的
    lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] ## 提取基因名称
    fit$subFeature <- lasso_geneids
    if(fit$subFeature[1] == "(Intercept)"){
      fit$subFeature <- fit$subFeature[2:length(fit$subFeature)]
    }
    RS <- as.numeric(predict(fit, as.matrix(training[, -1])))
    auc <- AUC(training,RS)
    fit$AUC <- auc
    
    name2 <- paste0(name, "[", i, "]")
    auc <- c(auc)
    for (j in 1:length(val)) {
      data <- val[[j]]
      RS <- as.numeric(predict(fit, as.matrix(data[, -1])))
      au <- AUC(data,RS)
      auc <- c(auc, au)
    }
    fit$name <- name2
    fit$auc <- auc
    
    models[[i*10]] <- fit
  }
  return(models)
}




#GBM
my_GBM <- function(training, val, name){
  y <- data.matrix(training[, 1])
  x <- training[, -1]
  
  fit <- gbm(formula = y ~ .,
             data = as.data.frame(x),
             n.trees = 5000,
             interaction.depth = 5,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  best <- which.min(fit$cv.error)
  fit <- gbm(formula = y ~ .,
             data = as.data.frame(x),
             n.trees = best,
             interaction.depth = 5,
             n.minobsinnode = 10,
             shrinkage = 0.001, n.cores = 8)
  fit$subFeature = rownames(summary.gbm(fit, plotit = FALSE))[summary.gbm(fit, plotit = FALSE)$rel.inf > 0]
  
  RS <- as.numeric(predict(fit, as.data.frame(x)))
  auc <- AUC(training,RS)
  fit$AUC <- auc
  
  auc <- c(auc)
  for (i in 1:length(val)) {
    vals <- val[[i]]
    y <- data.matrix(vals[, 1])
    x <- vals[, -1]
    RS <- as.numeric(predict(fit, as.data.frame(x)))
    au <- AUC(vals,RS)
    auc <- c(auc, au)
  }
  
  fit$name <- name
  fit$auc <- auc
  return(list(fit))
}






#KNN
my_knn <- function(training, val, name){
  # 设置10折交叉训练
  control <- trainControl(method = 'cv',number = 10)
  # knn模型训练
  fit <- train(Y~.,training,
                 method = 'knn',
                 preProcess = c('center','scale'),
                 trControl = control,
                 tuneLength = 2)
  
  fit$subFeature <- fit[["coefnames"]]
  RS <- predict(fit,newdata = training)
  auc <- AUC(training,RS)
  fit$AUC <- auc
  
  auc <- c(auc)
  for (i in 1:length(val)) {
    data <- val[[i]]
    RS <- predict(fit,newdata = data)
    au <- AUC(data,RS)
    auc <- c(auc, au)
  }
  fit$name <- name
  fit$auc <- auc
  return(list(fit))
}










SimpleHeatmap <- function(Cindex_mat = NULL, 
                          avg_Cindex = NULL, 
                          CohortCol = NULL, 
                          barCol = NULL,
                          col = c("#4195C1", "#FFFFFF", "#CB5746"), # 红蓝配色
                          cellwidth = 1, cellheight = 0.5, 
                          cluster_columns, cluster_rows, title){
  col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                            col = list("Cohort" = CohortCol),
                            show_annotation_name = F)
  
  row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                            gp = gpar(fill = barCol, col = NA),
                                            add_numbers = T, numbers_offset = unit(-10, "mm"),
                                            axis_param = list("labels_rot" = 0),
                                            numbers_gp = gpar(fontsize = 9, col = "white"),
                                            width = unit(3, "cm")),
                         show_annotation_name = F)
  
  Heatmap(as.matrix(Cindex_mat), name = "AUC",
          row_title = title,
          
          
          right_annotation = row_ha, 
          top_annotation = col_ha,
          col = col, 
          rect_gp = gpar(col = "black", lwd = 1), # 边框设置为黑色
          cluster_columns = cluster_columns, cluster_rows = cluster_rows, # 不进行聚类，无意义
          show_column_names = FALSE, 
          show_row_names = TRUE,
          row_names_side = "left",
          width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
          height = unit(cellheight * nrow(Cindex_mat), "cm"),
          column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
          column_title = NULL,
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                      x, y, gp = gpar(fontsize = 10))
          }
  )
}













