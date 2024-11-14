data <- read.table("testdata.txt", sep = "\t", header = T, row.names = 1)



colnames(data)[1] <- "Y"
# 划分数据集为训练集和测试集
set.seed(123) # 设置随机种子以获得可重复的结果

data <- data[sample(nrow(data)), ]
subsets <- split(data, 1:5)





m <- c("RF","lasso", "stepAIC",
       "svm", "xgBoost", "Rideg",
       "Enet", "knn", "GBM")
#m <- c("Rideg","RSF")
#m <- c("lasso", "superpc")



cindexs <- as.data.frame(matrix(nrow = 2, ncol = 1))
colnames(cindexs) <- m[1]
rownames(cindexs) <- c("train", "test")
cindexs <- list(cindexs,cindexs,cindexs,cindexs,cindexs)


source("ML.R")
library(RColorBrewer)
library(ComplexHeatmap)


for (k in 1:5) {
  training <- subsets[-k]
  training <- rbind(training[[1]], training[[2]],
                    training[[3]], training[[4]])
  val <- list(subsets[[k]])
  
  for(i in m){
    name <- i
    print(name)
    eval(parse(text = paste0("fit1 <- my_", i, "(training, val, name)")))
    for (j in 1:length(fit1)) {
      fit <- fit1[[j]] 
      gene <- fit$subFeature
      for (o in 1:length(gene)) {
        gene[o] <- gsub('[`"`]', '', gene[o])
      }
      name1 <- fit$name
      cindexs[[k]][, fit$name] <- fit$auc
    }
  }
}


tmp <- cindexs[[1]]
tmp[1, ] <- (cindexs[[1]][1, ] + 
  cindexs[[2]][1, ] +
  cindexs[[3]][1, ] +
  cindexs[[4]][1, ] +
  cindexs[[5]][1, ])/5
tmp[2, ] <- (cindexs[[1]][2, ] + 
               cindexs[[2]][2, ] +
               cindexs[[3]][2, ] +
               cindexs[[4]][2, ] +
               cindexs[[5]][2, ])/5
cindexs <- tmp








#save(models, file = "models.rda")
cindexs <- t(cindexs)
write.csv(cindexs, "cindexs_G2S2.csv")

Cindex_mat <- read.csv("cindexs_G2S2.csv", check.names = F, row.names = 1, header = T)
avg_Cindex <- sort(apply(Cindex_mat, 1, mean), decreasing = T) # 计算每种算法在所有队列中平均C-index，并降序排列
Cindex_mat <- Cindex_mat[names(avg_Cindex), ] # 对C-index矩阵排序
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # 保留三位小数
#fea_sel <- fea_list[[rownames(Cindex_mat)[1]]] # 最优模型（即测试集[或者训练集+测试集]C指数均值最大）所筛选的特征


CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired")[c(1,3)]# 设置绘图时的队列颜色
names(CohortCol) <- colnames(Cindex_mat)

# 调用简易绘图函数
cellwidth = 1; cellheight = 0.5
source("ML.R")
hm <- SimpleHeatmap(Cindex_mat = Cindex_mat, # 主矩阵0..
                    avg_Cindex = avg_Cindex, # 侧边柱状图
                    CohortCol = CohortCol, # 列标签颜色
                    barCol = "steelblue", # 右侧柱状图颜色
                    col = c("#1CB8B2", "#FFFFFF", "#EEB849"), # 热图颜色
                    cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F,
                    title = ">=G2/S2") # 是否对行列进行聚类
p1 <- hm
pdf(file.path("G2S2_heatmap of cindex.pdf"), width = cellwidth * ncol(Cindex_mat) + 7, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right") # 热图注释均放在右侧
invisible(dev.off())






