library(ape)
library("seqinr")
library(micropan) 
library(phangorn)
library('ggtree')
library('treeio')
library(gtable)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
###############################小型矩阵#############################
data <- read.table("path/allhunhegnewpdb.tab",sep = "\t", header=T)
head(data)
nrow(data)
ncol(data)
countData <- data[2:49]
rownames(countData) <-data$id
head(countData)
data1 <- as.matrix(countData)


cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(
    round(data1[i, j], 1), 
    x, y,
    gp = gpar(
      fontsize = 3
    ))
}

#rect_gp = gpar(type = "none")
cell_fun
col_fun = colorRamp2(c(0, 8, 25), c("red", "white", "blue"))
g<-Heatmap(
  data1, name = "distance", 
  col = col_fun,
  #cell_fun = cell_fun,
  width = unit(20, "cm"),
  height = unit(20, "cm"),
  #column_title = "sample",
  #column_title_side = "bottom",
  #column_title_rot = 0,
  #row_title = "Genes",
  #row_title_side = "left",
  #row_title_rot = 90,
  column_title_gp = gpar(
    col = "black",
    fontsize = 16,
    fontface = "italic",
    fill = "grey",
    border = "black"),
  #row_title_gp = gpar(
  #col = "black",
  #fontsize = 16,
  #fontface = "italic",
  #fill = "grey",
  #border = "black"),
  row_names_rot = 45,
  column_names_rot = 45,
  row_names_gp = gpar(
    col = "black",
    fontsize = 8
  ),
  column_names_gp= gpar(
    col = "black",
    fontsize = 8
  )
  
)

g
####################################小型矩阵热图################

a<-Heatmap(data1, name = "mat", rect_gp = gpar(type = "none"),
           row_dend_width = unit(30, "cm"),
           show_column_dend = FALSE,
           show_column_names = FALSE,
           show_heatmap_legend = FALSE,
           column_title = "wen100")
a
# 计算列之间的距离和聚类
dist_mat = dist(t(data1), method = "euclidean")
hclust_col = hclust(dist_mat, method = "ward.D")

# 创建列树图
col_dend = as.dendrogram(hclust_col)

#dist_mat = dist(mat)
#hclust_tree = hclust(dist_mat)

# 将聚类结果转换为dendrogram对象

df_tree <- as.phylo(col_dend)# 将聚类结果转成系统发育格式
write.tree(phy=df_tree, file="te.nwk") 

#write.tree(phy=a, file="heatmaptree.nwk") # 输出newick格式文件
############################################单独树图##############

MMu<-upgma(data1)
MMw<-wpgma(data1)
MMnj<-NJ(data1)
MMnb<-neighborNet(data1)



par(mfrow=c(1,3),plt=c(0, 0.8, 0.1, 0.9),mar=c(1, 2, 1, 1) + 0.1)

#plot(MMu,main="UPGMA",cex=0.5)#UPGMA法构建树
#plot(MMu,main="UPGMA",cex=0.8)#UPGMA法构建树
plot(MMu,main="UPGMA",cex=1)#UPGMA法构建树
plot(MMw,main="WPGMA",cex=1)#WPGMA法构建树
plot(MMnj,main="NJ",cex=1)#NJ法构建树



df_tree <- as.phylo(MMu)# 将聚类结果转成系统发育格式
write.tree(phy=df_tree, file="3506upgmatree.nwk") # 输出newick格式文件
df_tree <- as.phylo(MMw)# 将聚类结果转成系统发育格式
write.tree(phy=df_tree, file="3506wpgmatree.nwk") # 输出newick格式文件
df_tree <- as.phylo(MMnj)# 将聚类结果转成系统发育格式
write.tree(phy=df_tree, file="3506NJtree.nwk") # 输出newick格式文件




