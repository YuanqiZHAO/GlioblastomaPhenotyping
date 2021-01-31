# source_file = list.files('/Users/yuanqizhao/Desktop/SNFtool/R')
# print(sessionInfo())
# install.packages('parallel')
# library(parallel)
source('standardNormalization.R')
source('dist2.R')
source('affinityMatrix.R')
source('calNMI.R')
source('displayClusters.r')
source('estimateNumberOfClustersGivenGraph.R')
source('groupPredict.r')
source('SNF.R')
source('spectralClustering.r')
source('internal.R')
### 25 .5 10 0.0820928696879084
### 30 .5 30 "0.0912106312356116"
### 30 .5 20 "0.0943999579800279"
K = 30;	  ### test 10-30
alpha = 0.5;  
T = 20; 

exp = as.matrix(read.table('/Users/yuanqizhao/Desktop/CompMed_Project/data/expression.txt',header=TRUE,sep='\t',row.names=1))
methy = as.matrix(read.table('/Users/yuanqizhao/Desktop/CompMed_Project/data/methylation.txt',header=TRUE,sep='\t',row.names=1))
cnv = as.matrix(read.table('/Users/yuanqizhao/Desktop/CompMed_Project/data/cnv.txt',header=TRUE,sep='\t',row.names=1))
truelabel = read.table('/Users/yuanqizhao/Desktop/CompMed_Project/data/true_subtype.txt',header=TRUE,sep='\t',row.names=1)
truelabel = truelabel[,1]
# print(length(truelabel))
# ###truelabel = c(matrix(1,100,1),matrix(2,100,1)); ##the ground truth of the simulated data;

# #### Test
# # exp = t(exp[1:6,1:6])
# # methy = t(methy[1:6,1:6])
# # cnv = t(cnv[1:6,1:6])
# # truelabel = truelabel[1:6]
# # print(c('Input',class(exp),class(methy),class(cnv)))

norm_exp = standardNormalization(exp)
write.table(norm_exp, file="../Result/norm_exp.txt", row.names=FALSE, col.names=FALSE)
# norm_methy = standardNormalization(methy)
# write.table(norm_methy, file="../Result/norm_methy.txt", row.names=FALSE, col.names=FALSE)
norm_cnv = standardNormalization(cnv)
write.table(norm_cnv, file="../Result/norm_cnv.txt", row.names=FALSE, col.names=FALSE)
print('Normalization Complete')
# print(norm_exp[1:6,1:6])
# print(methy[1:6,1:6])
# print(norm_cnv[1:6,1:6])


Dist1 = dist2(as.matrix(norm_exp),as.matrix(norm_exp));
write.table(Dist1, file="../Result/Dist1.txt", row.names=FALSE, col.names=FALSE)
Dist2 = dist2(as.matrix(methy),as.matrix(methy));
write.table(Dist2, file="../Result/Dist2.txt", row.names=FALSE, col.names=FALSE)
Dist3 = dist2(as.matrix(norm_cnv),as.matrix(norm_cnv));
write.table(Dist3, file="../Result/Dist3.txt", row.names=FALSE, col.names=FALSE)
print('Distance Complete')
# print(Dist1)
# print(Dist2)
# print(Dist3)

W1 = affinityMatrix(Dist1, K, alpha)
write.table(W1, file="../Result/W1.txt", row.names=FALSE, col.names=FALSE)
# print(c('W1',class(W1)))

W2 = affinityMatrix(Dist2, K, alpha)
# print(c('W2',class(W2)))
write.table(W2, file="../Result/W2.txt", row.names=FALSE, col.names=FALSE)

W3 = affinityMatrix(Dist3, K, alpha)
write.table(W3, file="../Result/W3.txt", row.names=FALSE, col.names=FALSE)
print('Graph Complete')

# W1 = as.matrix(read.table('/Users/yuanqizhao/Desktop/SNFtool/Result/W1.txt',header=FALSE,sep=' '))
# W2 = as.matrix(read.table('/Users/yuanqizhao/Desktop/SNFtool/Result/W2.txt',header=FALSE,sep=' '))
# W3 = as.matrix(read.table('/Users/yuanqizhao/Desktop/SNFtool/Result/W3.txt',header=FALSE,sep=' '))
# print(c(dim(W1),dim(W2),dim(W3)))


displayClusters(W1,truelabel);
displayClusters(W2,truelabel);
displayClusters(W3,truelabel);

W = SNF(list(W1,W2,W3), K, T)
write.table(W3, file="../Result/W.txt", row.names=FALSE, col.names=FALSE)
print(c('W', dim(W)))
print('Fused Graph')


### Given number of clusters
C = 3			# number of clusters
group = spectralClustering(W, C); 	# the final subtypes information

write.table(group, file="../Result/group.txt")
# print('Subtyped')
# print(group)
## you can evaluate the goodness of the obtained clustering results by calculate Normalized mutual information (NMI): if NMI is close to 1, it indicates that the obtained clustering is very close to the "true" cluster information; if NMI is close to 0, it indicates the obtained clustering is not similar to the "true" cluster information.

# W_show = (W*1000)^6
# print(W)
displayClusters(W, group);

SNFNMI = calNMI(group, truelabel)
print(c('Mutual Info', SNFNMI))
# print(truelabel)
# print(group)

# ## you can also find the concordance between each individual network and the fused network

# ConcordanceMatrix = concordanceNetworkNMI(list(W, W1,W2));



# # ### finding cluster numbers by itself
# # ## Here we provide two ways to estimate the number of clusters. Note that,
# # ## these two methods cannot guarantee the accuracy of esstimated number of
# # ## clusters, but just to offer two insights about the datasets.
# # estimationResult = estimateNumberOfClustersGivenGraph(W, 2:5)