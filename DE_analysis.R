library(limma)
library(Biobase)
library(beadarray)

setwd("D:/20Fall/02-718/project/")
exp_df <- read.delim("multi-omics_data/HT_HG-U133A_filtered.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cluster <- read.delim("SNF_subtype.txt", header = FALSE, row.names  = 1)
summary(exp_df[,c(1:5)])

pheno <- new("AnnotatedDataFrame", data = cluster)

eset_non_norm <- new("ExpressionSet", exprs = as.matrix(exp_df), phenoData=pheno)

boxplot(exprs(eset_non_norm), main = "Before normalization", las=2)
# Log transformation
# exprs(eset_non_norm) <- log2(exprs(eset_non_norm))
#head(exprs(eset_non_norm))
# boxplot(exprs(eset_non_norm), main = "After log transformation",las=2)
eset_norm <- normaliseIllumina(eset_non_norm, method="quantile", transform="none")
# Quantile normalization
boxplot(exprs(eset_norm), main = "After normalization", las=2)

design <- model.matrix(~0 + as.factor(cluster$V2))
colnames(design) <- levels(as.factor(cluster$V2))

# Fit linear model
fit<-limma::lmFit(exprs(eset_norm),design)

contrasts_group <- c("G1-G2","G2-G3","G3-G1")
cont.matrix <- makeContrasts(contrasts =  contrasts_group, levels = design)


##########Extract linear model fit for the contrasts###################                                                  
fit_1  <- contrasts.fit(fit, cont.matrix)
fit_1  <- eBayes(fit_1)


# Write DE results for each comparison group
for (i in 1:length(contrasts_group)) {
  write.table(topTable(fit_1,n = Inf, coef = i, sort.by = 'logFC', p.value = 0.05), paste0("DE_contrast_",i,".txt"),
              quote = FALSE, sep = "\t")
}
