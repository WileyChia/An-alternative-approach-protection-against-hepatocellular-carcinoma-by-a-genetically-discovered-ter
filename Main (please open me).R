# ================================================================
# Source file description
# This source file is the main source file for the project
# ================================================================

# Checking Dependencies
if (!require('BiocManager')) install.packages('BiocManager')
BiocManager::install(c("ggplot2", "survival", "randomForestSRC", "dplyr", 
                       "survminer", "edgeR", "statmod", "EnhancedVolcano", "clusterProfiler", 
                       "ReactomePA", "pathview", "enrichplot", "randomForest", "Boruta", 
                       "caret", "pROC", "ggpubr", "ggthemes", "pheatmap", "org.Hs.eg.db", 
                       "DOSE", "impute", "limma", "Hmisc", "forestplot", "glmnet", "foreign", 
                       "ggrisk", "rms", "survivalROC", "ggplotify", "grid", "devtools"), update = F)
devtools::install_github("Tong-Chen/ImageGP", quiet = T)

library(ImageGP)
library(ggplot2)
library(survival)
library(randomForestSRC)

source("Data.R")
source("Function.R")

# This data file records expression data for all TCGA patients
load("Data/XENA.TCGA.ASY.SYMBOL.Counts.RData")

# Plot OS, RFS, DFS for TLS+ and TLS- groups
p = data.frame(row.names = TLS.Sample)
p = FillSurv(p)
p$Group = AssignGroup(p)

OS(p, "OS for TLS+ & TLS-", 1)
OS(p, "OS for TLS+ & TLS-", 2)
OS(p, "OS for TLS+ & TLS-", 5)

RFS(p, "RFS for TLS+ & TLS-", 1)
RFS(p, "RFS for TLS+ & TLS-", 2)
RFS(p, "RFS for TLS+ & TLS-", 5)

DFS(p, "DFS for TLS+ & TLS-", 1)
DFS(p, "DFS for TLS+ & TLS-", 2)
DFS(p, "DFS for TLS+ & TLS-", 5)

p$Time = p$Status = NULL

MixCox(p, timeSpan = 1)
MixCox(p, timeSpan = 2)
MixCox(p, timeSpan = 5)
MixCox(p, timeSpan = 1, use = "RFS")
MixCox(p, timeSpan = 2, use = "RFS")
MixCox(p, timeSpan = 5, use = "RFS")
MixCox(p, timeSpan = 1, use = "DFS")
MixCox(p, timeSpan = 2, use = "DFS")
MixCox(p, timeSpan = 5, use = "DFS")

# Perform variance analysis
p = data.frame(row.names = TLS.Sample)
p$Group = AssignGroup(p)
gene.diff = DEG(data.frame(row.names = rownames(p),Group = p$Group), title = "TLS+ and TLS-")
gene.diff = subset(gene.diff, gene.diff$FDR < .05)

# Prepare data
gene.diff$RV = abs(gene.diff$logFC)
gene.diff = gene.diff[order(gene.diff$RV, decreasing = T), ]
gene.diff = subset(gene.diff, gene.diff$FDR < .05)

# KEGG enrichment
GOKEGG(gene.diff)

#### random forest, this step is VERY time consuming
expression.forest = SubSetRowAndCol(DATA.ASY.SYMBOL, row = rownames(gene.diff), col = TLS.Sample)
expression.forest = data.frame(t(expression.forest))
result.forest = RandomForest(expression.forest, AssignGroup(expression.forest))

sp_boxplot(result.forest$boxplot, melted=T, xvariable = "Variable", yvariable = "Importance",
           legend_variable = "finalDecision", legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Confirmed"),
           xtics_angle = 90)
caret::featurePlot(result.forest$featurePlot, result.forest$featureGroup, plot="box")
plot(result.forest$mtryplot)
dotPlot(result.forest$dotplot)
plot(result.forest$roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres = TRUE, main = "ROC of RF result")
gene.top20 = result.forest$model$importance
gene.top20 = gene.top20[order(gene.top20, decreasing = T),]
gene.top20 = names(gene.top20[1:20]) 

# Further multivariate logistic regression using the top 20 genes of random forest importance
expression.sig = SubSetRowAndCol(DATA.ASY.SYMBOL,  row = gene.top20, col = TLS.Sample)
expression.sig = data.frame(t(expression.sig))
name = rownames(expression.sig)
expression.sig = lapply(expression.sig, function(x) {
  med = median(x)
  x = as.factor(ifelse(x > med, "High", "Low"))
})
expression.sig = data.frame(expression.sig)
rownames(expression.sig) = name
expression.sig$Group = AssignGroup(expression.sig)
ROC(expression.sig, gene.top20, "ROC curve of forecast TLS" )

# Significant gene correlation plots
expression.sig = SubSetRowAndCol(DATA.ASY.SYMBOL, row = gene.top20, col = TLS.Sample)
expression.sig = data.frame(t(expression.sig))
expression.sig$Group = AssignGroup(expression.sig)
Violin(expression.sig)

# Draw heatmap
expression.sig = SubSetRowAndCol(DATA.ASY.SYMBOL, row = gene.top20)
expression.sig = data.frame(t(expression.sig))
expression.sig = subset(expression.sig, substr(rownames(expression.sig),14,14) == "0")
expression.sig$TLS = AssignGroup(expression.sig)
expression.sig$TLS = factor(
  expression.sig$TLS,
  levels = c("TLS-", "TLS+"))
expression.sig = subset(expression.sig, !is.na(expression.sig$TLS))
expression.sig = expression.sig[order(expression.sig$TLS),]
annotation_col = data.frame(row.names = rownames(expression.sig))
annotation_col$TLS = AssignGroup(annotation_col)
annotation_col$TLS = factor(
  annotation_col$TLS,
  levels = c("TLS-", "TLS+"))
colors = list(TLS = c("TLS-" = "blue", "TLS+" = "yellow"))
expression.sig$TLS = NULL
expression.sig = log2(expression.sig + 1)

AnnotationHeatmap(
  expression.sig,
  annotation_col,
  colors, 
  c(table(annotation_col$TLS)[1]))

# GO enrichment analysis
options(digits = 5)
GOBP(gene.top20, type = "symbol")
GOCC(gene.top20, type = "symbol")
GOMF(gene.top20, type = "symbol")
GODO(gene.top20, type = "symbol")
options(digits = 22)

# Immune cell infiltration
result.c = read.table("Data/CIBERSORT-Results.txt", sep = "\t", header = T, row.names = "Mixture")
result.c$P.value = result.c$Correlation = result.c$RMSE = NULL
expression.topgene = SubSetRowAndCol(DATA.ASY.SYMBOL, row = gene.top20, col = rownames(result.c))
expression.topgene = data.frame(t(expression.topgene))
CorGroup(result.c, expression.topgene, method = "spearman", cluster_cols = T, title = "cor between immscore and topgene")

# Feature genes correlated with TLS feature genes
expression.tls9 = SubSetRowAndCol(DATA.ASY.SYMBOL, row = GENE.TLS9.SYMBOL, col = rownames(result.c))
expression.tls9 = data.frame(t(expression.tls9))
CorGroup(expression.topgene, expression.tls9, method = "spearman", cluster_cols = T, title = "cor between immscore and TLS9")
expression.tls12 = SubSetRowAndCol(DATA.ASY.SYMBOL, row = GENE.TLS.SYMBOL, col = rownames(result.c))
expression.tls12 = data.frame(t(expression.tls12))
CorGroup(expression.topgene, expression.tls12, method = "spearman", cluster_cols = T, title = "cor between immscore and TLS12")
expression.tls40 = SubSetRowAndCol(DATA.ASY.SYMBOL, row = GENE.TLS40.SYMBOL, col = rownames(result.c))
expression.tls40 = data.frame(t(expression.tls40))
CorGroup(expression.topgene, expression.tls40, method = "spearman", cluster_cols = T, title = "cor between immscore and TLS40")
expression.tls50 = SubSetRowAndCol(DATA.ASY.SYMBOL, row = GENE.TLS50.SYMBOL, col = rownames(result.c))
expression.tls50 = data.frame(t(expression.tls50))
CorGroup(expression.topgene, expression.tls50, method = "spearman", cluster_cols = T, title = "cor between immscore and TLS50")

# Drug sensitivity
DrugSensitivity(gene.top20)

# Prognostic model
result.lasso = LASSO(gene.top20, type = "symbol")
rs = result.lasso$RiskModel
rs = FillSurv(rs)
rs$Group = ifelse(rs$IsLowRisk, "LowRisk", "HighRisk")
OS(rs,"20-gene lasso")
RFS(rs,"20-gene lasso")
DFS(rs,"20-gene lasso")
colnames(rs)[2] = "Factor"
TimeROC(rs, "20-gene lasso", method = "KM")

rs$Time = rs$Status = rs$IsLowRisk = rs$Group = NULL
colnames(rs)[1] = "Risk Score"
MixCox(rs)

# External validation of ICGC data
load("Data/ICGC.JP.ASY.RData")
load("Data/ICGC.JP.COL.RData")

# SMIM3 used to be named C5orf62, this data table in ICGC uses the former name, correct it
colnames(ICGC.JP.ASY)[ colnames(ICGC.JP.ASY) == "C5orf62" ] = "SMIM3"
x = SubSetRowAndCol(ICGC.JP.ASY, col = result.lasso$FiltedGenes, row = rownames(ICGC.JP.COL[ICGC.JP.COL$Type == "T",]))
x = log2(x + 1)
x = na.omit(x)
result.icgc = data.frame(row.names = rownames(x))
result.icgc$RiskScore = eval(parse(text = result.lasso$Formula))
cutoff = eval(parse(text = result.lasso$Cutoff))
result.icgc$IsLowRisk = result.icgc$RiskScore < cutoff
# Merge survival information
result.icgc = merge(result.icgc, ICGC.JP.COL, by = "row.names")
rownames(result.icgc) = result.icgc$Row.names
result.icgc$Row.names = result.icgc$Type = NULL
result.icgc$Group = ifelse(result.icgc$IsLowRisk, "LowRisk", "HighRisk")
result.icgc$Time = as.numeric(result.icgc$Time)
result.icgc$Status = as.numeric(result.icgc$Status)
result.icgc = na.omit(result.icgc)
OS(result.icgc, "OS in ICGC cohort")
result.icgc$Factor = result.icgc$RiskScore
TimeROC(result.icgc, "ROC in ICGC cohort", method = "KM")

result.icgc$IsLowRisk = result.icgc$Factor = result.icgc$Group = NULL
result.icgc$Gender = ifelse(result.icgc$Gender == "male", "Male", "Female")
result.icgc$Stage = ifelse(result.icgc$Stage == 1 , "I", "II-IV")
result.icgc = data.frame(
  Age = result.icgc$Age, 
  Gender = result.icgc$Gender,
  Stage = result.icgc$Stage, 
  RiskScore = result.icgc$RiskScore,
  Time = result.icgc$Time,
  Status = result.icgc$Status)
MixCox(result.icgc, needSurv = F)

# Risk model heatmap
expression.lasso = SubSetRowAndCol(DATA.ASY.SYMBOL, row = result.lasso$FiltedGenes, col = rownames(result.lasso$RiskModel))
expression.lasso = data.frame(t(expression.lasso))
expression.lasso = merge(expression.lasso, result.lasso$RiskModel, by = "row.names")
expression.lasso = expression.lasso[order(expression.lasso$RiskScore),]
rownames(expression.lasso) = expression.lasso$Row.names
expression.lasso$Row.names = expression.lasso$RiskScore = expression.lasso$IsLowRisk = NULL
annotation_col = data.frame(row.names = rownames(expression.lasso))
annotation_col = FillSurvTable(annotation_col)
annotation_col$NewTumorTime = as.numeric(annotation_col$NewTumorTime) / 365
annotation_col$Time = as.numeric(annotation_col$Time) / 365
annotation_col$NewTumorEvent = factor(ifelse(annotation_col$NewTumorEvent == "NONE", "", annotation_col$NewTumorEvent), levels = c("YES", "NO"))
annotation_col$TLS = AssignGroup(annotation_col, til = "TIL")
annotation_col$TLS = ifelse(annotation_col$TLS == "TLS+", 3, ifelse(annotation_col$TLS == "TIL", 2, ifelse(annotation_col$TLS == "TLS-", 1, 0)))
annotation_col$Status = ifelse(annotation_col$Status == 0, "Alive", "Deceased")
annotation_col$T = ifelse(annotation_col$T == "T1", 1,
                          ifelse(annotation_col$T == "T2", 2,
                                 ifelse(annotation_col$T == "T3", 3, 4)))
annotation_col$N = ifelse(annotation_col$N == "NX", 0,
                          ifelse(annotation_col$N == "N0", 1, 2))
annotation_col$M = ifelse(annotation_col$M == "MX", 0,
                          ifelse(annotation_col$M == "M0", 1, 2))
annotation_col$Grade = ifelse(annotation_col$Grade == "G1", 1,
                          ifelse(annotation_col$Grade == "G2", 2,
                                 ifelse(annotation_col$Grade == "G3", 3, 4)))
annotation_col$Stage = ifelse(annotation_col$Stage == "I", 1,
                              ifelse(annotation_col$Stage == "II", 2,
                                     ifelse(annotation_col$Stage == "III", 3, 4)))
annotation_col = merge(annotation_col, result.lasso$RiskModel, by = "row.names")
rownames(annotation_col) = annotation_col$Row.names
annotation_col$IsLowRisk = annotation_col$Row.names = NULL
expression.lasso = SubSetRowAndCol(expression.lasso, row = rownames(annotation_col))
AnnotationHeatmap(expression.lasso, annotation_col)