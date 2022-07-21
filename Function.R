# ================================================================
# Source file description
# This source file contains all the predefined functions used in Main.
# ================================================================

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Group the TCGA samples with reference to the TLS data that has been observed under the microscope, based on the row names of the data parameter
# ----------------------------------------------------------------
# Parameter description.
# ----------------------------------------------------------------
# data: the sample to be grouped
# isC: whether the sample is of type c()
# tls: when the sample belongs to the TLS group, the content returned by the function
# notls: when the sample does not belong to the TLS, the function returns the content
# notls: when the sample does not belong to the TLS, the function returns the content
# til: when the sample belongs to the TIL group, the function returns the content

# ================================================================
AssignGroup = function(data, isC = F, tls = "TLS+", notls = "TLS-", til = "TLS+"){
  if (isC) {
    return(ifelse(data %in% TLS.None, notls, ifelse(data %in% TLS.TIL, til, ifelse(data %in% TLS.TLS, tls, ""))))
  }
  return(ifelse(rownames(data) %in% TLS.None, notls, ifelse(rownames(data) %in% TLS.TIL, til,  ifelse(rownames(data) %in% TLS.TLS, tls, ""))))
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Plot the patient's survival curve
# ----------------------------------------------------------------
# Parameter description.
# ----------------------------------------------------------------
# data: the sample data waiting to be plotted for OS, the data needs to contain Group, Time and Status
# title: title of the image
# time: the time range to be calculated, in years
# ================================================================
OS <- function(data, title, time = 10) {
  library(dplyr)
  library(survival)
  library(survminer)
  data <- data.frame(data)
  if (max(data$Time) > 1000) {
    data$Time = data$Time / 365
  }
  if (time != 10) {
    data$Status = ifelse(data$Time > time, ifelse(data$Status == 1, 0 ,data$Status), data$Status)
    data$Time = ifelse(data$Time > time, time , data$Time)
  }
  yt <- Surv(data$Time, data$Status)
  fitt <- surv_fit(yt ~ Group, data = data )
  ptpv <- as.character(survminer::surv_pvalue(fitt)$pval)
  Fit <<- fitt
  Yt <<- yt
  lab = names(table(data$Group))
  pt <- ggsurvplot(fitt,
                   data = data,
                   conf.int = TRUE,
                   pval = TRUE,
                   pval.method = TRUE,
                   legend.title = "Group",
                   legend.labs = c(lab[1], lab[2]),
                   xlab = "Time (year)",
                   surv.median.line = "hv",
                   risk.table = TRUE,
                   title = title) 
  print(paste("O S p:",ptpv))
  rm(data)
  rm(ptpv)
  pt
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Plot the patient's disease-free survival curve
# ----------------------------------------------------------------
# Parameter description.
# ----------------------------------------------------------------
# data: sample data waiting to be plotted for DFS, data needs to contain Group, will read other data from RData file
# title: title of the image
# time: the time range to be calculated, in years
# ================================================================
DFS <- function(data, title, time = 10) {
  library(dplyr)
  library(survival)
  library(survminer)
  data <- data.frame(data)
  data$Time = data$Status = NULL
  surv <- GetSurvTable()
  data$PID = substr(rownames(data),1,12)
  surv$PID = substr(rownames(surv),1,12)
  surv = subset(surv, surv$NewTumorEvent != "NONE")
  data <- merge(data, surv, by = "PID")
  data$PID = NULL
  data$NewTumorTime = as.numeric(data$NewTumorTime)
  data$Time = ifelse(data$NewTumorEvent == "YES", data$NewTumorTime, data$Time)
  data$Time = data$Time / 365
  data = data[nchar(data$NewTumorEvent) != 4,]
  data$Status <- ifelse(data$Status == 1, 1, ifelse(data$NewTumorEvent == "YES", 1, 0))
  data$Status = ifelse(data$Time > time, ifelse(data$Status == 1, 0 ,data$Status), data$Status)
  data$Time = ifelse(data$Time > time, time , data$Time)
  rm(surv)
  yt <- Surv(data$Time, data$Status)
  fitt <- surv_fit(yt ~ Group, data = data )
  ptpv <- as.character(survminer::surv_pvalue(fitt)$pval)
  Fit <<- fitt
  Yt <<- yt
  lab = names(table(data$Group))
  pt <- ggsurvplot(fitt,
                   data = data,
                   conf.int = TRUE,
                   pval = TRUE,
                   pval.method = TRUE,
                   legend.title = "Group",
                   legend.labs = c(lab[1], lab[2]),
                   xlab = "Time (year)",
                   surv.median.line = "hv",
                   risk.table = TRUE,
                   title = title)
  print(paste("DFS p:",ptpv))
  rm(data)
  rm(ptpv)
  pt
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Plot the patient's early relapse-free survival curve
# ----------------------------------------------------------------
# Parameter description.
# ----------------------------------------------------------------
# data: sample data waiting to be plotted for RFS, data needs to contain Group, will read other data from RData file
# title: title of the image
# time: the time range to be calculated, in years
# ================================================================
RFS <- function(data, title, time = 10) {
  library(dplyr)
  library(survival)
  library(survminer)
  data <- data.frame(data)
  data$Time = data$Status = NULL
  surv <- GetSurvTable()
  data$PID = substr(rownames(data),1,12)
  surv$PID = substr(rownames(surv),1,12)
  data <- merge(data, surv, by = "PID")
  data$PID = NULL
  data$NewTumorTime = as.numeric(data$NewTumorTime)
  data$Time = data$NewTumorTime
  data$Time = data$Time / 365
  data = data[nchar(data$NewTumorEvent) != 4,]
  data$Status <- ifelse(data$NewTumorEvent == "YES", 1, 0)
  data$Status = ifelse(data$Time > time, ifelse(data$Status == 1, 0 ,data$Status), data$Status)
  data$Time = ifelse(data$Time > time, time , data$Time)
  rm(surv)
  yt <- Surv(data$Time, data$Status)
  fitt <- surv_fit(yt ~ Group, data = data )
  ptpv <- as.character(survminer::surv_pvalue(fitt)$pval)
  Fit <<- fitt
  Yt <<- yt
  lab = names(table(data$Group))
  pt <- ggsurvplot(fitt,
                   data = data,
                   conf.int = TRUE,
                   pval = TRUE,
                   pval.method = TRUE,
                   legend.title = "Group",
                   legend.labs = c(lab[1], lab[2]),
                   xlab = "Time (year)",
                   surv.median.line = "hv",
                   risk.table = TRUE,
                   break.x.by = .5,
                   xlim = c(0,2),
                   title = title)
  print(paste("RFS p:",ptpv))
  rm(data)
  rm(ptpv)
  pt
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Get some of the patient's survival information data
# ================================================================
GetSurv = function(){
  load(file = "Data/XENA.TCGA.COL.CLEAN.RData")
  surv <- data.frame(row.names = rownames(DATA.COL.CLEAN),
                     Status = c(as.numeric(sub("LIVING", "0", sub("DECEASED", "1", DATA.COL.CLEAN$Status)))),
                     Time = c(as.numeric(DATA.COL.CLEAN$Time)))
  surv = surv[complete.cases(surv),]
  surv[surv==''] <- NA
  return(unique(na.omit(surv))) 
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Will add partial survival data to the parameter data and return the added data
# ================================================================
FillSurv = function(data){
  data = data.frame(data)
  surv = GetSurv()
  data$PID = substr(rownames(data),1,12)
  data$ID = rownames(data)
  surv$PID = substr(rownames(surv),1,12)
  data = merge(data, surv, by = "PID")
  rm(surv)
  rownames(data) = data$ID
  data$ID = data$PID = NULL
  return(data)
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Get the patient's full survival information data
# ================================================================
GetSurvTable = function(){
  load(file = "Data/XENA.TCGA.COL.CLEAN.RData")
  surv <- data.frame(row.names = rownames(DATA.COL.CLEAN),
                     Age = c(as.numeric(DATA.COL.CLEAN$Age)),
                     Gender = c(DATA.COL.CLEAN$Gender),
                     Stage =c(DATA.COL.CLEAN$Stage) ,
                     T = c(DATA.COL.CLEAN$T),
                     N = c(DATA.COL.CLEAN$N),
                     M = c(DATA.COL.CLEAN$M),
                     Grade = c(DATA.COL.CLEAN$Grade),
                     NewTumorTime = c(ifelse(is.na(DATA.COL.CLEAN$NewTumorTime),"NA",DATA.COL.CLEAN$NewTumorTime)),
                     NewTumorEvent = c(ifelse(DATA.COL.CLEAN$NewTumorEvent != "YES" & DATA.COL.CLEAN$NewTumorEvent != "NO","NONE",DATA.COL.CLEAN$NewTumorEvent) ),
                     Status = c(as.numeric(sub("LIVING", "0", sub("DECEASED", "1", DATA.COL.CLEAN$Status)))),
                     Time = c(as.numeric(DATA.COL.CLEAN$Time)))
  surv = surv[complete.cases(surv),]
  surv[surv==''] <- NA
  return(unique(na.omit(surv))) 
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Will add all the survival data to the parameter data and return the added data
# ================================================================
FillSurvTable = function(data){
  data = data.frame(data)
  surv = GetSurvTable()
  data$PID = substr(rownames(data),1,12)
  data$ID = rownames(data)
  surv$PID = substr(rownames(surv),1,12)
  data = merge(data, surv, by = "PID")
  rm(surv)
  rownames(data) = data$ID
  data$ID = data$PID = NULL
  return(data)
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Gene Difference Analysis
# ----------------------------------------------------------------
# Parameter description.
# ----------------------------------------------------------------
#group: dataframe, rowname is the sample name, which needs to contain the Group column
#title: image title
# ================================================================
DEG = function(group, title = "DEG"){
  library(edgeR)
  load("Data/XENA.TCGA.ASY.SYMBOL.Counts.RData")
  countsData <- DATA.ASY.SYMBOL
  countsData <- countsData[which(rowSums(countsData == 0) == 0), ]
  countsData <- data.frame(t(countsData))
  countsData <- subset(countsData, rownames(countsData) %in% rownames(group))
  countsData <- data.frame(countsData)
  countsData <- merge(countsData, group, by = "row.names")
  rownames(countsData) = countsData$Row.names
  countsData$Row.names = NULL
  group <- as.factor(countsData$Group)
  countsData$Group <- NULL
  countsData <- t(countsData)
  countsData = data.frame(countsData)
  countsData = as.matrix(2^countsData - 1)
  exp = apply(countsData, 2, as.integer)
  rownames(exp) = rownames(countsData)
  colnames(exp) = colnames(countsData)
  countsData = exp
  library(statmod)
  dgelist <- DGEList(counts = countsData, group = group)
  keep <- rowSums(cpm(dgelist) > 1) >= 2
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
  dgelist_norm <- calcNormFactors(dgelist, method = "TMM")
  design <- model.matrix(~group)
  dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
  fit <- glmFit(dge, design, robust = TRUE)
  lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
  gene_diff <- lrt
  gene_diff = data.frame(gene_diff)
  gene_diff <- gene_diff[order(gene_diff$FDR, gene_diff$logFC, decreasing = c(FALSE, TRUE)), ]
  library(EnhancedVolcano)
  gene_diff <- data.frame(gene_diff)
  p = EnhancedVolcano(
    gene_diff,
    lab = rownames(gene_diff),
    x = "logFC",
    y = "PValue",
    title = "Differential expression analysis",
    subtitle = title,
    pCutoff = .05,
    FCcutoff = 0,
    col = c("black", "blue", "green", "red"),
    colAlpha = 1,
    legendPosition = "right",
    legendLabSize = 14,
    legendIconSize = 5.0)
  plot.volcano <<- p
  return(gene_diff)
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# KEGG enrichment analysis
# ----------------------------------------------------------------
# Parameter description.
# ----------------------------------------------------------------
# deg: results of the difference analysis
# genes: set of genes to be explored
# ================================================================
GOKEGG = function(deg, genes = NULL){
  library(clusterProfiler)
  library(ReactomePA)
  if (!is.null(genes)) {
    deg = subset(deg, rownames(deg) %in% genes)
  }
  deg = deg[order(deg$logFC, decreasing = T),]
  logFC = deg$logFC
  names(logFC) = rownames(deg)
  gene <- names(logFC)
  gene = bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  rownames(gene) = gene$SYMBOL
  names(logFC) = gene[names(logFC),]$ENTREZID
  k = gseKEGG(
    logFC,
    keyType  = 'ncbi-geneid',
    organism = 'hsa',
    pvalueCutoff = .001,
    pAdjustMethod = "none")
  library(pathview)
  library(enrichplot)
  plot(gseaplot2(k,
                 1:10,
                 pvalue_table = T)
  )
  kk = gsePathway(logFC)
  plot(gseaplot2(kk,
                 1:10,
                 pvalue_table = T))
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Crop the rows and columns of a dataframe
# ================================================================
SubSetRowAndCol = function(data, row = NULL, col = NULL){
  if (!is.null(row)) {
    data = subset(data, rownames(data) %in% row)
  }
  if (!is.null(col)) {
    data = t(data)
    data = t(subset(data, rownames(data) %in% col))
  }
  return(data.frame(data))
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Random forest algorithm
# ----------------------------------------------------------------
# Parameter description.
# ----------------------------------------------------------------
# data: expression matrix, listed as gene, behavior sample
# group: numeric - do regression, factorial - do classification
# ================================================================
RandomForest = function(data, group){
  library(ggplot2)
  library(randomForest)
  library(Boruta)
  library(ImageGP)
  data = data.frame(data)
  # Splitting the training and validation sets
  data$Group = as.factor(group)
  set.seed(666666)
  data.Train = sample(nrow(data), 0.7 * nrow(data))
  data.Test =  data[-data.Train,]
  data.Train=  data[data.Train,]
  trang = data.Train$Group
  testg = data.Test$Group
  data.Train$Group = data.Test$Group = NULL
  data.Train = scale(as.matrix(data.Train))
  data.Test = scale(as.matrix(data.Test))
  data.Train =  data.frame(data.Train)
  data.Test = data.frame(data.Test)
  data.Train$Group = trang
  data.Test$Group = testg
  set.seed(666666)
  boruta <- Boruta(Group ~ ., data = data.Train, pValue = 0.05, mcAdj = T, maxRuns = 300)
  # Print important genes
  table(boruta$finalDecision)
  library(dplyr)
  boruta.imp <- function(x){
    imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
    colnames(imp) <- c("Variable","Importance")
    imp <- imp[is.finite(imp$Importance),]
    variableGrp <- data.frame(Variable=names(x$finalDecision), 
                              finalDecision=x$finalDecision)
    showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                          finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
    variableGrp <- rbind(variableGrp, showGrp)
    boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
    sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>% 
      summarise(median=median(Importance)) %>% arrange(median)
    sortedVariable <- as.vector(sortedVariable$Variable)
    boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
    invisible(boruta.variable.imp)
  }
  boruta.variable.imp <- boruta.imp(boruta)
  # Draw the important genes
  boxplot = boruta.variable.imp
  boruta.finalVarsWithTentative <- 
    data.frame(Item=getSelectedAttributes(boruta, withTentative = T), Type="Boruta_with_tentative")
  featurePlot = data.Train[,boruta.finalVarsWithTentative$Item]
  featureGroup = data.Train$Group
  generateTestVariableSet <- function(num_toal_variable){
    max_power <- ceiling(log10(num_toal_variable))
    tmp_subset <- c(unlist(sapply(1:max_power, function(x) (1:10)^x, simplify = F)), ceiling(max_power/3))
    base::unique(sort(tmp_subset[tmp_subset<num_toal_variable]))
  }
  boruta_train_data <- data.Train[, boruta.finalVarsWithTentative$Item]
  boruta_mtry <- generateTestVariableSet(ncol(boruta_train_data))
  library(caret)
  # Create model with default parameters
  set.seed(666666)
  trControl <- trainControl(method="repeatedcv", 
                            savePredictions = "all",
                            number = 10, 
                            verboseIter = TRUE,
                            repeats = 5)
  set.seed(666666)
  tuneGrid <- expand.grid(mtry = boruta_mtry)
  borutaConfirmed_rf_default <- train(Group ~ ., data = data.Train, method="rf", 
                                      tuneGrid = tuneGrid, # 
                                      metric="Accuracy",
                                      num.threads = 8,
                                      trControl = trControl)
  print(borutaConfirmed_rf_default)
  plotmtry = borutaConfirmed_rf_default
  plotdot = varImp(borutaConfirmed_rf_default)
  borutaConfirmed_rf_default_finalmodel <- borutaConfirmed_rf_default$finalModel
  # Evaluate the classification effect of the model on the training set
  g = data.Train$Group
  dt = data.Train
  dt$Group = NULL
  predictions_train <- predict(borutaConfirmed_rf_default_finalmodel, newdata = dt)
  confusionMatrix(predictions_train, g)
  # Blind evaluation to assess the effectiveness of the model when applied to the test set
  g = data.Test$Group
  dt = data.Test
  dt$Group = NULL
  prediction_prob <- predict(borutaConfirmed_rf_default_finalmodel, newdata = dt, type = "prob")
  library(pROC)
  roc_curve <- roc(g, prediction_prob[,1])
  print(roc_curve)
  return(list(model = borutaConfirmed_rf_default_finalmodel, 
              boxplot = boxplot, 
              featurePlot = featurePlot, 
              featureGroup = featureGroup,
              mtryplot = plotmtry,
              dotplot = plotdot,
              roc = roc_curve))
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# For plotting ROC curves
# ================================================================
ROC = function(data, factors, title){
  library(pROC)
  group = data$Group
  data$Group = NULL
  exp = parse(text = paste0("lm(group == 'TLS-'~", paste0(factors,collapse = "+") ,", data=data)"))
  model1 <- eval(exp)
  pre <- predict(model1, type='response')
  modelroc <- roc(group, pre)
  p = plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
           grid.col=c("green", "red"), max.auc.polygon=TRUE,
           auc.polygon.col="skyblue", print.thres = TRUE, main = title)
  return(p)
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# For plotting violin, sample rows in data, columns for genes, data need to contain Group
# ================================================================
Violin = function(data){
  library(ggpubr)
  data = data.frame(data)
  genes = colnames(data)
  data = cbind(data,Level = "")
  data = cbind(data,GENE = "")
  dxDATA = data
  data = data[0,]
  genes = subset(genes , genes != "Group")
  for (gene in genes) {
    dxDATA$GENE = gene
    exp = parse(text = paste0("dxDATA$",gene))
    dxDATA$Level = eval(exp)
    data = rbind( data , dxDATA) 
  }
  library(ggthemes)
  p <- ggviolin(data, x="Group", y="Level", color = "Group", fill = "Group", palette = "jco", facet.by = "GENE")
  p = p+stat_compare_means(aes(label=..p.signif..), label.x = 1.5, method = "wilcox.test")
  p=p+geom_boxplot(aes(x = Group,y = Level),width = .1,fill="black",color="gray")
  p=p+stat_summary(fun="mean",geom="point",shape=23,size=1.5,fill="white")
  p=p+scale_color_tableau()
  p
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# For drawing a heat map with comments
# ================================================================
AnnotationHeatmap <- function(data, annotation_col, annotation_colors = NULL, gaps_col = NULL, title = "Heatmap") {
  library(pheatmap)
  library(ggplot2)
  p <- pheatmap(
    t(data),
    scale = "row",
    trace = "none",
    gaps_col = gaps_col,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    color = colorRampPalette(colors = c("blue","white","red"))(100),
    border = F,
    cluster_row = T,
    cluster_col = F,
    show_colnames = F,
    show_rownames = T,
    silent = T,
    main = title
  )
  require(ggplotify)
  p = as.ggplot(p)
  return(p)
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# GO Biological Process
# ================================================================
GOBP <- function(genes, type = "ensg") {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  type = ifelse(type == "ensg", "ENSEMBL", "SYMBOL")
  # go生物过程
  genes <- substr(genes, 0, 15)
  goBP <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = type,
    ont = "BP",
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5
  )
  return(enrichplot::dotplot(goBP))
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# GO cell components
# ================================================================
GOCC <- function(genes, type = "ensg") {
  library(clusterProfiler)
  genes <- substr(genes, 0, 15)
  type = ifelse(type == "ensg", "ENSEMBL", "SYMBOL")
  # go细胞组分  基因产物位于哪个细胞器起作用
  goCC <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = type,
    ont = "CC",
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5
  )
  return(enrichplot::dotplot(goCC))  # barplot
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# GO molecule function
# ================================================================
GOMF <- function(genes, type = "ensg") {
  library(clusterProfiler)
  type = ifelse(type == "ensg", "ENSEMBL", "SYMBOL")
  # go分子功能
  genes <- substr(genes, 0, 15)
  goMF <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = type,
    ont = "MF",
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5
  )
  enrichplot::dotplot(goMF)
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Disease ontology
# ================================================================
GODO <- function(genes, type = "ensg") {
  library(DOSE)
  library(org.Hs.eg.db)
  type = ifelse(type == "ensg", "ENSEMBL", "SYMBOL")
  genes <- substr(genes, 0, 15)
  geneid <- mapIds(
    x = org.Hs.eg.db,
    keys = genes,
    keytype = type,
    column = "ENTREZID"
  )
  # GO疾病
  goDO <- enrichDO(
    gene = geneid,
    ont = "DO",
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5
  )
  enrichplot::dotplot(goDO)
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Drug sensitivity
# ================================================================
DrugSensitivity = function(genes, orderby = "p"){
  library(impute)
  library(limma)
  load(file = "Data/CellMiner.DTP.NCI60.RData")
  load(file = "Data/CellMiner.DTP.RNASeq.RData")
  DTP.NCI60 = as.matrix(DTP.NCI60)
  exp = DTP.RNASeq[genes,]
  outData = data.frame()
  for (gene in row.names(exp)) {
    x = as.numeric(exp[gene,])
    for (drug in row.names(DTP.NCI60)) {
      y = as.numeric(DTP.NCI60[drug,])
      corT = cor.test(x, y, method = "pearson")
      cor = corT$estimate
      pvalue = corT$p.value
      if (!is.na(pvalue)) {
        if (pvalue < .05) {
          outV = cbind(gene, drug, cor, pvalue)
          outData = rbind(outData, outV)
        }}
    }
  }
  outData$RC = abs(as.numeric(outData$cor))
  outData = outData[order(outData$RC, decreasing = T),]
  outData = outData[order(outData$pvalue),]
  library(ggplot2)
  library(ggpubr)
  plotList = list()
  plotListB = list()
  plotNum = 16
  if (nrow(outData) < plotNum) {
    plotNum = nrow(outData)
  }
  for (i in 1:plotNum) {
    g = outData[i, 1]
    d = outData[i, 2]
    x = as.numeric(exp[g, ])
    y = as.numeric(DTP.NCI60[d, ])
    cor = sprintf("%.03f", as.numeric(outData[i, 3]))
    p = 0
    if (as.numeric(outData[i, 4] < .001)) {
      p = "p < 0.001"
    } else {
      p = paste0("p = ", sprintf("%.03f", as.numeric(outData[i, 4])))
    }
    df = as.data.frame(cbind(x, y))
    pl = ggplot(data = df, aes(x = x, y = y)) +
      geom_point(size = 1) +
      stat_smooth(method = "lm", formula = y ~ x) + 
      labs(x = "Expression", y = "IC50" , title = paste0(d, ", ", g), subtitle = paste0("Cor=", cor, ", ", p)) + 
      theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text = element_blank()) + 
      theme_bw()
    plotList[[i]] = pl
    
    colnames(df)[2] = "IC50"
    df$group = ifelse(df$x > median(df$x), "high", "low")
    
    c = list(c("low", "high"))
    pl = ggboxplot(df,
                   x = "group",
                   y = "IC50",
                   fill = "group",
                   add = "jitter", size = .5,
                   xlab = paste0("the_expression_of_", g),
                   ylab = paste0("IC50_of_",d)) +
      stat_compare_means(comparisons = c,
                         method = "wilcox.test",
                         symnum.args = list(cutpoints = c(0, .001, .01, .05, 1),
                                            symbols = c("***", "**", "*", "ns")))
    plotListB[[i]] = pl
  }
  nr = ceiling(sqrt(plotNum))
  nc = ceiling(plotNum / nr)
  plot(ggarrange(plotlist = plotList, nrow = nr, ncol = nc))
  plot(ggarrange(plotlist = plotListB, nrow = nr, ncol = nc))
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Plot the correlation between two sets of data
# ================================================================
CorGroup = function(data1, data2, method = "pearson",
                    title = "cor",
                    cluster_cols = F, 
                    cluster_rows = T,
                    display_sig = T,
                    show_rownames = T, 
                    show_colnames = T){
  library(Hmisc)
  library(pheatmap)
  data1 = data.frame(data1)
  data2 = data.frame(data2)
  res1 <- rcorr(as.matrix(data1), as.matrix(data2), type = method) 
  corV = res1$r
  corV = corV[,length(data1) + 1 : length(data2)]
  if (is.null(dim(corV))) {
    corV = corV[1:length(data1)]
  } else{
    corV = corV[1:length(data1),]
  }
  pV = res1$P
  pV = pV[,length(data1) + 1 : length(data2)]
  if (is.null(dim(corV))) {
    pV = pV[1:length(data1)]
  } else{
    pV = pV[1:length(data1),]
  }
  if (!display_sig) {
    if (is.null(dim(corV))) {
      corV = as.matrix(corV)
      pV = as.matrix(pV)
      pheatmap(corV, 
               main = title,
               show_rownames = show_rownames,
               show_colnames = show_colnames,
               cluster_cols = cluster_cols,
               cluster_rows = cluster_rows,
      )
    }else{
      pheatmap(corV, 
               main = title,
               show_rownames = show_rownames,
               show_colnames = show_colnames,
               cluster_cols = cluster_cols,
               cluster_rows = cluster_rows,
      )
    }
    return(list(cor = corV, p = pV))
  }
  if (is.null(dim(corV))) {
    corV = as.matrix(corV)
    pV = as.matrix(pV)
    pheatmap(corV, 
             main = title,
             show_rownames = show_rownames,
             show_colnames = show_colnames,
             cluster_cols = cluster_cols,
             cluster_rows = cluster_rows,
             display_numbers = 
               ifelse(corV == 1, 
                      "ns", 
                      ifelse(pV < 0.001, 
                             "***", 
                             ifelse(pV < 0.01, 
                                    "**",
                                    ifelse(pV < 0.05, "*", "")))))
  }else{
    pheatmap(corV, 
             main = title,
             show_rownames = show_rownames,
             show_colnames = show_colnames,
             cluster_cols = cluster_cols,
             cluster_rows = cluster_rows,
             display_numbers = ifelse(corV == 1, "ns", ifelse(pV < 0.001, "***", ifelse(pV < 0.01, "**", ifelse(pV < 0.05, "*", "")))))
  }
  return(list(cor = corV, p = pV))
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Univariate Cox regression
# ================================================================
UnivariateCox <- function(data, pValue = 0.05, plot = F, returnAll = F) {
  library("survival")
  library("survminer")
  options(scipen = 200, digits = 3)
  surv <- GetSurv()
  data = data.frame(data)
  data$PID = substr(rownames(data),1,12)
  surv$PID = substr(rownames(surv),1,12)
  data = merge(data, surv, by = "PID")
  data$PID = NULL
  time <- as.numeric(data$Time)
  status <- as.numeric(data$Status)
  data$Time = data$Status = NULL
  univ_formulas <- sapply(colnames(data),
                          function(x)
                            as.formula(paste(
                              "Surv(time, status) ~ data$", x, ""
                            )))
  univ_models <- lapply(univ_formulas, function(x) {
    coxph(x, data = data)
  })
  univ_results <- lapply(univ_models,
                         function(x) {
                           x <- summary(x)
                           p.value <- signif(x$wald["pvalue"], digits = 2)
                           HR <- signif(x$coef[2], digits = 2)
                           HR.confint.lower <- signif(x$conf.int[, "lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[, "upper .95"], 2)
                           HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
                           res <- c(p.value, HR)
                           names(res) <- c("p.value", "HR (95% CI for HR)")
                           return(res)
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  res <- as.data.frame(res)
  HR <- gsub("[\\(\\)]", "", res$`HR (95% CI for HR)`)
  HR <- gsub("-", " ", HR)
  HR <-
    as.data.frame(do.call(cbind, strsplit(HR, " ")), stringsAsFactors = F)
  names(HR) <- rownames(res)
  HR <- data.frame(t(HR))
  hLow <- as.numeric(HR$X2)
  hHigh <- as.numeric(HR$X3)
  hCenter <- as.numeric(HR$X1)
  if (returnAll) {
    return(res)
  }
  topGenes <- subset(res, p.value <= pValue)
  cat(paste0("\np <", pValue, " factores: (", length(rownames(topGenes)), ") \n"))
  cat(paste(rownames(topGenes), collapse = " \\ "))
  cat("\n")
  res$p.value <-
    ifelse(as.numeric(res$p.value) < 0.001, "< 0.001", res$p.value)
  tabletext <- cbind(
    c("Item", "\n", rownames(res)),
    c("HR (95% CI for HR)", "\n", res$HR),
    c("P Value", "\n", res$p.value)
  )
  library(forestplot)
  pic <- forestplot(
    title = "Hazard Ratio Plot",
    txt_gp = fpTxtGp(
      label = gpar(cex = 1.25),
      ticks = gpar(cex = 1.1),
      xlab = gpar(cex = 1.2),
      title = gpar(cex = 1.2)
    ),
    col = fpColors(
      box = "#1c61b6",
      lines = "#1c61b6",
      zero = "gray50"
    ),
    zero = 1,
    xticks = c(.5, .81, 1, 1.25, 1.5, 1.9, 2),
    lwd.ci = 2,
    colgap = unit(8, "mm"),
    cex = 0.9,
    lineheight = "auto",
    ci.vertices = TRUE,
    ci.vertices.height = 0.4,
    labeltext = tabletext,
    graph.pos = 3,
    mean = c(NA, NA, hCenter),
    lower = c(NA, NA, hLow),
    upper = c(NA, NA, hHigh),
    boxsize = 0.4
  )
  if (plot) {
    plot(pic)
  }
  return(rownames(topGenes))
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# LASSO regression
# ================================================================
LASSO<- function(genes,type= "ensg",dat = NULL){
  library(glmnet)
  library(foreign)
  library(ggrisk)
  library(rms)
  library(survival)
  library(survminer)
  library(survivalROC)
  set.seed(666666)
  if (is.null(dat)) {
    if(type == "ensg"){
      dat <- data.frame(t(DATA.ASY))
    }else{
      dat <- data.frame(t(DATA.ASY.SYMBOL))
    }
  } else{
    dat = data.frame(dat)
  }
  load(file = "Data/XENA.TCGA.COL.Survival.RData")
  surv <- data.frame(t(DATA.SURV))
  dat$ID <- rownames(dat)
  surv$ID <- rownames(surv)
  dat <- merge(dat, surv, by = "ID")
  dat$Time <- as.numeric(dat$Time)
  dat$Status <- ifelse(dat$Status == "Alive", 0, 1)
  rownames(dat) =   dat$ID
  dat$ID <- NULL
  rm(surv)
  dat = subset(dat, Time > 0 )
  time = as.numeric(dat$Time)
  status = as.numeric(dat$Status)
  dat$Time = dat$Status = NULL
  dat = data.frame(t(dat))
  dat$SYMBOL = rownames(dat)
  dat = subset(dat, dat$SYMBOL %in% genes)
  dat$SYMBOL = NULL
  dat = t(dat)
  x <- dat
  y <-  Surv(time, status)
  fit <- glmnet(x, y, family = "cox")
  png(file = "Picture\\LASSO.png")
  plot(fit, label = T)
  plot(fit, xvar = "lambda", label = T)
  dev.off()
  x <- sapply(data.frame(x), as.numeric)
  lasso_fit <- cv.glmnet(data.matrix(x), y, family = "cox")
  png(file = "Picture\\LASSO_lambda.png")
  plot(lasso_fit)
  coef(lasso_fit, s = "lambda.min")
  dev.off()
  coefficient <- coef(lasso_fit, s = lasso_fit$lambda.min)
  Active.Index <- which(as.numeric(coefficient) != 0)
  active.coefficients <- as.numeric(coefficient)[Active.Index]
  sig_gene_multi_cox <- rownames(coefficient)[Active.Index]
  sigGenes = data.frame(Gene = c(rownames(coefficient)[Active.Index]), Value = c(as.numeric(coefficient)[Active.Index]))
  print(sigGenes)
  cat("Filted genes: ")
  cat(sig_gene_multi_cox)
  cat("\n")
  data <- data.frame(dat)
  x <- data.frame(x)
  formula = paste("c(", paste(paste0(sigGenes$Value," * x$",sigGenes$Gene),collapse = "+"),")")
  print(paste("Risk score formula: ", formula))
  exp = parse(text = paste("c(", paste(paste0(sigGenes$Value," * x$",sigGenes$Gene),collapse = "+"),")"))   
  data$RS <- eval(exp)
  cutoffFormula = paste(paste0(sigGenes$Value," * median(x$",sigGenes$Gene,")"),collapse = "+")
  print(paste("Cutoff formula: ", cutoffFormula))
  exp = parse(text = cutoffFormula)   
  riskScore <- eval(exp)
  cat(paste("Risk score:", riskScore, "\n"))
  RiskModel = data.frame(
    row.names = c(rownames(data)), 
    IsLowRisk = c(data$RS < riskScore),
    RiskScore = c(data$RS))
  return(list(RiskModel = RiskModel, RiskScore = riskScore,FiltedGenes = sig_gene_multi_cox, Formula = formula, Cutoff = cutoffFormula))
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Plot the time-dependent ROC curve
# ================================================================
TimeROC = function(data, title, method = "NNE"){
  library(survivalROC)
  data = data.frame(data)
  if (max(data$Time) > 1000) {
    data$Time = data$Time / 365
  }
  nobs <- NROW(data)
  roc <- survivalROC(
    Stime = data$Time,
    status = data$Status,
    marker = data$Factor,
    predict.time = 1,
    span = 0.01*nobs^(-0.20),
    method = method)
  plot(roc$FP, roc$TP,
       type = "l", col = "green", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False positive rate",
       ylab = "True positive rate",
       thresholds="best",
       print.thres="best",
       main = title)
  library(ggplotify)
  abline(a = 0, b = 1, lty = 3)
  roc3 = survivalROC(
    Stime = data$Time,
    status = data$Status,
    marker = data$Factor,
    predict.time = 3,
    span = 0.01*nobs^(-0.20),
    method = method)
  lines(roc3$FP, roc3$TP, type="l",col="yellow",xlim=c(0,1), ylim=c(0,1))
  roc5 = survivalROC(
    Stime = data$Time,
    status = data$Status,
    marker = data$Factor,
    predict.time = 5,
    span = 0.01*nobs^(-0.20),
    method = method)
  lines(roc5$FP, roc5$TP, type="l",col="red",xlim=c(0,1), ylim=c(0,1))
  cat(paste("AUC of 1 year =",round(roc$AUC, 3), "\n"))
  cat(paste("AUC of 3 years =",round(roc3$AUC,3), "\n"))
  cat(paste("AUC of 5 years =",round(roc5$AUC,3), "\n"))
  legend(0.6,0.2,c(paste("AUC of 1 year =\t",round(roc$AUC, 3)),
                   paste("AUC of 3 years =",round(roc3$AUC,3)),
                   paste("AUC of 5 years =",round(roc5$AUC,3))),
         x.intersp=1, y.intersp=0.8,
         lty= 1 ,lwd= 2,
         col=c("green","yellow","red"),
         bty = "n",
         seg.len=1,cex=0.8)
}

# ================================================================
# Function Description.
# ----------------------------------------------------------------
# Univariate + multivariate Cox regression
# ================================================================
MixCox = function(data, use = "OS", timeSpan = 10, needSurv = T){
  library("survival")
  library("survminer")
  options(scipen = 200, digits = 3)
  data = data.frame(data)
  if (needSurv) {
    surv <- GetSurvTable()
    data$PID = substr(rownames(data),1,12)
    surv$PID = substr(rownames(surv),1,12)
    surv = data.frame(
      row.names = rownames(surv), 
      Age = as.numeric(surv$Age),
      Gender = surv$Gender,
      Time = surv$Time, 
      Status = surv$Status, 
      NewTumorEvent = surv$NewTumorEvent, 
      NewTumorTime = surv$NewTumorTime,
      Stage = ifelse(surv$Stage == "I" | surv$Stage == "II", "I - II", "III - IV"),
      T = ifelse(surv$T == "T1" | surv$T == "T2", "T1 - T2", "T3 - T4"),
      N = ifelse(surv$N == "N0", "N0", "N1 - Nx"),
      M = ifelse(surv$M == "M0", "M0", "M1 - Mx"),
      Grade = ifelse(surv$Grade == "G1" | surv$Grade == "G2", "G1 - G2", "G3 - G4"),
      PID = surv$PID)
    surv = subset(surv, surv$NewTumorEvent != "NONE")
    data = merge(surv, data, by = "PID", all.y = F)
    data$PID = NULL
    if (use == "DFS") {
      data$NewTumorTime = as.numeric(data$NewTumorTime)
      data$Time = ifelse(data$NewTumorEvent == "YES", data$NewTumorTime, data$Time)
      data$Status <- ifelse(data$Status == 1, 1, ifelse(data$NewTumorEvent == "YES", 1, 0))
    } else if (use == "RFS") {
      data$Time = as.numeric(data$NewTumorTime)
      data$Status <- ifelse(data$NewTumorEvent == "YES", 1, 0)
    }
  }
  time <- as.numeric(data$Time)
  status <- as.numeric(data$Status)
  if (timeSpan != 10) {
    timeSpan = timeSpan * 365
    status = ifelse(time > timeSpan, ifelse(status == 1, 0 ,status), status)
    time = ifelse(time > timeSpan, timeSpan , time)
  }
  data$Time = data$Status = data$NewTumorEvent = data$NewTumorTime = NULL
  data = lapply(data, function(x){
    if (class(x) == "character") {
      x = as.factor(x)
    }
    return(x)
  })
  data = data.frame(data)
  formatNum = function(x, showmin = T){
    return(
      ifelse(showmin,
             ifelse(x < .001, "< 0.001", sprintf("%0.3f", x)),
             ifelse(x < .001, "0.001", sprintf("%0.3f", x))))
  }
  univ_formulas <- sapply(colnames(data), function(x) as.formula(paste("Surv(time, status) ~ data$", x)))
  univ_models <- lapply(univ_formulas, function(x) { coxph(x, data = data) })
  univ_results <- lapply(univ_models, function(x) {
    levels = x$xlevels[[1]]
    name = x$terms[[3]][[3]]
    x <- summary(x)
    p.total <- formatNum(x$wald["pvalue"])
    if (is.null(levels)) {
      HR <- formatNum(x$coef[2], F)
      HR.confint.lower <- formatNum(x$conf.int[, "lower .95"], F)
      HR.confint.upper <- formatNum(x$conf.int[, "upper .95"], F)
      res <- c(as.character(name), "", p.total, HR, HR.confint.lower, HR.confint.upper)
      names(res) <- c("name", "faturename", "p.value", "HR", "lHR", "uHR")
      return(list(res))
    }
    HR.CI95 = x[["conf.int"]]
    P = x[["coefficients"]]
    res <- c(as.character(name), "", p.total, "", "", "")
    names(res) <- c("name", "faturename", "p.value", "HR", "lHR", "uHR")
    sublevels = list()
    sublevels = append(sublevels, list(res))
    i = 1
    for (variable in levels) {
      if (i == 1) {
        res <- c("", variable, "", "", "", "")
        names(res) <- c("name", "faturename", "p.value", "HR", "lHR", "uHR")
        sublevels = append(sublevels, list(res))
      } else {
        res <- c("", variable, 
                 formatNum(P[, "Pr(>|z|)"][i - 1]), 
                 formatNum(HR.CI95[, "exp(coef)"][i - 1], F), 
                 formatNum(HR.CI95[, "lower .95"][i - 1], F), 
                 formatNum(HR.CI95[, "upper .95"][i - 1], F))
        names(res) <- c("name", "faturename", "p.value", "HR", "lHR", "uHR")
        sublevels = append(sublevels, list(res))
      }
      i = i + 1
    }
    return(sublevels)
  })
  res = list()
  needmuti = c()
  for (item in univ_results) {
    res = append(res, item)
    if (item[[1]][["p.value"]] != "") {
      if (item[[1]][["p.value"]] == "< 0.001") {
        needmuti = append(needmuti, item[[1]][["name"]])
      } else if (as.numeric(item[[1]][["p.value"]]) < .05) {
        needmuti = append(needmuti, item[[1]][["name"]])
      }
    }
  }
  res = data.frame(t(data.frame(res)))
  rownames(res) = c(1:length(rownames(res)))
  hCenter = as.numeric(res$HR)
  hLow = as.numeric(res$lHR)
  hHigh = as.numeric(res$uHR)
  CI95 = c()
  for (i in 1:length(rownames(res))) {
    r = res[i,]
    if (r["p.value"] == "") {
      CI95 = append(CI95, "Reference")
    } else if(r["faturename"] == "" & r["HR"] == ""){
      CI95 = append(CI95, "")
    } else{
      CI95 = append(CI95, paste(r["lHR"], "-", r["uHR"]))
    }
  }
  tabletext <- cbind(
    c("Item", res$name),
    c("", res$faturename),
    c("HR (95% CI for HR)", CI95),
    c("P Value", res$p.value)
  )
  library(forestplot)
  tit = ifelse(timeSpan == 10, "Univariate analysis",  paste("Univariate analysis of", use, "(", timeSpan/365 ,"years )"))
  picu <- forestplot(
    title = tit,
    xticks = seq(from = 0, to = 7, by = 1),  
    clip = c(0 , 8), 
    hrzl_lines = gpar(col="#444444"),   
    is.summary=c(TRUE,rep(FALSE,length(CI95))),  
    mean = c(NA, hCenter),
    lower = c(NA, hLow),
    upper = c(NA, hHigh),
    zero = 1,
    labeltext = tabletext,
    graph.pos = 3,
    graphwidth=unit(80,"mm"),
    vertices = TRUE,  
    grid = TRUE,      
    boxsize = 0.2,     
    new_page = FALSE
  )
  sig = needmuti
  exp = parse(text = paste("coxph(Surv(time, status) ~ ", paste(needmuti, collapse = "+"), ", data = data)"))
  result = eval(exp)
  allfactor = colnames(data)
  factors = result[["xlevels"]]
  result = summary(result)
  P = result[["coefficients"]]
  rownames(P) = sig
  HR = result[["conf.int"]]
  rownames(HR) = sig
  sublevels = list()
  i = 1
  for (f in allfactor) {
    if (f %in% sig) {
      levels = factors[[f]]
      p.total = formatNum(P[f,]["Pr(>|z|)"], T)
      if (is.null(levels)) {
        HRv <- formatNum(HR[f,]["exp(coef)"], F)
        HR.confint.lower <- formatNum(HR[f,]["lower .95"], F)
        HR.confint.upper <- formatNum(HR[f,]["upper .95"], F)
        res <- c(f, "", p.total, HRv, HR.confint.lower, HR.confint.upper)
        names(res) <- c("name", "faturename", "p.value", "HR", "lHR", "uHR")
        sublevels = append(sublevels, list(res))
      } else {
        res <- c(f, "", "", "", "", "")
        names(res) <- c("name", "faturename", "p.value", "HR", "lHR", "uHR")
        sublevels = append(sublevels, list(res))
        res <- c(f, levels[1], "", "", "", "")
        names(res) <- c("name", "faturename", "p.value", "HR", "lHR", "uHR")
        sublevels = append(sublevels, list(res))
        HRv <- formatNum(HR[f,]["exp(coef)"], F)
        HR.confint.lower <- formatNum(HR[f,]["lower .95"], F)
        HR.confint.upper <- formatNum(HR[f,]["upper .95"], F)
        res <- c(f, levels[2], p.total, HRv, HR.confint.lower, HR.confint.upper)
        names(res) <- c("name", "faturename", "p.value", "HR", "lHR", "uHR")
        sublevels = append(sublevels, list(res))
      }
    }
  }
  res = data.frame(t(data.frame(sublevels)))
  rownames(res) = c(1:length(rownames(res)))
  hCenter = c(rep(NA,length(tabletext[,1])))    
  hLow = c(rep(NA,length(tabletext[,1])))    
  hHigh = c(rep(NA,length(tabletext[,1])))    
  index = c(rep(NA,length(tabletext[,1])))
  CI95 = c(rep(NA,length(tabletext[,1])))
  p = c(rep(NA,length(tabletext[,1])))
  i = 0
  for (item in tabletext[,1]) {
    i = i + 1
    if (item != "") {
      if (item %in% res[["name"]]) {
        index[i] = item
        subD = res[res["name"] == item,]
        n = 0
        rn = length(subD[,1]) > 1
        while (length(subD[,1]) > 0) {
          if (subD[1,][3] == "") {
            p [i + n] = "-"
          } else {
            p [i + n] = subD[1,][3]
          }
          if (subD[1,][3] == "") {
            if (n != 0) {
              CI95[i + n] = "Reference"
            }
          } else {
            CI95[i + n] = paste(subD[1,][5], "-", subD[1,][6])
          }
          if (is.null(subD[1,][4])) {
            hCenter[i + n] = NA
          } else {
            hCenter[i + n] = subD[1,][4]
          }
          if (is.null(subD[1,][5])) {
            hLow[i + n] = NA
          } else {
            hLow[i + n] = subD[1,][5]
          }
          if (is.null(subD[1,][6])) {
            hHigh[i + n] = NA
          } else {
            hHigh[i + n] = subD[1,][6]
          }
          subD = subD[-1,]
          n = n + 1
        }      
        for (x in 1:n) {
          res = res[-1,]
        }
      }
    }
  }
  app = function(x){
    if (is.null(x)) {
      return(NA)
    }
    return(x)
  }
  hCenter = lapply(hCenter, app)
  hCenter = as.numeric(hCenter)
  hLow = lapply(hLow, app)
  hLow = as.numeric(hLow)
  hHigh = lapply(hHigh, app)
  hHigh = as.numeric(hHigh)
  hCenter = hCenter[2:length(hCenter)]
  hLow = hLow[2:length(hLow)]
  hHigh = hHigh[2:length(hHigh)]
  p = lapply(p, app)
  CI95 = lapply(CI95, app)
  p = unlist(p)
  CI95 = unlist(CI95)
  CI95[1] = "HR (95% CI for HR)"
  p[1] = "P Value"
  tabletext2 <- cbind(c(CI95), c(p))
  tit = ifelse(timeSpan == 10, "Multivariate analysis",  paste("Multivariate analysis of", use, "(", timeSpan/365 ,"years )"))
  picm <- forestplot(
    title = tit,
    xticks = seq(from = 0, to = 7, by = 1), 
    clip = c(0 , 8), 
    hrzl_lines = gpar(col="#444444"),   
    is.summary=c(TRUE,rep(FALSE,length(CI95))),  
    mean = c(NA, hCenter),
    lower = c(NA, hLow),
    upper = c(NA, hHigh),
    zero = 1,
    labeltext = tabletext2,
    graph.pos = 1,
    graphwidth=unit(80,"mm"),  
    vertices = TRUE,  
    grid = TRUE,      
    boxsize = 0.2,     
    new_page = FALSE   
  )
  library(grid)
  grid.newpage()
  borderWidth <- unit(4, "pt")
  width <- unit(convertX(unit(1, "npc") - borderWidth, unitTo = "npc", valueOnly = TRUE)/2, "npc")
  pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  plot(picu)
  upViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  plot(picm)
  upViewport()
}