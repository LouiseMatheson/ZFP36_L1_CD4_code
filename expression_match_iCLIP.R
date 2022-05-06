#!/usr/bin/env Rscript
Args = commandArgs(trailingOnly=TRUE)

library(tidyverse)

gene_list <- read.table(Args[grep(".txt",Args)], stringsAsFactors = F)$V1

load(Args[grep("RData",Args)])
print(ls())


#changing this so it can take CDS lengths and hits for expression matching instead - simplest way is to see if the CDS data objects exist and reassign to the UTR names
# so can only supply one or the other-  will preferentially use CDS if both were in the RData file
if(exists("iCLIP_CDS_hits")) { iCLIP_3UTR_hits <- iCLIP_CDS_hits}
if(exists("CDS_lengths")) { UTR_lengths <- CDS_lengths }

# filter gene list so that only genes for which we have both a 3'UTR length and an expression value are included (no point in considering genes without an annotated 3'UTR)
# also ensure we have the same set of genes for both of these.
UTR_lengths <- UTR_lengths[rownames(UTR_lengths) %in% rownames(Expression_norm_per_kb),]
Expression_norm_per_kb <- Expression_norm_per_kb[rownames(Expression_norm_per_kb) %in% rownames(UTR_lengths),]
gene_list <- gene_list[gene_list %in% rownames(UTR_lengths)]

geneSetName <- Args[grep(".txt",Args)]
geneSetName <-gsub("(.*).txt","\\1", geneSetName)
geneSetName <-gsub(".*/(.*?)$","\\1", geneSetName)

Celltype <- Args[grep(".RData",Args)]
Celltype <- gsub("(.*).RData", "\\1", Celltype)
Celltype <- gsub(".*/(.*?)$", "\\1", Celltype)
                                  
print(geneSetName)
print(Celltype)


### Define functions for expression/3'UTR length matching
#takes counts data (should only contain counts unless baseMeanGrep specified which only selects columns with counts), some form of processed data (resData) which includes a column (resCol) containing the assignment of genes into gene sets. 
#Will then identify expression matched genes for required gene sets, from the remaining genes that are not in any of the gene sets specified
# Will assume that gene names are rownames (of both countData - which must be the case - and resData), unless geneCol != NA in which case will look for column entitled geneCol specified
# Instead of providing a column in resData containing assigment to gene sets, can provide the genes within each set - in this case should be provided as a list and resCol left as NA
# will calculate base mean for countData across
# Default for new column containing expression matched (plus original gene set) assignment is "Expression_matched_genes"
# also adding option to consider 3'UTR length: will look for all genes no more than 25% different; or if this is less than 1000 genes will look for closest 1000 (but with min cutoff of 10bp)
# UTRdata should have rownames and a single column containing UTR lengths
ExpressionMatch <- function(countData, resData, geneSets, geneSetNames = NA, resCol = NA, geneCol = NA, nRandom = 100, 
                            baseMeanGrep = "", newCol = "Expression_matched_genes", avoid_genes = NULL, UTR3p = F, UTRdata = NULL, 
                            UTRthreshold = 10){
  countData <- as.data.frame(countData)
  resData <- as.data.frame(resData)
  countData$baseMean <- rowMeans(countData[,grepl(baseMeanGrep, colnames(countData))])
  geneSet_counts <- list()
  if(is.na(geneCol)){
    countData <- countData[rownames(countData) %in% rownames(resData),]
  }
  else {
    countData <- countData[rownames(countData) %in% resData$geneCol,]
  }
  if(is.na(resCol)){
    resData[,newCol] <- NA
    if(is.na(geneSetNames)){
      geneSetNames <- paste("geneSet",1:length(geneSets),sep = "_")
    }
    for(i in 1:length(geneSets)){
      if(is.na(geneCol)){
        resData[rownames(resData) %in% geneSets[[i]],newCol] <- geneSetNames[i]
      }
      else {
        resData[resData[,geneCol] %in% geneSets[[i]],newCol] <- geneSetNames[i]
      }
    }
    for(i in 1:length(geneSets)){
      genes_in_set <- geneSets[[i]]
      geneSet_counts[[i]] <- countData[rownames(countData) %in% genes_in_set,]
      countData <- countData[!rownames(countData) %in% genes_in_set,]
    }
  }
  else {
    resData[,newCol] <- as.character(resData[,resCol])
    if(is.na(geneSetNames)){
      geneSetNames <- geneSets
    }
    for(i in 1:length(geneSets)){
      resData[,newCol] <- gsub(geneSets[i], geneSetNames[i], resData[,newCol])
      if(is.na(geneCol)){
        genes_in_set <- rownames(resData)[resData[,resCol] == geneSets[i]]
      }
      else {
        genes_in_set <- resData[resData[,resCol] == geneSets[i], geneCol]
      }
      geneSet_counts[[i]] <- countData[rownames(countData) %in% genes_in_set,]
      countData <- countData[!rownames(countData) %in% genes_in_set,]
    }
  }
  countData <- countData[!rownames(countData) %in% avoid_genes,]
  if(UTR3p == T){
    geneSet_UTRs <- list()
    for(i in 1:length(geneSets)){
      geneSet_UTRs[[i]] <- UTRdata[match(rownames(geneSet_counts[[i]]),rownames(UTRdata)),1]
    }
    UTRdata <- na.omit(UTRdata)
    UTRdata$threshold <- UTRdata[,1] > UTRthreshold
    UTRdata <- UTRdata[UTRdata$threshold == T,]
    countData <- countData[rownames(countData) %in% rownames(UTRdata),]
  }
  for(i in 1:length(geneSets)){
    expression_matched_genes <- character()
    for(j in 1:nrow(geneSet_counts[[i]])){
      if(UTR3p==T){
        UTRdata1 <- UTRdata[rownames(UTRdata) %in% rownames(countData),]
        utr_length <- geneSet_UTRs[[i]][j]
        utr_min <- floor(0.75*utr_length)
        utr_max <- ceiling(1.25*utr_length)
        genes_to_keep <- rownames(UTRdata1)[UTRdata1[,1] > utr_min & UTRdata1[,1] < utr_max]
        if(length(genes_to_keep) < 1000){
          UTRdata1$diff <- abs(UTRdata1[,1] - utr_length)
          genes_to_keep <- rownames(UTRdata1[order(UTRdata1$diff),])[1:1000]
        }
        countData1 <- countData[rownames(countData) %in% genes_to_keep,]
      }
      else {
        countData1 <- countData
      }
      closest_gene_index <- sample(order(abs(countData1$baseMean - geneSet_counts[[i]]$baseMean[j]))[1:nRandom],1)
      expression_matched_genes[j] <- rownames(countData1)[closest_gene_index]
      countData <- countData[!rownames(countData) == rownames(countData1)[closest_gene_index],]
    }
    if(is.na(geneCol)){
      resData[rownames(resData) %in% expression_matched_genes,newCol] <- paste("Expression matched to",geneSetNames[i])
    }
    else {
      resData[resData[,geneCol] %in% expression_matched_genes,newCol] <- paste("Expression matched to",geneSetNames[i])
    }
  }
  return(resData)
}

# To generate multiple set of expression matched genes:
ExpressionMatch_multi <- function(nSets, countData, resData, geneSets, geneSetNames = NA, resCol = NA, geneCol = NA, nRandom = 100, 
                                  baseMeanGrep = "", newCol = "Expression_matched_genes", avoid_genes = NULL, UTR3p = F, UTRdata = NULL, 
                                  UTRthreshold = 10){
  if(is.na(geneSetNames)){
    if(is.na(resCol)){
      geneSetNames <- paste("geneSet",1:length(geneSets),sep = "_")
    }
    else {
      geneSetNames <- geneSets
    }
  }
  expression_matched_sets <- list()
  for(s in 1:nSets){
    expression_matched_sets[[s]] <- list()
    eM <- ExpressionMatch(countData = countData, resData = resData, geneSets = geneSets, geneSetNames = geneSetNames, resCol = resCol, 
                          geneCol = geneCol, nRandom = nRandom, baseMeanGrep = baseMeanGrep, newCol = newCol, avoid_genes = avoid_genes, 
                          UTR3p = UTR3p, UTRdata = UTRdata, UTRthreshold = UTRthreshold)
    eM_levels <- levels(na.omit(as.factor(eM[,newCol])))
    eM_levels <- eM_levels[!eM_levels %in% geneSetNames]
    for(l in eM_levels){
      if(is.na(geneCol)){
        expression_matched_sets[[s]][[as.character(l)]] <- as.character(na.omit(rownames(eM)[eM[,newCol]== l]))
      }
    }
  }
  return(expression_matched_sets)
}


expression_matched <-ExpressionMatch_multi(nSets = 100, countData = Expression_norm_per_kb, resData = Expression_norm_per_kb, geneSets = list(gene_list), geneSetNames = geneSetName, nRandom = 100, UTR3p = T, UTRdata = UTR_lengths)


exMatch_iCLIP3UTR_overlap <- sapply(expression_matched, function(x) sum(x[[1]] %in% iCLIP_3UTR_hits))


iCLIP_hits <- data.frame(Gene = gene_list)
iCLIP_hits$iCLIP_hit <- iCLIP_hits$Gene %in% iCLIP_3UTR_hits


exMatch_overlaps <- data.frame(Overlap = exMatch_iCLIP3UTR_overlap)
exMatch_overlaps$upper_95percent <- rep(quantile(exMatch_iCLIP3UTR_overlap, 0.95),length(exMatch_iCLIP3UTR_overlap))
exMatch_overlaps$lower_5percent <- rep(quantile(exMatch_iCLIP3UTR_overlap, 0.05),length(exMatch_iCLIP3UTR_overlap))

pdf(paste0(geneSetName, "_", Celltype, ".pdf"))
ggplot(iCLIP_hits, aes(x = geneSetName))+geom_bar(aes(fill = iCLIP_hit)) + scale_fill_manual(values =c("grey60","mediumturquoise")) + stat_summary(data = exMatch_overlaps, fun = "median", geom = "bar", width=0.7, fill = "darkcyan", aes(y = Overlap))+geom_errorbar(data = exMatch_overlaps[1,],aes(ymin= lower_5percent, ymax =upper_95percent), width=0.5)+ labs(y = "Number of genes", title = paste0(geneSetName,"\n",Celltype))
dev.off()

Expression_norm_per_kb$Mean <- rowMeans(Expression_norm_per_kb)
exMatch_expression <- data.frame(Gene = gene_list, Expression = Expression_norm_per_kb$Mean[match(gene_list, rownames(Expression_norm_per_kb))], Gene_set = geneSetName, Type = geneSetName)

for(i in 1:length(expression_matched)){
  exMatch_expression <- rbind(exMatch_expression, data.frame(Gene = expression_matched[[i]][[1]], Expression = Expression_norm_per_kb$Mean[match(expression_matched[[i]][[1]], rownames(Expression_norm_per_kb))], Gene_set = paste0("EM_", i), Type ="Expression matched"))
}

for(i in 1:10){
  random_genes <- sample(rownames(Expression_norm_per_kb),180)
  exMatch_expression <- rbind(exMatch_expression, data.frame(Gene = random_genes, Expression = Expression_norm_per_kb$Mean[match(random_genes, rownames(Expression_norm_per_kb))], Gene_set = paste0("Random_", i), Type = "Random"))
}

exMatch_expression$Type <- factor(exMatch_expression$Type, levels = c(geneSetName,"Expression matched", "Random"))
exMatch_expression$Gene_set <- as.factor(exMatch_expression$Gene_set)
exMatch_expression$Gene_set <- relevel(exMatch_expression$Gene_set, geneSetName)


### How about 3'UTR length?
exMatch_expression$UTRlength <- UTR_lengths[match(exMatch_expression$Gene, rownames(UTR_lengths)),1]

pdf(paste0(geneSetName, "_", Celltype, "_expression_UTR_match.pdf"), width=9)
ggplot(exMatch_expression, aes(x = Gene_set, y = Expression, fill = Type)) + geom_boxplot(outlier.size =0) + scale_fill_manual(values = c("firebrick","lightblue1", "lightgoldenrod"))
ggplot(exMatch_expression, aes(x = Gene_set, y = UTRlength, fill = Type)) + geom_boxplot(outlier.size =0)+scale_y_log10() + scale_fill_manual(values = c("firebrick","lightblue1", "lightgoldenrod"))
dev.off()

exMatch_expression$iCLIP_hit_FDR0.05 <- exMatch_expression$Gene %in% iCLIP_3UTR_hits

write.table(exMatch_expression, file = paste0(geneSetName, "_",Celltype,"_expression_matched.txt"), sep = "\t", quote = F, row.names =F)

exMatch_processed_data <- list(iCLIP_hits=iCLIP_hits, exMatch_overlaps= exMatch_overlaps, exMatch_expression=exMatch_expression)
saveRDS(exMatch_processed_data, paste0(geneSetName, "_", Celltype, "_expressionMatch_processed_data.rds"))
