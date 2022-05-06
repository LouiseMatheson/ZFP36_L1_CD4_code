#!/usr/bin/env Rscript

# Expects a single gtf file as an argument (if >1 will use first listed), and 
# one or more tab delimited text files containing 'Gene' column containing 
# Ensembl gene IDs corresponding to the GTF provided. 
Args = commandArgs(trailingOnly=TRUE)

library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
genome <- BSgenome.Mmusculus.UCSC.mm10

gtf_file <- Args[grepl("gtf$|gtf.gz$", Args)][1] # if more than one provided, use first
Args <- Args[!grepl("gtf$|gtf.gz$", Args)]

gtf <- import(gtf_file)

# define function for finding locations of sequence motifs
findMotif <- function(seqData, motif, baseCol ="base", densWindow = 25, dist = "overall", non_overlapping = T){
  seq <- paste(seqData[,baseCol], collapse = "")
  motif_indices <- gregexpr(paste0("(?=",motif,")"), seq, perl=T)[[1]]
  seqData[,motif] <- F
  nc <- nchar(motif)
  if(motif_indices[1]>0){
    for(i in 1:length(motif_indices)){
      seqData[motif_indices[i],motif] <- TRUE
    }
  }
  if(length(motif_indices)>1 & !(dist == "none")){
    seqRLE <- rle(seqData[,motif])
    rleT <- which(seqRLE$values == T) # this gives indices of seqRLE that are T (ie where motif is), but if there were 2xT in a row (only poss for homonuc motif) wouldn't be apparent
    rleI <- integer()
    for(t in 1:length(rleT)){
      rleI <- c(rleI, rep(rleT[t], seqRLE$lengths[rleT[t]])) #should give same number of indices for RLE as there are motif indices
    }
    if(dist == "both"|dist == "upstream") {
      seqData[,paste0(motif,"_dist_upstream")] <- NA
      
      if(rleI[length(rleI)] == rleI[length(rleI)-1]) {
        seqData[motif_indices[length(motif_indices)],paste0(motif,"_dist_upstream")] <- 1-nc
      } else {
        seqData[motif_indices[length(motif_indices)],paste0(motif,"_dist_upstream")] <- seqRLE$lengths[rleI[length(rleI)] - 1]+1-nc
      }
      if(non_overlapping==T){ 
        seqData[,paste0(motif,"_dist_upstream_pos")] <- NA
        if(seqData[motif_indices[length(motif_indices)],paste0(motif,"_dist_upstream")] >= 0){
          seqData[motif_indices[length(motif_indices)],paste0(motif,"_dist_upstream_pos")] <- seqData[motif_indices[length(motif_indices)],paste0(motif,"_dist_upstream")]
        } else {
          seqData[motif_indices[length(motif_indices)],paste0(motif,"_dist_upstream_pos")] <- min(((motif_indices[length(motif_indices)]-nc) - motif_indices[1:length(motif_indices)-1])[((motif_indices[length(motif_indices)]-nc) - motif_indices[1:length(motif_indices)-1]) >= 0])
        }
      }
    }
    if(dist == "both"|dist == "downstream"){
      seqData[,paste0(motif,"_dist_downstream")] <- NA
      if(rleI[2] == rleI[1]) {
        seqData[motif_indices[1],paste0(motif,"_dist_downstream")] <- 1-nc
      }
      else {
        seqData[motif_indices[1],paste0(motif,"_dist_downstream")] <- seqRLE$lengths[rleI[1] + 1]+1-nc
      }
      if(non_overlapping==T){ 
        seqData[,paste0(motif,"_dist_downstream_pos")] <- NA
        if(seqData[motif_indices[1],paste0(motif,"_dist_downstream")] >= 0){
          seqData[motif_indices[1],paste0(motif,"_dist_downstream_pos")] <- seqData[motif_indices[1],paste0(motif,"_dist_downstream")]
        } else {
          seqData[motif_indices[1],paste0(motif,"_dist_downstream_pos")] <- min(((motif_indices[2:length(motif_indices)]-motif_indices[1])-nc)[((motif_indices[2:length(motif_indices)]-motif_indices[1])-nc) >= 0] )
        }
      }
    }
    if(dist == "overall"){
      seqData[,paste0(motif,"_dist_overall")] <- NA
      if(rleI[length(rleI)] == rleI[length(rleI)-1]) {
        seqData[motif_indices[length(motif_indices)],paste0(motif,"_dist_overall")] <- 1-nc
      }
      else {
        seqData[motif_indices[length(motif_indices)],paste0(motif,"_dist_overall")] <- seqRLE$lengths[rleI[length(rleI)] - 1]+1-nc
      }
      if(non_overlapping==T){ 
        if(seqData[motif_indices[length(motif_indices)],paste0(motif,"_dist_overall")] >= 0){
          seqData[motif_indices[length(motif_indices)],paste0(motif,"_dist_overall_pos")] <- seqData[motif_indices[length(motif_indices)],paste0(motif,"_dist_overall")]
        } else {
          seqData[motif_indices[length(motif_indices)],paste0(motif,"_dist_overall_pos")] <- min(((motif_indices[length(motif_indices)]-nc) - motif_indices[1:length(motif_indices)-1])[((motif_indices[length(motif_indices)]-nc) - motif_indices[1:length(motif_indices)-1]) >= 0])
        }
      }
      if(rleI[2] == rleI[1]) {
        seqData[motif_indices[1],paste0(motif,"_dist_overall")] <- 1-nc
      }
      else {
        seqData[motif_indices[1],paste0(motif,"_dist_overall")] <- seqRLE$lengths[rleI[1] + 1]+1-nc
      }
      if(non_overlapping==T){ 
        if(seqData[motif_indices[1],paste0(motif,"_dist_overall")] >= 0){
          seqData[motif_indices[1],paste0(motif,"_dist_overall_pos")] <- seqData[motif_indices[1],paste0(motif,"_dist_overall")]
        } else {
          seqData[motif_indices[1],paste0(motif,"_dist_overall_pos")] <- min(((motif_indices[2:length(motif_indices)]-motif_indices[1])-nc)[((motif_indices[2:length(motif_indices)]-motif_indices[1])-nc) >= 0] )
        }
      }
    }
    if(length(motif_indices)>2){
      for(i in 2:(length(motif_indices)-1)){
        if(rleI[i] == rleI[i-1]){
          upstream_dist <- 1-nc
        }
        else {
          upstream_dist <- seqRLE$lengths[rleI[i] - 1]+1-nc
        }
        if(rleI[i] == rleI[i+1]){
          downstream_dist <- 1-nc
        }
        else {
          downstream_dist <- seqRLE$lengths[rleI[i] + 1]+1-nc
        }
        if(dist == "both"|dist == "upstream") {
          seqData[motif_indices[i],paste0(motif,"_dist_upstream")] <- upstream_dist
        }
        if(dist == "both"|dist == "downstream") {
          seqData[motif_indices[i],paste0(motif,"_dist_downstream")] <- downstream_dist
        }
        if(dist == "overall") {
          seqData[motif_indices[i],paste0(motif,"_dist_overall")] <- min(downstream_dist, upstream_dist)
        }
        if(non_overlapping==T){ 
          if(downstream_dist < 0){
            downstream_dist <- min(((motif_indices[(i+1):length(motif_indices)]-motif_indices[i])-nc)[((motif_indices[(i+1):length(motif_indices)]-motif_indices[i])-nc) >= 0] )
          }
          if(upstream_dist < 0){
            upstream_dist <- min(((motif_indices[i]-nc) - motif_indices[1:(i-1)])[((motif_indices[i]-nc) - motif_indices[1:(i-1)]) >= 0])
          }
          if(dist == "both"|dist == "upstream") {
            seqData[motif_indices[i],paste0(motif,"_dist_upstream_pos")] <- upstream_dist
          }
          if(dist == "both"|dist == "downstream") {
            seqData[motif_indices[i],paste0(motif,"_dist_downstream_pos")] <- downstream_dist
          }
          if(dist == "overall") {
            seqData[motif_indices[i],paste0(motif,"_dist_overall_pos")] <- min(downstream_dist, upstream_dist)
          }
        }
      }
    }
  }
  else {
    seqData[,paste0(motif,"_dist_",dist)] <- NA
    if(non_overlapping==T){ seqData[,paste0(motif,"_dist_",dist, "_pos")] <- NA }
  }
  if(length(densWindow) > 0 & motif_indices[1]>0){
    for(W in densWindow){
      for(i in 1:nrow(seqData)){
        seqData[i,paste(motif,"dens",as.character(W), sep ="_")] <- sum(seqData[,motif][max(1,(((i-W)-(nc-1)))):min((i+W),nrow(seqData))])
      }
    }
  }
  else {
    for(W in densWindow){
      seqData[,paste0(motif,"_dens_",as.character(W))] <- 0
    }
  }
  if(motif_indices[1]>0){
    for(i in 1:(nc-1)){
      MI <- motif_indices + i
      for(mi in MI){
        if(seqData[mi,motif] == F){
          seqData[mi,grep(paste0(motif, "_dist_"), colnames(seqData))] <- seqData[(mi-i),grep(paste0(motif, "_dist_"), colnames(seqData))] 
        }
      }
    }
  }
  return(seqData)
}


for(a in Args){
  df <- read_tsv(a)
  
  df$nTATT <- 0
  df$nTATTTATT <- 0
  df$TATT_minDist <- NA
  df$TATTTATT_minDist <- NA
  df$transcript_id <- NA
  
  genes_without_3UTR <- character()

  for(gene in unique(df$Gene)){
    utrCoord <- gtf[gtf$gene_id == gene & gtf$type == "three_prime_utr"]
    
    if(length(utrCoord) == 0){
      print(paste("No 3'UTR coordinates not found for", gene, "\nSkipping to next gene"))
      genes_without_3UTR <- c(genes_without_3UTR,gene)
      next
    }
    
    # we ideally only want one 3'UTR to analyse, since includes counting number of motifs etc. 
    # select based on whether the isoform is protein coding (vs eg NMD), whether it has a CCDS ID, and then prioritise best TSL.
    # if tied, then check whether overlapping, analyse the longest, or if more than one that don't overlap analyse both separately and take highest (or 'best' eg lowest min distance) to report
    as_tibble(utrCoord) %>%
      group_by(transcript_id) %>%
      summarise(utrLength = sum(width), protein_coding = gene_biotype[1] == "protein_coding", TSL = as.numeric(sub(" .*$", "", transcript_support_level[1])), CCDS = !is.na(ccds_id[1])) -> transcripts_with_3UTR
    
    if(TRUE %in% transcripts_with_3UTR$protein_coding) {
      transcripts_with_3UTR %>%
        filter(protein_coding == T) -> transcripts_with_3UTR
    }
    if(TRUE %in% transcripts_with_3UTR$CCDS) {
      transcripts_with_3UTR %>%
        filter(CCDS == T) -> transcripts_with_3UTR
    }
    if(sum(is.na(transcripts_with_3UTR$TSL)) < nrow(transcripts_with_3UTR)) { # if all are NA, just keep all
      transcripts_with_3UTR %>%
        filter(!is.na(TSL)) %>%
        filter(TSL == min(TSL)) -> transcripts_with_3UTR
    }
    
    if(nrow(transcripts_with_3UTR) > 1) {
      transcripts_with_3UTR %>%
        arrange(desc(utrLength)) %>%
        add_column(eliminate = F) -> transcripts_with_3UTR
      for(i in 1:(nrow(transcripts_with_3UTR)-1)) {
        if(transcripts_with_3UTR$eliminate[i] == F){
          for(j in (i+1):nrow(transcripts_with_3UTR)){
            tid_i <- transcripts_with_3UTR$transcript_id[i]
            tid_j <- transcripts_with_3UTR$transcript_id[j]
            if(length(findOverlaps(utrCoord[utrCoord$transcript_id == tid_i], utrCoord[utrCoord$transcript_id == tid_j])) > 0) {
              transcripts_with_3UTR$eliminate[j] <- T
            }
          }
        }
      }
      transcripts_with_3UTR <- filter(transcripts_with_3UTR, eliminate == F)
    }
    
    as_tibble(utrCoord) %>%
      arrange(start) %>%
      filter(transcript_id == transcripts_with_3UTR$transcript_id[1]) -> utrCoord_subset
    chr <- paste0("chr",utrCoord_subset$seqnames[1])
    strand <- utrCoord_subset$strand[1]
    if(strand == "-") {
      utrCoord_subset %>%
        arrange(desc(start)) -> utrCoord_subset
    }
    utrSeq <- character()
    utrPositions <- integer()
    for(r in 1:nrow(utrCoord_subset)){
      utrSeq <- paste0(utrSeq, getSeq(genome, chr, utrCoord_subset$start[r], utrCoord_subset$end[r], strand = strand, as.character = T))
      if(strand == "+") {
        utrPositions <- c(utrPositions, c(utrCoord_subset$start[r]:utrCoord_subset$end[r]))
      } else {
        utrPositions <- c(utrPositions, c(utrCoord_subset$end[r]:utrCoord_subset$start[r]))
      }
    }
    utr3_seq <- data.frame(position = utrPositions, base = strsplit(utrSeq, "")[[1]])
    utr3_seq$base <- factor(utr3_seq$base, levels = c("A","C","G","T","N"))
    utr3_seq <- findMotif(utr3_seq, motif = "TATT")
    utr3_seq <- findMotif(utr3_seq,motif = "TATTTATT")
    df$nTATT[df$Gene == gene] <- sum(utr3_seq$TATT)
    df$nTATTTATT[df$Gene == gene] <- sum(utr3_seq$TATTTATT)
    df$TATT_minDist[df$Gene == gene] <- min(na.omit(utr3_seq$TATT_dist_overall_pos))
    df$TATTTATT_minDist[df$Gene == gene] <- min(na.omit(utr3_seq$TATTTATT_dist_overall_pos))
    df$transcript_id[df$Gene == gene] <- transcripts_with_3UTR$transcript_id[1]
    
    if(nrow(transcripts_with_3UTR) > 1) {
      for(i in 2:nrow(transcripts_with_3UTR)) {
        as_tibble(utrCoord) %>%
          arrange(start) %>%
          filter(transcript_id == transcripts_with_3UTR$transcript_id[1]) -> utrCoord_subset
        chr <- paste0("chr",utrCoord_subset$seqnames[1])
        strand <- utrCoord_subset$strand[1]
        if(strand == "-") {
          utrCoord_subset %>%
            arrange(desc(start)) -> utrCoord_subset
        }
        utrSeq <- character()
        utrPositions <- integer()
        for(r in 1:nrow(utrCoord_subset)){
          utrSeq <- paste0(utrSeq, getSeq(genome, chr, utrCoord_subset$start[r], utrCoord_subset$end[r], strand = strand, as.character = T))
          if(strand == "+") {
            utrPositions <- c(utrPositions, c(utrCoord_subset$start[r]:utrCoord_subset$end[r]))
          } else {
            utrPositions <- c(utrPositions, c(utrCoord_subset$end[r]:utrCoord_subset$start[r]))
          }
        }
        utr3_seq <- data.frame(position = utrPositions, base = strsplit(utrSeq, "")[[1]])
        utr3_seq$base <- factor(utr3_seq$base, levels = c("A","C","G","T","N"))
        utr3_seq <- findMotif(utr3_seq, motif = "TATT")
        utr3_seq <- findMotif(utr3_seq,motif = "TATTTATT")
        nTATT <- sum(utr3_seq$TATT)
        nTATTTATT <- sum(utr3_seq$TATTTATT)
        minDist_TATT <- min(na.omit(utr3_seq$TATT_dist_overall_pos))
        minDist_TATTTATT <- min(na.omit(utr3_seq$TATTTATT_dist_overall_pos))
        
        # if these results are 'better' than with the previous isoform tested, replace in the df
        if(nTATTTATT > df$nTATTTATT[df$Gene == gene] | is.na(df$TATT_minDist[df$Gene == gene])) {
          df$nTATT[df$Gene == gene] <- nTATT
          df$nTATTTATT[df$Gene == gene] <- nTATTTATT
          df$TATT_minDist[df$Gene == gene] <- minDist_TATT
          df$TATTTATT_minDist[df$Gene == gene] <- minDist_TATTTATT
          df$transcript_id[df$Gene == gene] <- transcript_ids[i]
        } else if(!(is.na(df$TATT_minDist[df$Gene == gene]) | is.na(minDist_TATT))) {
          if(minDist_TATT < df$TATT_minDist[df$Gene == gene]) {
            df$nTATT[df$Gene == gene] <- nTATT
            df$nTATTTATT[df$Gene == gene] <- nTATTTATT
            df$TATT_minDist[df$Gene == gene] <- minDist_TATT
            df$TATTTATT_minDist[df$Gene == gene] <- minDist_TATTTATT
            df$transcript_id[df$Gene == gene] <- transcript_ids[i]
          }
        }
      }
    }
  }
    
  print("Genes for which no 3'UTR sequence found:")
  print(genes_without_3UTR)
  write.table(df, file = paste0("TATT_",a), row.names = F, quote =F, sep = "\t")
}
