library(rBLAST)
library(rtracklayer)
library(Rsamtools)

########################################################################

#functions

########################################################################

orderStartEnd <- function(hg, startCol, endCol){
  
  check <- sapply(1:nrow(hg), function(x){
    return((hg[x,startCol] > hg[x,endCol]))
  })
  
  if (TRUE %in% check){
    new_start <- rep(0,nrow(hg))
    new_end <- rep(0, nrow(hg))
    
    for (i in 1:nrow(hg)){
      if (hg[i,startCol] > hg[i,endCol]){
        new_start[i] <- hg[i,endCol]
        new_end[i] <- hg[i,startCol]
      } else{
        new_start[i] <- hg[i,startCol]
        new_end[i] <- hg[i,endCol]
      }
    }
    
    hg[,startCol] <- new_start
    hg[,endCol] <- new_end
  }
  
  return(hg)
}

stitch <- function(sub){
  
  #find what hits could be combined into 
  #a candidate homologous gene
  
  #assuming that the genes arent overlapping
  
  #special case where there is only one row
  if (nrow(sub) == 1){
    group <- paste(as.character(sub$SubjectID[1]),'1',sep=':')
    sub$group <- group
    
    return(sub)
  }
  
  #first check if the sub dataframe contains a warning column
  #if no than add it in
  
  sub <- sub[order(sub$S.start,decreasing = F),]
  
  #check if there are any patterns going down the sequence
  
  groupD <- rep(0,nrow(sub))
  index = 1
  
  q_s <- sub$Q.start[index]
  s_s <- sub$S.start[index]
  groupD[1] <- index
  
  for (i in 2:nrow(sub)){
    
    if ((sub$Q.start[i] > q_s) & (sub$S.start[i] > s_s)){
      groupD[i] <- index
    } else{
      index <- index + 1
      q_s <- sub$Q.start[i]
      s_s <- sub$S.start[i]
      groupD[i] <- index
    }
  }
  
  uD <- unique(groupD)
  lensD <- rep(0,length(uD))
  for (i in 1:length(uD)){
    lensD[i] <- sum(sub$Alignment.Length[which(groupD == uD[i])])
  }
  
  
  #check if there are any patterns going up the sequence
  
  groupU <- rep(0,nrow(sub))
  i = nrow(sub)
  index <- index + 1
  
  q_s <- sub$Q.start[i]
  s_s <- sub$S.start[i]
  groupU[i] <- index
  
  for (i in (nrow(sub)-1):1 ){
    if ((sub$Q.start[i] > q_s) & (sub$S.start[i] < s_s)){
      groupU[i] <- index
    } else{
      index <- index + 1
      q_s <- sub$Q.start[i]
      s_s <- sub$S.start[i]
      groupU[i] <- index
    }
  }
  
  uU <- unique(groupU)
  lensU <- rep(0,length(uU))
  for (i in 1:length(uU)){
    lensU[i] <- sum(sub$Alignment.Length[which(groupU == uU[i])])
  }
  
  #choose the group that makes the segment part of a larger gene
  
  #make sure that it is true for all regions in the group...
  #there may be overlap in some genes
  
  group <- rep(0,nrow(sub))
  
  for (i in 1:nrow(sub)){
    valueU <- lensU[which(uU == groupU[i])]
    valueD <- lensD[which(uD == groupD[i])]
    if (valueU > valueD){
      group[i] <- groupU[i]
    }else{
      group[i] <- groupD[i]
    }
  }
  
  #get the coverage
  group <- paste(as.character(sub$SubjectID[1]),group,sep=':')
  sub$group <- group
  
  return(sub)
}

predictHomologs <- function(seq, bl_subject){
  
  #blast the sequences against genome
  
  #Run BLAST query
  hits <- predict(bl_subject, seq)
  
  #get the sequences from the genome 
  #and realign them?
  
  hits <- orderStartEnd(hits, 'S.start', 'S.end')
  
  #get the length of the query coverage per sequence and any overlaps
  
  #for each sequence with hits 
  sequences_included <- as.character(unique(hits$SubjectID)) #all the sequences with hits
  #coverage <- rep(0,length(sequences_included))
  
  
  h <- 1
  sub <- hits[which(hits$SubjectID == sequences_included[h]),]
  
  sub <- stitch(sub)#find what hits could be combined into a candidate homologous gene
  tot_sub <- sub
  
  if (length(sequences_included) > 1){
    for (h in 2:length(sequences_included)){
      sub <- hits[which(hits$SubjectID == sequences_included[h]),]
      
      sub <- stitch(sub)#find what hits could be combined into a candidate homologous gene
      tot_sub <- rbind(tot_sub,sub)
    }
  }
  
  return(tot_sub)
}

getCDSseq <- function(homologs, fa_subject){
  
  groups <- unique(homologs$group)
  
  seqs <- sapply(groups, function(h){
    
    sub_sub <- homologs[which(homologs$group == h),]
    
    strands <- sapply(1:nrow(sub_sub), function(x){
      if (sub_sub$S.start[x] <= sub_sub$S.end[x]){
        return('+')
      }else{
        return('-')
      }
    })
    
    gr_sub <- GRanges(
      seqnames = sub_sub$SubjectID,
      ranges = IRanges(sub_sub$S.start, end = sub_sub$S.end),
      strand = strands)
    
    seq1 <- getSeq(fa_subject, gr_sub)
    seq1 <- DNAStringSet(paste(as.character(seq1), sep='', collapse = ''))
    names(seq1) <- h
    
    return(seq1)
  })
  
  return(seqs)
  
}


getCDSseqGFF <- function(hg, fa_ref){
  
  gr <- GRanges(
    seqnames = hg$seqid,
    ranges = IRanges(hg$start, end = hg$end),
    strand = hg$strand)
  
  seq <- getSeq(fa_ref, gr)
  #names(seq) <- as.character(1:length(seq)) #rename for testing
  seq <- DNAStringSet(paste(as.character(seq), sep='', collapse = ''))
  names(seq) <- 'query'
  
  return(seq)
}

compressStitches <- function(homologs){
  
  groups <- unique(homologs$group)
  
  Num_Seg <- sapply(1:length(groups),function(x){
    return(length(which(homologs$group == groups[x])))
  })
  
  SubjectID <- sapply(1:length(groups),function(x){
    indexes <- which(homologs$group == groups[x])
    return(unique(homologs$SubjectID[indexes]))
  })
  
  S.start <- sapply(1:length(groups),function(x){
    indexes <- which(homologs$group == groups[x])
    return(min(homologs$S.start[indexes]))
  })
  
  S.end <- sapply(1:length(groups),function(x){
    indexes <- which(homologs$group == groups[x])
    return(max(homologs$S.end[indexes]))
  })
  
  Total_Mismatches <- sapply(1:length(groups),function(x){
    indexes <- which(homologs$group == groups[x])
    return(sum(homologs$Mismatches[indexes]))
  })
  
  Total_Gap_Openings <- sapply(1:length(groups),function(x){
    indexes <- which(homologs$group == groups[x])
    return(sum(homologs$Mismatches[indexes]))
  })
  
  Total_Alignment_Length <- sapply(1:length(groups),function(x){
    indexes <- which(homologs$group == groups[x])
    return(sum(homologs$Alignment.Length[indexes]))
  })
  
  Total_Percent_Ident <- sapply(1:length(groups),function(x){
    len <- Total_Alignment_Length[x]
    return( (len - Total_Mismatches[x] - Total_Gap_Openings[x])/len )
  })
  
  c_homologs <- data.frame(groups, SubjectID, Num_Seg, Total_Percent_Ident, Total_Alignment_Length, Total_Mismatches, Total_Gap_Openings, S.start, S.end)
  return(c_homologs)
}

getOverlappingFeatures <- function(c_homologs, gene_pos){
  query <- GRanges(seqnames = c_homologs$SubjectID, ranges = IRanges(c_homologs$S.start,c_homologs$S.end))  
  subject <- GRanges(seqnames = gene_pos$seqid, ranges = IRanges(gene_pos$start,gene_pos$end))
  
  ov <- findOverlaps(query,subject)
  osv <- subjectHits(ov)
  oqv <- queryHits(ov)
  
  c_homologs$hits <- rep('',nrow(c_homologs))
  
  for (i in 1:nrow(c_homologs)){
    indexes <- which(oqv == i)
    c_homologs$hits[i] <- paste(as.character(gene_pos$gene[osv[indexes]]),sep='', collapse=', ')
  }
  
  return(c_homologs)
}

checkStats <- function(seqs, seq, c_homologs, ratio_min = 0.75, ratio_max = 1.25){
  
  c_homologs$size_diff <- sapply(1:length(seqs), function(h){
    ratio <- abs(width(seqs[[h]]) - width(seq))
    return(ratio)
  })
  
  # - similar size?
  c_homologs$size_check <- sapply(1:length(seqs), function(h){
    ratio <- width(seqs[[h]])/width(seq)
    return((ratio > ratio_min) & (ratio < ratio_max))
  })
  
  #get a score
  c_homologs$score <- sapply(1:length(seqs), function(h){
    len <- width(seq)
    size_diff_score <- (len - c_homologs$size_diff[h])/len
    return(sum(c(c_homologs$Total_Percent_Ident[h],size_diff_score))/2)
  })
  
  return(c_homologs)
}

checkBackBlast <- function(top_homologs, top_c_homologs, hg, fa_subject, bl_reference){
  
  top_seqs <- getCDSseq(top_homologs, fa_subject)
  
  back_check <- sapply(1:length(top_seqs), function(h){
    cl <- predict(bl_reference, top_seqs[[h]])
    cl <- orderStartEnd(cl, 'S.start', 'S.end')
    query <- GRanges(seqnames = cl$SubjectID, ranges = IRanges(cl$S.start,cl$S.end))  
    subject <- GRanges(seqnames = hg$seqid, ranges = IRanges(hg$start,hg$end))
    
    return(any(query %over% subject))
  })
  
  #this will most likely not be needed as all 
  #the checks should be on the c_homologs table
  
  #back_check_homologs <- rep(FALSE,nrow(top_homologs))
  #for (h in names(top_seqs)){
  #  indexes <- which(top_homologs$group == h)
  #  back_check_homologs[indexes] <- back_check
  #}
  
  #top_homologs$check_back <- back_check_homologs
  
  #this is needed as a precaution just in case 
  #the names get mixed up
  
  indexes <- sapply(top_c_homologs$groups, function(h){
    return(which(names(top_seqs) == h))
  })
  top_c_homologs$back_check <- back_check[indexes]
  return(top_c_homologs)
  
}

########################################################################
