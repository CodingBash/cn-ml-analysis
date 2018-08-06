#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

library(GenomicRanges)
library(ggplot2) 
library(reshape) 
library(reshape2)
library(made4)
library(cluster)
library(spatstat) # "im" function 

retrieveCores <- function(dir){
  return(read.table(dir, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
}

retrieveOrganoidSlices <- function(dir){
  return(read.table(dir, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
}

retrieveGeneListTrainingSet <- function(loaded_samples, gene_list, impute = FALSE, reference, res_dir){
  matrix_training_set <- data.frame(stringsAsFactors = FALSE)
  for(sample in loaded_samples){
    print(sample)
    facets_segment_data <- retrieveFacetsSegments(sample, sample_subdir = "/", reference = reference, dir = res_dir)
    facets_segment_data <- segmentsToBedFormat(facets_segment_data, median.clust = FALSE)
    gene_score <- data.frame(stringsAsFactors = FALSE)
    for(gene.i in seq(nrow(gene_list))){
      gene <- gene_list[gene.i, ]
      segments <- facets_segment_data[facets_segment_data$chrom == gene$chrom & ((facets_segment_data$start <= gene$start & gene$start <= facets_segment_data$end) | (facets_segment_data$start <= gene$end & gene$end <= facets_segment_data$end)), ]
      # TODO: genes may intersect across two different segments - will create logic to determine value
      value <- ifelse(nrow(segments) > 0, segments$value, NA)
      gene_score <- rbind(gene_score, cbind(gene, value))
    }
    gene_score <- gene_score[,c(1,2,3,5,4)]
    colnames(gene_score)[4] <- "value"
    feature.entry <- t(data.frame(gene_score[,"value"], row.names = gene_score$name))
    row.names(feature.entry) <- sample
    matrix_training_set <- rbind(matrix_training_set, feature.entry)  
  }
  
  # TODO: Impute NAs with column mean
  for(i in 1:ncol(matrix_training_set)){
    matrix_training_set[is.na(matrix_training_set[,i]), i] <- median(matrix_training_set[,i], na.rm = TRUE)
  }
  
  
  melted_training_set <- do.call(rbind, lapply(seq(nrow(matrix_training_set)), function(index){
    return(do.call(rbind, lapply(colnames(matrix_training_set[index, ]), function(featureId, index){
      featureEntry <- data.frame(score = matrix_training_set[index, featureId], featureId = featureId, sampleId = rownames(matrix_training_set[index, ])[1])
      return(featureEntry)
    }, index)))
  }))
  
  
  return(list(melted=melted_training_set, matrix=matrix_training_set))
}

retrieveSlicingTrainingSet <- function(loaded_samples, Acores, Dcores, ADcores, organoidSlicesFile){
  melted_training_set <- data.frame(stringsAsFactors = FALSE)  
  matrix_training_set <- data.frame(stringsAsFactors = FALSE)
  
  organoidSlices <- retrieveOrganoidSlices(organoidSlicesFile)
  for(sample in loaded_samples){
    #
    # Retrieve necessary data for feature value calculation (bins with cnlr and COREs)
    #
    final_incidence_table <- data.frame(sample= sample)
    if(!missing(ADcores)){
      print("Adding ADcores")
      incidence.input.core <- ADcores[,c(6,7)] # TODO: CORE should have headers so that I can call based on colname instead of index
      colnames(incidence.input.core) <- c("start", "end")
      incidence.input.events <- organoidSlices[organoidSlices$profID == sample, c("gstart", "gend")] 
      colnames(incidence.input.events) <- c("start", "end")
      incidence.output.table <- incidence(incidence.input.core, incidence.input.events, dropevents="Greedy",assoc="I")
      final_incidence_table <- cbind(final_incidence_table, t(as.data.frame(incidence.output.table)))
    } else {
      if(!missing(Acores)){
        print("Adding Acores")
        incidence.input.core <- Acores[,c(6,7)] # TODO: CORE should have headers so that I can call based on colname instead of index
        colnames(incidence.input.core) <- c("start", "end")
        incidence.input.events <- organoidSlices[organoidSlices$profID == sample & organoidSlices$lesion.type == 1, c("gstart", "gend")] 
        colnames(incidence.input.events) <- c("start", "end")
        incidence.output.table <- incidence(incidence.input.core, incidence.input.events, dropevents="Greedy",assoc="I")
        final_incidence_table <- cbind(final_incidence_table, t(as.data.frame(incidence.output.table)))
      }
      if(!missing(Dcores)){
        print("Adding Dcores")
        incidence.input.core <- Dcores[,c(6,7)] # TODO: CORE should have headers so that I can call based on colname instead of index
        colnames(incidence.input.core) <- c("start", "end")
        incidence.input.events <- organoidSlices[organoidSlices$profID == sample & organoidSlices$lesion.type == 0, c("gstart", "gend")] 
        colnames(incidence.input.events) <- c("start", "end")
        incidence.output.table <- incidence(incidence.input.core, incidence.input.events, dropevents="Greedy",assoc="I")
        final_incidence_table <- cbind(final_incidence_table, t(as.data.frame(incidence.output.table)))
      } 
    }
    rownames(final_incidence_table) <- sample
    final_incidence_table <- final_incidence_table[, -c(1)]
    matrix_training_set <- rbind(matrix_training_set, final_incidence_table)
  }
  
  #
  # Convert melt the matrix_training_set
  #
  melted_training_set <- do.call(rbind, lapply(seq(nrow(matrix_training_set)), function(index){
    return(do.call(rbind, lapply(colnames(matrix_training_set[index, ]), function(featureId, index){
      featureEntry <- data.frame(score = matrix_training_set[index, featureId], featureId = featureId, sampleId = rownames(matrix_training_set[index, ])[1])
      return(featureEntry)
    }, index)))
  }))
  
  return(list(melted=melted_training_set, matrix=matrix_training_set))
}

retrieveSegmentTrainingSet <- function(loaded_samples, Acores, Dcores, ADcores, sample_subdir = "/analysis/structural_variants/", reference = "NA12878", dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/", removeNAs = TRUE){
  melted_training_set <- data.frame(stringsAsFactors = FALSE)  
  matrix_training_set <- data.frame(stringsAsFactors = FALSE)
  for(sample in samples){
    #
    # Retrieve necessary data for feature value calculation (bins with cnlr and COREs)
    #
    facets_snp_data <- retrieveFacetsSnps(sample, sample_subdir = sample_subdir, reference = reference, dir = dir)
    snp_bed <- snpsToBedFormat(facets_snp_data)
    
    #
    # Preprocess COREs
    #
    gr <- NA
    if(missing(ADcores)){
      seqnames_A <- table(Acores[[1]])
      gr_A <- GRanges(
        seqnames = Rle(names(seqnames_A), as.vector(seqnames_A)),
        ranges = IRanges(Acores[[2]], Acores[[3]], names = row.names(Acores)),
        event=rep("A", nrow(Acores))
      )
      seqnames_D <- table(Dcores[[1]])
      gr_D <- GRanges(
        seqnames = Rle(names(seqnames_D), as.vector(seqnames_D)),
        ranges = IRanges(Dcores[[2]], Dcores[[3]], names = row.names(Dcores)),
        event=rep("D", nrow(Dcores))
      )
      # TODO: Do operations on GRanges object to simplify, if need be (i.e. reduction)
      gr <- c(gr_A, gr_D)
    } else {
      seqnames_AD <- table(ADcores[[1]])
      gr_AD <- GRanges(
        seqnames = Rle(names(seqnames_AD), as.vector(seqnames_AD)),
        ranges = IRanges(ADcores[[2]], ADcores[[3]], names = row.names(ADcores)),
        event=rep("AD", nrow(ADcores))
      )
      gr <- gr_AD
    }
    coreDf <- as.data.frame(gr, row.names = seq_along(gr$event))
    coreDf$score <- NA
    
    #
    # Determine value of each CORE feature (cnlr)
    #
    coreDf[, 7] <- sapply(seq(1, nrow(coreDf)), function(core.index){
      # Get snps with same chromosome as core, and in between the core region
      snps <- snp_bed[snp_bed[[1]] == coreDf[core.index, 1] & snp_bed[[2]] >= coreDf[core.index, 2] & snp_bed[[3]] <= coreDf[core.index, 3], ]
      # Calculate and assign median
      return(median(snps$value))
    } )
    
    # Convert back to Granges object - this is our final feature object
    featureSet <- makeGRangesFromDataFrame(coreDf, keep.extra.columns = TRUE)
    
    sampleTrainingSet <- data.frame(score = coreDf[,7])
    sampleTrainingSet$featureId <- rownames(coreDf)
    sampleTrainingSet$sampleId <- sample
    
    melted_training_set <- rbind(melted_training_set, sampleTrainingSet)
    
    matrix_training_entry <- t(as.data.frame(coreDf[,7]))
    rownames(matrix_training_entry) <- c(sample)
    colnames(matrix_training_entry) <- rownames(coreDf)
    colnames(matrix_training_entry) <- rownames(coreDf)
    matrix_training_set <- rbind(matrix_training_set, matrix_training_entry) 
  }
  if(removeNAs == TRUE){
    melted_training_set <- melted_training_set[-which(is.na(melted_training_set$score)),]
    matrix_training_set <- matrix_training_set[sapply(matrix_training_set, function(x) !any(is.na(x)))] 
  }
  return(list(melted=melted_training_set, matrix=matrix_training_set))
}

#Compute incidence of cores in a set of events using inclusion metric. If
#the events matrix represents events in a single profile, the result is a row in
#the incidence table.
#
#Input: 
#cores and events are both two-column matrices, with columns named "start" and 
#"end".
#dropevents (character) is a switch for the event-core matching rule. If
#dropevents="Drop" (default), best match is found
#among all pairs of events and cores and the incidence of the best matching core
#is set to that value. Both the event and the core are dropped, and the process
#is repeated until there are no cores left. If dropevents is not "Drop" or 
#"Greedy", only the core is dropped and the process is repeated until there 
#are no cores left. If dropevents="Greedy", the process essentially mimicks the
#greedy implementation of the CORE algorithm, with the order of the cores as 
#specified by the input. But the score for a core is maximized, not summed,
#over all events.
#
#Value: a numerical vector of length equal to nrow(cores), containing matching
#values for the cores.
#
incidence<-function(cores,events,dropevents="Drop",assoc="I"){
  fulltruth<-matrix(ncol=nrow(cores),nrow=nrow(events),data=0)
  if(assoc=="I")for(i in 1:nrow(cores))fulltruth[,i]<-
      (events[,"start"]<=cores[i,"start"]&events[,"end"]>=cores[i,"end"])*
      (cores[i,"end"]-cores[i,"start"]+1)/(events[,"end"]-events[,"start"]+1)
  else if(assoc=="J")for(i in 1:nrow(cores))fulltruth[,i]<-
      pmax(((pmin(events[,"end"],cores[i,"end"])-pmax(events[,"start"],cores[i,"start"])+1)/
              (pmax(events[,"end"],cores[i,"end"])-pmin(events[,"start"],cores[i,"start"])+1)),0)
  else stop("Undefined association")
  coretruth<-rep(0,nrow(cores))
  if(nrow(events)==0)return(coretruth)
  if(dropevents!="Greedy"){
    corelabel<-1:nrow(cores)
    for(i in 1:nrow(cores)){
      bestmatch<-which.max(fulltruth)
      whichcore<-(bestmatch-1)%/%nrow(fulltruth)+1
      whichevent<-(bestmatch-1)%%nrow(fulltruth)+1
      coretruth[corelabel[whichcore]]<-fulltruth[bestmatch]
      if(dropevents=="Drop")fulltruth<-fulltruth[-whichevent,,drop=F]
      fulltruth<-fulltruth[,-whichcore,drop=F]
      corelabel<-corelabel[-whichcore]
    }
  }
  else{ 
    weights<-rep(1,nrow(events))
    for(i in 1:nrow(cores)){
      bestmatch<-which.max(fulltruth[,i]*weights)
      coretruth[i]<-fulltruth[bestmatch,i]*weights[bestmatch]
      weights[bestmatch]<-weights[bestmatch]*(1-fulltruth[bestmatch,i])
    }
  }
  return(coretruth)
}


attachLabelsToSet <- function(matrix_training_set, labelData){
  sampleList <- rownames(matrix_training_set)
  labelLists <- lapply(names(labelData), function(label){
    aucList <- unlist(sapply(sampleList, function(sample, label){
      labelMatrix <- labelData[[label]]  
      auc <- c(labelMatrix[labelMatrix$SampleId == sample, ]$AUC, NA)[1]
      return(auc)
    }, label))
    return(aucList)
  })
  names(labelLists) <- names(labelData)
  labelDataframe <- do.call(cbind.data.frame, labelLists)
  labeled_matrix_training_set <- cbind(labelDataframe, matrix_training_set)
  return(labeled_matrix_training_set)
}