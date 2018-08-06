#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

library(GenomicRanges)
library(ggplot2) 
library(reshape) 
library(reshape2)
library(made4)
library(cluster)
library(spatstat) # "im" function 


visualizeUnclusteredHeatmap <- function(melted_training_set){
  ggplot(data = melted_training_set, aes(x = featureId, y = sampleId)) + 
    geom_tile(aes(fill = score), color = "white", size = 1) + 
    scale_fill_gradient2(low = "blue", mid="white", high = "tomato") + 
    xlab("core ID") + 
    theme_grey(base_size = 10) + 
    ggtitle("Heatmap (ggplot)") + 
    theme(axis.ticks = element_blank(), 
          panel.background = element_blank(), 
          plot.title = element_text(size = 12, colour = "gray50")) 
}

clusterTrainingSet <- function(melted_training_set, visualize = FALSE){
  # Unmelt training set for correlation analysis
  training_set_matrix <- dcast(data = melted_training_set,formula = sampleId~featureId,fun.aggregate = sum,value.var = "score")
  sampleIds <- training_set_matrix$sampleId
  training_set_matrix <- training_set_matrix[,-c(1)]
  training_set_matrix <- t(training_set_matrix)
  colnames(training_set_matrix) <- sampleIds
  training_set_matrix <- as.data.frame(training_set_matrix)
  
  # Remove samples with 0 variance
  nonzero_variance_samples <- unlist(lapply(colnames(training_set_matrix), function(sample){
    if(var(na.omit(training_set_matrix[, sample])) != 0){
      return(sample)
    }
  }))
  training_set_matrix <- training_set_matrix[,c(nonzero_variance_samples)]
  
  # Calculate distance matrix
  corRaw <- cor(training_set_matrix)
  dissimilarity <- 1 - corRaw
  distance.sample <- as.dist(dissimilarity)
  distance.core <- as.dist(t(dissimilarity))
  
  # Run hierarchical clustering
  hc.sample <- hclust(distance.sample)
  hc.core <- hclust(distance.core)
  
  if(visualize == TRUE){
    color.palette  <- colorRampPalette(c("blue4", "white", "red4"))(n=600)
    #    heatmap.2(t(training_set_matrix),  trace="none", dendrogram="row", Rowv = as.dendrogram(hc.sample), Colv = as.dendrogram(hc.core), density.info = 'none', scale='none', col = color.palette)
    heatmap.2(t(training_set_matrix), trace = "none", density="none", scale = "row", hclust=function(x) hclust(x, method="complete"), distfun=function(x) as.dist((1-cor(t(x)))/2), col = color.palette)
  }
  
  return(hc.sample)  
}