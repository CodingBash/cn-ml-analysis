
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
setwd("~/Documents/Git-Projects/Git-Research-Projects/cn-ml-analysis")
source("helperFunctions.r")
source("featureMatrixAssignment.r")
source("unsupervisedLearning.r")

#
# Load sample to retrieve feature set for
# TODO: There seems to be a scope conflict - samples is getting overwritten
#
samples <- load_samples(classes = c("T","F", "M"), sampleList = "./resources/sampleList.csv")

#
# Retrieve CORE features
#
setwd("~/Documents/Git-Projects/Git-Research-Projects/cnprep_cores")
ADcores <- retrieveCores("./hT_output/prev_run_7_27_2018_8_2/selectedCores/ADselectedCoresBP.bed") # BED file of recurrent regions
Acores <- retrieveCores("./hT_output/prev_run_7_27_2018_8_2/selectedCores/AselectedCoresBP.bed") # BED file of recurrent regions
Dcores <- retrieveCores("./hT_output/prev_run_7_27_2018_8_2/selectedCores/DselectedCoresBP.bed") # BED file of recurrent regions


head(ADcores)
head(Acores)
head(Dcores)

setwd("~/Documents/Git-Projects/Git-Research-Projects/cn-ml-analysis")
aucData <- readRDS("./resources/listSampleTESAUC.RDS")

head(aucData$Gemcitabine)
head(aucData$Paclitaxel)
head(aucData$`SN-38`)
head(aucData$`5-FU`)
head(aucData$Oxaliplatin)

#
# Retrieve training set
#
setwd("~/Documents/Git-Projects/Git-Research-Projects/FACETS_write_files")
segment_training_set <- retrieveSegmentTrainingSet(loaded_samples = samples, Acores = Acores, Dcores = Dcores, sample_subdir = "/", reference = "hN30", dir = "output/FACETS_Reference_hN30_8_2_18_1/")
segment_training_set$matrix <- attachLabelsToSet(matrix_training_set = segment_training_set$matrix, labelData = aucData)

reference <- "hN30"
res_dir <- "output/FACETS_Reference_hN30_8_2_18_1/" # Determine FACETS reference to use _ hN30

setwd("~/Documents/Git-Projects/Git-Research-Projects/gene-set-feature-scoring/")
gene_list <- read.table("./resources/genes_hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gene_list <- gene_list[, c(7,8,9,3)]
colnames(gene_list) <- c("chrom", "start", "end", "name")

setwd("~/Documents/Git-Projects/Git-Research-Projects/FACETS_write_files/")
gene_training_set <- retrieveGeneListTrainingSet(samples, gene_list, impute = TRUE, reference, res_dir)
gene_training_set$matrix <- attachLabelsToSet(matrix_training_set = gene_training_set$matrix, labelData = aucData)

setwd("~/Documents/Git-Projects/Git-Research-Projects/CNprep-Slicing-CORE-Analysis/")
Acores <- retrieveCores("./output/coresResults/prev_run_8_2_2018_2/selectedCores/AselectedCoresBP.bed") # BED file of amplification recurrent regions
Dcores <- retrieveCores("./output/coresResults/prev_run_8_2_2018_2/selectedCores/DselectedCoresBP.bed") # BED file of deletion recurrent regions
ADcores <- retrieveCores("./output/coresResults/prev_run_8_2_2018_2/selectedCores/ADselectedCoresBP.bed") # BED file of both recurrent regions

slicing_training_set <- retrieveSlicingTrainingSet(samples, Acores, Dcores, organoidSlicesFile = "./resources/slicingOutput/table/prev_run_8_2_2018_3/organoidSlices.txt")
slicing_training_set$matrix <- attachLabelsToSet(matrix_training_set = slicing_training_set$matrix, labelData = aucData)

options(repr.plot.width=15, repr.plot.height=15)
visualizeUnclusteredHeatmap(segment_training_set$melted)
visualizeUnclusteredHeatmap(slicing_training_set$melted)
visualizeUnclusteredHeatmap(gene_training_set$melted)

head(segment_training_set$melted)

options(repr.plot.width=15, repr.plot.height=15)
hc_segment <- clusterTrainingSet(segment_training_set$melted, visualize = TRUE)
hc_slicing <- clusterTrainingSet(slicing_training_set$melted, visualize = TRUE)
hc_gene <- clusterTrainingSet(gene_training_set$melted, visualize = TRUE)

options(repr.plot.width=15, repr.plot.height=7)
plot(hc_segment)
plot(hc_slicing)
plot(hc_gene)

options(repr.plot.width=15, repr.plot.height=15)


for(label.i in seq(6)){
    labelData <- aucData[[label.i]]
    test_set <- cbind(training_set$melted)
    
    test_set$sampleId <- sapply(test_set$sampleId, 
                                function(sample){
                                    return(paste0(sample, "-", label.i, "-", labelData[labelData$SampleId == sample, ]$AUC))
                                })
    clusterTrainingSet(test_set, visualize = TRUE)
}



#setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction/mlOutput")
#write.csv(training_set$matrix, file ="coreTrainingSet_7_31_2018_1.csv")
