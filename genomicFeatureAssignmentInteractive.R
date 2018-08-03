#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

setwd("~/Git-Projects/Git-Research-Projects/cn-ml-analysis")
source("genomicFeatureAssignment.R")

#
# Load sample to retrieve feature set for
# TODO: There seems to be a scope conflict - samples is getting overwritten
#
samples <- load_samples(classes = c("T","F", "M"), sampleList = "./resources/sampleList.csv")

#
# Retrieve CORE features
#
setwd("~/Git-Projects/Git-Research-Projects/cnprep_cores")
ADcores <- retrieveCores("./hT_output/prev_run_7_27_2018_8_2/selectedCores/ADselectedCoresBP.bed") # BED file of recurrent regions
Acores <- retrieveCores("./hT_output/prev_run_7_27_2018_8_2/selectedCores/AselectedCoresBP.bed") # BED file of recurrent regions
Dcores <- retrieveCores("./hT_output/prev_run_7_27_2018_8_2/selectedCores/DselectedCoresBP.bed") # BED file of recurrent regions

setwd("~/Git-Projects/Git-Research-Projects/cn-ml-analysis")
aucData <- readRDS("./resources/listSampleTESAUC.RDS")


#
# Retrieve training set
#
setwd("~/Git-Projects/Git-Research-Projects/FACETS_write_files")
#training_set <- retrieveTrainingSet(loaded_samples = samples, ADcores = ADcores, sample_subdir = "/", reference = "hN31", dir = "output/FACETS_Reference_hN31_7_28_18_2/")
training_set <- retrieveTrainingSet(loaded_samples = samples, Acores = Acores, Dcores = Dcores, sample_subdir = "/", reference = "hN31", dir = "output/FACETS_Reference_hN31_7_28_18_2/")
training_set$matrix <- attachLabelsToSet(matrix_training_set = training_set$matrix, labelData = aucData)

visualizeUnclusteredHeatmap(training_set$melted)
hc <- clusterTrainingSet(training_set$melted, visualize = TRUE)
plot(hc)

setwd("~/Git-Projects/Git-Research-Projects/cn-ml-analysis")
write.csv(training_set$matrix, file ="mlOutput/coreTrainingSet_8_2_2018_2.csv")
