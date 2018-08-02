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

cd_local("resources")
aucData <- readRDS("listSampleTESAUC.RDS")


#
# Retrieve training set
#
setwd("~/Git-Projects/Git-Research-Projects/FACETS_write_files")
training_set <- retrieveTrainingSet(loaded_samples = samples, ADcores = ADcores, sample_subdir = "/", reference = "hN31", dir = "output/FACETS_Reference_hN31_7_28_18_2/")
training_set$matrix <- attachLabelsToSet(matrix_training_set = training_set$matrix, labelData = aucData)

visualizeUnclusteredHeatmap(training_set$melted)
hc <- clusterTrainingSet(training_set$melted, visualize = TRUE)
plot(hc)

cd_local("mlOutput")
write.csv(training_set$matrix, file ="coreTrainingSet_7_31_2018_1.csv")
