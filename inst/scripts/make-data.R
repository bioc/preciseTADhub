### =========================================================================
### preciseTADhub is an ExperimentHub package that stores pre-trained random
### forest models that can be leveraged to predict TAD and/or chromatin loop
### boundaries
### -------------------------------------------------------------------------
###

# preciseTADhub contains 84 lists as RDS objects

# Each list contains 2 objects: 1) a train object from \code{caret} with RF
# model information, and 2) a data.frame of variable importance for each feature
# included in the model

# The file names are stuctured as follows: i_j_k_l.rds
# where i denotes the chromosome that was used as a holdout (i.e. for
#       testing), meaning all other chromosomes were used for training;
#       i = {CHR1, CHR2, ..., CHR21, CHR22}
#
#       j denotes the cell line;
#       j = {GM12878, K562}
#
#       k denotes the resolution (size of genomic bins);
#       k = {5kb, 10kb}
#
#       l denotes the TAD/loop caller used to define ground truth boundaries;
#       l = {Arrowhead, Peakachu}

# For example the file named "CHR1_GM12878_5kb_Arrowhead.rds" is a list whose
# first item is a RF model that was built on data for chromosomes 2-22
# (omitting CHR9; see https://doi.org/10.1101/2020.09.03.282186), binned using
# 5 kb bins, ground truth TAD boundaries were identified using the Arrowhead
# TAD caller at 5 kb on GM12878. All models included the same number of
# predictors including CTCF, RAD21, SMC3, and ZNF143. The second item in the
# list is a data.frame with variable importances for CTCF, RAD21, SMC3, and
# ZNF143.

# The pre-trained models set up users to apply them to predict their own
# boundaries on chromosomes that were heldout, per the framework in the
# preciseTAD paper (https://doi.org/10.1101/2020.09.03.282186).

# The following is an example script describing the steps involved in making
# one of the RDS objects stored in preciseTADhub:
# "CHR1_GM12878_5kb_Arrowhead.rds"

# This example follows from the preciseTAD R package (see the detailed vignette
# at https://bioconductor.org/packages/devel/bioc/html/preciseTAD.html)

library(preciseTAD)

# Read in ARROWHEAD-called TADs at 5kb
data(arrowhead_gm12878_5kb)

# Extract unique boundaries
bounds.GR <- extractBoundaries(domains.mat = arrowhead_gm12878_5kb,
                               preprocess = FALSE,
                               CHR = paste0("CHR", c(2:8,10:22)),
                               resolution = 5000)

# Read in GRangesList of 26 TFBS
data(tfbsList)
tfbsList_filt <- tfbsList[names(tfbsList) %in% c("Gm12878-Ctcf-Broad",
                                                 "Gm12878-Rad21-Haib",
                                                 "Gm12878-Smc3-Sydh",
                                                 "Gm12878-Znf143-Sydh")]

# Create the binned data matrix for CHR2-22 (training)
# using 5 kb binning, distance-type predictors from 4 TFBS from
# the GM12878 cell line, and random under-sampling
tadData <- createTADdata(bounds.GR = bounds.GR,
                         resolution = 5000,
                         genomicElements.GR = tfbsList_filt,
                         featureType = "distance",
                         resampling = "rus",
                         trainCHR = paste0("CHR", c(2:8,10:22)),
                         predictCHR = NULL)

# Perform random forest using TADrandomForest by tuning mtry over 10 values
# using 3-fold CV
CHR1_GM12878_5kb_Arrowhead <- TADrandomForest(trainData = tadData[[1]],
                            testData = NULL,
                            tuneParams = list(mtry = 2,
                                              ntree = 500,
                                              nodesize = 1),
                            cvFolds = 3,
                             cvMetric = "Accuracy",
                             verbose = TRUE,
                             model = TRUE,
                             importances = TRUE,
                             impMeasure = "MDA",
                             performances = FALSE)
