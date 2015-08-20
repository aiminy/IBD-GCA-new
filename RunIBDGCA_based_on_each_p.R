#Usage: R CMD BATCH --vanilla --slave RunIBDGCA_based_on_each_p.R& 
sink("CalculationBasedOnEachParent.txt")
source("IBDGCA_based_on_each_p.R")
sink()
