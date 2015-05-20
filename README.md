**A Generalized General Combining Ability model for genomic selection**

**Abstract**

In maize breeding, there are several previous studies that breeders pool the A/* and */B biparental crosses to form a training set to estimate the marker effects, and to use these effect estimations to predict the performance of certain trait for A/B cross. This method is referred as the General Combining Ability(GCA) model. The pooling method in GCA model is based on identifying the crosses that one of their parents is shared by A/B crosses. It is obvious that GCA model will not work for the cases that there are no crosses sharing parents with A/B cross. In this paper, we use Genotype-By-Sequencing data to develop a Generalized GCA(GGCA) model to address the shortcoming of the pedigree based GCA model by identifying Identical-By-Descent(IBD) regions on the genomes between A/B and other parents. We show GCA model is a special case of GGCA model. The application of GGCA model to several biparental crosses demonstrates  GGCA model has its special advantage over GCA model. 

**Method**

To run GGCA model, 5 input files are needed:

The file to define population and their parents:
Population_Parents_Data.RData

The file to identify the IBD segments from Beagle and with our modification to include the segments between parent itself: 
IBD_All_Snp.RData

The file to map the index from Beagle to SNP index: 
Snp_Index_4_all_Snp.RData

The file to have the marker effect estimation from rrBlup within each population:
CombinedMarkerEffectEstimation.RData

The file to include the all SNP and plant height data:
Pop_19_snp_plant_heigth_3_25_2015.RData


After you have these 5 files, make a directory named as data in your current directory, then you can start R, and type source("IBDGCA.R"): 
You will find the results in "prediction.acc.each.pop"

If you want to have these 5 files to have a try, please email me at ay247@cornell.edu

**Results**





