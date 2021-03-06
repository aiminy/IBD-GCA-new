**A Generalized General Combining Ability model for genomic selection**

**Abstract**

In maize breeding, there are several previous studies that breeders pool the A/* and */B biparental crosses to form a training set to estimate the marker effects, and to use these effect estimations to predict the performance of certain trait for A/B cross. This method is referred as the General Combining Ability(GCA) model. The pooling method in GCA model is based on identifying the crosses that one of their parents is shared by A/B crosses. It is obvious that GCA model will not work for the cases that there are no crosses sharing parents with A/B cross. In this paper, we use Genotype-By-Sequencing data to develop a Generalized GCA(GGCA) model to address the shortcoming of the GCA model by identifying Identical-By-Descent(IBD) regions on the genomes between A/B and other parents. We show GCA model is a special case of GGCA model. The application of GGCA model to several biparental crosses demonstrates  GGCA model has its special advantage over GCA model. 

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

In the IBD-GCA model, we used BEAGLE to identify the IBD segments between parents of target populations and that of training populations, and assigned weight 1 for IBD segments and weight 0 to NON-IBD segments. Considering the possibilities that some IBD segments could not be identified by BEAGLE, we relaxed weight assignment for NON-IBD segments to explore the effect of weight relaxation to NON-IBD segmenst on predction accuracies. To perform this weight relaxation, we fixed weight assignment for IBD segments to 1, and assigned a weight value for NON-IBD segments by sampling from uniform distribution[0,1] 10 times. The following Figure shows the effects of relaxing weight assignment to NON-IBD segments on prediction accuracires for different crosses.

![Image of Yaktocat](/Relax_NON_IBD_based_on_mapping_depth_2.png)
The above Figure shows three situations related to the effects of relaxing weight assignments to NON-IBD segments on prediction accuracies for these 19 crosses: 
* No effects:1008,1015,1018,1021,1115,1116,1117,1118. 
  * For these 8 crosses, increasing weights to NON-IBD segments does not have siginificant effects on prediction accuracuies, which indicates IBD segments we used have included the information for prediction.  
* Prediction accuracies decrease as the weights to NON-IBD segments increase: 1016,1019,1028,1114,1120
  * We observed the decreased prediction accuracies as the weights to NON-IBD segments increase for these 5 crosses. This situation could be due that borrowing information from NON-IBD segments could bring some noises, and these noises  could decrease prediction power.
* Prediction accuracies increase as the weight to NON-IBD segments increase:1017,1020,1023,1119,1121,1122
  * These 6 crosses shows the increased prediction accuracies as the weight to NON-IBD segments increase. This means by borrowing some information in NON-IBD segemnsts that BEAGLE failed to identify we retain some prediction power.

So far, we assigned weights based on whether the segment is in IBD or not only. We observed most of test crosses share at least one parent with other training crosses, so we combined the pedigree and IBD relationship together,and came up with a combined weight assignment shema as the following:

![Image of Yaktocat](/Ped_Plus_IBD.png)

In this combined weight assignment shema, we assigned more weight to the SNPs for the test crosses that have at least one parent shared with the training cross. For example, in above Figure we gave weight 2 to the SNPs that their effects are estimated from 1015,1016,1017,1018 crosses since these 4 crosses share one parent with test cross 1008, and they are also in IBD. For 1019, we gave weight 1 to the SNPs that their effect are estimated from 1019 cross since 1019 does not share partent with 1008, but these SNPs are in the IBD segments between parents of 1019 and 1008. Likewisely, we  performed weight relaxation for NON-IBD segments by sampling from uniform distribution[0,1] 10 times for this combined weight assignment schema. The following Figure shows the effects of relaxing weight assignment to NON-IBD segments on prediction accuracires for different crosses when using this combined weight assignment schema.

![Image of Yaktocat](/Results_Ped_Plus_IBD.png)

The above Figure shows we got the similar results as we performed weight relaxtion based on IBD relationship only.

These results demonstrate that IBD-GCA model can be seen as a filter that includes the related information between parents of tartget population and that of training population and excludes the unrelated noise between them. The suitable balancing between information inclusion and noise exclusion can help to obtain prediction power for the genomic selection of target population.
