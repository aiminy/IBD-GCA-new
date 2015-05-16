# IBD-GCA-new
This is a much improved IBD-GCA model for genomic prediction:
To run this program, 5 input files are needed:

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


After you have these 5 files in your current directory, then you can start R, and type source("IBDGCA.R"). 
You will find the results in "prediction.acc.each.pop"

If you want to have these 5 files to have a try, please email me at ay247@cornell.edu 


