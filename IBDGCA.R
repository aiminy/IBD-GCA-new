#Usage: start R, source("IBDGCA.R") 
#
#This program needs 5 input files. What each file included is described in readme file

load("Population_Parents_Data.RData")
#load("New_IBD_data.RData")
#load("New_IBD_data_2.RData")

load("IBD_All_Snp.RData")
#load("Snp_Index_After_filter.RData")

load("Snp_Index_4_all_Snp.RData")
#load("~/Marker_effect_estimation/CombinedMarkerEffectEstimation.RData")
load("CombinedMarkerEffectEstimation.RData")
load("Pop_19_snp_plant_heigth_3_25_2015.RData")


FindOtherPopulationIbdwithThisPopulation<-function(population.parents.name.4,wema.ibd.data.2,this.pop)
{

cat(this.pop,"\n")

maped.parents.4.p1<-intersect(unique(as.character(population.parents.name.4[which(as.character(population.parents.name.4[,1])==this.pop&as.character(population.parents.name.4[,3])=="p1"),2])),unique(as.character(wema.ibd.data.2[,1])))

maped.parents.4.p2<-intersect(unique(as.character(population.parents.name.4[which(as.character(population.parents.name.4[,1])==this.pop&as.character(population.parents.name.4[,3])=="p2"),2])),unique(as.character(wema.ibd.data.2[,1])))

maped.parents<-c(maped.parents.4.p1,maped.parents.4.p2)


wema.mapped.index<-which(as.character(wema.ibd.data.2[,1]) %in% maped.parents)

population.parents.name.5<-population.parents.name.4[which(as.character(population.parents.name.4[,1])!=this.pop),]

mapped.parents.2<-intersect(unique(as.character(wema.ibd.data.2[wema.mapped.index,5])),unique(as.character(population.parents.name.5[,2])))

mapped.population<-unique(as.character(population.parents.name.5[which(as.character(population.parents.name.5[,2]) %in% mapped.parents.2),1]))

all.pop.pair.start.end<-data.frame()

for(j in 1:length(mapped.population)){

sub.maped.parents<-intersect(unique(as.character(population.parents.name.5[which(as.character(population.parents.name.5[,1])==mapped.population[j]),2])),unique(as.character(wema.ibd.data.2[wema.mapped.index,5])))

sub.wema.mapped.index<-which(as.character(wema.ibd.data.2[wema.mapped.index,5]) %in% sub.maped.parents)


pop.pair.start.end<-cbind(rep(paste(this.pop,"at_c1",mapped.population[j],sep="_"),dim(wema.ibd.data.2[wema.mapped.index,][sub.wema.mapped.index,9:10])[1]),wema.ibd.data.2[wema.mapped.index,][sub.wema.mapped.index,9:10],rep(this.pop,dim(wema.ibd.data.2[wema.mapped.index,][sub.wema.mapped.index,9:10])[1]),rep(mapped.population[j],dim(wema.ibd.data.2[wema.mapped.index,][sub.wema.mapped.index,9:10])[1]))

colnames(pop.pair.start.end)<-c("pop.pair","start","end","target_pop","matched_pop")

all.pop.pair.start.end<-rbind(all.pop.pair.start.end,pop.pair.start.end)


}

return(all.pop.pair.start.end)

}


CallFindOtherPopulationIbdwithThisPopulation<-function(populatio.parents.name.4,wema.ibd.data.2){

pop.name<-unique(as.character(population.parents.name.4[,1]))

all.mapped.pop<-data.frame()


for(i in 1:length(pop.name)){

mapped.pop<-FindOtherPopulationIbdwithThisPopulation(population.parents.name.4,wema.ibd.data.2,pop.name[i])

all.mapped.pop<-rbind(all.mapped.pop,mapped.pop)


}

return(all.mapped.pop)

}


RefineIBDsegments<-function(pop.map){

all.pop<-unique(c(unique(as.character(pop.map[,4])),unique(as.character(pop.map[,5]))))

all.sorted.pop.map.c5.c6.binding<-list()

for(i in 1:length(all.pop)){

pop.map.c5<-pop.map[which(pop.map[,4]==all.pop[i]),]

pop.map.c6<-pop.map[which(pop.map[,5]==all.pop[i]),]
pop.map.c6.reformated<-cbind(pop.map.c6[,1:3],pop.map.c6[,5],pop.map.c6[,4])
colnames(pop.map.c6.reformated)[4:5]=c("target_pop","matched_pop")

pop.map.c5.c6.binding<-rbind(pop.map.c5,pop.map.c6.reformated)

sorted.pop.map.c5.c6.binding<-pop.map.c5.c6.binding[order(pop.map.c5.c6.binding[,5]),]

all.sorted.pop.map.c5.c6.binding[[i]]<-sorted.pop.map.c5.c6.binding
}

return(all.sorted.pop.map.c5.c6.binding)

}

Re<-CallFindOtherPopulationIbdwithThisPopulation(populatio.parents.name.4,data.ibd)
Re.refined<-RefineIBDsegments(Re)

Map2OneIndex<-function(i,sorted_ca_1021){

target.pop<-unique(as.character(sorted_ca_1021[,4]))

mapped.segment<-sorted_ca_1021[which(sorted_ca_1021[,2]<=i&sorted_ca_1021[,3]>=i),]

num=dim(mapped.segment)[1]

if(num!=0){

mapped.pop<-mapped.segment[,5]

pos.count.target.mapped<-list(pos=i,maping_count=length(unique(as.character(mapped.pop))),target_pop=unique(as.character(target.pop)),mapped_pop=unique(as.character(mapped.pop)))
}else{

pos.count.target.mapped<-list(pos=i,mapping_count=0,target_pop=unique(as.character(target.pop)),mapped_pop="No_mapping")

}

return(pos.count.target.mapped)

}

Map2OneIndex4OnePop<-function(ca.1021){

sorted.ca.1021<-ca.1021[order(ca.1021[,2]),]

genome.index<-seq(min(sorted.ca.1021[,2]),max(sorted.ca.1021[,3]))
#sorted_ca_1021<-sorted.ca.1021

Re.genome.index<-t(sapply(genome.index,Map2OneIndex,sorted.ca.1021))

return(Re.genome.index)

}

Re.refined.mapping.new<-lapply(Re.refined,Map2OneIndex4OnePop)

IBDGCAMarkerEffectEstimation4One<-function(i,Index_mapping,pop_coverage,combined_marker_estimation){

dim(combined_marker_estimation)

matched_pos<-as.numeric(Index_mapping[i,2])
matched_pop<-pop_coverage[i,4]$mapped_pop

target.pop.marker.index<-matched_pos

mapped.pop.index<-which(names(combined_marker_estimation[matched_pos,]) %in% matched_pop)

target.pop.marker.estimation<-mean(combined_marker_estimation[matched_pos,mapped.pop.index])


target.pop.marker.index.estimation<-cbind(target.pop.marker.index,target.pop.marker.estimation)

return(target.pop.marker.index.estimation)

}

Index.mapping<-cbind(seq(0,(length(snp.index.after.filter)-1)),snp.index.after.filter)

IBDGCAMarkerEffectEstimation4OnePop<-function(Re_refined_mapping_new_1){

Index.mapping.2<-Index.mapping[which(Index.mapping[,1] %in% Re_refined_mapping_new_1[,1]),]

Re_refined_mapping_new_1_estimation<-t(sapply(seq(1,length(Index.mapping.2[,1])),IBDGCAMarkerEffectEstimation4One,Index.mapping.2,Re_refined_mapping_new_1,combined.marker.estimation))

return(Re_refined_mapping_new_1_estimation)

}

#Re.refined.mapping.new.1.estimation<-IBDGCAMarkerEffectEstimation4OnePop(Re.refined.mapping.new[[1]])
Re.refined.mapping.new.estimation<-lapply(Re.refined.mapping.new,IBDGCAMarkerEffectEstimation4OnePop)


#load("Pop_19_snp_plant_heigth_3_25_2015.RData")

ImputatePhenotype<-function(Phenotype.Y.2){

NA.index=which(is.na(Phenotype.Y.2))

if(length(NA.index)!=0){
imputation.value<-mean(Phenotype.Y.2[-which(is.na(Phenotype.Y.2))])
Phenotype.Y.2[NA.index]<-imputation.value
}

return(Phenotype.Y.2)

}


IBDGCAPrediction4OnePop<-function(i){

target.pop<-unique(Re.refined.mapping.new[[i]][,3])

data.snp<-pop.acc.height.snp.16[,4:dim(pop.acc.height.snp.16)[2]]

pop.index<-which(pop.acc.height.snp.16[,1] %in% target.pop)
snp.index<-which(as.integer(names(data.snp)) %in% Re.refined.mapping.new.estimation[[i]][,1])

data.snp.after.qc<-data.snp[pop.index,snp.index]

Phenotype.Y<-ImputatePhenotype(pop.acc.height.snp.16[which(pop.acc.height.snp.16[,1] %in% target.pop),3])

Re.refined.mapping.new.estimation[[i]][which(is.na(Re.refined.mapping.new.estimation[[i]][,2])),2]<-0


target.pop.ibd.gca.trait.prediction<-as.matrix(data.snp.after.qc)%*%as.matrix(Re.refined.mapping.new.estimation[[i]][,2])

#cat(cor(target.pop.ibd.gca.trait.prediction,Phenotype.Y),"\n")
return(cbind(target.pop,cor(target.pop.ibd.gca.trait.prediction,Phenotype.Y)))

}

prediction.acc.each.pop<-t(sapply(1:19,IBDGCAPrediction4OnePop))

