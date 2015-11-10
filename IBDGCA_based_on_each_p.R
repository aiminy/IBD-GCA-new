#Usage: start R, source("IBDGCA4p1.R") 
#
#This program needs 5 input files. What each file included is described in readme file

#Load all functions needed for calculation
source("./WeightedIBDGCA_function.R")

#load population,their parents information
#pop.name parents.name  parents.type
#1008 BREAD_2754  p1
#1008 CML505  p1
#1008 IMAS-87 p1
#1008 KE_Maize19  p1
#1008 BREAD_2765  p2
#1008 IMAS-23 p2
load("./data/Population_Parents_Data.RData")

load("./data/IBD_All_Snp.RData")

load("./data/Snp_Index_4_all_Snp.RData")

#There are 19 cross populations from WEMA. we use each population as a training set to estimate marker effects,
#combine all estimation together to form a 955690*19 matrix. Row name of this matrix is SNP index, column name of
#this matrix is population name
load("./data/CombinedMarkerEffectEstimation.RData")

#Load phenotype and SNP data from 19 populations. The 1st column is population name. The 2nd column is accession name.
#The 3rd column is phenotype data. SNP data starts from 4th column. there are 3402 lines, and each line has 
#955690 SNP markers.
#The dimension of this data matrix:3402*955693
#pop.name acc height  1 2 3
#1008 WEMA_6x1008_MARS-WEMA_269939_tester_CML395_CML444 244.25  -9  -9  -9
#1008 WEMA_6x1008_MARS-WEMA_269940_tester_CML395_CML444 235.75  -9  -9  -9
#1008 WEMA_6x1008_MARS-WEMA_269941_tester_CML395_CML444 224.50  -9  -9  -9
load("./data/Pop_19_snp_plant_heigth_3_25_2015.RData")

FindOtherPopulationIbdwithThisPopulation<-function(population.parents.name.4,wema.ibd.data.2,this.pop,t_p,m_p)
{

cat(this.pop,"\t",t_p,"\t",m_p,"\n")
target.population.parent<-unique(as.character(population.parents.name.4[which(as.character(population.parents.name.4[,1])==this.pop&as.character(population.parents.name.4[,3])==as.character(t_p)),2]))


mapped.population.0<-population.parents.name.4[which(as.character(population.parents.name.4[,1])!=this.pop&as.character(population.parents.name.4[,3])==as.character(m_p)),]


mapped.population.00<-unique(as.character(mapped.population.0[which(mapped.population.0[,2] %in% c(target.population.parent)),1]))

cat("Populations_share_at_least_one_parent:",mapped.population.00,"\n")


mapped.population.parent<-unique(as.character(population.parents.name.4[which(as.character(population.parents.name.4[,1])!=this.pop&as.character(population.parents.name.4[,3])==as.character(m_p)),2]))


wema.mapped.index<-which(as.character(wema.ibd.data.2[,1]) %in% target.population.parent)

wema.ibd.data.3<-wema.ibd.data.2[wema.mapped.index,]


mapped.parents.based.on.ibd.4.target.population.parents<-unique(as.character(wema.ibd.data.3[,5]))

mapped.population.based.on.ibd.4.target.population.parents<-unique(as.character(mapped.population.0[which(mapped.population.0[,2] %in% c(mapped.parents.based.on.ibd.4.target.population.parents)),1]))

population.parents.name.6<-mapped.population.0

mapped.parents.2<-intersect(unique(as.character(wema.ibd.data.3[,5])),unique(as.character(population.parents.name.6[,2])))

mapped.population<-mapped.population.based.on.ibd.4.target.population.parents

cat("Mapped_Populations_Based_On_IBD:",mapped.population,"\n")

all.pop.pair.start.end<-data.frame()

for(j in 1:length(mapped.population)){

sub.maped.parents<-intersect(unique(as.character(population.parents.name.6[which(as.character(population.parents.name.6[,1])==mapped.population[j]),2])),unique(as.character(wema.ibd.data.3[,5])))

sub.wema.mapped.index<-which(as.character(wema.ibd.data.3[,5]) %in% sub.maped.parents)

n_mapping=dim(wema.ibd.data.3[sub.wema.mapped.index,9:10])[1]

pop.pair.start.end<-cbind(rep(paste(this.pop,"at_c1",mapped.population[j],sep="_"),n_mapping),wema.ibd.data.3[sub.wema.mapped.index,9:10],rep(this.pop,n_mapping),rep(mapped.population[j],n_mapping))

colnames(pop.pair.start.end)<-c("pop.pair","start","end","target_pop","matched_pop")

all.pop.pair.start.end<-rbind(all.pop.pair.start.end,pop.pair.start.end)

}

return(all.pop.pair.start.end)

}


CallFindOtherPopulationIbdwithThisPopulation<-function(populatio.parents.name.4,wema.ibd.data.2,t_p,m_p){
pop.name<-unique(as.character(population.parents.name.4[,1]))
all.mapped.pop<-data.frame()
for(i in 1:length(pop.name))
{
mapped.pop<-FindOtherPopulationIbdwithThisPopulation(population.parents.name.4,wema.ibd.data.2,pop.name[i],t_p,m_p)
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

Re.p1.p1<-CallFindOtherPopulationIbdwithThisPopulation(population.parents.name.4,data.ibd,"p1","p1")
Re.p1.p2<-CallFindOtherPopulationIbdwithThisPopulation(population.parents.name.4,data.ibd,"p1","p2")
Re.p2.p1<-CallFindOtherPopulationIbdwithThisPopulation(population.parents.name.4,data.ibd,"p2","p1")
Re.p2.p2<-CallFindOtherPopulationIbdwithThisPopulation(population.parents.name.4,data.ibd,"p2","p2")

Re.refined.p1.p1<-RefineIBDsegments(Re.p1.p1)
Re.refined.p1.p2<-RefineIBDsegments(Re.p1.p2)
Re.refined.p2.p1<-RefineIBDsegments(Re.p2.p1)
Re.refined.p2.p2<-RefineIBDsegments(Re.p2.p2)


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

SetStartAndEndIndex<-function(p1p1,p1p2,p2p1,p2p2){
all_pp<-c(p1p1[,2],p1p1[,3],p1p2[,2],p1p2[,3],p2p1[,2],p2p1[,3],p2p2[,2],p2p2[,3])
start_index=min(all_pp)
end_index=max(all_pp)
start_end_index<-list(s=start_index,e=end_index)
return(start_end_index)
}


CallSetStartAndEndIndex<-function(populatio.parents.name.4){
pop.name<-unique(as.character(populatio.parents.name.4[,1]))
Re.refined.mapping.new.p1.p1<-list()
Re.refined.mapping.new.p1.p2<-list()
Re.refined.mapping.new.p2.p1<-list()
Re.refined.mapping.new.p2.p2<-list()
j=1;
for(i in 1:length(pop.name)){
p11_index<-seq_along(Re.refined.p1.p1)[sapply(Re.refined.p1.p1, FUN=function(X) pop.name[i] %in% X[,4])]
p12_index<-seq_along(Re.refined.p1.p2)[sapply(Re.refined.p1.p2, FUN=function(X) pop.name[i] %in% X[,4])]
p21_index<-seq_along(Re.refined.p2.p1)[sapply(Re.refined.p2.p1, FUN=function(X) pop.name[i] %in% X[,4])]
p22_index<-seq_along(Re.refined.p2.p2)[sapply(Re.refined.p2.p2, FUN=function(X) pop.name[i] %in% X[,4])]
if(length(p11_index)!=0)
{
cat(pop.name[i],"\t",p11_index,"\t",p12_index,"\t",p21_index,"\t",p22_index,"\n")
se<-SetStartAndEndIndex(Re.refined.p1.p1[[p11_index]],Re.refined.p1.p2[[p12_index]],Re.refined.p2.p1[[p21_index]],Re.refined.p2.p2[[p22_index]])
ss=se$s
ee=se$e
Re.refined.mapping.new.p1.p1[[j]]<-Map2OneIndex4OnePopWholeIndex(Re.refined.p1.p1[[p11_index]],ss,ee)
Re.refined.mapping.new.p1.p2[[j]]<-Map2OneIndex4OnePopWholeIndex(Re.refined.p1.p2[[p12_index]],ss,ee)
Re.refined.mapping.new.p2.p1[[j]]<-Map2OneIndex4OnePopWholeIndex(Re.refined.p2.p1[[p21_index]],ss,ee)
Re.refined.mapping.new.p2.p2[[j]]<-Map2OneIndex4OnePopWholeIndex(Re.refined.p2.p2[[p22_index]],ss,ee)
j=j+1
}
}
Re.refined.mapping.new<-list
Re.refined.mapping.new<-list(p11=Re.refined.mapping.new.p1.p1,p12=Re.refined.mapping.new.p1.p2,p21=Re.refined.mapping.new.p2.p1,p22=Re.refined.mapping.new.p2.p2)
return(Re.refined.mapping.new)
}

Map2OneIndex4OnePopWholeIndex<-function(ca.1021,sss,eee){
sorted.ca.1021<-ca.1021[order(ca.1021[,2]),]
genome.index<-seq(sss,eee)
Re.genome.index<-t(sapply(genome.index,Map2OneIndex,sorted.ca.1021))
return(Re.genome.index)
}

Re.refined.mapping.new.with.parent<-CallSetStartAndEndIndex(population.parents.name.4)

mapping.matrix.p1.p1<-lapply(Re.refined.mapping.new.with.parent[[1]],ConvertMappingData2Matrix,combined.marker.estimation,as.character(1),as.character(0))
mapping.matrix.p1.p2<-lapply(Re.refined.mapping.new.with.parent[[2]],ConvertMappingData2Matrix,combined.marker.estimation,as.character(1),as.character(0))
mapping.matrix.p2.p1<-lapply(Re.refined.mapping.new.with.parent[[3]],ConvertMappingData2Matrix,combined.marker.estimation,as.character(1),as.character(0))
mapping.matrix.p2.p2<-lapply(Re.refined.mapping.new.with.parent[[4]],ConvertMappingData2Matrix,combined.marker.estimation,as.character(1),as.character(0))

j=1
mapping.matrix.ibd.all.parent.pair<-list()
for(i in 1:19){
mapping.matrix.ibd.all.parent.pair[[j]]<-mapping.matrix.p1.p1[[i]]+mapping.matrix.p1.p2[[i]]+mapping.matrix.p2.p1[[i]]+mapping.matrix.p2.p2[[i]]
j=j+1
}

matrix.based.on.pedigree.p1.p1<-lapply(Re.refined.mapping.new.with.parent[[1]],GenerateLabelMatrixBasedOnPedigreeOneParent,population.parents.name.4,combined.marker.estimation,"p1","p1")
matrix.based.on.pedigree.p1.p2<-lapply(Re.refined.mapping.new.with.parent[[2]],GenerateLabelMatrixBasedOnPedigreeOneParent,population.parents.name.4,combined.marker.estimation,"p1","p2")
matrix.based.on.pedigree.p2.p1<-lapply(Re.refined.mapping.new.with.parent[[3]],GenerateLabelMatrixBasedOnPedigreeOneParent,population.parents.name.4,combined.marker.estimation,"p2","p1")
matrix.based.on.pedigree.p2.p2<-lapply(Re.refined.mapping.new.with.parent[[4]],GenerateLabelMatrixBasedOnPedigreeOneParent,population.parents.name.4,combined.marker.estimation,"p2","p2")

j<-1
matrix.based.on.pedigree<-list()
for (i in 1:19){
cat(i,dim(matrix.based.on.pedigree.p1.p1[[i]]),"\n")
cat(i,dim(matrix.based.on.pedigree.p1.p2[[i]]),"\n")
cat(i,dim(matrix.based.on.pedigree.p2.p1[[i]]),"\n")
cat(i,dim(matrix.based.on.pedigree.p2.p2[[i]]),"\n")
matrix.based.on.pedigree[[j]]<-matrix.based.on.pedigree.p1.p1[[i]]+matrix.based.on.pedigree.p1.p2[[i]]+matrix.based.on.pedigree.p2.p1[[i]]+matrix.based.on.pedigree.p2.p2[[i]]
j=j+1
}

j<-1
matrix.based.on.pedigree.ibd<-list()
for(i in 1:19){
matrix.based.on.pedigree.ibd[[j]]<-matrix.based.on.pedigree[[i]]+mapping.matrix.ibd.all.parent.pair[[i]]
j=j+1
}

j<-1
matrix.based.on.pedigree.ibd.with.snp.index<-list()
for(i in 1:19){
matrix.based.on.pedigree.ibd.with.snp.index[[j]]<-cbind(Re.refined.mapping.new.with.parent[[1]][[i]][,c(1,3)],matrix.based.on.pedigree.ibd[[i]])
j=j+1
}

WeightIBDGCA5<-function(LL_matrix_pedigree_ibd,Estimation_Data,Pop_Snp_Data,w_ibd,w_no_ibd){
LL.matrix.pedigree.ibd<-matrix(unlist(LL_matrix_pedigree_ibd[,3:dim(LL_matrix_pedigree_ibd)[2]]),dim(LL_matrix_pedigree_ibd)[1],dim(LL_matrix_pedigree_ibd)[2]-2)
LL.matrix.data.label.6<-matrix(apply(LL.matrix.pedigree.ibd,1,sum),dim(LL.matrix.pedigree.ibd)[1],1)
LL.matrix.data.label.7=1/LL.matrix.data.label.6
LL.matrix.data.label.7[which(LL.matrix.data.label.7==Inf)]<-0
effect.estimation.combine.label<-Estimation_Data[c(unlist(LL_matrix_pedigree_ibd[,1])+1),-which(colnames(Estimation_Data) %in% c(unique(unlist(LL_matrix_pedigree_ibd[,2]))))]*LL.matrix.pedigree.ibd
Pop_Snp_Data_2<-Pop_Snp_Data[which(Pop_Snp_Data[,1]==unique(unlist(LL_matrix_pedigree_ibd[,2]))),4:dim(Pop_Snp_Data)[2]]
Pop_Snp_Data_3<-Pop_Snp_Data_2[,c(unlist(LL_matrix_pedigree_ibd[,1])+1)]
effect.estimation.combine.label.3<-matrix(apply(effect.estimation.combine.label,1,sum),dim(effect.estimation.combine.label)[1],1)
effect.estimation.combine.label.4=effect.estimation.combine.label.3*LL.matrix.data.label.7
prediction.y.2<-as.matrix(Pop_Snp_Data_3)%*%effect.estimation.combine.label.4
Imputated.Y<-ImputatePhenotype(as.matrix(Pop_Snp_Data[which(Pop_Snp_Data[,1] %in% c(unique(unlist(LL_matrix_pedigree_ibd[,2])))),3]))
cc2=cor(prediction.y.2,Imputated.Y)
Re2<-cbind(unique(unlist(LL_matrix_pedigree_ibd[,2])),w_ibd,w_no_ibd,cc2)
return(Re2)
}

matrix.based.on.pedigree.ibd.with.snp.index.prediction.re<-lapply(matrix.based.on.pedigree.ibd.with.snp.index,WeightIBDGCA5,combined.marker.estimation,pop.acc.height.snp.16,1,0)
matrix.based.on.pedigree.ibd.with.snp.index.prediction.re<-matrix(unlist(matrix.based.on.pedigree.ibd.with.snp.index.prediction.re),19,4,byrow=T)

save.image("Based_on_pedigree_ibd_prediction.RData")
savehistory("Based_on_pedigree_ibd_prediction.Rhistory")
q("yes")
