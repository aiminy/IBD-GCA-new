#load("./data/CombinedMarkerEffectEstimation.RData")

SplitSNP2IbdNoibd<-function(S,w_ibd,w_no_ibd){

n_ibd=S$maping_count
n_no_ibd=18-n_ibd

A=matrix(c(S$pos),1,length(c(S$pos)))
B=matrix(c(rep(w_ibd,n_ibd)),1,length(c(rep(w_ibd,n_ibd))))
C=matrix(c(rep(w_no_ibd,n_no_ibd)),1,length(c(rep(w_no_ibd,n_no_ibd))))

re<-cbind(A,B,C)

return(re)

}

ApplySplitSNP2IbdNoibd2OnePop<-function(Pop,w_ibd,w_no_ibd){

temp<-apply(Pop,1,SplitSNP2IbdNoibd,w_ibd,w_no_ibd)
temp<-t(temp)
return(temp)

}

WeightIBDGCA<-function(Mapping_Data,Estimation_Data,Pop_Snp_Data,w_ibd,w_no_ibd){

label.ibd.no.ibd<-ApplySplitSNP2IbdNoibd2OnePop(Mapping_Data,w_ibd,w_no_ibd)

label.ibd.no.ibd.2<-matrix(as.numeric(label.ibd.no.ibd),955690,19)

effect.estimation.combine.label.1008<-Estimation_Data[,2:dim(Estimation_Data)[2]]*label.ibd.no.ibd.2[,2:dim(label.ibd.no.ibd.2)[2]]

prediction.y<-as.matrix(Pop_Snp_Data[which(Pop_Snp_Data[,1]=="1008"),4:dim(Pop_Snp_Data)[2]])%*%as.matrix(apply(effect.estimation.combine.label.1008,1,mean))

cc=cor(prediction.y,as.matrix(Pop_Snp_Data[which(Pop_Snp_Data[,1]=="1008"),3]))

cat(w_ibd,w_no_ibd,cc,"\n")

}

GenerateLabelMatrix<-function(num_row,SS){
label.matrix<-matrix(rep(0,num_row*19),num_row,19)
colnames(label.matrix)<-SS
return(label.matrix)
}


ReplaceWithLabel4OneRow<-function(A,p,q){
B=A
B[5:length(B)][which(names(B[5:length(B)]) %in% B$mapped_pop)]<-p
B[5:length(B)][-which(names(B[5:length(B)]) %in% c(B$target_pop,B$mapped_pop))]<-q
B<-t(B)
return(B)
}

ImputatePhenotype<-function(Phenotype.Y.2){
NA.index=which(is.na(Phenotype.Y.2))
if(length(NA.index)!=0){
imputation.value<-mean(Phenotype.Y.2[-which(is.na(Phenotype.Y.2))])
Phenotype.Y.2[NA.index]<-imputation.value
}
return(Phenotype.Y.2)
}


WeightIBDGCA2<-function(Mapping_Data,Estimation_Data,Pop_Snp_Data,w_ibd,w_no_ibd){
LL.matrix<-GenerateLabelMatrix(dim(Mapping_Data)[1],colnames(Estimation_Data))
LL.matrix.data.label<-cbind(Mapping_Data,LL.matrix)
LL.matrix.data.label.2<-apply(LL.matrix.data.label,1,ReplaceWithLabel4OneRow,w_ibd,w_no_ibd)
LL.matrix.data.label.3 <- do.call(rbind,LL.matrix.data.label.2)
LL.matrix.data.label.4<-matrix(as.numeric(LL.matrix.data.label.3[,5:23]),length(as.numeric(LL.matrix.data.label.3[,5:23]))/19,19)
effect.estimation.combine.label<-Estimation_Data[c(unlist(Mapping_Data[,1])+1),]*LL.matrix.data.label.4
effect.estimation.combine.label.2<-effect.estimation.combine.label[,-which(colnames(effect.estimation.combine.label) %in% c(unique(Mapping_Data[,3])[[1]]))]
Pop_Snp_Data_2<-Pop_Snp_Data[which(Pop_Snp_Data[,1]==unique(Mapping_Data[,3])[[1]]),4:dim(Pop_Snp_Data)[2]]
Pop_Snp_Data_3<-Pop_Snp_Data_2[,c(unlist(Mapping_Data[,1])+1)]
prediction.y<-as.matrix(Pop_Snp_Data_3)%*%as.matrix(apply(effect.estimation.combine.label.2,1,mean))
Imputated.Y<-ImputatePhenotype(as.matrix(Pop_Snp_Data[which(Pop_Snp_Data[,1] %in% c(unique(Mapping_Data[,3])[[1]])),3]))
cc=cor(prediction.y,Imputated.Y)
Re<-cbind(unique(Mapping_Data[,3])[[1]],w_ibd,w_no_ibd,cc)
return(Re)
}


GiveRandomWeight<-function(i){
w_ibd=runif(1)
w_no_ibd=runif(1)
Re.all.pop.2<-do.call(rbind,lapply(Re.refined.mapping.new,WeightIBDGCA2,combined.marker.estimation,pop.acc.height.snp.16,as.character(w_ibd),as.character(w_no_ibd)))
return(Re.all.pop.2)
}

WeightIBDGCA3<-function(Mapping_Data,Estimation_Data,Pop_Snp_Data,w_ibd,w_no_ibd){

LL.matrix<-GenerateLabelMatrix(dim(Mapping_Data)[1],colnames(Estimation_Data))

LL.matrix.data.label<-cbind(Mapping_Data,LL.matrix)


LL.matrix.data.label.2<-apply(LL.matrix.data.label,1,ReplaceWithLabel4OneRow,w_ibd,w_no_ibd)

#print(head(LL.matrix.data.label.2))

LL.matrix.data.label.3 <- do.call(rbind,LL.matrix.data.label.2)

print(head(LL.matrix.data.label.3))


LL.matrix.data.label.4<-matrix(as.numeric(LL.matrix.data.label.3[,5:23]),length(as.numeric(LL.matrix.data.label.3[,5:23]))/19,19)

#print(head(LL.matrix.data.label.4))

colnames(LL.matrix.data.label.4)=colnames(LL.matrix.data.label.3)[5:23]

#print(head(LL.matrix.data.label.4))

#effect.estimation.combine.label<-combined.marker.estimation*LL.matrix.data.label.4
LL.matrix.data.label.5<-LL.matrix.data.label.4[,-which(colnames(LL.matrix.data.label.4) %in% c(unique(Mapping_Data[,3])[[1]]))]

#print(head(LL.matrix.data.label.5))

LL.matrix.data.label.6<-matrix(apply(LL.matrix.data.label.5,1,sum),dim(LL.matrix.data.label.5)[1],1)

#print(head(LL.matrix.data.label.6))

LL.matrix.data.label.7=1/LL.matrix.data.label.6

LL.matrix.data.label.7[which(LL.matrix.data.label.7==Inf)]<-0

#cat(dim(LL.matrix.data.label.7))
#print(head(LL.matrix.data.label.7))

effect.estimation.combine.label<-Estimation_Data[c(unlist(Mapping_Data[,1])+1),]*LL.matrix.data.label.4

#unique(Re.refined.mapping.new[[1]][,3])[[1]]

#pop.acc.height.snp.16

#which(colnames(effect.estimation.combine.label) %in% c(unique(Re.refined.mapping.new[[1]][,3])[[1]]))

effect.estimation.combine.label.2<-effect.estimation.combine.label[,-which(colnames(effect.estimation.combine.label) %in% c(unique(Mapping_Data[,3])[[1]]))]

Pop_Snp_Data_2<-Pop_Snp_Data[which(Pop_Snp_Data[,1]==unique(Mapping_Data[,3])[[1]]),4:dim(Pop_Snp_Data)[2]]
Pop_Snp_Data_3<-Pop_Snp_Data_2[,c(unlist(Mapping_Data[,1])+1)]


effect.estimation.combine.label.3<-matrix(apply(effect.estimation.combine.label.2,1,sum),dim(effect.estimation.combine.label.2)[1],1)

effect.estimation.combine.label.4=effect.estimation.combine.label.3*LL.matrix.data.label.7

#prediction.y<-as.matrix(Pop_Snp_Data_3)%*%as.matrix(apply(effect.estimation.combine.label.2,1,mean))

prediction.y.2<-as.matrix(Pop_Snp_Data_3)%*%effect.estimation.combine.label.4

Imputated.Y<-ImputatePhenotype(as.matrix(Pop_Snp_Data[which(Pop_Snp_Data[,1] %in% c(unique(Mapping_Data[,3])[[1]])),3]))

#cc=cor(prediction.y,Imputated.Y)
cc2=cor(prediction.y.2,Imputated.Y)

#Re<-cbind(unique(Mapping_Data[,3])[[1]],w_ibd,w_no_ibd,cc)
Re2<-cbind(unique(Mapping_Data[,3])[[1]],w_ibd,w_no_ibd,cc2)


Re3<-list(snp=Pop_Snp_Data_3,snp2=Pop_Snp_Data_2,Pre=prediction.y.2,Real=Imputated.Y,est2=LL.matrix.data.label.6,est3=LL.matrix.data.label.7,est=effect.estimation.combine.label.4,Divide_mapping_depth=Re2)

return(Re3)

#return(Re2)

}


GenerateLabelMatrixBasedOnPedigree<-function(Mapping_Data,population_parents_name,Estimation_Data){

s1<-as.character(population_parents_name[which(population_parents_name[,1] %in% unique(Mapping_Data[,3][[1]])),2])
 
s2<-as.character(population_parents_name[which(population_parents_name[,1] %in% c(colnames(Estimation_Data))),2])

s3<-intersect(s1,s2)

matached.pop.name<-unique(as.character(population_parents_name[which(population_parents_name[,2] %in% c(s3)),1]))

s4<-colnames(Estimation_Data)[which(colnames(Estimation_Data) %in% c(matached.pop.name))]

LL.matrix<-GenerateLabelMatrix(dim(Mapping_Data)[1],colnames(Estimation_Data))

#LL.matrix[,which(colnames(LL.matrix) %in% c(s4))]<-1

#print(head(LL.matrix))

s5<-s4[-which(s4 %in% c(unique(Mapping_Data[,3][[1]])))]

LL.matrix[,which(colnames(LL.matrix) %in% c(s5))]<-1

LL.matrix.2<-LL.matrix[,-which(colnames(LL.matrix) %in% c(unique(Mapping_Data[,3][[1]])))]

#print(head(LL.matrix.2))

#print(tail(LL.matrix.2))

#cat("target pop:",c(unique(Mapping_Data[,3][[1]])),"\n")

#cat(s5,"\n")

return(LL.matrix.2)

}

WeightIBDGCA4<-function(Mapping_Data,population_parents_name,Estimation_Data,Pop_Snp_Data,w_ibd,w_no_ibd){

LL.matrix<-GenerateLabelMatrix(dim(Mapping_Data)[1],colnames(Estimation_Data))

#LL.matrix.data.label<-cbind(Mapping_Data,LL.matrix)

LL.matrix.data.label<-cbind(Mapping_Data,LL.matrix)

LL.matrix.data.label.2<-apply(LL.matrix.data.label,1,ReplaceWithLabel4OneRow,w_ibd,w_no_ibd)

#print(head(LL.matrix.data.label.2))

LL.matrix.data.label.3 <- do.call(rbind,LL.matrix.data.label.2)

#print(head(LL.matrix.data.label.3))

LL.matrix.based.on.pedigree<-GenerateLabelMatrixBasedOnPedigree(Mapping_Data,population_parents_name,Estimation_Data)

#print(head(LL.matrix.based.on.pedigree))

#print(tail(LL.matrix.based.on.pedigree))

LL.matrix.data.label.4<-matrix(as.numeric(LL.matrix.data.label.3[,5:23]),length(as.numeric(LL.matrix.data.label.3[,5:23]))/19,19)

#print(head(LL.matrix.data.label.4))

colnames(LL.matrix.data.label.4)=colnames(LL.matrix.data.label.3)[5:23]

#print(head(LL.matrix.data.label.4))


#effect.estimation.combine.label<-combined.marker.estimation*LL.matrix.data.label.4
LL.matrix.data.label.5<-LL.matrix.data.label.4[,-which(colnames(LL.matrix.data.label.4) %in% c(unique(Mapping_Data[,3])[[1]]))]

#print(head(LL.matrix.data.label.5))

LL.matrix.pedigree.ibd<-LL.matrix.based.on.pedigree+LL.matrix.data.label.5

#LL.matrix.data.label.6<-matrix(apply(LL.matrix.data.label.5,1,sum),dim(LL.matrix.data.label.5)[1],1)

LL.matrix.data.label.6<-matrix(apply(LL.matrix.pedigree.ibd,1,sum),dim(LL.matrix.pedigree.ibd)[1],1)

#print(head(LL.matrix.data.label.6))

LL.matrix.data.label.7=1/LL.matrix.data.label.6

LL.matrix.data.label.7[which(LL.matrix.data.label.7==Inf)]<-0

#cat(dim(LL.matrix.data.label.7))
#print(head(LL.matrix.data.label.7))

#effect.estimation.combine.label<-Estimation_Data[c(unlist(Mapping_Data[,1])+1),]*LL.matrix.data.label.4

effect.estimation.combine.label<-Estimation_Data[c(unlist(Mapping_Data[,1])+1),-which(colnames(Estimation_Data) %in% c(unique(Mapping_Data[,3])[[1]]))]*LL.matrix.pedigree.ibd

#unique(Re.refined.mapping.new[[1]][,3])[[1]]

#pop.acc.height.snp.16

#which(colnames(effect.estimation.combine.label) %in% c(unique(Re.refined.mapping.new[[1]][,3])[[1]]))

#effect.estimation.combine.label.2<-effect.estimation.combine.label[,-which(colnames(effect.estimation.combine.label) %in% c(unique(Mapping_Data[,3])[[1]]))]

Pop_Snp_Data_2<-Pop_Snp_Data[which(Pop_Snp_Data[,1]==unique(Mapping_Data[,3])[[1]]),4:dim(Pop_Snp_Data)[2]]
Pop_Snp_Data_3<-Pop_Snp_Data_2[,c(unlist(Mapping_Data[,1])+1)]


effect.estimation.combine.label.3<-matrix(apply(effect.estimation.combine.label,1,sum),dim(effect.estimation.combine.label)[1],1)

effect.estimation.combine.label.4=effect.estimation.combine.label.3*LL.matrix.data.label.7

#prediction.y<-as.matrix(Pop_Snp_Data_3)%*%as.matrix(apply(effect.estimation.combine.label.2,1,mean))

prediction.y.2<-as.matrix(Pop_Snp_Data_3)%*%effect.estimation.combine.label.4

Imputated.Y<-ImputatePhenotype(as.matrix(Pop_Snp_Data[which(Pop_Snp_Data[,1] %in% c(unique(Mapping_Data[,3])[[1]])),3]))

#cc=cor(prediction.y,Imputated.Y)
cc2=cor(prediction.y.2,Imputated.Y)

#Re<-cbind(unique(Mapping_Data[,3])[[1]],w_ibd,w_no_ibd,cc)
Re2<-cbind(unique(Mapping_Data[,3])[[1]],w_ibd,w_no_ibd,cc2)


#Re3<-list(snp=Pop_Snp_Data_3,snp2=Pop_Snp_Data_2,Pre=prediction.y.2,Real=Imputated.Y,est2=LL.matrix.data.label.6,est3=LL.matrix.data.label.7,est=effect.estimation.combine.label.4,Divide_mapping_depth=Re2)

#Re4<-list(m1=LL.matrix.based.on.pedigree,m2=LL.matrix.data.label.5,m3=LL.matrix.pedigree.ibd)

#return(Re3)

#return(Re4)

return(Re2)
}