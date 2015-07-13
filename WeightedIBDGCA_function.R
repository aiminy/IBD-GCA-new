#load("./data/CombinedMarkerEffectEstimation.RData")

SplitSNP2IbdNoibd<-function(S,w_ibd,w_no_ibd){

#cat(S$pos,"\t",S$maping_count,"\t",S$target_pop,"\t",S$mapped_pop,"\n")

n_ibd=S$maping_count
n_no_ibd=18-n_ibd

#cat(length(S$pos),"\t",length(rep("1",n_ibd)),length(rep("0",n_no_ibd)),"\n")

A=matrix(c(S$pos),1,length(c(S$pos)))
B=matrix(c(rep(w_ibd,n_ibd)),1,length(c(rep(w_ibd,n_ibd))))
C=matrix(c(rep(w_no_ibd,n_no_ibd)),1,length(c(rep(w_no_ibd,n_no_ibd))))

#matrix(data = NA, nrow = 1, ncol = 1, byrow = FALSE,dimnames = NULL)

re<-cbind(A,B,C)

return(re)

}

ApplySplitSNP2IbdNoibd2OnePop<-function(Pop,w_ibd,w_no_ibd){

#temp<-matrix()

temp<-apply(Pop,1,SplitSNP2IbdNoibd,w_ibd,w_no_ibd)
temp<-t(temp)
return(temp)

}

WeightIBDGCA<-function(Mapping_Data,Estimation_Data,Pop_Snp_Data,w_ibd,w_no_ibd){

label.ibd.no.ibd<-ApplySplitSNP2IbdNoibd2OnePop(Mapping_Data,w_ibd,w_no_ibd)

label.ibd.no.ibd.2<-matrix(as.numeric(label.ibd.no.ibd),955690,19)

#cat(dim(label.ibd.no.ibd.2),"\n")
#print(label.ibd.no.ibd.2[1:20,])

#print(Estimation_Data[1:20,2:19])

effect.estimation.combine.label.1008<-Estimation_Data[,2:dim(Estimation_Data)[2]]*label.ibd.no.ibd.2[,2:dim(label.ibd.no.ibd.2)[2]]

#cat(dim(effect.estimation.combine.label.1008),"\n")

#print(effect.estimation.combine.label.1008[1:20,])


prediction.y<-as.matrix(Pop_Snp_Data[which(Pop_Snp_Data[,1]=="1008"),4:dim(Pop_Snp_Data)[2]])%*%as.matrix(apply(effect.estimation.combine.label.1008,1,mean))

cc=cor(prediction.y,as.matrix(Pop_Snp_Data[which(Pop_Snp_Data[,1]=="1008"),3]))

cat(w_ibd,w_no_ibd,cc,"\n")

}

#WeightIBDGCA(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,"0","1")
#WeightIBDGCA(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,"0.1","0.9")
#WeightIBDGCA(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,"0.2","0.8")
#WeightIBDGCA(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,"0.3","0.7")
#WeightIBDGCA(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,"0.4","0.6")
#WeightIBDGCA(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,"0.5","0.5")
#WeightIBDGCA(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,"0.6","0.4")
#WeightIBDGCA(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,"0.7","0.3")
#WeightIBDGCA(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,"0.8","0.2")
#WeightIBDGCA(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,"0.9","0.1")
#WeightIBDGCA(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,"1.0","0")


#for(i in 1:100){

#w_ibd=runif(1)
#w_no_ibd=runif(1)

#WeightIBDGCA(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,as.character(w_ibd),as.character(w_no_ibd))

#}

GenerateLabelMatrix<-function(num_row,SS){
label.matrix<-matrix(rep(0,num_row*19),num_row,19)
colnames(label.matrix)<-SS
return(label.matrix)
}

#LL.matrix<-GenerateLabelMatrix(colnames(combined.marker.estimation))
#LL.matrix.2<-LL.matrix[,-which(colnames(LL.matrix) %in% Re.refined.mapping.new[[1]][1,3]$target_pop)]

ReplaceWithLabel4OneRow<-function(A,p,q){

B=A

#names(LL.matrix.data.label[1,5:length(LL.matrix.data.label[1,])])

B[5:length(B)][which(names(B[5:length(B)]) %in% B$mapped_pop)]<-p

B[5:length(B)][-which(names(B[5:length(B)]) %in% c(B$target_pop,B$mapped_pop))]<-q

B<-t(B)

return(B)

}

#LL.matrix.3<-t(apply(LL.matrix.2,1,ReplaceWithLabel4OneRow,0.7,0.4))
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

#ReplaceWithLabel4OneRow(LL.matrix.data.label[1,],0.7,0.4)

LL.matrix.data.label.2<-apply(LL.matrix.data.label,1,ReplaceWithLabel4OneRow,w_ibd,w_no_ibd)
LL.matrix.data.label.3 <- do.call(rbind,LL.matrix.data.label.2)

#LL.matrix.data.label.4<-matrix(as.numeric(LL.matrix.data.label.3[,5:23]),955690,19)

LL.matrix.data.label.4<-matrix(as.numeric(LL.matrix.data.label.3[,5:23]),length(as.numeric(LL.matrix.data.label.3[,5:23]))/19,19)

#effect.estimation.combine.label<-combined.marker.estimation*LL.matrix.data.label.4

effect.estimation.combine.label<-Estimation_Data[c(unlist(Mapping_Data[,1])+1),]*LL.matrix.data.label.4

#unique(Re.refined.mapping.new[[1]][,3])[[1]]

#pop.acc.height.snp.16

#which(colnames(effect.estimation.combine.label) %in% c(unique(Re.refined.mapping.new[[1]][,3])[[1]]))

effect.estimation.combine.label.2<-effect.estimation.combine.label[,-which(colnames(effect.estimation.combine.label) %in% c(unique(Mapping_Data[,3])[[1]]))]
Pop_Snp_Data_2<-Pop_Snp_Data[which(Pop_Snp_Data[,1]==unique(Mapping_Data[,3])[[1]]),4:dim(Pop_Snp_Data)[2]]
Pop_Snp_Data_3<-Pop_Snp_Data_2[,c(unlist(Mapping_Data[,1])+1)]

prediction.y<-as.matrix(Pop_Snp_Data_3)%*%as.matrix(apply(effect.estimation.combine.label.2,1,mean))

Imputated.Y<-ImputatePhenotype(as.matrix(Pop_Snp_Data[which(Pop_Snp_Data[,1] %in% c(unique(Mapping_Data[,3])[[1]])),3]))

cc=cor(prediction.y,Imputated.Y)

Re<-cbind(unique(Mapping_Data[,3])[[1]],w_ibd,w_no_ibd,cc)

return(Re)

}

#WeightIBDGCA2(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,"0.7","0.4")

#Re.all.pop<-do.call(rbind,lapply(Re.refined.mapping.new,WeightIBDGCA2,combined.marker.estimation,pop.acc.height.snp.16,"0.7","0.4"))

#for(i in 1:100){
#
#w_ibd=runif(1)
#w_no_ibd=runif(1)
#
#WeightIBDGCA(Re.refined.mapping.new[[1]],combined.marker.estimation,pop.acc.height.snp.16,as.character(w_ibd),as.character(w_no_ibd))
#Re.all.pop.2<-do.call(rbind,lapply(Re.refined.mapping.new,WeightIBDGCA2,combined.marker.estimation,pop.acc.height.snp.16,as.character(w_ibd),as.character(w_no_ibd)))
#
#}

GiveRandomWeight<-function(i){
w_ibd=runif(1)
w_no_ibd=runif(1)

#cat(w_ibd,w_no_ibd,"\n")

#Re.all.pop.2<-cbind(w_ibd,w_no_ibd)

Re.all.pop.2<-do.call(rbind,lapply(Re.refined.mapping.new,WeightIBDGCA2,combined.marker.estimation,pop.acc.height.snp.16,as.character(w_ibd),as.character(w_no_ibd)))
return(Re.all.pop.2)
}



WeightIBDGCA3<-function(Mapping_Data,Estimation_Data,Pop_Snp_Data,w_ibd,w_no_ibd){

LL.matrix<-GenerateLabelMatrix(dim(Mapping_Data)[1],colnames(Estimation_Data))

LL.matrix.data.label<-cbind(Mapping_Data,LL.matrix)

#ReplaceWithLabel4OneRow(LL.matrix.data.label[1,],0.7,0.4)

LL.matrix.data.label.2<-apply(LL.matrix.data.label,1,ReplaceWithLabel4OneRow,w_ibd,w_no_ibd)

#print(head(LL.matrix.data.label.2))

LL.matrix.data.label.3 <- do.call(rbind,LL.matrix.data.label.2)

#print(head(LL.matrix.data.label.3))



#LL.matrix.data.label.4<-matrix(as.numeric(LL.matrix.data.label.3[,5:23]),955690,19)

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

#Random_index<-seq(1,1,1)

#Results4GiveRandomWeight<-sapply(Random_index,GiveRandomWeight)

#Re.all.pop.combined<-rbind(Re.all.pop.combined,do.call(rbind,lapply(Re.refined.mapping.new,WeightIBDGCA2,combined.marker.estimation,pop.acc.height.snp.16,0.5,0.5)))

#for(i in 1:20){

#w_ibd=runif(1)
#w_no_ibd=runif(1)
#
#Re.all.pop.combined<-rbind(Re.all.pop.combined,do.call(rbind,lapply(Re.refined.mapping.new,WeightIBDGCA2,combined.marker.estimation,pop.acc.height.snp.16,as.character(w_ibd),as.character(w_no_ibd))))
#
#}

#for(i in 1:20){

#w_ibd=runif(1)
#w_ibd=1
#w_no_ibd=runif(1)
#
#Re.all.pop.combined<-rbind(Re.all.pop.combined,do.call(rbind,lapply(Re.refined.mapping.new,WeightIBDGCA2,combined.marker.estimation,pop.acc.height.snp.16,as.character(w_ibd),as.character(w_no_ibd))))
#
#}

