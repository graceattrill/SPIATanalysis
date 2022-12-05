library(stringr)
library(data.table)

all.files<-list.files(path="/Users/gracie/Desktop/Immune SPIAT data", pattern="_object_results.csv", full.names=T, recursive=T, include.dirs=T)

for (i in all.files){
newdat<-str_replace(i,pattern=".csv", replacement="_tableI.csv")
print(paste0("Formatting", i))
tableI<-fread(i,data.table=F)
#get rid of XMin XMax YMin YMax, (& Analysis Region & Layer - T & I exported from different HALO versions)
tableI<- subset(tableI, select = -c(2:3,9:12))
#rename XY columns, HLA-ABC (& title - different HALO versions)
colnames(tableI)<-str_replace(colnames(tableI),pattern="Immune Registered XMin",replacement="XMin")
colnames(tableI)<-str_replace(colnames(tableI),pattern="Immune Registered XMax",replacement="XMax")
colnames(tableI)<-str_replace(colnames(tableI),pattern="Immune Registered YMin",replacement="YMin")
colnames(tableI)<-str_replace(colnames(tableI),pattern="Immune Registered YMax",replacement="YMax")
colnames(tableI)<-str_replace(colnames(tableI), pattern="HLA-ABC\\+ Melanoma", replacement="HLAABCpMelanoma")
colnames(tableI)<-str_replace(colnames(tableI), pattern="HLA-ABC- Melanoma", replacement="HLAABCnMelanoma")
colnames(tableI)<-str_replace(colnames(tableI),pattern="Image Location",replacement="Image File Name")

#Replace (um2)  
colnames(tableI)<-str_replace_all(colnames(tableI),pattern="\\(µm²\\)",replacement="")
colnames(tableI)<-str_replace_all(colnames(tableI),pattern="\\(µm\\)",replacement="")

fwrite(tableI, file = newdat,quote=F,sep=",")
}