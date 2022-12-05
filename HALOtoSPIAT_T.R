library(stringr)
library(data.table)

all.files<-list.files(inputfiles, full.names=T, recursive=T, include.dirs=T)

#Files are in alphabetical order so all.files will run through Immune file then T cell file for same patient
for (i in all.files){
    newdat<-str_replace(i,pattern=".csv", replacement="_tableT.csv")
    print(paste0("Formatting", i))
    tableT<-fread(i,data.table=F)
    #get rid of Xmin Xmax Ymin Ymax
    tableT<- subset(tableT, select = -c(2:3, 9:12))
    #rename XY columns
    colnames(tableT)<-str_replace(colnames(tableT),pattern="Immune Registered XMin",replacement="XMin")
    colnames(tableT)<-str_replace(colnames(tableT),pattern="Immune Registered XMax",replacement="XMax")
    colnames(tableT)<-str_replace(colnames(tableT),pattern="Immune Registered YMin",replacement="YMin")
    colnames(tableT)<-str_replace(colnames(tableT),pattern="Immune Registered YMax",replacement="YMax")
    #rename P1-8 columns
    colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39\\+CD103\\+PD-1\\+",replacement="P1")
    colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39\\+CD103\\+PD-1-",replacement="P2")
    colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39\\+CD103-PD-1\\+",replacement="P3")
    colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39\\+CD103-PD-1-",replacement="P4")
    colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39-CD103\\+PD-1\\+",replacement="P5")
    colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39-CD103\\+PD-1-",replacement="P6")
    colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39-CD103-PD-1\\+",replacement="P7")
    colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39-CD103-PD-1-",replacement="P8")
    #Replace um2  
    colnames(tableT)<-str_replace_all(colnames(tableT),pattern="\\(µm²\\)",replacement="")
    colnames(tableT)<-str_replace_all(colnames(tableT),pattern="\\(µm\\)",replacement="")
    
    #Add CD4
    tableT$CD4<-0
    tableT$CD4[tableT$CD8..T.cells==0 & tableT$T.cells==1]<-1
    
    
    fwrite(tableT, file = newdat,quote=F,sep=",")
}