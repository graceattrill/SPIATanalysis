library(stringr)
library(data.table)

#This HALOtoSPIAT uses the Immune panel's Melanoma
all.files<-list.files(path=inputfiles, pattern=".csv", full.names=T, recursive=T, include.dirs=T)

#These 3 vectors are so here so that the columns match up between the two tables in the final rbind
T.colsp<-c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "CD8+ T cells", "T cells")
T.colsm<-c("CD39 (Opal 520) Positive Classification", "CD39 (Opal 520) Positive Cytoplasm Classification", "CD39 (Opal 520) Cytoplasm Intensity", "CD39 (Opal 520) Cell Intensity", 
           "PD-1 (Opal 540) Positive Classification", "PD-1 (Opal 540) Positive Cytoplasm Classification", "PD-1 (Opal 540) Cytoplasm Intensity", "PD-1 (Opal 540) Cell Intensity",
           "CD3 (Opal 570) Positive Classification", "CD3 (Opal 570) Positive Cytoplasm Classification", "CD3 (Opal 570) Cytoplasm Intensity", "CD3 (Opal 570) Cell Intensity",
           "CD103 (Opal 620) Positive Classification", "CD103 (Opal 620) Positive Cytoplasm Classification", "CD103 (Opal 620) Cytoplasm Intensity", "CD103 (Opal 620) Cell Intensity",
           "CD8 (Opal 650) Positive Classification", "CD8 (Opal 650) Positive Cytoplasm Classification", "CD8 (Opal 650) Cytoplasm Intensity", "CD8 (Opal 650) Cell Intensity")
I.cols<-c("B cell", "NK cell", "Langerhans cell", "HLAABCpMelanoma", "HLAABCnMelanoma",
          "CD20 (Opal 520) Positive Classification",	"CD20 (Opal 520) Positive Cytoplasm Classification",	"CD20 (Opal 520) Cytoplasm Intensity",	"CD20 (Opal 520) Cell Intensity",
          "CD1a (Opal 540) Positive Classification",	"CD1a (Opal 540) Positive Cytoplasm Classification",	"CD1a (Opal 540) Cytoplasm Intensity",	"CD1a (Opal 540) Cell Intensity",
          "CD56 (Opal 570) Positive Classification",	"CD56 (Opal 570) Positive Cytoplasm Classification",	"CD56 (Opal 570) Cytoplasm Intensity",	"CD56 (Opal 570) Cell Intensity",
          "HLA-ABC (Opal 620) Positive Classification",	"HLA-ABC (Opal 620) Positive Cytoplasm Classification",	"HLA-ABC (Opal 620) Cytoplasm Intensity",	"HLA-ABC (Opal 620) Cell Intensity",
          "Langerin (Opal 650) Positive Classification",	"Langerin (Opal 650) Positive Cytoplasm Classification",	"Langerin (Opal 650) Cytoplasm Intensity",	"Langerin (Opal 650) Cell Intensity")


#Files are in alphabetical order so all.files will run through Immune file then T cell file for same patient
#Code is written so that Immune table is formatted before T cell file - checkti and ifelse runs corresponding T cell/immune code on file depending on whether "Immune" is in files name
for (i in all.files){
  checkti<- str_detect(i, "Immune")
  print(checkti)
  if (checkti == TRUE){
    melpin<-str_replace(i,pattern="tableI.csv", replacement="_combitable.csv")
    print(paste0("Formatting", i))
    tableIp<-fread(i,data.table=F)
    #add xy columns to tableI
    tableI<-tableIp[,1:6]
    #add T cell phenotype columns to tableI
    out1<-matrix(0,nrow=(nrow(tableI)),ncol=10)
    colnames(out1)<-c(T.colsp)
    tableI<-cbind(tableI,out1)
    #Add melanoma column to tableI
    tableI$Melanoma<-tableIp[,"Melanoma"]
    #add T cell marker columns
    out2<-matrix(0,nrow=(nrow(tableI)),ncol=20)
    colnames(out2)<-c(T.colsm)
    tableI<-cbind(tableI,out2)
    #add SOX10,DAPI & cell parameters
    tableI<-cbind(tableI,tableIp[,33:46])
    #add rest of Immune columns to tableI
    tableI<-cbind(tableI, tableIp[,c(7:9, 11:32)])
    
    #Replace (um2)  
    colnames(tableI)<-str_replace_all(colnames(tableI),pattern="\\(µm²\\)",replacement="")
    colnames(tableI)<-str_replace_all(colnames(tableI),pattern="\\(µm\\)",replacement="")
  }
    #If "Immune" not in file name, T cell code will run
  else{
    print(paste0("Formatting", i))
    print("replace")
    tableTp<-fread(i,data.table=F)
    print(dim(tableTp))
    #Remove T cell melanoma
    tableTp<-subset(tableTp, select = -c(17))
    #Add columns back to table T
    tableT<-tableTp[,1:16]
    out3<-matrix(0,nrow=(nrow(tableT)), ncol = 1)
    colnames(out3)<-"Melanoma"
    tableT<-cbind(tableT,out3)
    tableT<-cbind(tableT, tableTp[,17:50])
    #add columns from immune with 0 in all rows
    out4<-matrix(0,nrow=(nrow(tableT)),ncol=25)
    colnames(out4)<-c(I.cols)
    print(dim(out4))
    tableT<-cbind(tableT,out4)
    print(dim(tableT))
    #Replace um2  
    colnames(tableT)<-str_replace_all(colnames(tableT),pattern="\\(µm²\\)",replacement="")
    colnames(tableT)<-str_replace_all(colnames(tableT),pattern="\\(µm\\)",replacement="")
    #Set max T cell object Id +1 as minimum object Id for Immune panel
    objectidn<-max(tableT$`Object Id`)+1
    #Add object ID number to immune cells
    tableI$`Object Id`<-tableI$`Object Id` + objectidn
    
    combitable<-rbind(tableT,tableI)
    print(dim(combitable))
    
    #Add tumour & stroma columns to combitable for clusters
    combitable$Tumour<-0
    combitable$Tumour[combitable$Classifier.Label=="Tumour"]<-1
    combitable$Stroma<-0
    combitable$Stroma[combitable$Classifier.Label=="Stroma/Epidermis"]<-1
    
    
    fwrite(combitable, file = paste(melpin),quote=F,sep=",")
  }
}