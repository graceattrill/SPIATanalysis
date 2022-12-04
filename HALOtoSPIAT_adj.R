##Adjuvant cohort file conversion

library(stringr)
library(data.table)

setwd("output")

all.files<-list.files(path="adjuvantinput", full.names=T, recursive=T, include.dirs=T)

iteration<-1

for (i in all.files){
  melpin<-word(i, sep = "/",-1)
  melpin<-word(melpin,sep="_job",1)
  melpin<-word(melpin,sep="_CD39",1)
  melpin<-word(melpin,sep=".tif",1)
  melpin<-str_remove(melpin,"T cell")
  if (iteration>3 & iteration !=50){
    melpin<-word(melpin,sep="_",1)
  } else{
    melpin<-str_remove(melpin,"_")
  }
  
  print(paste(iteration, "Formatting", melpin))
  tableT<-fread(i,data.table=F)
  
  #Replace um2  
  colnames(tableT)<-str_remove(colnames(tableT),pattern="\\(µm²\\)")
  colnames(tableT)<-str_remove(colnames(tableT),pattern="\\(µm\\)")
  
  #Formatting depends on whether phenotyping was performed in HALO & naming of phenotypes in HALO
  #rename P1-8 columns
  
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39\\+CD103\\+PD-1\\+",replacement="P1")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39\\+CD103\\+PD-1-",replacement="P2")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39\\+CD103-PD-1\\+",replacement="P3")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39\\+CD103-PD-1-",replacement="P4")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39-CD103\\+PD-1\\+",replacement="P5")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39-CD103\\+PD-1-",replacement="P6")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39-CD103-PD-1\\+",replacement="P7")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD39-CD103-PD-1-",replacement="P8")
  
  
  colnames(tableT)<-str_replace(colnames(tableT),pattern="P1.CD8..T.cell",replacement="P1")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="P2.CD8..T.cell",replacement="P2")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="P3.CD8..T.cell",replacement="P3")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="P4.CD8..T.cell",replacement="P4")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="P5.CD8..T.cell",replacement="P5")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="P6.CD8..T.cell",replacement="P6")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="P7.CD8..T.cell",replacement="P7")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="P8.CD8..T.cell",replacement="P8")
  
  colnames(tableT)<-str_replace(colnames(tableT),pattern="CD8..T.cell",replacement="CD8+ T cell")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="T.cell",replacement="T cell")
  
  colnames(tableT)<-str_replace(colnames(tableT),pattern="DAPI Positive Classification",replacement="DAPI Positive")
  colnames(tableT)<-str_replace(colnames(tableT),pattern="DAPI.Positive",replacement="DAPI Positive")
  
  
  if (any(str_detect(colnames(tableT),"T cell"))!=TRUE & any(str_detect(colnames(tableT),"T.cell"))!=TRUE){
    print(any(str_detect(colnames(tableT),"T cell")))
    print(any(str_detect(colnames(tableT),"T.cell")))

    #Add T cell
    tableT$`T cell`<-0
    tableT$`T cell`[tableT$`DAPI Positive`==1 & tableT$`CD3 (Opal 570) Positive`==1 & tableT$`SOX10 (Opal 690) Positive Nucleus`==0]<-1
    #Add CD8
    tableT$`CD8+ T cell`<-0
    tableT$`CD8+ T cell`[tableT$`T cell`==1 & tableT$`CD8 (Opal 650) Positive`==1]<-1

    #Add P1-P8
    tableT$P1<-0
    tableT$P1[tableT$`CD8+ T cell`==1 & tableT$`CD39 (Opal 520) Positive`==1 & tableT$`CD103 (Opal 620) Positive`==1 & tableT$`PD-1 (Opal 540) Positive`==1]<-1
    
    tableT$P2<-0
    tableT$P2[tableT$`CD8+ T cell`==1 & tableT$`CD39 (Opal 520) Positive`==1 & tableT$`CD103 (Opal 620) Positive`==1 & tableT$`PD-1 (Opal 540) Positive`==0]<-1
    
    tableT$P3<-0
    tableT$P3[tableT$`CD8+ T cell`==1 & tableT$`CD39 (Opal 520) Positive`==1 & tableT$`CD103 (Opal 620) Positive`==0 & tableT$`PD-1 (Opal 540) Positive`==1]<-1
    
    tableT$P4<-0
    tableT$P4[tableT$`CD8+ T cell`==1 & tableT$`CD39 (Opal 520) Positive`==1 & tableT$`CD103 (Opal 620) Positive`==0 & tableT$`PD-1 (Opal 540) Positive`==0]<-1
    
    tableT$P5<-0
    tableT$P5[tableT$`CD8+ T cell`==1 & tableT$`CD39 (Opal 520) Positive`==0 & tableT$`CD103 (Opal 620) Positive`==1 & tableT$`PD-1 (Opal 540) Positive`==1]<-1
    
    tableT$P6<-0
    tableT$P6[tableT$`CD8+ T cell`==1 & tableT$`CD39 (Opal 520) Positive`==0 & tableT$`CD103 (Opal 620) Positive`==1 & tableT$`PD-1 (Opal 540) Positive`==0]<-1
    
    tableT$P7<-0
    tableT$P7[tableT$`CD8+ T cell`==1 & tableT$`CD39 (Opal 520) Positive`==0 & tableT$`CD103 (Opal 620) Positive`==0 & tableT$`PD-1 (Opal 540) Positive`==1]<-1
    
    tableT$P8<-0
    tableT$P8[tableT$`CD8+ T cell`==1 & tableT$`CD39 (Opal 520) Positive`==0 & tableT$`CD103 (Opal 620) Positive`==0 & tableT$`PD-1 (Opal 540) Positive`==0]<-1

    #Add Melanoma
    tableT$Melanoma<-0
    tableT$Melanoma[tableT$`SOX10 (Opal 690) Positive Nucleus`==1 & tableT$`DAPI Positive`==1]<-1
  }
  
  
  #Add CD4
  tableT$CD4<-0
  tableT$CD4[tableT$`CD8+ T cell`==0 & tableT$`T cell`==1]<-1
  
  #Check if image contains melanoma - if image does not contain melanoma, do not write
  melcheck<-any(str_detect(tableT$Melanoma,pattern="1"))
  print(melcheck)
  
  #Report missing CD8 phenotypes in missingphen table
  if (iteration==1){
    missingphen<-data.frame(melpin=character(),phenotype=character())
  }
  
  if(any(str_detect(tableT$`CD8+ T cell`,pattern="1"))==FALSE){
    missing<-data.frame(melpin=melpin,phenotype="CD8")
    missingphen<-rbind(missing,missingphen)}
  if(any(str_detect(tableT$P1,pattern="1"))==FALSE){
    missing<-data.frame(melpin=melpin,phenotype="P1")
    missingphen<-rbind(missing,missingphen)}
if(any(str_detect(tableT$P2,pattern="1"))==FALSE){
  missing<-data.frame(melpin=melpin,phenotype="P2")
  missingphen<-rbind(missing,missingphen)}
if(any(str_detect(tableT$P3,pattern="1"))==FALSE){
  missing<-data.frame(melpin=melpin,phenotype="P3")
  missingphen<-rbind(missing,missingphen)}
if(any(str_detect(tableT$P4,pattern="1"))==FALSE){
  missing<-data.frame(melpin=melpin,phenotype="P4")
  missingphen<-rbind(missing,missingphen)}
if(any(str_detect(tableT$P5,pattern="1"))==FALSE){
  missing<-data.frame(melpin=melpin,phenotype="P5")
  missingphen<-rbind(missing,missingphen)}
if(any(str_detect(tableT$P6,pattern="1"))==FALSE){
  missing<-data.frame(melpin=melpin,phenotype="P6")
  missingphen<-rbind(missing,missingphen)}
if(any(str_detect(tableT$P7,pattern="1"))==FALSE){
  missing<-data.frame(melpin=melpin,phenotype="P7")
  missingphen<-rbind(missing,missingphen)}
if(any(str_detect(tableT$P8,pattern="1"))==FALSE){
  missing<-data.frame(melpin=melpin,phenotype="P8")
  missingphen<-rbind(missing,missingphen)}

  if(melcheck==TRUE){
  fwrite(tableT, file = melpin,quote=F,sep=",")
  }
  iteration<-iteration+1
}