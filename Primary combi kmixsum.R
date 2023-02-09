library(data.table)
library(stringr)
library(spatstat)#for Kcross
library(Rphenograph) 
library(SPIAT)
setwd("/Users/gracie/Desktop/R Projects/SPIATprimary/CD4 output")

combi.files<-list.files(path="/Users/gracie/Desktop/R Projects/SPIATprimary/CD4 combifiles", pattern="cluster.csv", full.names=T, recursive=T, include.dirs=T)
iteration<-1

for(i in combi.files){

#define columns for SPIAT to read from
dye_columns_interest <- c("P1",
                          "P2", 
                          "P3", 
                          "P4", 
                          "P5", 
                          "P6", 
                          "P7", 
                          "P8",
                          "DAPI Positive Classification",
                          "B cell",
                          "NK cell",
                          "Langerhans cell",
                          "HLAABCpMelanoma",
                          "HLAABCnMelanoma",
                          "CD4")

#Setting dye_columns_interest & intensity_columns_interest as the same columns allows manually HALO-assigned phenotypes to be used
intensity_columns_interest <- dye_columns_interest

markers <- c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "DAPI", "B cell",
             "NK cell", "Langerhans cell", "HLAABCpMelanoma", "HLAABCnMelanoma","CD4")

#Convert file name to melpin
melpin<-str_remove(i,pattern = "extrainfilename")
melpin<-word(melpin, sep = "CD4")
print(paste(iteration,melpin))

#Set immune_phenotypes vector as all phenotypes except for Melanoma - simplifies later analysis
#Adjust for missing phenotypes for each patient, else:
immune_phenotypes<-c("P1","P2","P3","P4","P5","P6","P7", "P8", "CD4","B cell", "NK cell", "Langerhans cell")
T.phenotypes<-c("P1","P2","P3","P4","P5","P6","P7", "P8", "CD4")
phenotypes_of_interest<-c("P1","P2","P3","P4","P5","P6","P7", "P8", "CD4","B cell", "NK cell", "Langerhans cell","HLAABCpMelanoma", "HLAABCnMelanoma")
colour_vector<-c("chartreuse", "cyan", "blue", "darkorange", "deeppink", "chocolate4", "gold", "seagreen", "coral2","tan4", "violet", "azure2","yellow", "tan1")

formatted_image <- format_image_to_sce(format="HALO",
                                       path= i,
                                       markers=markers,
                                       dye_columns_interest=dye_columns_interest,
                                       intensity_columns_interest=intensity_columns_interest) 
print_column(formatted_image, column= "Phenotype")

#calculate cell proportions
p_cells <- calculate_cell_proportions(formatted_image, reference_celltypes=c("Total"))

#Remove any double positive B cell/NK cell or LC/NK cell and check phenotypes again
formatted_image<-select_phenotypes(formatted_image, keep=FALSE, phenotypes = c("B cell,NK cell", "NK cell,Langerhans cell"))
print_column(formatted_image, column = "Phenotype")

print("summary...")
#Calculate total summary distances
summary_distances<-calculate_summary_distances_between_cell_types(formatted_image,all_marker_combinations=TRUE)

if(iteration==1){
  all.summary.distances<-summary_distances
} else {
  all.summary.distances<-merge(all.summary.distances, summary_distances, by = c("Reference", "Nearest"), all = TRUE)
}

print("mixing...")
#Kcross & mixing scores with HLA & T cells

#Calculate mixing score - HLAABCp Melanoma to T - take both normalised and not normalised
mixscoretab<-data.frame(phenotype=character(),nmixingscore=numeric(),mixingscore=numeric(),Melpin=character(),HLA=character())
for(h in T.phenotypes){
  mixscorenum<-mixing_score_summary(formatted_image,reference_marker="HLAABCpMelanoma",target_marker=h,radius=20, column = "Phenotype")
  looptab<-data.frame(phenotype=paste0(h,"mix"), nmixingscore=mixscorenum$Normalised_mixing_score,mixingscore=mixscorenum$Mixing_score,Melpin=melpin,HLA="HLAp")
  mixscoretab<-rbind(mixscoretab,looptab)
}
if(iteration==1){
  all.mixing.scores<-mixscoretab
}else{
  all.mixing.scores<-rbind(all.mixing.scores,mixscoretab)
}

#Calculate mixing score - HLAABCn Melanoma to T
mixscoretab<-data.frame(phenotype=character(),nmixingscore=numeric(),mixingscore=numeric(),Melpin=character(),HLA=character())
for(k in T.phenotypes){
  mixscorenum<-mixing_score_summary(formatted_image,reference_marker="HLAABCnMelanoma",target_marker=k,radius=20, column = "Phenotype")
  looptab<-data.frame(phenotype=paste0(k,"mix"), nmixingscore=mixscorenum$Normalised_mixing_score,mixingscore=mixscorenum$Mixing_score,Melpin=melpin,HLA="HLAn")
  mixscoretab<-rbind(mixscoretab,looptab)
}
all.mixing.scores<-rbind(all.mixing.scores,mixscoretab)


print("crossing...")
#Kcross - HLAABCp Melanoma to T
Ktab<-data.frame(phenotype=character(),Kscore=numeric(),Melpin=character(),HLA=character())
for (l in T.phenotypes){
  df_cross<-bivariate_dependence(formatted_image, method = "Kcross", phenotypes = c(l, "HLAABCpMelanoma"), column = "Phenotype")
  pdf(paste0(melpin,l, "Kcross.pdf"), onefile = T, height = 10, width = 10)
  dev.off()
  Kcrossn<-AUC_of_cross_function(df_cross)
  Klooptab<-data.frame(phenotype=paste0(l, "Kcross"), Kscore = Kcrossn,Melpin=melpin,HLA="HLAp")
  Ktab<-rbind(Ktab, Klooptab)
  dev.off()
}

if(iteration==1){
  all.Kcross<-Ktab
} else{
  all.Kcross<-rbind(all.Kcross,Ktab)
}

#Kcross - HLAABCn Melanoma to T cell
Ktab<-data.frame(phenotype=character(),Kscore=numeric(),Melpin=character(),HLA=character())
#remove pts/T cell phenotypes as necessary
for (m in T.phenotypes){
  df_cross<-bivariate_dependence(formatted_image, method = "Kcross", phenotypes = c(m, "HLAABCnMelanoma"), column = "Phenotype")
  pdf(paste0(melpin,m, "Kcross.pdf"), onefile = T, height = 10, width = 10)
  dev.off()
  Kcrossn<-AUC_of_cross_function(df_cross)
  Klooptab<-data.frame(phenotype=paste0(m, "Kcross"), Kscore = Kcrossn,Melpin=melpin,HLA="HLAn")
  Ktab<-rbind(Ktab, Klooptab)
  dev.off()
}
all.Kcross<-rbind(all.Kcross,Ktab)


#Immune distance to nearest melanoma - re-do without HLA
  
  #define columns for SPIAT to read from
  dye_columns_interest <- c("P1",
                            "P2", 
                            "P3", 
                            "P4", 
                            "P5", 
                            "P6", 
                            "P7", 
                            "P8",
                            "DAPI Positive Classification",
                            "B cell",
                            "NK cell",
                            "Langerhans cell",
                            "Melanoma",
                            "CD4")
  
  #Setting dye_columns_interest & intensity_columns_interest as the same columns allows manually HALO-assigned phenotypes to be used
  intensity_columns_interest <- dye_columns_interest
  
  markers <- c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "DAPI", "B cell","NK cell", "Langerhans cell", "Melanoma","CD4")
  
#Adjust per patient based on missing phenotypes - else:
    immune_phenotypes<-c("P1","P2","P3","P4","P5","P6","P7", "P8", "B cell", "NK cell", "Langerhans cell","CD4")
    T.phenotypes<-c("P1","P2","P3","P4","P5","P6","P7", "P8","CD4")
    phenotypes_of_interest<-c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "B cell", "NK cell", "Langerhans cell","Melanoma","CD4")
    colour_vector<-c("chartreuse", "cyan", "blue", "darkorange", "deeppink", "chocolate4", "gold", "seagreen", "coral2","tan4", "violet", "azure2", "tan1")
  
  print(paste("Nearest neighbour analysis",iteration, melpin))
  
  formatted_image <- format_image_to_sce(format="HALO",
                                         path= i,
                                         markers=markers,
                                         dye_columns_interest=dye_columns_interest,
                                         intensity_columns_interest=intensity_columns_interest) 
  print_column(formatted_image, column= "Phenotype")
  #Remove any double positive B cell/NK cell or LC/NK cell and check phenotypes again
  formatted_image<-select_phenotypes(formatted_image, keep=FALSE, phenotypes = c("B cell,NK cell", "NK cell,Langerhans cell"))
  print_column(formatted_image, column = "Phenotype")
  
  avgdisttab<-data.frame(Reference=character(),Nearest=character(),Melpin=character(),Mean=numeric())
  for (n in immune_phenotypes){
    avgdist<-calculate_summary_distances_between_cell_types(formatted_image,column="Phenotype",all_marker_combinations=FALSE,combinations_of_interest=c(n,"Melanoma"))
    distlooptab<-data.frame(Reference="Melanoma",Nearest=n,Melpin=melpin,Mean=avgdist[2,3])
    if(n!="P1"){
      avgdisttab<-rbind(avgdisttab,distlooptab)
    } else{
      avgdisttab<-distlooptab
    }
  }
    if(iteration==1){
      all.avgdist<-avgdisttab
    } else {
  all.avgdist<-rbind(all.avgdist,avgdisttab)
    }
  
iteration<-iteration+1
}

#Prepare Kcross and mixing scores for joint export
col3<-rep(combi.files, each = 3)
col2<-rep(combi.files,each=2)
colnames(all.Kcross) = c("Phenotype", col2)
colnames(all.mixing.scores)<-c("Phenotype", col2)

Knmixprimarycombi<-rbind(all.Kcross,all.mixing.scores)
colnames(Knmixprimarycombi)<-str_remove(colnames(Knmixprimarycombi), pattern = "__combitable.csv")
colnames(Knmixprimarycombi)<-str_remove(colnames(Knmixprimarycombi), pattern = "/Volumes/Gracie/SPIAToutput/Combitable 1xmel/")
colnames(Knmixprimarycombi)<-word(colnames(Knmixprimarycombi), sep = " Immune.tif")
colnames(Knmixprimarycombi)<-word(colnames(Knmixprimarycombi), -1)
colnames(Knmixprimarycombi)[seq(2,ncol(Knmixprimarycombi),2)]<-paste0(colnames(Knmixprimarycombi)[seq(2,ncol(Knmixprimarycombi),2)]," HLAp")
colnames(Knmixprimarycombi)[seq(3,ncol(Knmixprimarycombi),2)]<-paste0(colnames(Knmixprimarycombi)[seq(3,ncol(Knmixprimarycombi),2)]," HLAn")
fwrite(Knmixprimarycombi,file="KmixprimarycombiCD4.csv", quote=F)

#Prepare summary distances for export
colnames(all.summary.distances) <- c("Target", "Nearest",col3)
colnames(all.summary.distances) <-str_remove(colnames(all.summary.distances), pattern = "__combitable.csv")
colnames(all.summary.distances)<-str_remove(colnames(all.summary.distances), pattern = "/Volumes/Gracie/SPIAToutput/Combitable 1xmel/")
colnames(all.summary.distances)<-word(colnames(all.summary.distances), sep = "Immune.tif")
col3num<-seq(3,length(col3),3)

#Take average of summary_distance data for each patient & put into single file
all.summary.distances$Allmean<-apply(all.summary.distances[,col3num],1,mean, na.rm=TRUE)
all.summary.distances$Allstddev<-apply(all.summary.distances[,col3num+1],1,mean, na.rm=TRUE)
all.summary.distances$Allmedian<-apply(all.summary.distances[,col3num+2,],1,mean, na.rm=TRUE)

#Last 3 columns [177:179] of avg.summ.dist contain averaged mean/stddev/median 
avg.summ.dist<-cbind(all.summary.distances[,1:2], all.summary.distances[,c(147:149)])
colnames(avg.summ.dist)<-c("Target", "Nearest", "Mean", "Std.Dev", "Median")

pdf(file = "primarysummarydistancescombiCD4.pdf", onefile = T, width = 15, height = 15)
plot_distance_heatmap(avg.summ.dist)
dev.off()

fwrite(avg.summ.dist,file="primarysummdistcombiCD4.csv", quote= F)


#Prepare average distances for export - convert pixels to um
all.avgdist$Mean<-all.avgdist$Mean/2
fwrite(all.avgdist,file="primaryminavgdistCD4.csv",quote=F)

#Edit files in excel to match up melpins & outcome data
good.summ.dist.all<-read.csv("editedgoodfile.csv",sep=",",header=T)
good.summ.dist.all<-as.data.frame(good.summ.dist.all)

col3num<-seq(3,length(colnames(good.summ.dist.all)),3)
good.summ.dist.all$Mean<-apply(good.summ.dist.all[,col3num],1,mean, na.rm=TRUE)
good.summ.dist.all$Std.Dev<-apply(good.summ.dist.all[,col3num+1],1,mean, na.rm=TRUE)
good.summ.dist.all$Median<-apply(good.summ.dist.all[,col3num+2,],1,mean, na.rm=TRUE)
head(good.summ.dist.all,3)
good.summ.dist<-cbind(good.summ.dist.all[,1:2],good.summ.dist.all[,93:95])
pdf(file = "primarysummarydistancescombigoodmean.pdf", onefile = T, width = 15, height = 15)
plot_distance_heatmap_GAscaled(good.summ.dist)
dev.off()

poor.summ.dist.all<-read.csv("editedpoorfile.csv",sep=",",header=T)
poor.summ.dist.all<-as.data.frame(poor.summ.dist.all)

col3num<-seq(3,length(colnames(poor.summ.dist.all)),3)
poor.summ.dist.all$Mean<-apply(poor.summ.dist.all[,col3num3],1,mean, na.rm=TRUE)
poor.summ.dist.all$Std.Dev<-apply(poor.summ.dist.all[,col3num+1],1,mean, na.rm=TRUE)
poor.summ.dist.all$Median<-apply(poor.summ.dist.all[,col3num+2,],1,mean, na.rm=TRUE)
head(poor.summ.dist.all,3)
poor.summ.dist<-cbind(poor.summ.dist.all[,1:2],poor.summ.dist.all[,72:74])
pdf(file = "primarysummarydistancescombipoormean.pdf", onefile = T, width = 15, height = 15)
plot_distance_heatmap_GAscaled(poor.summ.dist)
dev.off()

summ.dist.all<-read.csv("editedallfile.csv",sep=",",header=T)
summ.dist.all<-as.data.frame(summ.dist.all)

col3num<-seq(3,length(colnames(summ.dist.all)),3)
summ.dist.all$Mean<-apply(summ.dist.all[,col3num],1,mean, na.rm=TRUE)
summ.dist.all$Std.Dev<-apply(summ.dist.all[,col3num+1],1,mean, na.rm=TRUE)
summ.dist.all$Median<-apply(summ.dist.all[,col3num+2,],1,mean, na.rm=TRUE)
head(summ.dist.all,3)
summ.dist<-cbind(summ.dist.all[,1:2],summ.dist.all[,162:164])
pdf(file = "primarysummarydistancescombiallmean.pdf", onefile = T, width = 15, height = 15)
plot_distance_heatmap(summ.dist)
dev.off()

