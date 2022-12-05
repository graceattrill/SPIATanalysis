library(data.table)
library(stringr)
library(spatstat)#for Kcross
library(Rphenograph) 
library(SPIAT)

I.files<-list.files(path=inputfiles, pattern="tableI.csv", full.names=T, recursive=T, include.dirs=T)
iteration<-1

for (i in I.files){
  
  #Without HLA first
  dye_columns_interest <- c("B cell",
                            "NK cell",
                            "Langerhans cell",
                            "Melanoma",
                            "DAPI Positive Classification")
  
  intensity_columns_interest <- c("B cell",
                                  "NK cell",
                                  "Langerhans cell",
                                  "Melanoma",
                                  "DAPI Positive Classification")
  
  #Convert file name to melpin - check first step
  newdat<-str_replace(i, pattern = ".object_results_tableI.csv", replacement = "_I_SPIAT")
  newdat<-str_remove(newdat,pattern = "\\/Volumes\\/Gracie\\/SPIAToutput\\/Immune SPIAT data\\/")
  newdat<-word(newdat, sep = " Immune.tif")
  newdat<-word(newdat,-1)
  
  #Adjust input for missing phenotypes
    I.phenotypes<-c("B cell", "NK cell", "Langerhans cell")
    markers <- c("B cell", "NK cell", "Langerhans cell", "Melanoma", "DAPI")
    phenotypes_of_interest<-c("B cell", "NK cell", "Langerhans cell", "HLAABCpMelanoma", "HLAABCnMelanoma")
    colour_vector<-c("chartreuse", "cyan", "blue", "maroon","deeppink")
  
  formatted_image <- format_image_to_spe(format="HALO",
                                         path= i,
                                         markers=markers,
                                         dye_columns_interest=dye_columns_interest,
                                         intensity_columns_interest=intensity_columns_interest) 
  print_column(formatted_image, column= "Phenotype")
  formatted_image<-select_phenotypes(formatted_image, keep=FALSE, phenotypes = c("B cell,NK cell", "NK cell,Langerhans cell"))
  print_column(formatted_image, column = "Phenotype")
  
  
  #Kcross all Melanoma
  Ktab<-data.frame(phenotype=character(),Kscore=numeric())
  for (k in I.phenotypes){
    df_cross<-bivariate_dependence(formatted_image, method = "Kcross", phenotypes = c(k, "Melanoma"), column = "Phenotype")
    pdf(paste0(newdat,k, "Kcross.pdf"), onefile = T, height = 10, width = 10)
    dev.off()
    Kcrossn<-AUC_of_cross_function(df_cross)
    Klooptab<-data.frame(phenotype=paste0(k, "Kcross"), Kscore = Kcrossn)
    Ktab<-rbind(Ktab, Klooptab)
    dev.off()
  }
  
  #Calculate mixing score - Melanoma to immune
  mixscoretab<-data.frame(phenotype=character(),mixingscore=numeric())
  for(j in I.phenotypes){
    mixscorenum<-mixing_score_summary(formatted_image,reference_marker="Melanoma",target_marker=j,radius=20, column = "Phenotype")
    looptab<-data.frame(phenotype=j, score=mixscorenum$Normalised_mixing_score)
    mixscoretab<-rbind(mixscoretab,looptab)
  }
  if(iteration==1){
    all.mixing.scores<-mixscoretab
  }
  else{
    all.mixing.scores<-cbind(all.mixing.scores,mixscoretab[,2])
  }
  
  if(iteration==1){
    all.Kcross<-Ktab
  } else{
    all.Kcross<-cbind(all.Kcross,Ktab[,2])
  }
  
  #measure distances
  summary_distances<-calculate_summary_distances_between_phenotypes(formatted_image,all_marker_combinations=TRUE)
  
  if(iteration==1){
    all.summary.distances<-summary_distances
  } else {
    all.summary.distances<-merge(all.summary.distances, summary_distances, by = c("Target", "Nearest"), all = TRUE)
  }
  
  #Including HLA
  
  #Adjust input for missing phenotypes
  markers <- c("B cell", "NK cell", "Langerhans cell", "HLAABCpMelanoma", "HLAABCnMelanoma", "DAPI")
  
  dye_columns_interest <- c("B cell",
                            "NK cell",
                            "Langerhans cell",
                            "HLAABCpMelanoma",
                            "HLAABCnMelanoma",
                            "DAPI Positive Classification")
  
  intensity_columns_interest <- c("B cell",
                                  "NK cell",
                                  "Langerhans cell",
                                  "HLAABCpMelanoma",
                                  "HLAABCnMelanoma",
                                  "DAPI Positive Classification")
 
#upload file to R
formatted_image <- format_image_to_sce(format="HALO",
                                       path= i,
                                       markers=markers,
                                       dye_columns_interest=dye_columns_interest,
                                       intensity_columns_interest=intensity_columns_interest) 
print_column(formatted_image, column= "Phenotype")
formatted_image<-select_phenotypes(formatted_image, keep=FALSE, phenotypes = c("B cell,NK cell", "NK cell,Langerhans cell"))
print_column(formatted_image, column = "Phenotype")


#create plot with phenotypes
pdf(file=paste0(newdat, "Iplot.pdf"),onefile=T, height=15, width=15)
plot_cell_categories(formatted_image, phenotypes_of_interest, colour_vector)
dev.off()

#measure distances
summary_distances<-calculate_summary_distances_between_phenotypes(formatted_image,all_marker_combinations=TRUE)

if(iteration==1){
  all.summary.distances<-summary_distances
} else {
  all.summary.distances<-merge(all.summary.distances, summary_distances, by = c("Target", "Nearest"), all = TRUE)
}

#Calculate mixing score - HLAABCpMelanoma to immune
mixscoretabp<-data.frame(phenotype=character(),mixingscore=numeric())
for(j in I.phenotypes){
  mixscorenum<-mixing_score_summary(formatted_image,reference_marker="HLAABCpMelanoma",target_marker=j,radius=20, column = "Phenotype")
  looptab<-data.frame(phenotype=j, score=mixscorenum$Normalised_mixing_score)
  mixscoretabp<-rbind(mixscoretabp,looptab)
}

all.mixing.scores<-cbind(all.mixing.scores,mixscoretabp[,2])

#Calculate mixing score - HLAABCnMelanoma to immune
mixscoretabn<-data.frame(phenotype=character(),mixingscore=numeric())
for(j in I.phenotypes){
  mixscorenum<-mixing_score_summary(formatted_image,reference_marker="HLAABCnMelanoma",target_marker=j,radius=20, column = "Phenotype")
  looptab<-data.frame(phenotype=j, score=mixscorenum$Normalised_mixing_score)
  mixscoretabn<-rbind(mixscoretabn,looptab)
}

all.mixing.scores<-cbind(all.mixing.scores,mixscoretabn[,2])

#Kcross HLAABCpMelanmoma
Ktabp<-data.frame(phenotype=character(),Kscore=numeric())
for (k in I.phenotypes){
  df_cross<-bivariate_dependence(formatted_image, method = "Kcross", phenotypes = c(k, "HLAABCpMelanoma"), column = "Phenotype")
  pdf(paste0(newdat,k, "Kcross.pdf"), onefile = T, height = 10, width = 10)
  dev.off()
  Kcrossn<-AUC_of_cross_function(df_cross)
  Klooptab<-data.frame(phenotype=paste0(k, "Kcross"), Kscore = Kcrossn)
  Ktabp<-rbind(Ktabp, Klooptab)
  dev.off()
}

  all.Kcross<-cbind(all.Kcross,Ktabp[,2])


#Kcross HLAABCnMelanoma
Ktabn<-data.frame(phenotype=character(),Kscore=numeric())
for (k in I.phenotypes){
  df_cross<-bivariate_dependence(formatted_image, method = "Kcross", phenotypes = c(k, "HLAABCnMelanoma"), column = "Phenotype")
  pdf(paste0(newdat,k, "Kcross.pdf"), onefile = T, height = 10, width = 10)
  dev.off()
  Kcrossn<-AUC_of_cross_function(df_cross)
  Klooptab<-data.frame(phenotype=paste0(k, "Kcross"), Kscore = Kcrossn)
  Ktabn<-rbind(Ktabn, Klooptab)
  dev.off()
}

  all.Kcross<-cbind(all.Kcross,Ktabn[,2])
}


#Prepare Kcross and mixing scores for joint export
col3<-rep(I.files, each = 3)
colnames(all.Kcross) = c("Phenotype", col3)

colnames(all.mixing.scores)<-c("Phenotype", col3)

KnmixprimaryI<-rbind(all.Kcross,all.mixing.scores)
colnames(KnmixprimaryI)<-str_remove(colnames(KnmixprimaryI), pattern = "_object_results_tableI.csv")
colnames(KnmixprimaryI)<-str_remove(colnames(KnmixprimaryI), pattern = "/Volumes/Gracie/SPIAToutput/Immune SPIAT data/")
colnames(KnmixprimaryI)<-word(colnames(KnmixprimaryI), sep = " Immune.tif")
colnames(KnmixprimaryI)<-word(colnames(KnmixprimaryI), -1)
colnames(KnmixprimaryI)[seq(2,ncol(KnmixprimaryI),3)]<-paste0(colnames(KnmixprimaryI)[seq(2,ncol(KnmixprimaryI),3)]," Mel")
colnames(KnmixprimaryI)[seq(3,ncol(KnmixprimaryI),3)]<-paste0(colnames(KnmixprimaryI)[seq(3,ncol(KnmixprimaryI),3)]," HLAp")
colnames(KnmixprimaryI)[seq(4,ncol(KnmixprimaryI),3)]<-paste0(colnames(KnmixprimaryI)[seq(4,ncol(KnmixprimaryI),3)]," HLAn")
fwrite(KnmixprimaryI,file="KmixprimaryI.csv", quote=F)

#Prepare summary distances for export
colnames(all.summary.distances) <- c("Target", "Nearest",col3)
colnames(all.summary.distances) <-str_remove(colnames(all.summary.distances), pattern = ".object_results_tableI.csv")
colnames(all.summary.distances)<-str_remove(colnames(all.summary.distances), pattern = "/Volumes/Gracie/SPIAToutput/Immune SPIAT data/")
colnames(all.summary.distances)<-word(colnames(all.summary.distances), sep = "Immune.tif")
col3num<-seq(3,length(col3),3)

#Take average of summary_distance data for each patient & put into single file
all.summary.distances$Allmean<-apply(all.summary.distances[,col3num],1,mean, na.rm=TRUE)
all.summary.distances$Allstddev<-apply(all.summary.distances[,col3num+1],1,mean, na.rm=TRUE)
all.summary.distances$Allmedian<-apply(all.summary.distances[,col3num+2,],1,mean, na.rm=TRUE)

#Last 3 columns [177:179] of avg.summ.dist contain averaged mean/stddev/median 
avg.summ.dist<-cbind(all.summary.distances[,1:2], all.summary.distances[,177:179])
colnames(avg.summ.dist)<-c("Target", "Nearest", "Mean", "Std.Dev", "Median")

pdf(file = "primarysummarydistancesI.pdf", onefile = T, width = 15, height = 15)
plot_distance_heatmap(avg.summ.dist)
dev.off()

fwrite(all.summary.distances,file="primaryalldistI.csv", quote=F)
fwrite(avg.summ.dist,file="primarysummdistI.csv", quote= F)