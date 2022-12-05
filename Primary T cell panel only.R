library(data.table)
library(stringr)
library(spatstat)#for Kcross
library(Rphenograph) 
library(SPIAT)
setwd("primaryoutput")

T.files<-list.files(path="primaryinput",pattern="tableT.csv", full.names=T, recursive=T,include.dirs = T)

iteration<-1

for (i in T.files){
  
  #define columns for SPIAT to read from
  dye_columns_interest <- c("P1",
                            "P2", 
                            "P3", 
                            "P4", 
                            "P5", 
                            "P6", 
                            "P7", 
                            "P8", 
                            "Melanoma", 
                            "DAPI Positive Classification",
                            "CD4")
  
  intensity_columns_interest <- dye_columns_interest
  
  #Convert file name to melpin
  melpin<-word(i, sep = "/",-1)
  melpin<-word(melpin,sep=" T cell.tif")
  melpin<-word(melpin,-1)
  print(paste(iteration,melpin))
  
  #Adjust input if T cell populations are missing
    markers <- c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "Melanoma", "DAPI","CD4")
    T.phenotypes<-c("P1","P2","P3","P4","P5","P6","P7", "P8","CD4")
    phenotypes_of_interest<-c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "Melanoma","CD4")
    colour_vector<-c("chartreuse", "cyan", "blue", "darkorange", "deeppink", "chocolate4", "gold", "seagreen", "maroon","grey")  

#upload file to R
formatted_image <- format_image_to_spe(format="HALO",
                                       path= i,
                                       markers=markers,
                                       dye_columns_interest=dye_columns_interest,
                                       intensity_columns_interest=intensity_columns_interest) 
print_column(formatted_image, column= "Phenotype")

#create plot with phenotypes
pdf(file=paste0(melpin, "Tplot.pdf"),onefile=T, height=15, width=15)
plot_cell_categories_grey(formatted_image, phenotypes_of_interest, colour_vector)
dev.off()

#measure distances
print("summary...")
summary_distances<-calculate_summary_distances_between_cell_types(formatted_image,all_marker_combinations=TRUE)

if(iteration==1){
  all.summary.distances<-summary_distances
} else {
  all.summary.distances<-merge(all.summary.distances, summary_distances, by = c("Reference", "Nearest"), all = TRUE)
}

#Calculate mixing score - Melanoma to immune
print("mixing...")
mixscoretab<-data.frame(phenotype=character(),nmixingscore=numeric(),mixingscore=numeric(),Melpin=character())
for(j in T.phenotypes){
  mixscorenum<-mixing_score_summary(formatted_image,reference_marker="Melanoma",target_marker=j,radius=20, column = "Phenotype")
  looptab<-data.frame(phenotype=paste0(j,"mix"), nmixingscore=mixscorenum$Normalised_mixing_score,mixingscore=mixscorenum$Mixing_score,Melpin=melpin)
  mixscoretab<-rbind(mixscoretab,looptab)
}
if(iteration==1){
  all.mixing.scores<-mixscoretab
}
else{
  all.mixing.scores<-rbind(all.mixing.scores,mixscoretab)
}


#Kcross
print("crossing...")
Ktab<-data.frame(phenotype=character(),Kscore=numeric(),Melpin=character())
for (k in T.phenotypes){
  df_cross<-bivariate_dependence(formatted_image, method = "Kcross", phenotypes = c(k, "Melanoma"), column = "Phenotype")
  pdf(paste0(melpin,k, "Kcross.pdf"), onefile = T, height = 10, width = 10)
  dev.off()
  Kcrossn<-AUC_of_cross_function(df_cross)
  Klooptab<-data.frame(phenotype=paste0(k, "Kcross"), Kscore = Kcrossn,Melpin=melpin)
  Ktab<-rbind(Ktab, Klooptab)
  dev.off()
}

if(iteration==1){
  all.Kcross<-Ktab
} else{
  all.Kcross<-rbind(all.Kcross,Ktab)
}
  
print("CD8s")
##CD8+ T cell only
dye_columns_interest <- c("CD8+ T cells", "Melanoma", "DAPI Positive Classification")
intensity_columns_interest <- dye_columns_interest
markers<-c("CD8+ T cells", "Melanoma", "DAPI")
formatted_image <- format_image_to_spe(format="HALO",
                                       path= i,
                                       markers=markers,
                                       dye_columns_interest=dye_columns_interest,
                                       intensity_columns_interest=intensity_columns_interest) 
print_column(formatted_image, column= "Phenotype")

  mixscoretab<-data.frame(phenotype=character(),nmixingscore=numeric(),mixingscore=numeric(),Melpin=character())
  mixscorenum<-mixing_score_summary(formatted_image,reference_marker="Melanoma",target_marker="CD8+ T cells",radius=20, column = "Phenotype")
    looptab<-data.frame(phenotype=paste0("CD8+ T cell","mix"), nmixingscore=mixscorenum$Normalised_mixing_score,mixingscore=mixscorenum$Mixing_score,Melpin=melpin)
    mixscoretab<-rbind(mixscoretab,looptab)
    
    all.mixing.scores<-rbind(all.mixing.scores,mixscoretab)
    
    #Kcross
    Ktab<-data.frame(phenotype=character(),Kscore=numeric(),Melpin=character())
      df_cross<-bivariate_dependence(formatted_image, method = "Kcross", phenotypes = c("CD8+ T cells", "Melanoma"), column = "Phenotype")
      pdf(paste0(melpin,"CD8+ T cells", "Kcross.pdf"), onefile = T, height = 10, width = 10)
      dev.off()
      Kcrossn<-AUC_of_cross_function(df_cross)
      Klooptab<-data.frame(phenotype=paste0("CD8+ T cell", "Kcross"), Kscore = Kcrossn,Melpin=melpin)
      Ktab<-rbind(Ktab, Klooptab)
      dev.off()
      
      all.Kcross<-rbind(all.Kcross,Ktab)
      
      avgdisttab<-data.frame(Reference=character(),Nearest=character(),Melpin=character(),Mean=numeric())
        avgdist<-calculate_summary_distances_between_cell_types(formatted_image,column="Phenotype",all_marker_combinations=FALSE,combinations_of_interest=c("CD8+ T cells","Melanoma"))
        avgdisttab<-data.frame(Reference="Melanoma",Nearest="CD8+ T cell",Melpin=melpin,Mean=avgdist[2,3])
      if(iteration==1){
        all.avgdist<-avgdisttab
      } else {
        all.avgdist<-rbind(all.avgdist,avgdisttab)
      }

##T cell only 
print("T cells")

dye_columns_interest <- c("T cells", "Melanoma", "DAPI Positive Classification")
intensity_columns_interest <- dye_columns_interest
markers<-c("T cells", "Melanoma", "DAPI")
formatted_image <- format_image_to_spe(format="HALO",
                                       path= i,
                                       markers=markers,
                                       dye_columns_interest=dye_columns_interest,
                                       intensity_columns_interest=intensity_columns_interest) 
print_column(formatted_image, column= "Phenotype")

#calculate cell proportions
mixscoretab<-data.frame(phenotype=character(),nmixingscore=numeric(),mixingscore=numeric(),Melpin=character())
mixscorenum<-mixing_score_summary(formatted_image,reference_marker="Melanoma",target_marker="T cells",radius=20, column = "Phenotype")
looptab<-data.frame(phenotype=paste0("T cells","mix"), nmixingscore=mixscorenum$Normalised_mixing_score,mixingscore=mixscorenum$Mixing_score,Melpin=melpin)
mixscoretab<-rbind(mixscoretab,looptab)

all.mixing.scores<-rbind(all.mixing.scores,mixscoretab)

#Kcross
Ktab<-data.frame(phenotype=character(),Kscore=numeric(),Melpin=character())
  df_cross<-bivariate_dependence(formatted_image, method = "Kcross", phenotypes = c("T cells", "Melanoma"), column = "Phenotype")
  pdf(paste0(melpin,"T cells", "Kcross.pdf"), onefile = T, height = 10, width = 10)
  dev.off()
  Kcrossn<-AUC_of_cross_function(df_cross)
  Klooptab<-data.frame(phenotype=paste0("T cells", "Kcross"), Kscore = Kcrossn,Melpin=melpin)
  Ktab<-rbind(Ktab, Klooptab)
  dev.off()
  
  all.Kcross<-rbind(all.Kcross,Ktab)

  avgdist<-calculate_summary_distances_between_cell_types(formatted_image,column="Phenotype",all_marker_combinations=FALSE,combinations_of_interest=c("T cells","Melanoma"))
  distlooptab<-data.frame(Reference="Melanoma",Nearest="T cell",Melpin=melpin,Mean=avgdist[2,3])
  avgdisttab<-distlooptab
  all.avgdist<-rbind(all.avgdist,avgdisttab)

  iteration<-iteration+1
}


#Prepare Kcross and mixing scores for joint export

labels<-word(all.mixing.scores$phenotype,sep="mix")
all.mixing.scores$label<-labels
labels<-word(all.Kcross$phenotype,sep="Kcross")
all.Kcross$label<-labels
all.kmix<-merge(all.mixing.scores,all.Kcross,by=c("Melpin","label"))
fwrite(all.kmix,file="KmixprimarycombiCD4.csv", quote=F)

#Prepare summary distances for export
col3<-rep(T.files, each = 3)
col3num<-seq(3,length(col3),3)

colnames(all.summary.distances) <- c("Reference", "Nearest",col3)
colnames(all.summary.distances)<-word(colnames(all.summary.distances), sep = "/",-1)
colnames(all.summary.distances)<-word(colnames(all.summary.distances),sep=" T cell")
colnames(all.summary.distances)<-word(colnames(all.summary.distances),-1)


#Take average of summary_distance data for each patient & put into single file
all.summary.distances$Allmean<-apply(all.summary.distances[,col3num],1,mean, na.rm=TRUE)
all.summary.distances$Allstddev<-apply(all.summary.distances[,col3num+1],1,mean, na.rm=TRUE)
all.summary.distances$Allmedian<-apply(all.summary.distances[,col3num+2,],1,mean, na.rm=TRUE)

#Last 3 columns of avg.summ.dist contain averaged mean/stddev/median 
avg.summ.dist<-cbind(all.summary.distances[,1:2], all.summary.distances[,162:164])
colnames(avg.summ.dist)<-c("Reference", "Nearest", "Mean", "Std.Dev", "Median")

pdf(file = "primarysummarydistancescombiCD4.pdf", onefile = T, width = 15, height = 15)
plot_distance_heatmap(avg.summ.dist)
dev.off()

fwrite(all.summary.distances,file="primaryalldistT.csv", quote=F)
fwrite(avg.summ.dist,file="primarysummdistT.csv", quote= F)

#Export average distances
fwrite(all.avgdist,file="primaryCD8Tavgdist.csv",quote=F)

