library(data.table)
library(stringr)
library(spatstat)#for Kcross
library(Rphenograph) 
library(SPIAT)
setwd("adjuvantoutput")

T.files<-list.files(path="adjuvantinput", pattern="tableTadj.csv", full.names=T, recursive=T, include.dirs=T)

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
                            "DAPI Positive",
                            "CD4")
  
  intensity_columns_interest <- dye_columns_interest
  markers <- c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "Melanoma", "DAPI","CD4")
  
  
  #Convert file name to melpin
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
  
  #Adjust input if T cell populations are missing
    T.phenotypes<-c("P1","P2","P3","P4","P5","P6","P7", "P8","CD4")
    categories_of_interest<-c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "Melanoma","CD4")
    colour_vector<-c("chartreuse", "cyan", "blue", "darkorange", "deeppink", "chocolate4", "gold", "seagreen", "maroon","grey")
  
  immune.phenotypes<-c(T.phenotypes, "Melanoma")
  
  print(paste(iteration,melpin))
  
  #upload file to R
  formatted_image <- format_image_to_spe(format="HALO",
                                         path= i,
                                         markers=markers,
                                         dye_columns_interest=dye_columns_interest,
                                         intensity_columns_interest=intensity_columns_interest) 
  print(unique(formatted_image$Phenotype))
  
  #calculate cell proportions
  p_cells <- calculate_cell_proportions(formatted_image, reference_celltypes=c("Total"))
  
  #create plot with phenotypes
  tiff(file=paste0(melpin, "plotfigure.tif"), height=500, width=500)
  plot_cell_categories(formatted_image, categories_of_interest, colour_vector,"Phenotype")
  dev.off()
  
  #measure distances
  print("summary...")
  min_dist<-calculate_minimum_distances_between_celltypes(formatted_image,feature_colname="Phenotype",cell_types_of_interest=immune.phenotypes)
  summary_distances<-calculate_summary_distances_between_celltypes(min_dist)
  colnames(summary_distances)<-c("Pair",paste("Mean",melpin), paste("Min", melpin), paste("Max",melpin), paste("Std.Dev",melpin),paste("Median",melpin),"Reference","Target")
  
  if(iteration==1){
    all.summary.distances<-summary_distances
  } else {
    all.summary.distances<-merge(all.summary.distances, summary_distances, by = c("Pair", "Reference", "Target"), all = TRUE)
  }
  
  #Calculate mixing score - Melanoma to immune
  print("mixing...")
  mixscoretab<-data.frame(phenotype=character(),nmixingscore=numeric(),mixingscore=numeric(),Melpin=character())
  for(j in T.phenotypes){
    mixscorenum<-mixing_score_summary(formatted_image,reference_celltype="Melanoma",target_celltype=j,radius=20, feature_colname = "Phenotype")
    looptab<-data.frame(phenotype=paste0(j,"mix"), nmixingscore=mixscorenum$Normalised_mixing_score,mixingscore=mixscorenum$Mixing_score,Melpin=melpin)
    mixscoretab<-rbind(mixscoretab,looptab)
  }
  if(iteration==1){
    all.mixing.scores<-mixscoretab
  } else{
    all.mixing.scores<-rbind(all.mixing.scores,mixscoretab)
  }
  
  #Entropy
  print("Entropy...")
  ent<-calculate_entropy(formatted_image,T.phenotypes)
  enttab<-data.table(melpin=melpin,entropy=ent)
  if (iteration==1){
    all.entropy<-enttab
  } else {
    all.entropy<-rbind(all.entropy,enttab)
  }
  
  ##CD8+ T cell only
  print("CD8s")
  dye_columns_interest <- c("CD8+ T cell", "Melanoma", "DAPI Positive")
  intensity_columns_interest <- dye_columns_interest
  markers<-c("CD8+ T cell", "Melanoma", "DAPI")
  formatted_image <- format_image_to_spe(format="HALO",
                                         path= i,
                                         markers=markers,
                                         dye_columns_interest=dye_columns_interest,
                                         intensity_columns_interest=intensity_columns_interest) 
  
  CD8min_dist<-calculate_minimum_distances_between_celltypes(formatted_image,feature_colname="Phenotype",cell_types_of_interest=c("CD8+ T cell", "Melanoma"))
  CD8.dist<-calculate_summary_distances_between_celltypes(CD8min_dist)
  colnames(CD8.dist)<-colnames(summary_distances)
  
  if (iteration==1){
    all.CD8.dist<-CD8.dist
  } else{
    all.CD8.dist<-merge(all.CD8.dist, CD8.dist, by = c("Pair", "Reference", "Target"), all = TRUE)}
  
  
  mixscorenum<-mixing_score_summary(formatted_image,reference_celltype="Melanoma",target_celltype="CD8+ T cell",radius=20, feature_colname = "Phenotype")
  looptab<-data.frame(phenotype=paste0("CD8+ T cell","mix"), nmixingscore=mixscorenum$Normalised_mixing_score,mixingscore=mixscorenum$Mixing_score,Melpin=melpin)
  
  all.mixing.scores<-rbind(all.mixing.scores,looptab)
  
  ##T cell only 
  print("T cell")
  
  dye_columns_interest <- c("T cell", "Melanoma", "DAPI Positive")
  intensity_columns_interest <- dye_columns_interest
  markers<-c("T cell", "Melanoma", "DAPI")
  
  formatted_image <- format_image_to_spe(format="HALO",
                                         path= i,
                                         markers=markers,
                                         dye_columns_interest=dye_columns_interest,
                                         intensity_columns_interest=intensity_columns_interest) 
  
  #calculate cell proportions
  Tmin_dist<-calculate_minimum_distances_between_celltypes(formatted_image,feature_colname="Phenotype",cell_types_of_interest=c("T cell", "Melanoma"))
  T.dist<-calculate_summary_distances_between_celltypes(Tmin_dist)
  colnames(T.dist)<-colnames(summary_distances)
  
  if (iteration==1){
    all.T.dist<-T.dist
  } else{
    all.T.dist<-merge(all.T.dist, T.dist, by = c("Pair", "Reference", "Target"), all = TRUE)}
  
  mixscorenum<-mixing_score_summary(formatted_image,reference_celltype="Melanoma",target_celltype="T cell",radius=20, feature_colname = "Phenotype")
  looptab<-data.frame(phenotype=paste0(j,"mix"), nmixingscore=mixscorenum$Normalised_mixing_score,mixingscore=mixscorenum$Mixing_score,Melpin=melpin)
  all.mixing.scores<-rbind(all.mixing.scores,looptab)
  
  iteration<-iteration+1
} 

all.summary.distances<-rbind(all.summary.distances,all.CD8.dist,all.T.dist)

#Prepare Kcross and mixing scores for joint export
labels<-word(all.mixing.scores$phenotype,sep="mix")
all.mixing.scores$label<-labels
fwrite(all.mixing.scores,file="MixadjSPIAT.csv", quote=F)

#Prepare summary distances for export
col5<-rep(T.files, each = 5)
col5num<-seq(5,length(col3),5)

#Take average of summary_distance data for each patient & put into single file
all.summary.distances$Allmean<-apply(all.summary.distances[,col5num-1],1,mean, na.rm=TRUE)
all.summary.distances$Allstddev<-apply(all.summary.distances[,col5num+2],1,mean, na.rm=TRUE)
all.summary.distances$Allmedian<-apply(all.summary.distances[,col5num+3,],1,mean, na.rm=TRUE)

#Last 3 columns of avg.summ.dist contain averaged mean/stddev/median 
avg.summ.dist<-cbind(all.summary.distances[,1:3], all.summary.distances[,599:601])
colnames(avg.summ.dist)<-c("Pair","Reference","Target","Mean","Std.Dev","Median")
avg.summ.dist<-avg.summ.dist[avg.summ.dist$Reference != "T cell" & avg.summ.dist$Target != "T cell",]
avg.summ.dist<-avg.summ.dist[avg.summ.dist$Reference != "CD8+ T cell" & avg.summ.dist$Target != "CD8+ T cell",]

pdf(file = "adjuvantsummdist.pdf", onefile = T, width = 15, height = 15)
plot_distance_heatmap(avg.summ.dist)
dev.off()

fwrite(all.summary.distances,file="adjuvantalldist.csv", quote=F)
fwrite(avg.summ.dist,file="adjuvantsummdist.csv", quote= F)
fwrite(all.entropy,file="adjuvantallentropy.csv",quote=F,sep=",")


#Prepare good vs bad files externally to R, prep GAscaled function

rf.summ.dist.all<-read.csv("AdjuvantRFalldist.csv",sep=",",header=T)
rf.summ.dist.all<-as.data.frame(rf.summ.dist.all)
rf.summ.dist.all<-rf.summ.dist.all[3:96,]
rf.summ.dist.all<-rf.summ.dist.all[rf.summ.dist.all$Reference != "T cell" & rf.summ.dist.all$Nearest != "T cell",]
rf.summ.dist.all<-rf.summ.dist.all[rf.summ.dist.all$Reference != "CD8+ T cell" & rf.summ.dist.all$Nearest != "CD8+ T cell",]

col3num<-seq(3,length(colnames(rf.summ.dist.all)),3)
rf.summ.dist.all$AllMean<-apply(rf.summ.dist.all[,col3num],1,mean, na.rm=TRUE)
rf.summ.dist.all$AllStd.Dev<-apply(rf.summ.dist.all[,col3num+1],1,mean, na.rm=TRUE)
rf.summ.dist.all$AllMedian<-apply(rf.summ.dist.all[,col3num+2,],1,mean, na.rm=TRUE)
rf.summ.dist<-as.data.frame(cbind(rf.summ.dist.all[,1:2],rf.summ.dist.all[,234:236]))
colnames(rf.summ.dist)<-c("Reference","Nearest","Mean","Std.Dev","Median")
rf.summ.dist$Mean<-as.numeric(rf.summ.dist$Mean)
print(max(rf.summ.dist$Mean))

pdf(file = "adjuvantsummarydistancesrf.pdf", onefile = T, width = 15, height = 15)
plot_distance_heatmap_GAscaled(rf.summ.dist)
dev.off()


r.summ.dist.all<-read.csv("AdjuvantRalldist.csv",sep=",",header=T)
r.summ.dist.all<-as.data.frame(r.summ.dist.all)
r.summ.dist.all<-r.summ.dist.all[3:96,]
r.summ.dist.all<-r.summ.dist.all[r.summ.dist.all$Reference != "T cell" & r.summ.dist.all$Nearest != "T cell",]
r.summ.dist.all<-r.summ.dist.all[r.summ.dist.all$Reference != "CD8+ T cell" & r.summ.dist.all$Nearest != "CD8+ T cell",]

col3num<-seq(3,length(colnames(r.summ.dist.all)),3)
r.summ.dist.all$AllMean<-apply(r.summ.dist.all[,col3num],1,mean, na.rm=TRUE)
r.summ.dist.all$AllStd.Dev<-apply(r.summ.dist.all[,col3num+1],1,mean, na.rm=TRUE)
r.summ.dist.all$AllMedian<-apply(r.summ.dist.all[,col3num+2,],1,mean, na.rm=TRUE)
r.summ.dist<-as.data.frame(cbind(r.summ.dist.all[,1:2],r.summ.dist.all[,123:125]))
colnames(r.summ.dist)<-c("Reference","Nearest","Mean","Std.Dev","Median")
r.summ.dist$Mean<-as.numeric(r.summ.dist$Mean)
print(max(r.summ.dist$Mean))


pdf(file = "adjuvantsummarydistancesr.pdf", onefile = T, width = 15, height = 15)
plot_distance_heatmap_GAscaled(r.summ.dist)
dev.off()
