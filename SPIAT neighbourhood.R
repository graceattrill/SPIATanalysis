library(data.table)
library(stringr)
library(tidyr)
library(spatstat)
library(Rphenograph)
library(SPIAT)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
set.seed(25)

setwd(primaryoutput)

combi.files<-list.files(path="primaryinput", pattern="combitable.csv", full.names = T, recursive = T, include.dirs = T)
tumour_phenotypes<-c("P1,Tumour","P2,Tumour","P3,Tumour","P4,Tumour","P5,Tumour","P6,Tumour","P7,Tumour","P8,Tumour","B cell,Tumour",
                     "NK cell,Tumour","Langerhans cell,Tumour","CD4,Tumour")
stroma_phenotypes<-c("P1,Stroma","P2,Stroma","P3,Stroma","P4,Stroma","P5,Stroma","P6,Stroma","P7,Stroma","P8,Stroma","B cell,Stroma",
                     "NK cell,Stroma","Langerhans cell,Stroma", "CD4,Stroma")
melanoma_phenotypes<-c("Melanoma,Tumour","Melanoma,Stroma")
neighbour_phenotypes<-c(tumour_phenotypes,stroma_phenotypes)
immune_phenotypes<-c("P1","P2","P3","P4","P5","P6","P7","P8","B cell",
                     "NK cell","Langerhans cell","CD4")

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
                          "CD4",
                          "Tumour",
                          "Stroma")

intensity_columns_interest <- dye_columns_interest

markers <- c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "DAPI", 
             "B cell", "NK cell", "Langerhans cell", "Melanoma","CD4","Tumour","Stroma")

iteration<-1
iteration1<-1

for (i in combi.files) {
  melpin<-str_remove(i, pattern = "CD4_cluster.csv")
  melpin<-word(melpin,-1,sep=fixed("/"))
  print(paste("Formatting",iteration,melpin))
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
                            "CD4",
                            "Tumour",
                            "Stroma")
  
  intensity_columns_interest <- dye_columns_interest
  
  markers <- c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "DAPI", 
               "B cell", "NK cell", "Langerhans cell", "Melanoma","CD4","Tumour","Stroma")
  
  
  formatted_image <- format_halo_to_spe(path= i,
                                        markers=markers,
                                        dye_columns_interest=dye_columns_interest,
                                        intensity_columns_interest=intensity_columns_interest)
  print(unique(formatted_image$Phenotype))
  
  #remove double positive B cell/NK cell
  formatted_image<-select_celltypes(formatted_image,celltypes = c(tumour_phenotypes,stroma_phenotypes,melanoma_phenotypes), keep=TRUE)
  print(unique(formatted_image$Phenotype))
  
  #calculate cell proportions
  p_cells <- calculate_cell_proportions(formatted_image)
  print(p_cells)
  dev.off()
  
  #Entropy
  print("Entropy...")
  ent<-calculate_entropy(formatted_image,tumour_phenotypes)
  enttab<-data.table(melpin=melpin,entropy=ent)
  if (iteration==1){
    all.entropy<-enttab
  } else {
    all.entropy<-rbind(all.entropy,enttab)
  }
  
  avgmindist<-average_minimum_distance(formatted_image)*4.5
  print(avgmindist)
  
  disttab<-data.table(Melpin=melpin,Radius=avgmindist)
  if (iteration==1){
    all.disttab<-disttab
  } else {
    all.disttab<-rbind(all.disttab,disttab)}
  
  print("Clustering...")
  
  
  ###NEIGHBOURHOODS
  
  ##DBscan
  #Calculate & plot
  pdf(file=paste0(melpin,"ndmap.pdf"),height=15,width=15)
  neighbourd<-identify_neighborhoods(formatted_image, method="dbscan", cell_types_of_interest = neighbour_phenotypes,
                                     radius = avgmindist, feature_colname="Phenotype",no_pheno="OTHER")
  dev.off()
  print("Dbscan")
  #Heatmap
  nd_vis <- composition_of_neighborhoods(neighbourd, feature_colname="Phenotype")
  nd_vis <- nd_vis[nd_vis$Total_number_of_cells >=25,]
  pdf(file=paste0(melpin,"ndheatmap.pdf"), height=10, width=10,onefile=T)
  p =  plot_composition_heatmap(nd_vis, feature_colname="Phenotype")
  print(p)
  dev.off()
  nd_vis<-cbind(nd_vis,melpin)
  
  #Label clusters as tumour/stroma/marginal & possible TLS
  print("Labelling clusters...")
  #No melanoma cells in neighbourhoods - label as intratumoural/stroma/marginal clusters
  nd_vis<-nd_vis %>% separate(Phenotype,c("Phenotype","Classification"),sep=",")
  
  for (j in unique(nd_vis$Neighborhood)){
    nnum<-as.numeric(unique(nd_vis$Total_number_of_cells[nd_vis$Neighborhood==j]))
    
    #Determine what pct tumour vs stroma & make separate table
    tumoursum<-sum(as.numeric(nd_vis$Number_of_cells[nd_vis$Classification == "Tumour" & nd_vis$Neighborhood==j]))
    tumourpct<-(tumoursum/nnum)*100
    if (length(tumourpct)==0){
      tumourpct<-0
    }
    nd_vis$tumourpct[nd_vis$Neighborhood==j]<-tumourpct
    stromasum<-sum(as.numeric(nd_vis$Number_of_cells[nd_vis$Classification == "Stroma" & nd_vis$Neighborhood==j]))
    stromapct<-(stromasum/nnum)*100
    if (length(stromapct)==0){
      stromapct<-0
    }
    nd_vis$stromapct[nd_vis$Neighborhood==j]<-stromapct
    
    #Label cluster by tumour%
    if(min(nd_vis$tumourpct[nd_vis$Neighborhood==j])>72){
      nd_vis$Type[nd_vis$Neighborhood==j]<-"Intratumoural"
    } else if (min(nd_vis$stromapct[nd_vis$Neighborhood==j])>72){
      nd_vis$Type[nd_vis$Neighborhood==j]<-"Stromal"
    } else {
      nd_vis$Type[nd_vis$Neighborhood==j]<-"Marginal"
    }
    type<-unique(nd_vis$Type[nd_vis$Neighborhood==j])
    
    #Remove tumour/stroma annotation
    for (k in immune_phenotypes){
      cellsum<-with(nd_vis, sum(as.numeric(Number_of_cells[Phenotype==k & Neighborhood==j])))
      pctsum<-with(nd_vis, sum(as.numeric(Percentage[Phenotype==k & Neighborhood==j])))
      if (iteration1 ==1){
        phenvec<-c(k,NA,j,cellsum,nnum,pctsum, melpin,tumourpct,stromapct, type,NA)
      }else{
        phenvec<-c(k,NA,j,cellsum,nnum,pctsum, melpin,tumourpct,stromapct, type,NA,NA,NA)
      }
      nd_vis[nrow(nd_vis) + 1,]<-phenvec
    }
    
    #Determine what pct B cells
    bcellsum<-as.numeric(nd_vis$Number_of_cells[nd_vis$Phenotype=="B cell"& nd_vis$Neighborhood==j & is.na(nd_vis$Classification)])
    bcellpct<-(bcellsum/nnum)*100
    if (length(bcellpct)==0){
      bcellpct<-0
    }
    nd_vis$bcellpct[nd_vis$Neighborhood==j]<-bcellpct
    
    if (bcellpct>=10){
      nd_vis$TLS[nd_vis$Neighborhood==j]<-"Yes"
    } else{
      nd_vis$TLS[nd_vis$Neighborhood==j]<-"No"
    }
    
    #Note unclustered cell neighbourhoods are always Cluster 1
    if (iteration1==1){
      nd_vis$Freecells[nd_vis$Neighborhood==j]<-"Yes"
    } else{
      nd_vis$Freecells[nd_vis$Neighborhood==j]<-"No"
    }
    
    #Remove neighbourhoods with >=80% CD4s
    if(any(nd_vis$Percentage[nd_vis$Phenotype=="CD4" & nd_vis$Neighborhood==j]>=80)==TRUE){
      nd_vis<-nd_vis[nd_vis$Neighborhood!=j,]
    }
    iteration1<-iteration1+1
  }
  
  #Prepare data for heatmap - composition of clusters at each site
  #Isolate table with phenotypes only
  nd_visp<-nd_vis[is.na(nd_vis$Classification),]
  fwrite(nd_visp,file=paste0(melpin,"nd_visp.csv"),quote=F,sep=",")
  
  if (iteration==1){
    all_vis<-nd_visp
  } else{
    all_vis<-rbind(all_vis,nd_visp)
  }
  
  iteration<-iteration+1
  iteration1<-1
}

#Export intratumoural entropy
fwrite(all.entropy,file="primarytumourentropy.csv",quote=F,sep=",")

#rename clusters so nomatch
all_vis$Classification<-NULL
all_vis<-all_vis[all_vis$Percentage!=0,]
all_vis<-all_vis[all_vis$Freecells=="No",]

all_vis_combi_plot<-unite(all_vis,col="Neighborhood",c("Neighborhood","melpin"),sep="_",na.rm=TRUE,remove=FALSE)

#Remove 100% CD4 clusters from data
all_vis_combi_plot$Percentage<-as.numeric(all_vis_combi_plot$Percentage)
cd4ex<-subset(all_vis_combi_plot,Phenotype=="CD4" & Percentage >=80)
cd4exn<-unique(cd4ex$Neighborhood)
for (h in cd4exn){
  all_vis_combi_plot<-all_vis_combi_plot[all_vis_combi_plot$Neighborhood!=h,]
}

metadat<-as.data.frame(read.csv("clinicaldata.csv"))
outcomedat<-metadat[,c(1,5)]
colnames(outcomedat)<-str_remove(colnames(outcomedat),pattern="anythingextra")
all_vis_combi_plot<-merge(all_vis_combi_plot,outcomedat,by.x="melpin",by.y="Melpin",)

fwrite(all_vis_combi_plot, file="primaryallvisd.csv", quote=F, sep=",")
fwrite(all.disttab,file="primaryalldisttab.csv",quote=F,sep=",")

##Plot cluster Heatmaps including outcome,TLS,tissue localisation,cluster size

#Make cluster size for each image
clustersize <- unique(data.frame(Neighborhood = all_vis_combi_plot$Neighborhood, 
                                 Total_cells = all_vis_combi_plot$Total_number_of_cells))
rownames(clustersize) <- clustersize$Neighborhood
clustersize$Neighborhood <- NULL

#get data to matrix - need 3 variables Phenotype neighbourhood percentage
hmapdat <- all_vis_combi_plot[, c("Phenotype", "Neighborhood", 
                                  "Percentage")]
hmapdat <- reshape2::dcast(hmapdat, paste("Phenotype", 
                                          "~", "Neighborhood"), value.var = "Percentage")
rownames(hmapdat) <- hmapdat[, "Phenotype"]
hmapdat[, "Phenotype"] <- NULL
hmapdat[is.na(hmapdat)] <- -1
hmapdat<-as.matrix(hmapdat)

#Create heatmap annotations - include TLS, tissue localisation, cluster size, outcome for each cluster

annot<-unique(data.frame(Neighborhood=all_vis_combi_plot$Neighborhood))

#Outcome
ann.out <- unique(data.frame(Neighborhood = all_vis_combi_plot$Neighborhood, 
                             Outcome = all_vis_combi_plot$Outcome.group))

#TLS
ann.tls <- unique(data.frame(Neighborhood = all_vis_combi_plot$Neighborhood, 
                             TLS = all_vis_combi_plot$TLS))

#Tissue localisation
ann.tiss <- unique(data.frame(Neighborhood = all_vis_combi_plot$Neighborhood, 
                              Tissue = all_vis_combi_plot$Type))

ann.dat<-merge(ann.out,ann.tls,by="Neighborhood")%>%merge(ann.tiss,by="Neighborhood")
rownames(ann.dat)<-ann.dat$Neighborhood
ann.dat$Neighborhood<-NULL
fwrite(ann.dat,"neighbourannotationsdata.csv", row.names = T, sep=",")


#PUT ALL ANNOTATIONS IN THIS STEP
hann <- HeatmapAnnotation(df=ann.dat,
                          col = list(clustersize,       
                                     TLS = c("Yes" = "darkgoldenrod2", "No" = "floralwhite"),
                                     Outcome =c("B" = "darkorchid1", "A"= "forestgreen"),
                                     Tissue = c("Intratumoural" = "gold3", "Marginal" = "lightpink2", "Stromal" = "navyblue")))

map_cols <- circlize::colorRamp2(seq(min(hmapdat), max(hmapdat), length = 3), c("white","red", "purple"))

hmapdatnum<-apply(hmapdat,2,as.numeric)
rownames(hmapdatnum)<-rownames(hmapdat)

#Make heatmap
tiff(filename="allheatmap.tif",height=1000,width=5000)
ComplexHeatmap::Heatmap(hmapdatnum, col=map_cols, name = " ", 
                        cluster_columns = TRUE, cluster_rows = TRUE, row_names_gp = gpar(fontsize = 24), column_dend_height = unit(20, "mm"), row_dend_width = unit(20, "mm"),
                        show_column_names=FALSE,top_annotation = hann, border = TRUE, column_split=2, heatmap_legend_param = list(title_gp = gpar(fontsize = 20, fontface = "bold"),
                                                                                                                                  legend_height = unit(30, "mm"),
                                                                                                                                  legend_width = unit(30, "mm")))
dev.off()


pdf(file="allheatmapkmean11.pdf",height=10,width=50)
ComplexHeatmap::Heatmap(hmapdatnum, col=map_cols, name = " ", 
                        cluster_columns = TRUE, cluster_rows = TRUE, show_column_names=FALSE,top_annotation = hann, border = TRUE, column_km=11, column_km_repeats = 100)
dev.off()

#Find out what columns are in each cluster clusters

rcl.list <- column_order(HM)  #Extract metaclusters - which columns in which metacluster (output is a list)

lapply(rcl.list, function(x) length(x))  #check/confirm size clusters
class(hmapdat)
sum(duplicated(rownames(hmapdat)))

#Map metaclusters back to cluster & phenotype

for (i in 1:length(rcl.list)){
  metacluster<-paste("Metacluster",i)
  metcols<-as.vector(unlist(rcl.list[[i]]))
  clu <- hmapdat[,metcols]
  clu<-cbind(clu,metacluster)
  colnames(clu)[ncol(clu)]<-"Metacluster"
  clu <- cbind(rownames(clu), data.frame(clu, row.names=NULL))
  colnames(clu)[1]<-"Phenotype"
  if (i==1){
    metaclusters<-clu
  }else{
    metaclusters <- merge(metaclusters, clu,by=c("Metacluster","Phenotype"),all=TRUE)
  }
}

write.table(metaclusters, file= "metaclusterids.csv", sep=",", row.names=F,quote=F)

#Create heatmap data & annotations for TLS only
tls.ann.dat<-ann.dat[ann.dat$TLS=="Yes",]
tls.ann.dat$TLS<-NULL
tls.cols<-as.vector(rownames(tls.ann.dat))

tls.dat<-hmapdat[,tls.cols]

tann <- HeatmapAnnotation(df=tls.ann.dat,
                          col = list(clustersize,
                                     Outcome =c("B" = "darkorchid1", "A"= "forestgreen"),
                                     Tissue = c("Intratumoural" = "gold3", "Marginal" = "lightpink2", "Stromal" = "navyblue")))

map_cols <- circlize::colorRamp2(seq(min(tls.dat), max(tls.dat), length = 3), c("white","red", "purple"))

tls.datnum<-apply(tls.dat,2,as.numeric)
rownames(tls.datnum)<-rownames(tls.dat)

pdf(file="tlsheatmap.pdf",height=5,width=20)
ComplexHeatmap::Heatmap(tls.datnum, col=map_cols, name = " ", 
                        cluster_columns = TRUE, cluster_rows = TRUE, show_column_names=FALSE,top_annotation = tann, border = TRUE)
dev.off()

#Also map metaclusters back to outcome

#metaclustersid.csv file edited in excel to match pheatmap format
metcomp<-as.data.frame(read.csv("metaclusteridspheatmap.csv"))
metcomp<-metcomp[,c(4,9:20)]
metcompavg<-as.data.frame(metcomp %>% group_by(Metacluster) %>%  summarise_each(funs(round(mean(.),2))))
row.names(metcompavg)<-paste0("NMC",metcompavg$Metacluster)
metcompavg$Metacluster<-NULL
metcompavg<-as.matrix(metcompavg)

pdf(file="metaclustercompheatmapraw.pdf",width=10, height=10,onefile=T)
pheatmap(metcompavg,cluster_cols = F)
dev.off()

