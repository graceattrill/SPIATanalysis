library(data.table)
library(stringr)
library(spatstat)#for Kcross
library(Rphenograph) 
library(SPIAT)
setwd(primaryoutput) 

I.files<-list.files(path=primaryinput, pattern="tableI.csv", full.names=T, recursive=T, include.dirs=T)
I.phenotypes<-c("B cell", "NK cell", "Langerhans cell")

iteration<-1
for (i in I.files){
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
  
markers <- c("B cell", "NK cell", "Langerhans cell", "HLAABCpMelanoma", "HLAABCnMelanoma", "DAPI")
  
newdat<-str_replace(i, pattern = "tableI.csv", replacement = "I_SPIAT")

#upload file to R
formatted_image <- format_image_to_sce(format="HALO",
                                       path= i,
                                       markers=markers,
                                       dye_columns_interest=dye_columns_interest,
                                       intensity_columns_interest=intensity_columns_interest) 
print_column(formatted_image, column= "Phenotype")
formatted_image<-select_phenotypes(formatted_image, keep=FALSE, phenotypes = c("B cell,NK cell", "NK cell,Langerhans cell"))
print_column(formatted_image, column = "Phenotype")


#calculate cell proportions
p_cells <- calculate_cell_proportions(formatted_image, reference_celltypes=c("Total"))

#Calculate mixing score - Melanoma to immune
mixscoretab<-data.frame(phenotype=character(),mixingscore=numeric())
for(j in I.phenotypes){
  mixscorenum<-compute_mixing_score(formatted_image,reference_marker="HLAABCnMelanoma",target_marker=j,radius=20)
  looptab<-data.frame(phenotype=j, score=mixscorenum)
  mixscoretab<-rbind(mixscoretab,looptab)
}
if(iteration==1){
  all.mixing.scores<-mixscoretab
}
else{
  all.mixing.scores<-cbind(all.mixing.scores,mixscoretab[,2])
}


#Kcross
Ktab<-data.frame(phenotype=character(),Kscore=numeric())
for (k in I.phenotypes){
  df_cross<-bivariate_dependence(formatted_image, method = "Kcross", phenotypes = c(k, "HLAABCnMelanoma"), column = "Phenotype")
  Kcrossn<-AUC_of_cross_function(df_cross)
  Klooptab<-data.frame(phenotype=paste0(k, "Kcross"), Kscore = Kcrossn)
  Ktab<-rbind(Ktab, Klooptab)
  dev.off()
}

if(iteration == 1){
  all.Kcross<-Ktab
} else{
  all.Kcross<-cbind(all.Kcross,Ktab[,2])
}
  iteration<-iteration+1
}

colnames(all.mixing.scores)<-c("Phenotype", I.files)
colnames(all.Kcross) = c("Phenotype", c(I.files))

fwrite(all.mixing.scores, file = "SPIATprimarymixingscoresIHLAn.csv", quote=F)
fwrite(all.Kcross, file = "SPIATprimaryKcrossIHLAn.csv", quote=F)

