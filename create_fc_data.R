library(TinderMIX)
load("~/Google Drive (grecolab.fi@gmail.com)/shiny_BMDL/ttgate_dataset/rat_liver_invivo/common_genes.RData")
source("create_fc.R")

drug_list <- read.table("~/Google Drive (grecolab.fi@gmail.com)/shiny_BMDL/ttgate_dataset/rat_liver_invivo/drug_list.txt", sep="\t")

for(i in 1:length(drug_list)){
  drug_name = as.character(drug_list$V1)[i]
  print(drug_name)
  pd_path = paste("~/Google Drive (grecolab.fi@gmail.com)/shiny_BMDL/ttgate_dataset/rat_liver_invivo/datasets_old/", drug_name, "_pheno.txt",sep="")
  exp_path = paste("~/Google Drive (grecolab.fi@gmail.com)/shiny_BMDL/ttgate_dataset/rat_liver_invivo/datasets_old/", drug_name, "_exp.txt",sep="")
  out_path = "~/Google Drive (grecolab.fi@gmail.com)/shiny_BMDL/ttgate_dataset/rat_liver_invivo/dataset_FC/"
  
  pheno_data = read.delim(pd_path)
  rownames(pheno_data) = pheno_data$BARCODE
  colnames(pheno_data) = c("SampleID", "Dose","DoseLevel","DoseUnit","Time")
  exp_data = read.delim(exp_path)
  exp_data = exp_data[commonGenes,]
  
  dd = sort(unique(pheno_data$Dose))
  dataFC = compute_fc(exp_data, pheno_data, dose=dd[dd>0], time=sort(unique(pheno_data$Time)))
  
  dir.create(path = paste(out_path, drug_name, sep=""))
  save(dataFC, file = paste(out_path, drug_name, "/", drug_name, "FC_data.RData", sep=""))
  
  contour_res = suppressMessages(create_contour(dataFC$fc_data, dataFC$pdata, commonGenes,dose_index = 2,time_point_index =3 ,gridSize = 50))
  SST = contour_res$Stats
  GenesMap = contour_res$GenesMap
  save(contour_res, GenesMap, SST, file = paste(out_path, drug_name, "/", drug_name, "CMAP_data.RData", sep=""))
  
}
