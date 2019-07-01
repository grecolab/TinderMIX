# Create fold changes for each sample and corresponding controls (three FCs per sample) 

# load exp_data and pheno_data
#load("data/Exp_data_filtered_by_body_meth.RData")
# library(TinderMIX)
# data("WY14643")
# 
# dose = c(10,30,100)
# time = c(4, 8, 15, 29)
# 
# exp_data = WY14643$exp_data
# pheno_data = WY14643$pheno_data

compute_fc = function(exp_data, pheno_data, dose, time){
  fc_data = c()
  pdata = c()
  
  for (di in dose){
    print(di)
    for (tj in time){
      print(tj)
      idx = which(pheno_data$Time == tj & pheno_data$Dose == di)
      ctrl = which(pheno_data$Time == tj & pheno_data$Dose == 0)
      
      comb = expand.grid(idx, ctrl) # contains all combinations of samples and corresponding controls
      
      fm = c()
      for (k in 1:nrow(comb)){
        new_fc = exp_data[comb[k,1]]-exp_data[comb[k,2]]
        colnames(new_fc) = paste(colnames(exp_data)[comb[k,1]], colnames(exp_data)[comb[k,2]], sep =".")
        
        fm = cbind(fm, as.matrix(new_fc))
        
        pdata = rbind(pdata, c(colnames(new_fc), di, tj))
      }
      
    fc_data = cbind(fc_data, fm)
    }
  }
  
  colnames(pdata) = c("sample_id", "Dose", "Time")
  rownames(pdata) = pdata[,1]
  
  pdata = as.data.frame(pdata)
  pdata$Dose = as.integer(as.vector(pdata$Dose))
  pdata$Time = as.integer(as.vector(pdata$Time))
  
  return(list(fc_data=fc_data,pdata=pdata))
}


#save(fc_data, pdata, file = "data/FC_WY14643.RData")
