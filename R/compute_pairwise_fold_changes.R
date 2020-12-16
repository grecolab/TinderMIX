#'
#' This function starts from a phenodata and gene expression data matrix and compute all the possible pairwise foldchange values
#'
#' @param exp_data is the log2 expression matrix with genes on the rows and samples on the columns
#' @param pheno_data is a dataframe with phenodata informations. Samples are on the rows. The columns should include the dose and time point information. Doses of controls need to be indicated as 0
#' @param dose_index index of the column containiing the dose
#' @param time_index index of the column containing the time
#' @return a list containing two new matrices
#' \item{fc_data}{a matrix with all the possible fold_changes}
#' \item{pdata}{the new phenodata table}
#' @export
#'
compute_fc = function(exp_data, pheno_data, dose_index, time_index){
	
	dose = sort(unique(pheno_data[,dose_index]))
	dose = dose[dose>0]
	
	time= sort(unique(pheno_data[,time_index]))
	
	fc_data = c()
	pdata = c()
	
	for (di in dose){
		print(di)
		for (tj in time){
			print(tj)
			idx = which(pheno_data[,time_index] == tj & pheno_data[,dose_index] == di)
			ctrl = which(pheno_data[,time_index] == tj & pheno_data[,dose_index] == 0)
			
			if(length(idx)>0 & length(ctrl)>0){
				comb = expand.grid(idx, ctrl) # contains all combinations of samples and corresponding controls
				
				fm = c()
				for (k in 1:nrow(comb)){
					new_fc = exp_data[comb[k,1]]-exp_data[comb[k,2]]
					colnames(new_fc) = paste(colnames(exp_data)[comb[k,1]], colnames(exp_data)[comb[k,2]], sep =".")
					
					fm = cbind(fm, as.matrix(new_fc))
					
					pdata = rbind(pdata, c(names(new_fc), di, tj))
				}
				
				fc_data = cbind(fc_data, fm)
			}
			
			
		}
	}
	
	colnames(pdata) = c("SampleID", "Dose", "Time")
	rownames(pdata) = pdata[,1]
	
	pdata = as.data.frame(pdata)
	pdata$Dose = as.numeric(as.vector(pdata$Dose))
	pdata$Time = as.numeric(as.vector(pdata$Time))
	
	return(list(fc_data=fc_data,pdata=pdata))
}


#save(fc_data, pdata, file = "data/FC_WY14643.RData")
