#' This function create a table with the information on the dynamic-dose-dependent genes
#'
#' @import readxl
#' @param DDRGene is the results of the run_all_BMD_IC50 function
#' @param contour_res is the result of the create_contour function
#' @param nTimeInt number of time points
#' @param nDoseInt number of dose levels
#' @param biomart_dataset is a string specifying the dataset to use in the useEnsembl function. e.g rnorvegicus_gene_ensembl
#' @return a data frame 
#' @export
create_gene_table = function(DDRGene,contour_res, nTimeInt,nDoseInt,biomart_dataset = "rnorvegicus_gene_ensembl"){
	StatMat = cbind(DDRGene$MMA, contour_res$Stats[rownames(DDRGene$MMA),])
	gi = rownames(StatMat)
	
	print("converting ensembl")
	ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset=biomart_dataset)
	
	genedesc <- getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = gi, mart =ensembl)
	toRem = which(genedesc$external_gene_name %in% names(table(genedesc$external_gene_name)[table(genedesc$external_gene_name)>1]))
	toRem = toRem[seq_along(toRem) %%2 == 0]
	if(length(toRem)>0)genedesc = genedesc[-toRem,]
	rownames(genedesc) = genedesc$external_gene_name
	
	labels = sapply(gi, FUN = function(gg)colnames(StatMat)[which(StatMat[gg,1:(nDoseInt*nTimeInt)]!=0)])
	splitted_labels = do.call(rbind, strsplit(x = labels,split = "-"))
	gene_sign = sapply(gi, FUN = function(gg)StatMat[gg,which(StatMat[gg,1:(nDoseInt*nTimeInt)]!=0)])
	
	XX = cbind(gi, labels, splitted_labels, gene_sign)
	XX  = cbind(XX, StatMat[rownames(XX), ((nDoseInt*nTimeInt)+1):20])
	
	mmi = apply(XX[,6:8], 1, FUN = function(elem){
		if(elem[2]==1) dose_sign = "+" else dose_sign = "-"
		if(elem[1]==1) time_sign = "+" else time_sign = "-"
		if(elem[3]==1){
			dose_letter = "D"
			time_letter = "t"
		}
		if(elem[3]==2){
			dose_letter = "d"
			time_letter = "T"
		} 
		if(elem[3]==0){
			dose_letter = "D"
			time_letter = "T"
		}
		paste(dose_letter, dose_sign, time_letter, time_sign, sep="")
	})
	
	XX = cbind(mmi,genedesc[rownames(XX),2],XX)
	return(XX)
}