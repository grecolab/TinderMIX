task <- input[[1]]

tryCatch({
	source(file.path(getwd(), "server/scripts/tindermix/task_functions.R"), local = T)
	
	myInputs <- task$inputs
	
	plotHeight <- myInputs$height$value
	plotWidth <- myInputs$width$value
	
	contour_res_file <- myInputs$contour_res$value
	contour_res = readRDS(contour_res_file) 
	
	DDRGene_file <- myInputs$DDRGene$value
	DDRGene = readRDS(DDRGene_file) #load data from CREATE_COUNTOUR step
	
	biomart_dataset <- myInputs$biomart_dataset$value
	
	res = TinderMIX::create_gene_table(DDRGene,contour_res, biomart_dataset = biomart_dataset)
	
	gene_table_file <- task$outputs$gene_table_file$value
	
	write.table(res, file = gene_table_file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
	
},
error = function(e){
	addError('error99', e$message)
	end()
},
finally = {
	return(task)
})

