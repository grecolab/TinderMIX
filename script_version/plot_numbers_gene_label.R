task <- input[[1]]

tryCatch({
	source(file.path(getwd(), "server/scripts/tindermix/task_functions.R"), local = T)
	
	myInputs <- task$inputs
	
	plotHeight <- myInputs$height$value
	plotWidth <- myInputs$width$value
	
	drugName <- myInputs$drugName$value
	
	DDRGene_file <- myInputs$DDRGene$value
	DDRGene = readRDS(DDRGene_file) #load data from CREATE_COUNTOUR step
	
	gene_plot = TinderMIX::plot_number_genes_labels(DDRGene,drugName = drugName)
	
	nGenePlotFile <- task$outputs$nGenePlotFile$value
	
	pdf(file=nGenePlotFile, height=plotHeight, width=plotWidth)
	print(gene_plot)
	dev.off()
	
},
error = function(e){
	addError('error99', e$message)
	end()
},
finally = {
	return(task)
})

