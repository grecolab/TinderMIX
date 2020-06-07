task <- input[[1]]

tryCatch({
	source(file.path(getwd(), "server/scripts/tindermix/task_functions.R"), local = T)
	
	myInputs <- task$inputs
	
	plotHeight <- myInputs$height$value
	plotWidth <- myInputs$width$value
	
	DDRGene_file <- myInputs$DDRGene$value
	DDRGene = readRDS(DDRGene_file) #load data from CREATE_COUNTOUR step
	
	cake_diagram_plot = TinderMIX::plot_cake_diagrams_time_dose_effect(DDRGene)
	
	cake_diagram_plot_file <- task$outputs$cake_diagram_plot_file$value
	
	pdf(file=cake_diagram_plot_file, height=plotHeight, width=plotWidth)
	print(cake_diagram_plot)
	dev.off()
	
},
error = function(e){
	addError('error99', e$message)
	end()
},
finally = {
	return(task)
})

