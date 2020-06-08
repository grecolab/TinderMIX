task <- input[[1]]

tryCatch({
	source(file.path(getwd(), "server/scripts/tindermix/task_functions.R"), local = T)
	
	myInputs <- task$inputs
	
	fc_data_file <- myInputs$fc_data$value
	fc_data = readRDS(fc_data_file)
	
	dose_index <- myInputs$doseIDCol$value
	time_point_index <- myInputs$timeIDCol$value
	
	#parameter for fitting
	gridSize = myInputs$gridSize$value
	pvalFitting = myInputs$pvalFitting$value
	pvalFitting.adj.method = myInputs$pvalFitting.adj.method$value
	modelSelection = myInputs$modelSelection$value
	logScale = myInputs$logScale$value
	
	exp_data = fc_data$fc_data
	pheno_data = fc_data$pdata
	
	contour_res = suppressMessages(TinderMIX::create_contour(exp_data, pheno_data, 
																													 responsive_genes = rownames(exp_data)[1:100],
																													 dose_index = dose_index,
																													 time_point_index = time_point_index ,
																													 gridSize = gridSize,
																													 pvalFitting = pvalFitting,
																													 pvalFitting.adj.method = pvalFitting.adj.method,
																													 logScale = logScale, modelSelection = modelSelection))	
	contour_res_file <- task$outputs$contour_res$value
	saveRDS(object=contour_res, file=contour_res_file)
	
},
error = function(e){
	addError('error99', e$message)
	end()
},
finally = {
	return(task)
})

