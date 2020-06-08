task <- input[[1]]

tryCatch({
	source(file.path(getwd(), "server/scripts/tindermix/task_functions.R"), local = T)
	myInputs <- task$inputs
	
	contour_resFile <- myInputs$contour_res$value
	contour_res = readRDS(contour_resFile) #load data from CREATE_COUNTOUR step
	
	#parameter for fitting
	activity_threshold = myInputs$activity_threshold$value
	BMD_response_threshold = myInputs$BMD_response_threshold$value
	mode = myInputs$mode$value
	nDoseInt = myInputs$nDoseInt$value
	nTimeInt = myInputs$nTimeInt$value
	doseLabels = myInputs$doseLabels$value
	timeLabels = myInputs$timeLabels$value
	
	#saving parameters
	tosave = myInputs$tosave$value
	addLegend = myInputs$addLegend$value
	path = myInputs$path$value
	toPlot = myInputs$toPlot$value
	
	DDRGene = TinderMIX::run_all_BMD_IC50(contour_res = contour_res,
																				activity_threshold = activity_threshold,  
																				BMD_response_threshold = BMD_response_threshold, 
																				mode=mode, 
																				nDoseInt=nDoseInt, nTimeInt=nTimeInt, 
																				doseLabels = doseLabels, timeLabels = timeLabels,
																				tosave=FALSE, addLegend = FALSE, path = ".",
																				relGenes = contour_res$ggenes, toPlot = FALSE)
	
	DDRGene_res_file <- task$outputs$DDRGene$value
	saveRDS(object=DDRGene, file=DDRGene_res_file)
	
},
error = function(e){
	addError('error99', e$message)
	end()
},
finally = {
	return(task)
})

