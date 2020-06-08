task <- input[[1]]

tryCatch({
	source(file.path(getwd(), "server/scripts/tindermix/task_functions.R"), local = T)
	
	myInputs <- task$inputs
	
	exp_data_file <- myInputs$exp_data$value
	pheno_data_file <- myInputs$pheno_data$value
	dose_index <- myInputs$doseIDCol$value
	time_index <- myInputs$timeIDCol$value
	
	exp_data = readRDS(exp_data_file)
	pheno_data = readRDS(pheno_data_file)
	
	fc_data = TinderMIX::compute_fc(exp_data = exp_data,pheno_data = pheno_data,dose_index = dose_index,time_index = time_index)
	
	fc_data_file <- task$outputs$fc_data$value
	saveRDS(object=fc_data, file=fc_data_file)
	
},
error = function(e){
	addError('error99', e$message)
	end()
},
finally = {
	return(task)
})

