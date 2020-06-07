task <- input[[1]]

tryCatch({
	source(file.path(getwd(), "server/scripts/tindermix/task_functions.R"), local = T)
	
	myInputs <- task$inputs
	
	filename <- myInputs$filename$value

	exp_data = read.delim(filename)
	
	exp_file <- task$outputs$exp_data$value
	saveRDS(object=exp_data, file=exp_file)
	
},
error = function(e){
	addError('error99', e$message)
	end()
},
finally = {
	return(task)
})

