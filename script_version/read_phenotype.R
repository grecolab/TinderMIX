task <- input[[1]]

tryCatch({
	source(file.path(getwd(), "server/scripts/tindermix/task_functions.R"), local = T)
	
	myInputs <- task$inputs
	
	filename <- myInputs$filename$value

	pheno = read.delim(filename)

	pheno_data_file <- task$outputs$pheno_data$value
	saveRDS(object=pheno, file=pheno_data_file)
	
},
error = function(e){
	addError('error99', e$message)
	end()
},
finally = {
	return(task)
})

