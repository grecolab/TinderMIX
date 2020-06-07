suppressWarnings(suppressMessages(library(TinderMIX)))


# addoutput function from taskmanager
addOutput = function(name, type, value) {
	
	task$outputs[[name]] <<- list("type" = type, "value" = value)
}

# adderror function from taskmanager
addError = function(id, message) {
	
	error <- list(
		"id" = id,
		"message" = message
	)
	
	index <- length(task$errors) + 1
	task$errors[[index]] <<- error
}