#' read excel file as a list of dataframe
#'
#' @import readxl
#' @param filename is the path to the file
#' @param tibble boolean specifying if the content of each sheet should be read as tibble or dataframe
#' @return a data frame contained in the excel file
#' @export

read_excel_allsheets <- function(filename, tibble = FALSE) {
	sheets <- readxl::excel_sheets(filename)
	x <- lapply(sheets, function(X){
		y = readxl::read_excel(filename, sheet = X)
		if(nrow(y)==0){
			stop("empty sheet")
		}
		y[,1] = as.character(as.vector(y[[1]]))
		y
	})
	
	if(!tibble) x <- lapply(x, as.data.frame)
	names(x) <- sheets

	data = x[[1]]
	rownames(data) = data[,1]
	data = data[,-1]
	
	return(data)
	
}
