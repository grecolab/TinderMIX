#'
#' This function takes in input the result of the function create_contour and plot the dynamic dose responsive activation map of a specific gene
#'
#' @param contour_res is the result object from the create_contour function
#' @param geneName is the name of the gene
#' @return a ggplot object
#' @export
#'
plot_dynamic_dose_responsive_map = function(contour_res, geneName,activity_threshold,
																						BMD_response_threshold,mode,nTimeInt,nDoseInt,
																						timeLabels,doseLabels){

	immy = contour_res$RPGenes[[geneName]][[3]]
	coord = cbind(contour_res$RPGenes[[geneName]][[1]],contour_res$RPGenes[[geneName]][[2]])
	res2 = compute_BMD_IC50(immy = immy,coord = coord, geneName,
													activity_threshold = activity_threshold,
													BMD_response_threshold = BMD_response_threshold,
													mode = mode,
													nTimeInt = nTimeInt,nDoseInt=nDoseInt,
													timeLabels = timeLabels,
													doseLabels = doseLabels,toPlot = T, addLegend = T)
}