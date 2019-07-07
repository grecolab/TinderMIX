#'
#' This function fits a 3D regression model for every gene in the dataset and creates an N x N contour plot
#'
#' @importFrom reshape melt
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats predict lm pf
#' @param exp_data is the expression matrix with genes on the rows and samples on the columns
#' @param pheno_data is a dataframe with phenodata informations. Samples are on the rows. The columns should include the dose and time point information.
#' @param responsive_genes responsive_genes character vector with the genes statistically significant for the two-way anova
#' @param dose_index numeric value specifing the column of the phenodata table containing the doses
#' @param time_point_index numeric value specifing the column of the phenodata table containing the time points
#' @param gridSize numeric value specifing size of the z-grid
#' @return a list with list with estimated contour objects, 3D fitted objects, fitting statistics and feature values for time and dose
#' \item{GenesMap}{a matrix with the z-maps computed for each gene}
#' \item{RPGenes}{a list with the 3D fitted objects}
#' \item{Statis}{a matrix witht the fitting statistics: PValue,Adj.R.Square,RMSE}
#' \item{DFList}{a list with the data used for the fitting}
#'
#' @export
#'

#modelType string showing the model to fit. Possible options: poly2 and loess 

create_contour = function(exp_data, pheno_data, responsive_genes,dose_index, time_point_index, gridSize = 50){ #modelType = "poly2
  GenesMap = matrix(NA, ncol = length(responsive_genes), nrow = (gridSize*gridSize) )
  colnames(GenesMap) = responsive_genes
  RPGenes = list()
  DFList = list()

  Stats = matrix(0, nrow = length(responsive_genes), ncol = 3)
  rownames(Stats) = responsive_genes
  colnames(Stats) = c("PValue","Adj.R.Square","RMSE")
  index = 1
  pb = txtProgressBar(min = 1, max = length(responsive_genes), style = 3)
  for(g in responsive_genes){
    Exp = reshape::melt(exp_data[g,])
    Exp = cbind(rownames(Exp), Exp)
    colnames(Exp)[1] = "variable"
    Exp = cbind(Exp, pheno_data[as.character(Exp$variable),dose_index])
    Exp = cbind(Exp, pheno_data[as.character(Exp$variable),time_point_index])
    colnames(Exp) = c("Sample","Exp","Dose","Time")

    Exp = as.data.frame(Exp)
    Exp$Dose = as.numeric(as.vector(Exp$Dose))
    Exp$Time = as.numeric(as.vector(Exp$Time))
    model <- stats::lm(Exp ~ Dose * Time + I(Dose^2) + I(Time^2),data = Exp)
    
    # if(modelType == "poly2"){
    #   model <- stats::lm(Exp ~ Dose * Time + I(Dose^2) + I(Time^2),data = Exp)
    # }
    # if(modelType == "loess"){
    #   model = loess(Exp ~ Dose + Time,data = Exp)
    # }

    f <- summary(model)$fstatistic
    p <- stats::pf(f[1],f[2],f[3],lower.tail=F)
    adj.r.square = summary(model)$adj.r.squared
    RSS <- c(crossprod(model$residuals))
    MSE <- RSS / length(model$residuals)
    RMSE <- sqrt(MSE)
    Stats[index,] = c(p, adj.r.square,RMSE)

    x <-range(as.numeric(as.vector(Exp$Dose)))
    x <- seq(x[1], x[2], length.out=gridSize)
    y <- range(as.numeric(as.vector(Exp$Time)))
    y <- seq(y[1], y[2], length.out=gridSize)
    z <- outer(x,y, function(Dose,Time) stats::predict(model, data.frame(Dose,Time)))

    #GenesMap = cbind(GenesMap,as.vector(z)) # I want to compute the distance by the genes based on the mapping z
    GenesMap[,g] =  as.vector(z)
    RPGenes[[g]] = list(x,y,z)
    DFList[[g]] = Exp

    setTxtProgressBar(pb,index)
    index = index + 1
  }
  colnames(GenesMap) = responsive_genes#responsive_genes
  close(pb)

  return(list(GenesMap=GenesMap, RPGenes = RPGenes,Stats=Stats,DFList=DFList))
}

#'
#'This function plots the fitted 3d surface for the expression value of a gene
#'
#' @param toPlot is a list containing the predicted value for the x, y and z axis
#' @param DF is the data frame containing the information for the samples used in the fitting process
#' @return a plotly object
#'
#' @export
#'
#'
plot3d = function(toPlot = list(x,y,z), DF){
  p = plot_ly(x=~toPlot[[1]], y=~toPlot[[2]], z=~t(toPlot[[3]]),
              colors = c("#f5cb11", "#b31d83"),type="surface") %>%
    add_trace(data=DF, x=DF$Dose, y=DF$Time, z=DF$Exp, mode="markers", type="scatter3d",
              marker = list(opacity=0.7, symbol=105)) %>%
    layout(scene = list(
      aspectmode = "manual",
      aspectratio = list(x=1, y=1, z=1),
      xaxis = list(title = "Dose"),
      yaxis = list(title = "Time"),
      zaxis = list(title = "Exp")))
  return(p)
}

