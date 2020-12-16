model_fitting_and_validation = function(Exp, pvalFitting= 0.05, modelSelection=1:3){
  
  #TODO: if we decide to add different models, e.g. hill we can ask the user which method they want to use (e.g. with an extra parameter modelSelectionStrategy = c("anova","bic"))
  # or we can still choose the optimal across the polynomials with the anova strategy and then compare this optimal one against the other by comparing their akaike criterion values
  
  # fitting linear, poly2 and poly3 models
  ModList = list()
  null_model <- stats::lm(Exp ~ 1, data = Exp)
  lin_model <- stats::lm(Exp ~ Dose + Time,data = Exp)
  poly2_model <- stats::lm(Exp ~ I(Dose^2) + I(Time^2) + I(Dose * Time) + Dose + Time, data = Exp)
  poly3_model <- stats::lm(Exp ~ I(Dose^3) + I(Dose^2 * Time) + I(Dose * Time^2) + I(Time^3) + I(Dose^2) + I(Time^2) + I(Dose * Time) + Dose + Time  ,data = Exp)
  
  ModList[["linear"]] = lin_model
  ModList[["poly2"]] = poly2_model
  ModList[["poly3"]] = poly3_model
  
  #store fitted models in a list without names for do.call on anova
  MML = list(lin_model, poly2_model, poly3_model)
  
  # select only the ones in modelSelection parameter
  ModList = ModList[modelSelection]
  MML = MML[modelSelection]
  
  # computing model statistics and p-value
  SingleStats = c()
  for(model in ModList){
    f <- summary(model)$fstatistic
    p <- stats::pf(f[1],f[2],f[3],lower.tail=F)
    adj.r.square = summary(model)$adj.r.squared
    RSS <- c(crossprod(model$residuals))
    MSE <- RSS / length(model$residuals)
    RMSE <- sqrt(MSE)
    BICmod = stats::BIC(model)
    AICmod = stats::AIC(model)
    
    SingleStats = rbind(SingleStats, c(p, adj.r.square,RMSE, BICmod, AICmod))
  }
  
  colnames(SingleStats) = c("PValue","Adj.R.Square","RMSE","BIC","AIC")
  rownames(SingleStats) = names(ModList)
  
  # create a list including the null model and compute the anova between the null model, the linear, the poly2 and poly3 models
  MMListN = c(list(null_model), MML)
  anovaTest = do.call(anova, MMListN)
  #take anova pvalues
  anovaPVal = anovaTest$`Pr(>F)`
  #assign to pvalues model names
  names(anovaPVal) = c("Null",names(ModList))
  #set to 1 the pvalue for the null model
  anovaPVal[is.na(anovaPVal)] = 1
  
  #identify the smallest pvalue
  optIdx = which.min(anovaPVal)
  # get the statistics for the model associated to the smallest pvalue
  SS = SingleStats[(optIdx-1),]
  
  # the model associated with the smallest pvalue is returned. We don't filter for sign pvalues here, since we will probably apply a correction after
  optModel = ModList[[(optIdx-1)]]
  SS["PValue"] = anovaPVal[optIdx] #N.B. the pvalue returned for the model is the one of the anova with respect to the previous model in the nested hierarchy
  SS = c(SS, names(ModList)[(optIdx-1)])
  names(SS)[length(SS)] = "OptMod"
  toRet = list(optModel = optModel, stats = SS, modelsStats = SingleStats,anovaTest=anovaTest)
  return(toRet)
  
  # # if the model with smallest pvalue is significative with respect to his previous model in the anova and with respect to the null model, than that gene is selected as the otpimal one
  # if(anovaPVal[optIdx] < pvalFitting & SS["PValue"]<pvalFitting){
  #   optModel = ModList[[(optIdx-1)]]
  #   SS["PValue"] = anovaPVal[optIdx] #N.B. the pvalue returned for the model is the one of the anova with respect to the previous model in the nested hierarchy
  #   SS = c(SS, names(ModList)[(optIdx-1)])
  #   names(SS)[length(SS)] = "OptMod"
  #   toRet = list(optModel = optModel, stats = SS)
  #   return(toRet)
  # }else{
  #   #otherwise the gene is considered not fitted and the function return an empy list
  #   return(list(c(),SS*NA))
  # }
}


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
#' @param logScale boolean specifying if the fitting is performed by using the dose and time in log or linear scale
#' @param modelSelection is a vector of indices specifying which model to fit. 1:linear 2: poly2, 3: poly3
#' @return a list with list with estimated contour objects, 3D fitted objects, fitting statistics and feature values for time and dose
#' \item{GenesMap}{a matrix with the z-maps computed for each gene}
#' \item{RPGenes}{a list with the 3D fitted objects}
#' \item{Statis}{a matrix witht the fitting statistics: PValue,Adj.R.Square,RMSE}
#' \item{DFList}{a list with the data used for the fitting}
#' \item{ModList}{a list with the fitted models}
#' @export
#'

#modelType string showing the model to fit. Possible options: poly2 and loess 

create_contour = function(exp_data, pheno_data, responsive_genes,dose_index, time_point_index, gridSize = 50, 
                          pvalFitting.adj.method = "fdr",pvalFitting=0.05, logScale = FALSE,
                          modelSelection = c(1,2)){  #models is a vector of indices specifying which model to fit. 1:linear 2: poly2, 3: poly3
  
  GenesMap = matrix(NA, ncol = length(responsive_genes), nrow = (gridSize*gridSize) )
  colnames(GenesMap) = responsive_genes
  
  RPGenes = list()
  DFList = list()
  ModList = list()
  
  Stats = matrix(0, nrow = length(responsive_genes), ncol = 6)
  rownames(Stats) = responsive_genes
  colnames(Stats) = c("PValue","Adj.R.Square","RMSE","BIC","AIC", "OptMod")
  
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
    
    if(logScale){
      Exp$Dose = log(Exp$Dose)
      Exp$Time = log(Exp$Time)
      #fit three different models and store the results in a list
      
      fitted_models = model_fitting_and_validation(Exp, modelSelection=modelSelection)
      
      optModel = fitted_models$optModel
      Stats[index,] = fitted_models$stats
  
      x1 <-range(as.numeric(as.vector(exp(Exp$Dose))))
      x1 <- seq(x1[1], x1[2], length.out=gridSize)
      y1 <- range(as.numeric(as.vector(exp(Exp$Time))))
      y1 <- seq(y1[1], y1[2], length.out=gridSize)
      
      x <-range(as.numeric(as.vector(Exp$Dose)))
      x <- seq(x[1], x[2], length.out=gridSize)
      y <- range(as.numeric(as.vector(Exp$Time)))
      y <- seq(y[1], y[2], length.out=gridSize)
      
      z <- outer(x,y, function(Dose,Time) stats::predict(optModel, data.frame(Dose,Time)))
      List3d=list(x,y,z)
      #plot3d(toPlot = List3d,logScale = TRUE, DF = Exp)
      
    }else{
      fitted_models = model_fitting_and_validation(Exp, modelSelection=modelSelection)
      optModel = fitted_models$optModel
      Stats[index,] = fitted_models$stats
      
      x <-range(as.numeric(as.vector(Exp$Dose)))
      x <- seq(x[1], x[2], length.out=gridSize)
      y <- range(as.numeric(as.vector(Exp$Time)))
      y <- seq(y[1], y[2], length.out=gridSize)
      
      #Z is a matrix dose X time Z[1,1] = first time and first dose
      z <- outer(x,y, function(Dose,Time) stats::predict(optModel, data.frame(Dose,Time)))
      List3d=list(x,y,z)
      #plot3d(toPlot = List3d,logScale = FALSE, DF = Exp)
      
    }
    
    #GenesMap = cbind(GenesMap,as.vector(z)) # I want to compute the distance by the genes based on the mapping z
    GenesMap[,g] =  as.vector(z)
    RPGenes[[g]] = List3d
    DFList[[g]] = Exp
    ModList[[g]] = optModel
    
    # p = plot3d(toPlot = list(x,y,z),DF = Exp)
    # p

    setTxtProgressBar(pb,index)
    index = index + 1
  }
  colnames(GenesMap) = responsive_genes#responsive_genes
  close(pb)
  
  SST = as.data.frame(Stats)
  SST[,1] = as.numeric(as.vector(SST[,1]))
  SST[,2] = as.numeric(as.vector(SST[,2]))
  SST[,3] = as.numeric(as.vector(SST[,3]))
  
  adj.pval = p.adjust(SST[,1],method = pvalFitting.adj.method)
  SST = cbind(adj.pval,SST)
  ggenes = rownames(SST)[SST[,1]<pvalFitting]
  Stats = SST
  return(list(GenesMap=GenesMap, RPGenes = RPGenes,Stats=Stats,DFList=DFList,ggenes=ggenes, ModList=ModList))
}

plot_contour_plot = function(immy, coord,geneName){
  # geneName = "Hspbp1" #"Acmsd
  # immy = contour_res$RPGenes[[geneName]][[3]]
  # 
  # immy2 = clpr[[1]][[3]]
  # coord = cbind(contour_res$RPGenes[[geneName]][[1]],contour_res$RPGenes[[geneName]][[2]])
  
  res2 = compute_BMD_IC50(immy,coord, geneName,
                          activity_threshold = activity_threshold,
                          BMD_resonse_threhold = BMD_resonse_threhold,
                          mode = mode,
                          nTimeInt = nTimeInt,nDoseInt=nDoseInt,
                          timeLabels = timeLabels,
                          doseLabels = doseLabels, toPlot = TRUE, tosave = FALSE)
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
plot3d = function(toPlot = list(x,y,z), DF, logScale = FALSE){
  
  if(logScale){
    p = plot_ly(x=~toPlot[[1]], y=~toPlot[[2]], z=~t(toPlot[[3]]),
                colors = c("#f5cb11", "#b31d83"),type="surface") %>%
      add_trace(data=DF, x=exp(DF$Dose), y=exp(DF$Time), z=DF$Exp, mode="markers", type="scatter3d",
                marker = list(opacity=0.7, symbol=105)) %>%
      layout(scene = list(
        aspectmode = "manual",
        aspectratio = list(x=1, y=1, z=1),
        xaxis = list(title = "Dose"),
        yaxis = list(title = "Time"),
        zaxis = list(title = "Exp")))
  }else{
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
  }
 
  
  # p = plot_ly(x=~exp(x), y=~exp(y), z=~t(z),
  #             colors = c("#f5cb11", "#b31d83"),type="surface") %>%
  #   add_trace(data=pd, x=exp(pd$Dose), y=exp(pd$Day), z=pd$Counts, mode="markers", type="scatter3d",
  #             marker = list(opacity=0.7, symbol=105)) %>%
  #   layout(scene = list(
  #     aspectmode = "manual",
  #     aspectratio = list(x=1, y=1, z=1),
  #     xaxis = list(title = "Dose"),
  #     yaxis = list(title = "Time"),
  #     zaxis = list(title = "Counts")))
  
  return(p)
}

