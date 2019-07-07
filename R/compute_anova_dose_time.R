#'
#' This function computes a two way anova between dose and time for the expression value of every genes
#'
#' @importFrom stats lm aov p.adjust
#' @importFrom reshape melt
#'
#' @param exp_data is the expression matrix with genes on the rows and samples on the columns
#' @param pheno_data is a dataframe with phenodata informations. Samples are on the rows. The columns should include the dose and time point information.
#' @param dose_index numeric value specifing the column of the phenodata table containing the doses
#' @param time_point_index numeric value specifing the column of the phenodata table containing the time points
#' @return a matrix with pvalue associated to the dose, timepoint and the dose*timepoint effect
#' @export
#'
#'

compute_anova_dose_time  = function(exp_data, pheno_data, dose_index, time_point_index){
  PvalMat = c()
  anova_models = list()
  pb = txtProgressBar(min=1, max = nrow(exp_data), style = 3)
  for(i in 1:nrow(exp_data)){
    Exp = reshape::melt(exp_data[i,])
    Exp = cbind(rownames(Exp),Exp)
    colnames(Exp)[1] = "variable"
    Exp = cbind(Exp, pheno_data[as.character(Exp$variable),dose_index])
    Exp = cbind(Exp, pheno_data[as.character(Exp$variable),time_point_index])
    colnames(Exp) = c("Sample","Exp","Dose","Time")
    Exp$Dose = as.factor(Exp$Dose)
    Exp$Time = as.factor(Exp$Time)
    fit = aov(Exp ~ Dose * Time , data = Exp)
    sf = summary(fit)
    anova_models[[rownames(exp_data)[i]]] = fit
    PvalMat = rbind(PvalMat,sf[[1]][1:3,5])
    setTxtProgressBar(pb,i)
  }
  close(pb)

  colnames(PvalMat) = c("Dose","Time","DoseTime")
  rownames(PvalMat) = rownames(exp_data)
  PvalMatAdj = PvalMat
  PvalMatAdj[,1] = stats::p.adjust(PvalMatAdj[,1], method = "fdr")
  PvalMatAdj[,2] = stats::p.adjust(PvalMatAdj[,2], method = "fdr")
  PvalMatAdj[,3] = stats::p.adjust(PvalMatAdj[,3], method = "fdr")
  return(PvalMat)
}

#'
#' This function computes the venn diagram of the genes associated to time, dose or their interaction
#' @importFrom gplots venn
#'
#' @param PvalMat matrix with pvalue associated to the dose, timepoint and the dose*timepoint effect that is the output of the compute_anova_dose_time function
#' @param p.val.th is the threshold at which p.values are considered significant. Default = 0.01
#' @return a list containing the genes in each position of the venn diagram
#'
#'
#' @export
#'

build_items_list = function(PvalMat, p.val.th = 0.01){
  GL = list(Dose = rownames(PvalMat)[PvalMat[,1]<p.val.th],
            Time = rownames(PvalMat)[PvalMat[,2]<p.val.th],
            DoseTime = rownames(PvalMat)[PvalMat[,3]<p.val.th])
  ItemsList = gplots::venn(GL)
  ItemsList = attr(ItemsList,"intersections")
  return(ItemsList)
}
