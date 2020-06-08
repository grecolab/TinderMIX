
#'
#' This function takes in input the result of the function run_all_BMD_IC50 and plot a 3x3 
#' heatmap with the number of dose responsive genes for each label
#'
#' @param res is the result object from the run_all_BMD_IC50 function
#' @param drugName is the name of the drug that will be used in the title 
#' @return a ggplot object
#' @export
#'

plot_number_genes_labels = function(res, drugName ){
  library(ggplot2)
  
  M = res$MMA
  M = M[,c(3,6,9,2,5,8,1,4,7)]
  
  matrice = matrix(colSums(abs(M)),3,3)
  rownames(matrice) = c("Sensitive", "Intermediate","Resilient")
  colnames(matrice) = c("Early","Middle","Late")
  
  library(reshape2)
  melted_cormat <- melt(matrice)
  head(melted_cormat)
  colnames(melted_cormat)[1:2] = c("Dose","Time")
  melted_cormat$Dose = factor(melted_cormat$Dose, levels = rownames(matrice))
  melted_cormat$Time = factor(melted_cormat$Time, levels = colnames(matrice))
  
  p2 = ggplot(data = melted_cormat, aes(x=Dose, y=Time, fill=value)) + 
    geom_tile()+
    geom_text(aes(label = value), color = "white", size = 10) +
    ggtitle(drugName) + 
    theme(
      plot.title = element_text(size=14, face="bold.italic"),
      axis.title.x = element_text(size=14, face="bold"),
      axis.title.y = element_text(size=14, face="bold"),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.position = "none"
    )
  p2 
}

# this function takes in input the result of the function run_all_BMD_IC50 and plot a 3x3 multiplot with the number of dose responsive genes fir tge 12 segment of time and dose interaction
# the letters d and t (independently if they are capital or small ) stand for dose and time
# +/- indicate if the gene fc is increasing or decreasing with respect of dose and time
# capital letters are used to indicate which between dose and time has a stronger effect


#'
#' this function takes in input the result of the function run_all_BMD_IC50 and plot a 3x3 multiplot with the number of dose responsive genes fir tge 12 segment of time and dose interaction
#' the letters d and t (independently if they are capital or small ) stand for dose and time
#' +/- indicate if the gene fc is increasing or decreasing with respect of dose and time
#' capital letters are used to indicate which between dose and time has a stronger effect
#'
#' @param res is the result object from the run_all_BMD_IC50 function
#' @return a ggplot object
#' @export
#'

plot_cake_diagrams_time_dose_effect=function(res){
  M = res$MMA
  listPlot = list()
  
  library(ggplot2)
  library(easyGgplot2)
  
  index = 1
  for(i in c(1,4,7,2,5,8,3,6,9)){
    gi = rownames(M)[which(M[,i]!= 0)]
    
    if(length(gi)==0){
      df <- data.frame()
      p = ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
    }else{
      
      if(length(gi)==1){
        mi = matrix(M[gi, c("Dose","Time","Comparison(1Dose,2Time,0Both)")],nrow = 1, ncol = 3)
        rownames(mi) = gi
        colnames(mi) = c("Dose","Time","Comparison(1Dose,2Time,0Both)")
      }else{
        mi = M[gi, c("Dose","Time","Comparison(1Dose,2Time,0Both)")]
      }
      
      mmi = apply(mi, 1, FUN = function(elem){
        if(elem[1]==1) dose_sign = "+" else dose_sign = "-"
        if(elem[2]==1) time_sign = "+" else time_sign = "-"
        if(elem[3]==1){
          dose_letter = "D"
          time_letter = "t"
        }
        if(elem[3]==2){
          dose_letter = "d"
          time_letter = "T"
        } 
        if(elem[3]==0){
          dose_letter = "D"
          time_letter = "T"
        }
        paste(dose_letter, dose_sign, time_letter, time_sign, sep="")
      })
      
      #levels = c("Dose_1_Time_1_Comb_1","Dose_1_Time_1_Comb_0","Dose_1_Time_1_Comb_2",
      #                              "Dose_-1_Time_1_Comb_2","Dose_-1_Time_1_Comb_0","Dose_-1_Time_1_Comb_1",
      #                              "Dose_-1_Time_-1_Comb_1","Dose_-1_Time_-1_Comb_0","Dose_-1_Time_-1_Comb_2",
      #                              "Dose_1_Time_-1_Comb_2","Dose_1_Time_-1_Comb_0","Dose_1_Time_-1_Comb_1")
      
      # levels = c( "D(+)T(+)W(D)", "D(+)T(+)W(B)", "D(+)T(+)W(T)", 
      #             "D(-)T(+)W(T)", "D(-)T(+)W(B)", "D(-)T(+)W(D)",
      #             "D(-)T(-)W(D)", "D(-)T(-)W(B)", "D(-)T(-)W(T)",
      #             "D(+)T(-)W(T)", "D(+)T(-)W(B)", "D(+)T(-)W(D)")
      # 
      
      # the letters d and t (independently if they are capital or small ) stand for dose and time
      # +/- indicate if the gene fc is increasing or decreasing with respect of dose and time
      # capital letters are used to indicate which between dose and time has a stronger effect
      levels = c("D+t+", "D+T+", "d+T+", "d-T+", "D-T+", "D-t+", "D-t-", "D-T-", "d-T-", "d+T-", "D+T-", "D+t-")
      
      
      mmi = factor(mmi, levels = levels)
      ti = table(mmi)
      
      quadrant_names = factor(names(ti), levels = levels)
      df = data.frame(variable = quadrant_names, value = as.numeric(ti))
      #dose_time_lab = unlist(lapply(df$variable, FUN = function(elem)substr(x = elem, start = 1, stop = 8)))
      doseTimeEffect = tolower(df$variable)
      df = cbind(df, doseTimeEffect)
      
      p = ggplot(df, aes(x = variable, y = value, group =doseTimeEffect, fill=doseTimeEffect)) +  geom_bar(stat = "identity") #+ geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25)
      p = p + scale_y_continuous(1:12) + coord_polar(start=(-pi/2), direction = -1) + labs(x = "", y = "") 
      p = p + theme(legend.position = "none") + theme_minimal()+ ggtitle(colnames(M)[i])
      p = p + geom_vline(xintercept = 3.5) + geom_vline(xintercept = 6.5) + geom_vline(xintercept = 9.5) + geom_vline(xintercept = 12.5)
      p = p + theme(axis.title.y=element_blank()) + theme(legend.position = "none")
      p
      
    }
    
    listPlot[[index]] = p
    index = index + 1
  }
  
  ggplot2.multiplot(listPlot[[1]], listPlot[[2]], listPlot[[3]], listPlot[[4]],
                    listPlot[[5]], listPlot[[6]], listPlot[[7]], listPlot[[8]],
                    listPlot[[9]], cols=3)
  return(listPlot)
}


#'
#' this function takes in input the pathways enriched and retur a radar chart for each one of the gene label category
#'
#' @param Enriched_list is the list of dataframe resulting from enrichment for each gene category
#' @param n is the max number of pathways to plot in each radar plot
#' @param caxislabels is a character vector for center axis labels, overwriting values specified in axistype option. If NULL, the values specified by axistype option are used. Default is NULL.
#' @param vlcex is the size of the labels
#' @param kegg_level is the level of the kegg hierarchy to be considered in the plotting
#' @param mar are the margin settings for the plot
#' @return a ggplot object
#' @export
#'

plot_kegg_radar_chart = function(Enriched_list,  n = 5, vlcex = 1.5, kegg_level = 1,mar = c(2,1,1,1)){
  library(fmsb)
  
  par(mfrow = c(3,3), mar = mar)
  
  for(i in c(7,8,9,4,5,6,1,2,3)){
    EP_all = Enriched_list[[i]]
    if(kegg_level==1)ti = table(EP_all$Level1)
    else ti = table(EP_all$Level2)
    
    ti = sort(ti,decreasing = T)
    
    
    if(length(ti)>n) ti = ti[1:n]
    
    
    iid = which(names(ti) == "Environmental Information Processing")
    if(length(iid)>0){
      names(ti)[iid] = "Environmental \nInformation Processing"
    }
    
    iid = which(names(ti) == "Genetic Information Processing")
    if(length(iid)>0){
      names(ti)[iid] = "Genetic \nInformation Processing"
    }
    
    iid = which(names(ti) == "Cellular Processes")
    if(length(iid)>0){
      names(ti)[iid] = "Cellular \nProcesses"
    }
    
    iid = which(names(ti) == "Organismal Systems")
    if(length(iid)>0){
      names(ti)[iid] = "Organismal \nSystems"
    }
    
    iid = which(names(ti) == "Drug Development")
    if(length(iid)>0){
      names(ti)[iid] = "Drug \nDevelopment"
    }
    
    # iid = which(names(ti) == "Human Diseases")
    # if(length(iid)>0){
    #   names(ti)[iid] = "Human \nDiseases"
    # }
    
    pathways_levels = factor(names(ti))
    
    # Create data: note in High school for Jonathan:
    data <- as.data.frame(matrix( as.integer(ti) , ncol=length(ti)))
    colnames(data) <-pathways_levels
    
    # To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
    data <- rbind(rep(max(ti)+1,length(ti)) , rep(min(ti),length(ti)) , data)
    #data = apply(data,2,as.integer)
    
    # Custom the radarChart !
    radarchart( as.data.frame(data)  , axistype=1,
                
                #custom polygon
                pcol=rgb(0.2,0.5,0.5,0.9) , pfcol=rgb(0.2,0.5,0.5,0.5) , plwd=4 , 
                
                #custom the grid
                cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(from = 0, to = max(ti)+1, length.out = 5), cglwd=0.8,
                
                #custom labels
                vlcex=vlcex, title = names(Enriched_list)[i]
    )
    
  }
}
