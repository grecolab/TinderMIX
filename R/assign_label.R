#'
#' This function assigns a label to contour map
#'
#' @param map matrix containing the z-map for a specific gene or cluster prototype
#' @return a list with 9x9 matrices specifying if the gene is active at low, mid or high time points and dose levels
#'
#' @examples
#' data("FC_WY14643")
#' exp_data = fc_data
#' pheno_data = pdata
#'   PvalMat = compute_anova_dose_time(exp_data, pheno_data,dose_index = 2,time_point_index = 3)
#'   ItemsList = build_items_list(PvalMat)
#'   responsive_genes = unique(c(unlist(ItemsList$Dose),unlist(ItemsList$Time),unlist(ItemsList$`Dose:Time:DoseTime`),unlist(ItemsList$`Dose:Time`)))
#'   contour_res = create_contour(exp_data, pheno_data, responsive_genes,dose_index = 2,time_point_index =3 ,gridSize = 50)
#'   hls_res = hls_genes_clustering(contour_res$GenesMap,  nClust = c(5,10,15,20,25), method="pearson", hls.method = "ward")
#'   clpr = create_prototypes(clust_res = hls_res,summaryMat = hls_res$summaryMat,contour_res )
#'   create_tic_tac_toe(map = t(clpr$meanXYZ[[1]][[3]]))
#' @export
#'
create_tic_tac_toe = function(map = t(clpr$meanXYZ[[1]][[3]])){
  which(map == min(map),arr.ind = T)
  which(map == max(map),arr.ind = T)
  
  ranges = cut(seq(5, 20, length.out=50),3)
  
  low = which(ranges==levels(ranges)[1])
  middle = which(ranges==levels(ranges)[2])
  high = which(ranges==levels(ranges)[3])
  ids = list(low,middle,high)
  
  #ttt is a 3x3 matrix, with low-mid-high on the column and high-mid-low on the rows.
  ttt = matrix(0, 3,3)
  rownames(ttt) = c("High","Mid","Low")
  colnames(ttt) = c("Low","Mid","High")
  for(i in 1:3){
    for(j in 1:3){
      ttt[(4-i),j] = mean(mean(map[ids[[i]],ids[[j]]]))
    }
  }
  
  qq = quantile(ttt)
  if(abs(qq[5]) > abs(qq[1])){
    print("ciao")
    ttt_lab = ttt>qq[4]
  }else{
    print("mondo")
    ttt_lab= ttt<qq[2]
  }
  
  return(list(tic_tac_toe = ttt, ttt_label = ttt_lab))
  
}