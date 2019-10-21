#'
#' This function assigns a label to contour map
#'
#' @param map matrix containing the z-map for a specific gene or cluster prototype
#' @param BMD matrix containing the dose-response area
#' @param th a threshold to define the portion of dose-response area to be identified as labels for the gene.
#' @param nDoseInt number of dose related breaks in the gene label's table. default is 3
#' @param nTimeInt number of time related breaks in the gene label's table. default is 3
#' @param doseLabels vector of colnames (doses) for the gene label's table. default is  c("Sensitive","Intermediate","Resilient")
#' @param timeLabels vector of rownames (time points) for the gene label's table. default c("Late","Middle","Early")
#' @param mode is a character specifying when an area is called active. values can be "cumulative" or "presence". If presence, an area is called active if at least one of its pixel is on the BMD curve. If cumulative, the number of region needed to reach the th% of the cumulative of the number of pixel on the BMD curve is identified.
#' @param coord matrix with x and y coordinate. The first column contain the doses, while the second one the time points
#' @param nDoseInt number of dose related breaks in the gene label's table. default is 3
#' @param nTimeInt number of time related breaks in the gene label's table. default is 3
#' @param doseLabels vector of colnames (doses) for the gene label's table. default is  c("Sensitive","Intermediate","Resilient")
#' @param timeLabels vector of rownames (time points) for the gene label's table. default c("Late","Middle","Early")
#' @param myContour matrix with coordinate of bmd area border
#' @return a list with 9x9 matrices specifying if the gene is active at low, mid or high time points and dose levels
#'
#' @export
#'
label2DMap = function(map, BMD, coord, myContour, th=0.95, mode = "mix", nDoseInt=3, nTimeInt=3, doseLabels = c("Late","Middle","Early"), timeLabels = c("Sensitive","Intermediate","Resilient"), toplot = FALSE){
  
  mapSize = ncol(map)
  rangesDoses = cut(seq(5, 20, length.out=mapSize),nDoseInt)
  rangesTimes = cut(seq(5, 20, length.out=mapSize),nTimeInt)
  
  # print(rangesDoses)
  # print(rangesTimes)
  
  idsDoses = list()
  for(i in 1:nDoseInt){
    idsDoses[[i]] = which(rangesDoses==levels(rangesDoses)[i])
    if(toplot){
      print(exp(coord[idsDoses[[i]],1][length(idsDoses[[i]])]))
      graphics::abline(v = coord[idsDoses[[i]],1][length(idsDoses[[i]])], lty = 2)
    }
  }
  
  idsTimes = list()
  for(i in 1:nTimeInt){
    idsTimes[[i]] = which(rangesTimes==levels(rangesTimes)[i])
    if(toplot){
      print(exp(coord[idsTimes[[i]],2][length(idsTimes[[i]])]))
      graphics::abline(h = coord[idsTimes[[i]],2][length(idsTimes[[i]])], lty = 2)
    }
    
  }
  
  CM = matrix(0,mapSize,mapSize)
  for(i in 1:nrow(myContour)){
    CM[myContour[i,1], myContour[i,2]] = 1
  }
  
  
  #ttt is a 3x3 matrix, with low-mid-high on the column and high-mid-low on the rows.
  ttt = matrix(0, nrow = nTimeInt,ncol = nDoseInt)
  rownames(ttt) = timeLabels
  colnames(ttt) = doseLabels
  for(i in 1:nTimeInt){
    for(j in 1:nDoseInt){
      blockMat = matrix(0,mapSize,mapSize)
      blockMat[idsTimes[[i]],idsDoses[[j]]]=1
      ttt[i,j] = sum(blockMat * CM)/nrow(myContour)
    }
  }
  
  BMD[is.na(BMD)] = 0
  
  # percentage of coverage of each area by the BMD matrix
  ttt_area = matrix(0, nrow = nTimeInt,ncol = nDoseInt)
  rownames(ttt) = timeLabels
  colnames(ttt) = doseLabels
  for(i in 1:nTimeInt){
    for(j in 1:nDoseInt){
      ttt_area[i,j] = sum(BMD[idsTimes[[i]],idsDoses[[j]]])/(length(idsTimes[[i]]) * length(idsTimes[[j]]))
    }
  }
  
  if(mode == "cumulative"){
    ttt2 = ttt
    ttt_lab = ttt * 0
    cumulativa = 0
    while(cumulativa<th){
      idx = which(ttt2==max(ttt2), arr.ind = T)[1,]
      ttt_lab[idx[1], idx[2]] = TRUE
      cumulativa = cumulativa + ttt[idx[1], idx[2]]
      ttt2[idx[1], idx[2]]=0
    }
  }
  if(mode == "most_left"){ # identify the area that is closer to the (0,0) corner
    iid = which(myContour[,2] == min(myContour[,2]), arr.ind = T)
    myContour2 = matrix(0, nrow = length(iid), ncol = ncol(myContour))
    myContour2[,1] = myContour[iid,1]
    myContour2[,2] = myContour[iid,2]
    myContour = myContour2
    
    iid = which(myContour[,1] == max(myContour[,1]), arr.ind = T)[1]
    
    iid = myContour[iid,]
    
    ri = which(unlist(lapply(idsTimes,FUN = function(elem){
      iid[1] %in% elem
    })))
    ci = which(unlist(lapply(idsDoses,FUN = function(elem){
      iid[2] %in% elem
    })))
    ttt_lab = ttt * 0
    ttt_lab[ri, ci] = 1
  }
  if(mode == "presence"){
    ttt_lab = 0 * ttt
    ttt_lab[ttt>0] = 1
  }
  if(mode == "mix"){ # mix: takes into account of the number of points from the yellow label, the coverage of the area and the position of the area
    mm = matrix(c(7,8,9,4,5,6,1,2,3),3,3) / 9
    mm[ttt==0]=0
    ttt_sum = mm
    idx = which(ttt_sum == max(ttt_sum), arr.ind = T)
    ttt_lab = ttt_sum * 0
    ttt_lab[idx[1,1],idx[1,2]] = 1
  }
  
  
  # qq = quantile(ttt)
  # 
  # if(abs(qq[5]) >= abs(qq[1])){
  #   qq1 = quantile(ttt, th)
  #   
  #   if(qq1>0){
  #     ttt_lab = ttt>=qq1
  #   }else{
  #     ttt_lab = ttt>=quantile(ttt,1)
  #   }
  # }else{
  #   qq1 = quantile(ttt, 1-th)
  #   
  #   if(qq1>0){
  #     ttt_lab= ttt<=qq1
  #   }else{
  #     ttt_lab= ttt<=quantile(ttt,0)
  #   }
  # }
  # 
  # ttt = ttt
  
  return(list(tic_tac_toe = ttt, ttt_label = ttt_lab, ttt_area=ttt_area))
  
}

#   based on how many 1 I have in the tic-tac-toe blocks. the 1 are the BMD area.
#
#
# label2DMap = function(map, th=0.95, nDoseInt=3, nTimeInt=3, doseLabels = c("Low","Mid","High"), timeLabels = c("High","Mid","Low")){
#   mapSize = ncol(map)
#   rangesDoses = cut(seq(5, 20, length.out=mapSize),nDoseInt)
#   rangesTimes = cut(seq(5, 20, length.out=mapSize),nTimeInt)
#   
#   idsDoses = list()
#   for(i in 1:nDoseInt){
#     idsDoses[[i]] = which(rangesDoses==levels(rangesDoses)[i])
#   }
#   
#   idsTimes = list()
#   for(i in 1:nTimeInt){
#     idsTimes[[i]] = which(rangesTimes==levels(rangesTimes)[i])
#   }
#   
#   
#   #ttt is a 3x3 matrix, with low-mid-high on the column and high-mid-low on the rows.
#   ttt = matrix(0, nrow = nTimeInt,ncol = nDoseInt)
#   rownames(ttt) = timeLabels
#   colnames(ttt) = doseLabels
#   for(i in nTimeInt:1){
#     for(j in 1:nDoseInt){
#       ttt[i,j] = mean(mean(map[idsTimes[[i]],idsDoses[[j]]]))
#     }
#   }
#   
#   qq = quantile(ttt)
#   
#   if(abs(qq[5]) >= abs(qq[1])){
#     qq1 = quantile(ttt, th)
#     
#     if(qq1>0){
#       ttt_lab = ttt>=qq1
#     }else{
#       ttt_lab = ttt>=quantile(ttt,1)
#     }
#   }else{
#     qq1 = quantile(ttt, 1-th)
#     
#     if(qq1>0){
#       ttt_lab= ttt<=qq1
#     }else{
#       ttt_lab= ttt<=quantile(ttt,0)
#     }
#   }
#   
#   ttt = ttt
#   
#   return(list(tic_tac_toe = ttt, ttt_label = ttt_lab))
#   
# }


# it takes the first and last 1 in each row		
# it takes the first and last 1 in each column
bwtraceboundary= function(ternaryIMBMD){
  diffMat = 0 * ternaryIMBMD
  for(i in 1:nrow(ternaryIMBMD)){
    idx = which(ternaryIMBMD[i,]==1)
    diffMat[i, idx[1]] = 1
    diffMat[i, idx[length(idx)]] = 1
  }
  
  for(i in 1:ncol(ternaryIMBMD)){
    idx = which(ternaryIMBMD[,i]==1)
    diffMat[idx[1],i] = 1
    diffMat[idx[length(idx)],i] = 1
  }
  
  image_border = which(diffMat==1,arr.ind = 1)
  return(list(diffMat=diffMat,image_border=image_border))
}

rotate <- function(x) t(apply(x, 2, rev))
rad2deg <- function(rad) {(rad * 180) / (pi)}

#'
#' This function run the compute_BMD_IC50 for all genes and return a matrix with label associated to every gene
#'
#' @importFrom pracma gradient
#'
#' @param contour_res object resulting from the create_contour function
#' @param coord matrix with x and y coordinate. The first column contain the doses, while the second one the time points
#' @param geneName is a character string containing the gene name
#' @param activity_threshold threshold defining the responsive gene area. Eg. if the immy maps contains genes logFC, then an activity_threhdold = 0.58 means that the active area will be the one with an effect of 1.5 bigger or smaller than the controls
#' @param BMD_response_threshold a threshold to define the portion of dose-response area to be identified as labels for the gene.
#' @param nDoseInt number of dose related breaks in the gene label's table. default is 3
#' @param nTimeInt number of time related breaks in the gene label's table. default is 3
#' @param doseLabels vector of colnames (doses) for the gene label's table. default is  c("Sensitive","Intermediate","Resilient")
#' @param timeLabels vector of rownames (time points) for the gene label's table. default c("Late","Middle","Early")
#' @param tosave if true a png of the gene map is saved in path
#' @param path path of the folder where to save the gene map
#' @param relGenes vector of genes with signifincant pvalues from the fitting
#' @param mode is a character specifying when an area is called active. values can be "cumulative" or "presence". If presence, an area is called active if at least one of its pixel is on the BMD curve. If cumulative, the number of region needed to reach the th% of the cumulative of the number of pixel on the BMD curve is identified.
#' @return a list with two object: Mat is a matrix with genes on the rows and labels on the columns. GeneRes is a list of results from the compute_BMD_IC50 function, one for every gene
#'
#' @export
#'

run_all_BMD_IC50 = function(contour_res,activity_threshold = 0.1, BMD_response_threshold = 0.95, 
                            nDoseInt=3, nTimeInt=3, 
                            doseLabels = c("Late","Middle","Early"), timeLabels = c("Sensitive","Intermediate","Resilient"),
                            tosave=FALSE, toPlot = FALSE, addLegend = FALSE, path = ".",
                            relGenes, mode = "cumulative"){
  label_leg = c()
  for(i in doseLabels){
    for(j in timeLabels){
      label_leg = c(label_leg, paste(i,j,sep="-"))
    }
  }
  
  Map = list()
  goodGenes = c()
  GeneRes = list()
  geneDegree = c()
  geneDegreeLabel = c()
  
  pb = txtProgressBar(min = 1, max = length(contour_res$RPGenes), style = 3)
  for(i in 1:length(contour_res$RPGenes)){

    geneName = names(contour_res$RPGenes)[i]
    
    if((geneName %in% relGenes)==FALSE)
      next
    
    immy = contour_res$RPGenes[[geneName]][[3]]
    coord = cbind(contour_res$RPGenes[[geneName]][[1]],contour_res$RPGenes[[geneName]][[2]])
    tryCatch({
      res = compute_BMD_IC50(immy,coord, geneName,activity_threshold = activity_threshold,BMD_response_threshold=BMD_response_threshold,
                             nDoseInt=nDoseInt,nTimeInt=nTimeInt,doseLabels=doseLabels,timeLabels=timeLabels,tosave = tosave,path = path, toPlot = toPlot, addLegend = addLegend,mode=mode)
      GeneRes[[geneName]] = res
      if(is.null(res$verso)==FALSE){
        labels =  as.vector(res$label$ttt_label)
        Map[[geneName]] =  c(labels,res$verso)
        geneDegree = c(geneDegree,res$timedosedegree)
        geneDegreeLabel = c(geneDegreeLabel, res$degree_lab)
        goodGenes = c(goodGenes,geneName)
      }
    }, error = function(e) {
      print(NULL)
    })
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  names(geneDegree) = goodGenes
  names(geneDegreeLabel) = goodGenes
  
  Mat = do.call(rbind,Map)
  colnames(Mat) = c(label_leg,"Verso")
  verso = Mat[,"Verso"]
  Mat = Mat[,-ncol(Mat)]
  for(i in 1:nrow(Mat)){
    Mat[i,] = Mat[i,]* verso[i]
  }

  MMA = cbind(Mat, rowSums(abs(Mat)))
  
  labels = list()
  for(i in 1:nrow(Mat)){
    labels[[i]] = paste(colnames(Mat)[Mat[i,]!=0], collapse = " ")
  }
  

  return(list(Mat=Mat,MMA=MMA,GeneRes=GeneRes, geneDegree=geneDegree,geneDegreeLabel=geneDegreeLabel,labels=labels))
}

#'
#' This function identify the BMD area and the IC50 value in the time and dose maps 
#'
#' @importFrom pracma gradient
#' @importFrom pracma meshgrid
#' @importFrom raster contour
#'
#' @param immy z-maps of the fitted 3D model, with doses on the columns and time points on the rows
#' @param coord matrix with x and y coordinate. The first column contain the doses, while the second one the time points
#' @param geneName is a character string containing the gene name
#' @param activity_threshold threshold defining the responsive gene area. Eg. if the immy maps contains genes logFC, then an activity_threhdold = 0.58 means that the active area will be the one with an effect of 1.5 bigger or smaller than the controls
#' @param BMD_response_threshold a threshold to define the portion of dose-response area to be identified as labels for the gene.
#' @param nDoseInt number of dose related breaks in the gene label's table. default is 3
#' @param nTimeInt number of time related breaks in the gene label's table. default is 3
#' @param doseLabels vector of colnames (doses) for the gene label's table. default is  c("Sensitive","Intermediate","Resilient")
#' @param timeLabels vector of rownames (time points) for the gene label's table. default c("Late","Middle","Early")
#' @param toPlot it true the gene map is displayed
#' @param addLegend if true the legend will be added to the plot
#' @param tosave if true a png of the gene map is saved in path
#' @param path path of the folder where to save the gene map
#' @param mode is a character specifying when an area is called active. values can be "cumulative" or "presence". If presence, an area is called active if at least one of its pixel is on the BMD curve. If cumulative, the number of region needed to reach the th% of the cumulative of the number of pixel on the BMD curve is identified.
#' @return an object of class TinderMIX containing the fitted BMD object, the IC50 value. The function plot the map showing the responsive region.
#'
#' @export
#'

compute_BMD_IC50 = function(immy,coord, geneName,activity_threshold = 0.1, BMD_response_threshold = 0.95, nDoseInt=3, nTimeInt=3, 
                            doseLabels = c("Late","Middle","Early"), timeLabels = c("Sensitive","Intermediate","Resilient"), toPlot = TRUE,addLegend = TRUE, tosave=FALSE, path = ".", mode = "cumulative"){
  gridSize = nrow(immy)
  immy = rotate(rotate(rotate(immy)))
  activity_threshold = log2((10 + (10*activity_threshold))/10)
  
  if(max(abs(immy))<activity_threshold){
    print("No DE gene")
    
    if(toPlot){
      if(tosave){
        grDevices::png(filename = paste(path, geneName, ".png", sep=""))
      }
      graphics::image(coord[,1], coord[,2], rotate(immy), col = c("darkblue","darkgreen","brown"), xlab = "Dose",ylab = "Time", main = geneName)
      raster::contour(coord[,1], coord[,2], rotate(immy), add = TRUE,labcex = 1.3, col = "white")
    }
    
    #return(1)
    ans = list()
    ans$immy = immy
    ans$activity_threshold = activity_threshold
    ans$gradient = NULL
    ans$binaryIMBMD = NULL
    ans$verso = NULL
    ans$label = NULL
    ans$tracedarea = NULL
    ans$bmd = NULL
    ans$IC50 = NULL
    ans$geneName = geneName
    
    class(ans) = 'TinderMIX'
    return(ans)
  }
  
  binaryIMBMD = immy
  
  nonelig = which(abs(binaryIMBMD)<activity_threshold)		  
  elig = which(abs(binaryIMBMD)>=activity_threshold)
  
  binaryIMBMD[nonelig]=0; # non-eligible region
  binaryIMBMD[elig]=1; # eligible region
  
  gg = pracma::gradient(immy,coord[,1], coord[,2])
  gx = gg$X
  gy = gg$Y
  X <- pracma::meshgrid(coord[,1],coord[,1])$X
  Y <- pracma::meshgrid(coord[,2],coord[,2])$Y
  
  #IDENTIFY BMD REGION
  if(sum(abs(immy) < activity_threshold) == 0){  # % if non-eligible region in empty
    ternaryIMBMD = 0 * immy;
    ternaryIMBMD[abs(immy) >= activity_threshold & gx>=0]=1; # GREEN: gradient component in dose direction positive, i.e. dose increases
    ternaryIMBMD[abs(immy) >= activity_threshold & gx<0]=2; #% YELLOW: gradient component in dose direction negative, i.e. dose decreases
    #image(ternaryIMBMD)
  }else{
    ternaryIMBMD = 0 * immy;
    ternaryIMBMD[abs(immy)>=activity_threshold & gx>=0]=1; # GREEN: gradient component in dose direction positive, i.e. dose increases
    ternaryIMBMD[abs(immy)>=activity_threshold & gx<0]=2; # YELLOW: gradient component in dose direction negative, i.e. dose decreases
    ternaryIMBMD[abs(immy)<activity_threshold]=0; # BLUE: again non-eligible region
  }
  
  pixelsGreen = sum(ternaryIMBMD==1);
  pixelsYellow = sum(ternaryIMBMD==2);
  
  #how many pixel do the yellow and green area have on the biggest dose
  n50Green = sum(ternaryIMBMD[,gridSize]==1)
  n50Yellow=sum(ternaryIMBMD[,gridSize]==2)
  
  tabG = table(which(ternaryIMBMD==1,arr.ind = T)[,2])
  tabY = table(which(ternaryIMBMD==2,arr.ind = T)[,2])
  
  minG = min(as.numeric(names(tabG)))
  minY = min(as.numeric(names(tabY)))

  if((n50Yellow+n50Green)==0){ #none of the two areaa reach the highest dose
    print("No responsive area")
    if(toPlot){
      if(tosave){
        grDevices::png(filename = paste(path, geneName, ".png", sep=""))
      }
      graphics::image(coord[,1], coord[,2], rotate(immy), xlab = "Dose",ylab = "Time", main = geneName)
      raster::contour(coord[,1], coord[,2], rotate(immy), add = TRUE,labcex = 1.3, col = "black")
    }
    return(2)
  }
  
  if(n50Yellow==0){
    score_g  = 1000
    score_y = 0
  }else{
    if(n50Green==0){
    score_g  = 0
    score_y = 1000
    }else{
      gs50 = n50Green/(n50Yellow)
      gsArea = pixelsGreen/(pixelsYellow)
      
      ys50 = n50Yellow/(n50Green+1)
      ysArea = pixelsYellow/(pixelsGreen+1)
      
      score_g  = gs50 + gsArea - minG
      score_y = ys50 + ysArea - minY
    }
  }
  
  ternaryCopy = ternaryIMBMD
  
  if (score_g > score_y){
    badDigit = 2 #yellow
    ternaryIMBMD[ternaryIMBMD !=1]=0; #everything different than 1 is set to zero. only the green part is set to 1
    verso = 1; # crescente
  }else{
    badDigit = 1 #green
    ternaryIMBMD[ternaryIMBMD == 1] = 0; #set the green part to 0.
    ternaryIMBMD[ternaryIMBMD == 2] = 1; #set the yellow part to 1.
    verso = -1; # decrescente
  }
  
  # remove non BMD rows, cioe' rimuovi le righe, dove non ho tutti 1 a partire dalla dose BMD e la dose massima
  #toSetAsZero = c()
  for(i in 1:nrow(ternaryIMBMD)){
    nOnes = sum(ternaryIMBMD[i,])
    if(nOnes>0){
      numbers = which(ternaryIMBMD[i,]==1)
      consecNum = split(numbers, cumsum(c(1, diff(numbers) != 1)))
      iidx = which(unlist(lapply(consecNum, FUN = function(elem){gridSize %in% elem}))==FALSE)
      
      if(length(iidx)>0){
        ternaryIMBMD[i,unlist(consecNum[iidx])] = 0
      }
      
    }
  }
  
  if(sum(ternaryIMBMD)==0){
    print("No responsive area")
    return(2)
  }
  
  res = bwtraceboundary(ternaryIMBMD = ternaryIMBMD)
  myContour = res$image_border
  
  goodPoints = c() #sono i punti esterni alla regione individuata, ovvero il primo punto per ogni riga
  for(i in 1:gridSize){
    it = which(myContour[,1]==i) #riga
    if(length(it)>0){
      id = min(myContour[it,2]) # colonna
      goodPoints = c(goodPoints,which(myContour[,1]==i & myContour[,2]==id))
    }
  }
  
  if(length(goodPoints)>0) myContour = myContour[goodPoints,]
  
  if(is.null(nrow(myContour))){
    print("No dose resp gene")
    ans = list()
    ans$immy = immy
    ans$activity_threshold = activity_threshold
    ans$gradient = NULL
    ans$binaryIMBMD = NULL
    ans$verso = NULL
    ans$label = NULL
    ans$tracedarea = NULL
    ans$bmd = NULL
    ans$IC50 = NULL
    ans$geneName = geneName
    
    class(ans) = 'TinderMIX'
    return(ans)
  }
  
  #### If we have a gene that had an activity before the one in the dose response, with a different sign we remove it, since we accept only one sign in the dose responsiveness
  # e.g. if I have a gene where in the responsive area I first have no activity then positive, then negative and then positive again, even if the last positive are is part of the BMD I need to remove it since to be BMD a gene has to show the same activity always 
  toSetAsZero = c() 
  goodIdx = c()
  for(i in 1:nrow(myContour)){
    if(myContour[i,2]>1){
      if(badDigit %in% ternaryCopy[myContour[i,1],1:myContour[i,2]-1]){ #there was ternaryCopy here
        toSetAsZero = c(toSetAsZero,myContour[i,1])
        goodIdx=c(goodIdx,i)
      }
    }
  }
  
  if(length(toSetAsZero)>0){
    ternaryIMBMD[toSetAsZero,]=0
    myContour = myContour[-goodIdx,]
  }
  
  
  if(toPlot){
  if(tosave){
    grDevices::png(filename = paste(path, geneName, ".png", sep=""))
  }
    #there was ternaryCopy
  graphics::image(coord[,1], coord[,2], rotate(ternaryIMBMD),breaks = c(-1,0,1,2),col = c("darkblue","darkgreen","brown"), xlab = "Dose",ylab = "Time", main = geneName)
  raster::contour(coord[,1], coord[,2], rotate(immy), add = TRUE,labcex = 1.3, col = "white")
  if(addLegend){
    graphics::legend(graphics::grconvertX(30, "device"), graphics::grconvertY(1, "device"),
           c("Non responsive Area", "responsive Increasing","Responsive Decreasing",  "Dose-Response","IC50"),
           col =c(NA, NA,NA,"gold", "red"),
           lty = c(NA,NA,NA,1,1),
           fill = c("darkblue","darkgreen","brown",NA,NA),
           border = c("darkblue","darkgreen","brown",NA,NA),
           lwd = c(NA,NA,NA,3,3),
           xpd = NA, ncol = 2,box.lwd = 0,box.col = "white",bg = "white")
  }

  
  }
  if(sum(ternaryIMBMD)==0 || class(myContour)=="integer"){
   # print("No dose response gene")
    ans = list()
    ans$immy = immy
    ans$activity_threshold = activity_threshold
    ans$gradient = gg
    ans$binaryIMBMD = binaryIMBMD
    ans$verso = NULL
    ans$label = NULL
    ans$tracedarea = res
    ans$bmd = NULL
    ans$IC50 = NULL
    ans$geneName = geneName
    
    class(ans) = 'TinderMIX'
    if(toPlot){
      if(tosave) grDevices::dev.off()
    }
    return(ans)
  }
  
  
  BMD = ternaryIMBMD
  BMD[BMD==0] = NA
  if(toPlot){
    graphics::image(coord[,1], coord[,2],rotate(BMD), col = grDevices::adjustcolor( "yellow", alpha.f = 0.2), add = T)
  }
  
  bmd_area_idx = which(BMD==1,arr.ind = T)

  gradi = c()
  modL = c()
  for(i in 1:nrow(bmd_area_idx)){
    yc = gy[bmd_area_idx[i,1],bmd_area_idx[i,2]]
    xc = gx[bmd_area_idx[i,1],bmd_area_idx[i,2]]
    dg = rad2deg(atan2(yc, xc))
    md = sqrt(yc^2 + xc^2)
    gradi = c(gradi, dg)
    modL = c(modL,md)
  }
  
  tdd = (1/sum(modL)) * sum(gradi * modL)
  if(tdd>0 & tdd<90){
    degree_lab = "Time -"
  }else{
    if(tdd>=90 & tdd<=180){
      degree_lab = "Time +"
    }else{
      if(tdd <=0 & tdd>-90){
        degree_lab = "Time -"
      }else{
        degree_lab = "Time -"
      }
    }
  }
  
  restt0 = label2DMap(map = ternaryIMBMD,BMD = BMD, coord=coord,myContour = myContour,th = BMD_response_threshold,nDoseInt = nDoseInt,nTimeInt = nTimeInt,doseLabels = doseLabels,timeLabels = timeLabels,mode=mode,toplot = toPlot)
  
  
  ddd = cbind(coord[myContour[,2],1],coord[(gridSize + 1)-myContour[,1],2])
  colnames(ddd) = c("Dose","Time")
  
  if(toPlot){
    graphics::lines(ddd, col = "gold", lwd = 4, pch=16)
    graphics::points(ddd, col = "gold", lwd = 4, pch=16)
  }
  
  IC50 = c()
  for(i in 1:nrow(immy)){
    
    riga0 = immy[i,]
    
    riga = immy[i,is.na(BMD[i,])==FALSE]
    if(length(riga)==0){
      IC50 = c(IC50,NA)
    }
    if(length(riga)==1){
      IC50 = c(IC50,which(riga0 == riga))
    }else{
      riga2 = sort(riga)
      IC50=c(IC50,which(riga0==riga2[round(length(riga)/2)]))
    }

  }
  
  IC50 = IC50[myContour[,1]]
  
  ddd2= cbind(coord[IC50,1],coord[(gridSize + 1)-myContour[,1],2])
  colnames(ddd2) = c("Dose","Time")
  
  if(toPlot){
    graphics::lines(ddd2,col = "red", lwd = 4)
    graphics::points(ddd2,col = "red", pch = 16)
    if(tosave) grDevices::dev.off()
    
  }

  ans = list()
  ans$immy = immy
  ans$activity_threshold = activity_threshold
  ans$gradient = gg
  ans$binaryIMBMD = binaryIMBMD
  ans$verso = verso
  ans$label = restt0
  ans$tracedarea = res
  ans$bmd = ddd
  ans$IC50 = ddd2
  ans$geneName = geneName
  ans$timedosedegree = tdd
  ans$gradi = gradi
  ans$modL = modL
  ans$degree_lab = degree_lab
  class(ans) = 'TinderMIX'
  return(ans)
  
}

