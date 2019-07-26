#'
#' This function assigns a label to contour map
#'
#' @param map matrix containing the z-map for a specific gene or cluster prototype
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
label2DMap = function(map, coord, myContour, th=0.95, mode = "cumulative", nDoseInt=3, nTimeInt=3, doseLabels = c("Late","Middle","Early"), timeLabels = c("Sensitive","Intermediate","Resilient"), toplot = FALSE){
  mapSize = ncol(map)
  rangesDoses = cut(seq(5, 20, length.out=mapSize),nDoseInt)
  rangesTimes = cut(seq(5, 20, length.out=mapSize),nTimeInt)
  
  idsDoses = list()
  for(i in 1:nDoseInt){
    idsDoses[[i]] = which(rangesDoses==levels(rangesDoses)[i])
    if(toplot){
      graphics::abline(v = coord[idsDoses[[i]],1][length(idsDoses[[i]])], lty = 2)
    }
  }
  
  idsTimes = list()
  for(i in 1:nTimeInt){
    idsTimes[[i]] = which(rangesTimes==levels(rangesTimes)[i])
    if(toplot){
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
  }else{
    ttt_lab = 0 * ttt
    ttt_lab[ttt>0] = 1
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
  
  return(list(tic_tac_toe = ttt, ttt_label = ttt_lab))
  
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
        goodGenes = c(goodGenes,geneName)
      }
    }, error = function(e) {
      print(NULL)
    })
    setTxtProgressBar(pb,i)
  }
  close(pb)
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
  
  return(list(Mat=Mat,MMA=MMA,GeneRes=GeneRes, labels=labels))
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
  
  #activity_threshold = log2((10 + (10*0.01))/10)
  
  
  #activity_threshold = quantile(immy, activity_threshold)
  
  if(max(abs(immy))<activity_threshold){
    print("No DE gene")
    
    if(toPlot){
      if(tosave){
        grDevices::png(filename = paste(path, geneName, ".png", sep=""))
      }
      graphics::image(coord[,1], coord[,2], rotate(immy), col = c("darkblue","darkgreen","brown"), xlab = "Dose",ylab = "Time", main = geneName)
      raster::contour(coord[,1], coord[,2], rotate(immy), add = TRUE,labcex = 1.3, col = "white")
    }
    
    return(1)
  }
  
  binaryIMBMD = immy
  elig = which((abs(binaryIMBMD)<activity_threshold)==TRUE)
  nonelig = which((abs(binaryIMBMD)<activity_threshold)==FALSE)
  
  binaryIMBMD[elig]=0; # non-eligible region
  binaryIMBMD[nonelig]=1; # eligible region
  
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
     #image(coord[,1], coord[,2],rotate(ternaryIMBMD))
    # contour(coord[,1], coord[,2],immy, col="black", add = TRUE)
    # quiver(X,Y, gx, gy, scale = 0.5, col="blue")
  }
  
  # image(coord[,1], coord[,2],rotate(ternaryIMBMD))
  # contour(coord[,1], coord[,2],immy, col="black", add = TRUE)
  # quiver(X,Y, gx, gy, scale = 0.5, col="blue")
  
  pixelsGreen = sum(ternaryIMBMD==1);
  pixelsYellow =sum(ternaryIMBMD==2);
  
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
      graphics::image(coord[,1], coord[,2], rotate(immy), col = c("darkblue","darkgreen","brown"), xlab = "Dose",ylab = "Time", main = geneName)
      raster::contour(coord[,1], coord[,2], rotate(immy), add = TRUE,labcex = 1.3, col = "white")
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
    badDigit = 2
    ternaryIMBMD[ternaryIMBMD !=1]=0; #everything different than 1 is set to zero. only the green part is set to 1
    xx = which(ternaryIMBMD==1, arr.ind = T)[1,] # finds the first point on the contour in the image to initialize bwtraceboundary
    r = xx[1]
    c = xx[2]
    verso = 1; # crescente
  }else{
    badDigit = 1

    ternaryIMBMD[ternaryIMBMD == 1] = 0; #set the green part to 0.
    ternaryIMBMD[ternaryIMBMD == 2] = 1; #set the yellow part to 1.
    xx = which(ternaryIMBMD==1, arr.ind = T)[1,] #finds the first point on the contour in the image to initialize bwtraceboundary
    r = xx[1]
    c = xx[2] 
    verso = -1; # decrescente
  }
  
  # remove non BMD rows, cioe' rimuovi le righe, dove non ho tutti 1 a partire dalla dose BMD e la dose massima
  toSetAsZero = c()
  for(i in 1:nrow(ternaryIMBMD)){
    nOnes = sum(ternaryIMBMD[i,])
    if(nOnes>0){
      numbers = which(ternaryIMBMD[i,]==1)
      n = min(numbers)
      if(abs(sum(numbers - n:gridSize))>0){
        toSetAsZero=c(toSetAsZero, i)
      }
    }
  }
  
  if(length(toSetAsZero)>0){
    ternaryIMBMD[toSetAsZero,]=0
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
    id = min(myContour[it,2]) # colonna
    goodPoints = c(goodPoints,which(myContour[,1]==i & myContour[,2]==id))
  }
  
  if(length(goodPoints)>0) myContour = myContour[goodPoints,]
  
  toSetAsZero = c() #### 
  goodIdx = c()
  for(i in 1:nrow(myContour)){
    if(myContour[i,2]>1){
      if(badDigit %in% ternaryCopy[myContour[i,1],1:myContour[i,2]-1]){
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
  graphics::image(coord[,1], coord[,2], rotate(ternaryCopy), col = c("darkblue","brown","darkgreen"), xlab = "Dose",ylab = "Time", main = geneName)
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
  if(sum(ternaryIMBMD)==0){
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
  
  restt0 = label2DMap(map = ternaryIMBMD,coord=coord,myContour = myContour,th = BMD_response_threshold,nDoseInt = nDoseInt,nTimeInt = nTimeInt,doseLabels = doseLabels,timeLabels = timeLabels,mode=mode,toplot = toPlot)
  
  BMD = ternaryIMBMD
  BMD[BMD==0] = NA
  if(toPlot){
    graphics::image(coord[,1], coord[,2],rotate(BMD), col = grDevices::adjustcolor( "yellow", alpha.f = 0.2), add = T)
  }
  
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
  class(ans) = 'TinderMIX'
  return(ans)
  
}


# compute_BMD_IC50_first_version = function(immy,coord, geneName,BMD_threshold = 0.58){
# 
#   dosesSTDEpsLeft = 1e-10
#   dosesSTDEpsSouth = 1e-10
#   flagVerticalBMD = 0
# 
#   immy = rotate(rotate(rotate(immy)))
# 
#   if(max(abs(immy))<BMD_threshold){
#     print("No DE gene")
#     return(NULL)
#   }
# 
#   binaryIMBMD = immy
#   binaryIMBMD[abs(binaryIMBMD)>=BMD_threshold]=1; # eligible region
#   binaryIMBMD[abs(binaryIMBMD)<BMD_threshold]=0; # non-eligible region
# 
#   #image(rotate(binaryIMBMD)) # this image shows the eligible and non eligible region
#   gg = gradient(immy,coord[,1], coord[,2])
#   gx = gg$X
#   gy = gg$Y
#   X <- meshgrid(coord[,1],coord[,1])$X
#   Y <- meshgrid(coord[,2],coord[,2])$Y
# 
#   #IDENTIFY BMD REGION
#   if(sum(abs(binaryIMBMD) < BMD_threshold) == 0){  # % if non-eligible region in empty
#     ternaryIMBMD = 0 * immy;
#     ternaryIMBMD[abs(immy) >= BMD_threshold & gx>=0]=1; # GREEN: gradient component in dose direction positive, i.e. dose increases
#     ternaryIMBMD[abs(immy) >= BMD_threshold & gx<0]=0; #% YELLOW: gradient component in dose direction negative, i.e. dose decreases
#     #image(ternaryIMBMD)
#   }else{
#     ternaryIMBMD = immy;
#     ternaryIMBMD[abs(ternaryIMBMD)>=BMD_threshold & gx>=0]=1; # GREEN: gradient component in dose direction positive, i.e. dose increases
#     ternaryIMBMD[abs(ternaryIMBMD)>=BMD_threshold & gx<0]=2; # YELLOW: gradient component in dose direction negative, i.e. dose decreases
#     ternaryIMBMD[abs(ternaryIMBMD)<BMD_threshold]=0; # BLUE: again non-eligible region
#     # image(coord[,1], coord[,2],rotate(ternaryIMBMD))
#     # contour(coord[,1], coord[,2],immy, col="black", add = TRUE)
#     # quiver(X,Y, gx, gy, scale = 0.5, col="blue")
#   }
# 
#   pixelsGreen = sum(ternaryIMBMD==1);
#   pixelsYellow =sum(ternaryIMBMD==2);
# 
#   #ternaryIMBMD = apply(ternaryIMBMD,2,rev)
# 
#   if (pixelsGreen> pixelsYellow){
#     ternaryIMBMD[ternaryIMBMD !=1]=0; #everything different than 1 is set to zero. only the green part is set to 1
#     xx = which(ternaryIMBMD==1, arr.ind = T)[1,] # finds the first point on the contour in the image to initialize bwtraceboundary
#     r = xx[1]
#     c = xx[2]
#     verso = 1; # crescente
#   }else{
#     ternaryIMBMD[ternaryIMBMD == 1] = 0; #set the green part to 0.
#     ternaryIMBMD[ternaryIMBMD == 2] = 1; #set the yellow part to 1.
#     xx = which(ternaryIMBMD==1, arr.ind = T)[1,] #finds the first point on the contour in the image to initialize bwtraceboundary
#     r = xx[1]
#     c = xx[2]
#     verso = -1; # decrescente
#   }
# 
#   #res = create_tic_tac_toe(map = apply(ternaryIMBMD, 2, rev), verso)
#   restt0 = label2DMap(map = ternaryIMBMD, verso = verso)
# 
#   res = bwtraceboundary(ternaryIMBMD = ternaryIMBMD)
# 
#   myContour = res$image_border
#   #image(coord[,1], coord[,2], res$diffMat[50:1,]) #this plots the border
#   #image(coord[,1], coord[,2], t(res$diffMat[50:1,])) #this plots the border
# 
#   pixelDoses = myContour[,2]
#   pixelTimes = myContour[,1]
# 
#   indiciCol1 = which(pixelDoses==50)
#   indiciCol2 = which(pixelTimes==1)
#   indici = union(indiciCol1, indiciCol2)
# 
#   if(length(indici)>1){
#     #myContour = myContour[-indici,]
#     pixelDoses = pixelDoses[-indici]
#     pixelTimes = pixelTimes[-indici]
#   }
#   #plot(pixelDoses, 51- pixelTimes)
# 
#   kLeft = 1;
#   kSouth = 1;
#   sideLeft = c()
#   sideSouth = c()
# 
#   for(i in 1:length(pixelDoses)){
#     if ( (ternaryIMBMD[pixelTimes[i], pixelDoses[i]+1]==1 && pixelDoses[i]==1) ||
#          (ternaryIMBMD[pixelTimes[i], pixelDoses[i]+1]==1 &&
#           ternaryIMBMD[pixelTimes[i], pixelDoses[i]-1]!=1)){ #check left side
#       sideLeft = rbind(sideLeft,c(pixelDoses[i], pixelTimes[i]))
#       #points(sideLeft[kLeft,1], sideLeft[kLeft,2], col = "red")
#       kLeft= kLeft+1;
#     }else if((ternaryIMBMD[pixelTimes[i]-1, pixelDoses[i]]==1 && pixelTimes[i]==50) ||
#              (ternaryIMBMD[pixelTimes[i]-1, pixelDoses[i]]==1 && ternaryIMBMD[pixelTimes[i]+1, pixelDoses[i]]!=1)){
#       sideSouth = rbind(sideSouth,c(pixelDoses[i], pixelTimes[i]))
#       #points(sideSouth[kSouth,1], sideSouth[kSouth,2], col = "green")
#       kSouth = kSouth+1;
# 
#     }
#   }
# 
# 
#   dddLeft = cbind(coord[sideLeft[,1], 1], coord[51-sideLeft[,2], 2])
#   dddSouth = cbind(coord[sideSouth[,1], 1], coord[51-sideSouth[,2], 2])
# 
#   ddd = rbind(dddLeft,dddSouth)
#   timesBMD = ddd[,2]
#   dosesBMD = ddd[,1]
# 
#   timesBMDLeft = dddLeft[,2];
#   dosesBMDLeft = dddLeft[,1];
#   dosesMeanBMDLeft = mean(dosesBMDLeft)
#   dosesSTDBMDLeft = std(dosesBMDLeft)
# 
#   timesBMDSouth = dddSouth[,2]
#   dosesBMDSouth = dddSouth[,1]
#   dosesMeanBMDSouth = mean(timesBMDSouth)
#   dosesSTDBMDSouth = std(timesBMDSouth)
# 
#   image(coord[,1], coord[,2], rotate(ternaryIMBMD), col = c("darkblue","darkgreen"), xlab = "Dose",ylab = "Time", main = geneName)
# 
# 
#   #image(coord[,1], coord[,2], t(ternaryIMBMD[50:1,]), col = c("darkblue","darkgreen"), xlab = "Dose",ylab = "Time")
#   #points(dosesBMD, timesBMD)
#   fitBMDOut = NULL
# 
#   if (dosesSTDBMDLeft < dosesSTDEpsLeft && dosesSTDBMDSouth< dosesSTDEpsSouth){ #vertical and horizontal line
#     #vertical points to be fit, LSM doesn't work, so we considier the mean and draw a line
#     #horizontal points to be fit, so we considier the mean and draw a line
#     flagVerticalBMD = 1;
#     flagHorizontalBMD = 1;
#     abline(v = dosesMeanBMDLeft,col = "yellow", lwd = 3)
#     abline(h = dosesMeanBMDSouth,col = "yellow", lwd = 3)
#   }else if(dosesSTDBMDLeft>= dosesSTDEpsLeft && dosesSTDBMDSouth>= dosesSTDEpsSouth){ #fit both left and south
#     #non-vertical and non-horizontal points to be fit, LSM should work
#     res = optimal_fitting_by_r2(dosesBMD,timesBMD)
#     fitBMDOut = res$optMod
#     maxAdjRSquare = res$optr2
#     data = res$data
#     newdat = data.frame(doses = seq(min(dosesBMD), max(dosesBMD), length.out = 100))
#     newdat$pred = predict(fitBMDOut, newdata = newdat)
#     lines(x = newdat$doses, y = newdat$pred, col = "yellow", lwd = 3)
#     BMDx =  newdat$doses
#     BMDy = newdat$pred
#   }else if(dosesSTDBMDLeft< dosesSTDEpsLeft && dosesSTDBMDSouth >= dosesSTDEpsSouth){ #fit only sounth and draw a vline for the left
#     # vertical points to be fit, LSM doesn't work, so we considier the mean and draw a line
#     # non-horizontal points to be fit, LSM should work
#     flagVerticalBMD = 1;
#     res = optimal_fitting_by_r2(dosesBMDSouth,timesBMDSouth)
#     fitBMDOut = res$optMod
#     maxAdjRSquare = res$optr2
#     data = res$data
#     newdat = data.frame(doses = seq(min(dosesBMDSouth), max(dosesBMDSouth), length.out = 100))
#     newdat$pred = predict(fitBMDOut, newdata = newdat)
#     lines(x = newdat$doses, y = newdat$pred, col = "yellow", lwd = 3)
#     abline(v = dosesMeanBMDLeft,col = "yellow", lwd = 3)
#     BMDx =  newdat$doses
#     BMDy = newdat$pred
#   }else if(dosesSTDBMDLeft>= dosesSTDEpsLeft && dosesSTDBMDSouth< dosesSTDEpsSouth){
#     #non-vertical points to be fit, LSM should work
#     #horizontal points to be fit, so we considier the mean and draw a line
#     flagHorizontalBMD = 1;
#     res = optimal_fitting_by_r2(dosesBMDLeft,timesBMDLeft)
#     fitBMDOut = res$optMod
#     maxAdjRSquare = res$optr2
#     data = res$data
#     newdat = data.frame(doses = seq(min(dosesBMDLeft), max(dosesBMDLeft), length.out = 100))
#     newdat$pred = predict(fitBMDOut, newdata = newdat)
#     lines(x = newdat$doses, y = newdat$pred, col = "yellow", lwd = 3)
#     abline(h = dosesMeanBMDSouth)
#     BMDx =  newdat$doses
#     BMDy = newdat$pred
#   }
# 
#   #Identify IC50
#   maxDoseInMap = max(immy[ternaryIMBMD==1])
#   minDoseInMap = min(immy[ternaryIMBMD==1])
# 
#   meanDoseInMap = (maxDoseInMap + minDoseInMap)/2;
#   hh = hist(immy[ternaryIMBMD==1],plot=FALSE  )
#   nBreak = which((meanDoseInMap >= hh$breaks)==FALSE)[1]-1
#   perc = hh$counts[nBreak] / sum(hh$counts)
# 
#   #computing tolerance
#   if(perc>2){
#     IC50_mean_tol_perc = 0.025
#   }else{
#     IC50_mean_tol_perc = 0.05
#   }
# 
#   tol = abs(meanDoseInMap*IC50_mean_tol_perc)
# 
#   binaryIMIC50 = immy;
#   binaryIMIC50[binaryIMIC50<(meanDoseInMap - tol) | binaryIMIC50>(meanDoseInMap+tol)]=0
#   binaryIMIC50[binaryIMIC50>=(meanDoseInMap -tol) & binaryIMIC50<=(meanDoseInMap+tol)]=1
#   #image(coord[,1], coord[,2],rotate(binaryIMIC50)) #plot boundaries of the IC50 points
# 
#   binaryIMIC50 = binaryIMIC50 * ternaryIMBMD
# 
#   #binaryIMIC50 = rotate(rotate(rotate(binaryIMIC50)))
#   res = which(binaryIMIC50==1, arr.ind = T)
#   ccc = res[,1] #time
#   bbb = res[,2] #dose
# 
#   ddd = cbind(coord[bbb,1], coord[51- ccc,2])
#   timesIC50 = ddd[,2]
#   dosesIC50 = ddd[,1]
#   #points(dosesIC50, timesIC50,col = "red")
# 
#   dosesMeanIC50 = mean(dosesIC50)
#   dosesSTDIC50 = sd(dosesIC50)
#   resIC50 = NULL
#   if(dosesSTDIC50< dosesSTDEpsLeft){
#     #vertical points to be fit, LSM doesn't work, so we condier the mean
#     abline(v=dosesMeanIC50, col = "red", lwd = 3)
#   }else{
#     #non-vertical points to be fit, LSM should work
#     resIC50 = optimal_fitting_by_r2(doses = dosesIC50,times = timesIC50)
#     fitBMDOut = resIC50$optMod
#     maxAdjRSquare = resIC50$optr2
#     data = resIC50$data
#     i = resIC50$i
# 
#     newdat = data.frame(doses = seq(min(dosesIC50), max(dosesIC50),  length.out = 100))
#     newdat$pred = predict(fitBMDOut, newdata = newdat)
#     lines(x = newdat$doses, y = newdat$pred, col = "red", lwd = 3)
#     IC50x =  newdat$doses
#     IC50y = newdat$pred
#   }
# 
#   raster::contour(coord[,1], coord[,2], rotate(immy), add = TRUE,labcex = 1.3, col = "white")
#   legend(grconvertX(30, "device"), grconvertY(1, "device"),
#          c("Responsive Area", "Non responsive Area", "BMD","IC50"),
#          col =c(NA, NA,"yellow", "red"),
#          lty = c(NA,NA,1,1),
#          fill = c("darkgreen","blue",NA,NA),
#          border = c("darkgreen","blue",NA,NA),
#          lwd = c(NA,NA,3,3),
#          xpd = NA, ncol = 2,box.lwd = 0,box.col = "white",bg = "white")
# 
#   p = plot_ly(x = coord[,1],y=coord[,2], z = t(rotate(immy)), type = "contour")%>%
#     add_trace(x = IC50x, y = IC50y, type = "scatter", mode = "line", name = "IC50", width=10) %>%
#     add_trace(x = BMDx, y = BMDy, type = "scatter", mode= "line", name = "BMD",width=10) %>%
#     add_trace(x = coord[,1], y = coord[,2], z = t(rotate(ternaryIMBMD)), opacity = 0.2, showscale = FALSE)
# 
# 
#   ans = list()
#   ans$immy = immy
#   ans$BMD_threshold = BMD_threshold
#   ans$gradient = gg
#   ans$binaryIMBMD = binaryIMBMD
#   ans$verso = verso
#   ans$label = restt0
#   ans$tracedarea = res
#   ans$fittedBMD = fitBMDOut
#   ans$IC50 = meanDoseInMap
#   ans$binaryIMIC50=binaryIMIC50
#   ans$IC50fit = resIC50
#   ans$plotlyplot = p
#   class(ans) = 'TinderMIX'
#   return(ans)
# 
# }
# 
# 


# label2DMap_first_version = function(map, verso){
#   # which(map == min(map),arr.ind = T)
#   # which(map == max(map),arr.ind = T)
#   
#   ranges = cut(seq(5, 20, length.out=50),3)
#   
#   low = which(ranges==levels(ranges)[1])
#   middle = which(ranges==levels(ranges)[2])
#   high = which(ranges==levels(ranges)[3])
#   ids = list(low,middle,high)
#   
#   #ttt is a 3x3 matrix, with low-mid-high on the column and high-mid-low on the rows.
#   ttt = matrix(0, 3,3)
#   rownames(ttt) = c("High","Mid","Low")
#   colnames(ttt) = c("Low","Mid","High")
#   for(i in 3:1){
#     for(j in 1:3){
#       ttt[i,j] = mean(mean(map[ids[[i]],ids[[j]]]))
#     }
#   }
#   
#   qq = quantile(ttt)
#   #print(qq)
#   if(abs(qq[5]) >= abs(qq[1])){
#     if(qq[4]>0){
#       ttt_lab = ttt>=qq[4]
#     }else{
#       ttt_lab = ttt>=qq[5]
#     }
#   }else{
#     if(qq[2]>0){
#       ttt_lab= ttt<=qq[2]
#     }else{
#       ttt_lab= ttt<=qq[1]
#     }
#   }
#   
#   ttt = ttt*verso
#   
#   return(list(tic_tac_toe = ttt, ttt_label = ttt_lab))
#   
# } 


# optimal_fitting_by_r2 = function(doses, times){
#   n = 1:4
#   modelList = list()
#   adjustedR2List = c()
#   for(i in n){
#     model=lm( bquote( times ~ poly(doses,.(i)) ), data=data.frame(doses=doses, times = times)) 
#     
#     sm = summary(model)
#     modelList[[i]] = model
#     adjustedR2List = c(adjustedR2List,sm$adj.r.squared)
#   }
#   
#   optIdx = which.max(adjustedR2List)
#   fit = modelList[[optIdx]]
#   data = data.frame(doses,times)
#   
#   return(list(optMod = fit, optr2 = adjustedR2List[optIdx], data = data,i=optIdx))
#   
# }
