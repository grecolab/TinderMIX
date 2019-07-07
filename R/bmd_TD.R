#'
#' This function identify the BMD area and the IC50 value in the time and dose maps 
#'
#' @param immy z-maps of the fitted 3D model, with doses on the columns and time points on the rows
#' @param coord matrix with x and y coordinate. The first column contain the doses, while the second one the time points
#' @param geneName is a character string containing the gene name
#' @BMD_threshold threshold defining the responsive gene value
#' 
#' @return an object of class TinderMIX containing the fitted BMD object, the IC50 value. The function plot the map showing the responsive region.
#'
#' @examples
#' data("FC_WY14643")
#' exp_data = fc_data
#' pheno_data = pdata
#' PvalMat = compute_anova_dose_time(exp_data, pheno_data,dose_index = 2,time_point_index = 3)
#' ItemsList = build_items_list(PvalMat)
#' responsive_genes = unique(c(unlist(ItemsList$Dose),unlist(ItemsList$Time),unlist(ItemsList$`Dose:Time:DoseTime`),unlist(ItemsList$`Dose:Time`)))
#' contour_res = create_contour(exp_data, pheno_data, responsive_genes,dose_index = 2,time_point_index =3 ,gridSize = 50)
#' geneName = "Fam129a"
#' immy = contour_res$RPGenes[[geneName]][[3]]
#' coord = cbind(contour_res$RPGenes[[geneName]][[1]],contour_res$RPGenes[[geneName]][[2]])
#' res = compute_BMD_IC50(immy,coord, geneName,BMD_threshold = 0.58)
#' @export
#'

compute_BMD_IC50 = function(immy,coord, geneName,BMD_threshold = 0.58){
  
  dosesSTDEpsLeft = 1e-10
  dosesSTDEpsSouth = 1e-10
  flagVerticalBMD = 0
  
  immy = rotate(rotate(rotate(immy)))

  if(max(abs(immy))<BMD_threshold){
    print("No DE gene")
    return(NULL)
  }
  
  binaryIMBMD = immy
  binaryIMBMD[abs(binaryIMBMD)>=BMD_threshold]=1; # eligible region
  binaryIMBMD[abs(binaryIMBMD)<BMD_threshold]=0; # non-eligible region
  
  #image(rotate(binaryIMBMD)) # this image shows the eligible and non eligible region
  gg = gradient(immy,coord[,1], coord[,2])
  gx = gg$X
  gy = gg$Y
  X <- meshgrid(coord[,1],coord[,1])$X
  Y <- meshgrid(coord[,2],coord[,2])$Y
  
  #IDENTIFY BMD REGION
  if(sum(abs(binaryIMBMD) < BMD_threshold) == 0){  # % if non-eligible region in empty
    ternaryIMBMD = 0 * immy;
    ternaryIMBMD[abs(immy) >= BMD_threshold & gx>=0]=1; # GREEN: gradient component in dose direction positive, i.e. dose increases
    ternaryIMBMD[abs(immy) >= BMD_threshold & gx<0]=0; #% YELLOW: gradient component in dose direction negative, i.e. dose decreases
    #image(ternaryIMBMD)
  }else{
    ternaryIMBMD = immy;
    ternaryIMBMD[abs(ternaryIMBMD)>=BMD_threshold & gx>=0]=1; # GREEN: gradient component in dose direction positive, i.e. dose increases
    ternaryIMBMD[abs(ternaryIMBMD)>=BMD_threshold & gx<0]=2; # YELLOW: gradient component in dose direction negative, i.e. dose decreases
    ternaryIMBMD[abs(ternaryIMBMD)<BMD_threshold]=0; # BLUE: again non-eligible region
    # image(coord[,1], coord[,2],rotate(ternaryIMBMD))
    # contour(coord[,1], coord[,2],immy, col="black", add = TRUE)
    # quiver(X,Y, gx, gy, scale = 0.5, col="blue")
  }
  
  pixelsGreen = sum(ternaryIMBMD==1);
  pixelsYellow =sum(ternaryIMBMD==2);
  
  #ternaryIMBMD = apply(ternaryIMBMD,2,rev)
  
  if (pixelsGreen> pixelsYellow){
    ternaryIMBMD[ternaryIMBMD !=1]=0; #everything different than 1 is set to zero. only the green part is set to 1
    xx = which(ternaryIMBMD==1, arr.ind = T)[1,] # finds the first point on the contour in the image to initialize bwtraceboundary
    r = xx[1]
    c = xx[2]
    verso = 1; # crescente
  }else{
    ternaryIMBMD[ternaryIMBMD == 1] = 0; #set the green part to 0.
    ternaryIMBMD[ternaryIMBMD == 2] = 1; #set the yellow part to 1.
    xx = which(ternaryIMBMD==1, arr.ind = T)[1,] #finds the first point on the contour in the image to initialize bwtraceboundary
    r = xx[1]
    c = xx[2] 
    verso = -1; # decrescente
  }

  #res = create_tic_tac_toe(map = apply(ternaryIMBMD, 2, rev), verso)
  restt0 = create_tic_tac_toe(map = ternaryIMBMD, verso)
  
  library(imager)
  res = bwtraceboundary(ternaryIMBMD = ternaryIMBMD)
  
  myContour = res$image_border
  #image(coord[,1], coord[,2], res$diffMat[50:1,]) #this plots the border
  #image(coord[,1], coord[,2], t(res$diffMat[50:1,])) #this plots the border
  
  pixelDoses = myContour[,2]
  pixelTimes = myContour[,1]
  
  indiciCol1 = which(pixelDoses==50)
  indiciCol2 = which(pixelTimes==1)
  indici = union(indiciCol1, indiciCol2)
  
  if(length(indici)>1){
    #myContour = myContour[-indici,]
    pixelDoses = pixelDoses[-indici]
    pixelTimes = pixelTimes[-indici]
  } 
  #plot(pixelDoses, 51- pixelTimes)
  
  kLeft = 1;
  kSouth = 1;
  sideLeft = c()
  sideSouth = c()
  
  for(i in 1:length(pixelDoses)){
    if ( (ternaryIMBMD[pixelTimes[i], pixelDoses[i]+1]==1 && pixelDoses[i]==1) ||
         (ternaryIMBMD[pixelTimes[i], pixelDoses[i]+1]==1 && 
          ternaryIMBMD[pixelTimes[i], pixelDoses[i]-1]!=1)){ #check left side
      sideLeft = rbind(sideLeft,c(pixelDoses[i], pixelTimes[i]))
      #points(sideLeft[kLeft,1], sideLeft[kLeft,2], col = "red")
      kLeft= kLeft+1;
    }else if((ternaryIMBMD[pixelTimes[i]-1, pixelDoses[i]]==1 && pixelTimes[i]==50) || 
             (ternaryIMBMD[pixelTimes[i]-1, pixelDoses[i]]==1 && ternaryIMBMD[pixelTimes[i]+1, pixelDoses[i]]!=1)){
      sideSouth = rbind(sideSouth,c(pixelDoses[i], pixelTimes[i]))
      #points(sideSouth[kSouth,1], sideSouth[kSouth,2], col = "green")   
      kSouth = kSouth+1;
      
    }
  }
  
  
  dddLeft = cbind(coord[sideLeft[,1], 1], coord[51-sideLeft[,2], 2])
  dddSouth = cbind(coord[sideSouth[,1], 1], coord[51-sideSouth[,2], 2])
  
  ddd = rbind(dddLeft,dddSouth)
  timesBMD = ddd[,2]
  dosesBMD = ddd[,1]
  
  timesBMDLeft = dddLeft[,2];
  dosesBMDLeft = dddLeft[,1];
  dosesMeanBMDLeft = mean(dosesBMDLeft)
  dosesSTDBMDLeft = std(dosesBMDLeft)
  
  timesBMDSouth = dddSouth[,2]
  dosesBMDSouth = dddSouth[,1]
  dosesMeanBMDSouth = mean(timesBMDSouth)
  dosesSTDBMDSouth = std(timesBMDSouth)
  
  image(coord[,1], coord[,2], rotate(ternaryIMBMD), col = c("darkblue","darkgreen"), xlab = "Dose",ylab = "Time", main = geneName)
  
  #image(coord[,1], coord[,2], t(ternaryIMBMD[50:1,]), col = c("darkblue","darkgreen"), xlab = "Dose",ylab = "Time")
  #points(dosesBMD, timesBMD)
  fitBMDOut = NULL
  
  if (dosesSTDBMDLeft < dosesSTDEpsLeft && dosesSTDBMDSouth< dosesSTDEpsSouth){ #vertical and horizontal line
    #vertical points to be fit, LSM doesn't work, so we considier the mean and draw a line
    #horizontal points to be fit, so we considier the mean and draw a line
    flagVerticalBMD = 1;
    flagHorizontalBMD = 1;
    abline(v = dosesMeanBMDLeft,col = "yellow", lwd = 3)
    abline(h = dosesMeanBMDSouth,col = "yellow", lwd = 3)
  }else if(dosesSTDBMDLeft>= dosesSTDEpsLeft && dosesSTDBMDSouth>= dosesSTDEpsSouth){ #fit both left and south
    #non-vertical and non-horizontal points to be fit, LSM should work
    res = optimal_fitting_by_r2(dosesBMD,timesBMD)
    fitBMDOut = res$optMod
    maxAdjRSquare = res$optr2
    data = res$data
    newdat = data.frame(doses = seq(min(dosesBMD), max(dosesBMD), length.out = 100))
    newdat$pred = predict(fitBMDOut, newdata = newdat)
    lines(x = newdat$doses, y = newdat$pred, col = "yellow", lwd = 3)
    
  }else if(dosesSTDBMDLeft< dosesSTDEpsLeft && dosesSTDBMDSouth >= dosesSTDEpsSouth){ #fit only sounth and draw a vline for the left
    # vertical points to be fit, LSM doesn't work, so we considier the mean and draw a line
    # non-horizontal points to be fit, LSM should work
    flagVerticalBMD = 1;
    res = optimal_fitting_by_r2(dosesBMDSouth,timesBMDSouth)
    fitBMDOut = res$optMod
    maxAdjRSquare = res$optr2
    data = res$data
    newdat = data.frame(doses = seq(min(dosesBMDSouth), max(dosesBMDSouth), length.out = 100))
    newdat$pred = predict(fitBMDOut, newdata = newdat)
    lines(x = newdat$doses, y = newdat$pred, col = "yellow", lwd = 3)
    abline(v = dosesMeanBMDLeft,col = "yellow", lwd = 3)
    
  }else if(dosesSTDBMDLeft>= dosesSTDEpsLeft && dosesSTDBMDSouth< dosesSTDEpsSouth){
    #non-vertical points to be fit, LSM should work
    #horizontal points to be fit, so we considier the mean and draw a line
    flagHorizontalBMD = 1;
    res = optimal_fitting_by_r2(dosesBMDLeft,timesBMDLeft)
    fitBMDOut = res$optMod
    maxAdjRSquare = res$optr2
    data = res$data
    newdat = data.frame(doses = seq(min(dosesBMDLeft), max(dosesBMDLeft), length.out = 100))
    newdat$pred = predict(fitBMDOut, newdata = newdat)
    lines(x = newdat$doses, y = newdat$pred, col = "yellow", lwd = 3)
    abline(h = dosesMeanBMDSouth)
  }

  #Identify IC50
  maxDoseInMap = max(immy[ternaryIMBMD==1])
  minDoseInMap = min(immy[ternaryIMBMD==1])
  
  meanDoseInMap = (maxDoseInMap + minDoseInMap)/2;
  hh = hist(immy[ternaryIMBMD==1],plot=FALSE  )
  nBreak = which((meanDoseInMap >= hh$breaks)==FALSE)[1]-1
  perc = hh$counts[nBreak] / sum(hh$counts)
  
  #computing tolerance
  if(perc>2){
    IC50_mean_tol_perc = 0.025
  }else{
    IC50_mean_tol_perc = 0.05
  }
  
  tol = abs(meanDoseInMap*IC50_mean_tol_perc)
  
  binaryIMIC50 = immy;
  binaryIMIC50[binaryIMIC50<(meanDoseInMap - tol) | binaryIMIC50>(meanDoseInMap+tol)]=0
  binaryIMIC50[binaryIMIC50>=(meanDoseInMap -tol) & binaryIMIC50<=(meanDoseInMap+tol)]=1
  #image(coord[,1], coord[,2],rotate(binaryIMIC50)) #plot boundaries of the IC50 points
  
  binaryIMIC50 = binaryIMIC50 * ternaryIMBMD
  
  #binaryIMIC50 = rotate(rotate(rotate(binaryIMIC50)))
  res = which(binaryIMIC50==1, arr.ind = T)
  ccc = res[,1] #time
  bbb = res[,2] #dose
  
  ddd = cbind(coord[bbb,1], coord[51- ccc,2])
  timesIC50 = ddd[,2]
  dosesIC50 = ddd[,1]
  #points(dosesIC50, timesIC50,col = "red")
  
  dosesMeanIC50 = mean(dosesIC50)
  dosesSTDIC50 = sd(dosesIC50)
  resIC50 = NULL
  if(dosesSTDIC50< dosesSTDEpsLeft){
    #vertical points to be fit, LSM doesn't work, so we condier the mean
    abline(v=dosesMeanIC50, col = "red", lwd = 3)
  }else{
    #non-vertical points to be fit, LSM should work
    resIC50 = optimal_fitting_by_r2(doses = dosesIC50,times = timesIC50)
    fitBMDOut = resIC50$optMod
    maxAdjRSquare = resIC50$optr2
    data = resIC50$data
    i = resIC50$i
    
    newdat = data.frame(doses = seq(min(dosesIC50), max(dosesIC50),  length.out = 100))
    newdat$pred = predict(fitBMDOut, newdata = newdat)
    lines(x = newdat$doses, y = newdat$pred, col = "red", lwd = 3)
  }
  
  raster::contour(coord[,1], coord[,2], rotate(immy), add = TRUE,labcex = 1.3, col = "white")
  legend(grconvertX(30, "device"), grconvertY(1, "device"),
         c("Responsive Area", "Non responsive Area", "BMD","IC50"),
         col =c(NA, NA,"yellow", "red"),
         lty = c(NA,NA,1,1),
         fill = c("darkgreen","blue",NA,NA),
         border = c("darkgreen","blue",NA,NA),
         lwd = c(NA,NA,3,3),
         xpd = NA, ncol = 2,box.lwd = 0,box.col = "white",bg = "white")
  
  ans = list()
  ans$immy = immy
  ans$BMD_threshold = BMD_threshold
  ans$gradient = gg
  ans$binaryIMBMD = binaryIMBMD
  ans$verso = verso
  ans$label = restt0
  ans$tracedarea = res
  ans$fittedBMD = fitBMDOut
  ans$IC50 = meanDoseInMap
  ans$binaryIMIC50=binaryIMIC50
  ans$IC50fit = resIC50
  class(ans) = 'TinderMIX'
  return(ans)
  
}



create_tic_tac_toe = function(map, verso){
  # which(map == min(map),arr.ind = T)
  # which(map == max(map),arr.ind = T)
  
  ranges = cut(seq(5, 20, length.out=50),3)
  
  low = which(ranges==levels(ranges)[1])
  middle = which(ranges==levels(ranges)[2])
  high = which(ranges==levels(ranges)[3])
  ids = list(low,middle,high)
  
  #ttt is a 3x3 matrix, with low-mid-high on the column and high-mid-low on the rows.
  ttt = matrix(0, 3,3)
  rownames(ttt) = c("High","Mid","Low")
  colnames(ttt) = c("Low","Mid","High")
  for(i in 3:1){
    for(j in 1:3){
      ttt[i,j] = mean(mean(map[ids[[i]],ids[[j]]]))
    }
  }
  
  qq = quantile(ttt)
  print(qq)
  if(abs(qq[5]) >= abs(qq[1])){
    ttt_lab = ttt>=qq[4]
  }else{
    ttt_lab= ttt<=qq[2]
  }
  
  ttt = ttt*verso
  
  return(list(tic_tac_toe = ttt, ttt_label = ttt_lab))
  
}

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
  
  #image(diffMat)
  
  image_border = which(diffMat==1,arr.ind = 1)
  return(list(diffMat=diffMat,image_border=image_border))
}

optimal_fitting_by_r2 = function(doses, times){
  n = 1:4
  modelList = list()
  adjustedR2List = c()
  for(i in n){
    model=lm( bquote( times ~ poly(doses,.(i)) ), data=data.frame(doses=doses, times = times)) 
    
    # model <- lm(times ~ poly(doses,i))
    sm = summary(model)
    modelList[[i]] = model
    adjustedR2List = c(adjustedR2List,sm$adj.r.squared)
  }
  
  optIdx = which.max(adjustedR2List)
  fit = modelList[[optIdx]]
  data = data.frame(doses,times)
  
  return(list(optMod = fit, optr2 = adjustedR2List[optIdx], data = data,i=optIdx))
  
}

rotate <- function(x) t(apply(x, 2, rev))



