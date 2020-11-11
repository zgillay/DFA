# Version info 2020.10.22. 
# New features
#     Modified function: SensorContributionTongue: plots the calcFDA (ETpackage) daobject too
#     New function: plotFdaScores 
#     New function: plotFdaScores 
#     New function: plotFdaScores 
#     New function: getLdaVars 
#     New function: positionLegen: find the best corner for legend
#

# Version info 2020.10.12.
# New features
#     New function: SensorContributionTongue:  plots the contributions of the sensors for e-tongue

# Version info 2020.10.06.
# New features
#     New function: plotFDAmodel2: almost same as plotFDAmodel but with ellipses default confidenciaLevel=0.95 (needs to be improved)
#     New function: SensorContribution:  plots the contributions of the sensors
#     New function: daOptim: finds the optimal number of the latent variables for DA
#     New function: clc: clean screen
#     New fucntion: calcAveConfMatrix2: same as calcAveConfMatrix but it can be called with a list
#     Modified FDA_PCAmodels: it can be called with fullData instead of dadata
#     Modified plotFDAmodel: it can plot two groups



#Sample code
        # #DA (using fullData)
        #
        # GR=C_DAgroup  # Group column ID for DA   Needs to be set!!!

        # #DA Optimization
        # a=daOptima(fullData,GR,14)
        # dn<-which(a[,2]==max(a[,2]))[1]
        # 
        # #Plot cím
        # plotTitle="Title of the plot"          Needs to be set!!!
        # sub ="Sub Title"                      Needs to be set!!!

        # #DA
        # daObjectList=list()
        # for (i in 1:3) {
        #   TranInd1 = -seq(i+1, nrow(fullData$header), 3)
        #   daObjectList[[i]] <- FDA_PCAmodels(fullData,GR, NrPCs = 1:dn,  TranInd = TranInd1)
        #   sampleColor = fullData$colRep[[GR]]
        #   #plotFDAmodel(daObjectList[[i]], projCVres = TRUE, col = sampleColor[TranInd1], col2 = sampleColor[-TranInd1], main = Title)
        #   #legend("bottomright", legend= unique(fullData$heade[[GR]]), col = unique(sampleColor), pch = 16, cex = 0.8)
        #   plotFDAmodel2(daObjectList[[i]], groups= fullData$header[[GR]][TranInd1], projCVres = TRUE, col = sampleColor[TranInd1], col2 = sampleColor[-TranInd1], main = Title, sub=sub)
        # }
        # 
        # # NO Cross validation, the best modell
        # daObject <- FDA_PCAmodels(fullData,GR, NrPCs = 1:dn)
        # sampleColor = fullData$colRep[[GR]]
        # plotFDAmodel2(daObject, groups= fullData$header[[GR]], projCVres = FALSE, col = sampleColor,  main = plotTitle, sub=sub)
        # 
        # # Contribution plot
        # contribution=SensorContribution(fullData,daObject,plotTitle,corner)


clc <<- function() cat("\f")


FDA_PCAmodels <- function(fullData, ClassVar, TranInd = NA, NrPCs = NA){

  dataset <- fullData$header
  dataset$NIR <- fullData$NIR   #[,which(colnames(fullData$NIR)%in%c("X2703.5","X2601.84"))]
  dataset$col <- fullData$colRep
  
  #PCA elökészitése
  matrix.pca <- prcomp(dataset$NIR)
  PCAscores <- predict(matrix.pca)  
  
  #Data matrix készít
  i <- which(colnames(dataset) == ClassVar)
  if (any(is.na(NrPCs))){
	  D <- cbind(dataset[,i], data.frame(PCAscores))
  } else {
	  D <- cbind(dataset[,i], data.frame(PCAscores[,NrPCs]))
  }
  colnames(D)[1]<-'group'
  
  if (any(is.na(TranInd))){ #Nincs valicáció
    FDA_PCAmodel <- fda(group ~.,data = D, na.action = "na.omit")
    p <- predict(FDA_PCAmodel, newdata = D, type='variates')
    confusion <- confusion(FDA_PCAmodel)
    out <- list(FDA_PCAmodel = FDA_PCAmodel, scoors = p, confusion = confusion, matrix.pca = matrix.pca)
  } else {  #Van valicáció
    FDA_PCAmodel <- fda(group ~.,data = D[TranInd, ], na.action = "na.omit")
	  pTran <- predict(FDA_PCAmodel, D[TranInd, ], type='variates')
    pVal  <- predict(FDA_PCAmodel, D[-TranInd, ], type='variates')
    confusionTran <- confusion(FDA_PCAmodel)
    confusionVal <- confusion(FDA_PCAmodel, D[-TranInd, ])
    out <- list(FDA_PCAmodel = FDA_PCAmodel, scoors = pTran, scoorsVal = pVal, confusion = confusionTran, confusionVal = confusionVal, matrix.pca = matrix.pca)
  }
  return(out)
}#End of function


calcAveConfMatrix2 <- function(daObjectList) {
  c=list()
  for (i in 1:3){
    c[[i]] <- daObjectList[[i]]$confusion
  }
  ave_c <-c[[1]]
  for (i in 1:ncol(c[[1]])) {
    for (j in 1:nrow(c[[1]])) {
      ave_c[j, i] <- round(mean(c(c[[1]][j, i], c[[2]][j, i], c[[3]][j, i])), 2)
    }
    ave_c[, i] <- round(ave_c[, i]/sum(ave_c[, i])*100, 2)
  }
  cv=list()
  for (i in 1:3){
    cv[[i]] <- daObjectList[[i]]$confusionVal
  }
  ave_cv <-cv[[1]] 
  for (i in 1:ncol(cv[[1]])) {
    for (j in 1:nrow(cv[[1]])) {
      ave_cv[j, i] <- round(mean(c(cv[[1]][j, i], cv[[2]][j, i], cv[[3]][j, i])), 2)
    }
    ave_cv[, i] <- round(ave_cv[, i]/sum(ave_cv[, i])*100, 2)
  }
  accuracy_c <- sum(diag(ave_c))/sum(ave_c)*100
  accuracy_cv <-sum(diag(ave_cv))/sum(ave_cv)*100
  return(list(ave_c = ave_c, ave_cv = ave_cv, recognition = accuracy_c, prediction = accuracy_cv))
} #End of function


calcAveConfMatrix <- function(FDAmodel1, FDAmodel2, FDAmodel3) {
  ave_c <- c1 <- FDAmodel1$confusion
  c2 <- FDAmodel2$confusion
  c3 <- FDAmodel3$confusion
  for (i in 1:ncol(c1)) {
    for (j in 1:nrow(c1)) {
      ave_c[j, i] <- round(mean(c(c1[j, i], c2[j, i], c3[j, i])), 2)
    }
    ave_c[, i] <- round(ave_c[, i]/sum(ave_c[, i])*100, 2)
  }
  ave_cv <- c1 <- FDAmodel1$confusionVal
  c2 <- FDAmodel2$confusionVal
  c3 <- FDAmodel3$confusionVal
  for (i in 1:ncol(c1)) {
    for (j in 1:nrow(c1)) {
      ave_cv[j, i] <- round(mean(c(c1[j, i], c2[j, i], c3[j, i])), 2)
    }
    ave_cv[, i] <- round(ave_cv[, i]/sum(ave_cv[, i])*100, 2)
    
    
  }
  accuracy_c <- sum(diag(ave_c))/sum(ave_c)*100
  accuracy_cv <-sum(diag(ave_cv))/sum(ave_cv)*100
  return(list(ave_c = ave_c, ave_cv = ave_cv, recognition = accuracy_c, prediction = accuracy_cv))
} #End of function


getFdaVars <- function(FDA_PCAmodel, round = 2){
  #Vars <- round((LDAmodel$svd)^2/sum((LDAmodel$svd)^2)*100, round)
  round(c(FDA_PCAmodel$percent.explained[1], 
          (FDA_PCAmodel$percent.explained[-1] - FDA_PCAmodel$percent.explained[-length(FDA_PCAmodel$percent.explained)])), round)
} #End of function


plotFDAmodel <- function(FDA_PCAmodel, Roots = 1:2, col = 1, pch = 16, main = "", sub = "", projCVres = FALSE, 
                         col2 = 1, pch2 = "x", ...){
  #FDA_PCAmodel=daObject
  #Roots=1:2
  #projCVres = TRUE
  
  p <- FDA_PCAmodel$scoors
  vars <- getFdaVars(FDA_PCAmodel$FDA_PCAmodel)
  if (projCVres) {
	ylims <- apply(rbind(FDA_PCAmodel$scoors, FDA_PCAmodel$scoorsVal), 2, range)
	pCv <- FDA_PCAmodel$scoorsVal
  } else {
	ylims <- apply(FDA_PCAmodel$scoors, 2, range)
  }
	if (ncol(p)>1) {
		xlab=paste0('root ', Roots[1],' - ', vars[Roots[1]], "%")
		ylab=paste0('root ', Roots[2],' - ', vars[Roots[2]], "%")
  
		plot(p[,Roots[1]],p[,Roots[2]], xlab=xlab, ylab=ylab, col = col, pch = pch, main = main, sub = sub, xlim = ylims[,Roots[1]], ylim = ylims[,Roots[2]], ...)
		if (projCVres) {
			points(pCv[,Roots[1]], pCv[,Roots[2]], col = col2, pch = pch2, ...)
		}
	} else {
	  xlab="Samples"
	  ylab=paste0('root 1 ',  "100%")
		plot(p, col = col, pch = pch, main = main, sub = sub, ...)
	  if (projCVres) {
	    points(pCv, col = col2, pch = pch2, ...)
	  }
	}
  }#End of function

  
aveConfMatrix <- function(colConf){
	TrInd <- which(names(colConf) == "Training")
	VaInd <- which(names(colConf) == "Valid")
	
	AveTr <- 0
	for (i in TrInd){
		AveTr <- AveTr + (colConf[i]$Training)/colSums(colConf[i]$Training)*100
		#AveTr <- AveTr + t(t(colConf[i]$Training)/colSums(t(colConf[i]$Training))*100)
		
	}
	AveTr <- round(AveTr/3, 2)
	
	AveVa <- 0
	for (i in VaInd){
		 AveVa <- AveVa + (colConf[i]$Valid)/colSums(colConf[i]$Valid)*100
		#AveVa <- AveVa + t(t(colConf[i]$Valid)/colSums(t(colConf[i]$Valid))*100)
	}
	AveVa <- round(AveVa/3, 2)
	out <- list(Training = AveTr, Valid = AveVa)
}#End of function


plotFDAmodel2 <-function(FDA_PCAmodel, Roots = 1:2,col = 1, pch = 16,confidenciaLevel=0.95,
                         main = "", sub = "", projCVres = TRUE,col2 = 1, pch2 = 1, 
                         groups,legendFontSize = 0.8, corner=c(0,1), ...)  {
  # sample data
  
  #FDA_PCAmodel=dare
  #Roots=1:2
  #col = sampleColor[TranInd1]
  #groups=fullData$header[[GR]][TranInd1]
  
  
  x = FDA_PCAmodel$scoors[,Roots[1]]
  y = FDA_PCAmodel$scoors[,Roots[2]]
  
  if (projCVres) {
  xv = FDA_PCAmodel$scoorsVal[,Roots[1]]
  yv = FDA_PCAmodel$scoorsVal[,Roots[2]]
  }
  colorData = col

  ### prep legend ##
  legend0<-unique(groups)
  colorLegend0 <- c(unique(col)) # colors of the groups legend #
  

   df = data.frame(legend0, colorLegend0 )
   df <- df[order(df$legend0),]
   legend0<-as.character(df$legend0)
   colorLegend0<-as.character(df$colorLegend0)
  

  if (projCVres) {
    legendTextExt <- c(as.character(legend0),"Validation") # name of groups
    ylims <- apply(rbind(FDA_PCAmodel$scoors, FDA_PCAmodel$scoorsVal), 2, range)*1.4
    colorLegend<- c(colorLegend0,1)
    legendSimbols<-c(rep(pch,length(legend0)),pch2)
  }
  else {
    legendTextExt <- c(as.character(legend0))
    ylims <- apply(FDA_PCAmodel$scoors, 2, range)*1.4
    colorLegend<- colorLegend0
    legendSimbols<-c(rep(pch,length(legend0)))
  } 
  # size of legend
  
# Position of the legend  
   corner <- positionLegend(x,y,ylims,strOut = FALSE)


  trLegend1 <- list(border=TRUE, points=list(pch=legendSimbols, col=colorLegend), 
                    #lines=list(lty=1, col = unique(colorData)),
                    text=list(legendTextExt), background=0, alpha.background=1, 
                    cex=legendFontSize, columns=1)
  trAllLegends <- list(inside=list(fun=lattice::draw.key(trLegend1), corner=corner))
  
  ### prep legend END 
  
  mainText <- main
  
  if (!is.na(charmatch("FDA_PCAmodel", names(FDA_PCAmodel)))){
    v1<-FDA_PCAmodel$FDA_PCAmodel$percent.explained[["v1"]]
    v2<-FDA_PCAmodel$FDA_PCAmodel$percent.explained[["v2"]]
    xlab <-  paste0("root ",Roots[1]," - ", round(v1,2),"%")
    ylab <-  paste0("root ",Roots[2]," - ", round(v2-v1,2),"%")    
  } else if  (!is.na(charmatch("FDAmodel", names(FDA_PCAmodel)))){
    v1<-FDA_PCAmodel$FDAmodel$percent.explained[["v1"]]                  
    v2<-FDA_PCAmodel$FDAmodel$percent.explained[["v2"]]
    xlab <-  paste0("root ",Roots[1]," - ", round(v1,2),"%")
    ylab <-  paste0("root ",Roots[2]," - ", round(v2-v1,2),"%")
  } else {
    xlab <- paste0("root ",Roots[1])
    ylab <- paste0("root ",Roots[2])
  }
                  


  pch_data_el <- pch # pch of the points
  pch_center_el <- "x" # pch of group centers (set NULL to switch off)

  trRobust <- TRUE # if true calculates the conf level in robust way (less sensitive to outliers)
  
  trelPlot2 <- lattice::xyplot(y ~ x, main=mainText,sub=sub, xlab=xlab, ylab=ylab, col=colorData,
                               pch=pch_data_el, cex=0.85, legend=trAllLegends, xlim = ylims[,Roots[1]], ylim = ylims[,Roots[2]] , 
                               panel=function(x, y, ...) {
                                 
                                 lattice::panel.xyplot(x, y,   ...) # plot the data
                                 
                                 if (projCVres) {
                                   lattice::panel.xyplot(xv, yv, pch=pch2, col = col2 ,cex = 1  )
                                 }

                                 latticeExtra::panel.ellipse(x, y, col=colorLegend0 , center.pch=pch_center_el,   #col=unique(col)
                                                             robust=trRobust, groups=groups, scales="free", 
                                                             level=confidenciaLevel, lwd=1,lty=1, subscripts=TRUE)
                                 
                                  for (i in 1:length(legend0)){
                                      lattice::panel.text(mean(x[which(groups == legend0[i])]),mean(y[which(groups == legend0[i])]),
                                                          legend0[i],col=colorLegend0[i],cex=0.7)
                                  }
                               }
                               
  ) # end of call to trelPlot2
  print(trelPlot2)
}# end of function



daOptima<-function(fullData,GR,N=14){
  a=matrix(nrow = N,ncol = 3)
  a[1:3,1:3]=0 
  for (dn in 1:N){
    TranInd1 = -seq(1, nrow(fullData$header), 3)
    da1 <- FDA_PCAmodels(fullData,GR, NrPCs = 1:dn,  TranInd = TranInd1)
    
    TranInd2 = -seq(2, nrow(fullData$header), 3)
    da2 <- FDA_PCAmodels(fullData,GR, NrPCs = 1:dn,  TranInd = TranInd2)
    
    TranInd3 = -seq(3, nrow(fullData$header), 3)
    da3 <- FDA_PCAmodels(fullData,GR, NrPCs = 1:dn,  TranInd = TranInd3)
    
    conftable <- calcAveConfMatrix(da1,da2,da3)
    
    a[dn,1]=round(conftable[[3]])
    a[dn,2]=round(conftable[[4]])
    a[dn,3]=a[dn,1]-a[dn,2]
  }
  return(a)
}# End of function


SensorContributionTongue<-function(fullData,daObject,plotTitle="",corner="",isTable=FALSE){
 if (!is.na(charmatch("NIR", names(fullData)))){
        nCol=ncol(fullData$NIR)
        NIR<-matrix(nrow=(nCol+1),ncol=nCol,0)
        header<-matrix(nrow=(nCol+1),NA)
        for (i in 2:(nCol+1)) {
          NIR[i,i-1]<-1
        }
        colnames(NIR)<-colnames(fullData$NIR)
        PCAscoresHack <- predict(daObject$matrix.pca, newdata = NIR)
        D <- cbind(header, data.frame(PCAscoresHack))
        p <- predict(daObject$FDA_PCAmodel, newdata = D, type='variates')

    } else if (!is.na(charmatch("sensors", names(fullData)))){
          fullData$NIR <- fullData$sensors
          nCol=ncol(fullData$NIR)
          NIR<-matrix(nrow=(nCol+1),ncol=nCol,0)
          header<-matrix(nrow=(nCol+1),NA)
          
          for (i in 2:(nCol+1)) {
            NIR[i,i-1]<-sign(mean(fullData$NIR[,i-1]))
          }
          colnames(NIR)<-colnames(fullData$NIR)
          p <- predict(daObject$FDAmodel, newdata= NIR, type="variates")
          
        } else {print(paste0("fullData does not containt sensors or NIR"))}
  
stdev<-matrix(nrow=nCol,ncol=1)
for (i in 1:length(colnames(fullData$NIR))){
  stdev[i]<-sd(fullData$NIR[,i])
}

contribution<-matrix( nrow = nCol, ncol=(ncol(p)+1))
for (i in 1:ncol(p)){
  contribution[,i]<-(as.numeric(p[2:(nCol+1),i]-p[1,i]))*stdev
}

SUM=contribution[,1]^2
for (i in 2:ncol(p)){
  SUM=SUM+contribution[,i]^2
}
contribution[,(ncol(p)+1)]=sqrt(SUM)


Sensors<-colnames(fullData$NIR)
d<-data.frame(Sensors,stdev,contribution)
rownames(d)<-colnames(fullData$NIR)
 


colnames(d)=c(colnames(d)[1:(ncol(p)+2)],"Contribution")

d <- d[order(d[["Contribution"]],decreasing = TRUE),]

#ylims <- apply(cbind(d[["X1"]],d[["X2"]]), 2, range) #xlim = ylims[,1], ylim = ylims[,2],

# points(d[["X1"]],d[["X2"]],col=2,pch=16, xlab = "Root-1",ylab = "Root-2",cex.main=0.9,
#      main = paste("Sensor contribution in PCA-DA by Sensors \n",plotTitle),panel.first=grid())

if (!(corner=="")){
  legend(legend = c("All contributors","Major contriburiors"),col=c(2,1), pch = c(16,6), cex = 0.8, corner)
}



shift=(max(d[["X2"]])-min(d[["X2"]]))/40

for (i in 1:7){
  text(d$X1[i],d$X2[i]+shift,d$Sensors[i],cex=0.8,col=1)
  arrows(x=0, y=0, x1 = d$X1[i], y1 = d$X2[i], length = 0.15, angle = 15,
         code = 2, col = par("fg"), lty = par("lty"),  lwd = par("lwd"))
  
}
if (isTable){
  plot.new()
  testdf<-data.frame(Sensors=Sensors,Stdev=round(d$stdev,2),Contribution=round(d$Contribution,3))

  # show most of the options
  plotrix::addtable2plot(0. ,0.4,testdf,bty="o",display.rownames=TRUE,hlines=TRUE,
                vlines=TRUE,title=paste("Contribution Table, Roots: ", ncol(contribution)-3),xpad = 0.5,ypad = 1.5,)
}



return(d)
}# End of function


SensorContribution      <-function(fullData,daObject,plotTitle="",corner="",isTable=TRUE){
    # For NIR
    nCol=ncol(fullData$NIR)
    NIR<-matrix(nrow=(nCol+1),ncol=nCol)
    header<-matrix(nrow=(nCol+1),NA)
    
    NIR[1,]<-fullData$NIR[1,]
    for (i in 2:(nCol+1)) {
      NIR[i,]<-fullData$NIR[1,]
      NIR[i,i-1]<-NIR[i,i-1]+1
    }
    NIR=data.frame(NIR)
    colnames(NIR)<-colnames(fullData$NIR)
    
    PCAscoresHack <- predict(daObject$matrix.pca, newdata = NIR)
    D <- cbind(header, data.frame(PCAscoresHack))
    p <- predict(daObject$FDA_PCAmodel, newdata = D, type='variates')
 
  stdev<-matrix(nrow=nCol,ncol=1)
  for (i in 1:length(colnames(fullData$NIR))){
    stdev[i]<-sd(fullData$NIR[,i])
  }
  
  contribution<-matrix( nrow = nCol, ncol=(ncol(p)+1))
  for (i in 1:ncol(p)){
    contribution[,i]<-(as.numeric(p[2:(nCol+1),i]-p[1,i]))*stdev
  }
  
  SUM=contribution[,1]^2
  for (i in 2:ncol(p)){
    SUM=SUM+contribution[,i]^2
  }
  contribution[,(ncol(p)+1)]=sqrt(SUM)
  
  Sensors<-as.numeric(substr(colnames(NIR),2,100))
  
  d<-data.frame(Sensors,stdev,contribution)
  rownames(d)<-paste0("S",round(Sensors))
  
  colnames(d)=c(colnames(d)[1:(ncol(p)+2)],"Contribution")
  
  d <- d[order(d[["Contribution"]],decreasing = TRUE),]
  
  #ylims <- apply(cbind(d[["X1"]],d[["X2"]]), 2, range) #xlim = ylims[,1], ylim = ylims[,2],
  
  points(d[["X1"]],d[["X2"]],col=2,pch=16, xlab = "Root-1",ylab = "Root-2",cex.main=0.9,
       main = paste("Sensor contribution in PCA-DA by retention indeces \n",plotTitle),panel.first=grid())
  
    if (!(corner=="")){
      legend(legend = c("All contributors","Major contriburiors"),col=c(2,1), pch = c(16,6), cex = 0.8, corner)
    }

  
  
  numberMajorPoint=length(which(contribution[,(ncol(p)+1)]>max(contribution[,(ncol(p)+1)])/6))
  
  
  shift=(max(d[["X2"]])-min(d[["X2"]]))/40
  
  for (i in 1:numberMajorPoint){
    text(d$X1[i],d$X2[i]+shift,d$Sensors[i],cex=0.8,col=1)
    arrows(x=0, y=0, x1 = d$X1[i], y1 = d$X2[i], length = 0.15, angle = 15,
           code = 2, col = par("fg"), lty = par("lty"),  lwd = par("lwd"))
    
  }
  if (isTable){
  plot.new()
  
  
  testdf<-data.frame(Sensors=round(d$Sensors[1:15],2),Stdev=round(d$stdev[1:15],2),Contribution=round(d$Contribution[1:15],3))
  
  # show most of the options
  plotrix::addtable2plot(0. ,0.1,testdf,bty="o",display.rownames=TRUE,hlines=TRUE,
                         vlines=TRUE,title=paste("Contribution Table, Roots: ", ncol(contribution)-3),xpad = 0.5,ypad = 1.5,)
  }
  return(d)
}# End of function SensorContribution


calcFDA <- function(dataset, ClassVar, TranInd = NA){
  i <- which(colnames(dataset) == ClassVar)
  D <- cbind(dataset[,i], data.frame(dataset$sensors))
  
  colnames(D)[1]<-'group'
  if (any(is.na(TranInd))){
    FDAmodel <- mda::fda(group ~.,data = D, na.action = "na.omit")
    p <- predict(FDAmodel, newdata = D, type='variates')
    confusion <- mda::confusion(FDAmodel)
    out <- list(FDAmodel = FDAmodel, scoors = p, confusion = confusion, scores = p)
  } else {
    FDAmodel <- mda::fda(group ~.,data = D[TranInd, ], na.action = "na.omit")
    pTran <- predict(FDAmodel, D[TranInd, ], type='variates')
    pVal <- predict(FDAmodel, D[-TranInd, ], type='variates')
    confusionTran <- mda::confusion(FDAmodel)
    confusionVal <- mda::confusion(FDAmodel, D[-TranInd, ])
    out <- list(FDAmodel = FDAmodel, scores = pTran, scoresVal = pVal, confusion = confusionTran,
                confusionVal = confusionVal, scoors = pTran, scoorsVal = pVal)
  }
} # end of function calcFDA 


plotFdaScores <- function(FDAmodel, Roots = 1:2, col = 1, pch = 16, main = "", sub = "", projCVres = TRUE, 
                          col2 = 1, pch2 = "x", ...){

  if(projCVres) {
    projCVres <- any(FDAmodel$scoresVal>0)
  }
  p <- FDAmodel$scores
  vars <- getLdaVars(FDAmodel)
  if (projCVres) {
    ylims <- apply(rbind(FDAmodel$scores, FDAmodel$scoresVal), 2, range)
    pCv <- FDAmodel$scoresVal
  } else {
    ylims <- apply(FDAmodel$scores, 2, range)
  }
  if (ncol(p)>1) {
    xlab = vars[Roots[1]]
    ylab = vars[Roots[2]]
     
    
    plot(p[,Roots[1]],p[,Roots[2]], xlab=xlab, ylab=ylab, col = col, 
         pch = pch, main = main, sub = sub, xlim = ylims[,Roots[1]]*1.4, ylim = ylims[,Roots[2]]*1.4, ...)
    
    corner<- positionLegend(p[,Roots[1]],p[,Roots[2]],ylims)
    
    
    legend( corner,  legend = unique(fullDataET$group), col = unique(fullDataET$numRep$group), pch = 16)
    
    if (projCVres) {
      points(pCv[,Roots[1]], pCv[,Roots[2]], col = col2, pch = pch2, ...)
    }
  } else {
    plot(p, col = col, pch = pch, main = main, sub = sub, ...)
  }
} # end of function plotFdaScores


getLdaVars <- function(LDAmodel, round = 2){
  cumPerc <- LDAmodel$FDAmodel$percent.explained
  Vars <- round(c(cumPerc[1], diff(cumPerc)), round)
  paste0("root", 1:length(Vars), " - ", Vars, "%")
} # end of function getLdaVars 


positionLegend <- function(x,y,ylims,strOut=TRUE){
# Position of the legend  
epsx=(ylims[2,1]-ylims[1,1])/3
epsy=(ylims[2,2]-ylims[1,2])/4

a1=which( ylims[1,1]      <x & (x<ylims[1,1]+epsx)) #left
a2=which((ylims[2,1]-epsx)<x &  x<ylims[2,1]      ) # right

b1=which( ylims[1,2]      <y & (y<ylims[1,2]+epsy)) #down
b2=which((ylims[2,2]-epsy)<y &  y<ylims[2,2]      ) #up


leftdown  = length(intersect(a1,b1))
rightdown = length(intersect(a2,b1))
leftup    =length(intersect(a1,b2))
rightup   =length(intersect(a2,b2))

leng=data.frame(c(rightup,leftup,leftdown,rightdown),c(4,3,2,1))
colnames(leng)=c("points","Rank")
subleng=leng[which(leng$points==min(leng$points)),]
subleng=subleng[order(subleng$Rank),]

if (!strOut){
if (subleng$Rank[1]==4){corner= c(1,1)}
if (subleng$Rank[1]==3){corner=c(0,1)}
if (subleng$Rank[1]==2){corner=c(0,0)}
if (subleng$Rank[1]==1){corner=c(1,0)}
} else{
if (subleng$Rank[1]==4){corner= "topright"}
if (subleng$Rank[1]==3){corner="topleft"}
if (subleng$Rank[1]==2){corner="bottomleft"}
if (subleng$Rank[1]==1){corner="bottomright"}
}
return(corner)
}  # end of function positionLegend

2020.11.11.



