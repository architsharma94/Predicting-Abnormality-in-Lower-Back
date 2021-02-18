graphics.off() # Closing all of previous R's graphics windows.
rm(list=ls())  # CLearing R Memory
library(ggplot2)
library(ggpubr)
library(ks)
library(rjags)
library(runjags)
setwd("E:\\ABS - Applied Bayesian Stats\\z JAGS")
source("DBDA2E-utilities.R") 
#setwd("~/Documents/MATH2269_Bayesian/2020/presentations/Module 7/Application2/")


#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================

smryMCMC_HD = function(  codaSamples , compVal = NULL,  saveName=NULL) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  paramName = colnames(mcmcMat)
  for ( pName in paramName ) {
    if (pName %in% colnames(compVal)){
      if (!is.na(compVal[pName])) {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] , 
                                                          compVal = as.numeric(compVal[pName]) ))
      }
      else {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
      }
    } else {
      summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
    }
  }
  rownames(summaryInfo) = paramName
  
  # summaryInfo = rbind( summaryInfo , 
  #                      "tau" = summarizePost( mcmcMat[,"tau"] ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , data , xName="x" , yName="y" ,
                     showCurve=FALSE ,  pairsPlot=FALSE ,
                     saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]
  x = as.matrix(data[,xName])
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  guess = mcmcMat[,"guess"]
  zbeta0 = mcmcMat[,"zbeta0"]
  zbeta  = mcmcMat[,grep("^zbeta$|^zbeta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { zbeta = matrix( zbeta , ncol=1 ) }
  beta0 = mcmcMat[,"beta0"]
  beta  = mcmcMat[,grep("^beta$|^beta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { beta = matrix( beta , ncol=1 ) }
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph()
    nPtToPlot = 1000
    plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.5 ) # was cex=cex.cor*r
    }
    pairs( cbind( beta0 , beta , guess )[plotIdx,] ,
           labels=c( "beta[0]" , 
                     paste0("beta[",1:ncol(beta),"]\n",xName) , "guessing" ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Data with posterior predictive:
  # If only 1 predictor:
  if ( ncol(x)==1 ) {
    openGraph(width=7,height=6)
    par( mar=c(3.5,3.5,2,1) , mgp=c(2.0,0.7,0) )
    plot( x[,1] , y , xlab=xName[1] , ylab=yName , 
          cex=2.0 , cex.lab=1.5 , col="black" , main="Data with Post. Pred." )
    abline(h=0.5,lty="dotted")
    cVec = floor(seq(1,chainLength,length=30))
    xWid=max(x)-min(x)
    xComb = seq(min(x)-0.1*xWid,max(x)+0.1*xWid,length=201)
    for ( cIdx in cVec ) {
      lines( xComb , 
             guess[cIdx]*0.5 
             + (1.0-guess[cIdx])*1/(1+exp(-(beta0[cIdx]+beta[cIdx,1]*xComb ))) , 
             lwd=1.5 ,
             col="skyblue" )
      xInt = -beta0[cIdx]/beta[cIdx,1]
      arrows( xInt,0.5, xInt,-0.04, length=0.1 , col="skyblue" , lty="dashed" )
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"DataThresh",sep=""), type=saveType)
    }
  }
  # If only 2 predictors:
  if ( ncol(x)==2 ) {
    openGraph(width=7,height=7)
    par( mar=c(3.5,3.5,2,1) , mgp=c(2.0,0.7,0) )
    plot( x[,1] , x[,2] , pch=as.character(y) , xlab=xName[1] , ylab=xName[2] ,
          col="black" , main="Data with Post. Pred.")
    cVec = floor(seq(1,chainLength,length=30))
    for ( cIdx in cVec ) {
      abline( -beta0[cIdx]/beta[cIdx,2] , -beta[cIdx,1]/beta[cIdx,2] , col="skyblue" )
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"DataThresh",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  
  decideOpenGraph = function( panelCount , saveName , finished=FALSE , 
                              nRow=1 , nCol=3 ) {
    # If finishing a set:
    if ( finished==TRUE ) {
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,ceiling((panelCount-1)/(nRow*nCol))), 
                   type=saveType)
      }
      panelCount = 1 # re-set panelCount
      return(panelCount)
    } else {
      # If this is first panel of a graph:
      if ( ( panelCount %% (nRow*nCol) ) == 1 ) {
        # If previous graph was open, save previous one:
        if ( panelCount>1 & !is.null(saveName) ) {
          saveGraph( file=paste0(saveName,(panelCount%/%(nRow*nCol))), 
                     type=saveType)
        }
        # Open new graph
        openGraph(width=nCol*7.0/3,height=nRow*2.0)
        layout( matrix( 1:(nRow*nCol) , nrow=nRow, byrow=TRUE ) )
        par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
      }
      # Increment and return panel count:
      panelCount = panelCount+1
      return(panelCount)
    }
  }
  
  # Original scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(beta[0]) , main="Intercept" )
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
    histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx] )
  }
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( guess , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote("Mix. Coef.") , main="'guessing'" )
  panelCount = decideOpenGraph( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMarg") )
  
  # Standardized scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( zbeta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(z*beta[0]) , main="Intercept" )
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
    histInfo = plotPost( zbeta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(z*beta[.(bIdx)]) , main=xName[bIdx] )
  }
  panelCount = decideOpenGraph( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMargZ") )
  
  #-----------------------------------------------------------------------------
}

#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================

#Loading the data file
myData <- read.csv("train_spine.csv")
all_var = c("pelvic_incidence",
             "pelvic_tilt",
             "lumbar_lordosis_angle", 
             "sacral_slope",
             "pelvic_radius",
             "degree_spondylolisthesis",
             "pelvic_slope",
             "Direct_tilt",
             "thoracic_slope",
             "cervical_tilt",
             "sacrum_angle",
             "scoliosis_slope", "Class_att")
colnames(myData) <- all_var
head(myData)

#Converting the dependent variable to numeric
# Assigning 1 as Normal body structure and 0 as abnormal body structure
myData$Class_att <- as.numeric(as.factor(myData$Class_att)) - 1 # To get 0/1 instead of 1/2; Abnormal = 0; Normal = 1
head(myData)

# Summary Statistics
summary(myData)

#=== Descriptive look ===

# Distribution of Class Attrbute/ Target 
ggplot(myData, aes(x=factor(Class_att)))+
  geom_bar( width=0.7, fill="steelblue")+
  theme_minimal()+ labs(title="Frequency of Target Variable", 
                        x="Type of Body Structure (0:Abnormal, 1:Normal)", y = "No. of Observations")

# Scatter plots

myplots <- list()  # new empty list
for (i in 1:12) {
  p1 <- eval(substitute(
    ggplot(data=myData,aes(x=myData[ ,i], y = Class_att))+ 
      geom_point() +
      xlab(colnames(myData)[ i])
    ,list(i = i)))
  print(i)
  print(p1)
  myplots[[i]] <- p1  # add each plot into plot list
}


figure = ggarrange(myplots[[1]],myplots[[2]],myplots[[3]],myplots[[4]],myplots[[5]],myplots[[6]],myplots[[7]],myplots[[8]],myplots[[9]],myplots[[10]],myplots[[11]],myplots[[12]], ncol=3, nrow = 4)
figure


# THE DATA.
y = myData[,"Class_att"]
x = as.matrix(myData[,c(1:12)])

cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
show(round(cor(x),3) )
cat("\n")

x = as.matrix(myData[,c(2,4:6)])

# Loading the test data i.e. for evaluating the model
PredData = read.csv("test_spine.csv")
colnames(PredData) <- all_var
head(PredData)
PredData = PredData[,c(2,4:6,13)]
PredData$Class_att <- as.numeric(as.factor(PredData$Class_att)) - 1
summary(PredData)

xPred = as.matrix(PredData[,c(1:5)])

Nx = ncol(x)

PredData

# Specifying the data in a list, for later shipment to JAGS:
dataList <- list(
  x = x ,
  y = y ,
  xPred = xPred ,
  Ntotal = length(y),
  Nx = Nx, 
  Npred = nrow(xPred)
)


modelString = "
  # Standardize the data:
  data {
    for ( j in 1:Nx ) {
      xm[j]  <- mean(x[,j])
      xsd[j] <-   sd(x[,j])
      for ( i in 1:Ntotal ) {
        zx[i,j] <- ( x[i,j] - xm[j] ) / xsd[j]
      }
    }
  }
  # Specify the model for standardized data:
  model {
    for ( i in 1:Ntotal ) {
      # In JAGS, ilogit is logistic:
      y[i] ~ dbern( mu[i] )
      mu[i] <- ( guess*(1/2) 
                 + (1.0-guess)*ilogit(zbeta0+sum(zbeta[1:Nx]*zx[i,1:Nx])) )
    }
    # Priors vague on standardized scale:
    zbeta0 ~ dnorm( 0 , 1/2^2 )  
    for ( j in 1:Nx ) {
      zbeta[j] ~ dnorm( 0 , 1/2^2 )
    }
    guess ~ dbeta(1,9)
    # Transform to original scale:
    beta[1:Nx] <- zbeta[1:Nx] / xsd[1:Nx] 
    beta0 <- zbeta0 - sum( zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )
    
    # Compute predictions at every step of the MCMC
  for ( k in 1:Npred){
    pred[k] <- ilogit(beta0 + sum(beta[1:Nx] * xPred[k,1:Nx]))
  }
  }
  " # close quote for modelString
# Write out modelString to a text file
writeLines( modelString , con="TEMPmodel.txt" )

parameters = c( "zbeta0" , "beta0")
for ( i in 1:Nx){
  parameters = c(parameters, paste0("zbeta[",i,"]"))
}
for ( i in 1:Nx){
  parameters = c(parameters,  paste0("beta[",i,"]")) 
}
for ( i in 1:nrow(xPred)){
  parameters = c(parameters, paste0("pred[",i,"]")) 
}

# Specification of MCMC model runs parameters
parameters = c(parameters, "guess")
adaptSteps = 1000  # Number of steps to "tune" the samplers
burnInSteps = 1000
nChains = 4
thinSteps = 25
numSavedSteps = 2000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" ,
                        monitor=parameters  ,
                        data=dataList ,
                        # inits=initsList ,
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps ,
                        sample=numSavedSteps ,
                        thin=thinSteps , summarise=FALSE , plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )

#save.image(file="run-.RData")
#load(file="run-5.RData") # Load the results with 124,000 iterations

diagMCMC( codaSamples , parName="beta0" )
for ( i in 1:Nx){
  diagMCMC( codaSamples , parName=paste0("beta[",i,"]") )
}


compVal <- data.frame("beta0" = 15, "beta[1]" = 0, "beta[2]" = 0, "beta[3]" = 0, "beta[4]" =  0,  "beta[5]" =  0, 
                      "beta[6]" =  0, check.names=FALSE)

summaryInfo <- smryMCMC_HD( codaSamples = codaSamples , compVal = NULL )
print(summaryInfo)



plotMCMC(codaSamples = codaSamples , data = myData, xName=c("pelvic_tilt",
                                                            "sacral_slope",
                                                            "pelvic_radius",
                                                            "degree_spondylolisthesis"), yName="Class_att")




# Predictions for full records in training set
preds <- data.frame(PredProb = summaryInfo[12:42,3])

threshold <- 0.5# summaryInfo[427,3]
preds[which(preds$PredProb<threshold),2] <- 0
preds[which(preds$PredProb>threshold),2] <- 1

table(preds)

PredData[,5]
actual = PredData[,5]



# ============ Predictive check ============

confusionMatrix <- function(resp, pred){
  classRes <- data.frame(response = resp , predicted = pred)
  conf = xtabs(~ predicted + response, data = classRes)
  
  accuracy = sum(diag(conf))/sum(conf)
  accuracy
  precision = conf[1,1]/(conf[1,1]+conf[1,2])
  precision
  recall = conf[1,1]/(conf[1,1]+conf[2,1])
  recall
  Fscore = 2*((precision*recall)/(precision+recall))
  Fscore
  return(list(accuracy = accuracy, precision = precision, recall = recall, Fscore = Fscore, conf = conf))
}

confusionMatrix(resp = PredData[,5], pred = preds[,2])
