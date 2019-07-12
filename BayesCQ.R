library(BayesianTools)

piecewiseBayes <- function(data, fname = NULL){
  
  refPars <- matrix(NA,nrow=5,ncol=4)
  rnames <- c("a1","b1","delta","x0", "sigma")
  cnames <- c("best","sd", "lower", "upper")
  row.names(refPars) <- rnames
  colnames(refPars) <- cnames
  mnx <- mean(data$x,na.rm=TRUE)
  lm1 <-  lm(data$y~data$x, subset = data$x <= mnx)
  lm2 <-  lm(data$y~data$x, subset = data$x > mnx)
  
  #B0
  refPars[1,] <- c(coefficients(lm1)[1],
                   sd(data$y,na.rm=TRUE), 
                   -100,
                   100)
  #B1
  refPars[2,] <- c(coefficients(lm1)[2],10 , -100,100)
  #delta
  refPars[3,] <- c(coefficients(lm2)[2]-coefficients(lm1)[2],10 , -100,100)
  #x0
  refPars[4,] <- c(mean(data$x,na.rm=TRUE), 
                   sd(data$x,na.rm=TRUE), 
                   quantile(data$x,0.125),
                   quantile(data$x,0.875))
  refPars <- as.data.frame(refPars)
  
  refPars[5,] <- c(0.1, 
                   max(1,abs(sd(data$y,na.rm=TRUE))), 
                   0.001,
                   1000)
  refPars <- as.data.frame(refPars)
  parSel = c(1:5) #Parameters to select for calibration
  
  burnin <- 5000
  nsamples <- burnin*3
  
  #Function to Bayes fit a piecewise-linear (power-law) C-Q
  piecewise <- function(params,x = data$x){
    
    a1    <- params[1] #1st intercept
    b1    <- params[2] #1st slope
    delta <- params[3] #difference in slope
    x0    <- params[4] #breakpoint
    
    #y1 = a1 + b1 * x for x < X0
    #y2 = a2 + b2 * x for x>= x0
    #force lines to meet at x = X0
    #a1 + b1*X0 = a2 + b2*X0
    #a2 = a1+(b1-b2)*X0
    
    b2 <- b1 - delta #second slope by definition
    a2 <- a1 + delta*x0 #second intercept
    
    #y = a1 + b1 * x for x < X0
    #y = a1 + delta*x0 + (b1 - delta) * x for x>= x0
    
    #sigma <- params[5] #model error, not used
    
    pos <- x < x0
    
    ysim <- ifelse(pos,a1 + b1*x,a2 + b2*x)
    ysim
  }
  
  
  #likelihood
  ll <- function(par, sum = TRUE){
    params = refPars$best
    params[parSel] = par
    predicted <-   piecewise(params=params,x=data$x)
    diff <- c(predicted - data$y) 
    llValues <- dnorm(diff, sd = params[5], log = TRUE)  
    if (sum == FALSE) return(llValues)
    else return(sum(llValues))
  }
  #prior
  prior <- createTruncatedNormalPrior(
    mean = refPars$best[parSel],
    sd = refPars$sd[parSel],
    lower = refPars$lower[parSel], 
    upper = refPars$upper[parSel])
  
  setup <- createBayesianSetup(ll, prior,
                               names = rownames(refPars)[parSel])
  #Markov Chain Monte Carlo
  out <- runMCMC(bayesianSetup = setup, 
                 sampler = "DEzs", 
                 settings = list(iterations = nsamples,  
                                 message = FALSE, 
                                 nrChains = 3))
  
  #Sample from posterior
  post.parMatrix = getSample(out, start = 3000,end = burnin, thin = 5)
  
  #calc b2 so we can generate stats for it
  b2 <- post.parMatrix[,2] - post.parMatrix[,3]
  post.parMatrix <- cbind(post.parMatrix,b2)
  
  #Extract credible intervals
  CredInts <- getCredibleIntervals(post.parMatrix, quantiles = c(0.025, 0.5, 0.975))
  colnames(CredInts) <- c("a1","b1","delta","x0","sigma","b2")
  
  #get the predicted intervals 95%
  getPredInt <- function(parMatrix,model, data, numSamples = 1000,
                         quantiles = c(0.025, 0.975)){
    x2 <- sort(data$x,index.return = TRUE)
    y2 <- data$y[x2$ix]
    newdata = data.frame(x = x2$x, y = y2)
    sims <- apply(parMatrix, MARGIN = 1, FUN = function(pars) {
      model(c(pars),x = newdata$x)
    })
    ststs <- apply(sims, MARGIN = 1, function(x) quantile(x, c(0.025,0.5,0.975)))
    data.frame(x=x2$x,y=y2,p02.5 = ststs[1,],p50 = ststs[2,], p97.5 = ststs[3,])
  }
  
  gPI <- getPredInt(post.parMatrix[,1:5],piecewise, data)
  
  #do some plotting
  if(!is.null(fname)) pdf(file = fname)
  plot(data, type = "n", xlab = "", ylab = "")
  if (!is.null(data$id)) mtext(line = 1, side =3, data$id)
  if (is.null(data$nms)) {
    xnm <- "x"; 
    ynm <- "y";
  } else {
    xnm <- data$nms[1]
    ynm <- data$nms[2]
  }
  mtext(side = 1, line = 3, xnm)
  mtext(side = 2, line = 2.5, ynm)
  polygon(c(gPI$x,rev(gPI$x),gPI$x[1]),c(gPI[,"p02.5"],rev(gPI[,"p97.5"]),gPI[1,"p02.5"]), col = rgb(0,0,1,0.3), border =  rgb(0,0,1,0.3) )
  lines(gPI$x,gPI[,"p50"], col = rgb(0,0,1,0.7), lwd=2)
  points(data,pch=16)
  legend("bottom",
         col = c(rgb(0,0,1,0.3),rgb(0,0,1,0.7),"black"), 
         lwd = c(NA,2,1), 
         pch = c(NA,NA,16), 
         lty = c(NA,1,NA),
         fill = c(rgb(0,0,1,0.3), NA, NA),
         border = c(rgb(0,0,1,0.3), NA, NA),
         legend = c("95%","50%","Observed"),
         bty = "n")
  if(!is.null(fname)) dev.off()
  
  #output
  return(list(cedIntervals = CredInts,
              DIC = DIC(out),
              WAIC = WAIC(out)))
}

#Function to Bayes fit a linear (power-law) C-Q
lmBayes <- function(data, fname = NULL){
  
  refPars <- matrix(NA,nrow=3,ncol=4)
  rnames <- c("a1","b1","sigma")
  cnames <- c("best","sd", "lower", "upper")
  row.names(refPars) <- rnames
  colnames(refPars) <- cnames
  lm1 <-  lm(data$y~data$x)
  
  #B0
  refPars[1,] <- c(coefficients(lm1)[1],
                   sd(data$y,na.rm=TRUE), 
                   coefficients(lm1)[1] - 5*sd(data$y,na.rm=TRUE),
                   coefficients(lm1)[1] + 5*sd(data$y,na.rm=TRUE))
  #B1
  refPars[2,] <- c(coefficients(lm1)[2],10 , -100,100)
  #delta
  refPars[3,] <- c(0.1, 
                   max(1,abs(sd(data$y,na.rm=TRUE))), 
                   0.001,
                   100)
  
  refPars <- as.data.frame(refPars)
  
  refPars <- as.data.frame(refPars)
  parSel = c(1:3) #Parameters to select for calibration
  
  burnin <- 5000
  nsamples <- burnin*3
  #Model
  piecewise <- function(params,data = data){
    
    a1    <- params[1] #1st intercept
    b1    <- params[2] #1st slope
    #sigma <- params[3] #model error, not used
    
    ysim <- a1 + b1*data$x
    ysim
  }
  
  
  #likelihood
  ll <- function(par, sum = TRUE){
    x = refPars$best
    x[parSel] = par
    predicted <-   piecewise(x,data)
    diff <- c(predicted - data$y) 
    llValues <- dnorm(diff, sd = x[3], log = TRUE)  
    if (sum == FALSE) return(llValues)
    else return(sum(llValues))
  }
  
  #Prior distribution
  prior <- createTruncatedNormalPrior(
    mean = refPars$best[parSel],
    sd = refPars$sd[parSel],
    lower = refPars$lower[parSel], 
    upper = refPars$upper[parSel])
  
  setup <- createBayesianSetup(ll, prior,
                               names = rownames(refPars)[parSel])
  #Markov Chain Monte Carlo
  out <- runMCMC(bayesianSetup = setup, 
                 sampler = "DEzs", 
                 settings = list(iterations = nsamples,  
                                 message = FALSE, 
                                 nrChains = 3))
  
  #Sample from posterior
  post.parMatrix = getSample(out, start = 3000,end = burnin)
  
  #Credible interval
  CredInts <- getCredibleIntervals(post.parMatrix, quantiles = c(0.025, 0.5, 0.975))
  colnames(CredInts) <- c("a1","b1","sigma")
  
  #Predictive interval 95%
  getPredInt <- function(parMatrix,model, data, numSamples = 1000,
                         quantiles = c(0.025, 0.975)){
    x2 <- sort(data$x,index.return = TRUE)
    y2 <- data$y[x2$ix]
    newdata = data.frame(x = x2$x, y = y2)
    sims <- apply(parMatrix, MARGIN = 1, FUN = function(x) {
      model(c(x),data = newdata)
    })
    ststs <- apply(sims, MARGIN = 1, function(x) quantile(x, c(0.025,0.5,0.975)))
    data.frame(x=x2$x,y=y2,p02.5 = ststs[1,],p50 = ststs[2,], p97.5 = ststs[3,])
  }
  
  gPI <- getPredInt(post.parMatrix,piecewise, data)
  
              
  #Some plotting     
  if(!is.null(fname)) pdf(file = fname)
  plot(data, type = "n", xlab = "", ylab = "")
  if (!is.null(data$id)) mtext(line = 1, side =3, data$id)
  if (is.null(data$nms)) {
    xnm <- "x"; 
    ynm <- "y";
  } else {
    xnm <- data$nms[1]
    ynm <- data$nms[2]
  }
  mtext(side = 1, line = 3, xnm)
  mtext(side = 2, line = 2.5, ynm)
  
  polygon(c(gPI$x,rev(gPI$x),gPI$x[1]),c(gPI[,"p02.5"],rev(gPI[,"p97.5"]),gPI[1,"p02.5"]), col = rgb(0,0,1,0.3), border =  rgb(0,0,1,0.3) )
  lines(gPI$x,gPI[,"p50"], col = rgb(0,0,1,0.7), lwd=2)
  points(data,pch=16)
  legend("bottom",
         col = c(rgb(0,0,1,0.3),rgb(0,0,1,0.7),"black"), 
         lwd = c(NA,2,1), 
         pch = c(NA,NA,16), 
         lty = c(NA,1,NA),
         fill = c(rgb(0,0,1,0.3), NA, NA),
         border = c(rgb(0,0,1,0.3), NA, NA),
         legend = c("95%","50%","Observed"),
         bty = "n")
  if(!is.null(fname)) dev.off()
  return(list(credIntervals = CredInts,
              DIC = DIC(out),
              WAIC = WAIC(out))
  )
}

#Read the shapefile names from the databse to get a list of sites                   
sites <- list.files(path = "D:\\WIRT\\GIS",
                    pattern = "*.shp")
sites <- sites[grep("Basin",sites)]
SITE_REFs <- unlist(strsplit(sites,"Basin"))
SITE_REFs <- SITE_REFs[SITE_REFs != ""]
SITE_REFs <- unlist(strsplit(SITE_REFs,".shp"))

#Create a vector of rds file names saved from file xxxx                   
files <- paste0(paste0("D:\\WIRT\\Processed_CQDiscrete\\CQ_",SITE_REFs),".RDS")
CQExists <- file.exists(files)


setwd("D:\\WIRT\\Processed_CQDiscrete")
all.rds <- list.files(pattern="*.RDS")
all.rds <- all.rds[grep("CQ",all.rds)]
ulist <- list()
for (i in 1:length(all.rds)){
  print(i)
  
  #read the datafile for the site
  dat <- readRDS(all.rds[i])
  is.Q <- !is.na(dat$MeanDischarge_m3persec)
  counts <- apply(dat,MARGIN=2, FUN = function(x) { sum(!is.na(x) & is.Q)})
  pos.analytes <- grep("[|]",names(counts))
  counts <- counts[pos.analytes]
  pos <- grep("Water",names(counts))
  if (length(pos)>0) counts <- counts[-pos]
  pos <- grep("Start",names(counts))
  if (length(pos)>0) counts <- counts[-pos]
  pos <- grep("Lab",names(counts))
  if (length(pos)>0) counts <- counts[-pos]
  counts <- counts[counts>50]
  ulist[[all.rds[i]]] <- counts
}

defs <- c("PTM","Post Transition Metals",
          "TM", "Transition Metals",
          "MI","Major Ions",
          "ML", "Metaloids",
          "AEM","Alkaline Earth Metals",
          "ONM", "Other Nin Metals",
          "B","Biological",
          "H","Halogens",
          "E", "Earths")
unique.names <- unique(sort(unlist(lapply(ulist, FUN = function(x) names(x)))))
nameDB <- matrix(c("Acidity (tot) (CaCO3) | mg/L"        , "MI"          
                   , "Ag (tot) | mg/L"                     , "TM"         
                   , "Al (sol) | mg/L"                     , "PTM"          
                   , "Al (tot) | mg/L"                     , "PTM"         
                   , "Alkalinity (CO3-CO3) | mg/L"         , "MI"           
                   , "Alkalinity (HCO3-CaCO3) | mg/L"      , "MI" 
                   , "Alkalinity (HCO3-HCO3) | mg/L"       , "MI"  
                   , "Alkalinity (tot) (CaCO3) | mg/L"     , "MI" 
                   , "As (tot) | mg/L"                     , "ML"  
                   , "B (tot) | mg/L"                      , "ML"  
                   , "Ba (tot) | mg/L"                     , "AEM"   
                   , "C (sol inorg) {DIC} | mg/L"          , "ONM" 
                   , "C (sol org) {DOC, DOC as NPOC} | mg/L", "ONM" 
                   , "C (tot org) {TOC, TOC as NPOC} | mg/L", "ONM"
                   , "Ca (sol) | mg/L"                      , "AEM" 
                   , "Ca (tot) | mg/L"                      , "AEM"
                   , "Cd (tot) | mg/L"                      , "TM" 
                   , "Chlorophyll a (by vol) | mg/L"        ,"B"
                   , "Chlorophyll b (by vol) | mg/L"        ,"B" 
                   , "Chlorophyll c (by vol) | mg/L"        ,"B"
                   , "Chlorophyll sample volume | mL"       ,"B"          
                   , "Cl (sol) | mg/L"                      ,"H"                      
                   , "Co (tot) | mg/L"                      ,"TM" 
                   , "Cond calc 25 deg C | uS/cm"           ,"W" 
                   , "Cond comp 25 deg C (in situ) | uS/cm" ,"W"
                   , "Cond comp 25 deg C (lab) | uS/cm"     ,"W"
                   , "Cond uncomp (in situ) | uS/cm"        ,"W" 
                   , "Cond uncomp (lab) | uS/cm"            ,"TM" 
                   , "Cr (tot) | mg/L"                      ,"TM" 
                   , "Cu (tot) | mg/L"                      ,"TM"    
                   , "F (sol) | mg/L"                       ,"H"  
                   , "Fe (sol) | mg/L"                      ,"TM" 
                   , "Fe (tot) | mg/L"                      ,"TM"                
                   , "Hardness (tot) (CaCO3) {Ca+Mg} | mg/L", "MI"
                   , "K (sol) | mg/L"                       ,"AM" 
                   , "K (tot) | mg/L"                       ,"AM"
                   , "Mg (sol) | mg/L"                      ,"AEM"
                   , "Mg (tot) | mg/L"                      ,"AEM" 
                   , "Mn (sol) | mg/L"                      ,"TM"
                   , "Mn (tot) | mg/L"                      ,"TM" 
                   , "Mo (tot) | mg/L"                      ,"TM"
                   , "N (sum sol org) {DON} | mg/L"         ,"ONM" 
                   , "N (sum sol ox) {NOx-N, TON} | mg/L"   ,"ONM" 
                   , "N (tot kjel) {TKN} | mg/L"            ,"ONM"  
                   , "N (tot org) {TON} | mg/L"             ,"ONM" 
                   , "N (tot sol) {TN-filt} | mg/L"         ,"ONM"  
                   , "N (tot) {TN, pTN} | mg/L"             ,"ONM" 
                   , "Na (sol) | mg/L"                      ,"AM"  
                   , "Na (tot) | mg/L"                      ,"AM" 
                   , "NH3-N/NH4-N (sol) | mg/L"             ,"ONM" 
                   , "Ni (tot) | mg/L"                      ,"TM"
                   , "NO2-N (sol) | mg/L"                   ,"ONM" 
                   , "NO3-N (sol) | mg/L"                   ,"ONM"
                   , "NO3 (sol) | mg/L"                     ,"ONM" 
                   , "O - DO % | %"                         ,"ONM" 
                   , "O - DO (in situ) | mg/L"              ,"ONM"
                   , "O - DO | mg/L"                        ,"ONM" 
                   , "P (sol) | mg/L"                       ,"ONM"
                   , "P (tot org) {TOP} | mg/L"             ,"ONM" 
                   , "P (tot) {TP, pTP} | mg/L"             ,"ONM"
                   , "Pb (tot) | mg/L"                      ,"TM" 
                   , "pH (in situ) | no units"              ,"MI"
                   , "pH | no units"                       ,"MI"  
                   , "PO4-P (sol react) {SRP, FRP} | mg/L"  ,"ONM" 
                   , "S (tot) | mg/L"                       ,"ONM"
                   , "Salinity | mg/L"                      ,"MI" 
                   , "Sb (tot) | mg/L"                      ,"ML"
                   , "Se (tot) | mg/L"                      ,"ONM" 
                   , "Si (tot) | mg/L"                      ,"ML"
                   , "SiO2-Si (sol react) | mg/L"           ,"ML" 
                   , "SiO2 (sol react) | mg/L"              ,"ML"
                   , "Sn (tot) | mg/L"                      ,"PTM" 
                   , "SO4 (sol) | mg/L"                     ,"ONM"
                   , "SO4 (tot) | mg/L"                     ,"ONM" 
                   , "Sr (tot) | mg/L"                      ,"AEM"
                   , "Suspended solids (EDI) | mg/L"        ,"E" 
                   , "Suspended solids (gulp) | mg/L"       ,"E"
                   , "Suspended solids (pump) | mg/L"       ,"E" 
                   , "Suspended solids | mg/L"              ,"E"
                   , "Suspended solids <63u (EDI) | mg/L"   ,"E" 
                   , "Suspended solids <63u (gulp) | mg/L"  ,"E"
                   , "Suspended solids <63u (pump) | mg/L"  ,"E" 
                   , "Suspended solids >63u (EDI) | mg/L"   ,"E"
                   , "Suspended solids >63u (gulp) | mg/L"  ,"E" 
                   ,"Suspended solids >63u (pump) | mg/L"  ,"E"
                   , "TDSalts (sum of ions) | mg/L"         ,"MI" 
                   , "TDSolids (calc @180°C-by cond) | mg/L","E"
                   , "Ti (tot) | mg/L"                     ,"TM"  
                   , "TSS | mg/L"                          ,"MI" 
                   , "Turbidity (JTU) | JTU"               ,"E"
                   , "Turbidity (NTU) | NTU"               ,"E" 
                   , "V (tot) | mg/L"                       ,"TM"       
                   , "Zn (sol) | mg/L"                     ,"TM" 
                   , "Zn (tot) | mg/L"                     ,"TM"),
                 byrow = TRUE,ncol=2)

soluteSiteList <- vector("list",length=length(files))
for (j in 1:length(files)){
  print(files[j])
  solutes <- c("")
  if (file.exists(files[j])){ #Does a file exist
    print(file.exists(files[j]))
    CQ <- readRDS(file=files[j])  #may need a check here for any data
    siteref <- strsplit(files[j],split = "CQ_")[[1]][2]
    varsInUlist <- ulist[[paste0("CQ_",siteref)]]
    nmConc <- names(CQ)[names(CQ) %in% nameDB[,1]] #names(varsInUlist)
    dbColumns <- which(nameDB %in% nmConc)
    #length(nmConc)
    
    for (i in 1:length(nmConc)){
      #dates <- as.Date(substr(CQ$`Collect Date (Known Accuracy)`,1,11), format = "%d-%b-%Y")
      pos <- !is.na(as.numeric(c(CQ[,nmConc[i]])[[1]])) & is.finite(log10(as.numeric(CQ$MeanDischarge_m3persec) ))
      if (sum(pos)>50){
        #print(sum(pos))
       # meddate <- median(as.numeric(dates[pos ]))
        #print(dates[pos ][which(as.numeric(dates[pos ]) == meddate)])
        solutes <- c(solutes, nmConc[i])
      }
    }
  }
    soluteSiteList[[j]] <- solutes 
}

      
uniqueChems <- sort(unique(unlist(soluteSiteList)))
df <- data.frame(file = files)
df2 <- matrix(NA,nrow=length(files),ncol= length(uniqueChems)-1)
colnames(df2) <- uniqueChems[-1]
df3 <- cbind(df,df2)

for (j in 1:length(files)){
  n <- length(soluteSiteList[[j]])
  if (n>1){
    for (i in 2:n){
      pos <- which(colnames(df3) == soluteSiteList[[j]][i])
      df3[j,pos] <- 1
    }
  }
}



siteSolutes <- read.csv("D:\\WIRT\\Processed_CQDiscrete\\SiteSolutes.csv", skip = 1)

sites2process <- which(siteSolutes$Number>0)
bayesResults <- list()
library(stringr)
for (i in 1:length(sites2process)){
  print(paste0(siteSolutes[sites2process[i],3]))
  chems2process <- names(siteSolutes)[4:dim(siteSolutes)[2]][which(!is.na(siteSolutes[sites2process[i],4:dim(siteSolutes)[2]]) )]
  CQ <- readRDS(file=as.character(siteSolutes[sites2process[i],3]))
  names(CQ) <- str_replace_all(names(CQ) , "[[:punct:]]", ".")
  names(CQ) <- str_replace_all(names(CQ) , "\\| ", ".")
  names(CQ) <- gsub(" ",".",names(CQ))
  #names(CQ) <- gsub("(\\+|\\%|\\,|\\(|\\)|\\||\\/|\\{|\\}| > | < |\\-)",".",names(CQ))
  names(CQ) <- gsub("\\.\\.",".", names(CQ))  
  names(CQ) <- gsub("\\.\\.",".", names(CQ))  
  names(CQ) <- gsub("<",".", names(CQ)) 
  names(CQ) <- gsub(">",".", names(CQ))
  names(CQ) <- gsub("\\+",".", names(CQ))
  names(CQ) <- gsub("\\-",".", names(CQ))
  names(CQ) <- gsub("°",".", names(CQ))
  names(CQ)  <- gsub("\\.\\.",".", names(CQ) )  
  
  siteres <- list()
  for (j in chems2process){
    j2 <- gsub("\\.\\.",".",j)  
    j2 <- gsub("\\.\\.",".",j2)  
    j2 <- gsub("\\.\\.",".",j2) 
    #j2 <- paste0(gsub("\\.\\.",".",j2),".")
    print(paste(j2))
    concs <- CQ[,j2]
    pos.lor <- grep("<",concs)
    concs <- gsub("<","",concs)
    concs <- as.numeric(concs)
    if (length(pos.lor)>0) concs[pos.lor] <- concs[pos.lor]/2
    if(length(grep("pH",j2))>0){ 
      logC <- concs
       logC <- ifelse(logC > 14 | logC< 0, logC, NA)
    } else {
      logC <- log10(concs)
    }
    logQ <- log10(CQ[,"MeanDischarge.m3persec"])
    dates <- as.Date(substr(CQ[,"Collect.Date.Known.Accuracy."],1,11), format = "%d-%b-%Y")
    data <- list(x = logQ,
                 y = logC,
                 id = CQ$`Site.Ref`[1], 
                 nms = c("log10 Q","log10 C"),
                 date = dates
    )
    posx <- !is.finite(data$x)
    posy <- !is.finite(data$y)
    data$x[posx] <- NA
    data$y[posy] <- NA
    
    posxy <- which(!is.na(data$x) & !is.na(data$y))
    if (length(posxy)>100){
      data$x <- data$x[posxy]
      data$y <- data$y[posxy]
      
      if (sd(data$y)!=0){
      med.posxy <-  floor(median(1:length(posxy)))
      
      data$med.date <- data$date[posxy[med.posxy]]
      data$date <- data$date[posxy]
      data$start.date <- data$date[1]
      data$end.date <- data$date[length(data$date)]
      data$period <- (data$date <=  data$med.date)
      data$period[is.na(data$period)] <- ifelse(which(is.na(data$period))< max(which(data$period),na.rm=TRUE), TRUE,FALSE)
      #Whole Time series
      fname1 <- paste0(gsub("\\.","",paste0(CQ$`Site.Ref`[1], "_",j2)),"__All_PW.pdf")
      fname2 <- paste0(gsub("\\.","",paste0(CQ$`Site.Ref`[1], "_",j2)),"_All_LM.pdf")
      pW0 <- piecewiseBayes(data, fname = fname1)
      lm0 <- lmBayes(data, fname = fname2)
      
      #T1
      datat1 <- data
      datat1$x <- datat1$x[datat1$period]
      datat1$y <- datat1$y[datat1$period]
      datat1$date <- datat1$date[datat1$period]
      fname1 <- paste0(gsub("\\.","",paste0(CQ$`Site.Ref`[1], "_",j2)),"_T1_PW.pdf")
      fname2 <- paste0(gsub("\\.","",paste0(CQ$`Site.Ref`[1], "_",j2)),"_T1_LM.pdf")
      pW1 <- piecewiseBayes(datat1, fname = fname1)
      lm1 <- lmBayes(datat1, fname = fname2)
      
      #T2
      datat2 <- data
      datat2$x <- datat2$x[!data$period]
      datat2$y <- datat2$y[!data$period]
      datat2$date <- datat2$date[!data$period]
      fname1 <- paste0(gsub("\\.","",paste0(CQ$`Site.Ref`[1], "_",j2)),"_T2_PW.pdf")
      fname2 <- paste0(gsub("\\.","",paste0(CQ$`Site.Ref`[1], "_",j2)),"_T2_LM.pdf")
      pW2 <- piecewiseBayes(datat2, fname = fname1)
      lm2 <- lmBayes(datat2, fname = fname2)
      
      
      siteres[[j2]] <- list(pw0 = pW0, lm0 = lm0, 
                            pw1 = pW1, lm1 = lm1,
                            pw2 = pW2, lm2 = lm2,
                            data = data)
      }
    }
  }
  bayesResults[[CQ$`Site.Ref`[1]]] <- siteres
}
saveRDS(bayesResults,"bayesResults.RDS")

nmBayes <- names(bayesResults)
for (i in 1:length(bayesResults)) {
  print(nmBayes[i])
  res <- bayesResults[[nmBayes[i]]]
  for (j in names(res)){
  
    res2 <- res[[j]]
    data <- res2$data
    CVCQ <- matrix(c(
            mean(data$y),
            mean(data$x),
            sd(data$y)/mean(data$y),
            sd(data$x)/mean(data$x),
            sd(data$y)/mean(data$y) / (sd(data$x)/mean(data$x)),
            range(data$y),
            range(data$x),
            
            mean(data$y[data$period]),
            mean(data$x[data$period]),
            sd(data$y[data$period])/mean(data$y[data$period]),
            sd(data$x[data$period])/mean(data$x[data$period]),
            sd(data$y[data$period])/mean(data$y[data$period]) / ( sd(data$x[data$period])/mean(data$x[data$period])),
            range(data$y[data$period]),
            range(data$x[data$period]),
            
            mean(data$y[!data$period]),
            mean(data$x[!data$period]),
            sd(data$y[!data$period])/mean(data$y[!data$period]),
            sd(data$x[!data$period])/mean(data$x[!data$period]),
            sd(data$y[!data$period])/mean(data$y[!data$period]) / ( sd(data$x[!data$period])/mean(data$x[!data$period])),
            range(data$y[!data$period]),
            range(data$x[!data$period])
            
  ), nrow=1)
  colnames(CVCQ) <- c("meanlogC_all","meanlogQ_all", "CVC_all","CVQ_all","CVCCVQ_all","minC_all","maxC_all","minQ_all","maxQ_all",
                      "meanlogC_T1","meanlogQ_T1", "CVC_T1","CVQ_T1","CVCCVQ_T1","minC_T1","maxC_T1","minQ_T1","maxQ_T1",
                      "meanlogC_T2","meanlogQ_T2", "CVC_T2","CVQ_T2","CVCCVQ_T2","minC_T2","maxC_T2","minQ_T2","maxQ_T2")
  res2$CVCQ <- CVCQ
  bayesResults[[nmBayes[i]]][[j]] <- res2
  }
}
saveRDS(bayesResults,"bayesResults.RDS")            


setwd("D:\\WIRT\\Processed_CQDiscrete")
bayesResults <- readRDS("bayesResults.RDS")
#transform to a data frame
df.names <- c("site.id", "solute", "start.date","end.date","med.date",
              "n1","n2","nall",
              "pwall_a_025","pwall_a_500", "pwall_a_975", 
              "pwall_b1_025","pwall_b1_500", "pwall_b1_975",
              "pwall_delta_025","pwall_delta_500","pwall_delta_975", 
              "pwall_x0_025","pwall_x0_500","pwall_x0_975", 
              "pwall_sigma_025","pwall_sigma_500","pwall_sigma_975", 
              "pwall_b2_025","pwall_b2_500", "pwall_b2_975",
              "pwall_DIC", "pwall_IC", "pwall_pD", "pwall_pV",  "pwall_Dbar", 
              "pwall_Dhat", "pwall_WAIC1", "pwall_WAIC2", "pwall_lppd", "pwall_pWAIC1", "pwall_pWAIC2",
              "lmall_a_025", "lmall_a_500", "lmall_a_975", 
              "lmall_b_025", "lmall_b_500", "lmall_b_975",
              "lmall_sigma_025","lmall_sigma_500","lmall_sigma_975", 
              "lmall_DIC", "lmall_IC", "lmall_pD", "lmall_pV",  "lmall_Dbar", "lmall_Dhat", "lmall_WAIC1", "lmall_WAIC2", "lmall_lppd", "lmall_pWAIC1", "lmall_WAIC2",
              
              "pwt1_a_025","pwt1_a_500", "pwt1_a_975", 
              "pwt1_b1_025","pwt1_b1_500", "pwt1_b1_975",
              "pwt1_delta_025","pwt1_delta_500","pwt1_delta_975", 
              "pwt1_x0_025","pwt1_x0_500","pwt1_x0_975", 
              "pwt1_sigma_025","pwt1_sigma_500","pwt1_sigma_975", 
              "pwt1_b2_025","pwt1_b2_500", "pwt1_b2_975",
              "pwt1_DIC", "pwt1_IC", "pwt1_pD", "pwt1_pV",  "pwt1_Dbar", 
              "pwt1_Dhat", "pwt1_WAIC1", "pwt1_WAIC2", "pwt1_lppd", "pwt1_pWAIC1", "pwt1_pWAIC2",
              "lmt1_a_025", "lmt1_a_500", "lmt1_a_975", 
              "lmt1_b_025", "lmt1_b_500", "lmt1_b_975",
              "lmt1_sigma_025","lmt1_sigma_500","lmt1_sigma_975", 
              "lmt1_DIC", "lmt1_IC", "lmt1_pD", "lmt1_pV",  "lmt1_Dbar", "lmt1_Dhat", "lmt1_WAIC1", "lmt1_WAIC2", "lmt1_lppd", "lmt1_pWAIC1", "lmt1_WAIC2",
              
              "pwt2_a_025","pwt2_a_500", "pwt2_a_975", 
              "pwt2_b1_025","pwt2_b1_500", "pwt2_b1_975",
              "pwt2_delta_025","pwt2_delta_500","pwt2_delta_975", 
              "pwt2_x0_025","pwt2_x0_500","pwt2_x0_975", 
              "pwt2_sigma_025","pwt2_sigma_500","pwt2_sigma975", 
              "pwt2_b2_025","pwt2_b2_500", "pwt2_b2_975",
              "pwt2_DIC", "pwt2_IC", "pwt2_pD", "pwt2_pV",  "pwt2_Dbar", 
              "pwt2_Dhat", "pwt2_WAIC1", "pwt2_WAIC2", "pwt2_lppd", "pwt2_pWAIC1", "pwt2_pWAIC2",
              "lmt2_a_025", "lmt2_a_500", "lmt2_a_975", 
              "lmt2_b_025", "lmt2_b_500", "lmt2_b_975",
              "lmt2_sigma_025","lmt2_sigma_500","lmt2_sigma_975", 
              "lmt2_DIC", "lmt2_IC", "lmt2_pD", "lmt2_pV",  "lmt2_Dbar", "lmt2_Dhat", "lmt2_WAIC1", "lmt2_WAIC2", "lmt2_lppd", "lmt2_pWAIC1", "lmt2_WAIC2",
              
              "meanlogC_all","meanlogQ_all", "CVC_all","CVQ_all","CVCCVQ_all","minC_all","maxC_all","minQ_all","maxQ_all",
              "meanlogC_t1","meanlogQ_t1", "CVC_t1","CVQ_t1","CVCCVQ_t1","minC_t1","maxC_t1","minQ_t1","maxQ_t1",
              "meanlogC_t2","meanlogQ_t2", "CVC_t2","CVQ_t2","CVCCVQ_t2","minC_t2","maxC_t2","minQ_t2","maxQ_t2")
nrows <-  sum(unlist(lapply(bayesResults, FUN = function(x) length(x))))
ncols <- length(df.names)
df <- matrix(NA, nrow=nrows,ncol=ncols)
colnames(df) <- df.names
df <- data.frame(df, stringsAsFactors = FALSE)
counter <- 1
nmBayes <- names(bayesResults)
for (i in 1:length(nmBayes)){
  site.id <- nmBayes[i]
  results <- bayesResults[[site.id]]
  solutes <- names(results)
  for (j in 1:length(solutes)){
    print(paste(i,j))
    df[counter,"site.id"] <- site.id
    df[counter,"solute"] <- solutes[j]
    df[counter,"start.date"] <- as.character(results[[solutes[[j]]]]$data$start.date)
    df[counter,"end.date"] <- as.character(results[[solutes[[j]]]]$data$end.date)
    df[counter,"med.date"] <- as.character(results[[solutes[[j]]]]$data$med.date)
    df[counter,"n1"] <- as.character(results[[solutes[[j]]]]$data$med.date)
    df[counter,"n2"] <- as.character(results[[solutes[[j]]]]$data$med.date)
    df[counter,"nall"] <- as.character(results[[solutes[[j]]]]$data$med.date)
    
    
    df[counter,9:dim(df)[2]] <- c(
      unlist(results[[solutes[[j]]]][["pw0"]]),
      unlist(results[[solutes[[j]]]][["lm0"]]),
      unlist(results[[solutes[[j]]]][["pw1"]]),
      unlist(results[[solutes[[j]]]][["lm1"]]),
      unlist(results[[solutes[[j]]]][["pw2"]]),
      unlist(results[[solutes[[j]]]][["lm0"]]),
      unlist(results[[solutes[[j]]]][["CVCQ"]])
    )
    counter <- counter + 1
  }
}

#edit to fix up slight differences in solute names due to stop start above
newsols <- sub("\\.$","",df$solute)
df$solute <- newsols
saveRDS(df,file = "BayesianFits_WholeSeries.RDS")


setwd("D:\\WIRT\\Processed_CQDiscrete")
bf <- readRDS("BayesianFits_WholeSeries.RDS")



#plot of times
cqtimes <- bf[,c("start.date","med.date","end.date")]
cqtimes$start.date <- as.Date(cqtimes$start.date)
cqtimes$med.date <- as.Date(cqtimes$med.date)
cqtimes$end.date <- as.Date(cqtimes$end.date)

cqtimes <- cqtimes[sort(as.numeric(cqtimes$start.date),
                        index.return = TRUE)$ix,]

mindate <- which(cqtimes$start.date == min(cqtimes$start.date))

cqtimes2 <- order(cqtimes$start.date,cqtimes$med.date,cqtimes$end.date) 
cqtimes <- cqtimes[cqtimes2,]
cqtimes$index <- 1:length(cqtimes$start.date)
tmin <- min(cqtimes[,1])
tmax <- max(cqtimes[,3])
library(zoo)
dummy.zoo <- zoo(c(-1,-1),seq.Date(tmin,tmax,length.out = 2))
pdf(file = "SamplingDistribution.pdf",height = 8,width = 6)
par(oma=c(3,3,0.5,1),mar=c(0,0,0,0))
layout(cbind(c(2,1)),width = 1,heights =c(1,3))
plot(dummy.zoo,type="n",ylim=c(1,length(cqtimes$start.date)),
     xlab="",ylab="",yaxt = "n", yax.flip=TRUE)
axis(side=2,las=1,at=seq(0,1000,by=200))
for(i in 1:dim(cqtimes)[1]){
  lines(zoo(c(i,i),seq.Date(cqtimes[i,1],cqtimes[i,3],length.out = 2)))
  
}
for(i in 1:dim(cqtimes)[1]){
  points(zoo(c(i),cqtimes[i,2]),pch=19,col="#33A02C",cex=0.5)
}

sample.density <- seq.Date(as.Date("1950-01-01"),as.Date("2019-01-01"),by="year")
ntsperyr <- vector("numeric",length(sample.density))
for (i in 1:length(sample.density)){
  ntsperyr[i] <- length(which(cqtimes[,1] < sample.density[i] &
                         cqtimes[,3] > sample.density[i]))
}
samden <- zoo(ntsperyr,sample.density)
plot(samden,xaxt="n",xlab="",ylab = "",lwd=2,yaxt="n")
axis(side =2, at=seq(0,800,by=200),las=1)
dev.off()

#Classify all as PL or BSPL
#DIC: The idea is that models with smaller DIC should be preferred to models with larger DIC. Models are penalized both by the value of D ¯ {\displaystyle {\bar {D}}} \bar{D}, which favors a good fit, but also (similar to AIC) by the effective number of parameters p D {\displaystyle p_{D}} p_{D}
# IC: the model with the lowest BIC is preferred

#function input
#  dataframe with column numbers
#Step 1: Compare DIC (or IC) of a linear to piecewise model
#        Test1 = Piecewise Model Preferred = TRUE
#Setp 2: Check if in the piecewise model the prediction interval for delta spans zero
#        Test2 = delta not in zero = TRUE
#Step 3: Decide whether piecewise or linear
#        Test3 = Test1 and Test2
#Step 4: Classify the piecewise model type
          #E2 = pwb1025
          #F2 = pwb1975
#         Slope1_Type =IF(AND(E2<0,F2<0),"D",IF(AND(E2<0,F2>0),"C","A"))
#         AQ2 = Test3
#         AS2 = Slope1_Type
#         G2 = pwb2025
#         H2 = pwb2975
#         Slope2_Type ==IF(AQ2,IF((AVERAGE(E2:F2)-AVERAGE(G2:H2))<0,"D","A"),AS2)
#         PiecewiseModel_Type = =CONCATENATE(Slope1_Type, Slope2_Type)
#Step 5: Subclassify by types
#        Av2 = PiecewiseModel_Type == "AA"
#        subtype_AA=IF(AV2,IF(AND(G2<0,H2<0),"2",IF(AND(G2>0,H2>0),"1","")),"")
#        AX2 = PiecewiseModel_Type == "DD"
#        subtype_DD =IF(AX2,IF(AND(G2<0,G2<0),"1",IF(AND(G2>0,H2>0),"2","")),"")
#        Full_PWModel = CONCATENATE(PiecewiseModel_Type,subtype_AA,subtype_DD)
#Step 6: Linear Model Type
#         Z2 = lmb1025
#         AA2 = lmb1975
#        LinearModelType = IF(AND(Z2<0,AA2>0), "C",IF(AVERAGE(Z2:AA2)<0,"D","A"))
#Step 7: Final Model Classification
#         Model_Type =IF(Test3,PiecewiseModel_Type,LinearModelType)
#Step 8: Get Model Parameters
#         If Linear Model
#             lm_x500, lm_b500
#         If PW Model
#             pw_a500, pw_x500, pw_delta500, pw_b500

CQModelold <- function(df){
  #c(df$lm_DIC, df$pw_DIC, 
  #df$delta025, df$delta975, 
  #df$pwb025, df$pwb975,
  #df$lmb025,df$lmb975
  #df$lm_x500, df$lm_b500, 
  #df$pw_a500, df$pw_x500, df$pw_delta500, df$pw_b500)
  
  #Step 1: Compare DIC (or IC) of a linear to piecewise model
  #        Test1 = Piecewise Model Preferred = TRUE
  Test1 <- df$lm_DIC > df$pw_DIC
  #Setp 2: Check if in the piecewise model the prediction interval for delta spans zero
  #        Test2 = delta not in zero = TRUE
  Test2 <- !(df$pw_delta_025 < 0 & df$pw_delta_975 > 0)
  
  #Step 3: Decide whether piecewise or linear
  #        Test3 = Test1 and Test2
  Test3 <- Test1 & Test2
  
  
  #PW tests
  #Step 4: Classify the piecewise model type
  #E2 = pwb1025
  #F2 = pwb1975
  E2 <- df$pw_b1_025
  F2 <- df$pw_b1_975
  #         Slope1_Type =IF(AND(E2<0,F2<0),"D",IF(AND(E2<0,F2>0),"C","A"))
  Slope1_Type <- ifelse(E2<0 & F2<0,"D",ifelse(E2<0 & F2>0,"C","A"))
  #         AQ2 = Test3
  #         AS2 = Slope1_Type
  #         G2 = pwb2025
  #         H2 = pwb2975
  G2 <- df$pw_delta_025
  H2 <- df$pw_delta_975
  pw_d <- df$pw_delta_500
  pw_b <- df$pw_b1_500
  #         Slope2_Type ==IF(AQ2,IF((AVERAGE(E2:F2)-AVERAGE(G2:H2))<0,"D","A"),AS2)
  Slope2_Type <- ifelse(Test3,ifelse(df$pw_b2_500>0,"A","D"),Slope1_Type)
  #         PiecewiseModel_Type = =CONCATENATE(Slope1_Type, Slope2_Type)
  PiecewiseModel_Type <- paste0(Slope1_Type,Slope2_Type)
  #Step 5: Subclassify by types
  #        Av2 = PiecewiseModel_Type == "AA"
  AV2 <- PiecewiseModel_Type == "AA"
  #        subtype_AA=IF(AV2,IF(AND(G2<0,H2<0),"2",IF(AND(G2>0,H2>0),"1","")),"")
  subtype_AA <- ifelse(AV2,ifelse(G2<0 & H2<0,"2",ifelse(G2>0 & H2>0,"1","")),"")
  #        AX2 = PiecewiseModel_Type == "DD"
  AX2 <- PiecewiseModel_Type == "DD"
  #        subtype_DD =IF(AX2,IF(AND(G2<0,G2<0),"1",IF(AND(G2>0,H2>0),"2","")),"")
  subtype_DD <- ifelse(AX2,ifelse(G2<0 & H2<0,"1",ifelse(G2>0 & H2>0,"2","")),"")
  #        Full_PWModel = CONCATENATE(PiecewiseModel_Type,subtype_AA,subtype_DD)
  Full_PWModel <- paste0(PiecewiseModel_Type,subtype_AA,subtype_DD)
  #Step 6: Linear Model Type
  #         Z2 = lmb1025
  
  
  Z2 <- df$lm_b_025
  #         AA2 = lmb1975
  AA2 <- df$lm_b_975
  #        LinearModelType = IF(AND(Z2<0,AA2>0), "C",IF(AVERAGE(Z2:AA2)<0,"D","A"))
  LinearModelType <- ifelse(Z2<0 & AA2>0,"C",ifelse(Z2<0 & AA2<0,"D","A"))
  #Step 7: Final Model Classification
  #         Model_Type =IF(Test3,PiecewiseModel_Type,LinearModelType)
  Model_Type <- ifelse(Test3,Full_PWModel,LinearModelType)
  #Step 8: Get Model Parameters
  #         If Linear Model
  #             lm_x500, lm_b500
  if(!Test3){
    params <- c(df$lm_x_500, df$lm_b_500)
  } else{
    params <- c(df$pw_a_500, df$pw_x_500, df$pw_delta_500, df$pw_b1_500,df$pw_b2_500)
  }
  #         If PW Model
  #             pw_a500, pw_x500, pw_delta500, pw_b500
  list(mod = Model_Type,pars = params)
}


CQModel <- function(df){
  #c(df$lm_DIC, df$pw_DIC, 
  #df$delta025, df$delta975, 
  #df$pwb025, df$pwb975,
  #df$lmb025,df$lmb975
  #df$lm_x500, df$lm_b500, 
  #df$pw_a500, df$pw_x500, df$pw_delta500, df$pw_b500)
  
  #Step 1: Compare DIC (or IC) of a linear to piecewise model
  # Test1 = Piecewise Model Preferred = TRUE
  Test1 <- df$lm_DIC > df$pw_DIC
  #Setp 2: Check if in the piecewise model the prediction interval for delta spans zero
  # Test2 = delta not in zero = TRUE
  Test2 <- !(df$pw_delta_025 < 0 & df$pw_delta_975 > 0)
  
  #Step 3: Decide whether piecewise or linear
  # Test3 = Test1 and Test2
  Test3 <- Test1 & Test2
  
  
  #PW tests
  Test_b1_0 <- df$pw_b1_025<0 & df$pw_b1_975>0
  Test_b2_0 <- df$pw_b2_025<0 & df$pw_b2_975>0
  Test_b2_gt0 <- df$pw_b2_500 > 0
  Test_b1_gt0 <- df$pw_b1_500 > 0
  Test_delat0 <- df$pw_delta_025 < 0 & df$pw_delta_975 > 0
  Test_delatgt0 <- df$pw_delta_025 > 0 & df$pw_delta_975 > 0
  
  if (Test_b1_0) { 
    if (Test_b2_0){
      model <- "C"
    } else {
      if (Test_b2_gt0) {
        model <- "CA"
      } else {
        model <- "CD"
      }
    }
  } else {
    if (Test_b1_gt0){
      if (Test_b2_0){
        model <- "AC"
      } else {
        if (Test_b2_gt0 ) {
          if (Test_delat0){
            model <- "A"
          } else {
            if (Test_delatgt0) {
              model <- "AA1" 
            } else {
              model <- "AA2"
            }
          }
        } else {
          if (Test_delat0) {
            model <- "A"
          } else {
            model = "AD"
            }
        }
      }
    } else {
      if (Test_b2_0 ) {
        model <- "DC"
      } else {
        if (Test_b2_gt0){
          model <- "DA"
        } else {
          if (Test_delat0){
            model <- "D"
          } else {
            if (Test_delatgt0) {
              model <- "DD2"
            } else {
              model <- "DD1"
            }
          }
        }
      }
    }
  }
 
  if(!Test3){
    params <- c(df$lm_x_500, df$lm_b_500)
    model <- ifelse(df$lm_b_025<0 & df$lm_b_975>0,"C",ifelse(df$lm_b_500>0,"A","D"))
  } else{
    params <- c(df$pw_a_500, df$pw_x_500, df$pw_delta_500, df$pw_b1_500, df$pw_b2_500)
  }
  #         If PW Model
  #             pw_a500, pw_x500, pw_delta500, pw_b500
  list(mod = model,pars = params)
}
#c(df$lm_DIC, df$pw_DIC, 
#df$delta025, df$delta975, 
#df$pwb025, df$pwb975,
#df$lmb025,df$lmb975
#df$lm_x500, df$lm_b500, 
#df$pw_a500, df$pw_x500, df$pw_delta500, df$pw_b500)
df <- bf[,c("lmall_DIC","pwall_DIC",
            "pwall_delta_025","pwall_delta_975",
            "pwall_b1_025","pwall_b1_975",
            "pwall_b2_025","pwall_b2_975",
            "lmall_b_025","lmall_b_975",
            "lmall_a_500","lmall_b_500",
            "pwall_a_500","pwall_x0_500", "pwall_delta_500","pwall_b1_500","pwall_b2_500")]
names(df) <- c("lm_DIC", "pw_DIC", 
  "pw_delta_025", "pw_delta_975", 
  "pw_b1_025", "pw_b1_975",
  "pw_b2_025", "pw_b2_975",
  "lm_b_025","lm_b_975",
  "lm_x_500", "lm_b_500", 
  "pw_a_500", "pw_x_500", "pw_delta_500", "pw_b1_500","pw_b2_500")
  
 allres <- list();
 for (i in 1:dim(df)[1]){
   res <- CQModelold(df[i,])
   res$dates <- bf[i,c("start.date","end.date","med.date")]
   allres[[i]] <- res
 }
 
 #t1
 df <- bf[,c("lmt1_DIC","pwt1_DIC",
             "pwt1_delta_025","pwt1_delta_975",
             "pwt1_b1_025","pwt1_b1_975",
             "pwt1_b2_025","pwt1_b2_975",
             "lmt1_b_025","lmt1_b_975",
             "lmt1_a_500","lmt1_b_500",
             "pwt1_a_500","pwt1_x0_500", "pwt1_delta_500","pwt1_b1_500","pwt1_b2_500")]
 names(df) <- c("lm_DIC", "pw_DIC", 
                "pw_delta_025", "pw_delta_975", 
                "pw_b1_025", "pw_b1_975",
                "pw_b2_025", "pw_b2_975",
                "lm_b_025","lm_b_975",
                "lm_x_500", "lm_b_500", 
                "pw_a_500", "pw_x_500", "pw_delta_500", "pw_b1_500","pw_b2_500")
 
 t1res <- list();
 for (i in 1:dim(df)[1]){
   res <- CQModelold(df[i,])
   res$dates <- bf[i,c("start.date","end.date","med.date")]
   t1res[[i]] <- res
 }
 
 
 #t2
 df <- bf[,c("lmt2_DIC","pwt2_DIC",
             "pwt2_delta_025","pwt2_delta_975",
             "pwt2_b1_025","pwt2_b1_975",
             "pwt2_b2_025","pwt2_b2_975",
             "lmt2_b_025","lmt2_b_975",
             "lmt2_a_500","lmt2_b_500",
             "pwt2_a_500","pwt2_x0_500", "pwt2_delta_500","pwt2_b1_500","pwt2_b2_500")]
 names(df) <- c("lm_DIC", "pw_DIC", 
                "pw_delta_025", "pw_delta_975", 
                "pw_b1_025", "pw_b1_975",
                "pw_b2_025", "pw_b2_975",
                "lm_b_025","lm_b_975",
                "lm_x_500", "lm_b_500", 
                "pw_a_500", "pw_x_500", "pw_delta_500", "pw_b1_500","pw_b2_500")
 
 t2res <- list();
 for (i in 1:dim(df)[1]){
   res <- CQModelold(df[i,])
   res$dates <- bf[i,c("start.date","end.date","med.date")]
   t2res[[i]] <- res
 }
# names(bf)
# [1] "site.id"        "solute"         "start.date"     "end.date"       "med.date"      
# [6] "pwall_a025"     "pwall_a500"     "pwall_a975"     "pwall_b025"     "pwall_b500"    
# [11] "pwall_b975"     "pwall_delta025" "pwall_delta500" "pwall_delta975" "pwall_x0025"   
# [16] "pwall_x0500"    "pwall_x0975"    "pwall_sigma025" "pwall_sigma500" "pw_allsigma975"
# [21] "pwall_DIC"      "pwall_IC"       "pwall_pD"       "pwall_pV"       "pwall_Dbar"    
# [26] "pwall_Dhat"     "pwall_WAIC1"    "pwall_WAIC2"    "pwall_lppd"     "pwall_pWAIC1"  
# [31] "pwall_pWAIC2"   "lmall_a025"     "lmall_a500"     "lmall_a975"     "lmall_b025"    
# [36] "lmall_b500"     "lmall_b975"     "lmall_sigma025" "lmall_sigma500" "lmall_sigma975"
# [41] "lmall_DIC"      "lmall_IC"       "lmall_pD"       "lmall_pV"       "lmall_Dbar"    
# [46] "lmall_Dhat"     "lmall_WAIC1"    "lmall_WAIC2"    "lmall_lppd"     "lmall_pWAIC1"  
# [51] "lmall_WAIC2.1"  "pwt1_a025"      "pwt1_a500"      "pwt1_a975"      "pwt1_b025"     
# [56] "pwt1_b500"      "pwt1_b975"      "pwt1_delta025"  "pwt1_delta500"  "pwt1_delta975" 
# [61] "pwt1_x0025"     "pwt1_x0500"     "pwt1_x0975"     "pwt1_sigma025"  "pwt1_sigma500" 
# [66] "pw_t1sigma975"  "pwt1_DIC"       "pwt1_IC"        "pwt1_pD"        "pwt1_pV"       
# [71] "pwt1_Dbar"      "pwt1_Dhat"      "pwt1_WAIC1"     "pwt1_WAIC2"     "pwt1_lppd"     
# [76] "pwt1_pWAIC1"    "pwt1_pWAIC2"    "lmt1_a025"      "lmt1_a500"      "lmt1_a975"     
# [81] "lmt1_b025"      "lmt1_b500"      "lmt1_b975"      "lmt1_sigma025"  "lmt1_sigma500" 
# [86] "lmt1_sigma975"  "lmt1_DIC"       "lmt1_IC"        "lmt1_pD"        "lmt1_pV"       
# [91] "lmt1_Dbar"      "lmt1_Dhat"      "lmt1_WAIC1"     "lmt1_WAIC2"     "lmt1_lppd"     
# [96] "lmt1_pWAIC1"    "lmt1_WAIC2.1"   "pwt2_a025"      "pwt2_a500"      "pwt2_a975"     
# [101] "pwt2_b025"      "pwt2_b500"      "pwt2_b975"      "pwt2_delta025"  "pwt2_delta500" 
# [106] "pwt2_delta975"  "pwt2_x0025"     "pwt2_x0500"     "pwt2_x0975"     "pwt2_sigma025" 
# [111] "pwt2_sigma500"  "pw_t2sigma975"  "pwt2_DIC"       "pwt2_IC"        "pwt2_pD"       
# [116] "pwt2_pV"        "pwt2_Dbar"      "pwt2_Dhat"      "pwt2_WAIC1"     "pwt2_WAIC2"    
# [121] "pwt2_lppd"      "pwt2_pWAIC1"    "pwt2_pWAIC2"    "lmt2_a025"      "lmt2_a500"     
# [126] "lmt2_a975"      "lmt2_b025"      "lmt2_b500"      "lmt2_b975"      "lmt2_sigma025" 
# [131] "lmt2_sigma500"  "lmt2_sigma975"  "lmt2_DIC"       "lmt2_IC"        "lmt2_pD"       
# [136] "lmt2_pV"        "lmt2_Dbar"      "lmt2_Dhat"      "lmt2_WAIC1"     "lmt2_WAIC2"    
# [141] "lmt2_lppd"      "lmt2_pWAIC1"    "lmt2_WAIC2.1"   "meanlogC_all"   "meanlogQ_all"  
# [146] "CVC_all"        "CVQ_all"        "CVCCVQ_all"     "minC_all"       "maxC_all"      
# [151] "minQ_all"       "maxQ_all"       "meanlogC_t1"    "meanlogQ_t1"    "CVC_t1"        
# [156] "CVQ_t1"         "CVCCVQ_t1"      "minC_t1"        "maxC_t1"        "minQ_t1"       
# [161] "maxQ_t1"        "meanlogC_t2"    "meanlogQ_t2"    "CVC_t2"         "CVQ_t2"        
# [166] "CVCCVQ_t2"      "minC_t2"        "maxC_t2"        "minQ_t2"        "maxQ_t2"



modt1 <- unlist(lapply(t1res,FUN = function(x) x$mod))
modt2 <- unlist(lapply(t2res,FUN = function(x) x$mod))
modall <- unlist(lapply(allres,FUN = function(x) x$mod))
t2days <- unlist(lapply(t1res,FUN = function(x) as.Date(x$dates[1,2]) - as.Date(x$dates[1,3])))
t1days <- unlist(lapply(t1res,FUN = function(x) as.Date(x$dates[1,3]) - as.Date(x$dates[1,1])))
models <- data.frame(mod.t1 = modt1, mod.t2 = modt2, mod.all = modall , t1.days = t1days, t2.days = t2days)

bf <- cbind(bf,data.frame(mod.t1 = modt1, mod.t2 = modt2, mod.all = modall , t1days = t1days, t2days = t2days))
#bf <- cbind(bf,data.frame(t1days = t1days, t2days = t2days))
pos <- which(bf$mod.t1 %in% c("AA1","AA2","AD","DA","CA","AC","DD1","DD2","DC","CD") & bf$mod.t2 %in% c("AA1","AA2","AD","DA","CA","AC","DD1","DD2","DC","CD"))
plot(bf$pwt1_x0_500[pos],bf$pwt2_x0_500[pos]-bf$pwt1_x0_500[pos])
summary(lm((bf$pwt2_x0_500[pos]-bf$pwt1_x0_500[pos])~bf$pwt1_x0_500[pos]))
abline(coefficients(lm((bf$pwt2_x0_500[pos]-bf$pwt1_x0_500[pos])~bf$pwt1_x0_500[pos])))



#Test whether the proportion of chemostatic models was more than  just chance
mods <- c("A","D","C")
pos1 <- which(as.character(bf$mod.t1) %in% mods & as.character(bf$mod.t2) %in% mods )

pos <- which( as.character(bf$mod.t1) %in% mods & as.character(bf$mod.t2) %in% mods & (
  (bf$lmt1_b_500 < 0 & (bf$lmt2_b_500 - bf$lmt1_b_500) > 0 & (bf$lmt2_b_500 - bf$lmt1_b_500) < -2*bf$lmt1_b_500) |
    (bf$lmt1_b_500 >= 0 & (bf$lmt2_b_500 - bf$lmt1_b_500) < 0 & (bf$lmt2_b_500 - bf$lmt1_b_500) > -2*bf$lmt1_b_500)))
length(pos)
length(pos1)

binom.test(length(pos), length(pos1), p = 3/8,
           alternative ="greater",
           conf.level = 0.95)

#same but for segmented models

#first exponent
mods <- c("A","D","C","AA1","AA2","AD","DA","CA","AC","DD1","DD2","DC","CD")
pos1 <- which(bf$mod.t1 %in% mods & bf$mod.t2 %in% mods )

pos <- which( bf$mod.t1 %in% mods & bf$mod.t2 %in% mods & (
  (bf$pwt1_b1_500 < 0 & (bf$pwt2_b1_500 - bf$pwt1_b1_500) > 0 & (bf$pwt2_b1_500 - bf$pwt1_b1_500) < -2*bf$pwt1_b1_500) |
    (bf$pwt1_b1_500 >= 0 & (bf$pwt2_b1_500 - bf$pwt1_b1_500) < 0 & (bf$pwt2_b1_500 - bf$pwt1_b1_500) > -2*bf$pwt1_b1_500)))
length(pos)
length(pos1)

binom.test(length(pos), length(pos1), p = 3/8,
           alternative ="greater",
           conf.level = 0.95)

#both
mods <- c("A","D","C","AA1","AA2","AD","DA","CA","AC","DD1","DD2","DC","CD")
pos1 <- which(bf$mod.t1 %in% mods & bf$mod.t2 %in% mods )

pos <- which( bf$mod.t1 %in% mods & bf$mod.t2 %in% mods & (
  (bf$pwt1_b2_500 < 0 & (bf$pwt2_b2_500 - bf$pwt1_b2_500) > 0 & (bf$pwt2_b2_500 - bf$pwt1_b2_500) < -2*bf$pwt1_b2_500) |
    (bf$pwt1_b2_500 >= 0 & (bf$pwt2_b2_500 - bf$pwt1_b2_500) < 0 & (bf$pwt2_b2_500 - bf$pwt1_b2_500) > -2*bf$pwt1_b2_500)) &
    (
      (bf$pwt1_b1_500 < 0 & (bf$pwt2_b1_500 - bf$pwt1_b1_500) > 0 & (bf$pwt2_b1_500 - bf$pwt1_b1_500) < -2*bf$pwt1_b1_500) |
        (bf$pwt1_b1_500 >= 0 & (bf$pwt2_b1_500 - bf$pwt1_b1_500) < 0 & (bf$pwt2_b1_500 - bf$pwt1_b1_500) > -2*bf$pwt1_b1_500)))
length(pos)
length(pos1)

binom.test(length(pos), length(pos1), p = 3/8,
           alternative ="less",
           conf.level = 0.95)

#second
mods <- c("A","D","C","AA1","AA2","AD","DA","CA","AC","DD1","DD2","DC","CD")
pos1 <- which(bf$mod.t1 %in% mods & bf$mod.t2 %in% mods )

pos <- which( bf$mod.t1 %in% mods & bf$mod.t2 %in% mods & (
  (bf$pwt1_b2_500 < 0 & (bf$pwt2_b2_500 - bf$pwt1_b2_500) > 0 & (bf$pwt2_b2_500 - bf$pwt1_b2_500) < -2*bf$pwt1_b2_500) |
    (bf$pwt1_b2_500 >= 0 & (bf$pwt2_b2_500 - bf$pwt1_b2_500) < 0 & (bf$pwt2_b2_500 - bf$pwt1_b2_500) > -2*bf$pwt1_b2_500)))
length(pos)
length(pos1)

binom.test(length(pos), length(pos1), p = 3/8,
           alternative ="greater",
           conf.level = 0.95)



#change in db with change in meanlog Q
mods <- c("A","D","C")
pos <- which(bf$mod.t1 %in% mods & bf$mod.t2 %in% mods )
plot(c(bf$meanlogQ_t1[pos],bf$meanlogQ_t2[pos]),c(bf$lmt1_b_500[pos],bf$lmt2_b_500[pos]))
summary(lm((bf$lmt2_b_500[pos] - bf$lmt1_b_500[pos])~(bf$meanlogQ_t2[pos]-bf$meanlogQ_t1[pos])))

mod.types = levels(models$mod.t2)
combos <- expand.grid(1:length(mod.types),1:length(mod.types))
adj.list <- data.frame(from = mod.types[combos[,1]],to = mod.types[combos[,2]],value = rep(0,length(combos[,1])))

for (i in 1:nrow(adj.list)){
  adj.list[i,3] = length(which(as.character(models[,1]) == as.character(adj.list[i,1]) & as.character(models[,2]) == as.character(adj.list[i,2]) & 
                                 models$t1.days > 3650 & models$t2.days > 3650))
}
adj.list[,2] <- paste0(as.character(adj.list[,2]), ".t2")
adj.list[,1] <- paste0(as.character(adj.list[,1]), ".t1")

library(circlize)

library(ggsci)
#cols <- rev(add_transparency(
#         pal_simpsons("springfield")(length(mod.types)),  
#         transparency = 0.2))
library(RColorBrewer)
set.seed(13)
cols <- sample(brewer.pal(12, "Paired"), 12,replace = FALSE)[1:length(mod.types)]
#names(cols) <- mod.types
grid.col1 = cols
#names(grid.col1) <- mod.types

order = c(unique(adj.list[,1]),rev(unique(adj.list[,2])))
# "A.t1"   "C.t1"   "D.t1"   "AA1.t1" "AA2.t1" "DD1.t1" "DD2.t1" "CA.t1" 
# "CD.t1"  "AC.t1"  "AD.t1"  "DC.t1"  "DA.t1"  "DA.t2"  "DC.t2"  "AD.t2" 
# "AC.t2"  "CD.t2"  "CA.t2"  "DD2.t2" "DD1.t2" "AA2.t2" "AA1.t2" "D.t2"  
# "C.t2"   "A.t2"

#order = c(c("AA2.t1","A.t1","AA1.t1","AC.t1","AD.t1", "C.t1","CA.t1","CD.t1","DA.t1","DD1.t1","D.t1", "DC.t1", "DD2.t1"),rev(c("AA2.t2","A.t2","AA1.t2","AC.t2","AD.t2", "C.t2","CA.t2","CD.t2","DA.t2","DD1.t2","D.t2", "DC.t2", "DD2.t2")))
#color.order <- sapply(unlist(lapply(strsplit(order,".t"),FUN = function(x) x[1])), FUN = function(y) which(mod.types == y))
gridcols <- c(grid.col1,rev(grid.col1))
names(gridcols) <- order

chordDiagram(adj.list[,c(1,2,3)],order = order, 
             grid.col = gridcols,
             #col = grid.col1,
             small.gap = 3,
             big.gap = 15,
             link.sort = TRUE, link.decreasing = TRUE, 
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(unique(c(adj.list[,1:2])))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

#abline(h = 0, lty = 2, col = "#00000080")
#mtext(side = 1, line = -2, "Time Period 1")
#mtext(side = 3, line = -2, "Time Period 2")
circos.clear()


mod.types = levels(models$mod.t2) 
#c("A","C","D","AA1","AA2","DD1","DD2","CA","CD","AC","AD", "DC", "DA")
combos <- expand.grid(1:length(mod.types),1:length( sort(unique(bf$solute))))
model.solute.list <- data.frame(
  from =  sort(unique(bf$solute))[combos[,2]],
  to = mod.types[combos[,1]],
  value = rep(0,length(combos[,1])))

for (i in 1:nrow(model.solute.list)){
  model.solute.list[i,3] = length(which(as.character(bf$solute) == as.character(model.solute.list[i,1]) & as.character(models[,3]) == as.character( model.solute.list[i,2])))
}

#par(mfrow = c(1,3))
#Nutrients
oldSols <- c("Alkalinity.HCO3.HCO3.mg.L","Alkalinity.tot.CaCO3.mg.L","Ca.sol.mg.L","Cl.sol.mg.L","K.tot.mg.L","Mg.sol.mg.L","SiO2.Si.sol.react.mg.L","SiO2.sol.react.mg.L", "Hardness.tot.CaCO3.Ca.Mg.mg.L","Cond.comp.25.deg.C.in.situ.uS.cm")
newSols <- c("HCO3","CaCO3","Ca","Cl","K","Mg","Si","SiO2","Ca.Mg","EC")
oldNut <- c("PO4.P.sol.react.SRP.FRP.mg.L","P.tot.TP.pTP.mg.L",
            "N.sum.sol.org.DON.mg.L","N.sum.sol.ox.NOx.N.TON.mg.L","N.tot.kjel.TKN.mg.L",
            "N.tot.org.TON.mg.L","N.tot.TN.pTN.mg.L","NH3.N.NH4.N.sol.mg.L", "NO3.N.sol.mg.L","C.sol.org.DOC.DOC.as.NPOC.mg.L","C.tot.org.TOC.TOC.as.NPOC.mg.L")
newNut <- c("PO4.P","TP","DON","NOx","TKN","TON","TN","NH3","N03","DOC","TOC")
order <- c(mod.types,newNut)
pos.nutrients <- which(model.solute.list[,1] %in% oldNut)
pos.tracers <- which(model.solute.list[,1] %in% oldSols)
Nutrient.mat <- data.frame(
  from = as.character(model.solute.list[pos.nutrients,2]),
  to = as.character(model.solute.list[pos.nutrients,1]), 
  value = model.solute.list[pos.nutrients,3] ,
                           stringsAsFactors = FALSE)
for (j in 1:length(oldNut)){
  Nutrient.mat[Nutrient.mat[,2] == oldNut[j],2] <- newNut[j]
}
gridcols <- c(grid.col1,rep("grey",length(newNut)))
names(gridcols) <- c(order)
chordDiagram(Nutrient.mat,order = order, 
             link.sort = TRUE, link.decreasing = FALSE, 
             grid.col = gridcols,
             small.gap = 3,
             big.gap = 15,
             #col = gridcols[1:length(grid.col1)],
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(unique(c(Nutrient.mat[,1:2])))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()

#tracers
oldSols <- c("Alkalinity.HCO3.HCO3.mg.L","Alkalinity.tot.CaCO3.mg.L","Ca.sol.mg.L","Cl.sol.mg.L","K.tot.mg.L","Mg.sol.mg.L","SiO2.Si.sol.react.mg.L","SiO2.sol.react.mg.L", "Hardness.tot.CaCO3.Ca.Mg.mg.L","Cond.comp.25.deg.C.in.situ.uS.cm")
newSols <- c("HCO3","CaCO3","Ca","Cl","K","Mg","Si","SiO2","Ca.Mg","EC")
order = c(mod.types,newSols)
pos.tracers <- which(model.solute.list[,1] %in% oldSols)
Tracer.mat <- data.frame(
  from = as.character(model.solute.list[pos.tracers,2]),
  to = as.character(model.solute.list[pos.tracers,1]), 
  value = model.solute.list[pos.tracers,3] ,
                           stringsAsFactors = FALSE)
for (j in 1:length(oldSols)){
  Tracer.mat[Tracer.mat[,2] == oldSols[j],2] <- newSols[j]
}
gridcols <- c(grid.col1,rep("grey",length(newSols)))
names(gridcols) <- c(order)
chordDiagram(Tracer.mat,order = order, 
             link.sort = TRUE, link.decreasing = FALSE, 
             grid.col = gridcols,
             annotationTrack = "grid", 
             small.gap = 3,
             big.gap = 15,
             preAllocateTracks = list(track.height = max(strwidth(unlist(unique(c(Tracer.mat[,1:2])))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()

#Solids
oldSolids <- c("Suspended.solids.63u.gulp.mg.L","Suspended.solids.63u.pump.mg.L","Suspended.solids.pump.mg.L","TDSalts.sum.of.ions.mg.L","TDSolids.calc.180.C.by.cond.mg.L","TSS.mg.L","Turbidity.NTU.NTU")
newsolids <- c("SS63g","SS63p","SS","TDS","TDScal","TSS","Turb")
order = c(mod.types,newsolids)
pos.sed <- which(model.solute.list[,1] %in% oldSolids)
Solids.mat <- data.frame(
  from = as.character(model.solute.list[pos.sed,2]),
  to = as.character(model.solute.list[pos.sed,1]), 
  value = model.solute.list[pos.sed,3] ,
                           stringsAsFactors = FALSE)
for (j in 1:length(oldSolids)){
  Solids.mat[Solids.mat[,2] == oldSolids[j],2] <- newsolids[j]
}
gridcols <- c(grid.col1,rep("grey",length(newsolids)))
names(gridcols) <- c(order)
chordDiagram(Solids.mat,order = order, 
             link.sort = TRUE, link.decreasing = FALSE, 
             annotationTrack = "grid", 
             grid.col = gridcols,
             small.gap = 3,
             big.gap = 15,
             preAllocateTracks = list(track.height = max(strwidth(unlist(unique(c(Solids.mat[,1:2])))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()



#the db v b plot
#extract all the parameters for the C - C, A - A and D - D patterns
modt1 <- unlist(lapply(t1res,FUN = function(x) x$mod))
modt2 <- unlist(lapply(t2res,FUN = function(x) x$mod))
modall <- unlist(lapply(allres,FUN = function(x) x$mod))
t2days <- unlist(lapply(t1res,FUN = function(x) as.Date(x$dates[1,2]) - as.Date(x$dates[1,3])))
t1days <- unlist(lapply(t1res,FUN = function(x) as.Date(x$dates[1,3]) - as.Date(x$dates[1,1])))
models <- data.frame(mod.t1 = modt1, mod.t2 = modt2, mod.all = modall , t1.days = t1days, t2.days = t2days)
pos.same <- which((modt1  == "A" & modt2 == "A" | 
                  modt1  == "C" & modt2 == "C" | 
                  modt1  == "D" & modt2 == "D")  & t1days > 365*4 & t2days > 365*4)
pars.t1 <- t(sapply(pos.same, FUN = function(x) t1res[[x]]$pars))
pars.t2 <- t(sapply(pos.same, FUN = function(x) t2res[[x]]$pars))
solutes <- bf$solute[pos.same]
db <- data.frame(mod = models$mod.t1[pos.same], 
                 solute = solutes,
                 t1.days = models$t1.days[pos.same],
                 b1 = pars.t1[,2],b2 = pars.t2[,2], deltab = pars.t2[,2] - pars.t1[,2], stringsAsFactors = FALSE)
db$solClass <- ifelse(db$solute %in% oldNut, "Nutrients","Tracers")
db$solClass[db$solute %in% oldSolids] <- "Solids"

library(brms)

b1.0 <-
  brm(data = db, family = gaussian,
      deltab ~ 1 + b1,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 7)

db$solClass <- as.factor(db$solClass)
library(dplyr)
library(ggplot2)
library(ggthemes)
f <-
  fitted(b1.0, newdata = db) %>%  # we can use the same `nd` data from last time
  as_tibble() %>%
  bind_cols(db)

pdf(file = "dbvsb1.pdf", width = 6, height = 5)
par(oma=c(0,0,0,0),mar=c(4,4,1,1))
plot(db$b1,db$deltab, ylim=c(-0.4,0.4),xlim=c(-1,1), type="n", ylab = "", xlab = expression(paste(b[1])), las=1, xaxs="i",yaxs="i")
mtext(side = 2, line = 2.5, las= 1, expression(paste(Delta,"b")))
polygon(c(0,-1,-1,0),c(0,2,0,0),col = "dark grey")
polygon(c(0,1,1,0),c(0,-2,0,0),col = "dark grey")
lines(c(-1,1),c(0,0)) 
lines(c(0,0),c(-1,1)) 
lines(c(1,-1),c(-2,2)) 
box()
points(db$b1[db$solute %in% oldNut],db$deltab[db$solute %in% oldNut], pch=21,  bg=gridcols[3],col="black",cex = models$t1.days[pos.same]/365/8)
points(db$b1[db$solute %in% oldSolids],db$deltab[db$solute %in% oldSolids], pch=21,  col="black",bg=gridcols[6],cex = models$t1.days[pos.same]/365/8)
points(db$b1[db$solute %in% oldSols],db$deltab[db$solute %in% oldSols], pch=21, col="black",bg=gridcols[8],cex = models$t1.days[pos.same]/365/8)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))
}
b1mod <- f$b1
dbmod <- f$Estimate
bd025 <- f$Q2.5
bd975 <- f$Q97.5
pos.sort <- sort(b1mod,index.return = TRUE)$ix

polygon(c(b1mod[pos.sort],rev(b1mod[pos.sort]),b1mod[pos.sort][1]),
        c(bd025[pos.sort],rev(bd975[pos.sort]),bd025[pos.sort][1]),col=add.alpha(gridcols[4],0.8), border = NA)
lines(b1mod[pos.sort], dbmod[pos.sort], lty=2,lwd=2, col = "black")

legend(-0.95,-0,pch=21,pt.cex=c(1,1,1,0.5,1,2),col="black",pt.bg = gridcols[c(3,6,8,1,1,1)],legend = c("Nutrients","Solids","Tracers","4 yr","8 yr","16 yr"), bty="n",box.col=NA)

dev.off()



dt.test2 <- data.frame(dQ = bf$meanlogQ_t2[bf$solute %in% sol2]-bf$meanlogQ_t1[bf$solute %in% sol2], dC = bf$meanlogC_t2[bf$solute %in% sol2]-bf$meanlogC_t1[bf$solute %in% sol2], model = models$mod.all)



ggplot(data = as.data.frame(f),aes(x = b1, color = solClass)) +
  geom_smooth(
    aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
        fill = solClass),
    stat = "identity", 
    alpha = 1/4, size = 1/2) +
  geom_point(aes(y = deltab),
             size = 2) +
  scale_colour_pander() +
  scale_fill_pander() +
  scale_x_continuous("b1", expand = c(0, 0)) +
  ylab("deltab") +
  theme_pander() + 
  theme(text = element_text(family = "Times"),
        legend.position = "none")

b1.1 <-
  brm(data = db, family = gaussian,
      deltab ~ 1 + b1 + solClass,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 7)
f <-
  fitted(b1.1, newdata = db) %>%  # we can use the same `nd` data from last time
  as_tibble() %>%
  bind_cols(db)

ggplot(data = as.data.frame(f),aes(x = b1, color = solClass)) +
  geom_smooth(
    aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
        fill = solClass),
    stat = "identity", 
    alpha = 1/4, size = 1/2) +
  geom_point(aes(y = deltab),
             size = 2/3) +
  scale_colour_pander() +
  scale_fill_pander() +
  scale_x_continuous("b1", expand = c(0, 0)) +
  ylab("deltab") +
  theme_pander() + 
  theme(text = element_text(family = "Times"),
        legend.position = "none") +
  facet_wrap(~solClass)


b1.2 <-
  brm(data = db, family = gaussian,
      deltab ~ 1 + b1 * solClass,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 7)
library(tibble)

f <-
  fitted(b1.2, newdata = db) %>%  # we can use the same `nd` data from last time
  as_tibble() %>%
  bind_cols(db)

ggplot(data = as.data.frame(f),aes(x = b1, color = solClass)) +
  geom_smooth(
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = solClass),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  geom_point(aes(y = deltab),
             size = 1) +
  scale_colour_pander() +
  scale_fill_pander() +
  scale_x_continuous("b1", expand = c(0, 0)) +
  ylab("deltab") +
  theme_pander() + 
  theme(text = element_text(family = "Times"),
        legend.position = "none") +
  facet_wrap(~solClass)

db$t1d <- log10(db$t1.days ) - mean(log10(db$t1.days ))
b1.3 <-
  brm(data = db, family = gaussian,
      deltab ~ 1 + b1 + solClass + t1d,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 7,control = list(adapt_delta = 0.90, max_treedepth = 15))

f <-
  fitted(b1.3, newdata = db) %>%  # we can use the same `nd` data from last time
  as_tibble() %>%
  bind_cols(db)

ggplot(data = as.data.frame(f),aes(x = b1, color = solClass)) +
  geom_smooth(
    aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
        fill = solClass),
    stat = "identity", 
    alpha = 1/4, size = 1/2) +
  geom_point(aes(y = deltab),
             size = 1) +
  scale_colour_pander() +
  scale_fill_pander() +
  scale_x_continuous("b1", expand = c(0, 0)) +
  ylab("deltab") +
  theme_pander() + 
  theme(text = element_text(family = "Times"),
        legend.position = "none") + 
  facet_wrap(~solClass)

plot(marginal_effects(b1.3), points = T)

posterior_summary(b1.3) %>% round(digits = 2)


b1.4 <-
  brm(data = db, family = gaussian,
      deltab ~ 1 + b1 +  t1d,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 7,control = list(adapt_delta = 0.90, max_treedepth = 15))

f <-
  fitted(b1.4, newdata = db) %>%  # we can use the same `nd` data from last time
  as_tibble() %>%
  bind_cols(db)

ggplot(data = as.data.frame(f),aes(x = b1, color = solClass)) +
  geom_smooth(
    aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
        fill = solClass),
    stat = "identity", 
    alpha = 1/4, size = 1/2) +
  geom_point(aes(y = deltab),
             size = 1) +
  scale_colour_pander() +
  scale_fill_pander() +
  scale_x_continuous("b1", expand = c(0, 0)) +
  ylab("deltab") +
  theme_pander() + 
  theme(text = element_text(family = "Times"),
        legend.position = "none") + 
  facet_wrap(~solClass)

plot(marginal_effects(b1.4), points = T)

b1.5 <-
  brm(data = db, family = gaussian,
      deltab ~ 1 + b1 * t1d,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 7,control = list(adapt_delta = 0.90, max_treedepth = 15))

f <-
  fitted(b1.5, newdata = db) %>%  # we can use the same `nd` data from last time
  as_tibble() %>%
  bind_cols(db)

ggplot(data = as.data.frame(f),aes(x = b1, color = solClass)) +
  geom_smooth(
    aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
        fill = solClass),
    stat = "identity", 
    alpha = 1/4, size = 1/2) +
  geom_point(aes(y = deltab),
             size = 1) +
  scale_colour_pander() +
  scale_fill_pander() +
  scale_x_continuous("b1", expand = c(0, 0)) +
  ylab("deltab") +
  theme_pander() + 
  theme(text = element_text(family = "Times"),
        legend.position = "none")

plot(marginal_effects(b1.5), points = T)

b1.5 <- add_criterion(b1.5, c("loo", "waic"), reloo = TRUE)
b1.4 <- add_criterion(b1.4, c("loo", "waic"), reloo = TRUE)
b1.3 <- add_criterion(b1.3, c("loo", "waic"), reloo = TRUE)
b1.2 <- add_criterion(b1.2, c("loo", "waic"), reloo = TRUE)
b1.1 <- add_criterion(b1.1, c("loo", "waic"), reloo = TRUE)
b1.0 <- add_criterion(b1.0, c("loo", "waic"), reloo = TRUE)

loo_compare(b1.5,b1.4,b1.3, b1.2,b1.1,b1.0,
            criterion = "waic")
loo_compare(b1.5,b1.4,b1.3, b1.2,b1.1,b1.0,
            criterion = "loo")

model_weights(b1.5,b1.4,b1.3, b1.2,b1.1,b1.0, 
              weights = "waic") %>% 
  round(digits = 2)
model_weights(b1.5,b1.4,b1.3, b1.2,b1.1,b1.0, 
              weights = "loo") %>% 
  round(digits = 2)





library(dplyr)

db$solClass <- ifelse(db$solute %in% oldNut, "Nutrients","Tracers")
db$solClass[db$solute %in% oldSolids] <- "Solids"
modt1 <- unlist(lapply(t1res,FUN = function(x) x$mod))
modt2 <- unlist(lapply(t2res,FUN = function(x) x$mod))
modall <- unlist(lapply(allres,FUN = function(x) x$mod))
t2days <- unlist(lapply(t1res,FUN = function(x) as.Date(x$dates[1,2]) - as.Date(x$dates[1,3])))
t1days <- unlist(lapply(t1res,FUN = function(x) as.Date(x$dates[1,3]) - as.Date(x$dates[1,1])))
a.t1 <- unlist(lapply(t1res,FUN = function(x) {
  pars <- x$pars
  if (length(x$pars) == 2){
    par <- x$pars[1]
  } else {
    par <- x$pars[1]
  }
  par
  }));
x0.t1 <- unlist(lapply(t1res,FUN = function(x) {
  pars <- x$pars
  if (length(x$pars) == 2){
    par <- NA
  } else {
    par <- x$pars[2]
  }
  par
}));
delta.t1 <- unlist(lapply(t1res,FUN = function(x) {
  pars <- x$pars
  if (length(x$pars) == 2){
    par <- NA
  } else {
    par <- x$pars[3]
  }
  par
}));
b1.t1 <- unlist(lapply(t1res,FUN = function(x) {
  pars <- x$pars
  if (length(x$pars) == 2){
    par <- x$pars[2]
  } else {
    par <- x$pars[4]
  }
  par
}));
b2.t1 <- unlist(lapply(t1res,FUN = function(x) {
  pars <- x$pars
  if (length(x) == 2){
    par <- NA
  } else {
    par <- x$pars[5]
  }
  par
}));
a.t2 <- unlist(lapply(t2res,FUN = function(x) {
  pars <- x$pars
  if (length(x$pars) == 2){
    par <- x$pars[1]
  } else {
    par <- x$pars[1]
  }
  par
}));
x0.t2 <- unlist(lapply(t2res,FUN = function(x) {
  pars <- x$pars
  if (length(x$pars) == 2){
    par <- NA
  } else {
    par <- x$pars[2]
  }
  par
}));
delta.t2 <- unlist(lapply(t2res,FUN = function(x) {
  pars <- x$pars
  if (length(x$pars) == 2){
    par <- NA
  } else {
    par <- x$pars[3]
  }
  par
}));
b1.t2 <- unlist(lapply(t2res,FUN = function(x) {
  pars <- x$pars
  if (length(x$pars) == 2){
    par <- x$pars[2]
  } else {
    par <- x$pars[4]
  }
  par
}));
b2.t2 <- unlist(lapply(t2res,FUN = function(x) {
  pars <- x$pars
  if (length(x) == 2){
    par <- NA
  } else {
    par <- x$pars[5]
  }
  par
}));
a.tall <- unlist(lapply(allres,FUN = function(x) {
  pars <- x$pars
  if (length(x$pars) == 2){
    par <- x$pars[1]
  } else {
    par <- x$pars[1]
  }
  par
}));
x0.tall <- unlist(lapply(allres,FUN = function(x) {
  pars <- x$pars
  if (length(x$pars) == 2){
    par <- NA
  } else {
    par <- x$pars[2]
  }
  par
}));
delta.tall <- unlist(lapply(allres,FUN = function(x) {
  pars <- x$pars
  if (length(x$pars) == 2){
    par <- NA
  } else {
    par <- x$pars[3]
  }
  par
}));
b1.tall <- unlist(lapply(allres,FUN = function(x) {
  pars <- x$pars
  if (length(x$pars) == 2){
    par <- x$pars[2]
  } else {
    par <- x$pars[4]
  }
  par
}));
b2.tall <- unlist(lapply(allres,FUN = function(x) {
  pars <- x$pars
  if (length(x) == 2){
    par <- NA
  } else {
    par <- x$pars[5]
  }
  par
}));
models <- data.frame(mod.t1 = modt1, mod.t2 = modt2, mod.all = modall , t1.days = t1days, t2.days = t2days, a.t1 = a.t1, delta.t1 = delta.t1, x0.t1 = x0.t1, b1.t1 = b1.t1, b2.t1 = b2.t1,
                     a.t2 = a.t2, delta.t2 = delta.t2, x0.t2 = x0.t2, b1.t2 = b1.t2, b2.t2 = b2.t2,
                     a.tall = a.tall, delta.tall = delta.tall, x0.tall = x0.tall, b1.tall = b1.tall, b2.tall = b2.tall)
models$solute <- bf$solute
models$solClass <- ifelse(models$solute %in% oldNut, "Nutrients","Tracers")
models$solClass[models$solute %in% oldSolids] <- "Solids"
models$site <- bf$site.id
models$start.date <- bf$start.date
models$end.date <- bf$end.date
models$med.date <- bf$med.date

#add in the hydroclimate
monthly.clim <- readRDS(file = "D:/WIRT/ProcessedClimate/climDBmonthly.RDS")
models$PET.t1 <- rep(NA,nrow(models))
models$PET.t2 <- rep(NA,nrow(models))
models$PET.all <- rep(NA,nrow(models))
models$P.t1 <- rep(NA,nrow(models))
models$P.t2 <- rep(NA,nrow(models))
models$P.all <- rep(NA,nrow(models))
models$Q.t1 <- rep(NA,nrow(models))
models$Q.t2 <- rep(NA,nrow(models))
models$Q.all <- rep(NA,nrow(models))
models$EaonP.t1 <- rep(NA,nrow(models))
models$EaonP.t2 <- rep(NA,nrow(models))
models$EaonP.all <- rep(NA,nrow(models))
models$Aridity.t1 <- rep(NA,nrow(models))
models$Aridity.t2 <- rep(NA,nrow(models))
models$Aridity.all <- rep(NA,nrow(models))
models$dEonE <- rep(NA,nrow(models))
models$dPETonPET <- rep(NA,nrow(models))
models$dQonQ <- rep(NA,nrow(models))
for (i in 1:nrow(models)){
  print(i)
clim <- monthly.clim[models$site[i]][[models$site[[i]]]]
strt.yr <- as.numeric(format(as.Date(models$start.date[1]),"%Y"))
med.yr <-  as.numeric(format(as.Date(models$med.date[1]),"%Y"))
end.yr <-  as.numeric(format(as.Date(models$end.date[1]),"%Y"))
pos.t1.1 <- which(clim$year == strt.yr & clim$month == 3)
pos.t1.2 <- which(clim$year == max(strt.yr + 1,med.yr) & clim$month == 2)
clim$flowyear <- clim$year
clim$flowyear[clim$month<=3] <- clim$flowyear[clim$month<=3]-1

clim$PET <- abs(clim$PET)
Qannual <-  clim %>%
  group_by(flowyear) %>%
  summarise_at(c("PET","Precip","Discharge"),sum)

pos.t2.0 <- which(Qannual$flowyear== strt.yr)
pos.t2.1 <- which(Qannual$flowyear== med.yr)
pos.t2.2 <- which(Qannual$flowyear== end.yr)

PET.t1 <- mean(Qannual$PET[pos.t2.0:pos.t2.1],na.rm = TRUE)
PET.t2 <- mean(Qannual$PET[pos.t2.1:pos.t2.2],na.rm = TRUE)
PET.all <- mean(Qannual$PET[pos.t2.0:pos.t2.2],na.rm = TRUE)
models$PET.t1[i] <- PET.t1
models$PET.t2[i] <- PET.t2
models$PET.all[i] <- PET.all

P.t1 <- mean(Qannual$Precip[pos.t2.0:pos.t2.1],na.rm = TRUE)
P.t2 <- mean(Qannual$Precip[pos.t2.1:pos.t2.2],na.rm = TRUE)
P.all <- mean(Qannual$Precip[pos.t2.0:pos.t2.2],na.rm = TRUE)
models$P.t1[i] <- P.t1
models$P.t2[i] <- P.t2
models$P.all[i] <- P.all

Q.t1 <- mean(Qannual$Discharge[pos.t2.0:pos.t2.1],na.rm = TRUE)
Q.t2 <- mean(Qannual$Discharge[pos.t2.1:pos.t2.2],na.rm = TRUE)
Q.all <- mean(Qannual$Discharge[pos.t2.0:pos.t2.2],na.rm = TRUE)
models$Q.t1[i] <- Q.t1
models$Q.t2[i] <- Q.t2
models$Q.all[i] <- Q.all

EaonP.t1 = (P.t1 - Q.t1)/P.t1
EaonP.t2 = (P.t2 - Q.t2)/P.t2
EaonP.all = (P.all - Q.all)/P.all
models$EaonP.t1[i] <- EaonP.t1
models$EaonP.t2[i] <- EaonP.t2
models$EaonP.all[i] <- EaonP.all

Aridity.t1 <- PET.t1/P.t1
Aridity.t2 <- PET.t2/P.t2
Aridity.all <- PET.all/P.all
models$Aridity.t1[i] <- Aridity.t1
models$Aridity.t2[i] <- Aridity.t2
models$Aridity.all[i] <- Aridity.all

dPonP <- (P.t2 - P.t1)/P.all
dEonE <- ((P.t2 - Q.t2) - (P.t1 - Q.t1))/(P.all - Q.all)
dPETonPET <- ( PET.t2 -  PET.t1)/ PET.all
dQonQ <- (Q.t2 - Q.t1)/Q.all
models$dPonP[i] <- dPonP
models$dEonE[i] <- dEonE
models$dPETonPET[i] <- dPETonPET
models$dQonQ[i] <- dQonQ
}

models2 <- cbind(bf,models[,(1:dim(models)[2])[-which(names(models) %in% names(bf))]])



sol <- c("PO4.P.sol.react.SRP.FRP.mg.L","P.tot.TP.pTP.mg.L",
         "N.sum.sol.org.DON.mg.L","N.sum.sol.ox.NOx.N.TON.mg.L","N.tot.kjel.TKN.mg.L",
         "N.tot.org.TON.mg.L","N.tot.TN.pTN.mg.L","NH3.N.NH4.N.sol.mg.L", "NO3.N.sol.mg.L")
t.test(models2$meanlogC_t2[bf$solute %in% sol], models2$meanlogC_t1[bf$solute %in% sol],"less", paired = TRUE)
t.test(models2$meanlogQ_t2[bf$solute %in% sol], models2$meanlogQ_t1[bf$solute %in% sol],"less", paired = TRUE)

dt.test <- data.frame(dQ = models2$meanlogQ_t2[models2$solute %in% sol]-models2$meanlogQ_t1[models2$solute %in% sol], dC = models2$meanlogC_t2[models2$solute %in% sol]-models2$meanlogC_t1[models2$solute %in% sol], mod = models2$mod.all[models2$solute %in% sol])
dt.test$mod[dt.test$mod %in% c("A","AA1","AA2")] <- "A"
dt.test$mod[dt.test$mod %in% c("D","DD1","DD2")] <- "D"
dt.test$mod <- as.character(dt.test$mod)
dt.test$mod <- as.factor(dt.test$mod)
dt.test <- dt.test[dt.test$mod == "A",]

res <- lm(dC~dQ,data = dt.test)
library(mgcv)
res <- gam(dC~s(dQ),data = dt.test)



b1.x <-
  brm(data = dt.test, family = gaussian,
      dC ~ 1 + dQ,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 7)
f <-
  fitted(b1.x, newdata = dt.test) %>%  # we can use the same `nd` data from last time
  as_tibble() %>%
  bind_cols(dt.test)

plot(f$dQ,f$dC, xlab = "",ylab = "")
mtext(side=2,line=2.5,las=1,text = expression(paste(Delta,"C",sep="")))
mtext(side=1,line=2.5,las=1,text = expression(paste(Delta,"Q",sep="")))
pos.dq <- sort(f$dQ,index.return = TRUE)$ix
lines(f$dQ[pos.dq],f$Estimate[pos.dq])
polygon(c(f$dQ[pos.dq],rev(f$dQ[pos.dq]),f$dQ[pos.dq][1]),c(f$Q2.5[pos.dq],rev(f$Q97.5[pos.dq]),f$Q2.5[pos.dq][1]),col=add.alpha(gridcols[3],0.5))
abline(h=0,lty=2)
abline(v=0,lty=2)


unique(models2$solClass)
sol <- "Nutrients"
dt.test <- data.frame(dQ = models2$meanlogQ_t2[models2$solClass %in% sol]-models2$meanlogQ_t1[models2$solClass %in% sol], dC = models2$meanlogC_t2[models2$solClass %in% sol]-models2$meanlogC_t1[models2$solClass %in% sol], mod = models2$mod.all[models2$solClass %in% sol])
dt.test$mod[dt.test$mod %in% c("A","AA1","AA2")] <- "A"
dt.test$mod[dt.test$mod %in% c("D","DD1","DD2")] <- "D"
dt.test$mod <- as.character(dt.test$mod)
dt.test$mod <- as.factor(dt.test$mod)
#dt.test <- dt.test[dt.test$mod == "D",]

res <- lm(dC~dQ,data = dt.test)
plot(res)
library(mgcv)
res <- gam(dC~s(dQ),data = dt.test)
plot(res)


b1.x <-
  brm(data = dt.test, family = gaussian,
      dC ~ 1 + dQ,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 7)
f <-
  fitted(b1.x, newdata = dt.test) %>%  # we can use the same `nd` data from last time
  as_tibble() %>%
  bind_cols(dt.test)

plot(f$dQ,f$dC, xlab = "",ylab = "",col = gridcols[as.numeric(as.factor(models2$solute[models2$solClass %in% sol]))],pch=19,ylim=c(-1.5,1.5))
abline(h=0,lty=2)
abline(v=0,lty=2)
mtext(side=2,line=2.5,las=1,text = expression(paste(Delta,"C",sep="")))
mtext(side=1,line=2.5,las=1,text = expression(paste(Delta,"Q",sep="")))
pos.dq <- sort(f$dQ,index.return = TRUE)$ix
lines(f$dQ[pos.dq],f$Estimate[pos.dq])
polygon(c(f$dQ[pos.dq],rev(f$dQ[pos.dq]),f$dQ[pos.dq][1]),c(f$Q2.5[pos.dq],rev(f$Q97.5[pos.dq]),f$Q2.5[pos.dq][1]),col=add.alpha(gridcols[4],0.5))
legend("topleft",pch=19,col = gridcols[1:length(levels(as.factor(models2$solute[models2$solClass %in% sol])))],legend = levels(as.factor(models2$solute[models2$solClass %in% sol])),bty="n")




sol <- oldNut
#i.e. - "P.tot.TP.pTP.mg.L" and "NH3.N.NH4.N.sol.mg.L"  as these looked to be low
dt.test <- data.frame(dQ = models2$meanlogQ_t2[models2$solute %in% sol]-models2$meanlogQ_t1[models2$solute %in% sol], dC = models2$meanlogC_t2[models2$solute %in% sol]-models2$meanlogC_t1[models2$solute %in% sol], mod = models2$mod.all[models2$solute %in% sol], solute = models2$solute[models2$solute %in% sol])
dt.test$mod[dt.test$mod %in% c("A","AA1","AA2")] <- "A"
dt.test$mod[dt.test$mod %in% c("D","DD1","DD2")] <- "D"
dt.test$mod <- as.character(dt.test$mod)
dt.test$mod <- as.factor(dt.test$mod)

res <- lm(dC~dQ+solute,data = dt.test)
plot(res)
library(mgcv)
res <- gam(dC~s(dQ, by= solute),data = dt.test)
plot(res)
summary(res)

b1.x <-
  brm(data = dt.test, family = gaussian,
      dC ~ 1 + dQ+solute,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 7)
f <-
  fitted(b1.x, newdata = dt.test) %>%  # we can use the same `nd` data from last time
  as_tibble() %>%
  bind_cols(dt.test)

plot(f$dQ,f$dC, xlab = "",ylab = "",col = gridcols[as.numeric(as.factor(models2$solute[models2$solute %in% sol]))],pch=19)
abline(h=0,lty=2)
abline(v=0,lty=2)
mtext(side=2,line=2.5,las=1,text = expression(paste(Delta,"C",sep="")))
mtext(side=1,line=2.5,las=1,text = expression(paste(Delta,"Q",sep="")))
pos.dq <- sort(f$dQ,index.return = TRUE)$ix
lines(f$dQ[pos.dq],f$Estimate[pos.dq])
polygon(c(f$dQ[pos.dq],rev(f$dQ[pos.dq]),f$dQ[pos.dq][1]),c(f$Q2.5[pos.dq],rev(f$Q97.5[pos.dq]),f$Q2.5[pos.dq][1]),col=add.alpha(gridcols[4],0.5))
legend("topleft",pch=19,col = gridcols[1:length(levels(as.factor(models2$solute[models2$solute %in% sol])))],legend = levels(as.factor(models2$solute[models2$solute %in% sol])),bty="n")



sol <- oldNut[-c(2,8,10,11)]
#i.e. - "P.tot.TP.pTP.mg.L" and "NH3.N.NH4.N.sol.mg.L"  as these looked to be low, test below shows these were significant negative intercept effect
dt.test <- data.frame(dQ = models2$meanlogQ_t2[models2$solute %in% sol]-models2$meanlogQ_t1[models2$solute %in% sol], dC = models2$meanlogC_t2[models2$solute %in% sol]-models2$meanlogC_t1[models2$solute %in% sol], mod = models2$mod.all[models2$solute %in% sol])
dt.test$mod[dt.test$mod %in% c("A","AA1","AA2")] <- "A"
dt.test$mod[dt.test$mod %in% c("D","DD1","DD2")] <- "D"
dt.test$mod <- as.character(dt.test$mod)
dt.test$mod <- as.factor(dt.test$mod)
res <- lm(dC~dQ,data = dt.test)
summary(res)
library(mgcv)

res <- gam(dC~s(dQ),data = dt.test)
summary(res)


b1.x <-
  brm(data = dt.test, family = gaussian,
      dC ~ 1 + dQ,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 7)
f <-
  fitted(b1.x, newdata = dt.test) %>%  # we can use the same `nd` data from last time
  as_tibble() %>%
  bind_cols(dt.test)
print(b1.x)
f.nut <- f




sol <- "Tracers"
dt.test <- data.frame(dQ = models2$meanlogQ_t2[models2$solClass %in% sol]-models2$meanlogQ_t1[models2$solClass %in% sol], dC = models2$meanlogC_t2[models2$solClass %in% sol]-models2$meanlogC_t1[models2$solClass %in% sol], mod = models2$mod.all[models2$solClass %in% sol])
dt.test$mod[dt.test$mod %in% c("A","AA1","AA2")] <- "A"
dt.test$mod[dt.test$mod %in% c("D","DD1","DD2")] <- "D"
dt.test$mod <- as.character(dt.test$mod)
dt.test$mod <- as.factor(dt.test$mod)
#dt.test <- dt.test[dt.test$mod == "D",]

res <- lm(dC~dQ,data = dt.test)
summary(res)
library(mgcv)
res <- gam(dC~s(dQ),data = dt.test)
summary(res)


b1.x <-
  brm(data = dt.test, family = gaussian,
      dC ~ 1 + dQ,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 7)
f <-
  fitted(b1.x, newdata = dt.test) %>%  # we can use the same `nd` data from last time
  as_tibble() %>%
  bind_cols(dt.test)
f.trace <- f
print(b1.x)

sol <- "Solids"
dt.test <- data.frame(dQ = models2$meanlogQ_t2[models2$solClass %in% sol]-models2$meanlogQ_t1[models2$solClass %in% sol], dC = models2$meanlogC_t2[models2$solClass %in% sol]-models2$meanlogC_t1[models2$solClass %in% sol], mod = models2$mod.all[models2$solClass %in% sol])
dt.test$mod[dt.test$mod %in% c("A","AA1","AA2")] <- "A"
dt.test$mod[dt.test$mod %in% c("D","DD1","DD2")] <- "D"
dt.test$mod <- as.character(dt.test$mod)
dt.test$mod <- as.factor(dt.test$mod)
#dt.test <- dt.test[dt.test$mod == "A",]

res <- lm(dC~dQ,data = dt.test)
summary(res)
library(mgcv)
res <- gam(dC~s(dQ),data = dt.test)
summary(res)


b1.x <-
  brm(data = dt.test, family = gaussian,
      dC ~ 1 + dQ,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 10), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 7)
f <- fitted(b1.x, newdata = dt.test) %>%  # we can use the same `nd` data from last time
  as_tibble() %>%
  bind_cols(dt.test)
print(b1.x)
f.solid <- f


par(mfrow=c(1,3),mar=c(4,4,0.5,1),oma=c(0.5,1,0.5,0.5))
sol <- oldNut[-c(2,8,10,11)]
plot(f.nut$dQ,f.nut$dC, xlab = "",ylab = "",col = "black",bg= gridcols[3],pch=21, ylim=c(-0.6,0.6),xlim=c(-1.5,1.5))
abline(h=0,lty=2)
abline(v=0,lty=2)
mtext(side=2,line=2.5,text = expression(paste(Delta,mu,"log(C)",sep="")))
mtext(side=1,line=2.5,las=1,text = expression(paste(Delta,mu,"log(Q)",sep="")))
pos.dq <- sort(f.nut$dQ,index.return = TRUE)$ix
lines(f.nut$dQ[pos.dq],f.nut$Estimate[pos.dq])
polygon(c(f.nut$dQ[pos.dq],rev(f.nut$dQ[pos.dq]),f.nut$dQ[pos.dq][1]),c(f.nut$Q2.5[pos.dq],rev(f.nut$Q97.5[pos.dq]),f.nut$Q2.5[pos.dq][1]),col=add.alpha(gridcols[4],0.5))
legend("topleft",legend = "(a) Nutrients",bty="n")

sol <- "Tracers"
plot(f.trace$dQ,f.trace$dC, xlab = "",ylab = "",ylim=c(-0.6,0.6),xlim=c(-1.5,1.5),pch=21,col = "black",bg= gridcols[5])
abline(h=0,lty=2)
abline(v=0,lty=2)
mtext(side=2,line=2.5,text = expression(paste(Delta,mu,"log(C)",sep="")))
mtext(side=1,line=2.5,las=1,text = expression(paste(Delta,mu,"log(Q)",sep="")))
pos.dq <- sort(f.trace$dQ,index.return = TRUE)$ix
lines(f.trace$dQ[pos.dq],f.trace$Estimate[pos.dq])
polygon(c(f.trace$dQ[pos.dq],rev(f.trace$dQ[pos.dq]),f.trace$dQ[pos.dq][1]),c(f.trace$Q2.5[pos.dq],rev(f.trace$Q97.5[pos.dq]),f.trace$Q2.5[pos.dq][1]),col=add.alpha(gridcols[4],0.5))
legend("topleft",legend = "(b) Tracers",bty="n")

sol <- "Solids"
plot(f.solid$dQ,f.solid$dC, xlab = "",ylab = "",ylim=c(-0.6,0.6),xlim=c(-1.5,1.5),pch=21,col = "black",bg= gridcols[8])
abline(h=0,lty=2)
abline(v=0,lty=2)
mtext(side=2,line=2.5,text = expression(paste(Delta,mu,"log(C)",sep="")))
mtext(side=1,line=2.5,las=1,text = expression(paste(Delta,mu,"log(Q)",sep="")))
pos.dq <- sort(f.solid$dQ,index.return = TRUE)$ix
lines(f.solid$dQ[pos.dq],f.solid$Estimate[pos.dq])
polygon(c(f.solid$dQ[pos.dq],rev(f.solid$dQ[pos.dq]),f.solid$dQ[pos.dq][1]),c(f.solid$Q2.5[pos.dq],rev(f.solid$Q97.5[pos.dq]),f.solid$Q2.5[pos.dq][1]),col=add.alpha(gridcols[4],0.5))
legend("topleft",legend = "(c) Solids",bty="n")






ggplot(data = as.data.frame(f),aes(x = dC, color = mod)) +
  geom_smooth(
    aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
        fill = mod),
    stat = "identity", 
    alpha = 1/4, size = 1/2) +
  geom_point(aes(y = dC),
             size = 2/3) +
  scale_colour_pander() +
  scale_fill_pander() +
  scale_x_continuous("dQ", expand = c(0, 0)) +
  ylab("dC") +
  theme_pander() + 
  theme(text = element_text(family = "Times"),
        legend.position = "none") +
  facet_wrap(~mod)




plot(bf$meanlogQ_t2[bf$solute %in% sol]-bf$meanlogQ_t1[bf$solute %in% sol],bf$meanlogC_t2[bf$solute %in% sol]-bf$meanlogC_t1[bf$solute %in% sol])
lines(seq(-0.8,0.4,by=0.2),predict(res,newdata = data.frame(dQ = seq(-0.8,0.4,by=0.2))))
library(mgcv)
res <- gam(dC~s(dQ),data = dt.test)
res2 <- gam(dC~s(dQ),data = dt.test2)
sol2 <- c("Ca.sol.mg.L","Cl.sol.mg.L","K.tot.mg.L","Mg.sol.mg.L","SiO2.Si.sol.react.mg.L" ,"Hardness.tot.CaCO3.Ca.Mg.mg.L")
plot(bf$meanlogQ_t2[bf$solute %in% sol2]-bf$meanlogQ_t1[bf$solute %in% sol2],bf$meanlogC_t2[bf$solute %in% sol2]-bf$meanlogC_t1[bf$solute %in% sol2])
res2 <- lm(dC~dQ,data = dt.test2)
dt.test2 <- data.frame(dQ = bf$meanlogQ_t2[bf$solute %in% sol2]-bf$meanlogQ_t1[bf$solute %in% sol2], dC = bf$meanlogC_t2[bf$solute %in% sol2]-bf$meanlogC_t1[bf$solute %in% sol2], model = models$mod.all)



lmbudyko <- lm(dEonE~dQonQ+dPETonPET,data = models)
summary(lmbudyko)
pattern <- c("A","D","C") #c("CA","AC","CD","DC","AA1","AA2","DD1","DD2","AD","DA")
pos <- models$mod.t1 =="A" & models$mod.t2 == "A" |
  models$mod.t1 =="C" & models$mod.t2 == "C" | 
  models$mod.t1 =="D" & models$mod.t2 == "D"#&  models$solClass == "Nutrients" 
plot(models$dQonQ[pos],
     abs(models$b1.t2[pos] - models$b1.t1[pos]))

plot(models$b1.t2[models$mod.t1 == pattern] - models$b1.t1[models$mod.t1 == pattern], models$b2.t2[models$mod.t1 == pattern] - models$b2.t1[models$mod.t1 == pattern])

plot(models$b1.t2[models$mod.t1 == pattern] - models$b1.t1[models$mod.t1 == pattern], models$dEonE[models$mod.t1 == pattern])
plot(models$b1.t2[models$mod.t1 == pattern] - models$b1.t1[models$mod.t1 == pattern], models$dPonP[models$mod.t1 == pattern])
plot(models$b1.t2[models$mod.t1 == pattern] - models$b1.t1[models$mod.t1 == pattern], models$dPETonPET[models$mod.t1 == pattern])


par(mfrow=c(3,1),oma=c(0,0,0,0),mar=c(4,4,0.5,0.5))
n <- 2*365
plot(sort(unique(bf$meanlogQ_t2[models$t1.days>n])), seq(1/length(unique(bf$meanlogQ_t2[models$t1.days>n])),1, length.out = length(unique(bf$meanlogQ_t2[models$t1.days>n]))), type = "l", lty=2, xlab = "mean logQ",ylab = "CDF",las=1, xlim = c(-2.5,1.5))
lines(sort(unique(bf$meanlogQ_t1[models$t1.days>n])), seq(1/length(unique(bf$meanlogQ_t1[models$t1.days>n])),1,length.out = length(unique(bf$meanlogQ_t1[models$t1.days>n]))),lty=1)
ks.test(unique(bf$meanlogQ_t2[models$t1.days>n]), unique(bf$meanlogQ_t1[models$t1.days>n]), alternative = "greater")

n <- 7*365
plot(sort(unique(bf$meanlogQ_t2[models$t1.days>n])), seq(1/length(unique(bf$meanlogQ_t2[models$t1.days>n])),1,length.out = length(unique(bf$meanlogQ_t2[models$t1.days>n]))), type = "l", lty=2, xlab = "mean logQ",ylab = "CDF",las=1,xlim = c(-2.5,1.5))
lines(sort(unique(bf$meanlogQ_t1[models$t1.days>n])), seq(1/length(unique(bf$meanlogQ_t1[models$t1.days>n])),1,length.out = length(unique(bf$meanlogQ_t1[models$t1.days>n]))),lty=1)
ks.test(unique(bf$meanlogQ_t2[models$t1.days>n]), unique(bf$meanlogQ_t1[models$t1.days>n]), alternative = "greater")

n <- 14*365
plot(sort(unique(bf$meanlogQ_t2[models$t1.days>n])), seq(1/length(unique(bf$meanlogQ_t2[models$t1.days>n])),1,length.out = length(unique(bf$meanlogQ_t2[models$t1.days>n]))), type = "l", lty=2, xlab = "mean logQ",ylab = "CDF",las=1,xlim = c(-2.5,1.5))
lines(sort(unique(bf$meanlogQ_t1[models$t1.days>n])), seq(1/length(unique(bf$meanlogQ_t1[models$t1.days>n])),1,length.out = length(unique(bf$meanlogQ_t1[models$t1.days>n]))),lty=1)
ks.test(unique(bf$meanlogQ_t2[models$t1.days>n]),unique(bf$meanlogQ_t1[models$t1.days>n]), alternative = "greater")

par(mfrow=c(3,1),oma=c(0,0,0,0),mar=c(4,4,0.5,0.5))
n <- 1*365
plot(density(unique(bf$meanlogQ_t2[models$t1.days>n])), type = "l", lty=2, xlab = "mean logQ",ylab = "CDF",las=1, xlim = c(-2.5,1.5), main = "")
lines(density(unique(bf$meanlogQ_t1[models$t1.days>n])),lty=1)
ks.test(unique(bf$meanlogQ_t2[models$t1.days>n]), unique(bf$meanlogQ_t1[models$t1.days>n]), alternative = "greater")

n <- 7*365
plot(density(unique(bf$meanlogQ_t2[models$t1.days>n])), type = "l", lty=2, xlab = "mean logQ",ylab = "CDF",las=1,xlim = c(-2.5,1.5), main = "")
lines(density(unique(bf$meanlogQ_t1[models$t1.days>n])),lty=1)
ks.test(unique(bf$meanlogQ_t2[models$t1.days>n]), unique(bf$meanlogQ_t1[models$t1.days>n]), alternative = "greater")

n <- 14*365
plot(density(sort(unique(bf$meanlogQ_t2[models$t1.days>n]))), type = "l", lty=2, xlab = "mean logQ",ylab = "CDF",las=1,xlim = c(-2.5,1.5), main = "")
lines(density(unique(bf$meanlogQ_t1[models$t1.days>n])),lty=1)
ks.test(unique(bf$meanlogQ_t2[models$t1.days>n]),unique(bf$meanlogQ_t1[models$t1.days>n]), alternative = "greater")

pvals <- vector("numeric",length = 20)
for (i in 1:20){
n <- i*365
res <- ks.test(unique(bf$meanlogQ_t2[models$t1.days>n]),unique(bf$meanlogQ_t1[models$t1.days>n]), alternative = "greater")
pvals[i] <- res$p.value
}
dev.off()
plot(pvals)
abline(h=0.01)
abline(h=0.05,lty=2)



plot(models$b2.t2[models$mod.t1 == "CA" & models$mod.t2 == "CA"],models$b2.t1[models$mod.t1 == "CA" & models$mod.t2 == "CA"]-models$b2.t2[models$mod.t1 == "CA" & models$mod.t2 == "CA"], xlim=c(-1.2,1.2),ylim=c(-0.7,0.7))
mod <- "DA"
points(models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],models$b2.t1[models$mod.t1 == mod & models$mod.t2 == mod]-models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],pch=2,col="blue")
mod <- "CD"
points(models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],models$b2.t1[models$mod.t1 == mod & models$mod.t2 == mod]-models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],pch=2,col="blue")
mod <- "AD"
points(models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],models$b2.t1[models$mod.t1 == mod & models$mod.t2 == mod]-models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],pch=2,col="blue")
mod <- "AA1"
points(models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],models$b2.t1[models$mod.t1 == mod & models$mod.t2 == mod]-models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],pch=2,col="blue")
mod <- "AA2"
points(models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],models$b2.t1[models$mod.t1 == mod & models$mod.t2 == mod]-models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],pch=2,col="blue")
mod <- "DD1"
points(models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],models$b2.t1[models$mod.t1 == mod & models$mod.t2 == mod]-models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],pch=2,col="blue")
mod <- "DD2"
points(models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],models$b2.t1[models$mod.t1 == mod & models$mod.t2 == mod]-models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],pch=2,col="blue")
mod <- "AC"
points(models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],models$b2.t1[models$mod.t1 == mod & models$mod.t2 == mod]-models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],pch=2,col="blue")
mod <- "DC"
points(models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],models$b2.t1[models$mod.t1 == mod & models$mod.t2 == mod]-models$b2.t2[models$mod.t1 == mod & models$mod.t2 == mod],pch=2,col="blue")

par(mfrow=c(4,3), mar=c(1,1,1,1),oma=c(1,1,1,1))
mod <- "A"
plot(models$b2.t1[models$mod.t1 == mod],models$b2.t2[models$mod.t1 == mod]-models$b2.t1[models$mod.t1 == mod],pch=2,col="blue",xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "C"
plot(models$b2.t1[models$mod.t1 == mod],models$b2.t2[models$mod.t1 == mod]-models$b2.t1[models$mod.t1 == mod],pch=2,col="orange",xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "D"
plot(models$b2.t1[models$mod.t1 == mod],models$b2.t2[models$mod.t1 == mod]-models$b2.t1[models$mod.t1 == mod],pch=2,col="red",xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "AA1"
plot(models$b2.t1[models$mod.t1 == mod],models$b2.t2[models$mod.t1 == mod]-models$b2.t1[models$mod.t1 == mod],pch=19,col="blue",xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "AA2"
plot(models$b2.t1[models$mod.t1 == mod],models$b2.t2[models$mod.t1 == mod]-models$b2.t1[models$mod.t1 == mod],pch=20,col="blue",xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "DD1"
plot(models$b2.t1[models$mod.t1 == mod],models$b2.t2[models$mod.t1 == mod]-models$b2.t1[models$mod.t1 == mod],pch=19,col="red",xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "DD2"
plot(models$b2.t1[models$mod.t1 == mod],models$b2.t2[models$mod.t1 == mod]-models$b2.t1[models$mod.t1 == mod],pch=20,col="red",xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "CA"
plot(models$b2.t1[models$mod.t1 == mod],models$b2.t2[models$mod.t1 == mod]-models$b2.t1[models$mod.t1 == mod],pch=20,col="green",xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "AD"
plot(models$b2.t1[models$mod.t1 == mod],models$b2.t2[models$mod.t1 == mod]-models$b2.t1[models$mod.t1 == mod],pch=19,col="pink",xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "DA"
plot(models$b2.t1[models$mod.t1 == mod],models$b2.t2[models$mod.t1 == mod]-models$b2.t1[models$mod.t1 == mod],pch=1,col="pink",xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "CD"
plot(models$b2.t1[models$mod.t1 == mod],models$b2.t2[models$mod.t1 == mod]-models$b2.t1[models$mod.t1 == mod],pch=1,col="red",xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "DC"
plot(models$b2.t1[models$mod.t1 == mod],models$b2.t2[models$mod.t1 == mod]-models$b2.t1[models$mod.t1 == mod],pch=1,col="red",xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
#mod <- "AC"
#plot(models$b2.t2[models$mod.t1 == mod],models$b2.t2[models$mod.t1 == mod]-models$b2.t1[models$mod.t1 == mod],pch=19,col="green",xlim=c(-1.2,1.2),ylim=c(-0.7,0.7))
#abline(0,-1); abline(v=0); abline(h=0);

par(mfrow=c(4,3), mar=c(1,1,1,1),oma=c(1,1,1,1))
mod <- "A"
plot(models$b1.t1[models$mod.t1 == mod],models$b1.t2[models$mod.t1 == mod]-models$b1.t1[models$mod.t1 == mod],xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "C"
plot(models$b1.t1[models$mod.t1 == mod],models$b1.t2[models$mod.t1 == mod]-models$b1.t1[models$mod.t1 == mod],xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "D"
plot(models$b1.t1[models$mod.t1 == mod],models$b1.t2[models$mod.t1 == mod]-models$b1.t1[models$mod.t1 == mod],xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "AA1"
plot(models$b1.t1[models$mod.t1 == mod],models$b1.t2[models$mod.t1 == mod]-models$b1.t1[models$mod.t1 == mod],xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "AA2"
plot(models$b1.t1[models$mod.t1 == mod],models$b1.t2[models$mod.t1 == mod]-models$b1.t1[models$mod.t1 == mod],xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "DD1"
plot(models$b1.t1[models$mod.t1 == mod],models$b1.t2[models$mod.t1 == mod]-models$b1.t1[models$mod.t1 == mod],xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "DD2"
plot(models$b1.t1[models$mod.t1 == mod],models$b1.t2[models$mod.t1 == mod]-models$b1.t1[models$mod.t1 == mod],xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "CA"
plot(models$b1.t1[models$mod.t1 == mod],models$b1.t2[models$mod.t1 == mod]-models$b1.t1[models$mod.t1 == mod],xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "AD"
plot(models$b1.t1[models$mod.t1 == mod],models$b1.t2[models$mod.t1 == mod]-models$b1.t1[models$mod.t1 == mod],xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "DA"
plot(models$b1.t1[models$mod.t1 == mod],models$b1.t2[models$mod.t1 == mod]-models$b1.t1[models$mod.t1 == mod],xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "CD"
plot(models$b1.t1[models$mod.t1 == mod],models$b1.t2[models$mod.t1 == mod]-models$b1.t1[models$mod.t1 == mod],xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);
mod <- "DC"
plot(models$b1.t1[models$mod.t1 == mod],models$b1.t2[models$mod.t1 == mod]-models$b1.t1[models$mod.t1 == mod],xlim=c(-1.2,1.2),ylim=c(-0.7,0.7)); text(-1,-0.6,mod)
abline(0,-1); abline(v=0); abline(h=0);

plot(models$dEonE,predict(lmbudyko,newdata = models))
  
boxplot(b1.t1~mod.t1,data=models, ylim=c(-1.5,1.5))
boxplot(b2.t1~mod.t1,data=models,add=TRUE,col="blue")

boxplot(b1.t2~mod.t2,data=models, ylim=c(-1.5,1.5))
boxplot(b2.t2~mod.t2,data=models,add=TRUE,col="blue")

boxplot(b1.tall~mod.all,data=models, ylim=c(-1.5,1.5))
boxplot(b2.tall~mod.all,data=models,add=TRUE,col="blue")

pos.same <- which((modt1  == "A" & modt2 == "A" | 
                     modt1  == "C" & modt2 == "C" | 
                     modt1  == "D" & modt2 == "D")  & t1days > 365*4 & t2days > 365*4)
pars.t1 <- t(sapply(pos.same, FUN = function(x) t1res[[x]]$pars))
pars.t2 <- t(sapply(pos.same, FUN = function(x) t2res[[x]]$pars))
solutes <- bf$solute[pos.same]
db <- data.frame(mod = models$mod.t1[pos.same], 
                 solute = solutes,
                 t1.days = models$t1.days[pos.same],
                 b1 = pars.t1[,2],b2 = pars.t2[,2], deltab = pars.t2[,2] - pars.t1[,2], stringsAsFactors = FALSE)








library(plotrix)

plotRoseMod <- function(models,mod = "A", modtx = "mod.t1"){
  #mod <- "A"
  pos<- models[[modtx]] == mod & models$t1.days > 365*1 & models$t2.days > 365*1
  mydata <- data.frame(x = models$b1.t1[pos],
                       y = models$b1.t2[pos]- models$b1.t1[pos], 
                       mt1 = models$mod.t1[pos],
                       mt2 = models$mod.t2[pos])
  
  #calculate angles and speed from xy data
  angles <- atan2(mydata$y,mydata$x)*180/pi #degrees
  speed <- sqrt(mydata$x^2+mydata$y^2)
  mydata$ws <- speed
  mydata$wd <- angles
  mydata$wd[mydata$wd<0] <- 360 + mydata$wd[mydata$wd<0]
  angle <- 45/2
  angles <- seq(0,360,by=angle)
  mydata$sector <-  angle * (as.numeric(cut(mydata$wd, seq(0, 360, angle)))-0.5)
  
  mid.angles <- seq(angle/2,360,by=angle)
  
  #loop through each mid.angles count number of models
  mods <- levels(models$mod.t1)
  piedata <- vector("list",length(mid.angles))
  for (i in 1:length(piedata)){
    piedata[[i]] <- vector("numeric",length = length(mods))
  }
  for (i in 1:length(mid.angles)){
    vec <- piedata[[i]] 
    for (j in 1:length(mods)){
      vec[j] <- sum(as.character(mydata$mt2) == mods[j] & 
                      mydata$sector == mid.angles[i])
    }
    piedata[[i]] <- vec
  }
  for (i in 1:length(piedata)){
    piedata[[i]] <- c(0,cumsum(piedata[[i]])/length(mydata$sector)*100)
  }
  
  
  radial.extents <- piedata
  #radial.pie <-
  #function (radial.extents){
  sector.edges = NULL; 
  sector.colors = NULL; 
  cs1 = c(0, 1); 
  cs2 = c(0, 1); 
  cs3 = c(0, 1); 
  alpha = 1; 
  labels = NA; 
  label.pos = NULL;
  radlab = FALSE;
  start = 0; 
  clockwise = FALSE; 
  label.prop = 1.1;
  radial.lim = NULL;
  main = "";
  xlab = ""; 
  ylab = ""; 
  #mar = c(2, 2, 3, 2); 
  show.grid = TRUE; 
  show.grid.labels = 4; 
  show.radial.grid = TRUE; 
  rad.col = "black";
  grid.col = "gray";
  grid.bg = "transparent"; 
  grid.unit = NULL; 
  radial.labels = NULL;
  boxed.radial = TRUE; 
  add = FALSE;
  
  if (is.null(radial.lim)) radial.lim <- range(radial.extents)
  if (is.null(sector.edges)) {
    if (clockwise) {
      sector.edges <- seq(2 * pi + start, start, length.out = length(radial.extents) + 1)
    } else {
      sector.edges <- seq(start, 2 * pi + start, length.out = length(radial.extents) +  1)
    }
  }
  if (is.null(label.pos)) label.pos <- sector.edges[-length(sector.edges)] + diff(sector.edges)/2
  if (show.grid) {
    maxrad <- max(unlist(radial.extents))
    if (length(radial.lim) < 3) {
      grid.pos <- seq(0,maxrad,length.out = 6)
    } else {
      #grid.pos <- radial.lim
      grid.pos <- seq(0,maxrad,length.out = 6)
    }
    if (grid.pos[1] < radial.lim[1]) grid.pos <- grid.pos[-1]
    maxlength <- max(grid.pos - radial.lim[1])
  } else {
    grid.pos <- NA
    maxlength <- diff(radial.lim)
  }
  par(pty = "s")
  
  plot(0, xlim = c(-maxrad, maxrad), ylim = c(-maxrad, 
                                              maxrad), type = "n", axes = FALSE, xlab = xlab, ylab = ylab, asp=1)
  if (show.grid) {
    radial.grid(labels = rep("",length(c(0,pi/2,pi/2+pi/4,pi,3*pi/2,2*pi-pi/4))), label.pos = c(0,pi/2,pi/2+pi/4,pi,3*pi/2,2*pi-pi/4), 
                radlab = radlab, radial.lim = radial.lim, start = start, 
                rad.col = rad.col, 
                clockwise = clockwise, label.prop = label.prop, 
                grid.pos = grid.pos, grid.col = grid.col, grid.bg = grid.bg)
  }
  
  nsectors <- length(radial.extents)
  dtheta <- 2*pi/360
  set.seed(13)
  library(RColorBrewer)
  darkcols <- sample(brewer.pal(12, "Paired"), 12,replace = FALSE)
  if (is.list(radial.extents)) {
    
    for (sector in 1:nsectors) {
      annuli <- radial.extents[[sector]]
      
      for (annulus in 1:(length(annuli) - 1)) {
        
        drawSectorAnnulus(sector.edges[[sector]]+dtheta, 
                          sector.edges[[sector + 1]]-dtheta, 
                          annuli[annulus], 
                          annuli[annulus + 1], 
                          darkcols[annulus])
      }
    }
  } else {
    if (is.null(sector.colors)) 
      sector.colors <- rainbow(nsectors)
    for (sector in 1:nsectors) {
      drawSectorAnnulus(sector.edges[sector], sector.edges[sector + 
                                                             1], 0, radial.extents[sector], sector.colors[sector])
    }
  }
  
  if (show.grid) {
    lines(c(-maxrad, maxrad),c(0,0),col=rad.col,lwd=1)
    lines(c(0,0),c(-maxrad, maxrad),col=rad.col,lwd=1)
    lines(c(maxrad*cos(13*pi/8), maxrad*cos(5*pi/8)),c(maxrad*sin(13*pi/8), maxrad*sin(5*pi/8)),col=rad.col,lwd=1)
  }
}
models$b2.t1[models$mod.t1 == "A"] <- models$b1.t1[models$mod.t1 == "A"]
models$b2.t2[models$mod.t2 == "A"] <- models$b1.t2[models$mod.t2 == "A"]
models$b2.t1[models$mod.t1 == "C"] <- models$b1.t1[models$mod.t1 == "C"]
models$b2.t2[models$mod.t2 == "C"] <- models$b1.t2[models$mod.t2 == "C"]
models$b2.t1[models$mod.t1 == "D"] <- models$b1.t1[models$mod.t1 == "D"]
models$b2.t2[models$mod.t2 == "D"] <- models$b1.t2[models$mod.t2 == "D"]
plotRoseMod2 <- function(models, mod = "A", modtx = "mod.t1"){
  #mod <- "A"
  pos<- models[[modtx]] == mod & models$t1.days > 365*1 & models$t2.days > 365*1
  mydata <- data.frame(x = models$b2.t1[pos],
                       y = models$b2.t2[pos]- models$b2.t1[pos], 
                       mt1 = models$mod.t1[pos],
                       mt2 = models$mod.t2[pos])
  #calculate angles and speed from xy data
  angles <- atan2(mydata$y,mydata$x)*180/pi #degrees
  speed <- sqrt(mydata$x^2+mydata$y^2)
  mydata$ws <- speed
  mydata$wd <- angles
  mydata$wd[mydata$wd<0] <- 360 + mydata$wd[mydata$wd<0]
  angle <- 45/2
  angles <- seq(0,360,by=angle)
  mydata$sector <-  angle * (as.numeric(cut(mydata$wd, seq(0, 360, angle)))-0.5)
  
  mid.angles <- seq(angle/2,360,by=angle)
  
  #loop through each mid.angles count number of models
  mods <- levels(mydata$mt1)
  piedata <- vector("list",length(mid.angles))
  for (i in 1:length(piedata)){
    piedata[[i]] <- vector("numeric",length = length(mods))
  }
  for (i in 1:length(mid.angles)){
    vec <- piedata[[i]] 
    for (j in 1:length(mods)){
      vec[j] <- sum(as.character(mydata$mt2) == mods[j] & 
                      mydata$sector == mid.angles[i])
    }
    piedata[[i]] <- vec
  }
  for (i in 1:length(piedata)){
    piedata[[i]] <- c(0,cumsum(piedata[[i]])/length(mydata$sector)*100)
  }
  
  
  radial.extents <- piedata
  #radial.pie <-
  #function (radial.extents){
  sector.edges = NULL; 
  sector.colors = NULL; 
  cs1 = c(0, 1); 
  cs2 = c(0, 1); 
  cs3 = c(0, 1); 
  alpha = 1; 
  labels = NA; 
  label.pos = NULL;
  radlab = FALSE;
  start = 0; 
  clockwise = FALSE; 
  label.prop = 1.1;
  radial.lim = NULL;
  main = "";
  xlab = ""; 
  ylab = ""; 
  #mar = c(2, 2, 3, 2); 
  show.grid = TRUE; 
  show.grid.labels = 4; 
  show.radial.grid = TRUE; 
  rad.col = "black";
  grid.col = "gray";
  grid.bg = "transparent"; 
  grid.unit = NULL; 
  radial.labels = NULL;
  boxed.radial = TRUE; 
  add = FALSE;
  
  if (is.null(radial.lim)) radial.lim <- range(radial.extents)
  if (is.null(sector.edges)) {
    if (clockwise) {
      sector.edges <- seq(2 * pi + start, start, length.out = length(radial.extents) + 1)
    } else {
      sector.edges <- seq(start, 2 * pi + start, length.out = length(radial.extents) +  1)
    }
  }
  if (is.null(label.pos)) label.pos <- sector.edges[-length(sector.edges)] + diff(sector.edges)/2
  if (show.grid) {
    maxrad <- max(unlist(radial.extents))
    if (length(radial.lim) < 3) {
      grid.pos <- seq(0,maxrad,length.out = 6)
    } else {
      #grid.pos <- radial.lim
      grid.pos <- seq(0,maxrad,length.out = 6)
    }
    if (grid.pos[1] < radial.lim[1]) grid.pos <- grid.pos[-1]
    maxlength <- max(grid.pos - radial.lim[1])
  } else {
    grid.pos <- NA
    maxlength <- diff(radial.lim)
  }
  par(pty = "s")
  
  plot(0, xlim = c(-maxrad, maxrad), ylim = c(-maxrad, 
                                              maxrad), type = "n", axes = FALSE, xlab = xlab, ylab = ylab, asp=1)
  if (show.grid) {
    radial.grid(labels = rep("",length(c(0,pi/2,pi/2+pi/4,pi,3*pi/2,2*pi-pi/4))), label.pos = c(0,pi/2,pi/2+pi/4,pi,3*pi/2,2*pi-pi/4), 
                radlab = radlab, radial.lim = radial.lim, start = start, 
                rad.col = rad.col, 
                clockwise = clockwise, label.prop = label.prop, 
                grid.pos = grid.pos, grid.col = grid.col, grid.bg = grid.bg)
  }
  
  nsectors <- length(radial.extents)
  dtheta <- 2*pi/360
  set.seed(13)
  library(RColorBrewer)
  darkcols <- sample(brewer.pal(12, "Paired"), 12,replace = FALSE)
  if (is.list(radial.extents)) {
    
    for (sector in 1:nsectors) {
      annuli <- radial.extents[[sector]]
      
      for (annulus in 1:(length(annuli) - 1)) {
        
        drawSectorAnnulus(sector.edges[[sector]]+dtheta, 
                          sector.edges[[sector + 1]]-dtheta, 
                          annuli[annulus], 
                          annuli[annulus + 1], 
                          darkcols[annulus])
      }
    }
  } else {
    if (is.null(sector.colors)) 
      sector.colors <- rainbow(nsectors)
    for (sector in 1:nsectors) {
      drawSectorAnnulus(sector.edges[sector], sector.edges[sector + 
                                                             1], 0, radial.extents[sector], sector.colors[sector])
    }
  }
  
  if (show.grid) {
    lines(c(-maxrad, maxrad),c(0,0),col=rad.col,lwd=1)
    lines(c(0,0),c(-maxrad, maxrad),col=rad.col,lwd=1)
    lines(c(maxrad*cos(13*pi/8), maxrad*cos(5*pi/8)),c(maxrad*sin(13*pi/8), maxrad*sin(5*pi/8)),col=rad.col,lwd=1)
  }
}

pdf(file = "roseplot.pdf",width = 7, height = 4)
par(mfrow=c(4,7), oma=c(0.5,0.5,2,1),mar=c(0.5,0.5,0.5,0.5))
plotRoseMod(models,mod = "A")
mtext(side=3,line=-0.5,adj=0,"A",cex=0.8)
plotRoseMod(models,mod = "C")
mtext(side=3,line=-0.5,adj=0,"C",cex=0.8)
mtext(side=3,line=0.5,expression(paste(b[1])),cex=1.5)
plotRoseMod(models,mod = "D")
mtext(side=3,line=-0.5,adj=0,"D",cex=0.8)

plotRoseMod2(models,mod = "A")
mtext(side=3,line=-0.5,adj=0,"A",cex=0.8)
plotRoseMod2(models,mod = "C")
mtext(side=3,line=-0.5,adj=0,"C",cex=0.8)
mtext(side=3,line=0.5,expression(paste(b[2])),cex=1.5)
plotRoseMod2(models,mod = "D")
mtext(side=3,line=-0.5,adj=0,"D",cex=0.8)

plot(1,1,typ="n",axes=FALSE,xlab="",ylab="", xlim=c(0,12),ylim=c(0,12), xaxs="i",yaxs="i")
set.seed(13)
library(RColorBrewer)
darkcols <- sample(brewer.pal(12, "Paired"), 12,replace = FALSE)
points(c(1,1,1,1, 5,5,5, 9,9,9,9),
       c(2,5,8,11,5,8,11,2,5,8,11), pch = 22,col = "black", 
       bg = darkcols[c(4,3,2,1,7,6,5,9,11,10,8)],cex=2)
mtext(side=3,line = 0, expression(paste("C-Q ",t[2])))
text(1,2, paste(levels(models$mod.t1)[4]),cex=0.6, pos=4)
text(1,5, paste(levels(models$mod.t1)[3]),cex=0.6,pos=4)
text(1,8, paste(levels(models$mod.t1)[2]),cex=0.6,pos=4)
text(1,11,paste(levels(models$mod.t1)[1]),cex=0.6, pos=4)
text(5,5, paste(levels(models$mod.t1)[7]),cex=0.6,pos=4) #CD
text(5,8, paste(levels(models$mod.t1)[6]),cex=0.6,pos=4)
text(5,11,paste(levels(models$mod.t1)[5]),cex=0.6, pos=4)
text(9,2, paste(levels(models$mod.t1)[9]),cex=0.6,pos=4) #DA
text(9,5, paste(levels(models$mod.t1)[11]),cex=0.6,pos=4)
text(9,8, paste(levels(models$mod.t1)[10]),cex=0.6, pos=4)
text(9,11,paste(levels(models$mod.t1)[8]),cex=0.6,pos=4)
#text(9,11,paste(levels(models$mod.t1)[12]),cex=0.9,pos=4)


plotRoseMod(models,mod = "AA1")
mtext(side=3,line=-0.5,adj=0,"AA1",cex=0.8)
plotRoseMod(models,mod = "CA")
mtext(side=3,line=-0.5,adj=0,"CA",cex=0.8)
plotRoseMod(models,mod = "DD1")
mtext(side=3,line=-0.5,adj=0,"DD1",cex=0.8)

plotRoseMod2(models,mod = "AA1")
mtext(side=3,line=-0.5,adj=0,"AA1",cex=0.8)
plotRoseMod2(models,mod = "CA")
mtext(side=3,line=-0.5,adj=0,"CA",cex=0.8)
plotRoseMod2(models,mod = "DD1")
mtext(side=3,line=-0.5,adj=0,"DD1",cex=0.8)

plot(0, xlim = c(-1, 1), ylim = c(-1, 1), type = "n", axes = FALSE, xlab = expression(paste(b[1],"(",t[1],")")), ylab = expression(paste(Delta,"b")), asp=1)
radial.grid(labels = rep("",length(c(0,pi/2,pi/2+pi/4,pi,3*pi/2,2*pi-pi/4))), label.pos = c(0,pi/2,pi/2+pi/4,pi,3*pi/2,2*pi-pi/4), 
            radlab = FALSE, radial.lim = c(0,1), start = 0, 
            rad.col = "black", 
            clockwise = FALSE, label.prop = 1.1, 
            grid.pos = seq(0,1,length.out = 6), grid.col ="grey", grid.bg ="transparent")
lines(c(1*cos(13*pi/8), 1*cos(5*pi/8)),c(1*sin(13*pi/8), 1*sin(5*pi/8)),col="black",lwd=1)
mtext(side=1,line=0,expression(paste("-",Delta,"b")))
mtext(side=3,line=0,expression(paste("+",Delta,"b")))
mtext(side=2,line=0,paste("-b"),las=1)
mtext(side=4,line=0,paste("+b"),las=1)

plotRoseMod(models,mod = "AA2")
mtext(side=3,line=-0.5,adj=0,"AA2",cex=0.8)
plotRoseMod(models,mod = "CD")
mtext(side=3,line=-0.5,adj=0,"CD",cex=0.8)
plotRoseMod(models,mod = "DD2")
mtext(side=3,line=-0.5,adj=0,"DD2",cex=0.8)


plotRoseMod2(models,mod = "AA2")
mtext(side=3,line=-0.5,adj=0,"AA2",cex=0.8)
plotRoseMod2(models,mod = "CD")
mtext(side=3,line=-0.5,adj=0,"CD",cex=0.8)
plotRoseMod2(models,mod = "DD2")
mtext(side=3,line=-0.5,adj=0,"DD2",cex=0.8)
plot(0, xlim = c(-1, 1), ylim = c(-1, 1), type = "n", axes = FALSE,
     xlab = "", ylab = "")

plotRoseMod(models,mod = "AD")
mtext(side=3,line=-0.5,adj=0,"AD",cex=0.8)
plot(0, xlim = c(-1, 1), ylim = c(-1, 1), type = "n", axes = FALSE,
     xlab = "", ylab = "")
plotRoseMod(models,mod = "DA")
mtext(side=3,line=-0.5,adj=0,"DA",cex=0.8)


plotRoseMod2(models,mod = "AD")
mtext(side=3,line=-0.5,adj=0,"AD",cex=0.8)
plot(0, xlim = c(-1, 1), ylim = c(-1, 1), type = "n", axes = FALSE,
     xlab = "", ylab = "")
plotRoseMod2(models,mod = "DA")
mtext(side=3,line=-0.5,adj=0,"DA",cex=0.8)
plot(0, xlim = c(-1, 1), ylim = c(-1, 1), type = "n", axes = FALSE,
     xlab = "", ylab = "")
dev.off()
# > summary(lmsolids)
# 
# Call:
#   lm(formula = db$deltab[db$solute %in% oldSols] ~ db$b1[db$solute %in% 
#                                                            oldSols])
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.122879 -0.017079  0.002921  0.022940  0.109874 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                   0.013851   0.008396   1.650   0.1048  
# db$b1[db$solute %in% oldSols] 0.150884   0.059826   2.522   0.0146 *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.04067 on 54 degrees of freedom
# Multiple R-squared:  0.1054,	Adjusted R-squared:  0.08881 
# F-statistic: 6.361 on 1 and 54 DF,  p-value: 0.01465
# 
# > summary(lmnutrients)
# 
# Call:
#   lm(formula = db$deltab[db$solute %in% oldNut] ~ db$b1[db$solute %in% 
#                                                           oldNut])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.29044 -0.02113  0.01698  0.03529  0.15392 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                  -0.02551    0.00827  -3.085  0.00265 **
#   db$b1[db$solute %in% oldNut]  0.09401    0.03240   2.902  0.00459 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.06527 on 98 degrees of freedom
# Multiple R-squared:  0.07911,	Adjusted R-squared:  0.06971 
# F-statistic: 8.419 on 1 and 98 DF,  p-value: 0.004586
# 
# > summary(lmtracers)
# 
# Call:
#   lm(formula = db$deltab[db$solute %in% oldSolids] ~ db$b1[db$solute %in% 
#                                                              oldSolids])
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.163281 -0.025220 -0.004653  0.032191  0.110037 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                     -0.01144    0.01113  -1.027   0.3109  
# db$b1[db$solute %in% oldSolids]  0.11944    0.04821   2.477   0.0179 *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.05367 on 37 degrees of freedom
# Multiple R-squared:  0.1423,	Adjusted R-squared:  0.1191 
# F-statistic: 6.137 on 1 and 37 DF,  p-value: 0.01793


#Get data for piecewise regressions
#Extract second exponent and plot how this changes as a function of model type

pos.same <- which(models$mod.t1 %in% c("AA1","AA2","DD1","DD2","AC","AD","DC","DA","CA","CD") & models$mod.t2 %in% c("AA1","AA2","DD1","DD2","AC","AD","DC","DA","CA","CD") & (as.Date(bf$end.date) - as.Date(bf$start.date))/365 > 20)
pars.t1 <- t(sapply(pos.same, FUN = function(x) t1res[[x]]$pars))
pars.t2 <- t(sapply(pos.same, FUN = function(x) t2res[[x]]$pars))
solutes <- bf$solute[pos.same]
#c(df$pw_a500, df$pw_x500, df$pw_delta500, df$pw_b500)
db <- data.frame(mod = models$mod.t1[pos.same], 
                 solute = solutes,
                 a1 = pars.t1[,1],
                 x1 = pars.t1[,2], 
                 delta1 = pars.t1[,3], 
                 b11 = pars.t1[,4], 
                 b12 =  pars.t1[,4] - pars.t1[,3],
                 a2 = pars.t2[,1], 
                 x2 = pars.t2[,2], 
                 delta2 = pars.t2[,3], 
                 b21 = pars.t2[,4], 
                 b22 =  pars.t2[,4] - pars.t2[,3],
                 stringsAsFactors = FALSE)

#c("AA1","AA2","DD1","DD2","AC","AD","DC","DA","CA","CD")
par(mfrow=c(2,2))
cols <- gridcols[sapply(paste(mods), FUN = function(x) { which(names(gridcols) == x)})]
  mods <- models$mod.t1[pos.same]

plot(db$a1,db$a2, col = cols, pch = 19, xlab = "",ylab = "", xlim=c(-4,4),ylim=c(-4,4))
mtext(side=1,line=2.5,"a(t1)")
mtext(side=2,line=2.5,"a(t2)")
lines(c(-2,4),c(-2,4))

legend("topleft",pch = 19, 
       col = gridcols[which(gridcols %in% unique(cols))], 
       legend = names(gridcols)[which(gridcols %in% unique(cols))], bty = "n")
plot(db$b11,db$b21,col = cols, pch = 19, xlab = "",ylab = "", ylim = c(-1,1),xlim=c(-1,1))
mtext(side=1,line=2.5,"b(t1)")
mtext(side=2,line=2.5,"b(t2)")
abline(0,1)
plot(db$x1,db$x2,col = cols, pch = 19, xlab = "",ylab = "", ylim=c(-2.5,2.5),xlim=c(-2.5,2.5))
mtext(side=1,line=2.5,"X0(t1)")
mtext(side=2,line=2.5,"X0(t2)")
abline(0,1)
plot(db$delta1,db$delta2,col = cols, pch = 19, xlab = "",ylab = "", xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
mtext(side=1,line=2.5,expression(paste(delta,"(t1)",sep="")))
mtext(side=2,line=2.5,expression(paste(delta,"(t2)",sep="")))
abline(0,1)
#Check how delta used to classify models



#Spatial plots
library(rgdal)
library(sf)
sites <- list.files(path = "D:\\WIRT\\GIS",
                    pattern = "*.shp")
sites <- sites[grep("Basin",sites)]
SITE_REFs <- unlist(strsplit(sites,"Basin"))
SITE_REFs <- SITE_REFs[SITE_REFs != ""]
SITE_REFs <- unlist(strsplit(SITE_REFs,".shp"))
setwd("D:\\WIRT\\GIS")
catchments <- 
 lapply(sites, function(x) st_read("D:\\WIRT\\GIS", layer = gsub(".shp","",x)))
pos.read <- which(unlist(lapply(catchments, FUN = function(x) length(names(x)))) == 7)
#sort by area, largest to smallest

plot.order <- sort(unlist(lapply(catchments, FUN = st_area)), decreasing = TRUE,index.return = TRUE)$ix
plot.order <- plot.order[which(plot.order %in% pos.read)]
catchments <- do.call(rbind, lapply(sites[plot.order], function(x) st_read("D:\\WIRT\\GIS", layer = gsub(".shp","",x))))
names(catchments)
catchments$SITE_REF

catchments["dC_Cl"] <- NA
for (i in 1:length(catchments$SITE_REF)){
  sr <- as.character(catchments$SITE_REF[i])
  print(sr)
  pos <- which(models2$site.id == sr & models2$solute == "Cl.sol.mg.L")
  print(pos)
  if (length(pos)>0){
  catchments$dC_Cl[i] <- models2$meanlogC_t2[pos] - models2$meanlogC_t1[pos] 
  }
}


catchments["dC_Q"] <- NA
for (i in 1:length(catchments$SITE_REF)){
  sr <- as.character(catchments$SITE_REF[i])
  print(sr)
  pos <- which(models2$site.id == sr)
  print(pos)
  if (length(pos)>0){
    catchments$dC_Q[i] <- models2$meanlogQ_t2[pos[1]] - models2$meanlogQ_t1[pos[1]] 
  }
}
plot(catchments["dC_Q"])

catchments["PET"] <- NA
for (i in 1:length(catchments$SITE_REF)){
  sr <- as.character(catchments$SITE_REF[i])
  print(sr)
  pos <- which(models2$site.id == sr)
  print(pos)
  if (length(pos)>0){
    catchments$PET[i] <- models2$PET.all[pos[1]] 
  }
}
plot(catchments["PET"])

catchments["P"] <- NA
for (i in 1:length(catchments$SITE_REF)){
  sr <- as.character(catchments$SITE_REF[i])
  print(sr)
  pos <- which(models2$site.id == sr)
  print(pos)
  if (length(pos)>0){
    catchments$P[i] <- models2$P.all[pos[1]] 
  }
}
plot(catchments["P"])

Aust <- st_read("F:\\GIS\\WorldHighRes\\Australia.shp")

catchments["Aridity"] <- NA
for (i in 1:length(catchments$SITE_REF)){
  sr <- as.character(catchments$SITE_REF[i])
  print(sr)
  pos <- which(models2$site.id == sr)
  print(pos)
  if (length(pos)>0){
    catchments$Aridity[i] <- models2$Aridity.all[pos[1]] 
  }
}
plot(catchments["Aridity"])

catchments["db_nut"] <- NA
for (i in 1:length(catchments$SITE_REF)){
  sr <- as.character(catchments$SITE_REF[i])
  print(sr)
  pos <- which(models2$site.id == sr & models2$solClass == "Nutrients")
  print(pos)
  if (length(pos)>0){
    catchments$db_nut[i] <- mean(models2$b1.t2[pos] -  models2$b1.t1[pos],na.rm=TRUE)
  }
}
plot(catchments["db_nut"])

catchments["db_tra"] <- NA
for (i in 1:length(catchments$SITE_REF)){
  sr <- as.character(catchments$SITE_REF[i])
  print(sr)
  pos <- which(models2$site.id == sr & models2$solClass == "Tracers")
  print(pos)
  if (length(pos)>0){
    catchments$db_tra[i] <- mean(models2$b1.t2[pos] -  models2$b1.t1[pos],na.rm=TRUE)
  }
}
plot(catchments["db_tra"])


catchments["db_sol"] <- NA
for (i in 1:length(catchments$SITE_REF)){
  sr <- as.character(catchments$SITE_REF[i])
  print(sr)
  pos <- which(models2$site.id == sr & models2$solClass == "Solids")
  print(pos)
  if (length(pos)>0){
    catchments$db_sol[i] <- mean(models2$b1.t2[pos] -  models2$b1.t1[pos],na.rm=TRUE)
  }
}
plot(catchments["db_sol"])


catchments["db_turb"] <- NA
for (i in 1:length(catchments$SITE_REF)){
  sr <- as.character(catchments$SITE_REF[i])
  print(sr)
  pos <- which(models2$site.id == sr & models2$solute == "Turbidity.NTU.NTU" )
  print(pos)
  if (length(pos)>0){
    catchments$db_turb[i] <- mean(models2$b1.t2[pos] -  models2$b1.t1[pos],na.rm=TRUE)
  }
}
plot(catchments["db_turb"])


catchments["chemos"] <- NA
for (i in 1:length(catchments$SITE_REF)){
  sr <- as.character(catchments$SITE_REF[i])
  print(sr)
  pos <- which(models2$site.id == sr & ((models2$mod.t1 == "A" & models2$mod.t2 == "A") | (models2$mod.t1 == "C" & models2$mod.t2 == "C") | (models2$mod.t1 == "D" & models2$mod.t2 == "D")))
  print(pos)
  if (length(pos) > 0){
    catchments$chemos[i] <- sum(((models2$b1.t1[pos] < 0 & ((models2$b1.t2[pos] -  models2$b1.t1[pos]) >= 0) & ((models2$b1.t2[pos] -  models2$b1.t1[pos]) < -2*models2$b1.t1[pos])) | 
                                        (models2$b1.t1[pos] > 0 & ((models2$b1.t2[pos] -  models2$b1.t1[pos]) <= 0) & ((models2$b1.t2[pos] -  models2$b1.t1[pos]) > -2*models2$b1.t1[pos]))
                                      ),na.rm=TRUE)/length(pos)*100
  }
}
plot(catchments["chemos"])

catch.sp <- as_Spatial(catchments)
aus.sp <- as_Spatial(Aust)
library(raster)
aus.sp2 <- crop(aus.sp,extent(c(113.5,123,-35.5,-27.5)))

cols <- bpy.colors(n = 64, cutoff.tails = 0.1, alpha = 1.0) 



plot(reset = FALSE, catchments[-which(is.na(catchments["PET"])),"PET"], key.pos = 4,ylim=c(-35.5,-27.7),xlim=c(115,121.8))
plot(aus.sp2,add=TRUE)
plot(reset = FALSE, catchments[-which(is.na(catchments["P"])),"P"], key.pos = 4,ylim=c(-35.5,-27.7),xlim=c(115,121.8))
plot(aus.sp2,add=TRUE)
clim <- readRDS("D:\\WIRT\\ProcessedClimate\\climDBannual.RDS")
trends <- unlist(lapply(clim,FUN = function(x) coefficients(lm(x$Precip~x$year))[2]))
nms.trends <- unlist(strsplit(names(trends),split=".x\\$year"))
catchments$dP <- NA
for (i in 1:length(nms.trends)){
  pos <- which(as.character(catchments$SITE_REF) == nms.trends[i])
  if (length(pos)>0) catchments$dP[pos] <- trends[i]
}
plot(reset = FALSE, catchments["dP"], key.pos = 1,ylim=c(-35.5,-27.7),xlim=c(115,121.8), main = "Rainfall trend (mm/year)")
plot(aus.sp2,add=TRUE)
writeOGR(aus.sp2,"F:\\GIS\\WorldSimpleShapeFile\\SWWA.shp","SWWA",driver = "ESRI Shapefile")

SWWA <- readOGR("F:\\GIS\\WorldSimpleShapeFile\\SWWA_coast.shp")

par(mfrow=c(3,1))
plot(reset = FALSE, catchments["PET"], key.pos = 1,ylim=c(-35.5,-27.7),xlim=c(115,121.8), main = "Potential ET (mm/year)")
plot(SWWA,add=TRUE)

plot(reset = FALSE, catchments["dP"], key.pos = 1,ylim=c(-35.5,-27.7),xlim=c(115,121.8), main = "Rainfall (mm/year)")
plot(SWWA,add=TRUE)

SWWA <- st_read("F:\\GIS\\WorldSimpleShapeFile\\SWWA_coast.shp")

library(ggplot2)
catc <- ggplot() + geom_sf(data = Aust, color = "dark grey", fill = NA) + 
  geom_sf(data = catchments, color = NA, fill = "black") + 
  coord_sf(crs = st_crs(catchments),datum=NA) + 
  theme(panel.ontop = TRUE, 
        panel.grid = element_blank(), 
        line = element_blank(), 
        rect = element_blank()
  )
pet <- ggplot() + geom_sf(data = catchments, aes(fill = PET), lwd=0.1) + 
  scale_fill_viridis_c() +
  geom_sf(data = SWWA, color = "dark grey", fill = NA) + 
  coord_sf(crs = st_crs(catchments),datum=NA) + 
  theme(panel.ontop = TRUE, 
        panel.grid = element_blank(), 
        line = element_blank(), 
        rect = element_blank(),, 
        legend.title = element_blank(),
        legend.position = c(1.1,0.5),
        legend.text.align = 1
        ) + labs(title = "Potential ET (mm/year)")
p <- ggplot() + geom_sf(data = catchments, aes(fill = P), lwd=0.1) +
  scale_fill_viridis_c() +
  geom_sf(data = SWWA, color = "dark grey", fill = NA) +
  coord_sf(crs = st_crs(catchments),datum=NA) + 
  theme(panel.ontop = TRUE, 
        panel.grid = element_blank(), 
        line = element_blank(), 
        rect = element_blank(),, 
        legend.title = element_blank() ,
        legend.position = c(1.1,0.5),
        legend.text.align = 1) + 
  labs(title = "Precipitation (mm/year)")
dp <- ggplot() + geom_sf(data = catchments, aes(fill = dP), lwd=0.1) + 
  scale_fill_viridis_c() +
  geom_sf(data = SWWA, color = "dark grey", fill = NA) +
  coord_sf(crs = st_crs(catchments),datum=NA) + 
  theme(panel.ontop = TRUE, 
        panel.grid = element_blank(), 
        line = element_blank(), 
        rect = element_blank(), 
        legend.title = element_blank(),
        legend.position = c(1.1,0.5),
        legend.text.align = 1) + 
  labs(title = expression(paste("Precipitation Trend (mm/",year^2,")")))

par(oma=c(3,3,3,3),mar=c(0,0,0,0))
gridExtra::grid.arrange(catc,pet,p,dp,nrow=2)


plot(reset = FALSE, catchments["dP"], key.pos = 4,ylim=c(-35.5,-27.7),xlim=c(115,121.8), main = "Rainfall trend (mm/year)")
plot(SWWA,add=TRUE)

library(vioplot)
models2$b2.t1[is.na(models2$b2.t1)] <- models2$b1.t1[is.na(models2$b2.t1)]
models2$b2.t2[is.na(models2$b2.t2)] <- models2$b1.t2[is.na(models2$b2.t2)]
set.seed(13)
cols <- sample(brewer.pal(12, "Paired"), 12,replace = FALSE)[1:length(mod.types)]
par(mfrow=c(2,1), mar=c(3,4,0.5,0.5),oma=c(0,0,0,0))
vioplot(b1.t1~mod.t1,data=models2, col = cols, side = "left", ylim=c(-1,1.5), plotCentre = "line", colMed = NA,pchMed = NA,rectCol = NA)
vioplot(b1.t2~mod.t2,data=models2, col = cols, side = "right", add= TRUE, plotCentre = "line", colMed = "black", las=1,pchMed = 2,rectCol = NA)
mtext(side=2, las=1, line = 2.5, expression(paste(b[1])))
abline(h=0,lty=2)
legend("topleft",bty="n",legend = "(a)",adj=1)
vioplot(b2.t1~mod.t1,data=models2, col = cols, side = "left", ylim=c(-1,1.5), plotCentre = "line", colMed = "black",pchMed = 1,rectCol = NA)
vioplot(b2.t2~mod.t2,data=models2, col = cols, side = "right", add= TRUE, plotCentre = "line", colMed = "black", las=1,pchMed = 2,rectCol = NA)
abline(h=0,lty=2)
legend("topleft",bty="n",legend = "(b)",adj=1)
mtext(side=2, las=1, line = 2.5, expression(paste(b[2])))


#lONG TERM CHANGES IN c-q FOR SOME OF LARGEST DATASETS
setwd("D:\\WIRT\\Processed_CQDiscrete")
bayesResults <- readRDS("bayesResults.RDS")

#First identify lengths of each record
nmB <- names(bayesResults)
site.sol.list <- vector("list",length = length(nmB))
for (i in 1:length(nmB)){
  sols <- names(bayesResults[[nmB[i]]])
  if (length(sols)>0){
    solNumbers <- vector("numeric",length = length(sols))
  for (j in 1:length(sols)){
      solNumbers[j] <- 
        length(bayesResults[[nmB[i]]][[sols[j]]]$data$x)
  }
    site.sol.list[[i]] <- solNumbers
  }
}

#distribution of sample numers per C-Q
hist(log10(unlist( site.sol.list)))



