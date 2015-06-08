devtools::load_all();

FVEdata <- read.table("http://www.hsph.harvard.edu/fitzmaur/ala2e/fev1.txt", col.names=c('SubjectID', 'Height', 'Age', 'InitialHeight', 'InitialAge', 'LogFEV1'), skip=42  );

B = makePACEinputs(IDs= FVEdata$SubjectID, tVec=FVEdata$Age, yVec=FVEdata$LogFEV1);

y= B$Ly
t= B$Lt;

optns = CreateOptions()

 # FPCA checks the data validity for the PCA function. 
  if( CheckData(y,t) ){
    cat('FPCA has stopped.')
    return(FALSE);
  }  
 

  # FPCA sets the options structure that are still NULL
  optns = SetOptions(y, t, optns);

  
  # FPCA checks the options validity for the PCA function. 
  numOfCurves = length(y);
  if( CheckOptions(t, optns,numOfCurves) ){
    cat('FPCA has stopped.')
    return(FALSE);
  }


  # Bin the data (potentially):
  if ( optns$useBinnedData != 'OFF'){ 
      BinnedDataset <- GetBinnedDataset(y,t,optns)
      y = BinnedDataset$newy;
      t = BinnedDataset$newt; 
  }

  # Generate basic grids:
  # obsGrid:  the unique sorted pooled time points of the sample and the new data
  # regGrid: the grid of time points for which the smoothed covariance surface assumes values
  obsGrid = sort(unique( c(unlist(t), optns$newdata)));
  regGrid = seq(min(obsGrid), max(obsGrid),length.out = optns$nRegGrid);


  # Get the smoothed mean curve
  smcObj = GetSmoothedMeanCurve(y, t, obsGrid, regGrid, optns)
 

