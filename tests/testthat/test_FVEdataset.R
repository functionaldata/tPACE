devtools::load_all()




FVEdata <- read.table("http://www.hsph.harvard.edu/fitzmaur/ala2e/fev1.txt", col.names=c('SubjectID', 'Height', 'Age', 'InitialHeight', 'InitialAge', 'LogFEV1'), skip=42  );

mySample = makePACEinputs(IDs= FVEdata$SubjectID, tVec=FVEdata$Age, yVec=FVEdata$LogFEV1);

y= mySample$Ly
t= mySample$Lt

# optns = CreateOptions()
# system.time(tmp <- FPCA(y, t, optns))
# tmp$sigma2


optns1 <- CreateOptions(kernel='rect')
system.time(tmp1 <- FPCA(y, t, optns1))
plot(tmp1$phi[, 1]) # off
createCorrPlot(tmp1, 'Smoothed', TRUE)
