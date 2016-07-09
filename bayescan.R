source("/Users/Ying/Desktop/bin/BayeScan2.1/R functions/plot_R.r")
setwd("/Users/Ying/Desktop/analysis/greenbul/RAD/all-rerun/population-1-m0/")

#plot outliers with q-value lower than 5%
plot_bayescan("batch_4.bayescan.forest.eco_fst.txt",FDR=0.05)
results<-plot_bayescan("batch_4.bayesc_fst.txt",FDR=0.05)
results$outliers
results$nb_outliers

#plotting posterior distribution
mydata=read.table("batch_4.bayescan-fe.sel",colClasses="numeric")
parameter="Fst1"
plot(density(mydata[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution"))
parameter="alpha1"
parameter="logL"

library(boa)
boa.hpd(mydata[[parameter]],0.05)
