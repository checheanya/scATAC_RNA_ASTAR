#install.packages("feather")
install.packages("arrow")
setwd("~")

library(feather)
library(fitdistrplus)
library(MASS)
library(actuar)
library(stats)
library(sfsmisc)

path <- "rna_seq_x.feather"
df <- arrow::read_feather(path)

# set up for plots
colors <- c(rainbow(5),grey(0.6))
dist.names <- c("Lognormal", "Pareto", "Burr", "Loglogistic", "Weibull", "Gamma")
sample.list <- list()
# normal distr are not gonna work for single cell cause we have here not normalized data --> we can norm within or bw cells --> different analysis
# plus the matrix is super sparse --> classical distr won't fit 
# maybe:
# start from ATAC - select data on peaks --> select top N genes --> see noise for them and distr stuff

# interactive fitting of different distributions for each sample
for (i in 1:ncol(df)) {
  samplei.i_raw<-df[,i]
  samplei.i_raw<-samplei.i_raw[samplei.i_raw>0]
  samplei.i <- (samplei.i_raw - min(samplei.i_raw) + 0.0001) / (max(samplei.i_raw) - min(samplei.i_raw) + 0.0002)
  
  ln<-fitdist(samplei.i,"lnorm")
  p<-fitdist(samplei.i,"pareto",start = list(shape = 1, scale = 200),lower=c(0,0))
  b<-fitdist(samplei.i,"burr",start = list(shape1 = 1, shape2 = 1, rate = 1),lower=c(0,0,0))
  ll<-fitdist(samplei.i,"llogis",lower=c(0,0))
  wb<-fitdist(samplei.i,"weibull",lower=c(0,0))
  g<-fitdist(samplei.i,"gamma",lower=c(0, 0),start=list(scale=1,shape=1))
  
  samplei.i.list<-list(ln,p,b,ll,wb,g)
  sample.list[[i]] <- samplei.i.list
}

# AIC table construction
AIC_df <- as.data.frame(matrix(nrow=ncol(tpm.df),ncol=6))
colnames(AIC_df) <- dist.names
for(i in 1:ncol(tpm.df)){
  samplei.i.list <- sample.list[[i]]
  AIC_df[i,] <- gofstat(samplei.i.list)$aic
}

# get min from AIC table
for(i in 1:nrow(AIC_df)){
  AIC_df$min_AIC[i]<-colnames(AIC_df)[which.min(AIC_df[i,1:6])]
}
rownames(AIC_df) <- colnames(tpm.df)
AIC_df$cutoff_value <- cutoff_val

######################## RUN IN ONE GO #######################

#here you are saving the plots. Make sure to run from png to dev.off() in one go
png(paste0("./DISTR_RNA.png"), unit="cm", height=5, width=6, res=300) #you can modify the unit, height, width and resolution as you see fit
par(mfrow=c(2,2), mgp=c(5,0.4,0), oma = c (1, 2, 1, 0.5), cex=0.35) # you can modify mfrow to however many plots you need

# only get the first replicate to plot
for(i in seq_along(sample.list)){
  # only get the first replicate to plot
  samplei.i.list <- sample.list[[i]]
  if( i %in% c(1) ) {addlegend=T} else{addlegend=F}
  if(i %in% c(1,2)) {left_mar=2} else {left_mar=0}       # space for y ticks
  if(i %in% c(2,4)) {bottom_mar=2} else {bottom_mar=0} # space for x ticks
  
  par(mar=c(2, 2, 2, 0))
  cdfcomp(samplei.i.list,xlogscale=TRUE,ylogscale=TRUE,do.points=F,fitcol=colors, #lwd=1,
          ylim=c(10^-3,1),xlim=c(10^-2,10^5), addlegend = F, # addlegend
          ann=F) #main=cell_types[j]
  title(colnames(df_tpm)[i], line = 0.2)
  abline(v= cutoff_val, lty=2)
  
  if(i==1) legend("bottomright",bty="n" ,col=colors,legend=dist.names,lty=seq(length(dist.names)))
  # legend(x=10^0.2, y = 10^-0.7, bty="n",legend=cell_types[i])
}
mtext("Cumulative Distribution Function", side = 2, outer = T, cex=0.5) 
mtext("TPM", side = 1, outer = T, cex=0.5) 
dev.off()

################ END OF THE 1ST PLOT, START OF THE 2ND ##################

png(paste0("./DISTR_RNA_qq.png"), unit="cm", height=5, width=6, res=300)
par(mfrow=c(2,2), mgp=c(5,0.4,0), oma = c (1, 2, 1, 0.5), cex=0.35)  #bottom, left, top, right #location of label, tick mark label, tick mark
for (i in 1:ncol(tpm.df)) {  
  samplei.i.list <- sample.list[[i]]
  if( i %in% c(1) ) {addlegend=T} else{addlegend=F}
  if(i %in% c(1,2)) {left_mar=2} else {left_mar=0}       # space for y ticks
  if(i %in% c(2,4)) {bottom_mar=2} else {bottom_mar=0} # space for x ticks
  
  par(mar=c(2, 2, 2, 0))
  qqcomp(samplei.i.list, xlogscale=TRUE, ylogscale=TRUE, xlim=c(10^-1,10^3), ylim=c(10^-1,10^4),
         fitpch=20, fitcol=colors ,addlegend = F, main=cell_types[j], ann=FALSE) #, fitpch="."
  title(colnames(df_tpm)[i], line = 0.2)
  # if(i %% 2 == 1) eaxis(side=2, at = c(10^-2, 1, 10^2, 10^4, 10^6), cex=0.8)
  # if(i %in% c(3,4)) eaxis(side=1, at = c( 10^-2, 1, 10^2, 10^4, 10^6), cex=0.8)
  abline(v=3, lty=2)
  if(i==1) legend("bottomright",bty="n" ,col=colors,legend=dist.names,lty=seq(length(dist.names)))
  
  # legend(x=10^-1.5, y=10^4.3,bty="n",legend=paste0("rep. ", replicates[i],"\n",samples[i]))
}
mtext("Empirical quantiles", side = 2, outer = T, cex=0.5)
mtext("Theoretical quantiles", side = 1, outer = T, cex=0.5)
dev.off()




