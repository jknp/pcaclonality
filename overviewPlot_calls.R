#### Plotting of segmenetd values from calls object
overviewPlot_calls<-function(calls, scaling=1){
#	Get data from calls object
	Seg<-calls(calls)
#	Scaling of data
	Seg[which(Seg>scaling)]<- scaling
	Seg[which(Seg< -scaling)]<- -scaling
	
	color.palette = colorRampPalette(c("blue", "white", "red"))

	regChr<-chromosomes(calls)
	chrInd <- rep(0, length(regChr))
	chrInd[(regChr %% 2 == 0)] <- 1
	chrColor <- rep("black", length(regChr))
	chrColor[(regChr %% 2 == 0)] <- c("white")
	
	Y <- rep(FALSE, length(regChr))
	for (i in 2:length(regChr)){
		if ((regChr[i-1] != regChr[i])){ Y[i] <- TRUE }
	}
	Y[1] <- TRUE
	beginChr <- rep("", length(regChr))
	beginChr[Y] <- regChr[Y]
	
	def.par <- par
	fl <- layout(matrix(c(1,2,1,2,1,2), 2, 2, byrow = TRUE), width=c(1,9))
	# layout.show(fl)
	
	par(mar=c(3,2,4,0))
	image(z=matrix(chrInd, nrow=1), xaxt="n", yaxt="n", col=c("Black", "white"))
	axis(2, at=(which(Y)-1)/(length(Y)-1), labels=regChr[Y], tick=FALSE, las=1)
	par(mar=c(3,1,4,1))
	image(z=t(Seg), xaxt="n", yaxt="n", col=color.palette(25), main="calls values 0=neutral, 1=gain, -1=loss")
	par(def.par)
	return(invisible(NULL))
	gc()	
}
