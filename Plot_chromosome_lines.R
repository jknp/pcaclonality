
Chr.names<-function(Chromo){
uni.chr<-unique(Chromo)
uni.chr<-uni.chr[1:22]
temp<-rep(0,length(uni.chr))
temp2<-rep(0,length(uni.chr))
temp3<-rep(0,length(uni.chr))
for (i in 1:length(uni.chr)){
	temp[i]<-max(which(uni.chr[i]==Chromo))
	temp2[i]<-min(which(uni.chr[i]==Chromo))
	temp3[i]<-temp2[i]+((temp[i]-temp2[i])/2)
}
for (i in 1:length(temp)){
	abline(v=temp[i], col="black", lty="dashed")
	axis(1,at=temp3[i],labels=uni.chr[i],cex.axis=1.2)
}}
