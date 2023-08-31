PlotResp<-function(VIDSp,nplot,selval){
    library(ggplot2)
	x=1:100
	y=((100/x)*(1/sqrt(2*pi))*exp((-qnorm(1-(x/100))^2)/2))*c(sqrt(VIDSp))
	#y=(y-min(y))/max(y)
	forPlot=as.data.frame(cbind(x,y))
	names(forPlot)=c("Percent","ResponseToSelection")
	if(round(max(y)/10,0)<=0){limup=ceiling(max(y)/10)}else{limup=round(max(y)/10,0)}
	#p<-ggplot(forPlot, aes(Percent, ResponseToSelection)) + geom_point(colour = "blue") + scale_x_continuous(breaks=seq(0,100,5)[-1])
	p<-ggplot(forPlot, aes(Percent, ResponseToSelection)) + geom_point(colour = "blue") + scale_x_continuous(breaks=seq(0,100,5)[-1]) + scale_y_continuous(breaks=seq(0,max(y),limup))
	pdf(paste("RespSelPlot_",nplot,".pdf",sep=""))
	print(p)
	dev.off()
	jpeg(paste("RespSelPlot_",nplot,".jpg",sep=""),width = 680, height = 680)
	print(p)
	dev.off()	
	Rsel=y[selval]
	return(Rsel)
}