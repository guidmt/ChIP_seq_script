sampleinfo<-read.table(file="sampleinfo_myc_averageprofiles.txt",sep="\t",stringsAsFactors=F,header=T)

group_category<-unique(sampleinfo$Group)


for(i in group_category){

sub_sampleinfo<-sampleinfo[sampleinfo$Group==i,]

        output<-paste(i,"_homer_bin25_averageprofiles.norm.pdf",sep="")

        pdf(output,width=20,height=5)
        
        par(mfrow=c(1,4))

	for(r in 1:nrow(sub_sampleinfo)){	

	##read OHT file
	samplea<-read.table(file=sub_sampleinfo[r,1],sep="\t",stringsAsFactors=F,header=T)
	#samplea<-samplea[,-1]
	samplea_chip<-samplea[,5] - samplea[,2]
	##read LIF file
	sampleb<-read.table(file=sub_sampleinfo[r,2],sep="\t",stringsAsFactors=F,header=T)
	#sampleb<-sampleb[,-1]
	sampleb_chip<-sampleb[,5] - samplea[,2]
	
	label1<-sub_sampleinfo[r,3]
	label2<-sub_sampleinfo[r,4]
	title<-sub_sampleinfo[r,5]
	min.scale<-sub_sampleinfo[r,7]
	max.scale<-sub_sampleinfo[r,8]

#	mean.samplea<-apply(samplea_chip,2,mean)
#	mean.sampleb<-apply(sampleb_chip,2,mean)

#	matplot(mean.samplea,type="l",lty=1,xaxt="n",lwd=1.5,ylab="Normalized Reads Coverage",col="blue4",main=title,ylim=c(min(min.scale),max(max.scale)))

#	matplot(mean.sampleb,type="l",lty=1,xaxt="n",lwd=1.5,col="red4",add=T)
#	axis(1,at=c(1,ncol(samplea_chip)/2,ncol(samplea_chip)),labels=c("-3","0","3"))
#	legend("topleft",legend=c("Myc","LIF"),col=c("blue4","red4"),lty=c(1,1),bty = "n",lwd=2, cex=1.2)


	matplot(samplea_chip,type="l",lty=1,xaxt="n",lwd=1.5,ylab="Normalized Reads Coverage",col="#00C0AF",main=title,ylim=c(min(min.scale),max(max.scale)))
	matplot(sampleb_chip,type="l",lty=1,xaxt="n",lwd=1.5,col="#F8766D",add=T)
	axis(1,at=c(1,length(samplea_chip)/2,length(samplea_chip)),labels=c("-3","0","3"))
	legend("topleft",legend=c("LIF","Myc"),col=c("#F8766D","#00C0AF"),lty=c(1,1),bty = "n",lwd=2, cex=1.2)


}

dev.off()

}
