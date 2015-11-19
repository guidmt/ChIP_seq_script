#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
memory.limit(2000)
require(Cairo)

#outputoutput<-"heatmap_sort_SUZ12_MYC.png"
print(args[1])
print(args[2])

sampleinfo<- args[1]
output<-args[2]
fc<-args[3]
print(fc)

sampleinfo<-read.table(file=sampleinfo,sep="\t",stringsAsFactors=F,header=F)
###################################################################
##BEFORE I READ THE FILES THAT WILL USED TO RANK THE OTHERS HEATMAP
###################################################################
n<-1

#LIF
sample_a<-sampleinfo[n,1]
#OHT
sample_b<-sampleinfo[n,2]

#LIF
sample_a_label<-sampleinfo[n,3]
#OHT
sample_b_label<-sampleinfo[n,4]

color<-as.character(sampleinfo[n,5])

print(color)

#LIF
tab_sample_a<-read.table(file=as.character(sample_a),sep="\t",stringsAsFactors=F,header=T)

#I select only the sample with the ChIP signals, the number of columns change by the bin-size
tab_sample_a<-tab_sample_a[,c(803:ncol(tab_sample_a))]

#OHT
tab_sample_b<-read.table(file=as.character(sample_b),sep="\t",stringsAsFactors=F,header=T)

#I select only the sample with the ChIP signals, the number of columns change by the bin-size
tab_sample_b<-tab_sample_b[,c(803:ncol(tab_sample_b))]

###order always for the OHT

#order_genes<-tab_sample_b[order(tab_sample_b[,center_of_tss],decreasing=FALSE),1]
#in this case the center of tss is 401
##Qui ordino per OHT e utilizzando es suz12, questo ordinamento verrà utilizzato in tutto lo script
center_of_tss<-401
order_matrix<-order(tab_sample_b[,401],decreasing=FALSE)

tab_sample_a_order<-tab_sample_a[order_matrix,]
tab_sample_b_order<-tab_sample_b[order_matrix,]

rm(tab_sample_a)
rm(tab_sample_b)

range_colors<-data.frame(quantile(apply(tab_sample_b_order,2,mean)))

ranged_data<-c(min(tab_sample_b_order),range_colors[,1],max(tab_sample_b_order))

###
### Il numero dei ranges deve essere superiore al numero di colori di 1 
###

#if 100% of data are lower than 0.80

#pdf(file=output,height=10,width=30)

CairoPNG(filename=output, width=2000, height=1000)

par(mfrow=c(1,15),mar=c(4,0.8,2,0.8))
###
### Il numero dei ranges deve essere superiore al numero di colori di 1 
###
	#sample a
	range_colors<-data.frame(quantile(apply(tab_sample_a_order,2,mean)))

	ranged_data<-c(min(tab_sample_a_order),range_colors[,1],max(tab_sample_a_order))
	print(length(ranged_data))

	#if 100% of data are lower than 0.80
	if(ranged_data[6]<=0.80){

	color.scale<-c(c("white","white","white","white","white"), paste(color,3,sep="")) } else {

	color.scale<-c(c("white","white","white","white"), paste(color,seq(2:3),sep=""))

	}

	image(t(tab_sample_a_order), x=1:ncol(tab_sample_a_order), y=1:nrow(tab_sample_a_order),col=color.scale,breaks=ranged_data,axes=F,ylab="",xlab="",main=sample_a_label)
	axis(1, c(1, center_of_tss,ncol(tab_sample_a_order)), labels=c("-10", "0", "10"))

	#sample b
	range_colors<-data.frame(quantile(apply(tab_sample_b_order,2,mean)))

	ranged_data<-c(min(tab_sample_b_order),range_colors[,1],max(tab_sample_b_order))
	print(length(ranged_data))

	#if 100% of data are lower than 0.80
	if(ranged_data[6]<=0.80){

	color.scale<-c(c("white","white","white","white","white"), paste(color,3,sep="")) } else {

	color.scale<-c(c("white","white","white","white"), paste(color,seq(2:3),sep=""))

	}

	image(t(tab_sample_b_order), x=1:ncol(tab_sample_b_order), y=1:nrow(tab_sample_b_order),col=color.scale,breaks=ranged_data,axes=F,ylab="",xlab="",main=sample_b_label)
	axis(1, c(1, center_of_tss,ncol(tab_sample_b_order)), labels=c("-10", "0","10"))

	##estimation of fold-change
	difference_tab<-tab_sample_b_order-tab_sample_a_order

	difference_tab[difference_tab<fc]<-0

	main_diff<-paste(sample_b_label,sample_a_label,sep="/")

	image(t(difference_tab), x=1:ncol(difference_tab), y=1:nrow(difference_tab),col=color.scale,breaks=ranged_data,axes=F,ylab="",xlab="",main=main_diff)
	axis(1, c(1, center_of_tss,ncol(difference_tab)), labels=c("-10", "0","10"))

############################

n<-0
#i delete the reference
sampleinfo<-sampleinfo[-1,]

for(n in 1:nrow(sampleinfo)){

#LIF
sample_a<-sampleinfo[n,1]
#OHT
sample_b<-sampleinfo[n,2]
print(sample_a)
#LIF
sample_a_label<-sampleinfo[n,3]
#OHT
sample_b_label<-sampleinfo[n,4]

color<-as.character(sampleinfo[n,5])

#LIF
tab_sample_a<-read.table(file=as.character(sample_a),sep="\t",stringsAsFactors=F,header=T)
tab_sample_a<-tab_sample_a[,c(803:ncol(tab_sample_a))]

#OHT
tab_sample_b<-read.table(file=as.character(sample_b),sep="\t",stringsAsFactors=F,header=T)
tab_sample_b<-tab_sample_b[,c(803:ncol(tab_sample_b))]

tab_sample_a_order<-tab_sample_a[order_matrix,]
tab_sample_b_order<-tab_sample_b[order_matrix,]

	#sample a
	range_colors<-data.frame(quantile(apply(tab_sample_a_order,2,mean)))
	ranged_data<-c(min(tab_sample_a_order),range_colors[,1],max(tab_sample_a_order))

	if(ranged_data[6]<=0.80){

	color.scale<-c(c("white","white","white","white","white"), paste(color,3,sep="")) } else {

	color.scale<-c(c("white","white","white","white"), paste(color,seq(2:3),sep=""))

	}

	image(t(tab_sample_a_order), x=1:ncol(tab_sample_a_order), y=1:nrow(tab_sample_a_order),col=color.scale,breaks=ranged_data,axes=F,ylab="",xlab="",main=sample_a_label)
	axis(1, c(1, center_of_tss,ncol(tab_sample_a_order)), labels=c("-10", "0", "10"))

	#sample b
	range_colors<-data.frame(quantile(apply(tab_sample_b_order,2,mean)))

	ranged_data<-c(min(tab_sample_b_order),range_colors[,1],max(tab_sample_b_order))

	print(length(ranged_data))

	if(ranged_data[6]<=0.80){

	color.scale<-c(c("white","white","white","white","white"), paste(color,3,sep="")) } else {

	color.scale<-c(c("white","white","white","white"), paste(color,seq(2:3),sep=""))

	}


	image(t(tab_sample_b_order), x=1:ncol(tab_sample_b_order), y=1:nrow(tab_sample_b_order),col=color.scale,breaks=ranged_data,axes=F,ylab="",xlab="",main=sample_b_label)
	axis(1, c(1, center_of_tss,ncol(tab_sample_b_order)), labels=c("-10", "0","10"))
	#estimate fold-change
	difference_tab<-tab_sample_b_order-tab_sample_a_order

	#the values with a fold chagne <2 are represent with 0 
	difference_tab[difference_tab<fc]<-0

	range_colors<-data.frame(quantile(apply(difference_tab,2,mean)))
	ranged_data<-c(min(difference_tab),range_colors[,1],max(difference_tab))

	main_diff<-paste(sample_b_label,sample_a_label,sep="/")

	image(t(difference_tab), x=1:ncol(difference_tab), y=1:nrow(difference_tab),col=color.scale,breaks=ranged_data,axes=F,ylab="",xlab="",main=main_diff)
	axis(1, c(1, center_of_tss,ncol(difference_tab)), labels=c("-10", "0","10"))


}


dev.off()
