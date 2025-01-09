#load libraries needed

library(fitdistrplus)
library(mclust)
library(methods)
library(plyr)
library(selectspm)
library(spatstat)
library(spatstat.utils)
library(VGAM)
library(XML)
library(ggplot2)

source("stoloniferousFunctions.R")

ic.data<-list()

for(i in 1:20) 
{
	ic.data[[i]]<-read.delim(file.choose()) #repeat for 
}
#create ppp

ic.species.final<-list()
for(i in 1:20)
{
	plot(ic.data[[i]][,1:2])
	win.ic<-clickpoly(add=TRUE)
	ic.species.final[[i]]<-ppp(ic.data[[i]][,1],ic.data[[i]][,2],marks=ic.data[[i]][,3],window=win.ic)

}

names(ic.species.final)<-c("BC_Fractofusus" ,"D_Bradgatia","D_Pectinfrons","D_Fractofusus" ,"E_bradgatia","E_beothukis","E_CharniodiscusP" ,"E_CharniodiscusS" ,"E_Plum", "E_Primo2", "E_Fractofusus" , "LMP_Charniid1" ,"LMP_Charniid2" , "LMP_OstrichFeather" "G_Bradgatia" ,"StS_Charnia","BB_Charniodiscus","BB_Charnia", "BB_Plum","H14_Fractofusus") 

#calculate observed PCFs. 
#extract info needed for res.ic.f

get.ic.info.check<-function(ppp1)#need to ensure that input is a resonable e.g. m scale
{
	ic.info<-matrix(NA,1,9)
	colnames(ic.info)<-as.matrix(c("Count","MaxHeight","NoSizeclasses","PerInSmallest","IQR","LH*", "clusterRadius","ClusterSize/maxHeight","1-MinPCF"))
	
	#specimen numbers
	ic.info[,1]<-ppp1$n
	#maximum size
	ic.info[,2]<-max(marks(ppp1))
	#number of size classes
	mc.res<-Mclust(marks(ppp1))
	mc.class<-table(mc.res$classification)
	ic.info[,3]<-mc.res$G
	#percentage in the smallest one
	ic.info[,4]<-mc.class[1]/ppp1$n
	#IQR of isotrophy
	po.ppp1<-pairorient(ppp1,0,0.2,sigma=3,unit="degree") #where 0.2 is the max distance apart. 
	ic.info[,5]<-summary(po.ppp1$border)[5]-summary(po.ppp1$border)[2]
	#lhtests
	ic.info[,6]<-sum(Gest(ppp1)$han-Gest(ppp1)$theo)
	#cluster radius given by sigma
	ic.info[,7]<-kppm(unmark(ppp1)~1, clusters = "Thomas",rmax=0.2)$clustpar[2] #limits model fitting to the first 0.2m
	#max spcimen size/cluster size
	ic.info[,8]<-ic.info[,7]/ic.info[,2]
	#intra comp (1-pcf_min)
	test.env<-envelope(ppp1,pcf,nsim=19,nrank=1,savefuns = TRUE)
	lo1<-test.env$lo
	obs1<-test.env$obs
	ic.info[,9]<-1-min(obs1)
	return(ic.info)
	}


get.all.ic.info.check<-function(bp.species1)
{
	res<-list()
	for(i in 1:length(bp.species1))
	{
		res[[i]]<-get.ic.info.check(bp.species1[[i]])


	}
	res1<-do.call(rbind,res)
	row.names(res1)<-names(bp.species1)
	return(res1)
}

res.ic<-get.all.ic.info.check(ic.species.final)

####
#from plotting the PCFs (see below code) we visually identify intraspecific competition as follows:
ic.bi2<-c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,0,0)
  
######

#plot bar plot and plot best icq roses as per plot

p<-ggplot(as.data.frame(cbind(ic.bi2,res.ic)),aes(factor(ic.bi2),res.ic[,5]))
p+geom_boxplot()
p + geom_boxplot() + geom_jitter(height = 0, width = 0.1)

#example Directional plots
par(mfrow=c(5,1),mai=c(0.1,0.1,0.1,0.1))
rose(pairorient(ic.species.final$LMP_OstrichFeather,0,0.2,sigma=20,unit="degree"),col=col.order[14],main=names1[14])
rose(pairorient(ic.species.final$E_beothukis,0,0.2,sigma=20,unit="degree"),col=col.order[6],main=names1[6])
rose(pairorient(ic.species.final$D_Pectinfrons,0,0.2,sigma=20,unit="degree"),col=col.order[3],main=names1[3])
rose(pairorient(ic.species.final$E_bradgatia,0,0.2,sigma=20,unit="degree"),col=col.order[5],main=names1[5])
rose(pairorient(ic.species.final$H14_Fractofusus,0,0.2,sigma=40,unit="degree"),col=col.order[1],main=names1[20])

#Mann Whitney tests
wilcox.test(res.ic[,5]~ic.bi2)
wilcox.test(res.ic[,5]+res.ic[,8]~ic.bi2)

#regression analyses
summary(lm(res.ic.f[,14]~diff.tc.h[,3]+res.ic.f[,c(4,5,7,8)]))
AIC(lm(res.ic.f[,14]~diff.tc.h[,3]+res.ic.f[,c(4,5,7)]))
summary(lm(res.ic.f[,14]~diff.tc.h[,3]+res.ic.f[,c(4,5,7)]))
AIC(lm(res.ic.f[,14]~diff.tc.h[,3]+res.ic.f[,c(7)]+res.ic.f[,c(4)]))
summary(lm(res.ic.f[,14]~diff.tc.h[,3]+res.ic.f[,c(7)]+res.ic.f[,c(4)]))
AIC(lm(res.ic.f[,14]~diff.tc.h[,3]+res.ic.f[,c(5,7)]))
summary(lm(res.ic.f[,14]~diff.tc.h[,3]+res.ic.f[,c(5,7)]))
AIC(lm(res.ic.f[,14]~diff.tc.h[,3]+res.ic.f[,c(7)]))
summary(lm(res.ic.f[,14]~diff.tc.h[,3]+res.ic.f[,c(7)]))
AIC(lm(res.ic.f[,14]~diff.tc.h[,3]))
summary(lm(res.ic.f[,14]~diff.tc.h[,3]))
AIC(lm(res.ic.f[,14]~res.ic.f[,4]))
summary(lm(res.ic.f[,14]~res.ic.f[,4]))
AIC(lm(res.ic.f[,14]~res.ic.f[,7]))
summary(lm(res.ic.f[,14]~res.ic.f[,7]))


###

#plot of modelling

mypalette<-brewer.pal(11,"Spectral")
image(1:11,1,as.matrix(1:11),col=mypalette,xlab="Greens (sequential)",
       ylab="",xaxt="n",yaxt="n",bty="n")
col.order<-c(11,10,9,11,10,8,2,1,3,4,11,8,6,7,10,5,1,5,3)
pch.order<-c(5,15,15,15,16,16,16,16,16,16,16,17,17,17,6,2,18,18,18,1)

plot(res.ic.f[,14],simulate(best.ic.model)[,1],xlab="Min PCF",ylab="Simulated ",ylim=c(-0.5,1),xlim=c(-0.1,0.7),pch=pch.order,col=mypalette[col.order],xaxs="i", yaxs="i",cex=1.2)

legend(0.57,0.20 ,title="Taxa",   c("Beothukis","Bradgatia","Charniodiscus","Charnia","Culmofrons","Fractofusus","Ostrich Feather", "Pectinifrons","Primocandelabrum"), fill=mypalette[c(8,10,1,5,6,11,7,9,3)], horiz=FALSE, cex=0.7)

legend(0.47,0.20, title="Surfaces",   legend=c("Bed B","Bristy Cove","D surface","E surface","G surface","H14","LMP", "St Shotts"), pch=c(18,5,15,16,6,1,17,2),horiz=FALSE, cex=0.7)

abline(0,1,col="grey") #line for when simulated = observed



#plot pcf plots. The smoothed plots for Figure 1 are given by the #,funargs=c(stoyan=0.8)) argument, and the non-smoothed for Extended Figure 2, and the PCF values used for the analyses are not smoothed. . 

par(mfrow=c(3,3))

#plum, chaP

inter.e.plum.chaP<-superimpose(plum=unmark(ic.species.final$E_Plum),chaP=unmark(ic.species.final$E_CharniodiscusP))
pcf.e.plum.chaP<-envelope(inter.e.plum.chaP,pcfcross, nsim=999,nrank=50)#,funargs=c(stoyan=0.8))

plot(pcf.e.plum.chaP,main="E: Plumeropriscum and CharniodiscusP",ylim=c(0.5,2),legend=FALSE,xaxs="i",yaxs="i",xlab="Distance",ylab="PCF") #chaP primo1
plot(pcf(ic.species.final$E_Plum,stoyan=0.8),add=TRUE,col=2,legend=FALSE)
plot(pcf(ic.species.final$E_CharniodiscusP,stoyan=0.8),add=TRUE,col="orange",legend=FALSE)
#primo1, chaP

inter.e.primo.chaP<-superimpose(primo=unmark(ic.species.final$E_Primo2),chaP=unmark(ic.species.final$E_CharniodiscusP))
pcf.e.primo.chaP<-envelope(inter.bb.primo.chaP,pcfcross, nsim=999,nrank=50)

plot(pcf.e.primo.chaP,main="E: Primocandelabrum and CharniodiscusP ",ylim=c(0.5,2),legend=FALSE,xaxs="i",yaxs="i",xlab="Distance",ylab="PCF") #chaP primo1
plot(pcf(ic.species.final$E_Primo2),add=TRUE,col=2,legend=FALSE)
plot(pcf(ic.species.final$E_CharniodiscusP),add=TRUE,col="orange",legend=FALSE)


#plum, chaS
inter.e.plum.chaS<-superimpose(plum=unmark(ic.species.final$E_Plum),chaS=unmark(ic.species.final$E_CharniodiscusS))
pcf.e.plum.chaS<-envelope(inter.e.plum.chaS,pcfcross, nsim=999,nrank=50)

plot(pcf.e.plum.chaS,main="E: Plumeropriscum and CharniodiscusS ",ylim=c(0.5,2),legend=FALSE,xaxs="i",yaxs="i",xlab="Distance",ylab="PCF") #chaP primo1
plot(pcf(ic.species.final$E_Plum),add=TRUE,col=2,legend=FALSE)
plot(pcf(ic.species.final$E_CharniodiscusS),add=TRUE,col="orange",legend=FALSE)
 

#primo, chaS

inter.e.primo.chaS<-superimpose(primo=unmark(ic.species.final$E_Primo2),chaS=unmark(ic.species.final$E_CharniodiscusS))
pcf.e.primo.chaS<-envelope(inter.e.primo.chaS,pcfcross, nsim=999,nrank=50)

plot(pcf.e.primo.chaS,main="E: Primocandelabrum and CharniodiscusS ",ylim=c(0.5,2),legend=FALSE,xaxs="i",yaxs="i",xlab="Distance",ylab="PCF") #chaP primo1
plot(pcf(ic.species.final$E_Primo2),add=TRUE,col=2,legend=FALSE)
plot(pcf(ic.species.final$E_CharniodiscusS),add=TRUE,col="orange",legend=FALSE)

#fra, chaS

#inter.bb.fra.chaS<-superimpose(fra=unmark(ic.species.final$E_Fractofusus),chaS=unmark(ic.species.final$E_CharniodiscusS))
#pcf.bb.fra.chaS<-envelope(inter.bb.fra.chaS,pcfcross, nsim=999,nrank=50)

plot(pcf.bb.fra.chaS,main="E: Fractofusus and CharniodiscusS ",ylim=c(0.75,2),legend=FALSE,xaxs="i",yaxs="i",xlab="Distance",ylab="PCF") #chaP primo1
plot(pcf(ic.species.final$E_Fractofusus),add=TRUE,col=2,legend=FALSE)
plot(pcf(ic.species.final$E_CharniodiscusS),add=TRUE,col="orange",legend=FALSE)

#fra, chaP
#inter.bb.fra.chaP<-superimpose(fra=unmark(ic.species.final$E_Fractofusus),chaP=unmark(ic.species.final$E_CharniodiscusP))
#pcf.bb.fra.chaP<-envelope(inter.bb.fra.chaP,pcfcross, nsim=999,nrank=50)

plot(pcf.bb.fra.chaP,main="E: Fractofusus and CharniodiscusP ",ylim=c(0.75,2),legend=FALSE,xaxs="i",yaxs="i",xlab="Distance",ylab="PCF") #chaP primo1
plot(pcf(ic.species.final$E_Fractofusus),add=TRUE,col=2,legend=FALSE)
plot(pcf(ic.species.final$E_CharniodiscusP),add=TRUE,col="orange",legend=FALSE)

#Bed B Charnio Charnia ##shows heteromyopia. 

#inter.bb.char.cha<-superimpose(char=unmark(ic.species.final$BB_Charnia),cha=unmark(ic.species.final$BB_Charniodiscus))
#pcf.bb.char.cha<-envelope(inter.bb.char.cha,pcfcross, nsim=999,nrank=50)

 plot(pcf.bb.char.cha,main="BB: Charnia and Charniodiscus ",ylim=c(0.5,2),legend=FALSE,xaxs="i",yaxs="i",xlab="Distance",ylab="PCF") #chaP primo1
 plot(pcf(ic.species.final$BB_Charniodiscus),add=TRUE,col=2,legend=FALSE)
 plot(pcf(ic.species.final$BB_Charnia),add=TRUE,col="orange",legend=FALSE)
 

#LMP culmo ost


inter.culmo.ost<-superimpose(culmo=unmark(ic.species.final$LMP_Charniid2),ost=unmark(ic.species.final$LMP_OstrichFeather))
pcf.culmo.ost<-envelope(inter.culmo.ost,pcfcross, nsim=999,nrank=50)

plot(pcf.culmo.ost,main="LMP: Culmofrons - Ostrich Feathers",ylim=c(0.5,2),legend=FALSE,xaxs="i",yaxs="i",cex=5,xlab="Distance",ylab="PCF") #chaP fra
plot(pcf(ic.species.final$LMP_Charniid2),add=TRUE,col=2,legend=FALSE,cex=5)
plot(pcf(ic.species.final$LMP_OstrichFeather),add=TRUE,col="orange",legend=FALSE,cex=5)


#G
#pcf.g.bra<-envelope(ic.species.final$G_Bradgatia,pcf, nsim=999,nrank=50)

plot(pcf.g.bra,main="G: Bradgatia",ylim=c(0.25,2),legend=FALSE,xaxs="i",yaxs="i",cex=5,xlab="Distance",ylab="PCF") #chaP fra

############## plots for teh 





