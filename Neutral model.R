#Adam Burns - 2/10/2015
#aburns2@uoregon.edu
#From Burns et al. Contribution of neutral processes to the assembly of the gut microbial communities changes over host development
#Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics. Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity when stats=FALSE. For use in R.
#comun: A community table for communities of interest with local communities/samples as rows and taxa as columns. All samples must be rarefied to the same depth.
#If stats=TRUE the function will return fitting statistics.
#If stats=FALSE the function will return a table of observed and predicted values for each otu.

sncm.fit<- function(comun, stats=TRUE){
	require(minpack.lm)
	require(Hmisc)
	require(stats4)
	
	options(warn=-1)

	#Calculate the number of individuals per community
	N <- round(mean(apply(comun, 1, sum)))
	
	#Calculate the average relative abundance of each taxa across communities
	p.m <- apply(comun, 2, mean)
	p.m <- p.m[p.m != 0]
	p <- p.m/N
	
	#Calculate the occurrence frequency of each taxa across communities
	comun.bi <- 1*(comun>0)
	freq <- apply(comun.bi, 2, mean)
	freq <- freq[freq != 0]

	#Combine
	C <- merge(p, freq, by=0)
	C <- C[order(C[,2]),]
	C <- as.data.frame(C)
	C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
	p <- C.0[,2]
	freq <- C.0[,3]
	names(p) <- C.0[,1]
	names(freq) <- C.0[,1]

	#Calculate the limit of detection
	d = 1/N

	##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
	m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
	m.ci <- confint(m.fit, 'm', level=0.95)
	
	##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
	sncm.LL <- function(m, sigma){
		R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
		R = dnorm(R, 0, sigma)
		-sum(log(R))
	}
	m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
	
	##Calculate Akaike's Information Criterion (AIC)
	aic.fit <- AIC(m.mle, k=2)
	bic.fit <- BIC(m.mle)

	##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
	freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
	Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
	
	pred.ci <- binconf(freq.pred*nrow(comun), nrow(comun), alpha=0.05, method="wilson", return.df=TRUE)
	
	##Calculate AIC for binomial model
	bino.LL <- function(mu, sigma){
		R = freq - pbinom(d, N, p, lower.tail=FALSE)
		R = dnorm(R, mu, sigma)
		-sum(log(R))
	}
	bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
	
	aic.bino <- AIC(bino.mle, k=2)
	bic.bino <- BIC(bino.mle)
	
	##Goodness of fit for binomial model
	bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
	Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))

	bino.pred.ci <- binconf(bino.pred*nrow(comun), nrow(comun), alpha=0.05, method="wilson", return.df=TRUE)
	
		##Results
	if(stats==TRUE){
		fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), RMSE=numeric(), RMSE.bino=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
		fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, Rsqr, Rsqr.bino, RMSE, RMSE.bino, aic.fit, bic.fit, aic.bino, bic.bino, N, nrow(comun), length(p), d)
		return(fitstats)
	} else {
		A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
		A <- as.data.frame(A)
		colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
		B <- A[order(A[,1]),]
		return(B)
	}
}

#对模型中具体值的计算
a=read.csv("b_otu_ev.csv",row=1,header=T)
botu_A=a[1:27,]
botu_B=a[28:54,]
botu_C=a[55:81,]

b=read.csv("f_otu_ev.csv",row=1,header=T)
fotu_A=b[1:27,]
fotu_B=b[28:54,]
fotu_C=b[55:81,]

date=botu_A
N <- round(mean(apply(data, 1, sum)))
p.m <- apply(data, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
data.bi <- 1*(data>0)
freq <- apply(data.bi, 2, mean)
freq <- freq[freq != 0]

library(minpack.lm)
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] 
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))

#sncm.fit()为其上定义的中性模型函数
#comun为其上的spp,且为均一化(even)的绝对丰度（absolute）矩阵
b_A<-sncm.fit(botu_A,stats=TRUE)
b_A

a<-sncm.fit(fotu_A,stats=F)#stats=F
a

a1<-a[a$freq>a$pred.upr,]
a2<-a[a$pred.upr>a$freq,]
a2<-a2[a2$freq>a2$pred.lwr,]
a3<-a[a$pred.lwr>a$freq,]

dim(a1)
dim(a2)
dim(a3)

#对不同OTUs的定义
a[which(a$freq>a$pred.upr),'Type']="over"
a[which(a$pred.lwr>a$freq),'Type']="under"
a[which(a$pred.upr>a$freq & a$freq>a$pred.lwr,),'Type']="neutral"

#可视化
pdf(file="b_A.pdf")
plot(log(p),freq,cex.lab=1.5,cex.axis=1.5,col='white',xlab='log(Mean relative abundance)',ylab='Occurrence frequency')
points(log(a1$p),a1$freq,col='coral',pch=19,cex =0.5)
points(log(a2$p),a2$freq,col='gray50',pch=19,cex =0.5)
points(log(a3$p),a3$freq,col='purple',pch=19,cex =0.5)
lines(log(p),fitted(m.fit),col='red',lwd=3)
lines(log(p),a$pred.upr,lty=3,col="blue",lwd=2.5)
lines(log(p),a$pred.lwr,lty=3,col="blue",lwd=2.5)
text(-13,1,"R2=0.87)
text(-13,0.95,"m=0.71")
dev.off()

#model parametes 
b_A<-sncm.fit(botu_A,stats=TRUE)
b_B<-sncm.fit(botu_B,stats=TRUE)
b_C<-sncm.fit(botu_C,stats=TRUE)
f_A<-sncm.fit(fotu_A,stats=TRUE)
f_B<-sncm.fit(fotu_B,stats=TRUE)
f_C<-sncm.fit(fotu_C,stats=TRUE)

a=rbind(b_A,b_B,b_C,f_A,f_B,f_C)
rownames(a)=c("b_A","b_B","b_C","f_A","f_B","f_C")
a
write.csv(a,"model_parameters.csv")


##刘尧公号
library(Hmisc)
library(minpack.lm)
library(stats4)
 

spp=fotu_A
 
##将 Sloan 等（2006）的中性模型拟合到一个物种或分类群的丰度表，并返回几个拟合统计数据。或者，将根据它们在元群落中的丰度返回每个分类群的预测出现频率
#用非线性最小二乘法（Non-linear least squares，NLS）拟合模型参数
N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  #获取 m 值
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  #获取模型的 R2
 
#输出 3 个统计结果数据表，包括各物种或分类群的平均相对丰度（p.csv）、出现频率（freq.csv）和预测的出现频率（freq.pred.csv）
write.csv(p, file = "p.csv")
write.csv(freq, file = "freq.csv")
write.csv(freq.pred, file = "freq.pred.csv")
 
#p 是平均相对丰度（mean relative abundance）
#freq 是出现频率（occurrence frequency）的观测值
#freq.pred 是出现频率（occurrence frequency）的预测值，即中性模型的拟合值
 
points(log(a1$p),a1$freq,col='coral',pch=19,cex =0.5)
points(log(a2$p),a2$freq,col='gray50',pch=19,cex =0.5)
points(log(a3$p),a3$freq,col='purple',pch=19,cex =0.5)

#绘制统计图
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('gray50',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'purple'#出现频率低于中性群落模型预测的部分
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'coral'#出现频率高于中性群落模型预测的部分
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.5))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='red',lwd=2),default='native')
 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 
#grid.text(x=unit(0,'npc')-unit(-1,'lines'), y=unit(0,'npc')-unit(-15,'lines'),label='Mean Relative Abundance (log)', gp=gpar(fontface=2)) 
#grid.text(round(coef(m.fit)*N),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2)) 
#grid.text(label = "Nm=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2))
#grid.text(round(Rsqr,2),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
#grid.text(label = "Rsqr=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)





