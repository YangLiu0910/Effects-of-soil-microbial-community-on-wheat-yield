#不同群落多样性相关性分析
#alpha多样性
##corrplot相关性计算和可视化
library(corrplot)
a<-read.csv("total_alpha.csv",row=1,header=T)
a_bs=a[1:20,]
a_rs=a[21:40,]

#设置cor.mtest函数
cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat<-cor.mtest(a[,4:6]) #计算P-value 
env.spearman<-cor(a[,4:6],method = "spearman") #计算相关性

corrplot(env.spearman, add = F, tl.col="grey",type = "full",method = "number", number.cex=0.9,order = "original",
         tl.pos = "lt",p.mat = p.mat, sig.level = 0.05,insig = "blank")

#beta多样性
library(ggplot2)
library(ggpubr)
library(vegan)

bbray<- vegdist(botu)
fbray<- vegdist(fotu)
pbray<- vegdist(potu)


df<-data.frame(as.vector(fbray),as.vector(bbray))
names(df)<-c("bray1","bray2")

pdf("cor_f_b.pdf")
p=ggplot(data=df, aes(x=df$bray1,y=df$bray2))+xlab("Fungi")+ylab("Bacteria")
p+geom_point(color="red",alpha=0.5,size=2.8)+stat_cor(data=df, method = "spearman")+theme_bw()+stat_smooth(method="lm",color="black",se=T)
dev.off()



active=read.csv("active.csv",row=1,header=T)
act_dis=vegdist(scale(active), "euclid")

mantel(act_dis,pbray, method="spearman")







