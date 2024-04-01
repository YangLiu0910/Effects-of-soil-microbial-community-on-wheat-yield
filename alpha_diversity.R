library(vegan)
otu=read.csv("f_otu_ev.csv",row=1,header=T)
tree <- read.tree('otu_tree.tre')

alpha <- function(x, tree = NULL, base = exp(1)) {
       est <- estimateR(x)
       Richness <- est[1, ]
       Chao1 <- est[2, ]
       ACE <- est[4, ]
       Shannon <- diversity(x, index = 'shannon', base = base)
       Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
       Pielou <- Shannon / log(Richness, base)
       goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
       
       result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
       if (!is.null(tree)) {
              PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
              names(PD_whole_tree) <- 'PD_whole_tree'
              result <- cbind(result, PD_whole_tree)
       }
       result
}

#不包含谱系多样性，无需指定进化树；Shannon 公式的 log 底数我们使用 2
alpha_all <- alpha(otu, base = 2)##otu (samples as rows) and otu_even 
#包含谱系多样性时，指定进化树文件；Shannon 公式的 log 底数我们使用 2
alpha_all <- alpha(otu, tree, base = 2)
write.csv(alpha_all, 'f_alpha.csv')

#分组求均值和标准差
balpha=read.csv("b_alpha.csv",row=1,header=T)
falpha=read.csv("f_alpha.csv",row=1,header=T)


b_mean<-aggregate(balpha[,1:2],by=list(fa[,4]),FUN=mean)
b_sd<-aggregate(balpha[,1:2],by=list(fa[,4]),FUN=sd)

f_mean<-aggregate(falpha[,1:2],by=list(fa[,4]),FUN=mean)
f_sd<-aggregate(falpha[,1:2],by=list(fa[,4]),FUN=sd)


mean=cbind(f_mean,b_mean)
sd=cbind(f_sd,b_sd)
a=rbind(mean,sd)
write.csv(a,"alpha_mean+sd.csv")

##多样性的比较
library(ggplot2)
library(magrittr)
library(ggpubr)
df=talpha
p<-ggplot(df,aes(x=group,y=Shannon,fill="red"))+geom_bar(stat="identity",width=0.95)+theme_bw()
p+facet_grid(~type)

#boxplot
df$type<-factor(df$type,levels=c("P","F","B"))
p <- ggboxplot(df, x = "group", y = "Richness", fill = "Type", 
       facet.by="type",size=0.1,alpha=0.8)


+facet_grid(Split~Compartment)
p

#Kruskal wallis test(非参数)
library(agricolae)
group=fa[,4]
comparison<-with(falpha[,1:2],kruskal(Richness,group,group=TRUE, main="corn"))
comparison

comparison<-with(falpha[,1:2],kruskal(Shannon,group,group=TRUE, main="corn"))
comparison

#Two-way MANOVA 
Compartment=fa[,1]
Type=fa[,2]
fit <- manova(as.matrix(balpha[,1:2]) ~Compartment*Type,data=fa)
summary(fit)
summary.aov(fit)#与ANOVA的结果是一样的







