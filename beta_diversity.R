##加载程序包
library(vegan)#NMDS+多元方差分析
library(ggplot2)
library(ggrepel)
library(devtools)
library(ape)#PCoA

##距离计算和提取
botu=read.csv("b_otu_re.csv",row=1,header=T)
fotu=read.csv("f_otu_re.csv",row=1,header=T)

#Treatment group 
botu_ck=botu[1:27,]
botu_sm=botu[28:54,]
botu_pm=botu[55:81,]

fotu_ck=fotu[1:27,]
fotu_sm=fotu[55:81,]
fotu_pm=fotu[28:54,]


b_bray<-vegdist(botu)
f_bray<-vegdist(fotu)

b_bray_ck<-vegdist(botu_ck)
b_bray_sm<-vegdist(botu_sm)
b_bray_pm<-vegdist(botu_pm)

f_bray_ck<-vegdist(fotu_ck)
f_bray_sm<-vegdist(fotu_sm)
f_bray_pm<-vegdist(fotu_pm)

#Aggregate
botu_a=botu[1:27,]
botu_b=botu[28:54,]
botu_c=botu[55:81,]

fotu_a=fotu[1:27,]
fotu_b=fotu[28:54,]
fotu_c=fotu[55:81,]


b_bray_a<-vegdist(botu_a)
b_bray_b<-vegdist(botu_b)
b_bray_c<-vegdist(botu_c)

f_bray_a<-vegdist(fotu_a)
f_bray_b<-vegdist(fotu_b)
f_bray_c<-vegdist(fotu_c)


##样方分组
fa=read.csv("group.csv",row=1,header=T)
Treatment=fa[,1]
Layer=fa[,2]
Aggregate=fa[,3]

Treatment=fa[1:27,1]
Layer=fa[1:27,2]



fa$Treatment<-factor(fa$Treatment,levels=c("CK","SM","PM"))

#PCoA
res<- pcoa(as.dist(b_bray_b))##模型计算和选择
res$values$Relative_eig
si<-res$vectors[,1:2]##提取样方坐标

##可视化
PCoA1=round(res$values$Relative_eig[1],4)
PCoA2=round(res$values$Relative_eig[2],4)
si<-as.data.frame(si)
p<-ggplot(si,aes(x=Axis.1,y=Axis.2,color=Treatment,shape=Layer))
p_f_b=p+geom_point(size=3.5,alpha=0.65)+xlab(paste("PCoA1=",PCoA1*100,"%",sep=""))+ylab(paste("PCoA2=",PCoA2*100,"%",sep="")) + labs(title = "Bray-Curtis PCoA analysis")+theme(text=element_text(size=12),panel.grid=element_blank(),
panel.background=element_rect(fill='transparent', color='black')
,panel.border=element_rect(fill='transparent', color='transparent'))
p_f_b

#NMDS
m<-metaMDS(f_bray_a)##模型计算和选择
m$stress
si<-m$points[,1:2]#NMDS##提取样方坐标

si<-as.data.frame(si)
p<-ggplot(si,aes(x=MDS1,y=MDS2,color=Treatment,shape=Layer))
p_f_a=p+geom_point(size=4.5,alpha= 0.75)+xlab("MDS1")+ylab("MDS2")+ labs(title = paste("NMDS\n","Stress=",round(m$stress,3),sep="\t"))+theme(text=element_text(size=12),panel.grid=element_blank(),
panel.background=element_rect(fill='transparent', color='black')
,panel.border=element_rect(fill='transparent', color='transparent'))
p_f_a


#+stat_ellipse(linetype = 2,lwd=1)

library(patchwork)
bacteria=p_b_ck+p_b_sm+p_b_pm+plot_layout(guides = 'collect')
ggsave("bacterial_pcoa.pdf", bacteria, width = 12, height = 4)

bacteria=p_b_a+p_b_b+p_b_c+plot_layout(guides = 'collect')
ggsave("bacterial__arregate_nmds.pdf", bacteria, width = 12, height = 4)

fungi=p_f_ck+p_f_sm+p_f_pm+plot_layout(guides = 'collect')
ggsave("fungal_pcoa.pdf", fungi, width = 12, height = 4)

fungi=p_f_a+p_f_b+p_f_c+plot_layout(guides = 'collect')
ggsave("fungal_aggregate_nmds.pdf", fungi, width = 12, height = 4)

##adonis 全局检验
adonis_result=adonis(f_bray_b~Treatment*Layer)
summary(adonis_result)
#结果的输出
adonis_result<- data.frame(adonis_result$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
adonis_result
write.csv(adonis_result, file = 'f_b_adonis.csv')

##距离计算和主成分提取
bray<-vegdist(botu)
#PCoA
res<- pcoa(as.dist(bray))##模型计算和选择
res$values$Relative_eig
si<-res$vectors[,1:2]##提取样方坐标
b_pcoa=si
beta_pcoa=cbind(p_pcoa,f_pcoa,b_pcoa)
colnames(beta_pcoa)=c("p1","p2","f1","f2","b1","b2")
write.csv(beta_pcoa,"beta_pcoa.csv")





