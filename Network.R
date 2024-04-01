####数据筛选####
botu=read.csv("b_otu_re.csv",row=1,header=T)
fotu=read.csv("f_otu_re.csv",row=1,header=T)

b_A=botu[1:27,]
b_B=botu[28:54,]
b_C=botu[55:81,]

f_A=fotu[1:27,]
f_B=fotu[28:54,]
f_C=fotu[55:81,]


#根据OTU的相对丰度，建议筛选按丰度（0.01%）筛选
otu=f_A
otu <- otu[,colSums(otu)/nrow(otu)>=(0.01/100)]
genus1 <- t(otu)
genus1[genus1>0] <- 1
genus <- genus1[which(rowSums(genus1) >= 20), ] #共计20个样品，只保留在 10 个及以上样本中出现的OTU（MENA分析中有过滤步骤，建议筛选好后倒出作为构建网络的原始数据）
otu=otu[,rownames(genus)]
f_A=otu
dim(b_A)

otu_A=cbind(b_A,f_A)
otu_B=cbind(b_B,f_B)
otu_C=cbind(b_C,f_C)

#注释信息（构建网络所选节点的属性，例如：所属分类，指示物种，功能特性等）
gname<-read.csv("gname.csv",row.name=1,header=T)

####相关性计算####，获得邻接矩阵（R中选用两种方法计算，correlation and SpiecEasi，领接矩阵的名称进行标注与后期网络领接列表相对应）
###correlation 
library(Hmisc)
#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
corr <- rcorr(as.matrix(otu_C), type = 'spearman')
#阈值筛选
#将 spearman 相关系数低于 0.6 的关系剔除，即 r>=0.6
r <- corr$r
r[abs(r) < 0.8] <- 0

#选取显著性 p 值小于 0.01 的相关系数，即 p<0.01
p <- corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.01] <- -1
p[p<0.01 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'genus_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

####网络领接列表的获得####(对不同网络进行特别标注后，保存空间后期加载igraph可直接调取)
library(igraph)
#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
g

###网络图的精简（可选）
g<- simplify(g)#去掉自相关连线
is.simple(g)#去掉重复的连接
g<-delete.vertices(g,names(degree(g)[degree(g)==0]))##删除度为0的点（孤立点的删除）

###对g的节点的标注
rel<-colSums(otu_A)/nrow(otu_A)##注意对相应数据的替换
V(g)$rel<-rel#相对丰度
V(g)$phylum<-as.vector(gname[V(g)$name,2])
V(g)$type<-as.vector(gname[V(g)$name,1])

###对边正负相关的注释
E(g)$correlation <- E(g)$weight 
E(g)$type<- as.vector(E(g)$correlation>0)#positive
#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$weight <- abs(E(g)$weight)

###网络的输出进行可视化
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(g, 'C.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(g, 'network.gml', format = 'gml')

#输出边列表文件，用于NetShift分析，识别网络中的驱动节点
g=bs_g
edge <- data.frame(as_edgelist(g))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
    source = edge[[1]],
    target = edge[[2]],
    weight = E(g)$weight,
    correlation = E(g)$correlation
)
head(edge_list)

write.table(edge_list, 'network.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)



####网络结构的进一步分析和比较####（网络属性，节点属性，随机网络的比较，网络的度分布模式，网络节点的zp属性分析）

###网络整体属性的计算
# net_pro:自定义函数，提供需要计算网络性质的igraph对象，结果返回计算好的网络性质
igraph=C_g
igraph.weight=E(igraph)$correlation 
source("net_pro.R")
netpro_result<-net_pro(igraph)
write.csv(netpro_result,"C_network.pro.csv")

###网络节点属性的计算
# node_pro:自定义函数，提供需要计算网络性质的igraph对象，结果返回计算好的节点性质
igraph=C_g
source("node_pro.R")
nodepro_result<-node_pro(igraph)
write.csv(nodepro_result,"C_node.pro.csv")

###随机网络的属性比较
igraph=g
rand.g.netpro_result<-c()
for (i in 1:10){
  #random null model
  rand.g<- erdos.renyi.game(length(V(igraph)), length(E(igraph)),type = c("gnm"))
  tem_netpro_result<-net_pro(rand.g)
  rand.g.netpro_result<-cbind(rand.g.netpro_result,tem_netpro_result)
} 
write.csv(rand.g.netpro_result,"results/rand.g.1000.result.csv")

# 对随机矩阵结果求均值和sd值
result_summary<-cbind(rowMeans(rand.g.netpro_result),apply(rand.g.netpro_result,1,sd))
colnames(result_summary)<-c("Means","SD")
write.csv(result_summary,"results/wild.1000.result.summary.csv")


###网络度分布模式的计算
#计算节点度
igraph=bs_cult_g
V(igraph)$degree <- degree(igraph)
head(V(igraph)$degree)

#度分布统计
degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(names(degree_dist))
degree_count <- as.numeric(degree_dist)
dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)

#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(V(igraph)$degree, xlab = 'Degree', ylab = 'Frequency', 
    main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
    main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
    ylab = 'Log-count', main = 'Log-log degree distribution')

# 展示大随机网络度分布图
g_random <- erdos.renyi.game(1000, 10000,type = c("gnm"))
degree_distribution(g_random)
plot(degree_distribution(g_random, cumulative = FALSE),xlab="Degree",main="The distribution of degree for co-occurrence network")ggsave("Protist.pdf", p, width = 8, height =6)

##幂律分布
#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a =4, b = -1))
summary(mod)

##log 转化后线性关系的拟合，log函数向幂律函数的转化？
lm(log(count)~log(degree+1),data=dat)
log(count)~-1.010 *log(degree+1)+4.096

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b

#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2

#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
    dat_rand$count <- sample(dat_rand$count)
    SSre_rand <- sum((dat_rand$count-fit)^2)
    SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
    R2_rand <- 1 - SSre_rand/SStot_rand
    if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value

#作图展示
library(ggplot2)

p <- ggplot(dat, aes(x = degree, y = count)) +
geom_point(color = 'black') +
theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
stat_smooth(method = 'nls', formula = y ~ a*x^b,color="red",se = F) +
labs(x = 'Degree', y = 'Count')

#添加公式拟合的注释
label <- data.frame(
    formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
    R2 = sprintf('italic(R^2) == %.3f', R2),
    p_value = sprintf('italic(P) < %.3f', p_value)
)

p=p + geom_text(x = 60, y = 100, aes(label = formula), data = label, parse = TRUE, hjust = 0) +#注意根据实际坐标情况调整公式的位置
geom_text(x = 60, y = 90, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
geom_text(x = 60, y = 70, aes(label = p_value), data = label, parse = TRUE, hjust = 0)
p
ggsave("rs.pdf", p, width = 5, height =4.5)

##多个网络的分面图（且节点数可以不同）
#数据
dat2 <- read.csv('degree2.csv', stringsAsFactors = FALSE)#多个网络的度分布统计合并且分类

#构建计算函数
lm_labels <- function(dat) {
    #拟合
    mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
    a <- round(coef(mod)[1], 3)
    b <- round(coef(mod)[2], 3)
    
    #计算R2
    fit <- fitted(mod)
    SSre <- sum((dat$count-fit)^2)
    SStot <- sum((dat$count-mean(dat$count))^2)
    R2 <- round(1 - SSre/SStot, 3)
    
    #计算 p 值（置换检验原理）
    p_num <- 1
    dat_rand <- dat
    for (i in 1:999) {
        dat_rand$count <- sample(dat_rand$count)
        SSre_rand <- sum((dat_rand$count-fit)^2)
        SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
        R2_rand <- 1 - SSre_rand/SStot_rand
        if (R2_rand > R2) p_num <- p_num + 1
    }
    p_value <- p_num / (999+1)
    
    #方程式值的列表
    data.frame(formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
        R2 = sprintf('italic(R^2) == %.3f', R2),
        p_value = sprintf('italic(P) < %.3f', p_value))
}

#计算及作图
library(plyr)
library(ggplot2)

label <- ddply(dat2, 'Type', lm_labels)

p <- ggplot(dat2, aes(x = degree, y = count, color=Type)) +
geom_point() +
facet_wrap(~Type, nrow = 1) +
stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE, show.legend = FALSE) +
labs(x = 'Degree', y = 'Count')

p + geom_text(x = 40, y = 180, aes(label = formula), data = label, parse = TRUE, hjust = 0, color = 'black', show.legend = FALSE) +
geom_text(x = 40, y = 160, aes(label = R2), data = label, parse = TRUE, hjust = 0, color = 'black', show.legend = FALSE) +
geom_text(x = 40, y = 140, aes(label = p_value), data = label, parse = TRUE, hjust = 0, color = 'black', show.legend = FALSE)

###网络的模块化计算(Module)
igraph=g
# 模块性 modularity 注意模块化计算时对边的权重的转化，尤其是存在负的相关性时，weights =NULL 意为权重都是相等的
fc <- cluster_fast_greedy(igraph)# 聚类方法的选择cluster_walktrap,cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
modularity <- modularity(igraph,membership(fc))
modularity
comps <- membership(fc)# 每个OTU所属模块
max(comps)#模块数量

##按照某个模块构建子图
modules.g<-induced_subgraph(igraph,comps==1)
write.graph(modules.g,"2.graphml", format="graphml")
E(modules.g)$weight<-NA
plot(modules.g,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.width=1,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))


###网络节点z-p值的计算
igraph=C_g
#计算节点度
V(igraph)$degree <- degree(igraph)

#模块划分，详情 ?cluster_fast_greedy，有多种模型
set.seed(123)
E(igraph)$weight <- abs(E(igraph)$weight)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))#注意含有正负边的情况的模块的计算weight=null;或者E(igraph)$weight <- abs(E(igraph)$weight) 推荐用绝对值化的权重
#推荐用绝对值化的权重，两者的聚类结果有差异，因此对后面zp值的计算必有影响

#输出各节点（微生物 OTU）名称、节点度、及其所划分的模块的列表
nodes_list <- data.frame(
    nodes_id = V(igraph)$name, 
	degree = V(igraph)$degree, 
	modularity = V(igraph)$modularity
)
head(nodes_list)    #节点列表，包含节点名称、节点度、及其所划分的模块

write.table(nodes_list, 'nodes_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

##计算模块内连通度（Zi）和模块间连通度（Pi）
source('zi_pi.r')

#igraph 转邻接矩阵类型
adjacency_unweight=as.matrix(get.adjacency(igraph)) 

#节点属性列表，包含节点所划分的模块（注意此时第一列为行名）
nodes_list <- read.delim('nodes_list.txt', row.names = 1, sep = '\t', check.names = FALSE)

#两个文件的节点顺序要一致
nodes_list <- nodes_list[rownames(adjacency_unweight), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)


##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'
zi_pi
write.csv(zi_pi, 'C_zi_pi_result.csv')

library(ggplot2)
p=ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
geom_point(aes(color = type), alpha = 0.5, size = 2) +
scale_color_manual(values = c('gray','red','blue','purple'), 
    limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
    panel.background = element_blank(), legend.key = element_blank()) +
labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
geom_vline(xintercept = 0.62) +geom_text(aes(label=nodes_id), size=2)+
geom_hline(yintercept = 2.5)
p
ggsave("C1_zp.pdf", p, width = 5, height =4.5)


####网络抗毁性（鲁棒性/稳定性的检验）####
##微生物共发生网络的自然连通度(natural_connectivity)计算
#计算自然连通度
#adj_matrix 是网络邻接矩阵
nc <- function(adj_matrix) {
    #获取 0-1 矩阵，1 表示节点间存在边，0 表示不存在边
    adj_matrix <- as.matrix(adj_matrix)
    adj_matrix[abs(adj_matrix) != 0] <- 1
    
    #矩阵的特征分解，获取特征值 λ
    lambda <- eigen(adj_matrix, only.values = TRUE)$values
    lambda <- sort(lambda, decreasing = TRUE)
    
    #计算“平均特征根”，获得自然连通度
    lambda_sum <- 0
    N = length(lambda)
    for (i in 1:N) lambda_sum = lambda_sum + exp(lambda[i])
    lambda_average <- log(lambda_sum/N, base = exp(1))
    lambda_average
}


#计算自然连通度
igraph=A_g
#igraph 转邻接矩阵类型
adjacency_unweight=as.matrix(get.adjacency(igraph)) 
adj_matrix=adjacency_unweight#对应网络邻接列表的邻接矩阵
natural_connectivity <- nc(adj_matrix)
natural_connectivity

#模拟随机移除节点，就是随机在邻接矩阵中移除行和列中对应的 OTU
set.seed(123)
for (i in 1:300) {#移除的节点数
    
    #在邻接矩阵中随机移除 i 个节点
    remove_node <- sample(1:543, i)#总的节点数
    adj_matrix2 <- adj_matrix[-remove_node,-remove_node]
    
    #计算自然连通度
    natural_connectivity_remove <- nc(adj_matrix2)
    natural_connectivity <- c(natural_connectivity, natural_connectivity_remove)
}
dat <- data.frame(remove_node = 0:(length(natural_connectivity)-1), natural_connectivity = natural_connectivity)
write.csv(dat, 'dat_A1.csv')

#ggplot2 散点图+拟合线，展示随机移除节点后网络的自然连通度
library(ggplot2)
dat=read.csv("dat_A1.csv",row=1,header=T)
ggplot(dat, aes(remove_node, natural_connectivity,color=Type))+
geom_point() +theme_bw()+geom_smooth(se = FALSE)

p=ggplot(dat, aes(remove_node, natural_connectivity,color=Type)) +theme_bw()+
geom_line()+xlab("Number of remove node")+ylab("Natural connectivity")
p
ggsave("network_stability.pdf", p, width = 7, height =4.5)

###和给定微生物网络具有相同节点和边数量的随机网络的自然连通度计算
library(igraph)

#以广义随机图为例，首先需要定义图的集合
degree_dist <- table(degree(g))
degree_num <- as.numeric(names(degree_dist))
degree_count <- as.numeric(degree_dist)
names(degree_count) <- degree_num
degs <- rep(degree_num, degree_count)

#获得广义随机图（构建 3 个），并转换为邻接矩阵
set.seed(123)
g_rand1 <- degree.sequence.game(degs, method = 'simple')
adj_matrix_rand1 <- as.matrix(get.adjacency(g_rand1))
g_rand2 <- degree.sequence.game(degs, method = 'simple')
adj_matrix_rand2 <- as.matrix(get.adjacency(g_rand2))
g_rand3 <- degree.sequence.game(degs, method = 'simple')
adj_matrix_rand3 <- as.matrix(get.adjacency(g_rand3))

#随机移除节点，并计算自然连通度
natural_connectivity_rand1 <- nc(adj_matrix_rand1)
natural_connectivity_rand2 <- nc(adj_matrix_rand2)
natural_connectivity_rand3 <- nc(adj_matrix_rand3)

for (i in 1:150) {
    
    #在邻接矩阵中随机移除 i 个节点
    remove_node <- sample(1:190, i)
    adj_matrix_rand1_remove <- adj_matrix_rand1[-remove_node,-remove_node]
    adj_matrix_rand2_remove <- adj_matrix_rand2[-remove_node,-remove_node]
    adj_matrix_rand3_remove <- adj_matrix_rand3[-remove_node,-remove_node]
    
    #计算自然连通度
    natural_connectivity_rand1 <- c(natural_connectivity_rand1, nc(adj_matrix_rand1_remove))
    natural_connectivity_rand2 <- c(natural_connectivity_rand2, nc(adj_matrix_rand2_remove))
    natural_connectivity_rand3 <- c(natural_connectivity_rand3, nc(adj_matrix_rand3_remove))
}

#ggplot2 作图，微生物网络和随机网络一起，拟合线
library(ggplot2)

dat <- data.frame(remove_node = rep(0:150, 4), 
    natural_connectivity = c(natural_connectivity, natural_connectivity_rand1, natural_connectivity_rand2, natural_connectivity_rand3), 
    network = c(rep('microbial network', 151), rep('random network 1', 151), rep('random network 2', 151), rep('random network 3', 151)))
#write.csv(dat, 'dat.csv', row.names = FALSE, quote = FALSE)

ggplot(dat, aes(remove_node, natural_connectivity, color = network)) +
geom_smooth(se = FALSE)+theme_bw()

##3种或多种微生物网络鲁棒性评估
library(igraph)
library(ggplot2)

#读取网络邻接矩阵，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
adj1 <- read.csv('s.adj.csv', row.names = 1, check.names = FALSE)
adj2 <- read.csv('r.adj.csv', row.names = 1, check.names = FALSE)
adj3 <- read.csv('e.adj.csv', row.names = 1, check.names = FALSE)

#计算自然连通度
natural_connectivity1 <- nc(adj1)
natural_connectivity2 <- nc(adj2)
natural_connectivity3 <- nc(adj3)

#转化为 igraph 邻接列表，计算节点平均度
g1 <- graph_from_adjacency_matrix(as.matrix(adj1), mode = 'undirected', diag = FALSE)
g2 <- graph_from_adjacency_matrix(as.matrix(adj2), mode = 'undirected', diag = FALSE)
g3 <- graph_from_adjacency_matrix(as.matrix(adj3), mode = 'undirected', diag = FALSE)
average_degree1 <- mean(degree(g1))
average_degree2 <- mean(degree(g2))
average_degree3 <- mean(degree(g3))

#随机移除 200 个节点，并计算上述 3 种网络特征
for (i in 1:200) {
    
    #在邻接矩阵中随机移除 i 个节点
    remove_node <- sample(1:nrow(adj1), i)
    adj1_remove <- adj1[-remove_node,-remove_node]
    remove_node <- sample(1:nrow(adj2), i)
    adj2_remove <- adj2[-remove_node,-remove_node]
    remove_node <- sample(1:nrow(adj3), i)
    adj3_remove <- adj3[-remove_node,-remove_node]
    
    #计算自然连通度
    natural_connectivity1 <- c(natural_connectivity1, nc(adj1_remove))
    natural_connectivity2 <- c(natural_connectivity2, nc(adj2_remove))
    natural_connectivity3 <- c(natural_connectivity3, nc(adj3_remove))
    
    #计算节点平均度
    g1 <- graph_from_adjacency_matrix(as.matrix(adj1_remove), mode = 'undirected', diag = FALSE)
    g2 <- graph_from_adjacency_matrix(as.matrix(adj2_remove), mode = 'undirected', diag = FALSE)
    g3 <- graph_from_adjacency_matrix(as.matrix(adj3_remove), mode = 'undirected', diag = FALSE)
    average_degree1 <- c(average_degree1, mean(degree(g1)))
    average_degree2 <- c(average_degree2, mean(degree(g2)))
    average_degree3 <- c(average_degree3, mean(degree(g3)))
}

#ggplot2 作图，拟合线
dat <- data.frame(remove_node = rep(0:200, 6), 
    variable = c(rep(c(rep('natural_connectivity', 201), rep('average_degree', 201)), 3)), 
    values = c(natural_connectivity1, average_degree1, natural_connectivity2, average_degree2, natural_connectivity3, average_degree3), 
    network = c(rep('s', 201*2), rep('r', 201*2), rep('e', 201*2)))
#write.csv(dat, 'dat.csv', row.names = FALSE, quote = FALSE)

ggplot(dat, aes(remove_node, values, color = network)) +
geom_smooth(se = FALSE) +
facet_wrap(~variable, ncol = 2, scale = 'free_y')

#网络之间的比较
##不同网络图之间节点和边的简单比较
#node
nodenum_bs=names(V(bs_g))
nodenum_rs=names(V(rs_g))
nodenum_bs_c=names(V(bs_cult_g))
nodenum_rs_c=names(V(rs_cult_g))
nodenum_bs_w=names(V(bs_wild_g))
nodenum_rs_w=names(V(rs_wild_g))

#edges
edgenum_bs<-length(E(bs_g))
edgenum_rs<-length(E(rs_g))
edgenum_bs_c<-length(E(bs_cult_g))
edgenum_rs_c<-length(E(rs_cult_g))
edgenum_bs_w<-length(E(bs_wild_g))
edgenum_rs_w<-length(E(rs_wild_g))


# the numbers of shared edges of WT and OE network
nodenum_bs_rs<-intersect(nodenum_rs_w,nodenum_bs_w)
edgenum_bs_rs<-length(E(intersection(bs_g,rs_g)))

a=Reduce(intersect, list(nodenum_bs_c,nodenum_rs_c,nodenum_bs_w,nodenum_rs_w)) 

# 韦恩图展示特有边和共有边情况
library(VennDiagram)
#node
nodenum_bs=length(nodenum_bs)
nodenum_rs=length(nodenum_rs)
nodenum_bs_rs=length(nodenum_bs_rs)
pdf(file = paste0("nodenum_bs_rs_venn.pdf"),width = 13)
grid.newpage()
draw.pairwise.venn(area1=nodenum_bs,area2=nodenum_rs,cross.area=nodenum_bs_rs
                   ,category=c('BS','RS'),lwd=rep(1,1)
                   ,col=c('red','green'),fill=c('red','green')
                   ,cat.col=c('red','green')
                   ,rotation.degree=90)
dev.off()

#edge
pdf(file = paste0("edgenum_bs_rs_venn.pdf"),width = 13)
grid.newpage()
draw.pairwise.venn(area1=edgenum_bs,area2=edgenum_rs,cross.area=edgenum_bs_rs
                   ,category=c('BS','RS'),lwd=rep(1,1)
                   ,col=c('red','green'),fill=c('red','green')
                   ,cat.col=c('red','green')
                   ,rotation.degree=90)
dev.off()





