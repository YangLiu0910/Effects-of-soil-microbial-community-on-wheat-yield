##随机森林分类模型 判断解释变量对不同分组（类别型响应变量）的贡献度
library(randomForest)
set.seed(315)
rf = randomForest(soil_c[,4:14], soil_c$Treatment, importance=TRUE, proximity=TRUE, ntree = 1000)
print(rf)
summary(rf)
rf$predicted
imp= as.data.frame(rf$importance)
write.csv(imp,"soil_c_importance.csv")

#barplot
library(ggplot2)
imp=read.csv("soil_a_importance.csv",row=1,header=T)
imp=imp[order(imp[,6],decreasing = T),]
imp
p=ggplot(data=imp, mapping = aes(x=reorder(names,MeanDecreaseGini),y=MeanDecreaseGini,fill=Group)) + 
  geom_bar(stat="identity")+coord_flip()+theme_bw()
p
