save(meth_0.01, file = './data/meth_0.01.rda')
library(ggplot2)
meth_0.01$color = as.factor(meth_0.01$logFC>0)
#这里的logFC并不是取log后的差，而是甲基化水平直接做差
# 根据该变量，对基因进行排序
meth_order <- meth_0.01[meth_0.01$logFC %>% order(decreasing = T),]
rownames(meth_order) <- 1:34
meth_graph <- meth_order[c(1:20, 34:21),]
meth_graph$gene <- factor(meth_graph$gene, levels = meth_graph$gene)
g <- ggplot(data = meth_graph, aes(x = gene, y = logFC)) +  # colname 为字符串，可以用aes_string() 进行转换
  geom_bar(aes(fill = color), # 用到data中的数据，建立映射关系均要用到aes，fill表示填充，后面会设置具体的填充颜色，如不设置会采用系统默认，凡是表示分组都要转换为factor
           # colour="white",  # colour设置边框的颜色， width 设置
           stat='identity',      # stat=‘bin’ 表示系统自行计算频数作为y轴， stat = 'identity'表示使用自己的数据作为y轴    
           #width = 0,      # 两个bar的间隔，当x轴为连续数据时不起作用
           # binwidth = 0.3,  # 区段的长度， 计算该区段的频数
           na.rm = T        # 清除空值
  )  
g


g1 <- g +
  theme_bw() +  #设置主题， 将背景作为白色#996666
  scale_fill_manual(values = c('#f48c37', '#7030a0'), name = 'Type', labels = c('Increase with Cirrhosis to HCC','Decrease with Cirrhosis HCC')) + # 设置填充的方式，manual表示手动填充， 这里的‘1’和‘0’表示`fill=factor(flag)`中的levels， 那么为legend title， label 为legend标签
  theme(panel.grid.major = element_blank() # 通过设置主题将 主（major）次（minor）网格线去掉
        # panel.grid.minor = element_blank(),
        # legend.position=c(0.85, 0.85)
        )+      # 设置legend的位置，c(0,0)表示左下角， c(1,1)表示右上角
  # scale_x_continuous 改为 coord_cartesian 效果更好
  # scale_y_continuous(breaks = seq())+
  # scale_x_continuous(limits= c(-8, 3),   # 设置x轴， continuous表示数据为连续型，limits表示数据的范围
  # breaks = seq(-8, 3, 1.5),  # breaks 表示刻度的表示
  # expand = c(-0.05, 0))+     # 采用 stat=‘bin’会发现 左边，底边图不能对齐到轴线，采用expand可以进行设置，x轴设置了c(-0.05, 0),y轴设置了c(0, 0)这是做出的尝试， 并没有深究其原因
  scale_y_continuous(labels=percent, breaks = seq(-0.2, 0.3, 0.05)) + # 同 scale_x_continuous
  # coord_cartesian(xlim = c(-8, 3), ylim = c(0, 1500))+
  # xlab(title)+    # 设置x轴的标题
  # ylab('number of molecules')  + # y轴的标题
  # ggtitle('Distribution') + # 大标题
  theme( axis.title = element_text(size= 14),  # 设置轴标题，可以细致到 axis.title.x, axis.title.y
         axis.text = element_text(size = 10),      # 设置轴刻度
         legend.title = element_text(size = 12),   # 设置图例的标题
         legend.text = element_text(size = 12))    # 设置图例的项目
g1

ggsave('meth.png', dpi = 600, width = 27, height = 13, unit = 'cm')
