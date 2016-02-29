raw_data <- read_delim('./data/GO_pathway.txt', delim = '\t')
head(raw_data)
library(dplyr)
library(Cairo)

double_y_axis <- function(p1, p2){
  g1 <- ggplot_gtable(ggplot_build(p1))
  g2 <- ggplot_gtable(ggplot_build(p2))

  # overlap the panel of 2nd plot on that of 1st plot
  pp <- c(subset(g1$layout, name == "panel", se = t:r))
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)

  # axis tweaks
  ia <- which(g2$layout$name == "axis-l")
  ga <- g2$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)
  ax$grobs <- rev(ax$grobs)
  ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
  g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
  g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

  # draw it
  grid.draw(g)
  return(g)
}


split_fun <- function(x){
  tmp_split <- strsplit(x, '~')
  do.call('rbind', tmp_split)[,2]
}
go_data <- raw_data %>%
  dplyr::filter(grepl('GO',Term) & PValue<0.01) %>%
  select(c(1,2,3,5)) %>%
  dplyr::mutate(term = split_fun(Term)) %>%
  arrange(Category, PValue)
go_data$term <- factor(go_data$term, levels = unique(go_data$term))

# ggplot
library(ggplot2)
library(gtable)
library(grid)

grid.newpage()

# two plots
p1 <- ggplot(go_data)+
  geom_bar(aes(x = term, y=Count,  fill=Category), stat = 'identity') +
#   geom_line(aes(term, y=PValue, group=1), stat = 'identity', colour="#B00A13") +
#   geom_point(aes(term, y=PValue, group=1),colour="red", size=3, alpha=0.4)+
  theme_bw() +
  scale_fill_manual(values = c('#6f359d', '#F28C42'),
                    labels=c('Biological Process', 'Cellular Component'))+
  theme(axis.text.x = element_text(angle = 40, hjust=1, vjust=1),
        axis.ticks.x = element_blank(),
        legend.position='left'
       # axis.title.y=element_text(vjust=8)
        )+
  scale_y_continuous(expand = c(0,0), limits=c(0, 30))+
  xlab('Go Items')+
  ylab('The number of genes')

  # coord_fixed(ratio = 0.6)# 控制长宽比
p1
 p1 <- p1+coord_fixed(ratio = 0.6)# 控制长宽比


p2 <- ggplot(go_data, aes(term, PValue, group=1)) +
  geom_line( stat = 'identity', colour="#B00A13") +
  geom_point(colour="red", size=3, alpha=0.4)+
  scale_y_continuous(limits = c(0, 0.02))+
 # scale_y_reverse()+

  scale_x_discrete(breaks=NULL)+
  theme(axis.text.x = element_text(),
        axis.ticks.x= element_blank(),
  # panel.background = element_rect(fill = NA),
   # panel.grid = element_blank(),
  legend.position='top')+
  ylab('')+
  xlab('')
p2

  # extract gtable
g <- double_y_axis(p1,p2)



tiff("Plot600.tiff", type="cairo", width = 14, height = 8, units = 'in', res = 300)
grid.draw(g)
dev.off()


# 在右下角的框里调试好图片, 然后采用
# par('din') 获取宽高
CairoPNG('plotCairo.png', width = 13.58, height = 8.2, units='in', dpi=700)

grid.draw(g)
dev.off()

