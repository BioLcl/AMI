library(ggplot2)

mytheme<-theme(panel.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.border = element_rect(colour="black",fill=NA)) +
  theme(axis.text.x = element_text(size=15, color = "black",angle = 90),  # 设置x轴字体大小，以下同理
        axis.text.y = element_text(size=15, color = "black"), 
        axis.title.x = element_text(size=15, color = "black"), 
        axis.title.y= element_text(size=15,  color = "black"),
        legend.title = element_text(size=15,  color = "black"))+
  # ,
  # axis.text.x.bottom  = element_text(angle = 0, hjust = 0.5),
  # legend.text = element_text(size=13)) +# 分组条形组时，position = 'dodge'
  theme(panel.grid =element_blank())  +
  theme(axis.line = element_line(size=0.6, colour = "black")) +
  theme(panel.grid.major=element_line(colour=NA), panel.border = element_blank())+
  theme(axis.line = element_line(arrow = arrow(length = unit(0.35, 'cm'))))+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(15, "pt"),
        legend.key.height = unit(15, "pt"),
        legend.key.width = unit(15, "pt"),
        legend.text = element_text(size=13)
       )

setwd("D:/OneDrive/ACS/new_group2/random/importance")
data<-read.delim("Species_importance.txt",header = T)
data<-data[1:10,]
data$MeanDecreaseAccuracy<-data$MeanDecreaseAccuracy*100
data$ID <- factor(data$ID,levels=rev(data$ID),ordered = TRUE)
library(ggplot2)
colnames(data)
#rgb(224, 56, 55,max=255)
#rgb(90, 154, 211,max=255)
gg<-ggplot(data=data, aes(x=ID,y=MeanDecreaseAccuracy))+
  geom_hline(yintercept = c(0,0.5,1,1.5,2),color="gray")+
  geom_bar(stat="identity",fill=rgb(224, 56, 55,max=255),width = 0.5)+
  scale_y_continuous(limits=c(0, 2),expand = c(0.005,0))




exp_bar<-gg + theme(panel.grid.major = element_line(colour = "gray71", 
                                                    linetype = "blank"), 
                    axis.line = element_line(linetype = "blank"), 
                    axis.ticks = element_line(linetype = "blank"), 
                    panel.background = element_rect(fill = NA, 
                                                    linetype = "solid"),
                    panel.grid.minor = element_line(colour = NA, linetype = "blank"),
                    axis.text.x = element_text(colour = "black", vjust = 0.5,
                                               hjust = 1, angle = 90,size = 15),
                    axis.text.y = element_text(colour = "black",size=15), 
                    axis.title.x = element_text(size=17, color = "black"), 
                    axis.title.y= element_text(size=15,  color = "black"),
                    legend.text = element_text(color = "black", size = 15),
                    plot.title = element_text(hjust = 0.5,size = 20)) +
  labs(x = NULL, y = "Mean Decrease Accuray",title = "Species")+
  guides(fill=guide_legend(ncol=2,title = ""))+
  coord_flip()
exp_bar
ggsave("Species_imp_col_top10.png",exp_bar,width = 7,height =3.5 )
ggsave("Species_imp_col_top10.pdf",exp_bar,width = 7.2,height =3.5 )

exp_bar


