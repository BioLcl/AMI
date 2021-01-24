library(vegan)
library(ggplot2)
library(ggsignif)
library(ggpubr)
setwd("C:/Users/sunbo/OneDrive/FangCloudV2/Zheng lab/Zheng lab共享资料/刘成林/ACS/ACS_code/alpha_beta/16S")
species<-read.delim("asv_rarefied.txt",row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
species<-as.matrix(species) # raw data was counts for ASVs (each column is an independent sample)
species<-prop.table(species,2) # claculate the relative abundance (here 2 means divided by column sum)
species<-species*100 # transform to percentage
species<-as.data.frame(species)
colSums(species) # to check the result

group<-read.delim("sample-metadata.txt")
species<-species[,group$ID] # order



dist<-vegdist(t(species),"bray") # 计算bray curtis距离
distance<-as.matrix(dist)
colSums(distance)

### order 
distance<-distance[group$ID,]
distance<-distance[,group$ID]
rownames(distance)==colnames(distance)
rownames(distance)==group$ID

set.seed(1234)
######### anosim
anosim_result_all<-matrix(NA,1,3)
colnames(anosim_result_all)<-c("Group","p-vaule","statistic")

for (i in 2){
print (colnames(group)[i])
anosim_result<-anosim(distance,group[,i],permutations = 9999)
temp_anosim<-as.data.frame(t(matrix(c(colnames(group)[i],anosim_result$signif,anosim_result$statistic))))
colnames(temp_anosim)<-colnames(anosim_result_all)
anosim_result_all<-rbind(anosim_result_all,temp_anosim)
}

anosim_result_all<-anosim_result_all[-1,]

write.csv(anosim_result_all,"anosim_result_all.csv",row.names = F)
############ permanova


alpha_data<-read.delim("alpha.txt",row.names = 1) # alpha diversity was calculated previously
alpha_data<-alpha_data[group$ID,]
alpha_data<-cbind(alpha_data,group)



nmds_data<-metaMDS(t(species), distance = 'bray', k = 3)




nmds_point<-data.frame(nmds_data$point)
colnames(nmds_point)<-c("NMDS1","NMDS2","NMDS3")
#write.csv(nmds_point,"0-100nmds.csv")
nmds_point<-read.csv("0-100nmds.csv",row.names = 1)
rownames(nmds_point)==group$ID
nmds_plot_data<-cbind(nmds_point,group)
#######################################################
mytheme<-theme(panel.grid.major = element_line(linetype = "blank"), 
               panel.grid.minor = element_line(linetype = "blank"),
               panel.grid = element_line(color = 'white', linetype = 2, size = 0.1), 
               panel.background = element_rect(color = 'black', fill = 'transparent'),
               legend.key = element_rect(fill = 'transparent'),
               plot.title = element_text(hjust = 0.5,size = 15),
               axis.title = element_text(size = 12),
               axis.text.x = element_text(size = 9), 
               axis.text.y = element_text(size = 9),
               legend.text = element_text(size = 11)) 
#######################
alpha_mytheme<-theme(panel.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.border = element_rect(colour="black",fill=NA)) +
  theme(axis.text.x = element_text(size=11, color = "black"),  # ɨփx֡ז̥´󐡣¬ҔςͬÀ???        axis.text.y = element_text(size=15, color = "black"), 
        axis.title.x = element_text(size=15, color = "black"), 
        axis.title.y= element_text(size=12,  color = "black"),
        legend.title = element_text(size=15,  color = "black"))+
  # ,
  # axis.text.x.bottom  = element_text(angle = 0, hjust = 0.5),
  # legend.text = element_text(size=13)) +# ·֗錵Ўשʱ£¬position = 'dodge'
  theme(panel.grid =element_blank())  +
  theme(axis.line = element_line(size=0.3,color  = "black")) +
  theme(panel.grid.major=element_line(colour=NA), panel.border = element_blank())+
  theme(axis.line = element_line(arrow = arrow(length = unit(0.15, 'cm'))))+
  theme(axis.ticks.length=unit(-0.1, "cm"), 
       axis.text.x = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")),
        legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(15, "pt"),
        legend.key.height = unit(15, "pt"),
        legend.key.width = unit(15, "pt"),
        legend.text = element_text(size=13))





#######
nmds_plot_AMI<-ggplot(nmds_plot_data, aes(NMDS1 ,NMDS2, color=AMI,group = AMI))+
  scale_color_manual(values=c(Control="#487ebe",
                              Case="#e7352f"))+
  geom_point(aes(color = AMI), size = 1, alpha = 1)+
  stat_ellipse(type = "t", linetype = 1)+
  labs(title="",x = "NMDS1", y = "NMDS2")+
  guides(shape = guide_legend(title = ''),
         color = guide_legend(title = ''))+
  mytheme+
  theme(legend.position="top")+
  guides(color = guide_legend(reverse = T,title=NULL))
nmds_plot_AMI
ggsave("16s_AMI_NMDS.pdf", nmds_plot_AMI, width = 2.7, height = 3.6)
ggsave("16s_AMI_NMDS.png", nmds_plot_AMI, width = 2.7, height = 3.6)


AMI_my_comparisons <- list(c("Case", "Control"))

data<-alpha_data
AMI_fig_Shannon<-ggplot(data,aes(x=AMI,y=shannon,fill=AMI))+
  stat_boxplot(geom = "errorbar",width=0.15,color="black")+
  geom_boxplot(lwd=0.4,fatten = 1,outlier.size = 0.5)+labs(y="Shannon Index")+
  scale_fill_manual(values=c(rgb(231, 53, 47,maxColorValue = 255),rgb(72, 126, 190,maxColorValue = 255)))+
  scale_x_discrete(limits=c("Control","Case"))+
  alpha_mytheme+labs(fill ="",x="")+
  stat_compare_means(comparisons=AMI_my_comparisons,method="wilcox.test",label.y = 7.3)+
  guides(fill="none")+
  ylim(1,7.5)
AMI_fig_Shannon
ggsave("16s_AMI_shannon.pdf", AMI_fig_Shannon, width =2, height =3) 
ggsave("16s_AMI_shannon.png", AMI_fig_Shannon, width =2, height =3) 
max(data$simpson)
##
AMI_fig_Simpson<-ggplot(data,aes(x=AMI,y=simpson,fill=AMI))+
  stat_boxplot(geom = "errorbar",width=0.15,color="black")+
  geom_boxplot(lwd=0.4,fatten = 1,outlier.size = 0.5)+labs(y="Simpson Index")+
  scale_fill_manual(values=c(rgb(231, 53, 47,maxColorValue = 255),rgb(72, 126, 190,maxColorValue = 255)))+
  scale_x_discrete(limits=c("Control","Case"))+
  alpha_mytheme+labs(fill ="",x="")+
  stat_compare_means(comparisons=AMI_my_comparisons,method="wilcox.test",label.y = 1.08)+
  guides(fill="none")+
  ylim(0,1.1)
AMI_fig_Simpson
ggsave("16s_AMI_simpson.pdf", AMI_fig_Simpson, width =2, height =3) 
ggsave("16s_AMI_simpson.png", AMI_fig_Simpson, width =2, height =3) 
###
