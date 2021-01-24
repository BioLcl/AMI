library(vegan)
library(ggplot2)
library(ggsignif)
library(ggpubr)
setwd("D:/OneDrive/ACS/new_group2/medicine/Antibiotics")
species<-read.delim("Species.txt",row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
group<-read.delim("biaoxing.txt")
species<-species[,group$ID]



dist<-vegdist(t(species),"bray")
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

for (i in 2:6){
print (colnames(group)[i])
anosim_result<-anosim(distance,group[,i],permutations = 9999)
temp_anosim<-as.data.frame(t(matrix(c(colnames(group)[i],anosim_result$signif,anosim_result$statistic))))
colnames(temp_anosim)<-colnames(anosim_result_all)
anosim_result_all<-rbind(anosim_result_all,temp_anosim)
}
anosim_result_all<-anosim_result_all[-1,]
write.csv(anosim_result_all,"anosim_result_all.csv",row.names = F)
############ permanova
permanova_result_all<-matrix(NA,1,3)
colnames(permanova_result_all)<-c("Group","p-vaule","R2")
group<-group[,-10]

for (i in 2:13){
print (colnames(group)[i])
adonis_result<-adonis(distance~get(colnames(group)[i]),group,permutations = 9999)
adonis_result_detail<-as.data.frame(adonis_result$aov.tab)
temp_adonis<-as.data.frame(t(matrix(c(colnames(group)[i],adonis_result_detail[1,6],adonis_result_detail[1,5]))))
colnames(temp_adonis)<-colnames(permanova_result_all)
permanova_result_all<-rbind(permanova_result_all,temp_adonis)
}
permanova_result_all<-permanova_result_all[-1,]
write.csv(permanova_result_all,"permanova_result_all.csv",row.names = F)
permanova_result_all

alpha_data<-read.csv("alpha.csv",row.names = 1)
alpha_data<-t(alpha_data)
alpha_data<-alpha_data[group$ID,]
alpha_data<-cbind(alpha_data,group)
alpha_data$vfdb_all<-log10(alpha_data$vfdb_all)



nmds_data<-metaMDS(t(species), distance = 'bray', k = 3)


nmds_point<-data.frame(nmds_data$point)
colnames(nmds_point)<-c("NMDS1","NMDS2","NMDS3")
write.csv(nmds_point,"0-100nmds.csv")
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
  theme(axis.text.x = element_text(size=11, color = "black"),  # ɨփx֡ז̥´󐡣¬ҔςͬÀ�?        axis.text.y = element_text(size=15, color = "black"), 
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
ggsave("AMI_NMDS.pdf", nmds_plot_AMI, width = 2.7, height = 3.6)
ggsave("AMI_NMDS.png", nmds_plot_AMI, width = 2.7, height = 3.6)
####
nmds_plot_ST1TIN<-ggplot(nmds_plot_data, aes(NMDS1 ,NMDS2, color=ST1TIN,group = ST1TIN))+
  scale_color_manual(values=c("ST1TIN-"="#487ebe",
                              "ST1TIN+"="#e7352f"))+
  geom_point(aes(color = ST1TIN), size = 0.5, alpha = 1)+
  stat_ellipse(type = "t", linetype = 1)+
  labs(title="",x = "NMDS1", y = "NMDS2")+
  guides(shape = guide_legend(title = ''),
         color = guide_legend(title = ''))+
  mytheme+
  theme(legend.position="top")+
  guides(color = guide_legend(reverse = F,title=NULL))
nmds_plot_ST1TIN
ggsave("ST1TIN_NMDS.pdf", nmds_plot_ST1TIN, width = 2.7, height =3.4)
ggsave("ST1TIN_NMDS.png", nmds_plot_ST1TIN, width = 2.7, height =3.4)
########
nmds_plot_AHT<-ggplot(nmds_plot_data, aes(NMDS1 ,NMDS2, color=AHT,group = AHT))+
  scale_color_manual(values=c("AHT-"="#487ebe",
                              "AHT+"="#e7352f"))+
  geom_point(aes(color = AHT), size = 0.5, alpha = 1)+
  stat_ellipse(type = "t", linetype = 1)+
  labs(title="",x = "NMDS1", y = "NMDS2")+
  guides(shape = guide_legend(title = ''),
         color = guide_legend(title = ''))+
  mytheme+
  theme(legend.position="top")+
  guides(color = guide_legend(reverse = F,title=NULL))
nmds_plot_AHT
ggsave("AHT_NMDS.pdf", nmds_plot_AHT, width = 2.7, height =3.4)
ggsave("AHT_NMDS.png", nmds_plot_AHT,width = 2.7, height =3.4)
#######
nmds_plot_Antibiotics<-ggplot(nmds_plot_data, aes(NMDS1 ,NMDS2, color=Antibiotics,group = Antibiotics))+
  scale_color_manual(values=c("Antibiotics-"="#487ebe",
                              "Antibiotics+"="#e7352f"))+
  geom_point(aes(color = Antibiotics), size = 0.5, alpha = 1)+
  stat_ellipse(type = "t", linetype = 1)+
  labs(title="",x = "NMDS1", y = "NMDS2")+
  guides(shape = guide_legend(title = ''),
         color = guide_legend(title = ''))+
  mytheme+
  theme(legend.position="top")+
  guides(color = guide_legend(reverse = F,title=NULL))
nmds_plot_Antibiotics
ggsave("Antibiotics_NMDS.pdf", nmds_plot_Antibiotics, width = 2.7, height =3.4)
ggsave("Antibiotics_NMDS.png", nmds_plot_Antibiotics, width = 2.7, height =3.4)
######
# nmds_plot_PPI<-ggplot(nmds_plot_data, aes(NMDS1 ,NMDS2, color=PPI,group = PPI))+
#   scale_color_manual(values=c("PPI-"="#487ebe",
#                               "PPI+"="#e7352f"))+
#   geom_point(aes(color = PPI), size = 0.5, alpha = 1)+
#   stat_ellipse(type = "t", linetype = 1)+
#   labs(title="",x = "NMDS1", y = "NMDS2")+
#   guides(shape = guide_legend(title = ''),
#          color = guide_legend(title = ''))+
#   mytheme+
#   theme(legend.position="top")+
#   guides(color = guide_legend(reverse = F,title=NULL))
# nmds_plot_PPI
# ggsave("PPI_NMDS.pdf", nmds_plot_PPI, width = 2.5, height =3.3)
# ggsave("PPI_NMDS.png", nmds_plot_PPI, width = 2.5, height =3.3)

#########################
AMI_my_comparisons <- list(c("Case", "Control"))

data<-alpha_data
AMI_fig_Shannon<-ggplot(data,aes(x=AMI,y=shannon,fill=AMI))+
  stat_boxplot(geom = "errorbar",width=0.15,color="black")+
  geom_boxplot(lwd=0.4,fatten = 1,outlier.size = 0.5)+labs(y="Shannon Index")+
  scale_fill_manual(values=c(rgb(231, 53, 47,maxColorValue = 255),rgb(72, 126, 190,maxColorValue = 255)))+
  scale_x_discrete(limits=c("Control","Case"))+
  alpha_mytheme+labs(fill ="",x="")+
  stat_compare_means(comparisons=AMI_my_comparisons,method="wilcox.test",label.y = 5.65)+
  guides(fill="none")+
  ylim(1,5.75)
AMI_fig_Shannon
ggsave("AMI_shannon.pdf", AMI_fig_Shannon, width =2, height =3) 
ggsave("AMI_shannon.png", AMI_fig_Shannon, width =2, height =3) 

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

ggsave("AMI_simpson.pdf", AMI_fig_Simpson, width =2, height =3) 
ggsave("AMI_simpson.png", AMI_fig_Simpson, width =2, height =3) 
###
AMI_fig_VFs<-ggplot(data,aes(x=AMI,y=vfdb_all,fill=AMI))+
  stat_boxplot(geom = "errorbar",width=0.15,color="black")+
  geom_boxplot(lwd=0.4,fatten = 1,outlier.size = 0.5)+labs(y="Total VFs(log10)")+
  scale_fill_manual(values=c(rgb(231, 53, 47,maxColorValue = 255),rgb(72, 126, 190,maxColorValue = 255)))+
  scale_x_discrete(limits=c("Control","Case"))+
  alpha_mytheme+labs(fill ="",x="")+
  stat_compare_means(comparisons=AMI_my_comparisons,method="wilcox.test",label.y = 3.9)+
  guides(fill="none")+
  ylim(1,4)


ggsave("AMI_VFs.pdf", AMI_fig_VFs, width =2, height =3) 
ggsave("AMI_VFs.png", AMI_fig_VFs, width =2, height =3) 


#####.
ST1TIN_my_comparisons <- list(c("ST1TIN+", "ST1TIN-"))

data<-alpha_data
ST1TIN_fig_Shannon<-ggplot(data,aes(x=ST1TIN,y=shannon,fill=ST1TIN))+
  stat_boxplot(geom = "errorbar",width=0.15,color="black")+
  geom_boxplot(lwd=0.4,fatten = 1,outlier.size = 0.5)+labs(y="Shannon Index")+
  scale_fill_manual(values=c("ST1TIN-"="#487ebe",
                             "ST1TIN+"="#e7352f"))+
  scale_x_discrete(limits=c("ST1TIN-","ST1TIN+"))+
  alpha_mytheme+labs(fill ="",x="")+
  stat_compare_means(comparisons=ST1TIN_my_comparisons,method="wilcox.test",label.y = 5.65)+
  guides(fill="none")+
  ylim(1,5.75)

ST1TIN_fig_Shannon
ggsave("ST1TIN_shannon.pdf", ST1TIN_fig_Shannon, width =2.4, height =3) 
ggsave("ST1TIN_shannon.png", ST1TIN_fig_Shannon, width =2.4, height =3) 

##
ST1TIN_fig_Simpson<-ggplot(data,aes(x=ST1TIN,y=simpson,fill=ST1TIN))+
  stat_boxplot(geom = "errorbar",width=0.15,color="black")+
  geom_boxplot(lwd=0.4,fatten = 1,outlier.size = 0.5)+labs(y="Simpson Index")+
  scale_fill_manual(values=c("ST1TIN-"="#487ebe",
                             "ST1TIN+"="#e7352f"))+
  scale_x_discrete(limits=c("ST1TIN-","ST1TIN+"))+
  alpha_mytheme+labs(fill ="",x="")+
  stat_compare_means(comparisons=ST1TIN_my_comparisons,method="wilcox.test",label.y = 1.08)+
  guides(fill="none")+
  ylim(0,1.1)

ggsave("ST1TIN_simpson.pdf", ST1TIN_fig_Simpson, width =2.4, height =3) 
ggsave("ST1TIN_simpson.png", ST1TIN_fig_Simpson, width =2.4, height =3) 
#####
AHT_my_comparisons <- list(c("AHT+", "AHT-"))

data<-alpha_data
AHT_fig_Shannon<-ggplot(data,aes(x=AHT,y=shannon,fill=AHT))+
  stat_boxplot(geom = "errorbar",width=0.15,color="black")+
  geom_boxplot(lwd=0.4,fatten = 1,outlier.size = 0.5)+labs(y="Shannon Index")+
  scale_fill_manual(values=c("AHT-"="#487ebe",
                             "AHT+"="#e7352f"))+
  scale_x_discrete(limits=c("AHT-","AHT+"))+
  alpha_mytheme+labs(fill ="",x="")+
  stat_compare_means(comparisons=AHT_my_comparisons,method="wilcox.test",label.y = 5.65)+
  guides(fill="none")+
  ylim(1,5.75)

AHT_fig_Shannon
ggsave("AHT_shannon.pdf", AHT_fig_Shannon, width =2.4, height =3) 
ggsave("AHT_shannon.png", AHT_fig_Shannon,width =2.4, height =3) 

##
AHT_fig_Simpson<-ggplot(data,aes(x=AHT,y=simpson,fill=AHT))+
  stat_boxplot(geom = "errorbar",width=0.15,color="black")+
  geom_boxplot(lwd=0.4,fatten = 1,outlier.size = 0.5)+labs(y="Simpson Index")+
  scale_fill_manual(values=c("AHT-"="#487ebe",
                             "AHT+"="#e7352f"))+
  scale_x_discrete(limits=c("AHT-","AHT+"))+
  alpha_mytheme+labs(fill ="",x="")+
  stat_compare_means(comparisons=AHT_my_comparisons,method="wilcox.test",label.y = 1.08)+
  guides(fill="none")+
  ylim(0,1.1)

ggsave("AHT_simpson.pdf", AHT_fig_Simpson,width =2.4, height =3) 
ggsave("AHT_simpson.png", AHT_fig_Simpson,width =2.4, height =3) 
####
Antibiotics_my_comparisons <- list(c("Antibiotics+", "Antibiotics-"))

data<-alpha_data
Antibiotics_fig_Shannon<-ggplot(data,aes(x=Antibiotics,y=shannon,fill=Antibiotics))+
  stat_boxplot(geom = "errorbar",width=0.15,color="black")+
  geom_boxplot(lwd=0.4,fatten = 1,outlier.size = 0.5)+labs(y="Shannon Index")+
  scale_fill_manual(values=c("Antibiotics-"="#487ebe",
                             "Antibiotics+"="#e7352f"))+
  scale_x_discrete(limits=c("Antibiotics-","Antibiotics+"))+
  alpha_mytheme+labs(fill ="",x="")+
  stat_compare_means(comparisons=Antibiotics_my_comparisons,method="wilcox.test",label.y = 5.65)+
  guides(fill="none")+
  ylim(1,5.75)

Antibiotics_fig_Shannon
ggsave("Antibiotics_shannon.pdf", Antibiotics_fig_Shannon, width =2.4, height =3) 
ggsave("Antibiotics_shannon.png", Antibiotics_fig_Shannon, width =2.4, height =3) 

##
Antibiotics_fig_Simpson<-ggplot(data,aes(x=Antibiotics,y=simpson,fill=Antibiotics))+
  stat_boxplot(geom = "errorbar",width=0.15,color="black")+
  geom_boxplot(lwd=0.4,fatten = 1,outlier.size = 0.5)+labs(y="Simpson Index")+
  scale_fill_manual(values=c("Antibiotics-"="#487ebe",
                             "Antibiotics+"="#e7352f"))+
  scale_x_discrete(limits=c("Antibiotics-","Antibiotics+"))+
  alpha_mytheme+labs(fill ="",x="")+
  stat_compare_means(comparisons=Antibiotics_my_comparisons,method="wilcox.test",label.y = 1.08)+
  guides(fill="none")+
  ylim(0,1.1)

ggsave("Antibiotics_simpson.pdf", Antibiotics_fig_Simpson, width =2.4, height =3) 
ggsave("Antibiotics_simpson.png", Antibiotics_fig_Simpson, width =2.4, height =3) 
######
# PPI_my_comparisons <- list(c("PPI+", "PPI-"))
# 
# data<-alpha_data
# PPI_fig_Shannon<-ggplot(data,aes(x=PPI,y=shannon,fill=PPI))+
#   stat_boxplot(geom = "errorbar",width=0.15,color="black")+
#   geom_boxplot(lwd=0.4,fatten = 1,outlier.size = 0.5)+labs(y="Shannon Index")+
#   scale_fill_manual(values=c(rgb(231, 53, 47,maxColorValue = 255),rgb(72, 126, 190,maxColorValue = 255)))+
#   scale_x_discrete(limits=c("PPI-","PPI+"))+
#   alpha_mytheme+labs(fill ="",x="")+
#   stat_compare_means(comparisons=PPI_my_comparisons,method="wilcox.test",label.y = 5.65)+
#   guides(fill="none")+
#   ylim(1,5.75)
# 
# PPI_fig_Shannon
# ggsave("PPI_shannon.pdf", PPI_fig_Shannon,width =2.4, height =2.9)
# ggsave("PPI_shannon.png", PPI_fig_Shannon, width =2.4, height =2.9)
# 
# ##
# PPI_fig_Simpson<-ggplot(data,aes(x=PPI,y=simpson,fill=PPI))+
#   stat_boxplot(geom = "errorbar",width=0.15,color="black")+
#   geom_boxplot(lwd=0.4,fatten = 1,outlier.size = 0.5)+labs(y="Simpson Index")+
#   scale_fill_manual(values=c(rgb(231, 53, 47,maxColorValue = 255),rgb(72, 126, 190,maxColorValue = 255)))+
#   scale_x_discrete(limits=c("PPI-","PPI+"))+
#   alpha_mytheme+labs(fill ="",x="")+
#   stat_compare_means(comparisons=PPI_my_comparisons,method="wilcox.test",label.y = 1.08)+
#   guides(fill="none")+
#   ylim(0,1.1)
# 
# ggsave("PPI_simpson.pdf", PPI_fig_Simpson,width =2.4, height =2.9)
# ggsave("PPI_simpson.png", PPI_fig_Simpson, width =2.4, height =2.9)
