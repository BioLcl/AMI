
library(psych)
library(RColorBrewer)
library(ggplot2)
setwd("D:/FangCloudV2/Zheng lab/Zheng lab共享资料/刘成林/ACS/ACS_code/match_16s/species_16s_match_plot/")
rRNA<-read.delim("asv_match_profile.txt",sep="\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

metagenomics<-read.delim("Species.txt",sep="\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

both_ID<-intersect(colnames(rRNA),colnames(metagenomics))


metagenomics<-metagenomics[,both_ID]

metagenomics<-metagenomics[rowSums(metagenomics)!=0,]


metagenomics<-metagenomics/100
dir.create("compare_plot_diff_sep")

for (i in 1:nrow(rRNA)){

meta_ana<-metagenomics[unlist(strsplit(rownames(rRNA)[i],"[|]"))[2],]
rna_ana<-rRNA[i,]
ana_data<-t(rbind(meta_ana,rna_ana))

colnames(ana_data)[1]
colnames(ana_data)[2]
ana_data<-as.data.frame(ana_data)
asv_pic_id<-unlist(strsplit(rownames(rRNA)[i],"[|]"))[1]
asv_species<-unlist(strsplit(rownames(rRNA)[i],"[|]"))[2]
data_point<-ggplot(ana_data,aes(x=get(colnames(ana_data)[1]),y=get(colnames(ana_data)[2])))+
  geom_point(size=1.5)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(axis.text.x = element_text(size=15, color = "black"),  # 设置x轴字体大小，以下同理
        axis.text.y = element_text(size=15, color = "black"),
        axis.title.x = element_text(size=18, color = "black"),
        axis.title.y= element_text(size=18,  color = "black"),
        legend.title = element_text(size=15,  color = "black"))+
  theme(panel.grid =element_blank())  +
  theme(axis.line = element_line(size=0.8, colour = "black")) +
  theme(panel.grid.major=element_line(colour=NA), panel.border = element_blank())+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"),angle=90), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        legend.background = element_rect(fill="white"),
        legend.key.size = unit(15, "pt"),
        legend.key.height = unit(15, "pt"),
        legend.key.width = unit(15, "pt"))+
  labs(y=paste0(asv_pic_id,"\n",gsub("_"," ",asv_species),"\n","(16S Sequencing)"),x=paste0(sub("_"," ",colnames(ana_data)[1]),"\n","(Shotgun Sequencing)"))+
  geom_smooth(method = lm,colour="black")


pic_id_p<-gsub("[|]","_",colnames(ana_data)[2])
pic_id<-gsub(" ","_",pic_id_p)
pic_id<-paste0(asv_species,"_",asv_pic_id)
paste0(sub("_"," ",colnames(ana_data)[1]),"\n","shotgun")
ggsave(paste0("./compare_plot_diff_sep/",pic_id,".png"),data_point,width =7,height =7)
ggsave(paste0("./compare_plot_diff_sep/",pic_id,".pdf"),data_point,width =7,height =7)

}

