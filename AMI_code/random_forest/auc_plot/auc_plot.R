library(ggplot2)
setwd("D:/OneDrive/ACS/new_group2/random/plot")
species_auc<-read.table("auc_species_data.txt",header = T)
species_auc$specificities<-1-species_auc$specificities
species_auc$group<-"Species"

scfa_auc<-read.table("auc_scfa_mol_prop_data.txt",header = T)
scfa_auc$specificities<-1-scfa_auc$specificities
scfa_auc$group<-"SCFA Proportion"

TMAO_auc<-read.table("auc_TMAO_data.txt",header = T)
TMAO_auc$specificities<-1-TMAO_auc$specificities
TMAO_auc$group<-"TMAO related metabolites"

genus_auc<-read.table("auc_genus_data.txt",header = T)
genus_auc$specificities<-1-genus_auc$specificities
genus_auc$group<-"Genus"

# CAG_auc<-read.table("auc_CAG_data.txt",header = T)
# CAG_auc$specificities<-1-CAG_auc$specificities
# CAG_auc$group<-"CAG"

pathway_auc<-read.table("auc_pathway_data.txt",header = T)
pathway_auc$specificities<-1-pathway_auc$specificities
pathway_auc$group<-"Pathway"

auc_data<-rbind(species_auc,pathway_auc,genus_auc,scfa_auc,TMAO_auc)
auc_data$group<-factor(auc_data$group,levels=c("Species","Pathway","Genus","SCFA Proportion","TMAO related metabolites"))

species_auc$sensitivities


g <- ggplot(auc_data) + 
  geom_path(aes(x = specificities, y = sensitivities,color=group), size=0.5) + 
  labs(x = "False positive rate", y = "Ture positive rate",color="")  +
  scale_fill_discrete(breaks=c("Species","Genus","pathway","SCFA Proportion","TMAO related metabolites"))+
  scale_color_manual(values=c(rgb(224, 56, 55,max=255),
                              # rgb(106, 180, 45,max=255),
                              rgb(3, 172, 218,max=255),
                              "#ffd500",
                              "#bc6c25","#7209b7"),
                      labels=c("Species\nAUC=0.89 (0.81-0.97)\n",
                               # "CAG\nAUC=0.79 (0.67-0.91)\n",
                               "Pathway\nAUC=0.77 (0.64-0.90)\n",
                      "Genus\nAUC=0.87 (0.76-0.99)\n",
                      "SCFA Proportion\nAUC=0.81 (0.65-0.97)\n",
                      "TMAO related metabolites\nAUC=0.78 (0.66-0.91)"))+  geom_abline(intercept=0,slope=1,linetype=2 ,size=0.4)+

                              # labels=c("Species\n(AUC)=0.89 (0.81-0.97)\n",
                              # "Genus\n(AUC)=0.87\n",
                              # "CAG\n(AUC)=0.79\n",
                              # "Pathway\n(AUC)=0.77\n",
                              # "SCFA Proportion\n(AUC)=0.85\n",
                              # "TMAO related metabolites\n(AUC)=0.78"))+  geom_abline(intercept=0,slope=1,linetype=2 ,size=0.7)+
  theme(plot.title = element_text(face = 'bold',size=15))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        panel.background = element_rect(fill = NA, 
                                        colour = "black", linetype = "solid"), 
        plot.background = element_rect(colour = NA), 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA))+
  theme(panel.background = element_rect(size = 0.7), 
        plot.background = element_rect(size = 0.7))+
theme(axis.text.x = element_text(colour = "black"), 
      axis.text.y = element_text(colour = "black"),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 12))


ggsave("auc_95.png",g,width=5,height = 3)
ggsave("auc_95_scfa_mol.pdf",g,width=5,height = 3)

