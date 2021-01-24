get_corr_matrix_sig<-function(r_table,p_table,p_cutoff=0.05,r_sig=NA,c_sig=NA){
  if ((sum(p_table<p_cutoff))==0){
    output=FALSE
    print("no sig result")
    return(output)
  }else{
    for (i in (1:nrow(p_table))){
      for (j in (1:ncol(p_table))){
        if (p_table[i,j] < p_cutoff) {
          r_sig<-c(r_sig,i)
          c_sig<-c(c_sig,j)
        }
      }
    } 
    r_sig<-r_sig[-1]
    c_sig<-c_sig[-1]
    r_sig<-sort(as.numeric(unique(r_sig)))
    c_sig<-sort(as.numeric(unique(c_sig)))
    r_sig_table<-r_table[r_sig,c_sig]
    p_sig_table<-p_table[r_sig,c_sig]
    corr_sig_result<-list(r_sig_table,p_sig_table)
    
    names(corr_sig_result)<-c("r_sig_table","p_sig_table")
    return(corr_sig_result)
  }
}

setwd("D:/FangCloudV2/Zheng lab/Zheng lab共享资料/刘成林/ACS/ACS_code/random_forest/WGCNA/module_clinical_corelation")

SCFA<-read.delim("scfa_mol_ratio.txt",header=T,sep="\t",check.names = F)
pathway<-read.delim("pathway_module_abundance.txt",header=T,sep="\t",check.names = F)
species<-read.delim("Species_module_abundance.txt",header=T,sep="\t",check.names = F)
TMAO<-read.delim("TMAO.txt",header=T,sep="\t",check.names = F)
linchuang<-read.delim("linchuang.txt",header=T,sep="\t",check.names = F)
Genus<-read.delim("16s_genus_module.txt",header=T,sep="\t",check.names = F)


merge_data<-merge(species,Genus,by="ID",all=T)
merge_data<-merge(merge_data,pathway,by="ID",all=T)
merge_data<-merge(merge_data,linchuang,by="ID",all=T)
merge_data<-merge(merge_data,TMAO,by="ID",all=T)
merge_data<-merge(merge_data,SCFA,by="ID",all=T)
rownames(merge_data)<-merge_data$ID

library(psych)
ana_data<-merge_data[,-1]

#ana_list<-read.delim("ana_list.txt",sep="\t")

#ana_data<-ana_data[,intersect(ana_list$ID,colnames(ana_data))]

module_data<-ana_data[,1:9]
other_data<-ana_data[,10:22]
other_data_2<-ana_data[,23:ncol(ana_data)]


module_cor<-corr.test(module_data,module_data,method="spearman",adjust = "none")
module_cor_r<-module_cor$r
module_cor_p<-module_cor$p


trait_cor<-corr.test(module_data,other_data,method="spearman",adjust = "none")
trait_cor_r<-trait_cor$r
trait_cor_p<-trait_cor$p
newdata<-cbind(module_data,other_data)
trait_cor2<-corr.test(newdata,other_data_2,method="spearman",adjust = "none")
trait_cor2_r<-trait_cor2$r
trait_cor2_p<-trait_cor2$p

#cor_result<-corr.test(ana_data,ana_data,method="spearman",adjust = "none")
# cor_p<-cor_result$p
# cor_r<-cor_result$r
p_cut<-0.05
#corr_sig<-get_corr_matrix_sig(module_cor_r,module_cor_p,p_cut)

##
cor_r<-module_cor_r
cor_p<-module_cor_p


data_cor_index<-upper.tri(cor_r,diag = F)
cor_net_prep<-matrix(NA,1,4)
colnames(cor_net_prep)<-c("from","to","p","r")
for (i in 1:nrow(data_cor_index)){
  for (j in 1:ncol(data_cor_index)){
    if (data_cor_index[i,j]  ){
      new_res<-t(as.matrix(c(rownames(cor_p)[i],colnames(cor_p)[j],cor_p[i,j],cor_r[i,j])))
      cor_net_prep<-rbind(cor_net_prep,new_res)
   }
  }
}
###




cor_net_prep_trait<-matrix(NA,1,4)
colnames(cor_net_prep_trait)<-c("from","to","p","r")
for (i in 1:nrow(trait_cor_r)){
  for (j in 1:ncol(trait_cor_r)){
      new_res<-t(as.matrix(c(rownames(trait_cor_r)[i],colnames(trait_cor_r)[j],trait_cor_p[i,j],trait_cor_r[i,j])))
      cor_net_prep_trait<-rbind(cor_net_prep_trait,new_res)
  }
}
cor_net_prep_trait<-cor_net_prep_trait[-1,]
###
cor_net_prep_trait2<-matrix(NA,1,4)
colnames(cor_net_prep_trait2)<-c("from","to","p","r")
for (i in 1:nrow(trait_cor2_r)){
  for (j in 1:ncol(trait_cor2_r)){
    new_res<-t(as.matrix(c(rownames(trait_cor2_r)[i],colnames(trait_cor2_r)[j],trait_cor2_p[i,j],trait_cor2_r[i,j])))
    cor_net_prep_trait2<-rbind(cor_net_prep_trait2,new_res)
  }
}
cor_net_prep_trait2<-cor_net_prep_trait2[-1,]
####
cor_net_prep<-cor_net_prep[-1,]
cor_net_prep_raw<-cor_net_prep
cor_net_prep<-rbind(cor_net_prep,cor_net_prep_trait)
cor_net_prep<-rbind(cor_net_prep,cor_net_prep_trait2)

cor_net_prep<-as.data.frame(cor_net_prep)
cor_net_prep$padj<-p.adjust(cor_net_prep$p,method = "BH")
write.table(cor_net_prep,"cor_all_result_nopath.txt",row.names = F,sep="\t",quote = F)

cor_net_prep_sig<-cor_net_prep[cor_net_prep$padj<0.05,]





cor_net_prep_sig$abs_r<-abs(as.numeric(as.character(cor_net_prep_sig$r)))
table(cor_net_prep_sig$abs_r<=0.3)


cor_net_prep_sig$cor<-NA
cor_net_prep_sig$cor[as.numeric(as.character(cor_net_prep_sig$r))>0]<-"pos"
cor_net_prep_sig$cor[as.numeric(as.character(cor_net_prep_sig$r))<0]<-"neg"
cor_net_prep_sig$cor_large<-"small"
cor_net_prep_sig$cor_large[as.numeric(as.character(cor_net_prep_sig$r))>0.3]<-"large_pos"
cor_net_prep_sig$cor_large[as.numeric(as.character(cor_net_prep_sig$r))<(-0.3)]<-"large_neg"

write.table(cor_net_prep_sig,"edge_adj_nopath.txt",sep="\t",quote=F,row.names = F)
