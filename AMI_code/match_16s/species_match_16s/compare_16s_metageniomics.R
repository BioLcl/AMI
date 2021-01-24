get_corr_sig_list<-function(r_table,p_table,outputfile_name="./corr_sig_list.txt",
                            rowtable_name="table1",coltable_name="table2",
                            p_cutoff=0.05,p_method_name="p",relation_name="r"){
  if ((sum(p_table<p_cutoff,na.rm = T))==0){
    output="no sig result"
    print("no sig result")
    return(output)
  }else{
    corr_sig_list<-matrix(data=NA,nrow=1,ncol = 4)     
    
    colnames(corr_sig_list)<-c(rowtable_name,coltable_name,relation_name,p_method_name)
    
    for (i in (1:nrow(p_table))){
      for (j in (1:ncol(p_table))){
        if (is.na(p_table[i,j])){
          next
        } else if (p_table[i,j] < p_cutoff & r_table[i,j]>0) {
          new_res<-t(as.matrix(c(rownames(p_table)[i],colnames(p_table)[j],r_table[i,j],p_table[i,j])))
          corr_sig_list<-rbind(corr_sig_list,new_res)
          
        }
      }
    } 
    corr_sig_list<-corr_sig_list[-1,]
    corr_sig_list_output<-corr_sig_list
    write.table(corr_sig_list_output,outputfile_name,row.names = FALSE, sep = '\t', quote = FALSE, na = '')    
    return(corr_sig_list_output)
  }
}

impute_nodiag<-function(table,value){
  for (i in 1:nrow(table)){
    for (j in 1:ncol(table)){
      if (i!=j){
        table[i,j]=value
      }
    }
  }
  return(table)
}



get_corr_sig_label_matrix<-function(p_table,p_cutoff=0.05){
  sig_label_matrix<-matrix("",nrow(p_table),ncol(p_table))
  for (i in (1:nrow(p_table))){
    for (j in (1:ncol(p_table))){
      if (p_table[i,j] < p_cutoff) {
        sig_label_matrix[i,j]<-"*"
      }
    }
  } 
  return(sig_label_matrix)
}

heatmap_pdf_group<-function(r_table,outputfile_name="./corr.pdf",bk,width,hight,cluster=F,label_matrix=NA,
                            number_color="black",col_group=NA,ann_colors=NA){
  pdf(file =outputfile_name,width =width,height =hight )
  pheatmap(r_table, cluster_rows = F,
           scale="none",
           cluster_cols = cluster,
           #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
           color = c(colorRampPalette(colors =rdbu(100)[1:50])(length(bk)*2/12),colorRampPalette(colors =rdbu(100)[51:100])(length(bk)*10/12)),
           #color=rdbu(100),
           fontsize_row = 7.5,
           fontsize_col = 7.5,
           number_color=number_color,
           display_numbers=label_matrix,
           #legend_breaks=c(-0.5,-0.25,0,0.25,0.5,0.75,1),
           legend_breaks=c(-0.2,-0.1,0,0.25,0.5,0.75,1),
           annotation_col = col_group,
           annotation_colors =ann_colors,
           border_color=NA
  )
  dev.off()
}
library(psych)
library(base)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
setwd("D:/FangCloudV2/Zheng lab/Zheng lab????????/??????/ACS/ACS_code/match_16s/species_match_16s/")
rRNA<-read.csv("lefse_sig_asv.csv",stringsAsFactors = F,header = T,row.names = 1)

metagenomics<-read.delim("Species.txt",sep="\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
species_sig<-read.delim("species_lefse_sig_list.txt",header=F)
metagenomics<-metagenomics[species_sig$V1,]
both_ID<-intersect(colnames(rRNA),colnames(metagenomics))


metagenomics<-metagenomics[,both_ID]

metagenomics<-metagenomics[rowSums(metagenomics)!=0,]


metagenomics<-metagenomics/100

cor_data <- corr.test(t(metagenomics),t(rRNA),method = "spearman", adjust = "none")
cor_r<-cor_data$r
cor_p<-cor_data$p
table(cor_p<0.05)
corr_all_list<-get_corr_sig_list(cor_r,cor_p,"corr_all_sig_all_list_spearman_none_list.txt","species","asv",0.05,"p","r")





top_deter_matrix<-matrix(NA,1,4)
colnames(top_deter_matrix)<-c("metagenomics","16s","r","p")


for (i in 1:nrow(corr_all_list)){
  spe<-gsub("_"," ",corr_all_list[i,1])
  if (grepl(spe,corr_all_list[i,2])){
    top_deter_matrix<-rbind(top_deter_matrix,corr_all_list[i,])
  }
}
top_deter_matrix<-top_deter_matrix[-1,]
write.table(top_deter_matrix,"match_all_sig_spearman.txt",row.names = F,sep="\t")



grepl("Solobacterium",corr_all_list[1,2])


r_file<-list.files(path="./",pattern = "cor_r_*")
n=length(r_file)
merge.data<-read.csv(file=r_file[1],header = T,check.names = F)
colnames(merge.data)[1]="ID"

for (i in 2:n){
  newdata<-read.csv(file=r_file[i],header=T,check.names = F)
  colnames(newdata)[1]="ID"
  merge.data<-merge(merge.data,newdata,by="ID")
}
write.table(merge.data,"cor_r_all.txt",sep="\t",quote=F,row.names = F)


write.csv()
getwd()

cor_r_pic<-impute_nodiag(cor_r,NA)
cor_p_pic<-impute_nodiag(cor_p,NA)
corr_all_list<-get_corr_sig_list(cor_r_pic,cor_p_pic,"corr_all_list_pearson_BH_list.txt","genus","genus",1,"p","r")
corr_all_list<-as.data.frame(corr_all_list)
corr_all_list$BH<-p.adjust(as.numeric(corr_all_list$p),method = "BH")
write.table(corr_all_list,"corr_all_list_pearson_BH_list.txt",row.names = FALSE, sep = '\t', quote = FALSE, na = '')    




all_label<-get_corr_sig_label_matrix(cor_p_pic)
bk<-c(seq(-0.2,-0.01,by=0.01),seq(0,1,by=0.01))
bk


rdbu <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))
max(cor_r_pic,na.rm = T)
min(cor_r_pic,na.rm=T)
heatmap_pdf_group(cor_r_pic,"compare_raw_p.pdf",bk,15,15,F,all_label,"white")
dev.off()
######






