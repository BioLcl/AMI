library(psych)
library(base)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
rdbu <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))
heatmap_pdf_group<-function(r_table,outputfile_name="./corr.pdf",bk,width,hight,cluster=F,label_matrix=NA,number_color="black",col_group=NA,ann_colors=NA){
  pdf(file =outputfile_name,width =width,height =hight )
  pheatmap(r_table, cluster_rows = F,
           scale="none",
           #angle_col = 315,
           cluster_cols = cluster,
           #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
           color = c(colorRampPalette(colors =rdbu(100)[1:50])(length(bk)*5/11)[round(length(colorRampPalette(colors =rdbu(100)[1:50])(length(bk)*5/11))*1/5):round(length(colorRampPalette(colors =rdbu(100)[1:50])(length(bk)*5/11)))],
                     colorRampPalette(colors =rdbu(100)[51:100])(length(bk)*6/11)[1:round(length(colorRampPalette(colors =rdbu(100)[51:100])(length(bk)*6/11))*6/6)]),
           #color=rdbu(100),
           fontsize_row = 7.5,
           fontsize_col = 7.5,
           number_color=number_color,
           display_numbers=label_matrix,
           legend_breaks=-5:6,
           annotation_col = col_group,
           annotation_colors =ann_colors,
           border_color=NA
  )
  dev.off()
}

int <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}

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

get_corr_sig_list<-function(r_table,p_table,outputfile_name="./corr_sig_list.txt",rowtable_name="table1",coltable_name="table2",p_cutoff=0.05,p_method_name="p",relation_name="r"){
  if ((sum(p_table<p_cutoff))==0){
    output="no sig result"
    print("no sig result")
    return(output)
  }else{
    corr_sig_list<-matrix(data=NA,nrow=1,ncol = 4)  
    
    colnames(corr_sig_list)<-c(rowtable_name,coltable_name,relation_name,p_method_name)
    
    for (i in (1:nrow(p_table))){
      for (j in (1:ncol(p_table))){
        if (p_table[i,j] < p_cutoff) {
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
######
p_cut<-0.1
setwd("D:/AMI_code/Relat/clinical_species")
meta<-read.delim("meta.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
daixie_num<-ncol(meta)
meta<-as.data.frame(meta)
meta$ID<-rownames(meta)
genus<-read.delim("pathway_all.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

genus<-t(genus)
genus<-as.data.frame(genus)
colnames(genus)
#####
genus_o<-genus
genus_l<-read.delim("lefse_pathway_sig.txt",  sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,header=F)
genus_t<-genus[,match(genus_l[,1],colnames(genus))]



genus<-genus_t
#####
genus<-as.data.frame((genus))
genus$ID<-rownames(genus)
#####

biaoxing<-read.delim("biaoxing.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
biaoxing_num<-ncol(biaoxing)
biaoxing$ID<-rownames(biaoxing)

data_meta <- merge(biaoxing, meta, by = "ID",all=FALSE)

data_meta <- merge(data_meta, genus, by = "ID",all=FALSE)

rownames(data_meta)<-data_meta$ID
data_meta<-data_meta[order(data_meta$ID),]

daixie_num
biaoxing_num
out_daixie<-data_meta[,(2+biaoxing_num):(1+biaoxing_num+daixie_num)]
all_id<-data_meta$ID
# daixie_pos<-match(all_id,daixie$ID)
# out_daixie<-daixie[daixie_pos,-ncol(daixie)]
rm_daixie<-names(out_daixie)[colSums(out_daixie,na.rm = TRUE)==0]

print (paste("remove",rm_daixie,sep=": "))
out_daixie<-out_daixie[,!colSums(out_daixie,na.rm = TRUE)==0]
out_daixie<-out_daixie[order(rownames(out_daixie)),]


#genus
genus_pos<-match(all_id,genus$ID)
out_genus<-genus[genus_pos,-ncol(genus)]
#rmove abundance=0 genus

rm_genus<-names(genus)[colSums(out_genus)==0]
print (paste("remove tax",rm_genus,sep=": "))
out_genus<-out_genus[,!colSums(out_genus)==0]
out_genus<-out_genus[order(rownames(out_genus)),]
biaoxing<-data_meta[,2:(1+biaoxing_num)]

write.csv(out_daixie,"out_daixie.csv")
write.csv(out_genus,"out_genus.csv")
write.csv(data_meta,"data_meta.csv")

out_genus_o<-out_genus

######################### raw
out_genus_notarns<-out_genus
out_daixie_notrans<-out_daixie
impute_daixie<-function(x){
  x[is.na(x)]<-min(x,na.rm=T)/2
  return(x)
}

############################# all int
out_daixie<-out_daixie_notrans

out_genus<-out_genus_notarns
out_genus<-as.matrix(out_genus)

for (i in 1:ncol(out_genus)){
  out_genus[,i]<-int(out_genus[,i])
}
for (i in 1:ncol(out_daixie)){
  out_daixie[,i]<-int(out_daixie[,i])
}



lm_all_log<-cbind(biaoxing,out_daixie,out_genus)


fit_lmall_result_log<-matrix(NA,1,4)
colnames(fit_lmall_result_log)<-c("metabolome","genus","p","beta")
fit_lmall_result_log[1,1]<-"fit_lmall_result_log"
ncol(out_daixie)
##
lm_all_beta_matrix<-matrix(NA,ncol(out_daixie),ncol(out_genus))
colnames(lm_all_beta_matrix)<-colnames(out_genus)
rownames(lm_all_beta_matrix)<-colnames(out_daixie)
lm_all_padj_matrix<-lm_all_beta_matrix
#
library(broom)



cor_detail_result<-matrix(NA,1,6)
colnames(cor_detail_result)<-c("Metabolome","Species","beta","SE","t-value","p-value")

for (i in (biaoxing_num+1):(biaoxing_num+daixie_num)){
  for (j in (biaoxing_num+daixie_num+1):(ncol(lm_all_log))){
    fit<-lm(get(colnames(lm_all_log)[i])~get(colnames(lm_all_log)[j])+sex+age+BMI,data=lm_all_log)
    
    fit_data_detail<-as.data.frame(tidy(fit))

    lm_all_beta_matrix[i-biaoxing_num,j-biaoxing_num-daixie_num]<-fit_data_detail[2,4]
    lm_all_padj_matrix[i-biaoxing_num,j-biaoxing_num-daixie_num]<-fit_data_detail[2,5]
    fit_sig<-t(as.matrix(c(colnames(lm_all_log)[i],colnames(lm_all_log)[j],fit_data_detail[2,5],fit_data_detail[2,4])))
    fit_lmall_result_log<-rbind(fit_lmall_result_log,fit_sig)
    temp_stat_result<-t(as.matrix(c(colnames(lm_all_log)[i],colnames(lm_all_log)[j],fit_data_detail[2,2],fit_data_detail[2,3],fit_data_detail[2,4],fit_data_detail[2,5])))
    cor_detail_result<-rbind(cor_detail_result,temp_stat_result)
  }
}
cor_detail_result<-cor_detail_result[-1,]
cor_detail_result<-as.data.frame(cor_detail_result)
cor_detail_result$fdr<-p.adjust(as.numeric(cor_detail_result$`p-value`),method="BH")
write.table(cor_detail_result,"corlation_dtail_result.txt",sep="\t",row.names = F,quote=F)

write.csv(lm_all_beta_matrix,"species_linchuang_beta.csv")
write.csv(lm_all_padj_matrix,"species_scfa_p.csv")

#############3
################################################ aaa

lm_raw_p<-lm_all_padj_matrix
table(lm_raw_p<0.05)
fdr_padj<-matrix(p.adjust(lm_all_padj_matrix,method="BH"),nrow(lm_all_beta_matrix),ncol(lm_all_beta_matrix))
rownames(lm_all_padj_matrix)<-rownames(lm_all_beta_matrix)
colnames(lm_all_padj_matrix)<-colnames(lm_all_beta_matrix)

lm_all_padj_matrix<-fdr_padj
dev.off()
p_cut<-0.1
bk <- c(seq(-4,-0.01,by=0.01),seq(0,6,by=0.01))
if (sum(lm_all_padj_matrix<p_cut)>1){
  corr_sig_log<-get_corr_matrix_sig(lm_all_beta_matrix,lm_all_padj_matrix,p_cut,r_sig = c(NA,1:15))
  corr_sig_p_log<-corr_sig_log$p_sig_table
  corr_sig_r_log<-corr_sig_log$r_sig_table
  corr_sig_label<-get_corr_sig_label_matrix(corr_sig_p_log,p_cut)
  
  col_group_ord<-match(colnames(corr_sig_r_log),genus_l[,1])
  
  col_group<-genus_l[,3][col_group_ord]
  col_group<-as.data.frame(col_group)
  colnames(col_group)<-"Group"
  rownames(col_group)<-genus_l[,1][col_group_ord]

  ann_colors = list(
    Group = c(Case =	"#696969", Control =rgb(	220,220,220,maxColorValue = 255))
  )

  heatmap_sig_cluster=paste0("log_","sig_cluster_","heatmap.pdf")
  heatmap_pdf_group(corr_sig_r_log,heatmap_sig_cluster,bk,8.5,4.5,T,corr_sig_label,"white",col_group,ann_colors)
  dev.off()
}

