source("vcf_eval2.R")
library(circlize)
library(ComplexHeatmap)
library(ggpubr)
#####################
#Synthetic3 dataset #
#####################

setwd('/mnt/Data/benchmark_somatic/papier')
base = read.csv('data/NGV3/thruth_vcf/synthetic_challenge_set3_tumor.region.recode.vcf' ,sep="\t",header=F,comment.char="#")

#Dedup + BQSR
test_snp = read.csv('data/NGV3/SomaticAnalysis/tumor_NGv3_S1_vs_normal_NGv3_S1.somatic.snv.merge.filtered.vcf',sep="\t",header=F,comment.char="#")
test_indel = read.csv('data/NGV3/SomaticAnalysis/tumor_NGv3_S1_vs_normal_NGv3_S1.somatic.indel.merge.filtered.vcf',sep="\t",header=F,comment.char="#")

#Dedup
#test_snp = read.csv('data/NGV3_Dedup/SomaticAnalysis/tumor_NGv3_S1_vs_normal_NGv3_S1.somatic.snv.merge.filtered.vcf',sep="\t",header=F,comment.char="#")
#test_indel = read.csv('data/NGV3_Dedup/SomaticAnalysis/tumor_NGv3_S1_vs_normal_NGv3_S1.somatic.indel.merge.filtered.vcf',sep="\t",header=F,comment.char="#")

#BQSR
#test_snp = read.csv('data/NGV3_BQSR/SomaticAnalysis/tumor_NGv3_S1_vs_normal_NGv3_S1.somatic.snv.merge.filtered.vcf',sep="\t",header=F,comment.char="#")
#test_indel = read.csv('data/NGV3_BQSR/SomaticAnalysis/tumor_NGv3_S1_vs_normal_NGv3_S1.somatic.indel.merge.filtered.vcf',sep="\t",header=F,comment.char="#")

#None
#test_snp = read.csv('data/NGV3_sans_BQSR_Dedup/SomaticAnalysis/tumor_NGv3_S1_vs_normal_NGv3_S1.somatic.snv.merge.filtered.vcf',sep="\t",header=F,comment.char="#")
#test_indel = read.csv('data/NGV3_sans_BQSR_Dedup/SomaticAnalysis/tumor_NGv3_S1_vs_normal_NGv3_S1.somatic.indel.merge.filtered.vcf',sep="\t",header=F,comment.char="#")

base$ID = paste(base$V1,":",base$V2,":",base$V4,">",base$V5,sep="")
test_snp$ID = paste(test_snp$V1,":",test_snp$V2,":",test_snp$V4,">",test_snp$V5,sep="")
test_indel$ID = paste(test_indel$V1,":",test_indel$V2,":",test_indel$V4,">",test_indel$V5,sep="")
base_snp = base[which(base$V4 %in% c("A","T","G","C") & base$V5 %in% c("A","T","G","C")),]
base_indel = base[-which(base$V4 %in% c("A","T","G","C") & base$V5 %in% c("A","T","G","C")),]
base_snp = get_vaf(base_snp)
base_indel$VAF=1
mat_vc_snp = get_mat_vc(test_snp,base_snp)
roc_vc_snp = roc_vc(mat_vc_snp)
roc_vc_snp = format_roc_vc(roc_vc_snp)
roc_vc_cons_snp = best_roc_vote(mat_vc_snp,base_snp)
roc_vc_cons_snp_ngv3 = roc_vc_cons_snp

mat_vc_indel = get_mat_vc_indel(test_indel,base_indel)
roc_vc_indel = roc_vc(mat_vc_indel)
roc_vc_indel = format_roc_vc(roc_vc_indel)
roc_vc_cons_indel = best_roc_vote_indel(mat_vc_indel,base_indel)
roc_vc_cons_indel_ngv3 = roc_vc_cons_indel

test_TNScope_snv = read.csv('data/NGV3/TNScope/tumor_NGv3_S1_vs_normal_NGv3_S1.TNScope.somatic.snv.normalized.vcf', sep="\t", header=F,comment.char="#")
test_TNScope_snv$ID = paste(test_TNScope_snv$V1,":",test_TNScope_snv$V2,":",test_TNScope_snv$V4,">",test_TNScope_snv$V5,sep="")
test_TNScope_snv = test_TNScope_snv[which(test_TNScope_snv$V7=="PASS"),]
res_TNScope_snv = roc(base_snp,test_TNScope_snv)
#roc(base_snp,test_deepsomatic_snv[which(test_deepsomatic_snv$V6 >= 20),])
res_TNScope_snv = format_roc(res_TNScope_snv,tools="TNScope")

test_TNScope_indel = read.csv('data/NGV3/TNScope/tumor_NGv3_S1_vs_normal_NGv3_S1.TNScope.somatic.indel.normalized.vcf', sep="\t", header=F,comment.char="#")
test_TNScope_indel$ID = paste(test_TNScope_indel$V1,":",test_TNScope_indel$V2,":",test_TNScope_indel$V4,">",test_TNScope_indel$V5,sep="")
test_TNScope_indel = test_TNScope_indel[which(test_TNScope_indel$V7=="PASS" & test_TNScope_indel$V3=="."),]
res_TNScope_indel = roc(base_indel,test_TNScope_indel)
#roc(base_snp,test_deepsomatic_snv[which(test_deepsomatic_snv$V6 >= 20),])
res_TNScope_indel = format_roc(res_TNScope_indel,tools="TNScope")

test_dragen_snv = read.csv('data/NGV3/Dragen/normal_NGv3_S1.sort.RG.dedup.recall_tumor_NGv3_S1.sort.RG.dedup.recall.hard-filtered.snv.vcf', sep="\t", header=F,comment.char="#")
test_dragen_snv$ID = paste(test_dragen_snv$V1,":",test_dragen_snv$V2,":",test_dragen_snv$V4,">",test_dragen_snv$V5,sep="")
test_dragen_snv = test_dragen_snv[which(test_dragen_snv$V7=="PASS"),]
res_dragen_snv = roc(base_snp,test_dragen_snv)
res_dragen_snv = format_roc(res_dragen_snv,tools="Dragen")

test_dragen_indel = read.csv('data/NGV3/Dragen/normal_NGv3_S1.sort.RG.dedup.recall_tumor_NGv3_S1.sort.RG.dedup.recall.hard-filtered.indel.vcf', sep="\t", header=F,comment.char="#")
test_dragen_indel$ID = paste(test_dragen_indel$V1,":",test_dragen_indel$V2,":",test_dragen_indel$V4,">",test_dragen_indel$V5,sep="")
test_dragen_indel = test_dragen_indel[which(test_dragen_indel$V7=="PASS" & test_dragen_indel$V3=="."),]
res_dragen_indel = roc(base_indel,test_dragen_indel)
#roc(base_snp,test_deepsomatic_snv[which(test_deepsomatic_snv$V6 >= 20),])
res_dragen_indel = format_roc(res_dragen_indel,tools="Dragen")

res_ngv3_snv = rbind(roc_vc_snp,get_best_cons(roc_vc_cons_snp),res_TNScope_snv,res_dragen_snv)
res_ngv3_indel = rbind(roc_vc_indel,get_best_cons(roc_vc_cons_indel),res_TNScope_indel,res_dragen_indel)

res_ngv3_snv$Categ = c(rep("Classic",13),rep("Neural Network",3),rep("Comb",12),rep("Commercial",2))
res_ngv3_indel$Categ = c(rep("Classic",10),rep("Neural Network",3),rep("Comb",9),rep("Commercial",2))
#res_ngv3_snv$Categ = c(rep("Classic",13),rep("Neural Network",3),rep("Comb",15),rep("Commercial",2))
#res_ngv3_indel$Categ = c(rep("Classic",10),rep("Neural Network",3),rep("Comb",12),rep("Commercial",2))

res_ngv3_snv$Categ2 = "Individual"
res_ngv3_indel$Categ2 = "Individual"

res_ngv3_snv$Categ2[which(res_ngv3_snv$Categ=="Comb")] = "Comb"
res_ngv3_indel$Categ2[which(res_ngv3_indel$Categ=="Comb")] = "Comb"

###############################
#Seqc2 WES HC region dataset  #
###############################

base_snp = read.csv('data/Seqc2/thruth_vcf/high-confidence_sSNV_in_HC_regions_v1.2.hg38.exomeV6.recode.vcf',sep="\t",header=F,comment.char="#")
base_indel = read.csv('data/Seqc2/thruth_vcf/high-confidence_sINDEL_in_HC_regions_v1.2.hg38.exomev6.recode.vcf',sep="\t",header=F,comment.char="#")

#Dedup + BQSR
test_snp = read.csv('data/Seqc2/SomaticAnalysis/SRR7890883_S1_vs_SRR7890874_S1.somatic.snv.merge.HCregion.recode.vcf',sep="\t",header=F,comment.char="#")
test_indel = read.csv('data/Seqc2/SomaticAnalysis/SRR7890883_S1_vs_SRR7890874_S1.somatic.indel.merge.HCregion.recode.vcf',sep="\t",header=F,comment.char="#")

#Dedup
#test_snp = read.csv('data/Seqc2/SomaticAnalysis/SRR7890883_S1_vs_SRR7890874_S1.somatic.snv.merge.HCregion.recode.vcf',sep="\t",header=F,comment.char="#")
#test_indel = read.csv('data/Seqc2/SomaticAnalysis/SRR7890883_S1_vs_SRR7890874_S1.somatic.indel.merge.HCregion.recode.vcf',sep="\t",header=F,comment.char="#")

#BQSR
#test_snp = read.csv('data/Seqc2_BQSR/SomaticAnalysis/SRR7890883_S1_vs_SRR7890874_S1.somatic.snv.merge.HCregion.recode.vcf',sep="\t",header=F,comment.char="#")
#test_indel = read.csv('data/Seqc2_BQSR/SomaticAnalysis/SRR7890883_S1_vs_SRR7890874_S1.somatic.indel.merge.HCregion.recode.vcf',sep="\t",header=F,comment.char="#")

#None
#test_snp = read.csv('data/Seqc2_sans_BQSR_Dedup/SomaticAnalysis/SRR7890883_S1_vs_SRR7890874_S1.somatic.snv.merge.HCregion.recode.vcf',sep="\t",header=F,comment.char="#")
#test_indel = read.csv('data/Seqc2_sans_BQSR_Dedup/SomaticAnalysis/SRR7890883_S1_vs_SRR7890874_S1.somatic.indel.merge.HCregion.recode.vcf',sep="\t",header=F,comment.char="#")

base_snp$ID = paste(base_snp$V1,":",base_snp$V2,":",base_snp$V4,">",base_snp$V5,sep="")
base_indel$ID = paste(base_indel$V1,":",base_indel$V2,":",base_indel$V4,">",base_indel$V5,sep="")
test_snp$ID = paste(test_snp$V1,":",test_snp$V2,":",test_snp$V4,">",test_snp$V5,sep="")
test_indel$ID = paste(test_indel$V1,":",test_indel$V2,":",test_indel$V4,">",test_indel$V5,sep="")
base_snp = get_vaf(base_snp)
base_indel$VAF=1

mat_vc_snp = get_mat_vc(test_snp,base_snp)
roc_vc_snp = roc_vc(mat_vc_snp)
roc_vc_snp = format_roc_vc(roc_vc_snp)
#roc_vc_vote_snp = roc_vote(mat_vc_snp)
roc_vc_cons_snp = best_roc_vote(mat_vc_snp,base_snp)
roc_vc_cons_snp_seqc2 = roc_vc_cons_snp

mat_vc_indel = get_mat_vc_indel(test_indel,base_indel)
roc_vc_indel = roc_vc(mat_vc_indel)
roc_vc_indel = format_roc_vc(roc_vc_indel)
#roc_vc_vote_indel = roc_vote(mat_vc_indel)
roc_vc_cons_indel_seqc2 = best_roc_vote_indel(mat_vc_indel,base_indel)

test_TNScope_snv = read.csv('data/Seqc2/TNScope/SRR7890883_S1_vs_SRR7890874_S1.TNScope.somatic.snv.normalized.recode.vcf', sep="\t", header=F,comment.char="#")
test_TNScope_snv$ID = paste(test_TNScope_snv$V1,":",test_TNScope_snv$V2,":",test_TNScope_snv$V4,">",test_TNScope_snv$V5,sep="")
test_TNScope_snv = test_TNScope_snv[which(test_TNScope_snv$V7=="PASS"),]
res_TNScope_snv = roc(base_snp,test_TNScope_snv)
#roc(base_snp,test_deepsomatic_snv[which(test_deepsomatic_snv$V6 >= 20),])
res_TNScope_snv = format_roc(res_TNScope_snv,tools="TNScope")

test_TNScope_indel = read.csv('data/Seqc2/TNScope/SRR7890883_S1_vs_SRR7890874_S1.TNScope.somatic.indel.normalized.recode.vcf', sep="\t", header=F,comment.char="#")
test_TNScope_indel$ID = paste(test_TNScope_indel$V1,":",test_TNScope_indel$V2,":",test_TNScope_indel$V4,">",test_TNScope_indel$V5,sep="")
test_TNScope_indel = test_TNScope_indel[which(test_TNScope_indel$V7=="PASS" & test_TNScope_indel$V3=="."),]
res_TNScope_indel = roc(base_indel,test_TNScope_indel)
#roc(base_snp,test_deepsomatic_snv[which(test_deepsomatic_snv$V6 >= 20),])
res_TNScope_indel = format_roc(res_TNScope_indel,tools="TNScope")

test_dragen_snv = read.csv('data/Seqc2/Dragen/SRR7890874_S1.sort.RG.dedup.recall.reorder_SRR7890883_S1.sort.RG.dedup.recall.reorder.hard-filtered.snv.recode.vcf', sep="\t", header=F,comment.char="#")
test_dragen_snv$ID = paste(test_dragen_snv$V1,":",test_dragen_snv$V2,":",test_dragen_snv$V4,">",test_dragen_snv$V5,sep="")
test_dragen_snv = test_dragen_snv[which(test_dragen_snv$V7=="PASS"),]
res_dragen_snv = roc(base_snp,test_dragen_snv)
res_dragen_snv = format_roc(res_dragen_snv,tools="Dragen")

test_dragen_indel = read.csv('data/Seqc2/Dragen/SRR7890874_S1.sort.RG.dedup.recall.reorder_SRR7890883_S1.sort.RG.dedup.recall.reorder.hard-filtered.indel.recode.vcf', sep="\t", header=F,comment.char="#")
test_dragen_indel$ID = paste(test_dragen_indel$V1,":",test_dragen_indel$V2,":",test_dragen_indel$V4,">",test_dragen_indel$V5,sep="")
test_dragen_indel = test_dragen_indel[which(test_dragen_indel$V7=="PASS" & test_dragen_indel$V3=="."),]
res_dragen_indel = roc(base_indel,test_dragen_indel)
#roc(base_snp,test_deepsomatic_snv[which(test_deepsomatic_snv$V6 >= 20),])
res_dragen_indel = format_roc(res_dragen_indel,tools="Dragen")

res_seqc2_snv = rbind(roc_vc_snp,get_best_cons(roc_vc_cons_snp),res_TNScope_snv,res_dragen_snv)
res_seqc2_indel = rbind(roc_vc_indel,get_best_cons(roc_vc_cons_indel),res_TNScope_indel,res_dragen_indel)
res_seqc2_snv$Categ = c(rep("Classic",13),rep("Neural Network",3),rep("Comb",12),rep("Commercial",2))
res_seqc2_indel$Categ = c(rep("Classic",10),rep("Neural Network",3),rep("Comb",9),rep("Commercial",2))
#res_seqc2_snv$Categ = c(rep("Classic",13),rep("Neural Network",3),rep("Comb",15),rep("Commercial",2))
#res_seqc2_indel$Categ = c(rep("Classic",10),rep("Neural Network",3),rep("Comb",12),rep("Commercial",2))

res_seqc2_snv$Categ2 = "Individual"
res_seqc2_indel$Categ2 = "Individual"

res_seqc2_snv$Categ2[which(res_seqc2_snv$Categ=="Comb")] = "Comb"
res_seqc2_indel$Categ2[which(res_seqc2_indel$Categ=="Comb")] = "Comb"

############
#PERMED 2  #
############

base = read.csv('data/Permed/thruth_vcf/permed_merge2.truth.HC.vcf',sep="\t",header=F,comment.char="#")

#Dedup + BQSR
test_snp = read.csv( 'data/Permed/SomaticAnalysis/tumor_merge2_S1_vs_normal_merge2_S1.somatic.snv.merge.recode.vcf',sep="\t",header=F,comment.char="#")
test_indel = read.csv( 'data/Permed/SomaticAnalysis/tumor_merge2_S1_vs_normal_merge2_S1.somatic.indel.merge.recode.vcf',sep="\t",header=F,comment.char="#")

#Dedup
#test_snp = read.csv( 'data/Permed_Dedup/SomaticAnalysis/tumor_merge2_S1_vs_normal_merge2_S1.somatic.snv.merge.recode.vcf',sep="\t",header=F,comment.char="#")
#test_indel = read.csv( '/home/guille/NAS_IPC/data_Guille_A/benchmark/Permed_Dedup/SomaticAnalysis/tumor_merge2_S1_vs_normal_merge2_S1.somatic.indel.merge.recode.vcf',sep="\t",header=F,comment.char="#")

#BQSR
#test_snp = read.csv( 'data/Permed_BQSR/SomaticAnalysis/tumor_merge2_S1_vs_normal_merge2_S1.somatic.snv.merge.recode.vcf',sep="\t",header=F,comment.char="#")
#test_indel = read.csv( '/home/guille/NAS_IPC/data_Guille_A/benchmark/Permed_BQSR/SomaticAnalysis/tumor_merge2_S1_vs_normal_merge2_S1.somatic.indel.merge.recode.vcf',sep="\t",header=F,comment.char="#")

#None
#test_snp = read.csv( 'data/Permed_sans_BQSR_Dedup/SomaticAnalysis/tumor_merge2_S1_vs_normal_merge2_S1.somatic.snv.merge.recode.vcf',sep="\t",header=F,comment.char="#")
#test_indel = read.csv( '/home/guille/NAS_IPC/data_Guille_A/benchmark/Permed_sans_BQSR_Dedup/SomaticAnalysis/tumor_merge2_S1_vs_normal_merge2_S1.somatic.indel.merge.recode.vcf',sep="\t",header=F,comment.char="#")
base$ID = paste(base$V1,":",base$V2,":",base$V4,">",base$V5,sep="")
base_snp = base[which(base$V4 %in% c("A","T","G","C") & base$V5 %in% c("A","T","G","C")),]
base_indel = base[-which(base$V4 %in% c("A","T","G","C") & base$V5 %in% c("A","T","G","C")),]
test_snp$ID = paste(test_snp$V1,":",test_snp$V2,":",test_snp$V4,">",test_snp$V5,sep="")
test_indel$ID = paste(test_indel$V1,":",test_indel$V2,":",test_indel$V4,">",test_indel$V5,sep="")
base_snp = get_vaf(base_snp)
mat_vc_snp = get_mat_vc(test_snp,base_snp)
roc_vc_snp = roc_vc(mat_vc_snp)
roc_vc_snp = format_roc_vc(roc_vc_snp)
#roc_vc_vote_snp = roc_vote(mat_vc_snp)
roc_vc_cons_snp = best_roc_vote(mat_vc_snp,base_snp)
roc_vc_cons_snp_permed = roc_vc_cons_snp

base_indel$VAF=1
mat_vc_indel = get_mat_vc_indel(test_indel,base_indel)
roc_vc_indel = roc_vc(mat_vc_indel)
roc_vc_indel = format_roc_vc(roc_vc_indel)
#roc_vc_vote_indel = roc_vote(mat_vc_indel)
roc_vc_cons_indel = best_roc_vote_indel(mat_vc_indel,base_indel)
roc_vc_cons_indel_permed = roc_vc_cons_indel

test_TNScope_snv = read.csv('data/Permed/TNScope/tumor_merge2_S1_vs_normal_merge2_S1.TNScope.somatic.snv.normalized.recode.vcf' , sep="\t", header=F,comment.char="#")
test_TNScope_snv$ID = paste(test_TNScope_snv$V1,":",test_TNScope_snv$V2,":",test_TNScope_snv$V4,">",test_TNScope_snv$V5,sep="")
test_TNScope_snv = test_TNScope_snv[which(test_TNScope_snv$V7=="PASS"),]
res_TNScope_snv = roc(base_snp,test_TNScope_snv)
#roc(base_snp,test_deepsomatic_snv[which(test_deepsomatic_snv$V6 >= 20),])
res_TNScope_snv = format_roc(res_TNScope_snv,tools="TNScope")

test_TNScope_indel = read.csv('data/Permed/TNScope/tumor_merge2_S1_vs_normal_merge2_S1.TNScope.somatic.indel.normalized.recode.vcf' , sep="\t", header=F,comment.char="#")
test_TNScope_indel$ID = paste(test_TNScope_indel$V1,":",test_TNScope_indel$V2,":",test_TNScope_indel$V4,">",test_TNScope_indel$V5,sep="")
test_TNScope_indel = test_TNScope_indel[which(test_TNScope_indel$V7=="PASS" & test_TNScope_indel$V3=="."),]
res_TNScope_indel = roc(base_indel,test_TNScope_indel)
#roc(base_snp,test_deepsomatic_snv[which(test_deepsomatic_snv$V6 >= 20),])
res_TNScope_indel = format_roc(res_TNScope_indel,tools="TNScope")

test_dragen_snv = read.csv('data/Permed/Dragen/normal_merge2_S1.sort.RG.dedup.recall_tumor_merge2_S1.sort.RG.dedup.recall.hard-filtered.snv.recode.vcf' , sep="\t", header=F,comment.char="#")
test_dragen_snv$ID = paste(test_dragen_snv$V1,":",test_dragen_snv$V2,":",test_dragen_snv$V4,">",test_dragen_snv$V5,sep="")
test_dragen_snv = test_dragen_snv[which(test_dragen_snv$V7=="PASS"),]
res_dragen_snv = roc(base_snp,test_dragen_snv)
res_dragen_snv = format_roc(res_dragen_snv,tools="Dragen")

test_dragen_indel = read.csv('data/Permed/Dragen/normal_merge2_S1.sort.RG.dedup.recall_tumor_merge2_S1.sort.RG.dedup.recall.hard-filtered.indel.recode.vcf' , sep="\t", header=F,comment.char="#")
test_dragen_indel$ID = paste(test_dragen_indel$V1,":",test_dragen_indel$V2,":",test_dragen_indel$V4,">",test_dragen_indel$V5,sep="")
test_dragen_indel = test_dragen_indel[which(test_dragen_indel$V7=="PASS" & test_dragen_indel$V3=="."),]
res_dragen_indel = roc(base_indel,test_dragen_indel)
#roc(base_snp,test_deepsomatic_snv[which(test_deepsomatic_snv$V6 >= 20),])
res_dragen_indel = format_roc(res_dragen_indel,tools="Dragen")

res_permed_snv = rbind(roc_vc_snp,get_best_cons(roc_vc_cons_snp),res_TNScope_snv,res_dragen_snv)
res_permed_indel = rbind(roc_vc_indel,get_best_cons(roc_vc_cons_indel),res_TNScope_indel,res_dragen_indel)

res_permed_snv$Categ = c(rep("Classic",13),rep("Neural Network",3),rep("Comb",12),rep("Commercial",2))
res_permed_indel$Categ = c(rep("Classic",10),rep("Neural Network",3),rep("Comb",9),rep("Commercial",2))

res_permed_snv$Categ2 = "Individual"
res_permed_indel$Categ2 = "Individual"

res_permed_snv$Categ2[which(res_permed_snv$Categ=="Comb")] = "Comb"
res_permed_indel$Categ2[which(res_permed_indel$Categ=="Comb")] = "Comb"

###########
#HCC1143  #
###########

base = read.csv('data/HCC1143/thruth/dkfz.somatic.snv_mnv.recode.vcf',sep="\t",header=F,comment.char="#")
base = base[which(base$V7=="PASS"),]

#Dedup + BQSR
test_snp = read.csv('data/HCC1143/SomaticAnalysis/SRR6438473_S1_vs_SRR6438475_S1.somatic.snv.merge.filtered.vcf',sep="\t",header=F,comment.char="#")

#Dedup
#test_snp = read.csv('data/HCC1143_Dedup/SomaticAnalysis/SRR6438473_S1_vs_SRR6438475_S1.somatic.snv.merge.filtered.vcf',sep="\t",header=F,comment.char="#")

#BQSR
#test_snp = read.csv('data/HCC1143_BQSR/SomaticAnalysis/SRR6438473_S1_vs_SRR6438475_S1.somatic.snv.merge.filtered.vcf',sep="\t",header=F,comment.char="#")

#None
#test_snp = read.csv('data/HCC1143_sans_BQSR_Dedup/SomaticAnalysis/SRR6438473_S1_vs_SRR6438475_S1.somatic.snv.merge.filtered.vcf',sep="\t",header=F,comment.char="#")

base$ID = paste(base$V1,":",base$V2,":",base$V4,">",base$V5,sep="")
base_snp = base[which(base$V4 %in% c("A","T","G","C") & base$V5 %in% c("A","T","G","C")),]
test_snp$ID = paste(test_snp$V1,":",test_snp$V2,":",test_snp$V4,">",test_snp$V5,sep="")
base_snp = get_vaf_hcc1143(base_snp)
mat_vc_snp = get_mat_vc(test_snp,base_snp)
roc_vc_snp = roc_vc(mat_vc_snp)
roc_vc_snp = format_roc_vc(roc_vc_snp)
roc_vc_cons_snp = best_roc_vote(mat_vc_snp,base_snp)
roc_vc_cons_snp_hcc1143 = roc_vc_cons_snp

test_TNScope_snv = read.csv('data/HCC1143/TNScope/SRR6438473_S1_vs_SRR6438475_S1.TNScope.somatic.snv.normalized.vcf', sep="\t", header=F,comment.char="#")
test_TNScope_snv$ID = paste(test_TNScope_snv$V1,":",test_TNScope_snv$V2,":",test_TNScope_snv$V4,">",test_TNScope_snv$V5,sep="")
test_TNScope_snv = test_TNScope_snv[which(test_TNScope_snv$V7=="PASS"),]
res_TNScope_snv = roc(base_snp,test_TNScope_snv)
#roc(base_snp,test_deepsomatic_snv[which(test_deepsomatic_snv$V6 >= 20),])
res_TNScope_snv = format_roc(res_TNScope_snv,tools="TNScope")

test_dragen_snv = read.csv('data/HCC1143/Dragen/SRR6438475_S1.sort.RG.dedup.recall_SRR6438473_S1.sort.RG.dedup.recall.hard-filtered.vcf', sep="\t", header=F,comment.char="#")
test_dragen_snv$ID = paste(test_dragen_snv$V1,":",test_dragen_snv$V2,":",test_dragen_snv$V4,">",test_dragen_snv$V5,sep="")
test_dragen_snv = test_dragen_snv[which(test_dragen_snv$V7=="PASS"),]
res_dragen_snv = roc(base_snp,test_dragen_snv)
res_dragen_snv = format_roc(res_dragen_snv,tools="Dragen")

res_hcc1143_snv = rbind(roc_vc_snp,get_best_cons(roc_vc_cons_snp),res_TNScope_snv,res_dragen_snv)
res_hcc1143_snv$Categ = c(rep("Classic",13),rep("Neural Network",3),rep("Comb",12),rep("Commercial",2))

res_hcc1143_snv$Categ2 = "Individual"
res_hcc1143_snv$Categ2[which(res_hcc1143_snv$Categ=="Comb")] = "Comb"

###############################
#Validation SEQC2 - 3         #
###############################

base_snp = read.csv('data/Seqc2_3/thruth_vcf/high-confidence_sSNV_in_HC_regions_v1.2.hg38.exomeV6.recode.vcf',sep="\t",header=F,comment.char="#")
base_indel = read.csv('data/Seqc2_3/thruth_vcf/high-confidence_sINDEL_in_HC_regions_v1.2.hg38.exomev6.recode.vcf',sep="\t",header=F,comment.char="#")

test_snp = read.csv('data/Seqc2_3/SomaticAnalysis/SRR7890879_S1_vs_SRR7890880_S1.somatic.snv.merge.HCregion.recode.vcf',sep="\t",header=F,comment.char="#")
test_indel = read.csv('data/Seqc2_3/SomaticAnalysis/SRR7890879_S1_vs_SRR7890880_S1.somatic.indel.merge.HCregion.recode.vcf',sep="\t",header=F,comment.char="#")

base_snp$ID = paste(base_snp$V1,":",base_snp$V2,":",base_snp$V4,">",base_snp$V5,sep="")
base_indel$ID = paste(base_indel$V1,":",base_indel$V2,":",base_indel$V4,">",base_indel$V5,sep="")
test_snp$ID = paste(test_snp$V1,":",test_snp$V2,":",test_snp$V4,">",test_snp$V5,sep="")
test_indel$ID = paste(test_indel$V1,":",test_indel$V2,":",test_indel$V4,">",test_indel$V5,sep="")
base_snp = get_vaf(base_snp)
base_indel$VAF=1

mat_vc_snp = get_mat_vc(test_snp,base_snp)
roc_vc_snp = roc_vc(mat_vc_snp)
roc_vc_snp = format_roc_vc(roc_vc_snp)
#roc_vc_vote_snp = roc_vote(mat_vc_snp)
roc_vc_cons_snp = best_roc_vote(mat_vc_snp,base_snp)
roc_vc_cons_snp_seqc2_3 = roc_vc_cons_snp

mat_vc_indel = get_mat_vc_indel(test_indel,base_indel)
roc_vc_indel = roc_vc(mat_vc_indel)
roc_vc_indel = format_roc_vc(roc_vc_indel)
#roc_vc_vote_indel = roc_vote(mat_vc_indel)
roc_vc_cons_indel = best_roc_vote_indel(mat_vc_indel,base_indel)
roc_vc_cons_indel_seqc2_3 = roc_vc_cons_indel

best_comb_snv = "Lofreq-Muse-Mutect2-SomaticSniper-Strelka-Lancet"
best_comb_indel = "Mutect2-Strelka-Pindel-Varscan2"
best_comb_snv_indel = "Mutect2-Strelka-Varscan2"
best_tradeoff_snv = "Muse-Mutect2-Strelka"
best_tradeoff_indel = "Mutect2-Strelka-Varscan2"

res_seqc2_validation_snv = rbind(roc_vc_snp,roc_vc_cons_snp_seqc2_3[which(roc_vc_cons_snp_seqc2_3$Tools==best_comb_snv & roc_vc_cons_snp_seqc2_3$Nb_Vote==3),],roc_vc_cons_snp_seqc2_3[which(roc_vc_cons_snp_seqc2_3$Tools==best_tradeoff_snv & roc_vc_cons_snp_seqc2_3$Nb_Vote==2),])
res_seqc2_validation_snv$Categ = c(rep("Classic",13),rep("Neural Network",3),rep("Best-Comb-SNV",1),rep("Comb-Cost-Effective",1))

res_seqc2_validation_indel = rbind(roc_vc_indel,roc_vc_cons_indel_seqc2_3[which(roc_vc_cons_indel_seqc2_3$Tools==best_comb_indel & roc_vc_cons_indel_seqc2_3$Nb_Vote==2),],roc_vc_cons_indel_seqc2_3[which(roc_vc_cons_indel_seqc2_3$Tools==best_tradeoff_indel & roc_vc_cons_indel_seqc2_3$Nb_Vote==2),])
res_seqc2_validation_indel$Categ = c(rep("Classic",10),rep("Neural Network",3),rep("Best-Comb-INDEL",1),rep("Comb-Cost-Effective",1))

res_seqc2_snv$Categ2 = "Individual"
res_seqc2_indel$Categ2 = "Individual"

res_seqc2_snv$Categ2[which(res_seqc2_snv$Categ=="Comb")] = "Comb"
res_seqc2_indel$Categ2[which(res_seqc2_indel$Categ=="Comb")] = "Comb"

####################
# Figures and Tabs #
####################

#Tab with all the results for SNVs variant callers
ind_individual = which(res_ngv3_snv$Categ2=="Individual")
tab_results_ind_snv = data.frame(

    Tools = res_ngv3_snv$Tools[ind_individual],
    TPR_NGV3 = res_ngv3_snv$TPR[ind_individual],
    F1_NGV3 = res_ngv3_snv$F1[ind_individual],
    PPV_NGV3 = res_ngv3_snv$PPV[ind_individual],
    Rank_NGV3 = rank(1-as.numeric(as.character(res_ngv3_snv$F1[ind_individual]))),
    TPR_SEQC2 = res_seqc2_snv$TPR[ind_individual],
    F1_SEQC2 = res_seqc2_snv$F1[ind_individual],
    PPV_SEQC2 = res_seqc2_snv$PPV[ind_individual],
    Rank_SEQC2 = rank(1-as.numeric(as.character(res_seqc2_snv$F1[ind_individual]))),
    TPR_HCC1143 = res_hcc1143_snv$TPR[ind_individual],
    F1_HCC1143 = res_hcc1143_snv$F1[ind_individual],
    PPV_HCC1143 = res_hcc1143_snv$PPV[ind_individual],
    Rank_HCC1143 = rank(1-as.numeric(as.character(res_hcc1143_snv$F1[ind_individual]))),
    TPR_PERMED = res_permed_snv$TPR[ind_individual],
    F1_PERMED = res_permed_snv$F1[ind_individual],
    PPV_PERMED = res_permed_snv$PPV[ind_individual],
    Rank_PERMED = rank(1-as.numeric(as.character(res_permed_snv$F1[ind_individual])))
)
tab_results_ind_snv$Mean_F1 = apply(tab_results_ind_snv[,c("F1_NGV3","F1_SEQC2","F1_HCC1143","F1_PERMED")],1,function(x){ mean(as.numeric(as.character(x)))  })
tab_results_ind_snv$VarRank_snv = apply(tab_results_ind_snv[,c("Rank_NGV3","Rank_SEQC2","Rank_HCC1143","Rank_PERMED")],1,function(x){ var(as.numeric(as.character(x)))  })
tab_results_ind_snv$VarF1_snv = apply(tab_results_ind_snv[,c("F1_NGV3","F1_SEQC2","F1_HCC1143","F1_PERMED")],1,function(x){ var(as.numeric(as.character(x)))  })
#tab_results_ind_snv = tab_results_ind_snv[order(tab_results_ind_snv$Mean_F1,decreasing=T),]

#Tab with all the results for indels variant callers
ind_individual = which(res_ngv3_indel$Categ2=="Individual")
tab_results_ind_indel = data.frame(

    Tools = res_ngv3_indel$Tools[ind_individual],
    TPR_NGV3 = res_ngv3_indel$TPR[ind_individual],
    F1_NGV3 = res_ngv3_indel$F1[ind_individual],
    PPV_NGV3 = res_ngv3_indel$PPV[ind_individual],
    Rank_NGV3 = rank(1-as.numeric(as.character(res_ngv3_indel$F1[ind_individual]))),
    TPR_SEQC2 = res_seqc2_indel$TPR[ind_individual],
    F1_SEQC2 = res_seqc2_indel$F1[ind_individual],
    PPV_SEQC2 = res_seqc2_indel$PPV[ind_individual],
    Rank_SEQC2 = rank(1-as.numeric(as.character(res_seqc2_indel$F1[ind_individual]))),
    TPR_PERMED = res_permed_indel$TPR[ind_individual],
    F1_PERMED = res_permed_indel$F1[ind_individual],
    PPV_PERMED = res_permed_indel$PPV[ind_individual],
    Rank_PERMED = rank(1-as.numeric(as.character(res_permed_indel$F1[ind_individual])))
)
tab_results_ind_indel$Mean_F1 = apply(tab_results_ind_indel[,c("F1_NGV3","F1_SEQC2","F1_PERMED")],1,function(x){ mean(as.numeric(as.character(x)))  })
tab_results_ind_indel$VarRank_snv = apply(tab_results_ind_indel[,c("Rank_NGV3","Rank_SEQC2","Rank_PERMED")],1,function(x){ var(as.numeric(as.character(x)))  })


#Figure 1 (SNVs)
p1_A1=plot_roc_gg(res_ngv3_snv)
p1_A2=plot_roc_gg(res_seqc2_snv)
p1_A3=plot_roc_gg(res_permed_snv)
p1_A4=plot_roc_gg(res_hcc1143_snv)
ggarrange(p1_A1,p1_A2,p1_A3,p1_A4,nrow=2,ncol=2,labels=c("NGV3-SNVs","SEQC2-SNVs","PERMED-SNVs","HCC1143-SNVs"))

#Figure 1 (SNVs) pheatmap
mat_snp = tab_results_ind_snv[,c("F1_NGV3","F1_SEQC2","F1_PERMED","F1_HCC1143","Mean_F1")]
mat_snp = apply(mat_snp,2,function(x){ as.numeric(x) })
rownames(mat_snp)=tab_results_ind_snv$Tools
mat_snp = mat_snp[order(mat_snp[,"Mean_F1"],decreasing=T),]
col_fun = colorRamp2(c(0.5, 0.75, 1), c("deepskyblue3", "white", "firebrick3"))
p1_B=pheatmap(mat_snp,cluster_rows=FALSE,cluster_cols=FALSE,display_numbers=TRUE,col=col_fun,number_format="%.3f",fontsize=16) 


#Figure2 (indels)
p2_A1=plot_roc_gg(res_ngv3_indel)
p2_A2=plot_roc_gg(res_seqc2_indel)
p2_A3=plot_roc_gg(res_permed_indel)
ggarrange(p2_A1,p2_A2,p2_A3,nrow=2,ncol=2,labels=c("NGV3-indels","SEQC2-indels","PERMED-indels"))

#Figure 2 (indels) pheatmap
mat_snp = tab_results_ind_snv[,c("F1_NGV3","F1_SEQC2","F1_PERMED","F1_HCC1143","Mean_F1")]
mat_snp = apply(mat_snp,2,function(x){ as.numeric(x) })
rownames(mat_snp)=tab_results_ind_snv$Tools
mat_snp = mat_snp[order(mat_snp[,"Mean_F1"],decreasing=T),]
col_fun = colorRamp2(c(0.5, 0.75, 1), c("deepskyblue3", "white", "firebrick3"))
p1_B=pheatmap(mat_snp,cluster_rows=FALSE,cluster_cols=FALSE,display_numbers=TRUE,col=col_fun,number_format="%.3f",fontsize=16)


#Figure 3
tools_commun = intersect(tab_results_ind_snv,tab_results_ind_indels)

var_snv = tab_results_ind_snv$VarRank_snv
names(var_snv) = tab_results_ind_snv$Tools
var_indel = tab_results_ind_indel$VarRank_snv
names(var_indel) = tab_results_ind_indel$Tools

f1_snv = tab_results_ind_snv$Mean_F1
names(f1_snv) = tab_results_ind_snv$Tools
f1_indel = tab_results_ind_indel$Mean_F1
names(f1_indel) = tab_results_ind_indel$Tools

all_tools = unique(c(tab_results_ind_snv$Tools,tab_results_ind_indel$Tools))
df_var = data.frame(Tools=all_tools,Var_SNV=var_snv[all_tools],Var_indel=var_indel[all_tools])
rownames(df_var) = df_var[,1]
df_var = df_var[,-1]
par(mar=c(10,5,5,5))
barplot(t(as.matrix(df_var)),beside=T,las=2,xlab="",ylab="Variance")


#Figure 4
#Add cpu time
dat = read.csv("data/F1_and_cpu_time.csv",sep="\t",header=T)
dat_comb = dat[which(dat$Categ=="Comb"),]
dat_comb$Total_Time=NA
for(i in 1:nrow(dat_comb)){

    comb = strsplit(dat_comb$Tools[i],"-")[[1]]
    cpu_time=0
    for(j in 1:length(comb)){

        cpu_time = cpu_time + dat$CPU_Time[which(dat$Tools==comb[j])][1]

    }
    
    dat_comb$Total_Time[i] = cpu_time

}

p4_A1=ggplot(dat[which(dat$Class=="SNV"),], aes(x=CPU_Time, y=Mean_F1, size=Memory, fill=Categ, label=Labels)) +
    geom_point(alpha=0.5, shape=21, color="black") +
    scale_size(range = c(.1, 12), name="Memory") +
    geom_text_repel(aes(x=CPU_Time, y=Mean_F1, label=Labels), data=dat[which(dat$Class=="SNV"),],inherit.aes=FALSE,size=6,max.overlaps=15) +
    theme_linedraw(base_size=24) +
    theme(legend.position="bottom") +
    ylab("Mean F1 score") +
    xlab("CPU Time")


p4_A2=ggplot(dat[which(dat$Class=="INDEL"),], aes(x=CPU_Time, y=Mean_F1, size=Memory, fill=Categ, label=Labels)) +
    geom_point(alpha=0.5, shape=21, color="black") +
    scale_size(range = c(.1, 12), name="Memory") +
    geom_text_repel(aes(x=CPU_Time, y=Mean_F1, label=Labels), data=dat[which(dat$Class=="INDEL"),],inherit.aes=FALSE,size=6,max.overlaps=15) +
    theme_linedraw(base_size=24) +
    theme(legend.position="bottom") +
    ylab("Mean F1 score") +
    xlab("CPU Time")


p4_A3=ggplot(dat[which(dat$Class=="SNV, INDEL"),], aes(x=CPU_Time, y=Mean_F1, size=Memory, fill=Categ, label=Labels)) +
    geom_point(alpha=0.5, shape=21, color="black") +
    scale_size(range = c(.1, 12), name="Memory") +
    geom_text_repel(aes(x=CPU_Time, y=Mean_F1, label=Labels), data=dat[which(dat$Class=="SNV, INDEL"),],inherit.aes=FALSE,size=6,max.overlaps=15) +
    theme_linedraw(base_size=24) +
    theme(legend.position="bottom") +
    ylab("Mean F1 score") +
    xlab("CPU Time")


roc_vc_cons_snp_ngv3 = roc_vc_cons_snp_ngv3[order(roc_vc_cons_snp_ngv3$Tools,roc_vc_cons_snp_ngv3$Nb_Vote),]
roc_vc_cons_snp_seqc2 = roc_vc_cons_snp_seqc2[order(roc_vc_cons_snp_seqc2$Tools,roc_vc_cons_snp_seqc2$Nb_Vote),]
roc_vc_cons_snp_permed = roc_vc_cons_snp_permed[order(roc_vc_cons_snp_permed$Tools,roc_vc_cons_snp_permed$Nb_Vote),]
roc_vc_cons_snp_hcc1143 = roc_vc_cons_snp_hcc1143[order(roc_vc_cons_snp_hcc1143$Tools,roc_vc_cons_snp_hcc1143$Nb_Vote),]
#roc_vc_cons_snp_seqc2_ffpe = roc_vc_cons_snp_seqc2_ffpe[order(roc_vc_cons_snp_seqc2_ffpe$Tools,roc_vc_cons_snp_seqc2_ffpe$Nb_Vote),]
cons_snp = data.frame(
Tools=roc_vc_cons_snp_ngv3[,1],
Nb_Tools=roc_vc_cons_snp_ngv3[,10],
Nb_Vote=roc_vc_cons_snp_ngv3[,2],
Rank_NVG3=roc_vc_cons_snp_ngv3[,"Rank"],
Rank_Seqc2=roc_vc_cons_snp_seqc2[,"Rank"],
Rank_permed=roc_vc_cons_snp_permed[,"Rank"],
rank_hcc1143=roc_vc_cons_snp_hcc1143[,"Rank"],
TPR_NVG3=as.numeric(roc_vc_cons_snp_ngv3[,"TPR"]),
TPR_Seqc2=as.numeric(roc_vc_cons_snp_seqc2[,"TPR"]),
TPR_permed=as.numeric(roc_vc_cons_snp_permed[,"TPR"]),
TPR_hcc1143=as.numeric(roc_vc_cons_snp_hcc1143[,"TPR"]),
PPV_NVG3=as.numeric(roc_vc_cons_snp_ngv3[,"PPV"]),
PPV_Seqc2=as.numeric(roc_vc_cons_snp_seqc2[,"PPV"]),
PPV_permed=as.numeric(roc_vc_cons_snp_permed[,"PPV"]),
PPV_hcc1143=as.numeric(roc_vc_cons_snp_hcc1143[,"PPV"]),
F1_NGV3=roc_vc_cons_snp_ngv3[,"F1"],
F1_seqc2=roc_vc_cons_snp_seqc2[,"F1"],
F1_permed=roc_vc_cons_snp_permed[,"F1"],
F1_hcc1143=roc_vc_cons_snp_hcc1143[,"F1"]
)
cons_snp$Mean_rank = apply(cons_snp[,c(4,5,6,7)],1,mean)
cons_snp$Mean_TPR = apply(cons_snp[,c(8,9,10,11)],1,mean)
cons_snp$Mean_PPV = apply(cons_snp[,c(12,13,14,15)],1,mean)
cons_snp$Mean_F1 = apply(cons_snp[,c(16,17,18,19)],1,mean)
cons_snp$VarF1 = apply(cons_snp[,c(16,17,18,19)],1,var)
head(cons_snp[order(cons_snp$Mean_F1,decreasing=T),])

#summary table_cons snp
tab_summary_cons_snp=vector()
for(i in 2:13){

    best_i = cons_snp[which(cons_snp$Nb_Tools==i),]
    best_i = best_i[order(best_i$Mean_F1,decreasing=T),]
    best_i = best_i[1,]
    tab_summary_cons_snp = rbind(tab_summary_cons_snp,best_i)
}

tab_summary_cons_snp$Short_name = sapply(tab_summary_cons_snp$Tools,format_name)

p4_B1=plot_roc_ensemble_gg(tab_summary_cons_snp)



roc_vc_cons_indel_ngv3 = roc_vc_cons_indel_ngv3[order(roc_vc_cons_indel_ngv3$Tools,roc_vc_cons_indel_ngv3$Nb_Vote),]
roc_vc_cons_indel_seqc2 = roc_vc_cons_indel_seqc2[order(roc_vc_cons_indel_seqc2$Tools,roc_vc_cons_indel_seqc2$Nb_Vote),]
roc_vc_cons_indel_permed = roc_vc_cons_indel_permed[order(roc_vc_cons_indel_permed$Tools,roc_vc_cons_indel_permed$Nb_Vote),]
#roc_vc_cons_snp_seqc2_ffpe = roc_vc_cons_snp_seqc2_ffpe[order(roc_vc_cons_snp_seqc2_ffpe$Tools,roc_vc_cons_snp_seqc2_ffpe$Nb_Vote),]
cons_indel = data.frame(
Tools=roc_vc_cons_indel_ngv3[,1],
Nb_Tools=roc_vc_cons_indel_ngv3[,10],
Nb_Vote=roc_vc_cons_indel_ngv3[,2],
Rank_NVG3=roc_vc_cons_indel_ngv3[,"Rank"],
Rank_Seqc2=roc_vc_cons_indel_seqc2[,"Rank"],
Rank_permed=roc_vc_cons_indel_permed[,"Rank"],
TPR_NVG3=as.numeric(roc_vc_cons_indel_ngv3[,"TPR"]),
TPR_Seqc2=as.numeric(roc_vc_cons_indel_seqc2[,"TPR"]),
TPR_permed=as.numeric(roc_vc_cons_indel_permed[,"TPR"]),
PPV_NVG3=as.numeric(roc_vc_cons_indel_ngv3[,"PPV"]),
PPV_Seqc2=as.numeric(roc_vc_cons_indel_seqc2[,"PPV"]),
PPV_permed=as.numeric(roc_vc_cons_indel_permed[,"PPV"]),
F1_NGV3=roc_vc_cons_indel_ngv3[,"F1"],
F1_seqc2=roc_vc_cons_indel_seqc2[,"F1"],
F1_permed=roc_vc_cons_indel_permed[,"F1"])

cons_indel$Mean_rank = apply(cons_indel[,c(4,5,6)],1,mean)
cons_indel$Mean_TPR = apply(cons_indel[,c(7,8,9)],1,mean)
cons_indel$Mean_PPV = apply(cons_indel[,c(10,11,12)],1,mean)
cons_indel$Mean_F1 = apply(cons_indel[,c(13,14,15)],1,mean)
head(cons_indel[order(cons_indel$Mean_F1,decreasing=T),])

#summary table_cons indel
tab_summary_cons_indel=vector()
for(i in 2:10){

    best_i = cons_indel[which(cons_indel$Nb_Tools==i),]
    best_i = best_i[order(best_i$Mean_F1,decreasing=T),]
    best_i = best_i[1,]
    tab_summary_cons_indel = rbind(tab_summary_cons_indel,best_i)
}

tab_summary_cons_indel$Short_name = sapply(tab_summary_cons_indel$Tools,format_name)
p4_B2=plot_roc_ensemble_gg(tab_summary_cons_indel)


cons_snp$ID = paste(cons_snp$Tools,cons_snp$Nb_Vote,sep="_")
cons_indel$ID = paste(cons_indel$Tools,cons_indel$Nb_Vote,sep="_")
commun = intersect(cons_snp$ID,cons_indel$ID)
cons_snp_commun=cons_snp[which(cons_snp$ID %in% commun),]
cons_indel_commun=cons_indel[which(cons_indel$ID %in% commun),]
cons_merge = cbind(cons_snp_commun,cons_indel_commun)
cons_merge$TPR_Global = apply(cons_merge[,c(21,42)],1,mean)
cons_merge$PPV_Global = apply(cons_merge[,c(22,43)],1,mean)
cons_merge$F1_Global = apply(cons_merge[,c(23,44)],1,mean)
head(cons_merge[order(cons_merge$F1_Global,decreasing=T),])

#summary table_cons SNVs,indels
tab_summary_cons_merge=vector()
for(i in 2:8){

    best_i = cons_merge[which(cons_merge$Nb_Tools==i),]
    best_i = best_i[order(best_i$F1_Global,decreasing=T),]
    best_i = best_i[1,]
    tab_summary_cons_merge = rbind(tab_summary_cons_merge,best_i)
}
tab_summary_cons_merge$Short_name = sapply(tab_summary_cons_merge$Tools,format_name)
p4_B3=plot_roc_cons_merge_gg(tab_summary_cons_merge[,c("TPR_Global","PPV_Global","Short_name")])

ggarrange(p4_A1,p4_A2,p4_A3,p4_B1,p4_B2,p4_B3,nrow=2,ncol=3,labels=c("SNVs","indels","Both","SNVs","indels","Both"))
