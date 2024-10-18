library(stringr)
library(e1071)
library(ComplexHeatmap)
library(ggrepel)
library(future)
get_vaf = function(base){

	#vaf=apply(base,1,function(x){ as.numeric(gsub(".*bwaTVAF=([0-9]+.[0-9]+).*","\\1",x[8])) })
	vaf=apply(base,1,function(x){ as.numeric(gsub(".*VAF=([0-9]+.[0-9]+).*","\\1",x[8])) })
	base$VAF = vaf
	return(base)
}

get_vaf_hcc1143 = function(base){

	#vaf=apply(base,1,function(x){ as.numeric(gsub(".*bwaTVAF=([0-9]+.[0-9]+).*","\\1",x[8])) })
	vaf=apply(base,1,function(x){ as.numeric(gsub(".*AF=.*,([0-9]+.[0-9]+).*","\\1",x[8])) })
	base$VAF = vaf
	return(base)
}

get_mat_vc = function(test_snp,base_snp){

	variants_not_found = base_snp$ID[which(!base_snp$ID %in% test_snp$ID)]
	mat_vc = matrix(NA,ncol=16,nrow=nrow(test_snp)+length(variants_not_found))
	colnames(mat_vc) = c("FreeBayes","Lofreq","Muse","Mutect","Mutect2","SomaticSniper","Strelka","Vardict","Varscan2","Seurat","Lancet","Shimmer","Virmid","DeepSomatic","NeuSomatic","VarNet")
	rownames(mat_vc) = c(test_snp$ID,variants_not_found)
	for(i in 1:nrow(test_snp)){

		filt = strsplit(as.character(test_snp$V7[i]),"\\|")
		filt = filt[[1]]
		filt = filt[-1]
		call = strsplit(as.character(test_snp$V9[i]),":")
		call = call[[1]]
		
		if(c("LowVariantFreq") %in% filt){
			filt = filt[-which(filt == "LowVariantFreq")]		
		}
		
		
		if("ADP1" %in% call){		
			mat_vc[i,"FreeBayes"] = 0
			if("FreeBayes" %in% filt){
				mat_vc[i,"FreeBayes"] = 1
			}
		}

		if("ADLF" %in% call){		
			mat_vc[i,"Lofreq"] = 0
			if("Lofreq" %in% filt){
				mat_vc[i,"Lofreq"] = 1
			}
		}

		if("ADMC" %in% call){		
			mat_vc[i,"Mutect"] = 0
			if("Mutect" %in% filt){
				mat_vc[i,"Mutect"] = 1
			}
		}

		if("ADMU" %in% call){		
			mat_vc[i,"Muse"] = 0
			if("Muse" %in% filt){
				mat_vc[i,"Muse"] = 1
			}
		}

		if("ADM2" %in% call){		
			mat_vc[i,"Mutect2"] = 0
			if("Mutect2" %in% filt){
				mat_vc[i,"Mutect2"] = 1
			}
		}

		if("ADSK" %in% call){		
			mat_vc[i,"Strelka"] = 0
			if("Strelka" %in% filt){
				mat_vc[i,"Strelka"] = 1
			}
		}

		if("ADSN" %in% call){		
			mat_vc[i,"SomaticSniper"] = 0
			if("SomaticSniper" %in% filt){
				mat_vc[i,"SomaticSniper"] = 1
			}
		}

		if("ADVC" %in% call){		
			mat_vc[i,"Vardict"] = 0
			if("Vardict" %in% filt){
				mat_vc[i,"Vardict"] = 1
			}
		}

		if("ADVS2" %in% call){		
			mat_vc[i,"Varscan2"] = 0
			if("Varscan2" %in% filt){
				mat_vc[i,"Varscan2"] = 1
			}
		}

		if("ADSE" %in% call){		
			mat_vc[i,"Seurat"] = 0
			if("Seurat" %in% filt){
				mat_vc[i,"Seurat"] = 1
			}
		}

		if("ADLA" %in% call){		
			mat_vc[i,"Lancet"] = 0
			if("Lancet" %in% filt){
				mat_vc[i,"Lancet"] = 1
			}
		}

		if("ADSH" %in% call){		
			mat_vc[i,"Shimmer"] = 0
			if("Shimmer" %in% filt){
				mat_vc[i,"Shimmer"] = 1
			}
		}

		if("ADVI" %in% call){		
			mat_vc[i,"Virmid"] = 0
			if("Virmid" %in% filt){
				mat_vc[i,"Virmid"] = 1
			}
		}

        if("ADSS" %in% call){		
			mat_vc[i,"DeepSomatic"] = 0
			if("DeepSomatic" %in% filt){
				mat_vc[i,"DeepSomatic"] = 1
			}
		}
		
        if("ADNS" %in% call){		
			mat_vc[i,"NeuSomatic"] = 0
			if("NeuSomatic" %in% filt){
				mat_vc[i,"NeuSomatic"] = 1
			}
		}

        if("ADVN" %in% call){		
			mat_vc[i,"VarNet"] = 0
			if("VarNet" %in% filt){
				mat_vc[i,"VarNet"] = 1
			}
		}
        

	}
	mat_vc = as.data.frame(mat_vc)
	mat_vc$Vote = apply(mat_vc,1,sum,na.rm=T)
	mat_vc$Vote[which(is.na(mat_vc$Vote))] = 0
	mat_vc$Classif=0
	mat_vc$Classif[which(rownames(mat_vc) %in% base_snp$ID)] = 1
	mat_vc$VAF = NA
	ind_match = match(rownames(mat_vc),base_snp$ID)
	mat_vc$VAF[which(rownames(mat_vc) %in% base_snp$ID)] = as.numeric(as.character(base_snp$VAF[na.omit(ind_match)]))
	return(mat_vc)

}


get_mat_vc_indel = function(test_indel,base_indel){

	variants_not_found = base_indel$ID[which(!base_indel$ID %in% test_indel$ID)]
	mat_vc = matrix(NA,ncol=13,nrow=nrow(test_indel)+length(variants_not_found))
	colnames(mat_vc) = c("FreeBayes","Lofreq","Mutect2","Strelka","Vardict","Pindel","Scalpel","Varscan2","Seurat","Lancet","DeepSomatic","NeuSomatic","VarNet")
	rownames(mat_vc) = c(test_indel$ID,variants_not_found)
	
	for(i in 1:nrow(test_indel)){

		filt = strsplit(as.character(test_indel$V7[i]),"\\|")
		filt = filt[[1]]
		filt = filt[-1]
		call = strsplit(as.character(test_indel$V9[i]),":")
		call = call[[1]]
		
		if(c("LowVariantFreq") %in% filt){
			filt = filt[-which(filt == "LowVariantFreq")]		
		}
		
		
		if("ADP1" %in% call){		
			mat_vc[i,"FreeBayes"] = 0
			if("FreeBayes" %in% filt){
				mat_vc[i,"FreeBayes"] = 1
			}
		}

		if("ADLF" %in% call){		
			mat_vc[i,"Lofreq"] = 0
			if("Lofreq" %in% filt){
				mat_vc[i,"Lofreq"] = 1
			}
		}

		if("ADM2" %in% call){		
			mat_vc[i,"Mutect2"] = 0
			if("Mutect2" %in% filt){
				mat_vc[i,"Mutect2"] = 1
			}
		}

		if("ADPI" %in% call){		
			mat_vc[i,"Pindel"] = 0
			if("Pindel" %in% filt){
				mat_vc[i,"Pindel"] = 1
			}
		}

		if("ADSC" %in% call){		
			mat_vc[i,"Scalpel"] = 0
			if("Scalpel" %in% filt){
				mat_vc[i,"Scalpel"] = 1
			}
		}

		if("ADSK" %in% call){		
			mat_vc[i,"Strelka"] = 0
			if("Strelka" %in% filt){
				mat_vc[i,"Strelka"] = 1
			}
		}

		if("ADVC" %in% call){		
			mat_vc[i,"Vardict"] = 0
			if("Vardict" %in% filt){
				mat_vc[i,"Vardict"] = 1
			}
		}

		if("ADVS2" %in% call){		
			mat_vc[i,"Varscan2"] = 0
			if("Varscan2" %in% filt){
				mat_vc[i,"Varscan2"] = 1
			}
		}

		if("ADSE" %in% call){		
			mat_vc[i,"Seurat"] = 0
			if("Seurat" %in% filt){
				mat_vc[i,"Seurat"] = 1
			}
		}

		if("ADLA" %in% call){		
			mat_vc[i,"Lancet"] = 0
			if("Lancet" %in% filt){
				mat_vc[i,"Lancet"] = 1
			}
		}
		
        if("ADSS" %in% call){		
			mat_vc[i,"DeepSomatic"] = 0
			if("DeepSomatic" %in% filt){
				mat_vc[i,"DeepSomatic"] = 1
			}
		}
		
        if("ADNS" %in% call){		
			mat_vc[i,"NeuSomatic"] = 0
			if("NeuSomatic" %in% filt){
				mat_vc[i,"NeuSomatic"] = 1
			}
		}

        if("ADVN" %in% call){		
			mat_vc[i,"VarNet"] = 0
			if("VarNet" %in% filt){
				mat_vc[i,"VarNet"] = 1
			}
		}

	}
	mat_vc = as.data.frame(mat_vc)
	mat_vc$Vote = apply(mat_vc,1,sum,na.rm=T)
	mat_vc$Vote[which(is.na(mat_vc$Vote))] = 0
	mat_vc$Classif=0
	mat_vc$Classif[which(rownames(mat_vc) %in% base_indel$ID)] = 1
	mat_vc$VAF = NA
	ind_match = match(rownames(mat_vc),base_indel$ID)
	mat_vc$VAF[which(rownames(mat_vc) %in% base_indel$ID)] = as.numeric(as.character(base_indel$VAF[na.omit(ind_match)]))

	return(mat_vc)

}

roc_vote_n = function(mat_vc,comb,n){

	res = vector()
	
	for(i in 1:ncol(comb)){
		mat_vc$Vote2 = 0
		mat_vc$Vote2 = rowSums(mat_vc[,comb[,i]],na.rm=T)
		for(j in 1:n){
			mat_vc$Predict = 0
			mat_vc$Predict[which(mat_vc$Vote2 >= j)] = 1
			roc_ij = roc2(mat_vc)
			test = c(paste(comb[,i],collapse="-"),j)
			roc_ij = c(test,roc_ij)
			res = rbind(res,roc_ij)
		}
	}
	res = as.data.frame(res)
    res$Nb_Tools=n
	return(res)

}

roc_vc = function(mat_vc,vaf=1){

	variant_remove = which(mat_vc$VAF > vaf)
	cat("variants removed : ",length(variant_remove),"\n")
	if(length(variant_remove) > 0){
		mat_vc = mat_vc[-variant_remove,]	
	}
	res = vector()
	vc = colnames(mat_vc)
	vc = vc[-which(vc %in% c("Vote","Classif","Predict","VAF"))]
	for(i in 1:length(vc)){

		mat_vc$Predict = 0
		mat_vc$Predict[which(mat_vc[,vc[i]] == 1)] = 1
		roc_i = roc2(mat_vc)
		res = rbind(res,roc_i)	

	}

	res = cbind(vc,res)
	res = as.data.frame(res)
    res$Rank = rank(1-as.numeric(as.character(res$F1)))
	return(res)

}

plot_roc_gg = function(res_roc,title="",xlim=c(0,1),ylim=c(0,1),legend=TRUE){

    res_roc = res_roc[which(res_roc$Categ2=="Individual"),]
    res_roc = res_roc[order(res_roc$F1,decreasing=T),]
    res_roc = res_roc[1:10,]
    col_vc = c(
                "Mutect2"="#466791",
                "Lofreq"="#60bf37",
                "Strelka"="#953ada",
                "Muse"="#4fbe6c",
                "SomaticSniper"="#ce49d3",
                "Varscan2"="#a7b43d",
                "Mutect"="#5a51dc",
                "FreeBayes"="#d49f36",
                "Vardict"="#552095",
                "Pindel"="#507f2d",
                "Seurat"="#db37aa",
                "Lancet"="#84b67c",
                "Shimmer"="#a06fda",
                "Virmid"="#df462a",
                "Scalpel"="#5b83db",
                "Dragen"="#c76c2d",
                "TNScope"="#4f49a3",
                "DeepSomatic"="#82702d",
                "NeuSomatic"="#334c22",
                "VarNet"="#dd6bbb"              
                )    
    res_roc$PPV = as.numeric(as.character(res_roc$PPV))
    res_roc$TPR = as.numeric(as.character(res_roc$TPR))
	ggplot(res_roc, aes(PPV, TPR, colour=Tools, label=Tools)) + geom_point(size=3) + geom_text_repel(size=6, max.overlaps=15) + theme_linedraw(base_size=24) + theme(legend.position="none")+labs(y= "Mean recall (TPR)", x = "Mean precision (PPV)") + scale_color_manual(values=col_vc)
	#dev.off()
}

plot_roc_ensemble_gg = function(cons_summary,title="",xlim=c(0,1),ylim=c(0,1),legend=TRUE){

   col_vc = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")[1:nrow(cons_summary)]
	ggplot(cons_summary, aes(Mean_PPV, Mean_TPR, colour=col_vc, label=Short_name)) + geom_point(size=3) + geom_text_repel(size=6, max.overlaps=15) + theme_linedraw(base_size=24) + theme(legend.position="none") +labs(y= "Mean recall (TPR)", x = "Mean precision (PPV)")
	#dev.off()
}

plot_roc_cons_merge_gg = function(cons_summary,title="",xlim=c(0,1),ylim=c(0,1),legend=TRUE){
    
    col_vc = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")[1:nrow(cons_summary)]
	ggplot(cons_summary, aes(PPV_Global, TPR_Global, colour=col_vc, label=Short_name)) + geom_point(size=3) + geom_text_repel(size=6, max.overlaps=15) + theme_linedraw(base_size=24) + theme(legend.position="none") +labs(y= "Mean recall (TPR)", x = "Mean precision (PPV)")
	#dev.off()
}

best_roc_vote = function(mat_vc,base_snp){

    plan(multisession)

	mat_vc = as.data.frame(mat_vc)
	mat_vc$Classif = 0
	mat_vc$Classif[which(rownames(mat_vc) %in% base_snp$ID)]=1
	mat_vc$ID = rownames(mat_vc)
	
	#vc = c("FreeBayes","Lofreq","Muse","Mutect","Mutect2","SomaticSniper","Strelka","Vardict","Varscan2","Seurat","Lancet","Shimmer","Virmid","DeepSomatic","NeuSomatic","VarNet")
    vc = c("FreeBayes","Lofreq","Muse","Mutect","Mutect2","SomaticSniper","Strelka","Vardict","Varscan2","Seurat","Lancet","Shimmer","Virmid")
    comb2 = combn(vc,2)
	comb2_str = apply(comb2,2,function(x){ paste(x,collapse="-") })
	res_comb2 %<-% roc_vote_n(mat_vc,comb2,2)
	
	comb3 = combn(vc,3)
	comb3_str = apply(comb3,2,function(x){ paste(x,collapse="-") })
	res_comb3 %<-% roc_vote_n(mat_vc,comb3,3)
	
	comb4 = combn(vc,4)
	comb4_str = apply(comb4,2,function(x){ paste(x,collapse="-") })
	res_comb4 %<-% roc_vote_n(mat_vc,comb4,4)
	
    comb5 = combn(vc,5)
	comb5_str = apply(comb5,2,function(x){ paste(x,collapse="-") })
	res_comb5 %<-% roc_vote_n(mat_vc,comb5,5)

	comb6 = combn(vc,6)
	comb6_str = apply(comb6,2,function(x){ paste(x,collapse="-") })
	res_comb6 %<-% roc_vote_n(mat_vc,comb6,6)

	comb7 = combn(vc,7)
	comb7_str = apply(comb7,2,function(x){ paste(x,collapse="-") })
	res_comb7 %<-% roc_vote_n(mat_vc,comb7,7)

	comb8 = combn(vc,8)
	comb8_str = apply(comb8,2,function(x){ paste(x,collapse="-") })
	res_comb8 %<-% roc_vote_n(mat_vc,comb8,8)

	comb9 = combn(vc,9)
	comb9_str = apply(comb9,2,function(x){ paste(x,collapse="-") })
	res_comb9 %<-% roc_vote_n(mat_vc,comb9,9)

	comb10 = combn(vc,10)
	comb10_str = apply(comb10,2,function(x){ paste(x,collapse="-") })
	res_comb10 %<-% roc_vote_n(mat_vc,comb10,10)

    comb11 = combn(vc,11)
	comb11_str = apply(comb11,2,function(x){ paste(x,collapse="-") })
	res_comb11 %<-% roc_vote_n(mat_vc,comb11,11)

    comb12 = combn(vc,12)
	comb12_str = apply(comb12,2,function(x){ paste(x,collapse="-") })
	res_comb12 %<-% roc_vote_n(mat_vc,comb12,12)

    comb13 = combn(vc,13)
	comb13_str = apply(comb13,2,function(x){ paste(x,collapse="-") })
	res_comb13 %<-% roc_vote_n(mat_vc,comb13,13)

    #comb14 = combn(vc,14)
	#comb14_str = apply(comb14,2,function(x){ paste(x,collapse="-") })
	#res_comb14 %<-% roc_vote_n(mat_vc,comb14,14)

    #comb15 = combn(vc,15)
	#comb15_str = apply(comb15,2,function(x){ paste(x,collapse="-") })
	#res_comb15 %<-% roc_vote_n(mat_vc,comb15,15)

    #comb16 = combn(vc,16)
	#comb16_str = apply(comb16,2,function(x){ paste(x,collapse="-") })
	#res_comb16 %<-% roc_vote_n(mat_vc,comb16,16)
	
	res = list(res_comb2,res_comb3,res_comb4,res_comb5,res_comb6,res_comb7,res_comb8,res_comb9,res_comb10,res_comb11,res_comb12,res_comb13)
	#res = list(res_comb2,res_comb3,res_comb4,res_comb5,res_comb6,res_comb7,res_comb8,res_comb9,res_comb10,res_comb11,res_comb12,res_comb13,res_comb14,res_comb15,res_comb16)
	res = do.call(rbind,res)
	res = as.data.frame(res)
    colnames(res)[c(1,2)] = c("Tools","Nb_Vote")
	res$F1 = as.numeric(as.character(res$F1))
    res$Rank = rank(1-as.numeric(as.character(res$F1)))
	res = res[order(res$F1,decreasing=T),]
	return(res)
}

best_roc_vote_indel = function(mat_vc,base_indel){

    plan(multisession)

	mat_vc = as.data.frame(mat_vc)
	mat_vc$Classif = 0
	mat_vc$Classif[which(rownames(mat_vc) %in% base_indel$ID)]=1
	mat_vc$ID = rownames(mat_vc)
	
    vc = c("FreeBayes","Lofreq","Mutect2","Strelka","Vardict","Pindel","Scalpel","Varscan2","Seurat","Lancet")
	#vc = c("FreeBayes","Lofreq","Mutect2","Strelka","Vardict","Pindel","Scalpel","Varscan2","Seurat","Lancet","DeepSomatic","NeuSomatic","VarNet")
	comb2 = combn(vc,2)
	comb2_str = apply(comb2,2,function(x){ paste(x,collapse="-") })
	res_comb2 %<-% roc_vote_n(mat_vc,comb2,2)

	comb3 = combn(vc,3)
	comb3_str = apply(comb3,2,function(x){ paste(x,collapse="-") })
	res_comb3 %<-% roc_vote_n(mat_vc,comb3,3)
	
	comb4 = combn(vc,4)
	comb4_str = apply(comb4,2,function(x){ paste(x,collapse="-") })
	res_comb4 %<-% roc_vote_n(mat_vc,comb4,4)
	

	comb5 = combn(vc,5)
	comb5_str = apply(comb5,2,function(x){ paste(x,collapse="-") })
	res_comb5 %<-% roc_vote_n(mat_vc,comb5,5)

	comb6 = combn(vc,6)
	comb6_str = apply(comb6,2,function(x){ paste(x,collapse="-") })
	res_comb6 %<-% roc_vote_n(mat_vc,comb6,6)

	comb7 = combn(vc,7)
	comb7_str = apply(comb7,2,function(x){ paste(x,collapse="-") })
	res_comb7 %<-% roc_vote_n(mat_vc,comb7,7)

	comb7 = combn(vc,7)
	comb7_str = apply(comb7,2,function(x){ paste(x,collapse="-") })
	res_comb7 %<-% roc_vote_n(mat_vc,comb7,7)

	comb8 = combn(vc,8)
	comb8_str = apply(comb8,2,function(x){ paste(x,collapse="-") })
	res_comb8 %<-% roc_vote_n(mat_vc,comb8,8)

	comb9 = combn(vc,9)
	comb9_str = apply(comb9,2,function(x){ paste(x,collapse="-") })
	res_comb9 %<-% roc_vote_n(mat_vc,comb9,9)

	comb10 = combn(vc,10)
	comb10_str = apply(comb10,2,function(x){ paste(x,collapse="-") })
	res_comb10 %<-% roc_vote_n(mat_vc,comb10,10)

    #comb11 = combn(vc,11)
	#comb11_str = apply(comb11,2,function(x){ paste(x,collapse="-") })
	#res_comb11 %<-% roc_vote_n(mat_vc,comb11,11)

    #comb12 = combn(vc,12)
	#comb12_str = apply(comb12,2,function(x){ paste(x,collapse="-") })
	#res_comb12 %<-% roc_vote_n(mat_vc,comb12,12)

    #comb13 = combn(vc,13)
	#comb13_str = apply(comb13,2,function(x){ paste(x,collapse="-") })
	#res_comb13 %<-% roc_vote_n(mat_vc,comb13,13)
	
	res = list(res_comb2,res_comb3,res_comb4,res_comb5,res_comb6,res_comb7,res_comb8,res_comb9,res_comb10)
	#res = list(res_comb2,res_comb3,res_comb4,res_comb5,res_comb6,res_comb7,res_comb8,res_comb9,res_comb10,res_comb11,res_comb12,res_comb13)
	res = do.call(rbind,res)
	res = as.data.frame(res)
    colnames(res)[c(1,2)] = c("Tools","Nb_Vote")
	res$F1 = as.numeric(as.character(res$F1))
    res$Rank = rank(1-as.numeric(as.character(res$F1)))
	res = res[order(res$F1,decreasing=T),]
	return(res)
}


roc = function(base,test){

	region_size = 150603778

	P = nrow(base)
	N = region_size -nrow(base)	
	TP = sum(test$ID %in% base$ID)
	TN = N - nrow(test[-which(test$ID %in% base$ID),])
	FP = sum(!test$ID %in% base$ID)
	FN = sum(!base$ID %in% test$ID)

	TPR = TP/(TP+FN) #sensitivity
	FPR = FP/N #specificity
	F1 = (2*TP)/(2*TP + FP + FN)
	ACC = (TP+TN)/(TP+TN+FP+FN)
	PPV = TP/(TP+FP) #Predictive positiv value
	
	res = c(P,TP,FP,TPR,FPR,F1,ACC,PPV)
	names(res) = c("P","TP","FP","TPR","FPR","F1","ACC","PPV")

	return(res)

}

roc2 = function(mat_vc){

	region_size = 150603778
    #region_size = 40000000

	P = sum(mat_vc$Classif)
	N = region_size - P
	TP = sum(mat_vc$Predict==1 & mat_vc$Classif==1)
	#TN = N - sum(mat_vc$Predict==0 & mat_vc$Classif==0)
	FP = sum(mat_vc$Predict==1 & mat_vc$Classif==0)
	FN = sum(mat_vc$Predict==0 & mat_vc$Classif==1)

	TPR = TP/(TP+FN) #sensitivity
	FPR = FP/N #specificity
	F1 = (2*TP)/(2*TP + FP + FN)
	#ACC = (TP+TN)/(TP+TN+FP+FN)
	PPV = TP/(TP+FP) #Predictive positiv value

	res = c(P,TP,FP,TPR,FPR,F1,PPV)
	names(res) = c("P","TP","FP","TPR","FPR","F1","PPV")

	return(res)

}

format_roc = function(res_roc,tools=""){

    res=data.frame(Tools=tools, Nb_Vote=NA, P=res_roc[1], TP=res_roc[2], FP=res_roc[3], TPR=res_roc[4], FPR=res_roc[5],F1=res_roc[6],PPV=res_roc[8], Nb_Tools=1, Rank=NA)
    return(res)

}

format_roc_vc = function(res_roc_vc){

    res = data.frame(Tools=res_roc_vc[,1], Nb_Vote=NA, P=res_roc_vc[,2], TP=res_roc_vc[,3], FP=res_roc_vc[,4], TPR=res_roc_vc[,5], FPR=res_roc_vc[,6], F1=res_roc_vc[,7], PPV=res_roc_vc[,8], Nb_Tools=1, Rank=res_roc_vc[,9])
    return(res)

}

format_name = function(long_name){

    short_name = long_name

   names_sub = c(
                "Mutect2"="M2",
                "Lofreq"="LO",
                "Strelka"="SK",
                "Muse"="MU",
                "SomaticSniper"="SS",
                "Varscan2"="VC",
                "Mutect"="M1",
                "FreeBayes"="FB",
                "Vardict"="VD",
                "Pindel"="PL",
                "Seurat"="SR",
                "Lancet"="LC",
                "Shimmer"="SH",
                "Virmid"="VM",
                "Scalpel"="SC"                        
                )

    for(i in 1:length(names_sub)){

        short_name = str_replace_all(short_name,names(names_sub)[i],names_sub[i])        
    
    }

    return(short_name)

}

get_best_cons = function(res_cons){

    res_cons_split = split(res_cons,res_cons$Nb_Tools)
    res = do.call(rbind,lapply(res_cons_split,function(x){ x[1,] }))
    return(res)

}
