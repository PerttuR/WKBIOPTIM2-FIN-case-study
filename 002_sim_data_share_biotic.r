
# 2016-10-30: adaptacao do script para casos em que var2=""
# 2017-06-08: import script from WKSDO to WKBIOPTIM
# 2017-06-08: removed implementation of localMaxima that was giving some problems; new localMaxima2 solves problems of plateaus while ascending and plateaus in maxima 
# 2017-06-10: removed implementation of localMaxima that was giving some problems; new localMaxima2 solves problems of plateaus while ascending and plateaus in maxima 
# 2017-06-10: fixed error in sims sample size per strata
# 2017-06-11: added target_var3 [version V3]
# 2017-06-11: added sample_size dependency on number of length classes
# 2017-06-12: significant change to simulation algorithm: simplification and adaptation to two stage sampling [version V4] and remaining types of sampling [version V5]
# 2017-06-18: creations of faz_sim_sample and significant other changes and simplifications [version V7]
# 2017-06-19: reviewed code, added graphs; added median to summary
# 2017-06-19: adaptated to RDB CA
# 2017-06-21: fix one issue with mode determination during wkbioptim
# 2017-06-22: simplification, for cycle 
# 2017-06-28: ?? [version V9]
# 2017-08-XX: improved annotation and handling of results [version V10]
# 2017-08-XX: exploratory analysis of smooth and proportion for mode acceptance moved outside loop [version V10]
# 2017-08-XX: included parallel processing to speed up time ~ 5 times faster [version V11]
# 2017-09-15: open project Shrimp RCM MED [Maria Falciani] 
# 2017-09-15: added lookup conversion table for column names
# 2017-09-15: added check on var distribution of samples selected for simulation
# 2018-05-22: renamed to "002_sim_data.r" [previously "teste2_v11_share_pand_clean.r"]
# 2018-05-22: renamed var sample_id->sampId
# 2018-05-22: added check of consequences of min_n selection [now also in no of samples]
# 2018-05-22: added automatic detection of original_class_span
# 2018-05-27: added automatic determination of weight-length parameters
# 2018-05-28: added mode tests on categorical variables

		# Wishlist
			# couple sampling by weight [see fishPI algorithm...]
			# adapt to CA table
			# estratificacao desigual (e.g., 1 por class comp nos pequenos, todos os grandes)
			# improve format specification
			
			# add details of original design
			# add emdist
			# Kullback, S., Leibler, R.A. (1951). On information and sufficiency. Annals of Mathematical Statistics, 22: 79-86

# add simulation of population from stratified samples (via hans gerritsen var?)

	rm(list=ls())			

# read functions
	source("sample_level_funs1.R") # contains "expl.analysis.smooth.and.modes", sim_sample

# read packages
	library(rlist) # list.merge
	#library(reshape)
	library(parallel)		

# load input data
	load("000_Inputs\\input_data.rdata")

# read variable table
	variable_table <- read.csv2("000_Inputs\\variable_table_biotic.csv", as.is=TRUE)

	
# =================	
# Sampling design
# =================

		# set sampling design of sample data
		sampling_design <- list (stratified = FALSE, strata_var = "")
		# sampling_design <- list (stratified = TRUE, strata_var = "lenCls")
		# sampling_design <- list (stratified = TRUE, strata_var = "matStage")
		# sampling_design <- list (stratified = TRUE, strata_var = "sex")
	
	
# =================	
# Sim definitions 
# =================

	# setting of the minimum number of individuals considered representative
		min_n<-100

		table_select_samples<-table(df0$sampId)[table(df0$sampId)>=min_n]; 
		samples_to_analyze<-names(table_select_samples)		
		
		
# ======================
# Mode determination		
# ======================
	
	# wishlist:
		# include ouputs from mixdist package in LocalMaxima2
		# compare ouputs with Julia's find_mode()
		# rename "ls_auto_modes" to something more appropriate

		for (var1 in variable_table$variable)
			{
			source("sample_level_funs1.R")
			print(var1)
			if (var1==variable_table$variable[1])
				{
				ls_auto_modes<-func_detect_modes_in_samples(x = droplevels(df0[df0$sampId %in% samples_to_analyze,]), variable = var1, original_class_span = variable_table[variable_table$variable == var1, "original_class_span"], smooth_class_span = variable_table[variable_table$variable == var1, "smooth_class_span"], min_proportion_to_accept_mode = variable_table[variable_table$variable == var1, "min_proportion_to_accept_mode"])	
				} else {
						source("sample_level_funs1.R")	
						ls_auto_modes<-list.merge(ls_auto_modes, func_detect_modes_in_samples(x = droplevels(df0[df0$sampId %in% samples_to_analyze,]), variable = var1, original_class_span = variable_table[variable_table$variable == var1, "original_class_span"], smooth_class_span = variable_table[variable_table$variable == var1, "smooth_class_span"], min_proportion_to_accept_mode = variable_table[variable_table$variable == var1, "min_proportion_to_accept_mode"])	)
						#ls_auto_modes<-list.merge(ls_auto_modes, func_detect_modes_in_samples(x = droplevels(df0[df0$sampId %in% "2029_999",]), variable = var1, original_class_span = variable_table[variable_table$variable == var1, "original_class_span"], smooth_class_span = variable_table[variable_table$variable == var1, "smooth_class_span"], min_proportion_to_accept_mode = variable_table[variable_table$variable == var1, "min_proportion_to_accept_mode"])	)
						}
			}
	
	str(ls_auto_modes,3)	
	ls_auto_modes[["2001_999"]]$lenCls
	ls_auto_modes[["2001_999"]]$matStage
	ls_auto_modes[["2001_999"]]$age
	ls_auto_modes[["2001_999"]]$sex
	ls_auto_modes[["2001_999"]]$mature
	

	
# ======================
# Weight - Length Relationship		
# ======================

	coefs_weight_length<-coef(lm(log(df0$indWt)~log(df0$lenCls)))
	names(coefs_weight_length)<-c("a","b")
	# coefs_weight_length<- c(a= -7.40, b=3.04)
	
	
# =======================
# Simulations	
# =======================

		# if your data is a simple random sample of individuals from which all biological variables were sampled you will be able to test all sampling strategies
		# if your data was stratified with regards to one of the variables you have to take that into account
			# you can mantain the stratification and look into the consequences of a reduction or increase in the number of individuals samples per strata
			# you can attemp to de-stratify, creating a pseudo random sample and then test scenarios in it

	source("sample_level_funs1.R")			
														# check			

	
	# creates a storage object
	ls_DT_compiled<-sapply(samples_to_analyze, function(x) NULL)

	seed<-1
	set.seed(seed)
	ptc1<-Sys.time()	
for(sampId in samples_to_analyze)
	{

	print("==========================")
	print(sampId)
	print("==========================")
	
	# selects sample
		df1<-df0[df0$sampId == sampId,]
		# use this instead if you want to look at a specific sample
			#df1<-df0[df0$sampId == "2001_999",]

# ===============			
# Simulations of samples
# ===============

	# wishlist
		# # Sampling of different number of individuals without replacement (sample size dependent of size classes in the sample)
			# tmp.n_classes<-length(ls_auto_modes[["Length_class"]][["original_breaks"]])
			# sampling_options = list (n_sims = 100, stratified = FALSE, replacement=FALSE, sample_all_available = TRUE, sample_all_available_warning = TRUE, stages="one", samp_sizes = c(1:5*tmp.n_classes), strata_var = "none", vars_to_keep = c("Length_class", "Weight", "Age", "Sex", "Maturity_stage", "Mature"))
	
		sampling_options <- list (n_sims = 50, 
							stages="one", 																						# no of stages
								stratified = FALSE, strata_var = "", 															# stratification details
									stage1_samp_size=NA, samp_sizes = c(seq(30,300, by=30), nrow(df1)), 													# samp sizes
										replacement=FALSE, 	sample_all_available = TRUE, sample_all_available_warning = TRUE, 	# replacement options
											vars_to_keep = c(""))		
	
				ls_sims1<-faz_sim_sample(sampDes = sampling_design, sampOpt = sampling_options, df1o = df1)	
				
		
# ====================
# Building of sample statistics
# ====================


	# creates storage object
	ls_sims_stats<-lapply(sapply(as.character(sampling_options$samp_sizes), function(x) NULL), function(x) sapply(variable_table$variable, function(x) NULL)) 

	vars_numerical<-variable_table$variable[variable_table$type=="numerical"]
	vars_categorical<-variable_table$variable[variable_table$type=="categorical"]
	
	detectCores(all.tests = FALSE, logical = TRUE)
	
	# set the no of cores 
	# detectCores()
	cl <- makeCluster(5)
	#clusterExport(cl, varlist = c("make_summary_numeric","variable","df1","ls_auto_modes","coefs_weight_length","localMaxima2","localMaxima","vars_numerical","vars_categorical"))		
	clusterExport(cl, varlist = c("make_summary_numeric","df1","ls_auto_modes","coefs_weight_length","localMaxima2","localMaxima","vars_numerical","vars_categorical"))		
	
	ls_auto_modes_sample<-ls_auto_modes[[sampId]]

	
	
	for (j in 1:length(sampling_options$samp_sizes))
			{
			if(!sampling_options$stratified & (sampling_options$stages =="one" | sampling_options$stages =="two")) # not stratified, one or two stages
				{
				print(paste("Processing sample size", sampling_options$samp_sizes[j]))	
				if(sampling_options$stages == "one") w <- "1st_Stage" else w <- "2nd_Stage"
				 for (variable in vars_numerical)
					{	
					print(paste(".",variable, sep=""))
					source("sample_level_funs1.R")
					ls_sims_stats[[j]][[variable]]<-do.call("rbind",lapply(ls_sims1[[j]], function(x){make_summary_numeric(x[[w]], variable, a= coefs_weight_length[["a"]], b=coefs_weight_length[["b"]])}))
					#ls_sims_stats[[j]][[variable]]<-do.call("rbind",parLapply(cl, X = ls_sims1[[j]], function(x) make_summary_numeric(x, variable = variable, a= coefs_weight_length[["a"]], b=coefs_weight_length[["b"]])))				
					}
				 for (variable in vars_categorical)
					{	
					print(paste(".",variable, sep=""))
					source("sample_level_funs1.R")
					ls_sims_stats[[j]][[variable]]<-do.call("rbind",lapply(ls_sims1[[j]], function(x){make_summary_categorical(x[[w]], variable)}))
					#ls_sims_stats[[j]][[variable]]<-do.call("rbind",parLapply(cl, X = ls_sims1[[j]], function(x) make_summary_categorical(x[[w]], variable)))				
					}	
	
				# adds weight estimate from lenCls to all variables
					if ("lenCls" %in% vars_numerical)
						{
						for (variable in c(vars_numerical, vars_categorical)[c(vars_numerical, vars_categorical) != "lenCls"] )
							{
						ls_sims_stats[[j]][[variable]]$estim_weight<-ls_sims_stats[[j]][["lenCls"]]$estim_weight
							}
						}
			
				# ==================
				# models: 
				# ==================			
				
					print(".models")
					ls_sims_stats[[j]][["weight-length"]]<-do.call("rbind",lapply(ls_sims1[[j]], function(x, model){make_models_random(x[[w]], model="weight-length")}))
					ls_sims_stats[[j]][["sex-ratio"]]<-do.call("rbind",lapply(ls_sims1[[j]], function(x, model){make_models_random(x[[w]], model="sex-ratio")}))
					ls_sims_stats[[j]][["L50"]]<-do.call("rbind",lapply(ls_sims1[[j]], function(x, model){make_models_random(x[[w]], model="L50")}))
				}
			}
			
			#browser()	
				ls_DT_compiled[[sampId]]<-sapply(c(vars_categorical, vars_numerical), function(x) NULL)
				
				DT<-sapply(c(vars_categorical, vars_numerical), function(x) NULL)
				
				target_object = ls_sims_stats
				for (variable in c(vars_categorical, vars_numerical))
				{
				# compilation of results
				DT[[variable]]<-data.frame()
				for (i in names(target_object))
					{
					DT_sim<-data.frame(sampId = df1$sampId[1], sim = as.numeric(as.character(i)), target_object[[i]][[variable]])
					DT[[variable]]<-rbind(DT[[variable]], DT_sim)
					}
				ls_DT_compiled[[sampId]][[variable]]<-DT[[variable]]
				}
			stopCluster(cl)
}
		
str(ls_DT_compiled,3)
ls_DT_compiled[["20159001_1"]][["lenCls"]]		

head(ls_DT_compiled[["20159001_21"]][["lenCls"]])		

		
		ylimite_cv = c(0,6.5)
		ylimite_mwcv = c(0,85)
		variable = "lenCls"
		
		windows(15,7); par(mfcol=c(2,3))
		a<-	"20159001_1"
		DT2<-ls_DT_compiled[[a]][[variable]]	
		boxplot(cv~sim, data=DT2, main=paste("cv of the mean",a), xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5, ylim=ylimite_cv)
		boxplot(MWCV~sim, data=DT2, main=paste("MWCV of the sample",a), xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5, ylim=ylimite_mwcv)
		abline(h=tail(DT2)$MWCV, col="blue", lty=2, lwd=2)
		a<-	"20159001_11"
		DT2<-ls_DT_compiled[[a]][[variable]]	
		boxplot(cv~sim, data=DT2, main=paste("cv of the mean",a), xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5, ylim=ylimite_cv)
		boxplot(MWCV~sim, data=DT2, main=paste("MWCV of the sample",a), xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5, ylim=ylimite_mwcv)
		abline(h=tail(DT2)$MWCV, col="blue", lty=2, lwd=2)
		title(variable, outer=T, line=-1)
		a<-	"20159001_21"
		DT2<-ls_DT_compiled[[a]][[variable]]	
		boxplot(cv~sim, data=DT2, main=paste("cv of the mean",a), xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5, ylim=ylimite_cv)
		boxplot(MWCV~sim, data=DT2, main=paste("MWCV of the sample",a), xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5, ylim=ylimite_mwcv)
		abline(h=tail(DT2)$MWCV, col="blue", lty=2, lwd=2)
		title(variable, outer=T, line=-1)
#		
#		windows(15,7); par(mfcol=c(2,2))
#		a<-	"113_999"
#		DT2<-ls_DT_compiled[[a]][[variable]]	
#		boxplot(cv~sim, data=DT2, main=paste("cv of the mean",a), xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5, ylim=ylimite_cv)
#		boxplot(MWCV~sim, data=DT2, main=paste("MWCV of the sample",a), xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5, ylim=ylimite_mwcv)
#		abline(h=tail(DT2)$MWCV, col="blue", lty=2, lwd=2)
#		a<-	"248_999"
#		DT2<-ls_DT_compiled[[a]][[variable]]	
#		boxplot(cv~sim, data=DT2, main=paste("cv of the mean",a), xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5, ylim=ylimite_cv)
#		boxplot(MWCV~sim, data=DT2, main=paste("MWCV of the sample",a), xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5, ylim=ylimite_mwcv)
#		abline(h=tail(DT2)$MWCV, col="blue", lty=2, lwd=2)		
#		title(variable, outer=T, line=-1)
#
#
#		variable = "lenCls"
#		a<-	"20159001_21"
#		boxplot(n_modes_smooth~sim, data=DT2, main=paste("n_modes_smooth",a), xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5, ylim=ylimite_cv)
#		
#	
#	
#		mean(df0[df0$sampId=="88_999","age"])
#		median(df0[df0$sampId=="88_999","age"])
#		min(df0[df0$sampId=="88_999","age"])
#		max(df0[df0$sampId=="88_999","age"])
#	
#	str(ls_DT_compiled,1)
#	
#	
#	head(ls_DT_compiled[[2]])				
#	tail(ls_DT_compiled[[2]])				
#	head(ls_DT_compiled[[1]])				
#				
#				
#			DT2 <- ls_DT_compiled[["88_999"]][["age"]]	
#			DT2 <- ls_DT_compiled[["88_999"]][["lenCls"]]	
#			
#			# plot main outputs	(for one sample)
#			
#			# outputs one sample 1
#				windows(15,7); par(mfrow=c(2,2))
#					boxplot(mean~sim, data=DT2, main="mean of the sample", xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5)
#					plot(mean~sim, data=DT2, main="mean of the sample", xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5)
#					boxplot(cv~sim, data=DT2, main="cv of the mean", xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5)
#					boxplot(ttest_prob~sim, data=DT2, main="ttest_prob", las=2, cex.axis=1.1, cex.main=1.5)
#					abline(h=0.05, col="red", lty=2, lwd=2)
#					#title(outer=TRUE, main=paste(DT2$sampId[1]), line=-1, cex.main=1.5)			
#			
#			
#				# outputs one sample 2
#				windows(15,7); par(mfrow=c(2,2))
#					boxplot(min~sim, data=DT2, main="minimum of the sample", xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5)
#					boxplot(max~sim, data=DT2, main="maximum of the sample", xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5)
#					boxplot(median~sim, data=DT2, main="median of the sample", xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5)
#					boxplot(n_class_sampled~sim, data=DT2, main="No. length classes sampled", xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5)
#					#title(outer=TRUE, main=paste("sample",DT2$sample_id[1]), line=-1)
#		
#			# outputs one sample 3		
#				windows(15,7); par(mfrow=c(2,2))					
#					#boxplot(n_modes~sim, data=DT2, main="Number of modes (original)", xlab = "sample size")
#					boxplot(n_modes_smooth~sim, data=DT2, main="Number of modes (smooth)", xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5)
#					#boxplot(n_modes_correct~sim, data=DT2, main="Number of modes correct (original)", xlab = "sample size")
#					boxplot(ks_prob~sim, data=DT2, main="ks_prob", las=2, cex.axis=1.1, cex.main=1.5)
#					abline(h=0.05, col="red", lty=2, lwd=2)
#					boxplot(n_modes_correct_smooth~sim, data=DT2, main="Number of modes correct (smooth)", xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5)
#					boxplot(MWCV~sim, data=DT2, main="MWCV of the sample", xlab = "sample size", las=2, cex.axis=1.1, cex.main=1.5)
#					abline(h=tail(DT2)$MWCV, col="blue", lty=2, lwd=2)
#					#abline(h=20, col="red", lty=2, lwd=2)			
#			
