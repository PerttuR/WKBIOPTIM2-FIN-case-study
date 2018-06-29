# ====================	
# Exploratory Analysis 
# ====================	

load("000_Inputs\\input_data.rdata")
ls()
df0
# No of samples
sampId<-paste(df0$trpCode, df0$staNum)
cat("No. of samples:",length(unique(sampId)), "\n")

# quality checks

# Length
for(variable in c("lenCls","age","sex","matStage"))
{
cat("--------", "\n")
cat("-",variable,"-", "\n")
cat("--------", "\n")
cat("No. of NA in", variable,":",sum(is.na(df0[[variable]])), "\n")
cat("%. of NA in", variable, ":",sum(is.na(df0[[variable]]))/nrow(df0)*100, "\n")
cat("Max No. of fish in a sample:",max(table(sampId)), "\n")
cat("Min No. of fish in a sample:",min(table(sampId)), "\n")
cat("Max No. of NAs in a sample: ",max(tapply(df0[[variable]],sampId, function(x){sum(is.na(x))})), "\n")
cat("No samples with NAs: ",sum(tapply(df0[[variable]],sampId, function(x){sum(is.na(x))})>0), "\n")
cat("% samples with NAs: ",sum(tapply(df0[[variable]],sampId, function(x){sum(is.na(x))})>0)/length(unique(sampId))*100, "\n")
if(variable %in% c("lenCls","age"))
	{
	cat("Min",variable,":",min(df0[[variable]], na.rm=T), "\n")
	cat("Max",variable,":",max(df0[[variable]], na.rm=T), "\n")
	cat("width of", variable,"(as auto detected):",median(diff(sort(unique(df0[[variable]])))), "\n") # NEW
	}
}

# barplot of all variable
for(variable in c("lenCls","age","sex"))
{
if(variable %in% c("lenCls","age")) {niveis<-seq(min(df0[[variable]],na.rm=T), max(df0[[variable]],na.rm=T), by=median(diff(sort(unique(df0[[variable]])))))} else {niveis=unique(df0[[variable]])}
barplot(table(factor(df0[[variable]], levels=niveis)), las=2, cex.names=0.7)
savePlot(filename = paste("001_Exploratory_analysis\\001_Barplot_All_",variable,".png", sep=""), type = "png")
graphics.off()
}
	
# example of simulation on the entire dataset
for (variable in c("lenCls","age"))
{
if (variable %in% c("lenCls","age")) {niveis = seq(min(df0[[variable]], na.rm=T), max(df0[[variable]], na.rm=T), by=median(diff(sort(unique(df0[[variable]]))))) }
if (variable %in% c("sex")) {niveis<-c("F","M")}

	# Examples of possibilities of simulation/optimization
		# without replacement
		for( i in unique(sampId))
			{

			df2<-df0[sampId==i & !is.na(df0[[variable]]),]
			# sampling the lf with various sizes
				windows(15,7)
				par(mfrow=c(2,3))			
				sampsize<-nrow(df2)
				barplot(table(factor(df2[[variable]], levels=niveis)), las=2, cex.names=0.7, main=paste("original n (NAs excluded) =", nrow(df2)))	
				barplot(table(sample(factor(df2[[variable]], levels=niveis), size=sampsize, replace=FALSE)), las=2, cex.names=0.7,  main=paste("sampled",sampsize,"wor repl"))	
				for (j in c(200,150,100,50))
				if(sampsize>=j) 
					{
					barplot(table(sample(factor(df2[[variable]], levels=niveis), size=j, replace=FALSE)), las=2, cex.names=0.7,  main=paste("sampled",j,"wor repl"))
					} else {
							plot.new()
							}
			savePlot(filename = paste("001_Exploratory_analysis\\002_Simulation_",i,"_",variable,"_without_replacement.png", sep=""), type = "png")
			dev.off()
			}

		# with replacement	
		for( i in unique(sampId))
			{
			# Example of possibilities of simulation/optimization
			df2<-df0[sampId==i & !is.na(df0[[variable]]),]
			# sampling the lf with various sizes
				windows(15,7)
				par(mfrow=c(2,3))			
				sampsize<-nrow(df2)
				barplot(table(factor(df2[[variable]], levels=niveis)), las=2, cex.names=0.7, main=paste("original n (NAs excluded) =", nrow(df2)))	
				barplot(table(sample(factor(df2[[variable]], levels=niveis), size=sampsize, replace=TRUE)), las=2, cex.names=0.7,  main=paste("sampled",sampsize,"wr repl"))	
				for (j in c(200,150,100,50))
				if(sampsize>=j) 
					{
					barplot(table(sample(factor(df2[[variable]], levels=niveis), size=j, replace=TRUE)), las=2, cex.names=0.7,  main=paste("sampled",j,"wr repl"))
					} else {
							plot.new()
							}
			savePlot(filename = paste("001_Exploratory_analysis\\002_Simulation_",i,"_",variable,"_with_replacement.png", sep=""), type = "png")
			dev.off()
			}	
	
}	


	# =========================	
	# Select samples to analyse [adapt to your case]
	# =========================
		
		# explore tables below
		table(df0$trpCode)
		sort(table(df0$sampId))
	
		# define minimum number of individuals required for samples to be considered "representative" [adapt to your situation]
		# Note: 
			# if your protocol establishes the collection of a specific volume/weight than your sample size might be related to mean size of the individuals
			# under those circunstances it may be important to ensure that min_n is sufficiently low so that samples with larger individuals are not excluded
		# Cautionary note:
			# the exact implications of sampling with replacement from samples that have quite different sample sizes needs to be looked up
				# not sure how long into sample_size>>real_size we can go in different samples

			# set 		
			min_n<-100
			table_select_samples<-table(df0$sampId)[table(df0$sampId)>=min_n]; 
			# prints the number of samples being considered and the proportion of total samples they represent
			cat("No. selected samples:",length(table_select_samples),"\n")
			cat("% selected samples:",length(table_select_samples)/length(unique(df0$sampId))*100,"\n")
			samples_to_analyze<-names(table_select_samples)
			# check is representative of original data (if not, try adjusting min_n down and see if it improves)
				# note: you can adjust ctr_var1 and ctr_var2 to your needs
				#ctr_var1<-"quarter"
				#ctr_var2<-"foCatEu6"
			  ctr_var1<-"trpCode"
			  ctr_var2<-"trpCode"
				original_prop_samples<-round(prop.table(table(unique(df0[,c("sampId",ctr_var1,ctr_var2)])[,ctr_var1], unique(df0[,c("sampId",ctr_var1,ctr_var2)])[, ctr_var2]))*100,2)
				final_prop_samples<-round(prop.table(table(unique(df0[df0$sampId %in% samples_to_analyze,c("sampId",ctr_var1,ctr_var2)])[,ctr_var1], unique(df0[df0$sampId %in% samples_to_analyze,c("sampId",ctr_var1,ctr_var2)])[,ctr_var2]))*100,2)
				original_prop_indiv<-round(prop.table(table(df0[,ctr_var1], df0[, ctr_var2]))*100,2)
				final_prop_indiv<-round(prop.table(table(df0[df0$sampId %in% samples_to_analyze,ctr_var1], df0[df0$sampId %in% samples_to_analyze, ctr_var2]))*100,2)
				#differences in percent - if not, try adjusting min_n down and see if it improves
				cat("% dif in terms of samples: \n",original_prop_samples-final_prop_samples,"\n") 	
				cat("% dif in terms of indiv: \n",original_prop_indiv-final_prop_indiv,"\n") 	
				rm(ctr_var1, ctr_var2, original_prop_samples, final_prop_samples, original_prop_indiv, final_prop_indiv)

		
		
	# =======================
	# Set target variable, original_class_span, smooth_class_span, and threshold % for mode consideration [as proportion of individuals in sample]
	# =======================

	source("sample_level_funs1.R") # contains "expl.analysis.smooth.and.modes"

	
	variable = "lenCls"
	
	# setting of original class span for variable [default is automatic but you can choose manually]
		original_class_span <- median(diff(sort(unique(df0[[variable]])))) 	# automatic detection of the original length class
		#original_class_span <- 1 	# manual setting of interval of the original length class
		
		# exploratory analysis for determination of best smooth_class_span and min_proportion_to_accept_mode
			# outputs graphs that help determine the two parameters 
			expl.analysis.smooth.and.modes (df0 = df0, 
										samples_to_analyze = samples_to_analyze#[1:10]
										, 
										smooth_class_span = 2*original_class_span, 
										min_proportion_to_accept_mode = 0.01,
										save_plot = TRUE, 
										dir_save = "001_Exploratory_analysis\\", 
										file_root = paste("001_original_and_smoothed_data_",variable,"_", sep=""))			
				

				