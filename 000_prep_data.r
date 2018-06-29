#
# script to prepare data
#
rm(list=ls())

#Petri´s functions for generated CS data from IBAS IBTS surveys:




# load file
#option Nuno
#df0<-read.table("000_Inputs\\HER_SD_24_2014_subset.txt", header=TRUE, sep="\t")
#option all fin data
#df0<-read.table("000_Inputs\\CA_out.csv", header=TRUE, sep=";")
#option clened fin Herring data
df0<-read.table("000_Inputs\\CA_HER_cleaned_out.csv", header=TRUE, sep=";")
ref_tab<-read.csv2("000_Inputs\\col_names_conversion_table_FIN.csv")

colnames(df0)<-ref_tab$CA_Standard[match(tolower(colnames(df0)),tolower(ref_tab$Own_names))]

# removes columns not in accepted list
	df0[,!colnames(df0) %in% ref_tab$CA_Standard]<-NULL

# creates columns missing	
	for (i in ref_tab$CA_Standard)
		{if (!i %in% colnames(df0)) df0[i]<-"No info"} 

# Column prep [project specific]
	# tweak on Sex
		df0$sex<-as.character(df0$sex)
		df0$sex[df0$sex=="-" | is.na(df0$sex)]<-NA
		df0$sex<-factor(df0$sex, exclude=NULL)
	# Tweak on maturity
		df0$matStage<-as.character(df0$matStage)
		df0$matStage[df0$matStage=="#"]<-NA
		df0$matStage<-factor(df0$matStage)
	# creates mature
		df0$mature<-NA
		df0$mature[!is.na(df0$matStage) & df0$matStage %in% c(0,1,2)]<-0
		df0$mature[!is.na(df0$matStage) & !df0$matStage %in% c(0,1,2)]<-1
		df0$mature<-factor(df0$mature, levels=sort(unique(df0$mature)))
		df0$mature<-factor(df0$mature)

# creates sampID [adapt to your case]		
df0$sampId<-paste(df0$trpCode, df0$staNum, sep="_")
ls1<-split(df0, df0$sampId)
ls2<-lapply(ls1, function(x){x$indivId<-paste(x$sampId, 1:nrow(x),sep="_"); x})
df0<-do.call("rbind", ls2)
rownames(df0)<-NULL
		
save(df0, file="000_Inputs\\input_data.rdata")		


# create directories
	#dir.create("000_Inputs")
	#dir.create("001_Exploratory_Analysis")
	#dir.create("002_Simdata")
	#dir.create("003_Selected_Sample_Size")
