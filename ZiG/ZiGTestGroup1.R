library(metagenomeSeq)
library(readr)

# Pull in the feature table
IBD = load_meta(file = "~/Documents/R/IBDData/IBD_feature_table.txt")
# Pull in the taxonomy table for the feature table
taxaIBD = read.csv(file = "~/Documents/R/IBDData/feature-taxonomy-table.csv", stringsAsFactors = F)
# name the rows with the OTU (required for package)
j = taxaIBD[,1]
rownames(taxaIBD) = j
# Reorder the rows to match the feature table (required for package)
ord1 = match(rownames(IBD$counts),rownames(taxaIBD))
taxaIBD = taxaIBD[ord1, ]
# Pull in the phenotype table, this is a little different than the vignette because I am reading it straight from the CSV
clinIBD = read.csv(file = "~/Documents/R/IBDData/Pheno_table.csv", sep = '\t', row.names = 1, stringsAsFactors = F, na.strings=c("","-","NA"))
# Match the order of the phenotype table with the feature table
ord = match(colnames(IBD$counts),rownames(clinIBD))
clinIBD = clinIBD[ord, ]
# create an annotated dataframe of the phenotype data for input into the MRExperiment function later
phenDataIBD = AnnotatedDataFrame(clinIBD)
# create an annotated dataframe of the OTU data for input into the MRExperiment function
OTUIBD = AnnotatedDataFrame(taxaIBD)


#Create the experiment object
objIBD = newMRexperiment(IBD$counts,phenoData = phenDataIBD, featureData = OTUIBD)

# # Normalization
# q = metagenomeSeq::MRcounts(objIBD)
# # still trying to figure out normalization factor
# 
# # 
# objIBD1 = metagenomeSeq::cumNorm(objIBD,p = 0.5)
# 
# 
# rareFeat = which(rowSums(MRcounts(objIBD)>0)<10)
# objIBDTrim = objIBD[-rareFeat]
# ibdp = metagenomeSeq::cumNormStat(objIBDTrim, pFlag = T, main = "trimmed IBD data")
# 
# objIBDTrim = metagenomeSeq::cumNorm(objIBDTrim, p = ibdp)
# IBDStat = pData(objIBDTrim)$assignment
# 
# normFactor = metagenomeSeq::normFactors(objIBDTrim)
# normFactor = log2(normFactor/median(normFactor) + 1)
# mod = model.matrix(~IBDStat+normFactor)
# settings = metagenomeSeq::zigControl(maxit = 10, verbose = T)
# fit = metagenomeSeq::fitZig(obj = objIBDTrim, mod = mod, useCSSoffset = F, control = settings)


#######################
# Taking out the UC patients and only comparing Normal vs CD
uc = controls = grep("UC", pData(objIBD)$assignment)
objIBDTrim = objIBD[,-uc]
# Trim out any rare features with less than 5 samples
rareFeat = which(rowSums(MRcounts(objIBD)>0)<5)
objIBDTrim = objIBDTrim[-rareFeat,]
# Calculate proper percentile at which to normalize counts
ibdp = cumNormStat(objIBDTrim, pFlag = T, main = "Trimmed IBD data")
# Calculate scaling factors using the percentile calculated in the previous step
objIBDTrim = cumNorm(objIBDTrim, p = ibdp)
# Identify which samples have which disease status
IBDStat = pData(objIBDTrim)$assignment
# Identify which samples have which sex
sex = pData(objIBDTrim)$sex

#Calculate normalization factors
normFactor = normFactors(objIBDTrim)
normFactor = log2(normFactor/median(normFactor) + 1)
# Create the model matrix that will be used in the fit
mod = model.matrix(~IBDStat+normFactor)
# Create the settings that will be put into the ZiG fit
settings = zigControl(maxit = 10, verbose = T)
# Calculate the fit and put it in an object
fit1 = fitZig(obj = objIBDTrim, mod = mod, useCSSoffset = F, control = settings)
# Create a table with the fit and the data that we have from the ZiG
res_fit<-MRtable(fit1, number = length(fit1$taxa))
# pull out only those with p-values low enough
res_fit<-res_fit[res_fit$adjPvalues<0.05,]
# Calculate the number of samples needed to make this method effective and ammend that to the table we created
Min_effec_samp<-calculateEffectiveSamples(fit1)
Min_effec_samp<-Min_effec_samp[ names(Min_effec_samp)  %in% rownames(res_fit)]
res_fit$Min_sample<-Min_effec_samp
# Only pull out those features which have enough samples to be an effective measure
res_fit<-res_fit[res_fit$`+samples in group 0` >= Min_effec_samp & res_fit$`+samples in group 1` >= Min_effec_samp,]
# Return that object
return(res_fit)
# Return an object that shows which samples and taxon the differentially expressed features come from
library(plot3D)
# Set colors for heatmap
incol = c(jet.col(100))
# Pull out the predicted identities of the features of interest
taxons = taxaIBD[taxaIBD$Feature.ID %in% rownames(res_fit),]
# Creat the object that will turn into a heat map.
heater <- IBD$counts[row.names(IBD$counts) %in% taxons$Feature.ID,]
# Set the rownames to be the taxa instead of the feature
row.names(heater) = c("Veillonella dispar", "Prevotella melaninogenica 1", "Prevotella melaninogenica 2", "Streptococcus", "Prevotella", "Veillonella parvula", "Prevotella nigrescens")
#adjust the order to make sure they are grouped together
ord2 = match(row.names(clinIBD), colnames(heater))

dev.off()
# Set correct margins
par(mar = c(7,3,2,5))
# Use the image2D function from plot3D package to make a heat map
image2D(as.matrix(log10(heater+1)), yaxt = 'n', xaxt = 'n', xlab = '', ylab = '')
# Set the axis names to match the patient numbers and the predicted IDs
axis(2,at = seq(0,1,1/37),labels = colnames(heater), cex.axis = 0.5, las = 1)
axis(1, at = seq(0,1,1/6), labels = rownames(heater), cex.axis = 0.5, las = 2)
#heatmap(log10(as.matrix(heater+1)),col = incol, Rowv = NA, Colv = NA, cexRow = 0.75)
#which(as.matrix(heater)==0)

