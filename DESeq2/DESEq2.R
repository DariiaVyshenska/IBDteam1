# DESeq2 Analysis of Gastric Microbiota of IBD Patients

###
#setwd ("1. Spring 2018/1. MCB 599 Prob Solving BLDS/DESeq2/")

# bioconductor R install
#source("https://bioconductor.org/biocLite.R")

# biocLite()
#biocLite("phyloseq")
#biocLite("metagenomeSeq")
#biocLite("DESeq2")
#biocLite("bit")
#install.packages("pheatmap")

require (DESeq2)
require(phyloseq)
require(metagenomeSeq)
library(ape)
require(vegan)
library (ggplot2)
library(devtools)  # Load the devtools package
library(reshape2)
library ("pheatmap")

####
# Files
features <-read.table("IBD_feature_table.txt", row.names = 1, stringsAsFactors = F)
taxonomy<-read.csv("feature-taxonomy-table.csv", stringsAsFactors = F, row.names = 1)
metada<-read.table("sample-metadata_allFactors.tsv",stringsAsFactors = F ) # confounders data

# Feature Table
# Remove the first line with a # + remove the '#' of the second line
colnames(features) <-features[1,]
features <-features[-1,]

# Check if we have a matrix (required by phyloseq) 
str(features)  # we have integers
# Make integers to numeric
features_tmp<-apply(features, 2, function(x) x<-as.numeric(x))
rownames(features_tmp)<-rownames(features)
features<-features_tmp
features<-as.matrix(features)
features

# Taxonomy table
# Split strings by ;
# Apply: 1 is row and 2 is column
taxo_fixed_lines<-apply(taxonomy, 1, function(line_of_dataframe) strsplit(line_of_dataframe,";"))

final_taxo<-matrix(ncol=10,nrow=length(taxo_fixed_lines))
final_taxo<-as.data.frame(final_taxo)
rownames(final_taxo)<-rownames(taxonomy)
colnames(final_taxo)<-1:length(colnames(final_taxo))
for (i in 1:length(taxo_fixed_lines)){
  for (j in 1:length(taxo_fixed_lines[[i]][[1]]))
  {
    final_taxo[i,][j]<-taxo_fixed_lines[[i]][[1]][j]
  }
}
final_taxo<-final_taxo[,-c(8,9,10)]
colnames(final_taxo)<-c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

final_taxo<-as.matrix(final_taxo)
str(final_taxo)

# Metadata from QIIME2 
colnames(metada)<-metada[1,]
metada<-metada[-1,]
rownames(metada)<-metada[,1]
metada<-metada[,-1]
#metada<-as.factor(metada)

#we want first column to be patient id
metada$id_L_style<-rownames(metada) # saving original rownames in case we mess up
rownames(metada)<-metada$patient_id # make rownames the patient id

#Create a phyloseq object 
#official tutorial https://joey711.github.io/phyloseq/import-data.html

phylo_object <- phyloseq(otu_table(features,taxa_are_rows = TRUE), sample_data(metada), tax_table(final_taxo))
plot_bar(phylo_object, fill = "Class") #classes of the unnormalized data

#Maude's github
#https://github.com/walllab/ASD_microbiome16s_public
#ASD_microbiome_final.R #final results used 

###Run normalization of data only
deSeqNorm <- function(ps){
  library(DESeq2)
  ps_dds <- phyloseq_to_deseq2(ps, ~ disease )
  ps_dds <- estimateSizeFactors(ps_dds, type = "poscounts")
  ps_dds <- estimateDispersions(ps_dds)
  abund <- getVarianceStabilizedData(ps_dds)
  abund <- abund + abs(min(abund)) #don't allow deseq to return negative counts
  ps_deSeq <- phyloseq(otu_table(abund, taxa_are_rows = T), sample_data(ps), tax_table(ps))
  return(ps_deSeq)
}
norm_DeSeq_ps_sup5000<-deSeqNorm(phylo_object_sup5000)
plot_bar(norm_DeSeq_ps, fill = "Class")

#this function runs normalization and DESeq2 
#DO NOT FEED the normalized one (norm_DeSeq_ps)

###Run DESeq2
# looking at disease versus non disease
runDESeq <- function(ps){
  diagdds = phyloseq_to_deseq2(ps, ~ disease) 
  diagdds <- estimateSizeFactors(diagdds, type = "poscounts")
  diagdds <- DESeq(diagdds,fitType="parametric", betaPrior = FALSE) 
  res = results(diagdds, contrast = c("disease", "D", "N"))
  res$padj[is.na(res$padj)] = 1
  sig <- res[res$padj <.05,]
  sigtab <- data.frame(cbind(sig, tax_table(ps)[rownames(sig), ]))
  return(sigtab)
}

deseq_res<-runDESeq(phylo_object)
deseq_res
write.csv(deseq_res, file="deseq_res.csv")

# looking at crohns versus non disease
runDESeq1 <- function(ps){
  diagdds = phyloseq_to_deseq2(ps, ~ disease_state) 
  diagdds <- estimateSizeFactors(diagdds, type = "poscounts")
  diagdds <- DESeq(diagdds,fitType="parametric", betaPrior = FALSE) 
  res = results(diagdds, contrast = c("disease_state", "C", "N"))
  res$padj[is.na(res$padj)] = 1
  sig <- res[res$padj <.05,]
  sigtab <- data.frame(cbind(sig, tax_table(ps)[rownames(sig), ]))
  return(sigtab)
}

deseq_res1<-runDESeq1(phylo_object)
deseq_res1
write.csv(deseq_res1, file="deseq_res1.csv")

# looking at crohns versus non disease
runDESeq_pvalue <- function(ps){
  diagdds = phyloseq_to_deseq2(ps, ~ disease_state) 
  diagdds <- estimateSizeFactors(diagdds, type = "poscounts")
  diagdds <- DESeq(diagdds,fitType="parametric", betaPrior = FALSE) 
  res = results(diagdds, contrast = c("disease_state", "C", "N"))
  res$padj[is.na(res$padj)] = 1
  sig <- res[res$padj <.1,]
  sigtab <- data.frame(cbind(sig, tax_table(ps)[rownames(sig), ]))
  return(sigtab)
}

deseq_res2<-runDESeq_pvalue(phylo_object_sup5000)
deseq_res2
write.csv(deseq_res1, file="deseq_res1.csv")


# looking at UC versus non disease
runDESeq2 <- function(ps){
  diagdds = phyloseq_to_deseq2(ps, ~ disease_state) 
  diagdds <- estimateSizeFactors(diagdds, type = "poscounts")
  diagdds <- DESeq(diagdds,fitType="parametric", betaPrior = FALSE) 
  res = results(diagdds, contrast = c("disease_state", "U", "N"))
  res$padj[is.na(res$padj)] = 1
  sig <- res[res$padj <.05,]
  sigtab <- data.frame(cbind(sig, tax_table(ps)[rownames(sig), ]))
  return(sigtab)
}

deseq_res2<-runDESeq2(phylo_object)
deseq_res2
write.csv(deseq_res2, file="deseq_res2.csv")

# looking at UC versus non disease
runDESeq3 <- function(ps){
  diagdds = phyloseq_to_deseq2(ps, ~ disease_state) 
  diagdds <- estimateSizeFactors(diagdds, type = "poscounts")
  diagdds <- DESeq(diagdds,fitType="parametric", betaPrior = FALSE) 
  res = results(diagdds, contrast = c("disease_state", "C", "U"))
  res$padj[is.na(res$padj)] = 1
  sig <- res[res$padj <.05,]
  sigtab <- data.frame(cbind(sig, tax_table(ps)[rownames(sig), ]))
  return(sigtab)
}

deseq_res3<-runDESeq3(phylo_object)
deseq_res3
write.csv(deseq_res3, file="deseq_res3.csv")


# looking at male versus female
runbysex <- function(ps){
  diagdds = phyloseq_to_deseq2(ps, ~ sex) 
  diagdds <- estimateSizeFactors(diagdds, type = "poscounts")
  diagdds <- DESeq(diagdds,fitType="parametric", betaPrior = FALSE) 
  res = results(diagdds, contrast = c("sex", "m", "f"))
  res$padj[is.na(res$padj)] = 1
  sig <- res[res$padj <.05,]
  sigtab <- data.frame(cbind(sig, tax_table(ps)[rownames(sig), ]))
  return(sigtab)
}

bysex_res<-runbysex(phylo_object)
bysex_res
write.csv(bysex_res, file="bysex_res.csv")

#code final https://github.com/walllab/ASD_microbiome16s_public/blob/master/ASD_microbiome_final.R

plot_bar(phylo_object, fill = "Class")
plot_bar(phylo_object_sup5000, fill = "Class")

#remove the very shallow samples 
sort(colSums(features)) #=> remove any sample under 5000 reads L2S385 L2S380 L2S384
phylo_object_sup5000<-prune_samples(sample_sums(phylo_object)>=5000, phylo_object)

# plot new phyloseq oject that has removed samples with reads<5000
plot_bar(phylo_object_sup5000, fill = "Class")

deseq_res_sup5000<-runDESeq(phylo_object_sup5000)
deseq_res_sup5000
write.csv(deseq_res_sup5000, file="deseq_res_sup5000")

deseq_res_sup5000_cd<-runDESeq1(phylo_object_sup5000)
deseq_res_sup5000_cd
write.csv(deseq_res_sup5000_cd, file="deseq_res_sup5000_cd.csv")

deseq_res_sup5000_uc<-runDESeq2(phylo_object_sup5000)
deseq_res_sup5000_uc
write.csv(deseq_res_sup5000_uc, file="deseq_res_sup5000_uc.csv")

deseq_res3<-runDESeq3(phylo_object_sup5000)
deseq_res3
write.csv(deseq_res3, file="deseq_res_sup5000_cd_uc.csv")



#if you want you can Filter the taxa by prevalence 
filterTaxaByPrevolence <- function(ps, percentSamplesPresentIn){
  prevalenceThreshold <- percentSamplesPresentIn * nsamples(ps)
  toKeep <- apply(data.frame(otu_table(ps)), 1, function(taxa) return(sum(taxa > 0) > prevalenceThreshold))
  ps_filt <- prune_taxa(toKeep, ps)
  return(ps_filt)
}

phylo_object_sup5000_min0.03<-filterTaxaByPrevolence(phylo_object_sup5000, 0.03)

###Declare function to make boxplots
make_boxplots <- function(rsv_table, grouping, pvals, title = "", xlab = "", ylab = ""){
  # rsv_table should be a sample x feature table
  # grouping should match the samples and describe what condition the samples have
  # extraFeatureLabeling is a vector that matched the features in the rsv_table that describes extra info you want to display, in this case genus
  #pvals <- round(pvals, digits = 2)
  pvals[pvals < .01] = "***"
  df_grouped <- data.frame(cbind(rsv_table), grouping) # Create dataframe with RSV abundances and Group
  
  colnames(df_grouped) <- c(colnames(rsv_table), "Group") # Rename the columns so we can refer to them later
  
  grouped_stacked <- melt(df_grouped, id = "Group") # Put dataframe in form that's easy to use with ggplot
  
  # Include Genus name in dataframe for graph labelling
  #match_seq_to_extraInfo <- data.frame(rsv_name =colnames(rsv_table), extraInfo = extraFeatureLabeling) # Create little mapping dataframe for rsv_names to their genuses
  match_seq_to_pval <- data.frame(rsv_name = colnames(rsv_table), pval_adj = pvals) # Create little mapping dataframe for rsv_names to their genuses
  #grouped_stacked$extraInfo <- as.character(match_seq_to_extraInfo$extraInfo[match(grouped_stacked$variable, match_seq_to_extraInfo$rsv_name)]) # assign genus to each rsv in ggplot friendly format
  grouped_stacked$pval <- as.character(match_seq_to_pval$pval_adj[match(grouped_stacked$variable, match_seq_to_pval$rsv_name)]) # assign genus to each rsv in ggplot friendly format
  
  # Plot! The function facet_wrap will break graphs up by whatever variables you put after the '~' character. In this case, we want to break it up by RSV name AND genus name
  p <- ggplot(grouped_stacked, aes(x=Group, y = value)) + geom_boxplot(aes(fill = Group)) +
    geom_jitter(aes(x = Group, y = value), position=position_jitter(0.2), cex=1.5, color="gray44") + 
    facet_wrap(~ variable  + pval, scale = "free") + labs(title = title, x = xlab, y = ylab) + scale_y_log10() + theme_minimal()
  
  print(p)
}


deseq_box_sup<- make_boxplots(t(otu_table(features,taxa_are_rows = TRUE)[rownames(deseq_res_sup5000), ]), 
              sample_data(metada)$disease, paste(deseq_res_sup5000$Genus, deseq_res_sup5000$Species), 
               "Differential Abundance between IBD and Controls, by Features")
make_boxplots(t(otu_table(features,taxa_are_rows = TRUE)[rownames(deseq_res), ]), sample_data(metada)$disease_state, paste(deseq_res$Family, deseq_res$Genus, deseq_res$Species), deseq_res$padj)

deseq_box_sup_cd<- make_boxplots(t(otu_table(features,taxa_are_rows = TRUE)[colnames(deseq_res_sup5000_cd), ]), 
                              sample_data(metada)$disease_state, deseq_res_sup5000_cd$padj)


#extract the significant features and the counts for each sample
#cd versus control
# otu_abundance <- otu_table(norm_DeSeq_ps)
# otu_abundance <- as.data.frame(otu_abundance)
# otu_abundance <- otu_abundance[row.names(otu_abundance) %in% row.names(deseq_res_sup5000_cd),]
# library(plot3D)
# par(mar = c(7,3,2,3))
# incol = jet.col(100)
# image2D(as.matrix(otu_abundance), col = incol, xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
# axis(2,at = seq(0,1,1/37),labels = colnames(otu_abundance), cex.axis = 0.5, las = 1)
# axis(1, at = seq(0,1,0.5), labels = deseq_res_sup5000_cd$Genus, cex.axis = 0.75, las = 2)

dev.off()
#extract the significant features and the counts for each sample
#IBD versus control
otu_abundance <- otu_table(norm_DeSeq_ps_sup5000)
otu_abundance <- as.data.frame(otu_abundance)
otu_abundance <- otu_abundance[row.names(otu_abundance) %in% row.names(deseq_res_sup5000),]
par(mar = c(7,3,2,3))
incol = jet.col(100)
image2D(as.matrix(otu_abundance), col = incol, xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
axis(2,at = seq(0,1,1/34),labels = colnames(otu_abundance), cex.axis = 0.5, las = 1)
axis(1, at = seq(0,1,1), labels = deseq_res_sup5000$Genus, cex.axis = 0.75, las = 2)
mtext('Control', side = 2, line = 1.5, adj = 1)
mtext('CD', side = 2, line = 1.5, adj = 0)
par(xpd=TRUE,oma=c(0,0,0,6.15))
abline(h = 15.5/34, lwd = 3)
