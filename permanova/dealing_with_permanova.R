#bioconductor R install
#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("phyloseq")
biocLite("metagenomeSeq")
#biocLite("DESeq2")
require(phyloseq)
require(metagenomeSeq)
library(ape)
require(vegan)


#how to check where my work directory is (wd)
#getwd()

setwd("/Users/rmf/Documents/classes/2018 04 spring/MCB 599 Collaborative Problem Solving/MB_599")

#I need to edit the feature table: remove the first line with a # + remove the '#' of the second line
features<-read.table("feature_table.txt", row.names = 1, stringsAsFactors = F)
colnames(features)<-features[1,]
features<-features[-1,]
#to check if we have a matrix (required by phyloseq) 
str(features) #calling the function "structure"
#here it's characters => we need numbers 
features_tmp<-apply(features, 2, function(x) x<-as.numeric(x))
rownames(features_tmp)<-rownames(features)
features<-features_tmp
features<-as.matrix(features)


#I need to fix the import of the taxo file
taxonomy<-read.csv("feature-taxonomy-table.csv", stringsAsFactors = F, row.names = 1)
#Need to split the string using ;
#apply: 1 is row and 2 is column
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

#metada from qiime, can be edited by adding columns 
metada<-read.table("sample-metadata_allFactors.tsv",stringsAsFactors = F )
colnames(metada)<-metada[1,]
metada<-metada[-1,]
rownames(metada)<-metada[,1]
metada<-metada[,-1]
#metada<-as.factor(metada): sometimes weird things happen and pyou need factors

#Create a phyloseq object 
#official tutorial https://joey711.github.io/phyloseq/import-data.html

str(final_taxo)
phylo_object <- phyloseq(otu_table(features,taxa_are_rows = TRUE), sample_data(metada), tax_table(final_taxo))


plot_bar(ps, fill = "Class")

#remove the very shallow samples 
#google the statics on your phyloseq object 
sort(colSums(features)) #=> remove any sample under 5000 reads L2S385 L2S380 L2S384
phylo_object_sup5000<-prune_samples(sample_sums(phylo_object)>=5000, phylo_object)

#Normalization 
CSS_norm<-function(ps){
  library(metagenomeSeq)
  ps.metaG<-phyloseq_to_metagenomeSeq(ps)
  p_stat = cumNormStatFast(ps.metaG)
  ps.metaG = cumNorm(ps.metaG, p = p_stat)
  ps.metaG.norm <- MRcounts(ps.metaG, norm = T)
  ps_CSS<-phyloseq(otu_table(ps.metaG.norm, taxa_are_rows = T), sample_data(ps),tax_table(ps))
  return(ps_CSS)
}

#instead of the CSS 
deSeqNorm <- function(ps){
  library(DESeq2)
  ps_dds <- phyloseq_to_deseq2(ps, ~ Treatment )
  ps_dds <- estimateSizeFactors(ps_dds, type = "poscounts")
  ps_dds <- estimateDispersions(ps_dds)
  abund <- getVarianceStabilizedData(ps_dds)
  abund <- abund + abs(min(abund)) #don't allow deseq to return negative counts
  ps_deSeq <- phyloseq(otu_table(abund, taxa_are_rows = T), sample_data(ps), tax_table(ps), phy_tree(ps))
  return(ps_deSeq)
}

#testDESeq<-deSeqNorm(phylo_object_sup5000)

CSS_norm_ps_sup5000<-CSS_norm(phylo_object_sup5000)
plot_bar(CSS_norm_ps_sup5000, fill = "Class")

#if you want you can Filter the taxa by prevalence 
filterTaxaByPrevolence <- function(ps, percentSamplesPresentIn){
  prevalenceThreshold <- percentSamplesPresentIn * nsamples(ps)
  toKeep <- apply(data.frame(otu_table(ps)), 1, function(taxa) return(sum(taxa > 0) > prevalenceThreshold))
  ps_filt <- prune_taxa(toKeep, ps)
  return(ps_filt)
}

phylo_object_sup5000_min0.03<-filterTaxaByPrevolence(phylo_object_sup5000, 0.03)


#PERMANOVA QUESTION
# Is there any "factors" 
pval_factors=c()
for (i in 2:length(metada)){
  cat (i,"\t")
  #tmp_map<-metada[!is.na(metada[,i]),]
  #tmp_map<-metada[!metada[,i] == "",]
  tmp_map<-metada[,i]
  tmp_map<-as.data.frame(tmp_map)
  row.names(tmp_map)<-row.names(metada)
  ps.tmp<-phylo_object_sup5000_min0.03
  sample_data(ps.tmp) <- tmp_map
  as.data.frame(sample_data(ps.tmp))
  
  tmp_nb_samples<-dim(otu_table(ps.tmp))[2]
  OTU_tables_bray <- phyloseq::distance(ps.tmp, method = "bray")
  df_metada <- data.frame(sample_data(ps.tmp ))
  colnames(df_metada)<-colnames(metada)[i]
  form1<-as.formula(paste("OTU_tables_bray",colnames(metada)[i],sep="~"))
  tmp<-adonis(form1, data = df_metada, permutations = 9999) #adonis test read about it
  tmp<-tmp$aov.tab$`Pr(>F)`[1]
  pval_factors<-c(pval_factors,tmp)}

pval_factors<-p.adjust(pval_factors, method="fdr")

ps_pcoa <- ordinate(
  physeq = phylo_object_sup5000_min0.03, 
  method = "CAP", 
  distance = "bray",
  formula = ~ disease + disease_state + sex + IBD_degree + gastric_bx + nausea + dysphagia + reflux + GERD + reflux_GERD + dyspepsia + hematochezia + Barretts + esophageal_stricture + diarrhea + abdominal_pain + abdom_RUQ_epigastric_pain
)

plot(ps_pcoa)

#the page to play with the plotting 
potential_confounding_factor<-c("disease", "disease_state","sex", "IBD_degree", "gastric_bx", "nausea" , "dysphagia" , "reflux" , "GERD" , "reflux_GERD" , "dyspepsia" , "hematochezia" , "Barretts" , "esophageal_stricture" , "diarrhea" , "abdominal_pain" , "abdom_RUQ_epigastric_pain")

#how to save and read an R file 
saveRDS(phylo_object_sup5000_min0.03, file="phylo_object_sup5000_min0.03.RData")
test<-readRDS("phylo_object_sup5000_min0.03.RData")

otu_table_from_ps<-otu_table(phylo_object_sup5000_min0.03)
write.csv(otu_table_from_ps, file="otu_table_from_ps.csv")
