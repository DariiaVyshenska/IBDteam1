#bioconductor R install
#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("phyloseq")
biocLite("metagenomeSeq")
biocLite("DESeq2")
library(DESeq2)
require(phyloseq)
require(metagenomeSeq)
library(ape)
require(vegan)

#####
# you will need the following files:
# feature-taxonomy-table.csv
# sample-metadata_allFactors.tsv
# feature-taxonomy-table.csv
#####

#how to check where my work directory is (wd)
#getwd()
setwd("D:/GitHub/IBDteam1")

### functions

filterTaxaByPrevolence <- function(ps, percentSamplesPresentIn){
  prevalenceThreshold <- percentSamplesPresentIn * nsamples(ps)
  toKeep <- apply(data.frame(otu_table(ps)), 1, function(taxa) return(sum(taxa > 0) > prevalenceThreshold))
  ps_filt <- prune_taxa(toKeep, ps)
  return(ps_filt)
}

CSS_norm<-function(ps){
  library(metagenomeSeq)
  ps.metaG<-phyloseq_to_metagenomeSeq(ps)
  p_stat = cumNormStatFast(ps.metaG)
  ps.metaG = cumNorm(ps.metaG, p = p_stat)
  ps.metaG.norm <- MRcounts(ps.metaG, norm = T)
  ps_CSS<-phyloseq(otu_table(ps.metaG.norm, taxa_are_rows = T), sample_data(ps),tax_table(ps))
  return(ps_CSS)
}

cal_error_type_I_and_II <- function(ps_filt, MCSIM, treatment, control){
  #set.seed(100)
  otu_ps<-otu_table(ps_filt)
  otu_ps<-as.data.frame(otu_ps)
  totu_ps<-t(otu_ps)
  
  totu_ps_Aut<-totu_ps[sample_data(ps_filt)$disease_state == treatment, ]
  totu_ps_Aut<-as.data.frame(totu_ps_Aut)
  totu_ps_Aut<-as.matrix(totu_ps_Aut)
  
  totu_ps_Cont<-totu_ps[sample_data(ps_filt)$disease_state == control,]
  totu_ps_Cont<-as.data.frame(totu_ps_Cont)
  totu_ps_Cont<-as.matrix(totu_ps_Cont)
  
  #### Get a list of dirichlet-multinomial parameters for the data using DOM
  fit.totu_ps_Cont<-DM.MoM(totu_ps_Cont)
  fit.totu_ps_Aut<-DM.MoM(totu_ps_Aut)
  group.alpha<-rbind(fit.totu_ps_Aut$gamma, fit.totu_ps_Cont$gamma)
  
  ###Seeting up the groups
  nb_reads_totu_ps_Aut<-rowSums(totu_ps_Aut)
  nb_reads_totu_ps_Cont<-rowSums(totu_ps_Cont)
  group.Nrs <- list(nb_reads_totu_ps_Aut, nb_reads_totu_ps_Cont)
  
  #####Size and Power for the Several-Sample DM Parameter Test Comparison
  #Please set Monte-Carlo to be at least 1,000 (see HMP manual)
  type_II<-MC.Xdc.statistics(group.Nrs,  numMC = MCSIM,group.alpha) 
  type_I<-MC.Xdc.statistics(group.Nrs,  numMC =MCSIM,fit.totu_ps_Aut$gamma, type="hnull") 
  error_typeI_and_II<-list()
  error_typeI_and_II[[1]]<-type_I
  error_typeI_and_II[[2]]<-type_II
  names(error_typeI_and_II)<-c("error_typeI","error_typeII")
  return(error_typeI_and_II)
}

##### Formatting input, creating phyloseq object
## import of feature table
# need to edit the initial file with feature table: remove the first line with a # + remove the '#' of the second line
features<-read.table("feature_table.txt", row.names = 1, stringsAsFactors = F)
colnames(features)<-features[1,]
features<-features[-1,]
features_tmp<-apply(features, 2, function(x) x<-as.numeric(x))
rownames(features_tmp)<-rownames(features)
features<-features_tmp
features<-as.matrix(features)
## import of taxonomy file
taxonomy<-read.csv("feature-taxonomy-table.csv", stringsAsFactors = F, row.names = 1)
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
# import of metadata file 
metada<-read.table("sample-metadata_allFactors.tsv",stringsAsFactors = F )
colnames(metada)<-metada[1,]
metada<-metada[-1,]
rownames(metada)<-metada[,1]
metada<-metada[,-1]
##Create a phyloseq object 
#official tutorial https://joey711.github.io/phyloseq/import-data.html
str(final_taxo)
phylo_object <- phyloseq(otu_table(features,taxa_are_rows = TRUE), sample_data(metada), tax_table(final_taxo))
# remove the very shallow samples 

## filtering and normalizing phyloseq object 
sort(colSums(features)) #=> remove any sample under 5000 reads L2S385 L2S380 L2S384
phylo_object_sup5000<-prune_samples(sample_sums(phylo_object)>=5000, phylo_object)

#if you want you can Filter the taxa by prevalence 
phylo_object_sup5000_min0.03<-filterTaxaByPrevolence(phylo_object_sup5000, 0.03)

#Normalization 
CSS_norm_ps_sup5000<-CSS_norm(phylo_object_sup5000_min0.03)

## caluclating statistical power for CD vs N
# test the distribution (first - format, then - test)
otu_ps<-otu_table(CSS_norm_ps_sup5000)
otu_ps<-as.data.frame(otu_ps)
totu_ps<-t(otu_ps)
C.alpha.multinomial(totu_ps) #ok reject since zero ##pvalue = 0 => rejected => Dirichlet-Multinomial distribution

MCSIM <- 1000
CD_vs_N_error_type_I_II <- cal_error_type_I_and_II(CSS_norm_ps_sup5000, MCSIM, "C", "N")


U_vs_N_error_type_I_II <- cal_error_type_I_and_II(CSS_norm_ps_sup5000, MCSIM, "U", "N")
# weired behaviour for U vs N power calculation - type I always equals type II.
# for both CSS and by prevailance normalization - power around 0.6623377
