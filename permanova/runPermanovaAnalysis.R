#bioconductor R install
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("phyloseq")
biocLite("metagenomeSeq")
biocLite("DESeq2")
require(phyloseq)
require(metagenomeSeq)
library(ape)
require(vegan)

#getwd()
#setwd("/Users/rmf/Documents/classes/2018.04.spring/MCB.599.Collaborative.Problem.Solving/MB_599_permanova")

## feature table data wrangling
# remove the first line with a # + remove the '#' of the second line
features<-read.table("feature_table.txt", row.names = 1, stringsAsFactors = F)
colnames(features)<-features[1,]
features<-features[-1,]
# change from characters => numbers 
features_tmp<-apply(features, 2, function(x) x<-as.numeric(x))
rownames(features_tmp)<-rownames(features)
features<-features_tmp
features<-as.matrix(features)

# taxonomy table data wrangling
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

# metadata table data wrangling
#metada from qiime, can be edited by adding columns 
metada<-read.table("2018.9.27_sample-metadata_allFactors.txt",stringsAsFactors = F )
colnames(metada)<-metada[1,]
metada<-metada[-1,]
rownames(metada)<-metada[,1]
metada<-metada[,-1]
#metada<-as.factors(metada) # sometimes weird things happen and you need factors

# combine all 3 tables into a phyloseq object 
# official tutorial https://joey711.github.io/phyloseq/import-data.html
str(final_taxo)
phylo_object <- phyloseq(otu_table(features,taxa_are_rows = TRUE), sample_data(metada), tax_table(final_taxo))

# plot histogram of classes in phyloseq object
plot_bar(phylo_object, fill = "Class")

# data filtering (filtering on reads): remove the very shallow and very deep samples 
sort(colSums(features)) #=> remove any sample under 5000 reads L2S385 L2S380 L2S384, above 100,000
phylo_object_above100K<-prune_samples(sample_sums(phylo_object)<=100000, phylo_object)
phylo_object_sup5000<-prune_samples(sample_sums(phylo_object_above100K)>=5000, phylo_object_above100K)

# data normalization function
CSS_norm<-function(ps){
  library(metagenomeSeq)
  ps.metaG<-phyloseq_to_metagenomeSeq(ps)
  p_stat = cumNormStatFast(ps.metaG)
  ps.metaG = cumNorm(ps.metaG, p = p_stat)
  ps.metaG.norm <- MRcounts(ps.metaG, norm = T)
  ps_CSS<-phyloseq(otu_table(ps.metaG.norm, taxa_are_rows = T), sample_data(ps),tax_table(ps))
  return(ps_CSS)
}

# call CSS data normalization function
CSS_norm_ps_sup5000<-CSS_norm(phylo_object_sup5000)
plot_bar(CSS_norm_ps_sup5000, fill = "Class")

# create csv table for Dariia (for what??)
CSS_norm_ps_sup5000_table<-otu_table(CSS_norm_ps_sup5000)
CSS_norm_ps_sup5000_dataFrame<-as.data.frame(CSS_norm_ps_sup5000_table)
write.csv(CSS_norm_ps_sup5000_dataFrame, "CSS_norm_ps_sup5000_dataFrame.CSV",row.names = F)

# data filtering function (filter on prevalence of taxa in samples)
filterTaxaByPrevolence <- function(ps, percentSamplesPresentIn){
  prevalenceThreshold <- percentSamplesPresentIn * nsamples(ps)
  toKeep <- apply(data.frame(otu_table(ps)), 1, function(taxa) return(sum(taxa > 0) > prevalenceThreshold))
  ps_filt <- prune_taxa(toKeep, ps)
  return(ps_filt)
}

# call data filtering function (keeps only taxa present in 3% of samples or more)
# should remove taxa not present in more than 1 patient for this dataset
phylo_object_sup5000_min0.03<-filterTaxaByPrevolence(phylo_object_sup5000, 0.03)

# PERMANOVA TEST: are there any confounding factors 
pval_factors=c()
for (i in 2:length(metada)){
  cat (i,"\t")
  tmp_map<-metada[,i]
  tmp_map_factors<-as.factor(tmp_map)
  tmp_map<-as.data.frame(tmp_map)
  row.names(tmp_map)<-row.names(metada)
  ps.tmp<-phylo_object_sup5000_min0.03
  sample_data(ps.tmp) <- tmp_map
  as.data.frame(sample_data(ps.tmp))
  
  tmp_nb_samples<-dim(otu_table(ps.tmp))[2]
  OTU_tables_bray <- phyloseq::distance(ps.tmp, method = "bray")
  df_metada <- data.frame(sample_data(ps.tmp ))
  if (length(levels(df_metada$tmp_map)) > 1)
    {
    # commenting out: checking for levels earlier
    # tmp_map<-as.data.frame(tmp_map)
    # row.names(tmp_map)<-row.names(metada)
    # ps.tmp<-phylo_object_sup5000_min0.03
    # sample_data(ps.tmp) <- tmp_map
    # as.data.frame(sample_data(ps.tmp))
    # 
    # tmp_nb_samples<-dim(otu_table(ps.tmp))[2]
    # OTU_tables_bray <- phyloseq::distance(ps.tmp, method = "bray")
    # df_metada <- data.frame(sample_data(ps.tmp ))
    colnames(df_metada)<-colnames(metada)[i]
    form1<-as.formula(paste("OTU_tables_bray",colnames(metada)[i],sep="~"))
    tmp<-adonis(form1, data = df_metada, permutations = 9999) 
    tmp<-tmp$aov.tab$`Pr(>F)`[1]}
  else {tmp<-2} # assign 2 as pvalue (not possible, but will use for filtering out ones with only one factor later)
    pval_factors<-c(pval_factors,tmp)
    } 

# filtering on number of factors
# remove patient ID (no associated pval) and all factor names associated with pvals more than 1
names_pval_factors<-colnames(metada[2:length(metada)])[pval_factors<=1] # view to see which names to remove later
# remove all factors with level <= 1 (removes all pvalues with values greater than 1 (See above))
pval_factors<-pval_factors[pval_factors<=1]  # if less than or = 1 keep it

# calculate FDR-adjusted pvalue
names_pval_factors<-p.adjust(pval_factors, method="fdr")

# function to plot pcoa
ps_pcoa <- ordinate(
  physeq = phylo_object_sup5000_min0.03, 
  method = "CAP", 
  distance = "bray",
  formula = ~ disease + disease_state + sex + IBD_degree + gastric_bx + nausea + dysphagia + reflux + GERD + reflux_GERD + dyspepsia + hematochezia + esophageal_stricture + diarrhea + abdominal_pain + abdom_RUQ_epigastric_pain
)

### untested ###
# setting pcoa graph title parameters
title_prep<-c("Potential Confounding Factors")

# adding pcoa graph parameters
to_plot=list()
for (i in 1:length(metada)){
  to_plot[[i]]<-plot_ordination(
  physeq = phylo_object_sup5000_min0.03,
  ordination = ps_pcoa,
  color = metada[i],
  axes = c(1,2),
  title = title_prep[i]
  )
}

### end untested ###

# call plot the pcoa (without added parameters code, above)
plot(ps_pcoa)

# the page to play with the plotting 
# removed Barrett's as per names_pval_factors output
potential_confounding_factor<-c("disease", "disease_state","sex", "IBD_degree", "gastric_bx", "nausea" , "dysphagia" , "reflux" , "GERD" , "reflux_GERD" , "dyspepsia" , "hematochezia" , "esophageal_stricture" , "diarrhea" , "abdominal_pain" , "abdom_RUQ_epigastric_pain")

#how to save and read an R file 
saveRDS(phylo_object_sup5000_min0.03, file="phylo_object_sup5000_min0.03.RData")
test<-readRDS("phylo_object_sup5000_min0.03.RData")

otu_table_from_ps<-otu_table(phylo_object_sup5000_min0.03)
write.csv(otu_table_from_ps, file="otu_table_from_ps.csv")
