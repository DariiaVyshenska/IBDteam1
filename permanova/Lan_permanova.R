biocLite("metagenomeSeq")
require(metagenomeSeq)


#if you want you can Filter the taxa by prevalence 
filterTaxaByPrevolence <- function(ps, percentSamplesPresentIn){
  prevalenceThreshold <- percentSamplesPresentIn * nsamples(ps)
  toKeep <- apply(data.frame(otu_table(ps)), 1, function(taxa) return(sum(taxa > 0) > prevalenceThreshold))
  ps_filt <- prune_taxa(toKeep, ps)
  return(ps_filt)
}

phylo_object_sup5000_min0.03<-filterTaxaByPrevolence(phylo_object_sup5000, 0.03)

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
  # take out barretts
  formula = ~ disease + disease_state + sex + IBD_degree + gastric_bx + nausea + dysphagia + reflux + GERD + reflux_GERD + dyspepsia + hematochezia + esophageal_stricture + diarrhea + abdominal_pain + abdom_RUQ_epigastric_pain
)

plot(ps_pcoa)
# run this and we get the results?
ps_pcoa


