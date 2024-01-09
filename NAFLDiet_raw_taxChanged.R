source("NAFLDiet_sample.R")
library(DESeq2)

#### Trim the Number of Zeros in more than 10%
0.1*381
0.2*381
0.5*381

raw_abundance_MGS_data <- read_xlsx("all-tables.xlsx",
                                         sheet = "MGS-downsized-counts")

raw_abundance_MGS_data_127 <- raw_abundance_MGS_data %>% dplyr::select(rownames(sam127)) 

raw_abundance_MGS_data_127 %>% colSums()



sam1 <- sample_data(sam127)
tax1 <- tax_table(raw_tax_matrix_127)
ntaxa(tax1)

### transformation ###


otu.raw <- otu_table(raw_abundance_MGS_data_127, taxa_are_rows = TRUE)


random_tree = rtree(ntaxa(tax1), rooted=TRUE, tip.label=taxa_names(tax1))
ex1b.raw <- phyloseq(otu.raw, sam1, tax1, random_tree)

tax_table(ex1b.raw)

ex1b.raw %>% otu_table() %>% colSums()

flist<- filterfun(kOverA(38.1, 0))



### Use relative_abundance
microbiome::transform(re_ab_127, 'log10p') -> re_ab_127.log10p
microbiome::transform(re_ab_127, 'clr') -> re_ab_127.clr
microbiome::transform(re_ab_127, 'alr', shift=10, reference = 1) -> re_ab_127.alr


otu.clr <- otu_table(re_ab_127.clr, taxa_are_rows = TRUE)
otu.log10p <- otu_table(re_ab_127.log10p, taxa_are_rows = TRUE)
otu.alr <- otu_table(re_ab_127.alr, taxa_are_rows = TRUE)

ex1b.log10p <- phyloseq(otu.log10p, sam1, tax1, random_tree)
ex1b.clr <- phyloseq(otu.clr, sam1, tax1, random_tree)
ex1b.alr <- phyloseq(otu.alr, sam1, tax1, random_tree)


### Use raw
ent.trim <- filter_taxa(ex1b.raw, flist, TRUE)
ent.trim <- filter_taxa(ex1b.log10p, flist, TRUE)
ent.trim <- filter_taxa(ex1b.clr, flist, TRUE)
ent.trim <- filter_taxa(ex1b.alr, flist, TRUE)

## filter unclassified ones 
ent.trim <- subset_taxa(ent.trim, ! phylum %in% c("unclassified"))
ent.trim

### which is then filtered such that
### only OTUs with a mean greater than 10^-5 are kept.

ent.trim <- filter_taxa(ent.trim, function(x) mean(x) > 1e-5, TRUE)
ent.trim
####### at Phylum level
#In many biological settings, the set of all organisms from all samples are well-represented in the available taxonomic reference database. When (and only when) this is the case, it is reasonable or even advisable to filter taxonomic features for which a high-rank taxonomy could not be assigned. Such ambiguous features in this setting are almost always sequence artifacts that don’t exist in nature. It should be obvious that such a filter is not appropriate for samples from poorly characterized or novel specimens, at least until the possibility of taxonomic novelty can be satisfactorily rejected. Phylum is a useful taxonomic rank to consider using for this purpose, but others may work effectively for your data.

table(tax_table(ent.trim)[,"phylum"],exclude=NULL)


dds <- phyloseq_to_deseq2(ent.trim, design = ~ Month) 
dds = DESeq(dds, test="Wald", fitType = "local", betaPrior =  T)

library(MASS)

re_ab_127_deltaM120.fatFrac %>% as.data.frame() %>%
  dplyr::select(c(contains("00925"), contains("fat")))  
cor() %>% reshape2::melt() %>% dplyr::filter(! value == "1") %>%
  arrange(desc(abs(value)))

ml <- glm.nb()



ent.trim.G1 <- subset_samples(ent.trim, Group == 1)
ent.trim.G2 <- subset_samples(ent.trim, Group == 2)
ent.trim.G3 <- subset_samples(ent.trim, Group == 3)



###############################
#############
############ adjust for age, sex, BMI, and diabetes
#############
#############
#################################################

sample_data(ent.trim)
dds1 <- phyloseq_to_deseq2(ent.trim.G1, design = ~ Month + age + gender + BMI + Diabetes + Metformin) 
dds2 <- phyloseq_to_deseq2(ent.trim.G2, design = ~ Month + age + gender + BMI + Diabetes + Metformin) 
dds3 <- phyloseq_to_deseq2(ent.trim.G3, design = ~ Month + age + gender + BMI + Diabetes + Metformin) 


dds1 = DESeq(dds1, test="Wald", fitType = "local", betaPrior =  T)
dds2 = DESeq(dds2, test="Wald", fitType = "local", betaPrior =  T)
dds3 = DESeq(dds3, test="Wald", fitType = "local", betaPrior =  T)

resultsNames(dds1)
res1 = results(dds1, contrast = c("Month","M12","M0"))
res2 = results(dds2, contrast = c("Month","M12","M0"))
res3 = results(dds3, contrast = c("Month","M12","M0"))


rbind(cbind(res1, Group = rep(1)), cbind(res2, Group = rep(2))) %>%
  rbind(cbind(res3, Group = rep(3))) -> M12_M0_G123
M12_M0_G123 %>% as.data.frame %>% arrange(desc(padj))
M12_M0_G123[which(M12_M0_G123$pvalue < 0.05), ] %>% as.data.frame %>% group_by(Group) %>%
  summary()

M12_M0_G123[which(M12_M0_G123$pvalue < 0.05), ] %>% as.data.frame 

#######################################
###########
########### only different in M12 but not in M0
##########
#############################################

ent.trim.M0 <- subset_samples(ent.trim, Month == "M0")
ent.trim.M6 <- subset_samples(ent.trim, Month == "M6")
ent.trim.M12 <- subset_samples(ent.trim, Month == "M12")

ddsM0 <- phyloseq_to_deseq2(ent.trim.M0, design = ~ Group) 
ddsM12 <- phyloseq_to_deseq2(ent.trim.M12, design = ~ Group) 

ddsM0 = DESeq(ddsM0, test="Wald", fitType = "local", betaPrior =  T)
ddsM12 = DESeq(ddsM12, test="Wald", fitType = "local", betaPrior =  T)

resG13_M12 <- results(ddsM12)
resultsNames(ddsM12)
resG12_M12 <- results(ddsM12, contrast = c("Group", "2","1"))
resG32_M12 <- results(ddsM12, contrast = c("Group", "3","2"))

resG12_M12

############ only stastically significant
resG13_M12[which(resG13_M12$pvalue < 0.05 ), ] -> 
  sigG13M12
resG12_M12[which(resG12_M12$pvalue < 0.05 ), ] -> 
  sigG12M12
resG32_M12[which(resG32_M12$pvalue < 0.05 ), ] -> 
  sigG32M12
###############

resG13_M12[which(resG13_M12$padj < 0.05 & abs(resG13_M12$log2FoldChange) > 1), ] -> 
  sigG13M12
resG12_M12[which(resG12_M12$padj < 0.05 & abs(resG12_M12$log2FoldChange) > 1), ] -> 
  sigG12M12
resG32_M12[which(resG32_M12$padj < 0.05 & abs(resG32_M12$log2FoldChange) > 1), ] -> 
  sigG32M12

sigG13M12

c(rownames(sigG12M12), rownames(sigG13M12), rownames(sigG32M12)) %>% unique() -> M12.allchanged

intersect(intersect(rownames(sigG12M12), rownames(sigG13M12)), rownames(sigG32M12)) -> M12changed
M12changed

resG13_M0 <- results(ddsM0)
resultsNames(ddsM0)
resG12_M0 <- results(ddsM0, contrast = c("Group", "2","1"))
resG32_M0 <- results(ddsM0, contrast = c("Group", "3","2"))

########### only consider significantly stastically
resG13_M0[which(resG13_M0$pvalue < 0.05 ), ] -> 
  sigG13M0
resG12_M0[which(resG12_M0$pvalue < 0.05 ), ] -> 
  sigG12M0
resG32_M0[which(resG32_M0$pvalue < 0.05 ), ] -> 
  sigG32M0
##############

resG13_M0[which(resG13_M0$padj < 0.05 & abs(resG13_M0$log2FoldChange) > 1), ] -> 
  sigG13M0
resG12_M0[which(resG12_M0$padj < 0.05 & abs(resG12_M0$log2FoldChange) > 1), ] -> 
  sigG12M0
resG32_M0[which(resG32_M0$padj < 0.05 & abs(resG32_M0$log2FoldChange) > 1), ] -> 
  sigG32M0

c(rownames(sigG13M0), rownames(sigG12M0), rownames(sigG32M0)) %>% unique() -> M0changed


M0changed
intersect(M12.allchanged, M0changed) -> tmp0

M12.allchanged[!M12.allchanged %in% tmp0] -> changedList
#M12.allchanged[!M12.allchanged %in% tmp0] -> changedList

changedList %>% length()
########### Normalized count data #############
counts(ddsM0, normalized = TRUE) -> ddsM0.n.c
counts(ddsM12, normalized = TRUE) -> ddsM12.n.c

ddsM12.n.c[rownames(ddsM0.n.c) %in% changedList,] -> ddsM12.n.c.1

ddsM0.n.c[rownames(ddsM0.n.c) %in% rownames(ddsM12.n.c.1),] -> ddsM0.n.c.1


rowSums(ddsM12.n.c.1==0) %>% enframe() %>% arrange(desc(value))
rowSums(ddsM0.n.c.1==0) %>% enframe() %>% arrange(desc(value))

ddsM12.n.c.1 - ddsM0.n.c.1 -> sp_raw_deltaM120_normalized.count


save(sp_raw_deltaM120_normalized.count, file="sp_raw_deltaM120_normalized.count.rda")


############### Relative_abundance data #########

otu_table(ex1b) -> ex1b.table

ex1b.table[rownames(ex1b.table) %in% changedList] %>% as.data.frame() %>% dplyr::select(contains("M0")) -> ex1b.table.M0
ex1b.table[rownames(ex1b.table) %in% changedList] %>% as.data.frame() %>% dplyr::select(contains("M12")) -> ex1b.table.M12

ex1b.table.M0 - ex1b.table.M12 -> sp_re_ab_deltaM120

save(sp_re_ab_deltaM120, file = "sp_re_ab_deltaM120.rda")

#################################################
############### Relative_abundance data #########

load(file="M12.Kruskal.changedList.rda")
load(file="M0.Kruskal.changedList.rda")

M12.Kruskal.changedList[! M12.Kruskal.changedList$key%in% M0.Kruskal.changedList$key,] -> changedList.Kruskal

otu_table(ex1b) -> ex1b.table

ex1b.table[rownames(ex1b.table) %in% changedList.Kruskal$key] %>% as.data.frame() %>% dplyr::select(contains("M0")) -> ex1b.Kruskal.table.M0
ex1b.table[rownames(ex1b.table) %in% changedList.Kruskal$key] %>% as.data.frame() %>% dplyr::select(contains("M12")) -> ex1b.Kruskal.table.M12

ex1b.Kruskal.table.M12 - ex1b.Kruskal.table.M0 -> sp_Kruskal_re_ab_deltaM120

save(sp_Kruskal_re_ab_deltaM120, file = "sp_Kruskal_re_ab_deltaM120.rda")

#################################################
############### Relative_abundance data #########

load(file="M12.Kruskal.changedList.rda")
load(file="M0.Kruskal.changedList.rda")

rbind(M12.Kruskal.changedList, M0.Kruskal.changedList) %>% distinct()  -> changedList.Kruskal

otu_table(ex1b) -> ex1b.table

ex1b.table[rownames(ex1b.table) %in% changedList.Kruskal$key] %>% as.data.frame() %>% dplyr::select(contains("M0")) -> ex1b.Kruskal.table.M0
ex1b.table[rownames(ex1b.table) %in% changedList.Kruskal$key] %>% as.data.frame() %>% dplyr::select(contains("M12")) -> ex1b.Kruskal.table.M12

ex1b.Kruskal.table.M12 - ex1b.Kruskal.table.M0 -> sp_Kruskal_re_ab_changedBothM0M12_deltaM120

save(sp_Kruskal_re_ab_changedBothM0M12_deltaM120, file = "sp_Kruskal_re_ab_changedBothM0M12_deltaM120.rda")

########### ############################
