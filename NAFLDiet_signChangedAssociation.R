source("NAFLDiet_phyloseq_v2_127_sample.R")

## NAFLDiet_phyloseq_v2_127_signChangedAssociation.R
library("psych")
library(ppcor)
library(gplots)
library("ComplexHeatmap")
library(circlize)

LiverRelatedVar.1 <- c("Study.ID", "ID",
                       "A2.ALT.(ukat/L)","B.ALT.(ukat/L)","C2.ALT.(ukat/L)",
                       "A2.ApoB/ApoA1","B.ApoB/ApoA1","C2.ApoB/ApoA1",
                       "A2.AST.(ukat/L)","B.AST.(ukat/L)","C2.AST.(ukat/L)",
                       "A2.GGT.(ukat/L)","B.GGT.(ukat/L)","C2.GGT.(ukat/L)",
                       "A2.BMI.(kg/m2)","B.BMI.(kg/m2)","C2.BMI.(kg/m2)",
                       "A2.HDL.Cholesterol.(mmol/L)","B.HDL.Cholesterol.(mmol/L)","C2.HDL.Cholesterol.(mmol/L)",
                       "A2.LDL.cholesterol.(mmol/L)", "B.LDL.cholesterol.(mmol/L)","C2.LDL.cholesterol.(mmol/L)",
                       "A2.Triglycerides.(mmol/L)","B.Triglycerides.(mmol/L)","C2.Triglycerides.(mmol/L)")

LiverRelatedVar.2 <- c("M0.ALT.(ukat/L)","M6.ALT.(ukat/L)","M12.ALT.(ukat/L)",
                       "M0.ApoB/ApoA1","M6.ApoB/ApoA1","M12.ApoB/ApoA1",
                       "M0.AST.(ukat/L)","M6.AST.(ukat/L)","M12.AST.(ukat/L)",
                       "M0.GGT.(ukat/L)","M6.GGT.(ukat/L)","M12.GGT.(ukat/L)",
                       "M0.BMI.(kg/m2)","M6.BMI.(kg/m2)","M12.BMI.(kg/m2)",
                       "M0.HDL.Cholesterol.(mmol/L)","M6.HDL.Cholesterol.(mmol/L)","M12.HDL.Cholesterol.(mmol/L)",
                       "M0.LDL.cholesterol.(mmol/L)", "M6.LDL.cholesterol.(mmol/L)","M12.LDL.cholesterol.(mmol/L)",
                       "M0.Triglycerides.(mmol/L)","M6.Triglycerides.(mmol/L)","M12.Triglycerides.(mmol/L)"
)



corr_simple <- function(corr,sig=0.5){
  #convert data to numeric in order to run correlations
  #convert to factor first to keep the integrity of the data - each value will become a number rather than turn into NA
  #df_cor <- data %>% mutate_if(is.character, as.factor)
  #df_cor <- df_cor %>% mutate_if(is.factor, as.numeric)  #run a correlation and drop the insignificant ones
  #corr <- cor(df_cor)
  #prepare to drop duplicates and correlations of 1     
  corr[lower.tri(corr,diag=TRUE)] <- NA 
  #drop perfect correlations
  corr[corr == 1] <- NA   #turn into a 3-column table
  corr <- as.data.frame(as.table(corr))
  #remove the NA values from above 
  corr <- na.omit(corr)   #select significant values  
  corr <- subset(corr, abs(Freq) > sig) 
  #sort by highest correlation
  corr <- corr[order(-abs(corr$Freq)),]   #print table
  print(corr)  #turn corr back into matrix in order to plot with corrplot
  mtx_corr <- reshape2::acast(corr, Var1~Var2, value.var="Freq")
  
  #plot correlations visually
  corrplot(mtx_corr, is.corr=FALSE, tl.col="black", na.label=" ")
}

corr_simple()



########
######### only significantly changed vs. liver enzymes
#########
load(file = "sp_re_ab_deltaM120.rda")
load(file="sp_raw_deltaM120_normalized.count.rda")
load(file="sample.127.388.liver.delta.rda")
load(file = "LiverFatFra.4.rda")


sp_raw_deltaM120_normalized.count 

sample.127.full %>% dplyr::select(contains("lipids")) %>% as.data.frame()

sample.127.388.liver.delta %>% dplyr::select(contains(c("3","ID"))) -> sample.127.388.liver.delta.3


## use abundance deseq2 ##
sp_raw_deltaM120_normalized.count %>% t() %>% as.data.frame() %>% rownames_to_column() %>%
  separate_wider_delim(cols = rowname, delim = "-",names = c("cID", "Month")) %>%
  separate_wider_delim(cols = cID, delim = "_",names = c("Prefix", "ID")) %>%
  left_join(sample.127.388.liver.delta.3) %>%
  mutate(ID=as.character(ID)) %>%
  left_join(LiverFatFra.4)  %>%
  dplyr::select(-c(Prefix, ID, Month, Triglycerides1, Triglycerides2, 
                   fatfraction_v1, fatfraction_v2, ApoB.ApoA13)) -> df1



colnames(df1)
df1 %>% dim()

#A straightforward application of matrix algebra to 
#remove the effect of the variables in the y set from the x set. DESeq2 ## abundance difference ###########
df1 [1,] 
df1
df1[1,1:75]
partial.r(df1, 
          c(1:75,77:83),
          c(76),
          use="pairwise",
          method="spearman") -> partial.r.result


aggr(df1, plot = F, numbers = T, prop = T)

pcor.test.result <- data.frame(matrix(ncol=3, nrow=0))
colnames(pcor.test.result) <- c("var1","var2","value")

### 83 deltafat, without deltaFat 
for (i in 1:75){
  for (j in c(77:82)){
    var1 <- colnames(df1)[i]
    var2 <- colnames(df1)[j]
    value <-  pcor.test(df1[,i], df1[,j],  df1[,c(76)], method="spearman")
    rbind(data.frame(var1=var1, var2 = var2, value=value), pcor.test.result) -> pcor.test.result
  }
}

#########################

aggr(df1, plot = F, numbers = T, prop = T)

df1 %>% na.omit() -> df2
for (i in 1:75){
  for (j in c(77:83)){
    var1 <- colnames(df2)[i]
    var2 <- colnames(df2)[j]
    value <-  pcor.test(df2[,i], df2[,j],  df2[,c(76)], method="spearman")
    rbind(data.frame(var1=var1, var2 = var2, value=value), pcor.test.result) -> pcor.test.result
  }
}

############## 

pcor.test.result %>% arrange(desc(abs(value.estimate))) -> pcor.test.result.1

pcor.test.result.1 %>% 
  dplyr::select(value.p.value) %>% as.matrix() %>% p.adjust(method = "fdr") %>% as.data.frame() -> fdr.result

fdr.result %>% tail()


fdr.result %>% 
  mutate_if(is.numeric, ~ 1 * (. < 0.1))  -> fdr.result.1

fdr.result.1[fdr.result.1 == 0] <- NA

colnames(fdr.result.1) <- "fdr.res"

cbind(pcor.test.result.1, fdr.result.1) -> pcor.test.result.all

pcor.test.result.all

tax_table(ex1b) %>% as.data.frame() %>% mutate(nTax = paste0(MGS,": ",MainTax)) -> tax.table.dict


tax.table.dict[rownames(tax.table.dict) %in% pcor.test.result.all$var1,] %>% rownames_to_column() %>%
  left_join(pcor.test.result.all, by=c( "rowname" = "var1"), keep=T) -> pcor.test.result.all

pcor.test.result.all %>% dplyr::select(var1, nTax)

tax.table.dict[rownames(tax.table.dict) %in% pcor.test.result.all$var1,] %>% dplyr::select(-nTax) %>%
  arrange(phylum) -> taxa.75


phy_tree(ex1b)

M12.Kruskal.changedList

library(MicrobiotaProcess)
library(ggtree)
help(as.treedata)
tree <- as.treedata(taxa.75, include.rownames=TRUE)

taxa.75 %>% names()



######## Top 20 associated species ######################

pcor.test.result.all %>% group_by(var1) %>%   
  summarise(sumV = sum(value.estimate)) %>% arrange(desc(abs(sumV))) %>% 
  head(n=10) %>%
  dplyr::select(var1) %>% distinct() -> tmp.var1

pcor.test.result.all %>%
  arrange(desc(abs(value.estimate))) %>% head(n=20) %>%
  dplyr::select(var1) -> tmp.var1.2

rbind(tmp.var1, tmp.var1.2) %>% distinct() -> tmp.var1.3

pcor.test.result.all %>% 
  filter(var1 %in% tmp.var1.3$var1) %>%
  arrange(desc(abs(value.estimate))) -> df3

df3$var2 <-gsub("3","",as.character(df3$var2))
df3$var2 <-gsub(".cholesterol","",as.character(df3$var2))

df3 %>% dplyr::select(var1) 
df3 %>% dplyr::select(nTax) %>% distinct()

pdf("try.pdf", width = 12)
df3 %>% ggplot(aes(y=nTax, x=var2, fill=value.estimate)) +
  scale_fill_gradient2(midpoint=0, 
                       low="darkred", 
                       mid="white",
                       high="blue", 
                       space ="Lab", 
                       breaks = c(-0.4,-0.2,0,0.2,0.4))+
  geom_tile() + 
  geom_point(aes(size=fdr.res),color="white") +
  clean_background
dev.off()

fdr.result.1  %>% dim()

##############################

cbind(pcor.test.result.1, fdr.result) -> pcor.test.result.all

pcor.test.result.all

tax_table(ex1b) %>% as.data.frame() %>% mutate(nTax = paste0(MGS,": ",MainTax)) -> tax.table.dict
rownames(tax.table.dict)

tax.table.dict[rownames(tax.table.dict) %in% pcor.test.result.all$var1,] %>% rownames_to_column() %>%
  left_join(pcor.test.result.all, by=c( "rowname" = "var1"), keep=T) -> pcor.test.result.all

df5 <- pcor.test.result.all

aggr(df5, plot=F)

df5 %>% dim()


library(circlize)


### tree circle ### 
tax.table.dict[rownames(tax.table.dict) %in% tmp.df$sp,] %>% rownames_to_column(var = "sp1") %>%
  left_join(tmp.df, by=c("sp1" = "sp"), keep = T) -> tmp0.df

tmp0.df %>% select(1:11)  %>% as.matrix() -> tmp0.1.df
tax75 <- phyloseq::tax_table(tmp0.1.df)

tree75 = rtree(ntaxa(tax75), rooted=TRUE, tip.label=taxa_names(tax75))
as.treedata(tree75) -> tree
ggtree(tree, layout = "circular")

dend <- as.dendrogram(tree75)
dend

labels = tree75$tip.label
n = length(labels)
ct = cutree(dend, k=6)
ct
n
pdf("try.pdf")
circos.par(cell.padding = c(0, 0, 0, 0))
circos.initialize("a", xlim = c(0, 75)) # only one sector
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.3, 
             panel.fun = function(x, y) {
               for(i in seq_len(n)) {
                 circos.text(i-0.5, 0, labels[i], adj = c(0, 0.5), 
                             facing = "clockwise", niceFacing = TRUE,
                             col = ct[labels[i]], cex = 0.5)
               }
             })

suppressPackageStartupMessages(library(dendextend))
dend = color_branches(dend, k = 6, col = 1:6)
dend_height = attr(dend, "height")
circos.track(ylim = c(0, dend_height), bg.border = NA, 
             track.height = 0.4, panel.fun = function(x, y) {
               circos.dendrogram(dend)
             })
circos.clear()

dev.off()


as.data.frame.matrix(partial.r.result)[1:75,76:82] %>% rownames_to_column(var = "sp")-> tmp.df
tmp.df
tax.table.dict[rownames(tax.table.dict) %in% tmp.df$sp,] %>% rownames_to_column(var = "sp1") %>%
  left_join(tmp.df, by=c("sp1" = "sp"), keep = T) %>% dplyr::select(nTax, phylum, ALT3, AST3, GGT3,  
                                                                    LDL.cholesterol3,
                                                                    HDL.cholesterol3, 
                                                                    Triglycerides3, deltaFraction) %>%
  setNames(c("nTax", "phylum", "ALT", "AST", "GGT", "LDL","HDL", "Trig.", "FatFract.")) %>%
  column_to_rownames(var="nTax")  %>% mutate(phylum = as.factor(phylum)) -> mat

mat
mat %>% str()
tax.table.dict

split = mat$phylum

Heatmap(mat[,2:8], row_split = split)
split 

mat %>% dim()

dev.off()
pdf("try.pdf")
circos.clear()

circos.par(gap.degree=c(1,1,1,1,1,1,1,10), "start.degree" = 30)
circos.heatmap.initialize(mat, split = split)
col_fun1 = colorRamp2(c(-0.3, 0, 0.3), c("red", "white", "blue"))
circos.heatmap(mat[,2:8],  col = col_fun1)#,rownames.side = "inside")
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  get.cell.meta.data("ycenter")
  if(CELL_META$sector.numeric.index == 8) { # the last sector
    cn = colnames(mat[,2:8])
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[1]+1.2, n), 
                1:n, cn, 
                cex = 0.5, adj = c(0, 0), facing = "bending.inside")
  }
}, bg.border = NA)

dev.off()
circos.clear()

               
circos.track(track.index = get.current.track.index(),
             panel.fun = function(x, y) {
               if(CELL_META$sector.numeric.index == 1) { # the last sector
                 cn = colnames(mat)
                 n = length(cn)
                 circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                             1:n - 0.5, cn, 
                             cex = 0.5, adj = c(0, 0.5), facing = "inside")
               }
             }, bg.border = NA)



circos.clear()

dev.off()


set.seed(123)
mat1 = rbind(cbind(matrix(rnorm(50*5, mean = 1), nr = 50), 
                   matrix(rnorm(50*5, mean = -1), nr = 50)),
             cbind(matrix(rnorm(50*5, mean = -1), nr = 50), 
                   matrix(rnorm(50*5, mean = 1), nr = 50))
)
rownames(mat1) = paste0("R", 1:100)
colnames(mat1) = paste0("C", 1:10)
mat1 = mat1[sample(100, 100), ] # randomly permute rows
split = sample(letters[1:5], 100, replace = TRUE)
split = factor(split, levels = letters[1:5])

col_fun1 = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
circos.heatmap(mat1, split = split, col = col_fun1)

col_fun1






###############
##########
######### Global association 20 per
########
######################
load(file = "re_ab_127_20per.clinc.rda")


re_ab_127_20per.clinc %>% dplyr::select_if(is.numeric) %>%  
  cor(use="pairwise.complete.obs", method="spearman") -> per20_cor

per20_cor %>% dim()

sum(is.infinite(per20_cor))
per20_cor[ is.na(per20_cor)] <- 0 

aggr(per20_cor, plot=F)

help(cor)

heatmap.2(per20_cor, scale = "none", na.rm = TRUE, col = bluered(500), 
          trace = "none", density.info = "none")



###############
##########
######### Global association 2363
########
######################

load(file = "re_ab_127_2363.clinc.rda")


re_ab_127_2363.clinc[,colnames(re_ab_127_2363.clinc) %in% LiverRelatedVar.1]

re_ab_127_2363.clinc %>% dplyr::select_if(is.numeric) %>%  
  cor(use="pairwise.complete.obs", method="spearman") -> all2363_cor

all2363_cor %>% dim()

sum(is.infinite(all2363_cor))
all2363_cor[ is.na(all2363_cor)] <- 0 

aggr(all2363_cor, plot=F)


all2363_cor[lower.tri(all2363_cor,diag=TRUE)] <- NA 
#drop perfect correlations
all2363_cor[all2363_cor == 1] <- NA   #turn into a 3-column table
all2363_cor <- as.data.frame(as.table(all2363_cor))
#remove the NA values from above 
all2363_cor <- na.omit(all2363_cor)   #select significant values  
all2363_cor.1 <- subset(all2363_cor, abs(Freq) > 0.1) 
#sort by highest correlation
all2363_cor.1 <- all2363_cor.1[order(-abs(all2363_cor.1$Freq)),]  

all2363_cor.1[1:5,1:3] 

all2363_cor.1 %>% filter(grepl("Bacteroides ster", Var1)) %>% filter(!grepl("MGS", Var2)) %>% filter(Var2 %in% LiverRelatedVar.1)

re_ab_127_2363.clinc %>% dplyr::select(contains("Bacteroides stercoris")) -> Bacteroides.Stercoris

colSums(Bacteroides.Stercoris==0)

228/381

all2363_cor.1  %>% filter(! grepl("MGS", Var2)) %>%  filter( grepl("MGS", Var1)) 

pdf("global_association2363.pdf")
heatmap.2(all2363_cor, scale = "none", na.rm = TRUE, col = bluered(500), 
          trace = "none", density.info = "none")

dev.off()

#################################################
##########
########## change of MGS M12-M0 vs. Liver parameters
##########
##################################################

re_ab_127  %>% dplyr::select(ends_with("M0")) -> re_ab_127.M0
re_ab_127  %>% dplyr::select(ends_with("M6")) -> re_ab_127.M6
re_ab_127  %>% dplyr::select(ends_with("M12")) -> re_ab_127.M12


re_ab_127_deltaM120
re_ab_127_deltaM120 <- re_ab_127.M12 - re_ab_127.M0 %>% round(3)


LiverFatFra %>% separate_wider_delim(patient_id, delim = "_",
                                     names=c("prefix","ID")) %>%
  mutate(ID = as.numeric(ID)) -> LiverFatFra.1



ID127.new %>% filter(! Studie.ID %in% LiverFatFra.1$ID) %>% 
  dplyr::select(nID)

sam127 %>% dplyr::select(ID, Group) %>% distinct() -> sam127.tmp
LiverFatFra.1 %>% filter(ID %in% sam127$ID) -> LiverFatFra.2

LiverFatFra.2 %>% dplyr::select(-c(date_v1, date_v2, prefix)) %>%
  mutate(deltaFatFrac = fatfraction_v2-fatfraction_v1) %>% round(3) -> LiverFatFra.3

re_ab_127_deltaM120 %>% t() %>% as.data.frame() %>% rownames_to_column() %>%
  separate_wider_delim(cols = rowname, delim = "-",names = c("cID", "Month")) %>%
  separate_wider_delim(cols = cID, delim = "_",names = c("Prefix", "ID")) %>% mutate(ID=as.numeric(ID)) %>%
  left_join(LiverFatFra.3) %>% dplyr::select(-c(Prefix,  Month)) -> re_ab_127_deltaM120.fatFrac

re_ab_127_deltaM120.fatFrac 
re_ab_127_2363.clinc

re_ab_127_2363.clinc[,colnames(re_ab_127_2363.clinc) %in% LiverRelatedVar.1] -> re_ab_127_2363.clinic.liver

re_ab_127_2363.clinic.liver %>% as.data.frame %>% mutate(ALT1=`B.ALT.(ukat/L)`- `A2.ALT.(ukat/L)`,
                                                   ALT2=`C2.ALT.(ukat/L)`- `B.ALT.(ukat/L)`,
                                                   ALT3=`C2.ALT.(ukat/L)`- `A2.ALT.(ukat/L)`,
                                                   AST1=`B.AST.(ukat/L)`- `A2.AST.(ukat/L)`,
                                                   AST2=`C2.AST.(ukat/L)`- `B.AST.(ukat/L)`,
                                                   AST3=`C2.AST.(ukat/L)`- `A2.AST.(ukat/L)`,
                                                   GGT1=`B.GGT.(ukat/L)`- `A2.GGT.(ukat/L)`,
                                                   GGT2=`C2.GGT.(ukat/L)`- `B.GGT.(ukat/L)`,
                                                   GGT3=`C2.GGT.(ukat/L)`- `A2.GGT.(ukat/L)`,
                                                   ApoB.ApoA11=`B.ApoB/ApoA1`- `A2.ApoB/ApoA1`,
                                                   ApoB.ApoA12=`C2.ApoB/ApoA1`- `B.ApoB/ApoA1`,
                                                   ApoB.ApoA13=`C2.ApoB/ApoA1`- `A2.ApoB/ApoA1`,
                                                   LDL.cholesterol1=`B.LDL.cholesterol.(mmol/L)`- `A2.LDL.cholesterol.(mmol/L)`,
                                                   LDL.cholesterol2=`C2.LDL.cholesterol.(mmol/L)`- `B.LDL.cholesterol.(mmol/L)`,
                                                   LDL.cholesterol3=`C2.LDL.cholesterol.(mmol/L)`- `A2.LDL.cholesterol.(mmol/L)`,
                                                   HDL.cholesterol1=`B.HDL.Cholesterol.(mmol/L)`- `A2.HDL.Cholesterol.(mmol/L)`,
                                                   HDL.cholesterol2=`C2.HDL.Cholesterol.(mmol/L)`- `B.HDL.Cholesterol.(mmol/L)`,
                                                   HDL.cholesterol3=`C2.HDL.Cholesterol.(mmol/L)`- `A2.HDL.Cholesterol.(mmol/L)`,
                                                   Triglycerides1=`B.Triglycerides.(mmol/L)`- `A2.Triglycerides.(mmol/L)`,
                                                   Triglycerides2=`C2.Triglycerides.(mmol/L)`- `B.Triglycerides.(mmol/L)`,
                                                   Triglycerides3=`C2.Triglycerides.(mmol/L)`- `A2.Triglycerides.(mmol/L)`
                                                   
                                                   
) %>% dplyr::select(-c(2:25)) %>% distinct() -> re_ab_127_2363.clinic.liver.delta

re_ab_127_2363.clinic.liver.delta  %>% head()
re_ab_127_deltaM120.fatFrac %>% dim()
re_ab_127_2363.clinic.liver.delta %>% mutate(ID=as.numeric(ID)) %>% 
  left_join(as.data.frame(re_ab_127_deltaM120.fatFrac), copy = TRUE)


#######################################################
#####
#### Correlation of the sample
####
########################################################

sample.127.388v %>% 
  mutate(diet.group = as.factor(diet.group), 
         A2.Gender = as.factor(A2.Gender),
         A2.Smoking = as.factor(A2.Smoking),
         A2.Education = as.factor(A2.Education)
  ) %>% dplyr::select_if(is.numeric) %>%
  cor(use="pairwise.complete.obs", method="spearman") -> sample_cor


##################
sum(is.infinite(sample_cor))
sample_cor[ is.na(sample_cor)] <- 0 

aggr(sample_cor, plot=F)


pdf("sample_spearman_association.pdf")
heatmap.2(sample_cor, scale = "none", na.rm = TRUE, col = bluered(500), 
          trace = "none", density.info = "none")
dev.off()
#####################

sample.127.388v %>% 
  mutate(diet.group = as.factor(diet.group), 
         A2.Gender = as.factor(A2.Gender),
         A2.Smoking = as.factor(A2.Smoking),
         A2.Education = as.factor(A2.Education)
  ) %>% dplyr::select_if(is.numeric) %>%  dplyr::select(starts_with("B.")) %>%
  cor(use="pairwise.complete.obs", method="spearman") -> sample_cor

sample_cor[lower.tri(sample_cor, diag=TRUE)] <- NA
sample_cor[sample_cor == 1] <-NA

aggr(sample_cor, plot=F)

sample_cor %>% as.table() %>% as.data.frame() %>% na.omit() %>%
  subset(abs(Freq) > 0.4) %>% arrange(desc(abs(Freq))) -> sample.cor.2

sample.cor.2

#sample.cor.2 %>% reshape2::acast(sample.cor.2, Var1 ~ Var2, value.var = "Freq")

#######################################################
#####
#### Correlation of the spices
####
########################################################

re_ab_127.t.2363 %>% 
  cor(use="pairwise.complete.obs", method="spearman") -> species_cor

sum(is.infinite(species_cor))
sample_cor[ is.na(species_cor)] <- 0 

aggr(species_cor, plot=F)

species_cor[lower.tri(species_cor, diag=TRUE)] <- NA
species_cor[species_cor == 1] <-NA

aggr(species_cor, plot=F)

species_cor %>% as.table() %>% as.data.frame() %>% na.omit() %>%
  subset(abs(Freq) > 0.4) %>% arrange(desc(abs(Freq))) -> species_cor_2

species_cor_2

pdf("Global_species_spearman_association.pdf")
heatmap.2(species_cor, scale = "none", na.rm = TRUE, col = bluered(500), 
          trace = "none", density.info = "none")
dev.off()
species_cor %>% dim()


