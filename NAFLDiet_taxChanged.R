source("NAFLDiet_phyloseq_v2_127_sample.R")
library(DESeq2)


#### Trim the Number of Zeros in more than 10%
0.1*381
0.2*381
0.5*381

sam1 <- sample_data(sam127)
tax1 <- tax_table(raw_tax_matrix_127)
ntaxa(tax1)

### transformation ###

relative_abundance_MGS_data_127 %>% head()

otu1 <- otu_table(relative_abundance_MGS_data_127, taxa_are_rows = TRUE)


random_tree = rtree(ntaxa(tax1), rooted=TRUE, tip.label=taxa_names(tax1))
ex1b <- phyloseq(otu1, sam1, tax1, random_tree)




ex1b %>% otu_table() %>% colSums()



flist<- filterfun(kOverA(38.1, 0))

ent.trim <- filter_taxa(ex1b, flist, TRUE)
ent.trim 

## filter unclassified ones 
ent.trim <- subset_taxa(ent.trim, ! phylum %in% c("unclassified"))
ent.trim

### which is then filtered such that
### only OTUs with a mean greater than 10^-5 are kept.

ent.trim <- filter_taxa(ent.trim, function(x) mean(x) > 1e-5, TRUE)
ent.trim

rank_names(ex1b)
table(tax_table(ent.trim)[,"phylum"],exclude=NULL)


dds1
rs <- rowSums(counts(dds1))
counts(dds)
rmx <- apply(counts(dds1), 1, max)
plot(rs+1, rmx/rs, log="x")


dds <- phyloseq_to_deseq2(ent.trim, design = ~ Month) 
dds = DESeq(dds, test="Wald", fitType = "local", betaPrior =  T)



re_ab_127_deltaM120.fatFrac %>% as.data.frame() %>%
  dplyr::select(c(contains("00925"), contains("fat")))  
cor() %>% reshape2::melt() %>% dplyr::filter(! value == "1") %>%
  arrange(desc(abs(value)))



ent.trim

ent.trim.G1 <- subset_samples(ent.trim, Group == 1)
ent.trim.G2 <- subset_samples(ent.trim, Group == 2)
ent.trim.G3 <- subset_samples(ent.trim, Group == 3)

dds1 <- phyloseq_to_deseq2(ent.trim.G1, design = ~ Month) 
dds2 <- phyloseq_to_deseq2(ent.trim.G2, design = ~ Month) 
dds3 <- phyloseq_to_deseq2(ent.trim.G3, design = ~ Month) 


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

###########
###########
## only different in M12 but not in M0
##########
##########

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
changedList
counts(ddsM0, normalized = TRUE) -> ddsM0.n.c
counts(ddsM12, normalized = TRUE) -> ddsM12.n.c

ddsM12.n.c[rownames(ddsM0.n.c) %in% changedList,] -> ddsM12.n.c.1

ddsM0.n.c[rownames(ddsM0.n.c) %in% rownames(ddsM12.n.c.1),] -> ddsM0.n.c.1


rowSums(ddsM12.n.c.1==0) %>% enframe() %>% arrange(desc(value))
rowSums(ddsM0.n.c.1==0) %>% enframe() %>% arrange(desc(value))



ddsM12.n.c.1 - ddsM0.n.c.1 -> sp_raw_deltaM120

save(sp_raw_deltaM120, file="sp_raw_deltaM120.rda")



rbind(cbind(res1, Group = rep(1)), cbind(res2, Group = rep(2))) %>%
  rbind(cbind(res3, Group = rep(3))) -> M12_M0_G123
M12_M0_G123 %>% as.data.frame %>% arrange(desc(padj))
M12_M0_G123[which(M12_M0_G123$pvalue < 0.05), ] %>% as.data.frame %>% group_by(Group) %>%
  summary()

for (i in c("res1", "res2", "res3")){
  res = get(i)
  sigtab1 = res[which(res$padj < 0.05 ), ]
  dim(sigtab1) 
}


resResults = function(res){
  sigtab1 = res[which(res$padj < 0.05 ), ]
  sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(ent.raw.trim)[rownames(sigtab1), ], "matrix"))
  sigtab1 %>% dim()  
}


sigtab1 = res1[which(res1$pvalue < 0.05 ), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(ent.raw.trim)[rownames(sigtab1), ], "matrix"))
sigtab1 %>% dim()
sigtab1



rld <- rlog(dds1)
plotPCA(rld, intgroup="Month")
dev.copy2pdf(file = "try.pdf")


#### do variable stabilization ###
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
ddsvst = getVarianceStabilizedData(dds)
ddsvst %>% colSums()
dds
ent.raw.trim0 = ent.raw.trim
otu_table(ent.raw.trim) <- otu_table(ddsvst, taxa_are_rows = TRUE)
dds = DESeq(dds, test="Wald", fitType = "parametric")
res1 = results(dds, cooksCutoff = FALSE)
res1
resultsNames(dds)


sigtab1 = res1[which(res1$log2FoldChange > 1.5 ), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(ent.raw.trim)[rownames(sigtab1), ], "matrix"))
sigtab1 %>% dim()
sigtab1

sigtab1[rownames(sigtab1) %in% filtedOut20per, ] %>% dim()

Group.dds <- phyloseq_to_deseq2(ent.trim, ~ Group)
Group.dds
Group.dds = DESeq(Group.dds, test="Wald", fitType = "parametric")
res2 = results(Group.dds, cooksCutoff = FALSE)

sigtab2 = res2[which(res$log2FoldChange > 1 ), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(ent.trim)[rownames(sigtab2), ], "matrix"))
sigtab2 %>% dim()

sigtab2

res
res.1 = cbind(as(res, "data.frame"), as(tax_table(ent.trim)[rownames(res), ], "matrix"))
library(circlize)

# create matrix
mat1 <- matrix(runif(80), 10, 8)
mat2 <- matrix(runif(80), 10, 8)
rownames(mat1) <- rownames(mat2) <- paste0('a', 1:10)
colnames(mat1) <- colnames(mat2) <- paste0('b', 1:8)

# join together
matX <- cbind(mat1, mat2)
matX 
# set splits
split <- c(rep('a', 8), rep('b', 8))
split = factor(split, levels = unique(split))

library(Cairo)
help("circos.heatmap")
# create circular heatmap
col_fun1 = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
circos.heatmap(matX, split = split, col = col_fun1,  rownames.side = "outside")
circos.clear()

circos.heatmap(t(matX), split = split, col = col_fun1,  rownames.side = "inside")
set.seed(999)
n = 1000
df = data.frame(sectors = sample(letters[1:8], n, replace = TRUE),
                x = rnorm(n), y = runif(n))
df  %>% head()

circos.par("track.height" = 0.1)
circos.initialize(df$sectors, x = df$x)
circos.track(df$sectors, y = df$y,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(5), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })

col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(df$sectors, df$x, df$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "text", sector.index = "a", track.index = 1)
bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(df$sectors, df$x, bin.size = 0.2, bg.col = bgcol, col = NA)
circos.track(df$sectors, x = df$x, y = df$y,
             panel.fun = function(x, y) {
               ind = sample(length(x), 10)
               x2 = x[ind]
               y2 = y[ind]
               od = order(x2)
               circos.lines(x2[od], y2[od])
             })
circos.update(sector.index = "d", track.index = 2, 
              bg.col = "#FF8080", bg.border = "black")
circos.points(x = -2:2, y = rep(0.5, 5), col = "white")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "updated", col = "white")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  breaks = seq(xlim[1], xlim[2], by = 0.1)
  n_breaks = length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
              breaks[-1], rep(ylim[2], n_breaks - 1),
              col = rand_color(n_breaks), border = NA)
})



