source("NAFLDiet_sample.R")
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(microbiomeutil) # some utility tools
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling

#########################
########### TSNE ########
########################
library(Rtsne) # Load package

set.seed(423542)

method <- "tsne"
trans <- "hellinger"
distance <- "euclidean"


# Distance matrix for samples
ps <- microbiome::transform(re_ab_127,trans)
ps <- microbiome::transform(re_ab_127,'log10p')


# Calculate sample similarities
dm <- vegdist(otu_table(ps, taxa_are_rows = TRUE), distance)

# Run TSNE
tsne_out <- Rtsne(dm, dims = 2) 
proj <- tsne_out$Y
rownames(proj) <- rownames(otu_table(ps))

pdf("try.pdf")
p <- plot_landscape(proj, legend = T, size = 1) 
print(p)
dev.off()

ex1b.G1 <- subset_samples(ex1b,  Month != "M12")
ex1b.G2 <- subset_samples(ex1b, Group == "2" & Month != "M12")
ex1b.G3 <- subset_samples(ex1b, Group == "3" & Month != "M12")
otu <- abundances(ex1b.G1)
metadata <- meta(ex1b.G1)
metadata

rda.result <- vegan::rda(t(otu) ~ factor(metadata$Group),
                         na.action = na.fail, scale = TRUE)


plot(rda.result, choices = c(1,2), type = "points", pch = 15, scaling = 3, cex = 0.7, col = metadata$Group)
points(rda.result, choices = c(1,2), pch = 15, scaling = 3, cex = 0.7, col = metadata$Group)
pl <- ordihull(rda.result, metadata$Group, scaling = 3, label = TRUE)

##############
############## unweighted unifrac
############
#######
# if we remove OTUs that are detected at least 10 times in 5% of the samples
ex1b.rar.filtered <- core(ex1b, detection = 10, prevalence = 0.05)

summarize_phyloseq(ex1b.rar.filtered)
ex1b
ex1b.G1 <- subset_samples(ex1b, Group == "1" & Month != "M12")
ex1b.G2 <- subset_samples(ex1b, Group == "2" & Month != "M12")
ex1b.G3 <- subset_samples(ex1b, Group == "3" & Month != "M12")

ex1b.M60 <- subset_samples(ex1b, Month != "M12") %>% transform('log10p')
ex1b.M612 <- subset_samples(ex1b, Month != "M0")
ex1b.M012 <- subset_samples(ex1b, Month != "M6")

ordu.unwt.uni.M60 <- ordinate(ex1b.M60, "PCoA", "unifrac", weighted=F)
ordu.unwt.uni.M612 <- ordinate(ex1b.M612, "PCoA", "unifrac", weighted=F)
ordu.unwt.uni.M012 <- ordinate(ex1b.M012, "PCoA", "unifrac", weighted=F)


pdf("log10pTransformedWeightedM60.pdf")
ex1b.M60 <- subset_samples(ex1b, Month != "M12") %>% transform('log10p')
ordu.wt.uni.M60 <- ordinate(ex1b.M60, "PCoA", "unifrac", weighted=T)
wt.unifracM60 <- plot_ordination(ex1b.M60,
                                   ordu.wt.uni.M60,
                                   color="Group")
wt.unifracM60 <- wt.unifracM60 + theme_classic() + scale_color_brewer("Group", palette = "Set2")
print(wt.unifracM60 + stat_ellipse())
dev.off()





pdf("notTransformedWeightedM60.pdf")
ex1b.M60 <- subset_samples(ex1b, Month != "M12")
ordu.wt.uni.M60 <- ordinate(ex1b.M60, "PCoA", "unifrac", weighted=T)
wt.unifracM60 <- plot_ordination(ex1b.M60,
                                 ordu.wt.uni.M60,
                                 color="Group")
wt.unifracM60 <- wt.unifracM60 + theme_classic() + scale_color_brewer("Group", palette = "Set2")
print(wt.unifracM60 + stat_ellipse())
dev.off()

help(ordinate)

ordu.bray.M60 <- ordinate(ex1b.M60, "PCoA", "bray")

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA") #,"NMDS", "MDS", "PCoA")
ord_meths = c( "NMDS", "MDS")

llply(as.list(ord_meths), function(i){
  print(i)
})
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq, method=i, distance=dist)
  plot_ordination(physeq, ordi,  color="Group", shape="Group")
}, ex1b.M60, dist)

plist

names(plist) <- ord_meths

plist

pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"


p = ggplot(pdataframe, aes(Axis_1, Axis_2, color=Month, shape=Group, fill=Month))
p = p + geom_point(size=4) + geom_polygon()
p = p + facet_wrap(~method, scales="free")
p = p + scale_fill_brewer(type="qual", palette="Set1")
p = p + scale_colour_brewer(type="qual", palette="Set1")
p

# check for Eigen values 
barplot(ordu.bray.M60$values$Eigenvalues[1:10])

ordu.bray.M60 <- plot_ordination(ex1b.M60,
                                 ordu.bray.M60,
                                 color="Group")
ordu.bray.M60 <- ordu.bray.M60 + theme_classic() + scale_color_brewer("Group", palette = "Set2")
print(ordu.bray.M60 + stat_ellipse())


ordu.bray.M60$data%>% as.data.frame() %>% 
  pivot_wider(names_from = Month,
              values_from = c(Axis.1,Axis.2)) %>%
  mutate(sd=(Axis.1_M6-Axis.1_M0)*(Axis.1_M6-Axis.1_M0)+
           (Axis.2_M6-Axis.2_M0)*(Axis.2_M6-Axis.2_M0)) %>% summary()

pdf("Bray-Curtis_distance_M60.pdf")
ordu.bray.M60$data%>% as.data.frame() %>% 
  pivot_wider(names_from = Month,
              values_from = c(Axis.1,Axis.2)) %>%
  mutate(sd=(Axis.1_M6-Axis.1_M0)**2+
           (Axis.2_M6-Axis.2_M0)**2) %>%
  ggplot(aes(x=Axis.1_M0, y=Axis.2_M0,  label=ID)) +
  geom_point(shape=21, size=3, fill="pink") +
  #geom_text()+
  geom_point(aes(x=Axis.1_M6, y=Axis.2_M6), 
             fill="chocolate1", shape=23, size=3) +
  geom_segment(aes(x=Axis.1_M0, y=Axis.2_M0, 
                   xend=Axis.1_M6, yend=Axis.2_M6, color=sd),
               arrow=arrow(length = unit(0.2, "cm")),
               color="black")+
  scale_color_gradient2(midpoint=5.667e-03, 
                        low="gray90", 
                        mid="gray75",
                        high="gray45", 
                        space ="Lab" )+
  labs(
    x = "PCoA.1 [7.7%]",
    y = "PCoA.2 [5.5%]",
    title = "Bray-Curtis Distance (M0 -> M6)"
  )+
  facet_wrap(~Group, strip.position = "left", ncol = 1) +
  clean_background +
  theme(panel.background = element_rect(colour = "gray", size=1)) +
  guides(color="none")
dev.off()


ordu.bray.M612$data%>% as.data.frame() %>% 
  pivot_wider(names_from = Month,
              values_from = c(Axis.1,Axis.2)) %>%
  mutate(sd=(Axis.1_M6-Axis.1_M12)*(Axis.1_M6-Axis.1_M12)+
           (Axis.2_M6-Axis.2_M12)*(Axis.2_M6-Axis.2_M12)) %>% summary()

ordu.bray.M612 <- ordinate(ex1b.M612, "PCoA", "bray")
ordu.bray.M612 <- plot_ordination(ex1b.M612,
                                  ordu.bray.M612,
                                  color="Group")
ordu.bray.M612 <- ordu.bray.M612 + theme_classic() + scale_color_brewer("Group", palette = "Set2")
print(ordu.bray.M612 + stat_ellipse())


pdf("Bray-Curtis_distance_M612.pdf")

ordu.bray.M612$data%>% as.data.frame() %>% 
  pivot_wider(names_from = Month,
              values_from = c(Axis.1,Axis.2)) %>%
  mutate(sd=(Axis.1_M6-Axis.1_M12)*(Axis.1_M6-Axis.1_M12)+
           (Axis.2_M6-Axis.2_M12)*(Axis.2_M6-Axis.2_M12)) %>% 
  ggplot(aes(x=Axis.1_M6, y=Axis.2_M6,  label=ID)) +
  geom_point(shape=21, size=3, fill="pink") +
  #geom_text()+
  geom_point(aes(x=Axis.1_M12, y=Axis.2_M12), 
             fill="chocolate1", shape=23, size=3) +
  scale_color_gradient2(midpoint=5.273e-03, 
                        low="gray90", 
                        mid="gray75",
                        high="gray45", 
                        space ="Lab" )+
  geom_segment(aes(x=Axis.1_M6, y=Axis.2_M6, 
                   xend=Axis.1_M12, yend=Axis.2_M12, color=sd),
               arrow=arrow(length = unit(0.2, "cm")),
               color="black")+
  labs(
    x = "PCoA.1 [6.2%]",
    y = "PCoA.2 [5.1%]",
    title = "Bray-Curtis Distance (M6 -> M12)"
  )+
  facet_wrap(~Group, strip.position = "left", ncol=1) +
  clean_background +
  theme(panel.background = element_rect(colour = "gray", size=1)) +
  guides(color="none")
dev.off()





unwt.unifracM60 <- plot_ordination(ex1b.M60,
                                   ordu.unwt.uni.M60,
                                   color="Group")
unwt.unifracM60 <- unwt.unifracM60 + theme_classic() + scale_color_brewer("Group", palette = "Set2")
print(unwt.unifracM60 + stat_ellipse())

unwt.unifracM60$data %>% as.data.frame() %>% 
  pivot_wider(names_from = Month,
              values_from = c(Axis.1,Axis.2)) %>%
  mutate(sd=(Axis.1_M6-Axis.1_M0)*(Axis.1_M6-Axis.1_M0)+
           (Axis.2_M6-Axis.2_M0)*(Axis.2_M6-Axis.2_M0)) %>% summary()

pdf("unweighted_UniFrac_distance_M60.pdf")
unwt.unifracM60$data %>% as.data.frame() %>% 
  pivot_wider(names_from = Month,
              values_from = c(Axis.1,Axis.2)) %>%
  mutate(sd=(Axis.1_M6-Axis.1_M0)*(Axis.1_M6-Axis.1_M0)+
           (Axis.2_M6-Axis.2_M0)*(Axis.2_M6-Axis.2_M0)) %>%
  ggplot(aes(x=Axis.1_M0, y=Axis.2_M0,  label=ID)) +
  geom_point(shape=21, size=5, fill="pink") +
  #geom_text()+
  geom_point(aes(x=Axis.1_M6, y=Axis.2_M6), 
             color="chocolate1", shape=17, size=3) +
  geom_segment(aes(x=Axis.1_M0, y=Axis.2_M0, 
                   xend=Axis.1_M6, yend=Axis.2_M6, color=sd),
               arrow=arrow(length = unit(0.2, "cm")))+
  scale_color_gradient2(midpoint=1.143e-03, 
                        low="gray90", 
                        mid="gray75",
                        high="gray45", 
                        space ="Lab" )+
  labs(
    x = "PCoA.1 [11.2%]",
    y = "PCoA.2 [6%]",
    title = "Bray (M0 -> M6)",
    subtitle = "unweighted UniFrac distance"
          )+
  facet_wrap(~Group) +
  clean_background +
  theme(panel.background = element_rect(colour = "gray", size=1)) +
  guides(color="none")
dev.off()

ordu.unwt.uni <- ordinate(ex1b, "PCoA","unifrac",weighted=F)
ordu.wt.uni <- ordinate(ex1b, "PCoA","unifrac",weighted=T)

ordu.unwt.uni.filter <- ordinate(ex1b, "PCoA", "unifrac", weighted=F)

ordu.unwt.uni.G1 <- ordinate(ex1b.G1, "PCoA","unifrac",weighted=F)
ordu.unwt.uni.G2 <- ordinate(ex1b.G2, "PCoA","unifrac",weighted=F)
ordu.unwt.uni.G3 <- ordinate(ex1b.G3, "PCoA","unifrac",weighted=F)

ordu.unwt.uni.G1


# check for Eigen values 
barplot(ordu.unwt.uni$values$Eigenvalues[1:10])

unwt.unifrac <- plot_ordination(ex1b.G1, 
                                ordu.unwt.uni.G1, color="Month") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("Month", palette = "Set2")
print(unwt.unifrac + stat_ellipse())

unwt.unifrac2 <- plot_ordination(ex1b.G2, 
                                ordu.unwt.uni.G2, color="Month") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("Month", palette = "Set2")
print(unwt.unifrac + stat_ellipse())

unwt.unifrac3 <- plot_ordination(ex1b.G3, 
                                ordu.unwt.uni.G3, color="Month") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("Month", palette = "Set2")
print(unwt.unifrac + stat_ellipse())
require(gridExtra)
unwt.unifrac$data %>% ggplot(aes(x=Axis.1, y=Axis.2, color=Month, shape=Month))+
  geom_point() +
  clean_background
unwt.unifrac2
unwt.unifrac2$data %>% as.data.frame() %>% 
  pivot_wider(names_from = Month,
              values_from = c(Axis.1,Axis.2)) %>%
  mutate(sd=(Axis.1_M6-Axis.1_M0)*(Axis.1_M6-Axis.1_M0)+
           (Axis.2_M6-Axis.2_M0)*(Axis.2_M6-Axis.2_M0)) %>%
  ggplot(aes(x=Axis.1_M0, y=Axis.2_M0,  label=ID)) +
  geom_point(shape=21, size=5, fill="pink") +
  #geom_text()+
  geom_point(aes(x=Axis.1_M6, y=Axis.2_M6), 
             color="chocolate1", shape=17, size=3) +
  geom_segment(aes(x=Axis.1_M0, y=Axis.2_M0, 
                   xend=Axis.1_M6, yend=Axis.2_M6),
               arrow=arrow(length = unit(0.2, "cm")),
               color="black")+
  clean_background


help("pivot_wider")
### filtered
unwt.unifrac
# check for Eigen values 
barplot(ordu.unwt.uni.filter$values$Eigenvalues[1:10])

unwt.unifrac.rar.filtered <- plot_ordination(ex1b.rar.filtered, 
                                             ordu.unwt.uni.filter, color="Group") 
unwt.unifrac.rar.filtered <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac.rar.filtered <- unwt.unifrac + theme_classic() + scale_color_brewer("Group", palette = "Set2")
print(unwt.unifrac.rar.filtered)

###############  metaMDS
################
##################

re_ab_127_NMDS <- metaMDS(t(re_ab_127.M12), k=3)

re_ab_127_NMDS$points %>% as.data.frame() -> tmp_NMDS

cbind(tmp_NMDS,sam127) -> NMDS.2
NMDS.2
NMDS.2 %>% plot_ly(x = ~MDS1, y = ~MDS2, z = ~MDS3) %>%
  add_markers(color=~Group) 


ggplot(NMDS.2, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point()


