rm(list=ls()[])
library(vegan)
library(phyloseq)
library(tidyverse)
#library(mia)
library(ggplot2)
library(microbiome)
library(genefilter)
library(ggsignif)
library(superb)
library(car)
library(ggpubr)
#library(webr)
library(multcomp)
library(reshape2)
library(magrittr)
library(ade4)
library(factoextra)
#library(rayshader)
library(readxl)
library(mice)
library(mvnormtest)
library("ape")
library(VIM)
library(phylogram)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("tax_table", "phyloseq")
conflict_prefer("select", "dplyr")
conflicts_prefer(base::cbind)
conflicts_prefer(base::rbind)




relative_abundance_MGS_data <- read_xlsx("all-tables.xlsx",
                                         sheet = "MGS-downsized-relative-abundanc")



rownames(relative_abundance_MGS_data)
relative_abundance_MGS_data[,-c(1:11)] %>% 
  dplyr::select(c(1:415)) %>% 
  as.matrix() -> relative_abundance_matrix


relative_abundance_MGS_data[,c(2:11)] %>% as.matrix() -> raw_taxonomy_matrix


colnames(relative_abundance_MGS_data[,-c(1:11)]) %>% as.data.frame()-> tmp0
colnames(tmp0) <- c("ID")

data_grp  <- tmp0 %>% mutate(ID.c = ID)%>% filter( grepl("[-]", ID)) %>% 
  separate_wider_delim(cols = ID.c, delim = "-",names = c("nID", "Month")) %>%
  column_to_rownames(var="ID")

##### Select the 127 patients with M0, M6 and M12 data #####

data_grp %>% filter(Month == "M0") %>% dplyr::select(nID)  -> ID135.0
data_grp %>% filter(Month == "M6") %>% dplyr::select(nID)  -> ID135.6
data_grp %>% filter(Month == "M12") %>% dplyr::select(nID)  -> ID135.12
merge(merge(ID135.0, ID135.12),ID135.6) -> ID127
ID127 %>% mutate(Studie.ID = gsub("lunliv_","",nID)) -> ID127.new


relative_abundance_MGS_data[,-c(1:11)] %>% colnames() %>% enframe() %>%
  mutate(samples=value) %>%  filter( grepl("[-]", samples)) %>%
  separate_wider_delim(cols = samples, delim = "-",names = c("cID", "Month")) %>%
  separate_wider_delim(cols = cID, delim = "_",names = c("Prefix", "ID")) %>%
  dplyr::select(value, ID, Month) %>% column_to_rownames(var="value") -> sample.ID.Month




########## sample S2 #########


rowSums(is.na(sample_NAFLDiet)) %>% enframe() %>% cbind(sample_NAFLDiet$`Studie-ID`) -> sample.s0

sample_NAFLDiet %>% filter(`Studie-ID` %in% ID127.new$Studie.ID)  %>%
  mutate(`Studie-ID` = as.character(`Studie-ID`)) -> sample.127.full



colnames(sample_NAFLDiet) 

sample_NAFLDiet %>% dplyr::select(`Studie-ID`, Kostgrupp, `A2 Ålder (år)`,`A2 Kön`, 
                                  `A2 BMI (kg/m2)`,`A2 Diabetesdiagnos`,
                                  `A2 Metformin dos (mg)`) -> sample.s


sample_NAFLDiet %>% dplyr::select(`Studie-ID`, Kostgrupp, matches("A2")) -> sample.s
sample.s %>% str()

sample.s %>%  mutate(ID=as.character(`Studie-ID`), 
                     Group=as.factor(Kostgrupp),
                     age = as.numeric(`A2 Ålder (år)`),
                     BMI = as.numeric(`A2 BMI (kg/m2)`),
                     gender = as.factor(`A2 Kön`),
                     Diabetes =as.factor(`A2 Diabetesdiagnos`),
                    Metformin = as.factor(ifelse(is.na(`A2 Metformin dos (mg)`),0,1))) -> sample.s1
sample.s1 %>% dplyr::select(ID,Group,age, gender, BMI, Diabetes, Metformin) -> sample.s1
sample.s1

### replace the two NA age with mean value.


tempData <- mice(sample.s1[,1:7], m=5, maxit = 50, meth="pmm", seed=500)

sample.s2 <- complete(tempData,1)


sample.s2 %>% filter(ID %in% ID127.new$Studie.ID) -> sample.s127
sample.s %>% filter(`Studie-ID` %in% ID127.new$Studie.ID) %>% 
  mutate(ID=as.character(`Studie-ID`))-> sample.A2.127 


dplyr::left_join(sample.ID.Month, sample.s2) -> sam0
sam0
rownames(sam0) <- rownames(sample.ID.Month)


sample.ID.Month %>% rownames_to_column() %>% dplyr::left_join(sample.A2.127)%>%
  column_to_rownames(var="rowname") -> sample.A2.127.1

#### select the 127 individuals

sam0[sam0$ID %in% ID127.new$Studie.ID,] -> sam127

sample.A2.127.1[sample.A2.127.1$ID %in% ID127.new$Studie.ID,] -> sample.A2.127.2

relative_abundance_MGS_data_127 <- relative_abundance_MGS_data %>% dplyr::select(rownames(sam127)) 

cbind(relative_abundance_MGS_data[,c(2:11)], relative_abundance_MGS_data_127) %>% 
  dplyr::select(c(1:10)) %>% as.matrix()-> raw_tax_matrix_127

otu1 <- otu_table(relative_abundance_MGS_data_127, taxa_are_rows = TRUE)
sam1 <- sample_data(sam127)
tax1 <- phyloseq::tax_table(raw_tax_matrix_127)
ntaxa(tax1)

random_tree = rtree(ntaxa(tax1), rooted=TRUE, tip.label=taxa_names(tax1))

plot(random_tree)
dnd <- as.dendrogram(random_tree)


ex1b <- phyloseq(otu1, sam1, tax1, random_tree)

ex1b.M0 <- subset_samples(ex1b, Month == "M0")
ex1b.M6 <- subset_samples(ex1b, Month == "M6")
ex1b.M12 <- subset_samples(ex1b, Month == "M12")

relative_abundance_MGS_data[,rownames(sam127)] -> re_ab_127



ENcolnames <- read.table("variableENname.v1.csv", sep=",")
ENcolnames
colnames(sample.127.full) <- ENcolnames
colnames(sample_NAFLDiet) <- ENcolnames


###########
clean_background <- theme(plot.background = element_rect("white"),
                          panel.background = element_rect("white"),
                          panel.grid = element_line("white"),
                          axis.line = element_line("gray25"),
                          axis.text = element_text(size = 12, color = "gray25"),
                          axis.title = element_text(color = "gray25"),
                          legend.text = element_text(size = 12),
                          legend.key = element_rect("white"))
pal <- c("lightsalmon1", "gold1", "palegreen4")

