source("NAFLDiet_sample.R")

library(mediation)

sample_NAFLDiet %>% filter(`Studie-ID` %in% ID127.new$Studie.ID) -> sample.127.full

sample.127.full %>% colnames() 

ENcolnames <- read.table("variableENname.v1.csv", sep=",")
ENcolnames
colnames(sample.127.full) <- ENcolnames
colnames(sample_NAFLDiet) <- ENcolnames


####### alpha ########
div <- microbiome::alpha(ex1b, index="all")
div_shannon <- diversity((relative_abundance_matrix), index = "shannon")
div_shannon %>% filter(row.names(div_shannon) %in% rownames(sam127)) -> div_shannon

rownames(div_shannon) == rownames(sam127)

### vegdist(t(tmp1), binary = TRUE) is the Sørensen index of dissimilarity

relative_abundance_MGS_data[,rownames(sam127)] -> re_ab_127

re_ab_127  %>% dplyr::select(ends_with("M0")) -> re_ab_127.M0
re_ab_127  %>% dplyr::select(ends_with("M6")) -> re_ab_127.M6
re_ab_127  %>% dplyr::select(ends_with("M12")) -> re_ab_127.M12


re_ab_127[rowSums(re_ab_127) != 0,] -> re_ab_127.2363

re_ab_127.2363 %>% mutate_if(is.numeric, ~ 1 * (. != 0)) %>% colSums() %>% enframe() %>%
  setNames(c("nID","NoOfSpecies")) -> speciesNo127




