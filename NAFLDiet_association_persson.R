source("NAFLDiet_phyloseq_v2_127_sample.R")


library(Hmisc)


# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

vars = c("riglycerides", "BMI", "liver.fat", 
         "AST", "ALT", "GGT", "ApoB/ApoA1",
         "holesterol")

LiverRelatedVar <- c("A2.ALT.(ukat/L)", "A2.AST.(ukat/L)", "A2.GGT.(ukat/L)",
                     "A2.BMI.(kg/m2)", "A2.ApoB/ApoA1", 
                     "A2.HDL.Cholesterol.(mmol/L)",  "A.LDL.cholesterol.(mmol/L)",
                     "A2.Triglycerides.(mmol/L)", 
                     "B.ALT.(ukat/L)",  "B.AST.(ukat/L)", "B.GGT.(ukat/L)",
                     "B.BMI.(kg/m2)", "B.ApoB/ApoA1", 
                     "B.HDL.Cholesterol.(mmol/L)", "B.LDL.cholesterol.(mmol/L)",
                     "B.Triglycerides.(mmol/L)", 
                     "C2.ALT.(ukat/L)", "C2.AST.(ukat/L)", "C2.GGT.(ukat/L)", 
                     "C2.BMI.(kg/m2)", "C2.ApoB/ApoA1", 
                     "C2.HDL.Cholesterol.(mmol/L)",  "C2.LDL.cholesterol.(mmol/L)",
                     "C2.Triglycerides.(mmol/L)")

LiverRelatedVar.1 <- c("A2.ALT.(ukat/L)","B.ALT.(ukat/L)","C2.ALT.(ukat/L)",
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




LiverRelatedVar.1 %>% length()

sample.127.full %>% colnames() 




sample.127.full[colSums(is.na(sample.127.full)) < 3] -> sample.127.388v
 

rownames(relative_abundance_MGS_data) <- paste0(relative_abundance_MGS_data$MGS, 
                                                sep=": ",
                                                relative_abundance_MGS_data$MainTax
)

relative_abundance_MGS_data %>% dplyr::select(rownames(sam127)) -> re_ab_127.0
relative_abundance_MGS_data %>% dplyr::select(MGS, MainTax, rownames(sam127)) -> re_ab_127.mgs

re_ab_127.mgs %>% mutate(rn = paste0(MGS, ": ", MainTax)) %>% column_to_rownames(var="rn") %>%
  dplyr::select(-c(MGS, MainTax)) -> re_ab_127


rownames(re_ab_127)

t(re_ab_127) -> re_ab_127.t

re_ab_127.t[,colSums(re_ab_127.t) !=0] -> re_ab_127.t.2363



sample.ID.Month[rownames(sample.ID.Month) %in% rownames(re_ab_127.t.2363),] %>% 
  cbind(re_ab_127.t.2363) -> re_ab_127.1


### filter 20% ###
re_ab_127.t[,colSums(re_ab_127.t==0) <= 381*0.8] -> re_ab_127.t.20per

sample.127.388v$Study.ID
sample.ID.Month[rownames(sample.ID.Month) %in% rownames(re_ab_127.t.20per),] %>% 
  cbind(re_ab_127.t.20per) %>% 
  dplyr::select(-Month) %>% 
  right_join(sample.127.388v, by=c("ID" = "Study.ID"))  -> re_ab_127_20per.clinc

save(re_ab_127_20per.clinc, file = "re_ab_127_20per.clinc.rda")

sample.ID.Month[rownames(sample.ID.Month) %in% rownames(re_ab_127.t.2363),] %>% 
  cbind(re_ab_127.t.2363) %>% 
  dplyr::select(-Month) %>% 
  right_join(sample.127.388v, by=c("ID" = "Study.ID"))  -> re_ab_127_2363.clinc


save(re_ab_127_2363.clinc, file = "re_ab_127_2363.clinc.rda")

sample.ID.Month[rownames(sample.ID.Month) %in% rownames(re_ab_127.t.20per),] %>% 
  cbind(re_ab_127.t.20per) %>% mutate(ID = as.numeric(ID)) %>% 
  dplyr::select(-Month) %>% 
  right_join(sample.127.388v, join_by(ID == Study.ID))  -> re_ab_127_20per.clinc

#t(re_ab_127_20per.clinc) %>% na.omit %>% t() -> re_ab_127_20per.clinc

sample.127.388v %>% dim()
re_ab_127_20per.clinc %>% dim()

re_ab_127_20per.clinc[,colnames(re_ab_127_20per.clinc) %in% LiverRelatedVar.1] -> sample.127.388v.liver

re_ab_127.t.20per %>% cbind(sample.127.388v.liver)  -> re_ab_127_20per.liver
re_ab_127_20per.liver

re_ab_127_20per.liver%>% 
  cor(use="complete.obs") -> per20_cor_liver

per20_cor_liver %>% reshape2::melt() %>% filter(value !=1) %>% 
  dplyr::filter(grepl("MGS", Var1)) %>% arrange(desc(abs(value))) %>%
  filter(!grepl("MGS", Var2)) %>% head(n = 5)

re_ab_127.t.20per %>% cbind(sample.127.388v.liver)  %>%  as.data.frame() %>%
  dplyr::select(c(contains("hMGS.05548"), contains("hMGS.04310"))) %>% setNames(c("species1","species2")) -> tmp.df
cor.test(tmp.df$species1, tmp.df$species2)

#######################################################
#####
#### Correlation of the liver variable change vs. 20percent 
####
########################################################
sample.127.388v.liver
sample.127.388v.liver %>% as.data.frame %>% mutate(ALT1=`B.ALT.(ukat/L)`- `A2.ALT.(ukat/L)`,
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
                                
                                
) %>% dplyr::select(-c(1:25)) -> sample.127.388.liver.delta

re_ab_127_20per %>% cbind(sample.127.388.liver.delta)  %>%
  cor(use="complete.obs") -> per20_cor_liver.delta

per20_cor_liver.delta %>% reshape2::melt() %>% filter(value !=1) %>% arrange(desc(abs(value))) %>%
  dplyr::filter(grepl("MGS", Var1))  %>%
  filter(!grepl("MGS", Var2))  %>% head(n=45) %>% dplyr::select(Var1) %>% distinct()



pdf("per20_cor_liver_delta.pdf", width=10, height=12)
per20_cor_liver.delta %>% reshape2::melt() %>% filter(value !=1) %>% arrange(desc(abs(value))) %>%
  dplyr::filter(Var1 %in% species24var)  %>%
  filter(!grepl("MGS", Var2))  %>%
  ggplot(aes(x=Var1, y=Var2, color = value)) +
  scale_color_gradient2(midpoint=0, 
                        low="red", 
                        mid="white",
                        high="blue", 
                        space ="Lab" )+
  theme(axis.text.x = element_text(angle = 90)) +
  geom_point(size=11, shape=22) +
  geom_point(aes( alpha=abs(value)*100, size=5+abs(value)*100)) +
  clean_background
dev.off()


############## association with Fat Fraction #####


LiverFatFra %>% separate_wider_delim(patient_id, delim = "_",
                                     names=c("prefix","ID")) %>%
  mutate(ID = as.numeric(ID)) -> LiverFatFra.1



ID127.new %>% filter(! Studie.ID %in% LiverFatFra.1$ID) %>% 
  dplyr::select(nID)

sam127 %>% dplyr::select(ID, Group) %>% distinct() -> sam127.tmp
LiverFatFra.1 %>% filter(ID %in% sam127$ID) -> LiverFatFra.2

LiverFatFra.2 %>% dplyr::select(-c(date_v1, date_v2, prefix)) %>%
  mutate(deltaFatFrac = fatfraction_v2-fatfraction_v1) %>% round(3) -> LiverFatFra.3

re_ab_127.1 %>% mutate(ID = as.numeric(ID)) %>% 
  dplyr::select(-Month) %>%
  left_join(LiverFatFra.3, join_by(ID == ID)) %>% dplyr::select(-ID) %>% as.matrix() -> re_ab_127.LiverFatFrac

re_ab_127.LiverFatFrac %>% cor(use="complete.obs") -> re_ab_127_20per.fatFrac.cor
re_ab_127.LiverFatFrac %>% cor.test()
cor.test(1:5,2:6)

re_ab_127.LiverFatFrac %>% cor() -> re_ab_127_20per.fatFrac.cor

re_ab_127_20per.fatFrac.cor %>% reshape2::melt() %>% 
  filter(grepl("MGS", Var1)) %>% filter(!grepl("MGS", Var2)) %>%
  filter(! is.na(value)) %>%
  arrange(desc(abs(value))) %>% head()

re_ab_127  %>% dplyr::select(ends_with("M0")) -> re_ab_127.M0
re_ab_127  %>% dplyr::select(ends_with("M6")) -> re_ab_127.M6
re_ab_127  %>% dplyr::select(ends_with("M12")) -> re_ab_127.M12

re_ab_127_deltaM120 <- re_ab_127.M12 - re_ab_127.M0 %>% round(3)



re_ab_127_deltaM120 %>% t() %>% as.data.frame() %>% rownames_to_column() %>%
  separate_wider_delim(cols = rowname, delim = "-",names = c("cID", "Month")) %>%
  separate_wider_delim(cols = cID, delim = "_",names = c("Prefix", "ID")) %>%
  mutate(ID=as.numeric(ID)) %>%
  left_join(LiverFatFra.3, join_by(ID==ID)) %>% dplyr::select(-c(Prefix, ID, Month)) %>%
  as.matrix() -> re_ab_127_deltaM120.fatFrac

re_ab_127_deltaM120.fatFrac%>% cor(use="complete.obs") %>% reshape2::melt() %>%
  filter(grepl("MGS", Var1)) %>% filter(!grepl("MGS", Var2)) %>%
  arrange(desc(abs(value))) %>% filter(abs(value) > 0.4) %>% dplyr::select(Var1) %>% distinct()


raw_tax_matrix_127  
species_04_names  <- rev(species_04_names)
  
re_ab_127_deltaM120.fatFrac %>% as.data.frame() %>% dplyr::select(c(all_of(species_04_names),contains("fat"))) -> big04

re_ab_127.LiverFatFrac %>% na.omit() %>% cor() -> per20_M120_FatFrac_cor
re_ab_127.LiverFatFrac %>% cor(use="complete.obs") -> per20_M120_FatFrac_cor

big04 %>% cor(use="complete.obs") -> big04_cor

pdf("try.pdf", width=5, height=5)
big04_cor %>% reshape2::melt() %>% dplyr::filter(grepl("MGS",Var1))  %>%
  filter(grepl( "delta", Var2)) %>% na.omit() %>% arrange(desc((value))) %>%
  ggplot(aes(y=Var1, x=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient2(midpoint=0, 
                        low="red", 
                        mid="white",
                        high="blue", 
                        space ="Lab", 
                       breaks = c(-0.4,-0.2,0,0.2,0.4))+
  theme(axis.text.x = element_text(angle = 90)) +
  clean_background
dev.off()
heatmap.2(big04_cor, scale = "none", col = bluered(500), 
          trace = "none", density.info = "none")

######## lower correlation ###### 
per20_cor_liver %>% reshape2::melt() %>% filter(value !=1) %>% 
  dplyr::filter(Var1 %in% LiverRelatedVar.1) %>% arrange(desc(abs(value))) %>%
  filter(grepl("MGS", Var2)) %>% head(n = 5)
re_ab_127_20per %>% cbind(sample.127.388.liver)  %>%  
  dplyr::select(c(contains("A2.ALT"), contains("AF33"))) %>% setNames(c("ALT","AF33")) -> tmp.df
cor.test(tmp.df$ALT, tmp.df$AF33)

dev.off()


heatmap.2(per20_cor_liver, scale = "none", col = bluered(500), 
          trace = "none", density.info = "none")

res_20per_rcorr %>% arrange(desc(abs(cor))) %>% dplyr::filter(! grepl("MGS", row))

pdf("global_20per_cor.pdf", width=12)


re_ab_127_20per.clinc %>% cor() -> per20_cor
per20_cor %>% dim()
heatmap.2(per20_cor, scale = "none", col = bluered(500), 
          trace = "none", density.info = "none")
dev.off()

re_ab_127_20per.clinc %>% rcorr() -> global_20per_rcorr

flattenCorrMatrix(global_20per_rcorr$r, global_20per_rcorr$P) -> res_20per_rcorr

res_cor.var <- data.frame(matrix(nrow=0, ncol=4))
res_cor.var %>% dim()
colnames(res_cor.var) <- colnames(res_20per_rcorr)


res_cor.var %>% dplyr::select(row) %>% distinct()

#######################################################
#####
#### Correlation of the liver variable with 20percent 
####
########################################################

sample.127.388v %>% dplyr::select(c(Study.ID, any_of(LiverRelatedVar.1))) -> 
  sample.127.388.liver

t(re_ab_127.t.20per) %>% na.omit %>% t() -> re_ab_127_20per

re_ab_127_20per %>% dim()
sample.127.388.liver 
sample.ID.Month[rownames(sample.ID.Month) %in% rownames(re_ab_127.t.20per),] %>% 
  cbind(re_ab_127.t.20per) %>% mutate(ID = as.numeric(ID)) %>% 
  dplyr::select(-Month) %>%
  left_join(sample.127.388.liver, join_by(ID == Study.ID)) %>% 
  as.matrix() -> re_ab_127_20per.clinc


re_ab_127_20per.clinc %>% cor() -> per20_cor
heatmap.2(per20_cor, scale = "none", col = bluered(500), 
          trace = "none", density.info = "none")
dev.off()




 res_cor.var.2 %>% arrange(desc(abs(cor)))

############### Global ###################
re_ab_127.clinc %>% rcorr() -> 
  global_rcorr

re_ab_127.clinc %>% dim()

t(re_ab_127.clinc) %>% na.omit %>% t() -> re_ab_127.clinc

re_ab_127.clinc %>% cor() -> 
  global_cor

global_cor[global_cor > 0.7]

flattenCorrMatrix(global_rcorr$r,   global_rcorr$P)  %>% filter(! abs(cor) == Inf) %>%
  dplyr::arrange(desc(abs(cor))) %>% filter(abs(cor) > 0.4) 


library(ComplexHeatmap)

library("gplots")

pdf("global_cor.pdf")
heatmap.2(global_cor, scale = "none", col = bluered(500), 
          trace = "none", density.info = "none")

dev.off()




heatmap(global_cor)

res_cor <- rcorr(as.matrix(sample.127.388v))

res_cor$r %>% head()




flattenCorrMatrix(res_cor$r, res_cor$P) %>% 
  filter( ! (cor == 1)) %>% filter( ! cor == -1 ) %>% filter(! abs(cor) == Inf) %>%
  dplyr::arrange(desc(abs(cor))) %>% filter(abs(cor) > 0.4) -> res_cor.2

res_cor.2 %>% colnames()

res_cor.var <- data.frame(matrix(nrow=0, ncol=4))
res_cor.var %>% dim()
colnames(res_cor.var) <- colnames(res_cor.2)

for (var in vars){
res_cor.2[grep(var, res_cor.2$row),] %>% rbind(res_cor.var) -> res_cor.var
}

res_cor.var


fat_results <- data.frame(matrix(nrow=0, ncol=5))
fat_results
colnames(fat_results) <- c("Group","Month","mean","sd","var")

for(var in vars){
  sample.127.full %>% dplyr::select(c(contains(var), diet.group)) %>% group_by(diet.group) 
  
}

re_ab_127
