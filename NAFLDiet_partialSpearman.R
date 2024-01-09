source("NAFLDiet_sample.R")




library(ppcor)


LiverRelatedVar.1 <- c("Study.ID",
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




LiverFatFra %>% separate_wider_delim(patient_id, delim = "_",
                                     names=c("prefix","ID")) %>%
  mutate(ID = as.numeric(ID)) -> LiverFatFra.1



ID127.new %>% filter(! Studie.ID %in% LiverFatFra.1$ID) %>% 
  dplyr::select(nID)

sam127 %>% dplyr::select(ID, Group) %>% distinct() -> sam127.tmp
LiverFatFra.1 %>% filter(ID %in% sam127$ID) -> LiverFatFra.2

sam127 %>% dplyr::filter(Month=="M0") %>% mutate(ID=as.numeric(ID)) %>%
  dplyr::left_join( LiverFatFra.2) %>% mutate(ID=as.character(ID)) -> 
  LiverFatFra.3

relative_abundance_MGS_data %>% dplyr::select(rownames(sam127)) -> re_ab_127.0
relative_abundance_MGS_data %>% dplyr::select(MGS, MainTax, rownames(sam127)) -> re_ab_127.mgs

re_ab_127.mgs %>% mutate(rn = paste0(MGS, ": ", MainTax)) %>% column_to_rownames(var="rn") %>%
  dplyr::select(-c(MGS, MainTax)) -> re_ab_127

## without changes the species names ## don't run
re_ab_127 <- re_ab_127.0
##########


rownames(re_ab_127)

t(re_ab_127) -> re_ab_127.t

re_ab_127.t[,colSums(re_ab_127.t) !=0] -> re_ab_127.t.2363


sample.ID.Month[rownames(sample.ID.Month) %in% rownames(re_ab_127.t.2363),] %>% 
  cbind(re_ab_127.t.2363) -> re_ab_127.1

sam127 %>% dplyr::select(ID, Group) -> sam127.ID.Group
re_ab_127.1 %>% right_join(sam127.ID.Group) -> re_ab_127.2
library(lme4)
library(broom)
library(tidyr)

### filter 20% ###
re_ab_127.2 %>% setNames(paste0("sp",names(.))) -> re_ab_127.3
re_ab_127.3[,colSums(re_ab_127.3==0) <= 381*0.8] -> re_ab_127.3.20per

re_ab_127.3.20per %>%  filter(spMonth == "M0") %>% gather(key, value, -spGroup) %>% group_by(key) %>%
  do(tidy(kruskal.test(x= .$value, g = .$spGroup))) -> kruskal.test.res.20per.M0.df

kruskal.test.res.M0.df %>% filter(p.value < 0.05) %>% dim()
kruskal.test.res.20per.M0.df %>% filter(p.value < 0.05) %>% dim()

re_ab_127.3.20per %>%  filter(spMonth == "M12") %>% gather(key, value, -spGroup) %>% group_by(key) %>%
  do(tidy(kruskal.test(x= .$value, g = .$spGroup))) -> kruskal.test.res.20per.M12.df

kruskal.test.res.20per.M12.df %>% filter(p.value < 0.05) %>% dim()

kruskal.test.res.20per.M0.df %>% filter(p.value < 0.05) %>% dplyr::select(key) -> M0.Kruskal.changedList
kruskal.test.res.20per.M12.df %>% filter(p.value < 0.05) %>% dplyr::select(key) -> M12.Kruskal.changedList

dplyr::intersect(M0.Kruskal.changedList, M12.Kruskal.changedList)

save(M12.Kruskal.changedList, file="M12.Kruskal.changedList.rda")
save(M0.Kruskal.changedList, file="M0.Kruskal.changedList.rda")

### filter 20% ###
re_ab_127.t[,colSums(re_ab_127.t==0) <= 381*0.8] -> re_ab_127.t.20per

sample.127.full %>% str()

sample.ID.Month[rownames(sample.ID.Month) %in% rownames(re_ab_127.t.20per),] %>% 
  cbind(re_ab_127.t.20per)  %>% 
  dplyr::select(-Month) %>% 
  right_join(sample.127.full, by=c("ID" = "Study.ID"))  -> re_ab_127_20per.clinc


re_ab_127_20per.clinc
re_ab_127_20per.clinc[,colnames(re_ab_127_20per.clinc) %in% LiverRelatedVar.1] -> sample.127.388v.liver

sample.127.388v.liver %>% as.data.frame %>% mutate( BMI1=`B.BMI.(kg/m2)`- `A2.BMI.(kg/m2)`,
                                                    BMI2=`C2.BMI.(kg/m2)`- `B.BMI.(kg/m2)`,
                                                    BMI3=`C2.BMI.(kg/m2)`- `A2.BMI.(kg/m2)`,
                                                   ALT1=`B.ALT.(ukat/L)`- `A2.ALT.(ukat/L)`,
                                                   ALT2=`C2.ALT.(ukat/L)`- `B.ALT.(ukat/L)`,
                                                   ALT3=`C2.ALT.(ukat/L)`- `A2.ALT.(ukat/L)`,
                                                   ALT1=`B.ALT.(ukat/L)`- `A2.ALT.(ukat/L)`,
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
                                                   
                                                   
) %>% dplyr::select(-c(2:25)) -> sample.127.388.liver.delta

save(sample.127.388.liver.delta, file="sample.127.388.liver.delta.rda")

relative_abundance_MGS_data %>% dplyr::select(rownames(sam127)) -> re_ab_127.0

t(re_ab_127.0) -> re_ab_127.0.t

re_ab_127.0
re_ab_127.0.t[,colSums(re_ab_127.0.t) !=0] -> re_ab_127.0.t.2363

re_ab_127.0.t.2363 %>% cbind(sam127) -> re_ab_127.0.t_sam127

sample.127.full[colSums(is.na(sample.127.full)) < 1] -> sample.127.388v

sample.127.388v[,colnames(sample.127.388v) %in% LiverRelatedVar.1]
sample.127.388v
sample.127.388v[,colnames(sample.127.388v) %in% LiverRelatedVar.1] %>%
  mutate(Study.ID = as.character(Study.ID)) %>% 
  left_join(re_ab_127.0.t_sam127, by=c("Study.ID" = "ID")) -> re_ab_127.liver

re_ab_127.liver  %>% cbind(sample.127.388.liver.delta) %>% na.omit() -> re_ab_127.liver.clean
LiverFatFra.3 %>% mutate(deltaFraction = fatfraction_v2 - fatfraction_v1) %>%
  dplyr::select(ID, fatfraction_v1, fatfraction_v2, deltaFraction) -> 
  LiverFatFra.4

save(LiverFatFra.4, file = "LiverFatFra.4.rda")

re_ab_127.liver %>%  left_join(LiverFatFra.4, by=c("Study.ID" = "ID")) -> 
  re_ab_127.liver.1

re_ab_127.liver.1  %>% cbind(sample.127.388.liver.delta) -> re_ab_127.liver.delta




spearman.cor.liver.delta

#############
############ Association between the values
############# 

re_ab_127.liver.delta %>% select_if(is.numeric) -> re_ab_127.liver.delta.numeric
 
re_ab_127.liver.delta.numeric  

shapiro.test(re_ab_127.liver.delta.numeric$`17`)

library("psych")
#A straightforward application of matrix algebra to 
#remove the effect of the variables in the y set from the x set. 
re_ab_127.liver.delta.numeric %>% dplyr::select(2385) %>% head()
partial.r(re_ab_127.liver.delta.numeric, 
          c(2:2387),
          c(1,2385,2386,2387,2388),
          use="pairwise",
          method="spearman") -> partial.r.result

re_ab_127.liver.delta.numeric 
partial.r.result %>% reshape2::melt() %>% 
  dplyr::filter(Var1 %in% LiverRelatedVar.1) %>%
  dplyr::filter(! Var2 %in% LiverRelatedVar.1) %>%
  arrange(desc(abs(value)))

sample_NAFLDiet %>% colnames()

############### Association delta M12 vs M0 

re_ab_127  %>% dplyr::select(ends_with("M0")) -> re_ab_127.M0
re_ab_127  %>% dplyr::select(ends_with("M6")) -> re_ab_127.M6
re_ab_127  %>% dplyr::select(ends_with("M12")) -> re_ab_127.M12

re_ab_127_deltaM120 <- re_ab_127.M12 - re_ab_127.M0 %>% round(3)

cor_spear.par <- 0
for (i in 25:2387){
  
 cor_spear.par[i] <- pcor.test(re_ab_127.liver.delta.numeric$`A2.ALT.(ukat/L)`, 
          re_ab_127.liver.delta.numeric[i], 
           re_ab_127.liver.delta.numeric$`A2.BMI.(kg/m2)`,
          method="spearman")
}

cor_spear.par %>% unlist()
