source("NAFLDiet_sample.R")


##########################################
##############################
############################## gain and disappearance
#############################
#############################################

relative_abundance_MGS_data[,rownames(sam127)] -> re_ab_127


re_ab_127  %>% dplyr::select(ends_with("M0")) -> re_ab_127.M0
re_ab_127  %>% dplyr::select(ends_with("M6")) -> re_ab_127.M6
re_ab_127  %>% dplyr::select(ends_with("M12")) -> re_ab_127.M12




re_ab_127[rowSums(re_ab_127) != 0,] -> re_ab_127.2363

re_ab_127.2363  %>% dplyr::select(ends_with("M0")) -> re_ab_127.M0
re_ab_127.2363  %>% dplyr::select(ends_with("M6")) -> re_ab_127.M6
re_ab_127.2363  %>% dplyr::select(ends_with("M12")) -> re_ab_127.M12


pdf("matrix01.pdf")
pdf("try.pdf")
set.seed(1)
sample(rep(c(0,1),32), replace = TRUE) %>% matrix(nrow=8) %>% melt() -> tmp.matrix.draw
tmp.matrix.draw %>% ggplot(aes(x=Var1, y=Var2))+
  geom_point(aes(color=factor(value), alpha=0.3), size=12) +
  geom_text(aes(label = value), size=6)+
  scale_color_manual(values = c("red","blue"))+
  clean_background
dev.off()


re_ab_127.M0 %>% 
  mutate_if(is.numeric, ~ 1 * (. != 0)) -> re_ab_127.M0.01
re_ab_127.M6 %>% 
  mutate_if(is.numeric, ~ 1 * (. !=  0)) -> re_ab_127.M6.01
re_ab_127.M12 %>% 
  mutate_if(is.numeric, ~ 1 * (. !=  0)) -> re_ab_127.M12.01


re_ab_127.M6.01 - re_ab_127.M0.01  -> deltaM60.01
sam127

deltaM60.01 %>% t() %>% as.data.frame %>% 
  rownames_to_column() %>%  
  separate_wider_delim(cols = rowname, delim = "-",names = c("nID", "Month")) %>%
  separate_wider_delim(cols = nID, delim = "_",names = c("prefix", "ID")) %>%
  dplyr::select(-c(prefix, Month)) %>% column_to_rownames(var="ID") ->
  deltaM60.01_v1  

re_ab_127.M0.01 %>% rowSums() %>% as.data.frame() %>%
  mutate_if(is.numeric, ~1 * (. != 0)) %>% colSums()

re_ab_127.M6.01 %>% rowSums() %>% as.data.frame() %>%
  mutate_if(is.numeric, ~1 * (. != 0)) %>% colSums()

re_ab_127.M12.01 %>% rowSums() %>% as.data.frame() %>%
  mutate_if(is.numeric, ~1 * (. != 0)) %>% colSums()

data.frame(ID=sam127$ID, Group = sam127$Group) %>% distinct() -> tmp.dict
deltaM60.01_v1 %>% dim()

### unchanged 487
deltaM60.01_v1[1:20,1:10] -> tmp
deltaM60.01_v1 %>% select(where(~ all(. == 0))) %>% dim()

deltaM60.01_v1 %>% select(- where(~ all(. == 0))) %>% dim()

deltaM60.01_v1 %>% select( where(~ any(. == 1))) %>% 
  select(where(~any(. == -1))) %>% dim()

deltaM60.01_v1 %>% select( where(~ any(. == 1))) %>% 
  select(-where(~any(. == -1))) %>% dim()


t(re_ab_127.M0.01)[,colSums(t(re_ab_127.M0.01))!=0] %>% dim()
t(re_ab_127.M6.01)[,colSums(t(re_ab_127.M6.01))!=0] %>% dim()
t(re_ab_127.M12.01)[,colSums(t(re_ab_127.M12.01))!=0] %>% dim()

## M0 to M6
t(re_ab_127.M0.01) %>% 
  as.data.frame() -> t.M0.01

t.M0.01[,colSums(t.M0.01) == 0] %>% colnames() -> tmp.name
t.M0.01[,colSums(t.M0.01) != 0] %>% colnames() -> tmp.lose.name

t(re_ab_127.M6.01) %>% as.data.frame() -> t.M6.01

t.M6.01[,tmp.name] -> t.M6.01.tmp
t.M6.01[,tmp.lose.name] -> t.M6.01.tmp.lose
t.M6.01.tmp[,colSums(t.M6.01.tmp) !=0] %>% dim()
t.M6.01.tmp.lose[,colSums(t.M6.01.tmp.lose) == 0] %>% dim()

## M0 to M12
t(re_ab_127.M12.01) %>% as.data.frame() -> t.M12.01

t.M12.01[,tmp.name] -> t.M12.01.tmp
t.M12.01[,tmp.lose.name] -> t.M12.01.tmp.lose
t.M12.01.tmp[,colSums(t.M12.01.tmp) !=0] %>% dim()
t.M12.01.tmp.lose[,colSums(t.M12.01.tmp.lose) == 0] %>% dim()

t.M6.01.tmp[,colSums(t.M6.01.tmp) !=0] %>% colnames() -> tmp3.name
t.M12.01[,tmp3.name] -> t.M12.01.tmp3
t.M12.01.tmp3[colSums( t.M12.01.tmp3) == 0] %>% dim()

## M6 to M12 
t(re_ab_127.M6.01) %>% 
  as.data.frame() -> t.M6.01
t.M6.01[,colSums(t.M6.01) == 0] %>% colnames() -> tmp2.name
t.M6.01[,colSums(t.M6.01) != 0] %>% colnames() -> tmp2.lose.name
t(re_ab_127.M12.01) %>% as.data.frame() -> t.M12.01

t.M12.01[,tmp2.name] -> t.M12.01.tmp
t.M12.01[,tmp2.lose.name] -> t.M12.01.tmp.lose
t.M12.01.tmp[,colSums(t.M12.01.tmp) !=0] %>% dim()
t.M12.01.tmp.lose[,colSums(t.M12.01.tmp.lose) == 0] %>% dim()

t.M12.01.tmp[,colSums(t.M12.01.tmp) !=0] %>% colnames() -> tmp161.name

t.M6.01.tmp[,colSums(t.M6.01.tmp) !=0] %>% colnames() -> tmp179.name

table(colSums(t.M12.01[,tmp179.name]) != 0)

t.M6.01.tmp[,colSums(t.M6.01.tmp) !=0] -> tmp179 

t.M6.01.tmp.lose[,colSums(t.M6.01.tmp.lose) == 0] %>% rownames() -> tmp163.loss.name
t.M12.01.tmp[,colSums(t.M12.01.tmp) !=0] %>% rownames() -> tmp188.name

table(colSums(t.M0.01[,tmp161.name])==0 )

write.csv(tmp161.name, 'try1.csv')
write.csv(tmp179.name, 'try2.csv')
write.csv(tmp163.loss.name, 'try3.csv')
write.csv(tmp188.name, 'try4.csv')

## gain and lose, 1169
deltaM60.01_v1 %>% dplyr::select(where(~any(. ==1 ))) %>% dplyr::select(where(~any(. == -1 ))) %>% dim()

## loss but gained in other individual 1169
deltaM60.01_v1 %>% dplyr::select(where(~any(. < 0 ))) %>% dplyr::select( where(~any(. == 1 ))) %>% dim()

## unchanged 729
deltaM60.01_v1[,colSums(deltaM60.01_v1) == 0] %>% dim()

## complete gain, 387
deltaM126.01_v1 %>% dplyr::select(where(~any(. > 0 ))) %>% dplyr::select(- where(~any(. == -1 ))) %>% dim()
## gain and lose, 1133
deltaM126.01_v1 %>% dplyr::select(where(~any(. > 0 ))) %>% dplyr::select(where(~any(. == -1 ))) %>% dim()

## complete loss 331
deltaM126.01_v1 %>% dplyr::select(where(~any(. < 0 ))) %>% dplyr::select(- where(~any(. == 1 ))) %>% dim()
## loss but gained in other individual 1133
deltaM126.01_v1 %>% dplyr::select(where(~any(. < 0 ))) %>% dplyr::select( where(~any(. == 1 ))) %>% dim()

# unchaged 787
deltaM126.01_v1[,colSums(deltaM126.01_v1) == 0] %>% dim()

deltaM60.01_v1 %>% mutate(
  no_zeros = rowSums(. == 0),
  no_one = rowSums(. == 1),
  no_minusOne = -rowSums(. == -1),
  profit = no_one + no_minusOne
) %>% cbind(tmp.dict) -> deltaM60.01_v2

deltaM60.01_v2 %>% dplyr::select(-c(1:2363))

t(re_ab_127.M0.01) -> t01
t(re_ab_127.M6.01) -> t06
t(re_ab_127.M12.01) -> t012

t01 %>% rowSums() %>% as.data.frame() %>% head()
t06 %>% rowSums() %>% as.data.frame() %>% head()
t012 %>% rowSums() %>% as.data.frame() %>% head()

t06 %>% rowSums() %>% as.data.frame() %>% arrange(desc(.)) 

pdf("try.pdf")
t01 %>% rowSums() %>% melt() %>% ggplot(aes(x=value))+
  geom_histogram()
t06 %>% rowSums() %>% melt() %>% ggplot(aes(x=value))+
  geom_histogram()
t012 %>% rowSums() %>% melt() %>% ggplot(aes(x=value))+
  geom_histogram()
dev.off()


sam127 %>% filter(Month=="M0") %>% distinct()  -> tmp.dict.1

t01 %>% rowSums() %>% melt() %>% cbind(tmp.dict.1) %>%
  dplyr::select(ID, Group, value) %>% group_by(Group) %>%
  summarise(meanM0 = mean(value), sdM0 = sd(value)) -> 
  t0.summary

t06 %>% rowSums() %>% melt() %>% cbind(tmp.dict.1) %>%
  dplyr::select(ID, Group, value) %>% group_by(Group) %>%
  summarise(meanM6 = mean(value), sdM6 = sd(value)) -> 
  t6.summary

t012 %>% rowSums() %>% melt() %>% cbind(tmp.dict.1) %>%
  dplyr::select(ID, Group, value) %>% group_by(Group) %>%
  summarise(meanM12 = mean(value), sdM12 = sd(value)) -> 
  t12.summary

deltaM60.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% group_by(Group) %>%
  summarise(mean_1_M60=mean(no_one),
            sd_1_M60=sd(no_one),
            mean_minus1_M60 =  mean(no_minusOne),
            sd_minus1_M60 =  sd(no_minusOne)
  ) -> t60.summary

t60.summary


deltaM60.01_v2 %>% dplyr::select(-c(1:2363)) %>%
  cbind(tmp.dict.1) -> deltaM60.01_annova
t01 %>% rowSums() %>% melt() %>% cbind(deltaM60.01_annova) -> 
  deltaM60.01_annova 

######## Adjust for sex, age, BMI, and baseline 

ancova_model <- aov(no_one ~ Group + value + age + BMI + gender + Diabetes 
                    + Metformin, 
                    data=deltaM60.01_annova)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

ancova_model <- aov(no_minusOne ~ Group + value + age + BMI + gender + Diabetes 
                    + Metformin, 
                    data=deltaM60.01_annova)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)


deltaM60.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(no_minusOne)) %>% group_by(Group) 




pdf("try.pdf", width=20, height = 1.5)
deltaM60.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% group_by(Group) %>%
  summarise(mean_1=mean(no_one),
            mean_minus1 =  mean(no_minusOne)) %>%
  group_by(Group) %>%
  ggplot(aes(x=1, y=mean_1))+
  geom_col(aes(col="white", fill="red"), width=0.3) +
  geom_col(aes(x=1, y=mean_minus1, fill="blue"), width=0.3)+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray", size=1)+
  coord_flip() +
  facet_wrap(~Group) +
  clean_background

dev.off()

deltaM60.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange((no_minusOne)) %>% filter(Group == 1) 



pdf("Group1_lossGainSpecies.pdf")
deltaM60.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% filter(Group == 1) %>%
  ggplot(aes(x=1:43, y=no_one, fill="white"))+
  geom_col(aes(fill="red", alpha=0.5)) +
  geom_col(aes(x=1:43, y=no_minusOne, fill="blue", alpha=0.5))+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray", size=1.2) +
  geom_smooth(aes(x=1:43, y=profit),
              color="orange",fill="orange", span=1,
              level=0.999)+
  clean_background +
  xlim(0,47) +
  ylim(-230,100)+
  coord_flip()
dev.off()
deltaM60.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% filter(Group == 2) 

deltaM60.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% filter(Group == 2) 

pdf("Group2_lossGainSpecies.pdf")
pdf("try.pdf")
deltaM60.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% filter(Group == 2) %>%
  rowid_to_column() %>%
  ggplot(aes(x=rowid, y=no_one, fill="white"))+
  geom_col(aes(fill= "red",alpha=0.5)) +
  geom_col(aes(y=no_minusOne, fill="blue", alpha=0.5))+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray", size=1.2) +
  geom_smooth(aes( y=profit),
              color="orange",fill="orange", span=1,
              level=0.999)+
  clean_background +
  coord_flip() +
  xlim(0,47) +
  ylim(-230,100)
dev.off()

deltaM60.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% filter(Group == 3) %>% dim()

pdf("Group3_lossGainSpecies.pdf")
deltaM60.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% filter(Group == 3) %>%
  ggplot(aes(x=1:38, y=no_one, fill="white"))+
  geom_col(aes( fill="red", alpha=0.5)) +
  geom_col(aes(y=no_minusOne, fill="blue", alpha=0.5))+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray", size=1) +
  geom_smooth(aes(y=profit),
              color="orange",fill="orange", span=1,
              level=0.999)+
  clean_background +
  xlim(0,47) +
  ylim(-230,100)+
  coord_flip()
dev.off()

deltaM60.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(no_one))


deltaM60.01_v2 %>% dplyr::select(-c(1:2363)) %>% group_by(Group) %>%
  summarise(mean(no_zeros), mean(no_one), mean(no_minusOne))

###################################
##############
############## M12 to M6 ########
#############
####################################

re_ab_127.M12.01 - re_ab_127.M6.01  -> deltaM126.01

deltaM126.01 %>% t() %>% as.data.frame %>% 
  rownames_to_column() %>%  
  separate_wider_delim(cols = rowname, delim = "-",names = c("nID", "Month")) %>%
  separate_wider_delim(cols = nID, delim = "_",names = c("prefix", "ID")) %>%
  dplyr::select(-c(prefix, Month)) %>% column_to_rownames(var="ID") ->
  deltaM126.01_v1  

deltaM126.01_v1  %>% dim()

data.frame(ID=sam127$ID, Group = sam127$Group) %>% distinct() -> tmp.dict

deltaM126.01_v1 %>% mutate(
  no_zeros = rowSums(. == 0),
  no_one = rowSums(. == 1),
  no_minusOne = -rowSums(. == -1),
  profit = no_one-no_minusOne
) %>% cbind(tmp.dict) -> deltaM126.01_v2

deltaM126.01_v2 %>% dplyr::select(-c(1:2363))

t(re_ab_127.M0.01) -> t01
t(re_ab_127.M6.01) -> t06
t(re_ab_127.M12.01) -> t012


t01 %>% rowSums() %>% as.data.frame() %>% head()
t06 %>% rowSums() %>% as.data.frame() %>% head()
t012 %>% rowSums() %>% as.data.frame() %>% head()

pdf("try.pdf")
t01 %>% rowSums() %>% melt() %>% ggplot(aes(x=value))+
  geom_histogram()
t06 %>% rowSums() %>% melt() %>% ggplot(aes(x=value))+
  geom_histogram()
t012 %>% rowSums() %>% melt() %>% ggplot(aes(x=value))+
  geom_histogram()
dev.off()



sam127 %>% dplyr::filter(Month=="M0") %>% dplyr::select(-Month) -> tmp.dict.1
deltaM126.01_v2 %>% dplyr::select(-c(1:2363)) %>%
  cbind(tmp.dict.1) -> deltaM126.01_annova
t01 %>% rowSums() %>% melt() %>% cbind(deltaM126.01_annova) -> 
  deltaM126.01_annova 

######## Adjust for sex, age, BMI, and baseline 

ancova_model <- aov(no_one ~ Group + value + age + BMI + gender + Diabetes 
                    + Metformin, 
                    data=deltaM126.01_annova)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

ancova_model <- aov(no_minusOne ~ Group + value + age + BMI + gender + Diabetes 
                    + Metformin, 
                    data=deltaM126.01_annova)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

deltaM126.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(no_minusOne)) %>% head()

deltaM126.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% group_by(Group) %>%
  summarise(mean_one_126=mean(no_one),
            sd_one_126=sd(no_one),
            mean_minus_1_126 =  mean(no_minusOne),
            sd_minus1_1_126 =  sd(no_minusOne)
  ) -> t126.summary

cbind(t0.summary, t6.summary) %>% 
  cbind(t12.summary) %>% cbind(t60.summary) %>% dplyr::select(-Group) %>% 
  cbind(t126.summary) %>% mutate(across(where(is.numeric), round, digits=1)) ->
  summary.mean



###########
clean_background <- theme(plot.background = element_rect("white"),
                          panel.background = element_rect("white"),
                          panel.grid = element_line("white"),
                          axis.line = element_line("black"),
                          axis.text = element_text(size = 14, color = "black"),
                          axis.title = element_text(color = "black"),
                          legend.text = element_text(size = 14),
                          legend.key = element_rect("white"))

pdf("try.pdf", width=20, height=5)
summary.mean  %>% group_by(Group) %>%
  ggplot(aes(x=1:3, color=Group)) +
  geom_segment(aes(x=0.8, xend=1.2, y=meanM0, yend=meanM0), linewidth = 3)+
  geom_segment(aes(x=1.8, xend=2.2, y=meanM6, yend=meanM6), linewidth = 3)+
  geom_segment(aes(x=2.8, xend=3.2, y=meanM12, yend=meanM12), linewidth = 3)+
  geom_segment(aes(x=1.2, xend=1.4, y=meanM0, yend=meanM0+mean_1_M60), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=1.2, xend=1.4, y=meanM0, yend=meanM0+mean_minus1_M60), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=1.6, xend=1.8, y=meanM0+mean_1_M60, yend=meanM6), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=1.6, xend=1.8, y=meanM0+mean_minus1_M60, yend=meanM6), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=1.4, xend=1.6, y=meanM0 + mean_1_M60, 
                   yend=meanM0 + mean_1_M60), color = "lightblue", 
               linewidth = 1.8)+
  geom_segment(aes(x=1.4, xend=1.6, y=meanM0 + mean_minus1_M60, 
                   yend=meanM0 + mean_minus1_M60), color = "pink", 
               linewidth = 1.8)+
  geom_segment(aes(x=2.4, xend=2.6, y=meanM6 + mean_one_126 , 
                   yend=meanM6 + mean_one_126), color = "lightblue", 
               linewidth = 1.8)+
  geom_segment(aes(x=2.4, xend=2.6, y=meanM6 + mean_minus_1_126 , 
                   yend=meanM6 + mean_minus_1_126), color="pink", 
               linewidth = 1.8)+
  geom_segment(aes(x=2.2, xend=2.4, y=meanM6, yend=meanM6+mean_one_126), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=2.2, xend=2.4, y=meanM6, yend=meanM6+mean_minus_1_126), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=2.6, xend=2.8, y=meanM6+mean_one_126, yend=meanM12), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=2.6, xend=2.8, y=meanM6+mean_minus_1_126, yend=meanM12), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_label(aes(x=1, y=meanM0+6, label=meanM0),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=2, y=meanM6+6, label=meanM6),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=3, y=meanM12+6, label=meanM12),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=1.5, y=meanM0+ mean_1_M60 + 6, label=mean_1_M60),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=1.5, y=meanM0+mean_minus1_M60+6, label=mean_minus1_M60),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=2.5, y=meanM6+mean_one_126+6, label=mean_one_126),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=2.5, y=meanM6+mean_minus_1_126+6, label=mean_minus_1_126),color="black", 
             label.size = NA, fontface = "bold")+
  clean_background +
  facet_grid(~Group)

dev.off()



gain_lose <- data.frame(
  meanM0 = 2091,
  meanM6 = 2107,          
  meanM12 = 2124,
  mean_1_M60 = 179  ,    
  mean_minus1_M60 = -163,  
  mean_one_126    =  161,
  mean_minus_1_126 = -144,
  mean_one_120 = 188,
  mean_minus_1_120 = -155
)


pdf("try.pdf")
gain_lose %>% ggplot(aes(x=1:3)) +
  geom_segment(aes(x=0.8, xend=1.2, y=meanM0, yend=meanM0), linewidth = 3, color="red")+
  geom_segment(aes(x=1.8, xend=2.2, y=meanM6, yend=meanM6), linewidth = 3, color="red")+
  geom_segment(aes(x=2.8, xend=3.2, y=meanM12, yend=meanM12), linewidth = 3, color="red")+
  geom_segment(aes(x=1.2, xend=1.4, y=meanM0, yend=meanM0+mean_1_M60), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=1.2, xend=1.4, y=meanM0, yend=meanM0+mean_minus1_M60), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=1.6, xend=1.8, y=meanM0+mean_1_M60, yend=meanM6), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=1.6, xend=1.8, y=meanM0+mean_minus1_M60, yend=meanM6), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=1.2, xend=1.9, y=meanM0, yend=meanM0 + mean_one_120), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=1.2, xend=1.9, y=meanM0, yend=meanM0 +mean_minus_1_120), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=1.9, xend=2.1, y=meanM0 + mean_one_120, 
                   yend=meanM0 + mean_one_120), color = "lightblue", 
               linewidth = 1.8)+
  geom_segment(aes(x=1.9, xend=2.1, y=meanM0 + mean_minus_1_120, 
                   yend=meanM0 + mean_minus_1_120), color = "pink", 
               linewidth = 1.8)+
  geom_segment(aes(x=1.4, xend=1.6, y=meanM0 + mean_1_M60, 
                   yend=meanM0 + mean_1_M60), color = "lightblue", 
               linewidth = 1.8)+
  geom_segment(aes(x=1.4, xend=1.6, y=meanM0 + mean_minus1_M60, 
                   yend=meanM0 + mean_minus1_M60), color = "pink", 
               linewidth = 1.8)+
  geom_segment(aes(x=2.4, xend=2.6, y=meanM6 + mean_one_126 , 
                   yend=meanM6 + mean_one_126), color = "lightblue", 
               linewidth = 1.8)+
  geom_segment(aes(x=2.4, xend=2.6, y=meanM6 + mean_minus_1_126 , 
                   yend=meanM6 + mean_minus_1_126), color="pink", 
               linewidth = 1.8)+
  geom_segment(aes(x=2.2, xend=2.4, y=meanM6, yend=meanM6+mean_one_126), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=2.2, xend=2.4, y=meanM6, yend=meanM6+mean_minus_1_126), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=2.6, xend=2.8, y=meanM6+mean_one_126, yend=meanM12), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=2.6, xend=2.8, y=meanM6+mean_minus_1_126, yend=meanM12), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=2.1, xend=2.8, y=meanM0+mean_one_120, yend=meanM12), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_segment(aes(x=2.1, xend=2.8, y=meanM0+mean_minus_1_120, yend=meanM12), 
               linewidth = 1, linetype = 3, color = "gray")+
  geom_label(aes(x=1, y=meanM0+10, label=meanM0),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=2, y=meanM6+10, label=meanM6),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=3, y=meanM12+10, label=meanM12),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=1.5, y=meanM0+ mean_1_M60 + 10, label=mean_1_M60),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=1.5, y=meanM0+mean_minus1_M60+10, label=mean_minus1_M60),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=2.5, y=meanM6+mean_one_126+10, label=mean_one_126),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=2.5, y=meanM6+mean_minus_1_126+10, label=mean_minus_1_126),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=2, y=meanM0+mean_one_120+10, label=mean_one_120),color="black", 
             label.size = NA, fontface = "bold")+
  geom_label(aes(x=2, y=meanM0+mean_minus_1_120+10, label=mean_minus_1_120),color="black", 
             label.size = NA, fontface = "bold")+
  clean_background

dev.off()


pdf("try.pdf", width=20, height = 1.5)
deltaM126.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% group_by(Group) %>%
  summarise(mean_1=mean(no_one),
            mean_minus1 =  mean(no_minusOne)) %>%
  group_by(Group) %>%
  ggplot(aes(x=1, y=mean_1))+
  geom_col(aes(col="white", fill="blue"), width=0.3) +
  geom_col(aes(x=1, y=-(mean_minus1), fill="red"), width=0.3)+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray", size=1)+
  coord_flip() +
  facet_wrap(~Group) +
  clean_background

dev.off()

deltaM126.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% filter(Group == 1) %>% dim()


pdf("Group1_lossGainSpecies.pdf")
deltaM126.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% filter(Group == 1) %>%
  ggplot(aes(x=1:43, y=no_one, fill="white"))+
  geom_col(aes(fill="blue", alpha=0.5)) +
  geom_col(aes(x=1:43, y=-(no_minusOne), fill="red", alpha=0.5))+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray", size=1) +
  geom_smooth(aes(x=1:43, y=profit),
              color="orange",fill="orange", span=1,
              level=0.999)+
  clean_background +
  xlim(0,45) +
  ylim(-180,130)+
  coord_flip()
dev.off()

pdf("Group2_lossGainSpecies.pdf")
pdf("try.pdf")
deltaM126.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% filter(Group == 2) %>%
  ggplot(aes(x=1:46, y=no_one))+
  geom_col(aes(col="white", fill="red", alpha=0.5)) +
  geom_col(aes(x=1:46, y=-(no_minusOne), fill="blue", alpha=0.5))+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray", size=1) +
  geom_smooth(aes(x=1:46, y=profit),
              color="orange",fill="orange", span=1,
              level=0.999)+
  clean_background +
  coord_flip() +
  xlim(0,47) +
  ylim(-180,130)
dev.off()

deltaM126.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% filter(Group == 3) %>% dim()

pdf("Group3_lossGainSpecies.pdf")
deltaM126.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(profit)) %>% filter(Group == 3) %>%
  ggplot(aes(x=1:38, y=no_one))+
  geom_col(aes(col="white", fill="blue", alpha=0.5)) +
  geom_col(aes(x=1:38, y=-(no_minusOne), fill="red", alpha=0.5))+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray", size=1) +
  geom_smooth(aes(x=1:38, y=profit),
              color="orange",fill="orange", span=1,
              level=0.999)+
  clean_background +
  xlim(0,45) +
  ylim(-180,130)+
  coord_flip()
dev.off()

deltaM126.01_v2 %>% dplyr::select(-c(1:2363)) %>% 
  arrange(desc(no_one))


deltaM126.01_v2 %>% dplyr::select(-c(1:2363)) %>% group_by(Group) %>%
  summarise(mean(no_zeros), mean(no_one), mean(no_minusOne))
