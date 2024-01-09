source("NAFLDiet_sample.R")
library(FSA)
library(ggsci)


################################
########                #######
####### alpha diversity  (simpson, shannon, richness and faith's phylogy tree)
######
######                  ######
################################

div <- microbiome::alpha(ex1b, index="all")
head(div)

### 
library(picante)
data("phylocom")
phylocom$phylo %>% str()
phylocom$sample
phy_tree(ex1b) %>% str()
sample_data(ex1b)
prunedTree <- prune.sample(sam127, )


######## Faither's Phylogenic Diversity
pd(t(otu_table(ex1b)), phy_tree(ex1b), include.root = FALSE) -> pd.result
sample.A2.127.2 %>% rownames()
sam127
pdf("try.pdf")
pd.result %>% cbind(sam127) %>% dplyr::select(PD, SR, Month, Group) %>%
  group_by(Group) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>%
  mutate(Group=as.factor(Group)) %>%
  ggplot(aes(x=Month, y = PD,  group=interaction(Month, Group), color=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.2, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), alpha=0.2) +
  scale_fill_nejm() +
  scale_color_nejm() +
  #scale_color_manual(values = c("red","green","blue")) + 
  #showSignificance(c(0.7,1.3),1,-0.5,"NS")+
   showSignificance(c(1.9,2.3),1400,-0.5,"*")+
  #showSignificance(c(2.7,3.3),1,-0.5,"NS")+
  clean_background +
  labs(x = "Month",
       y = "Faith's Phylogenic Diversity",
       title = "Faith's Phylogenic diversity at three time points") +
  #scale_fill_discrete(name = "Groups", labels  = c("1: LC-PUFA \n (n = 129)", "2: HND \n (n = 138)", "3: NNR-recommend \n (n = 114)")) +
  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16)
  )

dev.off()



cbind(pd.result, sam127) %>% filter(Month=="M0") -> M0.pd
cbind(pd.result, sam127) %>% filter(Month=="M6") -> M6.pd
cbind(pd.result, sam127) %>% filter(Month=="M12") -> M12.pd

cbind(pd.result, sam127) %>% filter(Group=="1") %>% mutate(Month=as.factor(Month)) -> G1.pd
cbind(pd.result, sam127) %>% filter(Group=="2") -> G2.pd
cbind(pd.result, sam127) %>% filter(Group=="3") -> G3.pd

shapiro.test(M0.pd$PD)
shapiro.test(M6.pd$PD)
shapiro.test(M12.pd$PD)
shapiro.test(pd.result$PD)

## only group 1 is normal, G2, G3 not normal
shapiro.test(G3.pd$PD)

kruskal.test(PD ~ Month, data=G2.pd)
kruskal.test(PD ~ Month, data=G3.pd)

ancova_model <- aov(PD ~ Month, data=G1.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Month="Tukey"))
summary(postHocs)

########## M0 normal, M12 normal, M6 not normal
kruskal.test(PD ~ Group, data=M0.pd)
kruskal.test(PD ~ Group, data=M6.pd)
kruskal.test(PD ~ Group, data=M12.pd)





pairwise.wilcox.test(M6.pd$PD, M6.pd$Group, p.adjust.method = "bonferroni")

dunnTest(PD ~ Group, data=M6.pd, method="bonferroni")


ancova_model <- aov(PD ~ Group,
                    data=M0.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

ancova_model <- aov(PD ~ Group + age + BMI + gender + Diabetes + Metformin,
                    data=M0.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

ancova_model <- aov(PD ~ Group + age,
                    data=M0.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

ancova_model <- aov(PD ~ Group ,
                    data=M12.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)


ancova_model <- aov(PD ~ Group + age + BMI + gender + Diabetes + Metformin ,
                    data=M12.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)



  
###################
########
######## Richness
#######
###################
pd.result

pdf("meanRichness.pdf", width=540, height=480)
pd.result %>% cbind(sam127) %>% dplyr::select(PD, SR, Month, Group) %>%
  group_by(Group, Month) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>%
  mutate(Group=as.factor(Group)) %>%
  ggplot(aes(x=Month, y=SR, group=interaction(Month, Group))) +
  stat_summary(fun.y=mean, geom="line", position=position_dodge(width=0.2), size=1,
             aes(group=Group, color=Group))+
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge(width=0.2), size=0.1,
               aes(group=Group, color=Group))+
  stat_summary(fun=mean, geom="point", position=position_dodge(width=0.2), size=2.5, shape=21,
             aes(group=Group, colour=Group))+
  scale_fill_nejm() +
  scale_color_nejm() +  clean_background 
  

pdf("try.pdf")
pd.result %>% cbind(sam127) %>% dplyr::select(PD, SR, Month, Group) %>%
  group_by(Group) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>%
  mutate(Group=as.factor(Group)) %>%
  ggplot(aes(x=Month, y = SR,  group=interaction(Month, Group), color=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.2, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), alpha=0.2) +
  scale_fill_nejm() +
  scale_color_nejm() +
  #showSignificance(c(0.7,1.3),470,-0.5,"NS")+
  #showSignificance(c(1.7,2),470,-0.5,"NS")+
  showSignificance(c(1.9,2.3),465,-0.5,"*")+
  #showSignificance(c(2.7,3.3),470,-0.5,"NS")+
  clean_background +
  labs(x = "Month",
       y = "Richness") +
  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16)
  )

dev.off()

#### G1, G2 normal, G3 not normal
shapiro.test(G1.pd$SR)
shapiro.test(G2.pd$SR)
shapiro.test(G3.pd$SR)

shapiro.test(pd.result$SR)
cbind(pd.result, sam127) -> pd.SR
kruskal.test(pd.SR$SR ~ pd.SR$Month, data=pd.SR)

########## all normal
ancova_model <- aov(SR ~ Group,data=M0.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

ancova_model <- aov(SR ~ Group + age + BMI + gender + Diabetes + Metformin, data=M0.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

ancova_model <- aov(SR ~ Group + age, data=M0.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

ancova_model <- aov(SR ~ Group ,data=M12.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)


ancova_model <- aov(SR ~ Group + age + BMI + gender + Diabetes + Metformin, data=M12.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

ancova_model <- aov(SR ~ Group + age, data=M12.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

ancova_model <- aov(SR ~ Group ,data=M6.pd)
summary(ancova_model)
TukeyHSD(ancova_model, conf.level = 0.95)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)
confint(ancova_model)

ancova_model <- aov(SR ~ Group + age + BMI + gender + Diabetes + Metformin, data=M6.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

confint(ancova_model)


ancova_model <- aov(RS ~ Group + age, data=M6.pd)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

pairwise.t.test(M6.pd$SR, M6.pd$Group, p.adjust.method="bonferroni")




############# simpson ###################

div_simpson = vegan::diversity(t(relative_abundance_MGS_data_127), index = "simpson") 

div_simpson %>% enframe() %>% setNames(c("name","simpson")) -> div_simpson

cbind(div_simpson, sam127) %>% filter(Month=="M0") -> M0.simpson
cbind(div_simpson, sam127) %>% filter(Month=="M6") -> M6.simpson
cbind(div_simpson, sam127) %>% filter(Month=="M12") -> M12.simpson

div_simpson %>% str()

div_simpson$simpson
mshapiro.test(t(div_simpson$simpson))
mshapiro.test(t(div_shannon$shannon))


kruskal.test(simpson ~ Group, data = M0.simpson)
kruskal.test(simpson ~ Group, data = M6.simpson)
kruskal.test(simpson ~ Group, data = M12.simpson)


cbind(div_simpson, sam127) %>% mutate(value = simpson) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>% group_by(interaction(Group,Month)) %>%
  summarise(mean(simpson))

pdf("Simpson_diversity_relativeAbundance.pdf", width=8, height=5)
cbind(div_simpson, sam127) %>% mutate(value = simpson) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>%
  ggplot(aes(x=Month, y=value,  group=interaction(Month, Group), color=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.2, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), alpha=0.2) +
  #showSignificance(c(0.7,1.3),0.99,-0.5,"NS")+
  #showSignificance(c(1.7,2.3),0.99,-0.5,"NS")+
  #showSignificance(c(2.7,3.3),0.99,-0.5,"NS")+
  clean_background +
  scale_fill_nejm() +
  scale_color_nejm() +
  labs(x = "Month",
       y = "Simpson index",
       title = "Simpson diversity at three time points") +
  ylim(0.9,1) +
#  scale_fill_discrete(name = "Groups", labels  = c("1: LC-PUFA \n (n = 129)", "2: HND \n (n = 138)", "3: NNR-recommend \n (n = 114)")) +
  guides(color=FALSE, alpha = FALSE)+
  
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16)
  )
dev.off()

############ shannon ################

div_shannon <- diversity((relative_abundance_matrix), index = "shannon")
##div_shannon <- diversity(raw_abundance_matrix, index="shannon")

div_shannon %>% filter(row.names(div_shannon) %in% rownames(sam127)) -> div_shannon

cbind(div_shannon, sam127) %>% mutate(Month=factor(Month, levels=c("M0","M6","M12")))%>%
  group_by(Month, Group) %>% summarize(mean=mean(shannon), sd=sd(shannon)) 

cbind(div_shannon, sam127) %>% mutate(Month=factor(Month, levels=c("M0","M6","M12")))%>%
  group_by(Month) %>% summarize(mean=mean(shannon), sd=sd(shannon)) 

cbind(div_shannon, sam127) %>% mutate(Month=factor(Month, levels=c("M0","M6","M12")))%>%
  dplyr::select(shannon,Month, Group) %>%
  group_by(Month, Group) %>% summarize(mean=mean(shannon))



cbind(div_shannon, sam127) %>% filter(Month=="M0") -> M0.shannon
cbind(div_shannon, sam127) %>% filter(Month=="M6") -> M6.shannon
cbind(div_shannon, sam127) %>% filter(Month=="M12") -> M12.shannon


kruskal.test(shannon ~ Group, data=M0.shannon)
kruskal.test(shannon ~ Group, data=M6.shannon)
kruskal.test(shannon ~ Group, data=M12.shannon)

### 
car::Anova(aov(shannon ~ Group, data = M0.shannon))
car::Anova(aov(shannon ~ Group, data = M6.shannon))
car::Anova(aov(shannon ~ Group, data = M12.shannon))

summary(aov(shannon ~ Group, data = M0.shannon))
#The p-value is 0.327 that is greater than 0.05, 
#so the covariate shannon and the group 
#are independent to each other.

######## Adjust for sex, age, BMI, and baseline

M6.shannon %>% cbind(M0=M0.shannon$shannon)  -> M60.shannon

M60.shannon %>% str()


ancova_model <- aov(shannon ~ Group + age + BMI + gender + Diabetes + Metformin + M0,
                    data=M60.shannon)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

M12.shannon %>% cbind(M0.shannon$shannon)  -> M012.shannon
M012.shannon %>% rename( `M0.shannon$shannon` = "M0") -> M012.shannon

ancova_model <- aov(shannon ~ Group + age + BMI + gender + Diabetes + Metformin + M0,
                    data=M012.shannon)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)


cbind(div_shannon,sam127) %>% mutate(value = shannon) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>%
  group_by(Group) %>% summary()

cbind(div_shannon, sam127) %>% mutate(value = shannon) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>% ggplot(aes(x=Month, y=value,  group=interaction(Month, Group), color=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.2, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), alpha=0.2) +
  clean_background +
  scale_fill_nejm() +
  scale_color_nejm() +
  clean_background +
  labs(x = "Month",
       y = "Shannon index",
       title = "Shannon diversity at three time points") +
  #scale_fill_discrete(name = "Groups", labels  = c("1: LC-PUFA \n (n = 129)", "2: HND \n (n = 138)", "3: NNR-recommend \n (n = 114)")) +
  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16)
  )

pdf("Shannon_diversity_relativeAbundance.pdf", width=8, height=5)
cbind(div_shannon, sam127) %>% mutate(value = shannon) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>%
  ggplot(aes(x=Month, y=value,  group=interaction(Month, Group), alpha=0.1)) +
  geom_boxplot(aes(fill=Group, alpha=0.1))+
  geom_jitter() +
  stat_summary(fun.y=mean, geom="line", position=position_dodge(width=0.8), size=1,
               aes(group=Group, color=Group))+
  stat_summary(fun=mean, geom="point", position=position_dodge(width=0.8), size=2.5, 
               aes(group=Group, colour=Group))+
  scale_color_manual(values = c("red","green","blue")) + 
 # showSignificance(c(0.7,1.3),5.1,-0.05,"NS")+
#  showSignificance(c(1.7,2.3),5.1,-0.05,"NS")+
#  showSignificance(c(2.7,3.3),5.1,-0.05,"NS")+
  clean_background +
  labs(x = "Month",
       y = "MGS mean Shannon diversity",
       title = "Shannon diversity at three time points") +
  scale_fill_discrete(name = "Groups", labels  = c("1: LC-PUFA \n (n = 129)", "2: HND \n (n = 138)", "3: NNR-recommend \n (n = 114)")) +
  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
  )
dev.off()

pdf("Shannon_diversity_rawAbundance.pdf", width=8, height=5)
cbind(div_shannon, sam127) %>% mutate(value = shannon) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>%
  ggplot(aes(x=Month, y=value,  group=interaction(Month, Group), alpha=0.1)) +
  geom_boxplot(aes(fill=Group, alpha=0.1))+
  geom_jitter() +
  stat_summary(fun.y=mean, geom="line", position=position_dodge(width=0.8), size=1,
               aes(group=Group, color=Group))+
  stat_summary(fun=mean, geom="point", position=position_dodge(width=0.8), size=2.5, 
               aes(group=Group, colour=Group))+
  scale_color_manual(values = c("red","green","blue")) + 
  showSignificance(c(0.7,1.3),5.1,-0.05,"NS")+
  showSignificance(c(1.7,2.3),5.1,-0.05,"*")+
  showSignificance(c(2.7,3.3),5.1,-0.05,"NS")+
  clean_background +
  labs(x = "Months",
       y = "MGS mean Shannon diversity",
       title = "Shannon diversity at three time points") +
  scale_fill_discrete(name = "Groups", labels  = c("1: LC-PUFA \n (n = 129)", "2: HND \n (n = 138)", "3: NNR-recommend \n (n = 114)")) +
  guides(color=FALSE, alpha = FALSE)
dev.off()

div_shannon %>% cbind(sam127) %>% mutate(value = shannon)


div_shannon %>% cbind(sam127) %>% mutate(value = shannon) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>%
  dplyr::select(value, ID, Month, Group, age, gender, BMI, Diabetes, Metformin) %>% pivot_wider(names_from = Month, values_from = value) %>%
  mutate(delta1=M6-M0, 
         delta2=M12-M6, 
         delta3 = M12 -M0,
         delta_M06 = (M6-M0)/M0, 
         delta_M126 = (M12-M6)/M6,
         delta_M012 = (M12-M0)/M0) -> div_shannon_dis

div_shannon_dis %>% group_by(Group) %>% summary()

div_shannon_dis

######## Adjust for sex, age, BMI, 

ancova_model <- aov(delta1 ~ Group + age + BMI + gender + Diabetes 
                    + Metformin + M0, 
                    data=div_shannon_dis)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

ancova_model <- aov(delta2 ~ Group + age + BMI + gender + Diabetes 
                    + Metformin + M0, 
                    data=div_shannon_dis)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

ancova_model <- aov(delta_M06 ~ Group + age + BMI + gender + Diabetes 
                    + Metformin + M0, 
                    data=div_shannon_dis)
Anova(ancova_model, type = "III")
postHocs <- glht(ancova_model, linfct=mcp(Group="Tukey"))
summary(postHocs)

cbind(div_shannon,sam127) %>% mutate(value = shannon) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>%
  group_by(Group) %>% summary()


div_shannon_dis %>% dplyr::select(Group, delta1) %>%
  filter(! is.na(delta1)) %>% group_by(Group) %>% summarize(mean(delta1))

div_shannon_dis %>%  filter(! is.na(delta1)) %>% filter(! is.na(delta2))%>%
  group_by(Group) %>% summarize(m1 = mean(delta1), m2=mean(delta2), m3=mean(delta3)) %>%
  mutate(fc12o60=m1/m2) 

pdf("Shannon_diversity_change_M0-M6.pdf")
div_shannon_dis %>% dplyr::select(Group, delta1) %>% 
  filter(! is.na(delta1)) %>% group_by(Group) %>% 
  summarize(mean=mean(delta1)) -> mean_delta1
mean_delta1
div_shannon_dis %>% dplyr::select(Group, delta1) %>% 
  filter(! is.na(delta1)) %>% group_by(Group) %>% 
  ggplot(aes(x=delta1, fill=Group))+
  geom_density(alpha=0.4)+
  geom_vline(data=mean_delta1, aes(xintercept=mean, color=Group), size=1.2) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="black") +
  coord_flip() +
  facet_wrap(~Group) +
  theme(panel.spacing.x = unit(-0.4, "lines"))+
  clean_background  +
  labs(x="Shannon diversity change",
       y="The density of the inviduals") +
  ggtitle(expression(~Delta[M6-M0]))+
  scale_fill_discrete(name="Group", labels=c("1: LC-PUFA", "2: HND","3: NNR-recommend")) +
  guides(linetype="none", color="none") +
  theme(legend.position="bottom",
        axis.text = element_text(face="bold")) +
  xlim(-1.2, 1)

dev.off()

div_shannon_dis
pdf("Shannon_diversity_change_M6-M12.pdf")
div_shannon_dis %>% dplyr::select(Group, delta2) %>% 
  filter(! is.na(delta2)) %>% group_by(Group) %>% 
  summarize(mean=mean(delta2)) -> mean_delta2

div_shannon_dis %>% dplyr::select(Group, delta2) %>% 
  filter(! is.na(delta2)) %>% group_by(Group) %>% 
  ggplot(aes(x=delta2,fill=Group))+
  geom_density(alpha=0.4)+
  geom_vline(data=mean_delta2, aes(xintercept=mean, color=Group), size=1.2) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="black") +
  coord_flip() +
  facet_wrap(~Group) +
  theme(panel.spacing.x = unit(-0.4, "lines"))+
  clean_background  +
  labs(x="Shannon diversity change",
       y="The density of the inviduals") +
  ggtitle(expression(~Delta[M12-M6]))+
  scale_fill_discrete(name="Group", labels=c("1: LC-PUFA", "2: HND","3: NNR-recommend")) +
  guides(linetype="none", color="none") +
  theme(legend.position="bottom",
        axis.text = element_text(face="bold")) +
  xlim(-1.2, 1) 

dev.off()
div_shannon_dis

pdf("Shannon_diversity_change_M12-M0.pdf")
div_shannon_dis %>% dplyr::select(Group, delta3) %>% 
  filter(! is.na(delta3)) %>% group_by(Group) %>% 
  summarize(mean=mean(delta3)) -> mean_delta3

div_shannon_dis %>% dplyr::select(Group, delta3) %>% 
  filter(! is.na(delta3)) %>% group_by(Group) %>% 
  ggplot(aes(x=delta3,fill=Group))+
  geom_density(alpha=0.4)+
  geom_vline(data=mean_delta3, aes(xintercept=mean, color=Group), size=1.2) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="black") +
  coord_flip() +
  facet_wrap(~Group) +
  theme(panel.spacing.x = unit(-0.4, "lines"))+
  clean_background  +
  labs(x="Shannon diversity change",
       y="The density of the inviduals") +
  ggtitle(expression(~Delta[M12-M0]))+
  scale_fill_discrete(name="Group", labels=c("1: LC-PUFA", "2: HND","3: NNR-recommend")) +
  guides(linetype="none", color="none") +
  theme(legend.position="bottom",
        axis.text = element_text(face="bold")) +
  xlim(-1.2, 1) 

dev.off()


pdf("Shannon_diversity_change_M0-M6divM0.pdf")
div_shannon_dis %>% dplyr::select(Group, delta_M06) %>% 
  filter(! is.na(delta_M06)) %>% group_by(Group) %>% 
  summarize(mean=mean(delta_M06)) -> mean_deltaM06


div_shannon_dis %>% dplyr::select(Group, delta_M06) %>% 
  filter(! is.na(delta_M06)) %>% group_by(Group) %>% 
  ggplot(aes(x=delta_M06, fill=Group))+
  geom_density(alpha=0.4)+
  geom_vline(data=mean_deltaM06, aes(xintercept=mean, color=Group), size=1.2) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="black") +
  coord_flip() +
  facet_wrap(~Group) +
  theme(panel.spacing.x = unit(-0.2, "lines"))+
  clean_background  +
  labs(x="Shannon diversity change",
       y="The density of the inviduals") +
  ggtitle(expression(~Delta[M6-M0]/M0))+
  scale_fill_discrete(name="Group", labels=c("1: LC-PUFA", "2: HND","3: NNR-recommend")) +
  guides(linetype="none", color="none") +
  theme(legend.position="bottom",
        axis.text = element_text(face="bold")) 

dev.off()

div_shannon_dis %>% str()

pdf("Shannon_diversity_change_M12-M0divM0.pdf")
div_shannon_dis %>% dplyr::select(Group, delta_M012) %>% 
  filter(! is.na(delta_M012)) %>% group_by(Group) %>% 
  summarize(mean=mean(delta_M012)) -> mean_deltaM012


div_shannon_dis %>% dplyr::select(Group, delta_M012) %>% 
  filter(! is.na(delta_M012)) %>% group_by(Group) %>% 
  ggplot(aes(x=delta_M012, fill=Group))+
  geom_density(alpha=0.4)+
  geom_vline(data=mean_deltaM012, aes(xintercept=mean, color=Group), size=1.2) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="black") +
  coord_flip() +
  facet_wrap(~Group) +
  theme(panel.spacing.x = unit(-0.2, "lines"))+
  clean_background  +
  labs(x="Shannon diversity change",
       y="The density of the inviduals") +
  ggtitle(expression(~Delta[M12-M0]/M0))+
  scale_fill_discrete(name="Group", labels=c("1: LC-PUFA", "2: HND","3: NNR-recommend")) +
  guides(linetype="none", color="none") +
  theme(legend.position="bottom",
        axis.text = element_text(face="bold")) 

dev.off()

div_shannon_dis %>% str()
pdf("Shannon_diversity_change_M12-M6divM6.pdf")
div_shannon_dis %>% dplyr::select(Group, delta_M126) %>% 
  filter(! is.na(delta_M126)) %>% group_by(Group) %>% 
  summarize(mean=mean(delta_M126)) -> mean_deltaM126


div_shannon_dis %>% dplyr::select(Group, delta_M126) %>% 
  filter(! is.na(delta_M126)) %>% group_by(Group) %>% 
  ggplot(aes(x=delta_M126, fill=Group))+
  geom_density(alpha=0.4)+
  geom_vline(data=mean_deltaM126, aes(xintercept=mean, color=Group), size=1.2) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="black") +
  coord_flip() +
  facet_wrap(~Group) +
  theme(panel.spacing.x = unit(-0.2, "lines"))+
  clean_background  +
  labs(x="Shannon diversity change",
       y="The density of the inviduals") +
  ggtitle(expression(~Delta[M12-M6]/M6))+
  scale_fill_discrete(name="Group", labels=c("1: LC-PUFA", "2: HND","3: NNR-recommend")) +
  guides(linetype="none", color="none") +
  theme(legend.position="bottom",
        axis.text = element_text(face="bold")) 

dev.off()


############### Figure Summary ###############

p1 <-cbind(div_shannon, sam127) %>% mutate(value = shannon) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>% ggplot(aes(x=Month, y=value,  group=interaction(Month, Group), color=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.2, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), alpha=0.2) +
  clean_background +
  scale_fill_nejm() +
  scale_color_nejm() +
  clean_background +
  ylim(2.5,5) +
  labs(x = "Month",
       y = "Shannon index") +
  #scale_fill_discrete(name = "Groups", labels  = c("1: LC-PUFA \n (n = 129)", "2: HND \n (n = 138)", "3: NNR-recommend \n (n = 114)")) +
  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16)
  )

p2 <- cbind(div_simpson, sam127) %>% mutate(value = simpson) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>%
  ggplot(aes(x=Month, y=value,  group=interaction(Month, Group), color=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.2, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), alpha=0.2) +
  #showSignificance(c(0.7,1.3),0.99,-0.5,"NS")+
  #showSignificance(c(1.7,2.3),0.99,-0.5,"NS")+
  #showSignificance(c(2.7,3.3),0.99,-0.5,"NS")+
  clean_background +
  scale_fill_nejm() +
  scale_color_nejm() +
  labs(x = "Month",
       y = "Simpson index") +
  ylim(0.89,1) +
  #  scale_fill_discrete(name = "Groups", labels  = c("1: LC-PUFA \n (n = 129)", "2: HND \n (n = 138)", "3: NNR-recommend \n (n = 114)")) +
  guides(color=FALSE, alpha = FALSE)+
  
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16)
  )

p3 <- pd.result %>% cbind(sam127) %>% dplyr::select(PD, SR, Month, Group) %>%
  group_by(Group) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>%
  mutate(Group=as.factor(Group)) %>%
  ggplot(aes(x=Month, y = SR,  group=interaction(Month, Group), color=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.2, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), alpha=0.2) +
  scale_fill_nejm() +
  scale_color_nejm() +
  #showSignificance(c(0.7,1.3),470,-0.5,"NS")+
  #showSignificance(c(1.7,2),470,-0.5,"NS")+
  showSignificance(c(1.9,2.3),465,-0.5,"*")+
  #showSignificance(c(2.7,3.3),470,-0.5,"NS")+
  clean_background +
  labs(x = "Month",
       y = "Richness") +
  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16)
  )

p4 <- pd.result %>% cbind(sam127) %>% dplyr::select(PD, SR, Month, Group) %>%
  group_by(Group) %>% 
  mutate(Month=factor(Month, levels=c("M0","M6","M12"))) %>%
  mutate(Group=as.factor(Group)) %>%
  ggplot(aes(x=Month, y = PD,  group=interaction(Month, Group), color=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.2, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), alpha=0.2) +
  scale_fill_nejm() +
  scale_color_nejm() +
  #scale_color_manual(values = c("red","green","blue")) + 
  #showSignificance(c(0.7,1.3),1,-0.5,"NS")+
  showSignificance(c(1.9,2.3),1400,-0.5,"*")+
  #showSignificance(c(2.7,3.3),1,-0.5,"NS")+
  clean_background +
  labs(x = "Month",
       y = "Faith's Phylogenic Diversity") +
  scale_fill_discrete(name = "Groups", labels  = c("Group 1: LC-PUFA (n = 129)", "Group 2: HND (n = 138)", "Group 3: NNR-recommend (n = 114)")) +
  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16),
        legend.position="bottom"
  )

par(mfrow=c(1,2))
pdf("try.pdf", width=12, height=8)
p1
p2
p3
p4
dev.off()
