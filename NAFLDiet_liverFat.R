source("NAFLDiet_sample.R")

LiverFatFra <- read_xlsx("NAFLDIET_LiverfatFraction.xlsx")

LiverFatFra %>% separate_wider_delim(patient_id, delim = "_",
                                     names=c("prefix","ID")) %>%
  mutate(ID = as.numeric(ID)) -> LiverFatFra.1



ID127.new %>% filter(! Studie.ID %in% LiverFatFra.1$ID) %>% 
  dplyr::select(nID)

sam127 %>% dplyr::select(ID, Group) %>% distinct() -> sam127.tmp
LiverFatFra.1 %>% filter(ID %in% sam127$ID) -> LiverFatFra.2

sam127 %>% dplyr::filter(Month=="M0") %>% mutate(ID=as.numeric(ID)) %>%
  dplyr::left_join( LiverFatFra.2) -> LiverFatFra.3


sam127.tmp %>% filter(ID %in% LiverFatFra.2$ID) %>% mutate(ID=as.numeric(ID)) -> sam127.tmp

pdf("try.pdf")
dplyr::left_join(sam127.tmp, LiverFatFra.2) %>% 
  mutate(change=fatfraction_v2-fatfraction_v1) %>%
  arrange(desc(change)) %>% group_by(Group) %>% na.omit() %>% 
  ggplot(aes(x=1, y=fatfraction_v1,  fill=Group, alpha=0.1))+
  geom_boxplot(outlier.shape = NA) +
  geom_boxplot(aes(x=2, y=fatfraction_v2),outlier.shape = NA)+
  showSignificance(c(0.7,1.3),34,-0.05,"NS")+
  showSignificance(c(1.7,2),33,-0.05,"NS")+
  showSignificance(c(2,2.3),31,-0.05,"***")+
  showSignificance(c(1.7,2.3),35,-0.05,"*")+
  scale_color_manual(values = c("red","green","blue")) + 
  clean_background +
  labs(x="Month",
       y="fatfraction") +
#  scale_x_discrete(breaks=c("1","2"),labels=c("M0", "M12")) +
  guides(alpha="none") 

dev.off()

dplyr::left_join(sam127.tmp, LiverFatFra.2) %>% 
  mutate(change=fatfraction_v2-fatfraction_v1,
         percentage = change/fatfraction_v1*100) %>%
  arrange(desc(change)) %>% group_by(Group) %>% na.omit() %>% 
  summarise_all(list(mean,sd)) %>%
  dplyr::select(-c(ID_fn1, prefix_fn1, date_v1_fn1, 
                   date_v2_fn1,prefix_fn2, date_v1_fn2, date_v2_fn2 )) %>%
  as.data.frame()



dplyr::left_join(sam127.tmp, LiverFatFra.2) 
  
sam127 %>% dplyr::filter(Month=="M0") %>% mutate(ID=as.numeric(ID)) %>%
  dplyr::left_join( LiverFatFra.2) -> LiverFatFra.3

cor.test(LiverFatFra.3$fatfraction_v1, LiverFatFra.3$fatfraction_v2, method="spearman")


shapiro.test(LiverFatFra.4$fatfraction_v1)
shapiro.test(LiverFatFra.4$fatfraction_v2)

dplyr::left_join(sam127.tmp, LiverFatFra.2) %>%
  mutate(change=fatfraction_v2-fatfraction_v1) -> LiverFatFra.4


kruskal.test(fatfraction_v2 ~ Group, data=LiverFatFra.4)
kruskal.test(fatfraction_v1 ~ Group, data=LiverFatFra.4)
kruskal.test(change ~ Group, data=LiverFatFra.4)

pairwise.wilcox.test(LiverFatFra.4$fatfraction_v2, LiverFatFra.4$Group, p.adjust.method = "bonferroni")

pairwise.wilcox.test(LiverFatFra.4$change, LiverFatFra.4$Group, p.adjust.method = "bonferroni")

LiverFatFra.4 %>% dplyr::filter(Group == 1) -> LiverFatFra.4.G1
LiverFatFra.4 %>% dplyr::filter(Group == 2) -> LiverFatFra.4.G2
LiverFatFra.4 %>% dplyr::filter(Group == 3) -> LiverFatFra.4.G3

wilcox.test(LiverFatFra.4.G1$fatfraction_v1, LiverFatFra.4.G1$fatfraction_v2)
wilcox.test(LiverFatFra.4.G2$fatfraction_v1, LiverFatFra.4.G2$fatfraction_v2)
wilcox.test(LiverFatFra.4.G3$fatfraction_v1, LiverFatFra.4.G3$fatfraction_v2)



pdf("try.pdf")
dplyr::left_join(sam127.tmp, LiverFatFra.2) %>% mutate(change=fatfraction_v2-fatfraction_v1) %>%
  arrange(desc(change)) %>% filter(Group == 1) %>% 
  ggplot(aes(x=1:41, y=fatfraction_v1, fill="white"))+
  geom_col(aes(fill="blue", alpha=0.5)) +
  geom_col(aes(x=1:41, y=-fatfraction_v2, fill="red", alpha=0.5))+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray", size=1) +
  geom_smooth(aes(x=1:41, y=change),
              color="orange",fill="orange", span=1,
              level=0.999)+
  clean_background +
  xlim(0,45)+
  coord_flip()
dev.off()
  ggplot(aes(x=ID, y=change, color=Group))+
  geom_point()+
  geom_line()



vars = c("triglycerides", "BMI", "liver.fat", 
         "ASAT", "ALAT", "Glutamyltransferase",
         "Cholesterol")

fat_results <- data.frame(matrix(nrow=0, ncol=5))
fat_results
colnames(fat_results) <- c("Group","Month","mean","sd","var")

for(var in vars){
  sample.127.full %>% dplyr::select(c(contains(var), diet.group)) %>% group_by(diet.group) 
    
}

re_ab_127



for (var in vars) {
  
  
  sample.127.full %>% dplyr::select(c(contains(var), diet.group)) %>% group_by(diet.group) %>% 
    summarise_all(.funs = list(mean, sd)) %>% setNames(c("Group","M0mean","M6mean","M12mean",
                                                         "M0sd","M6sd","M12sd")) 
  
  sample.127.full %>% dplyr::select(c(contains(var), diet.group)) %>% group_by(diet.group) %>% 
    na.omit() %>%
    dplyr::summarise_all(mean) %>% pivot_longer(-diet.group) %>% setNames(c("Group","Month","mean")) %>%
    dplyr::left_join(
      sample.127.full %>% dplyr::select(c(contains(var), diet.group)) %>% group_by(diet.group) %>% 
        na.omit() %>%
        summarise_all(sd) %>% pivot_longer(-diet.group)  %>% setNames(c("Group","Month","sd"))
    )  %>% rbind(fat_results) -> fat_results
}

sample.127.full %>% dplyr::select(c(contains("ALAT"), diet.group))



fat_results %>% separate_wider_delim(Month, delim = ".", names_sep = ".", names = c("1","2"), too_many=c("merge")) %>%
  mutate(Group = as.factor(Group),
         Month.1 = factor(Month.1),
         Month.2 = factor(Month.2, 
                          levels=c("ALAT.(ukat/L)","ASAT.(ukat/L)", "Glutamyltransferase.(ukat/L)",
                                   "BMI.(kg/m2)","cholesterol.(mg)",
                                   "Cholesterol.(mmol/L)" , "HDL.Cholesterol.(mmol/L)",
                                   "LDL.cholesterol.(mmol/L)" ,   
                                   "Triglycerides.(mmol/L)"     ))
  ) %>% mutate(Month.2 = gsub("\\.", " ", Month.2)) -> fat_results.df



fat_results.df$Month.2 %>% as.factor() %>% levels() 
pd <- position_dodge(0.1) 


pdf("fatChange_v1.pdf")
fat_results.df %>% group_by(Group) %>%
  ggplot(aes(x=Month.1, y=mean, group=Group, color=Group)) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.05, position=pd, color="gray") +
  geom_point(position=pd, shape=21, fill="white", size=3) +
  geom_line(position=pd, size=1.5) +
  facet_wrap(~Month.2, scales = "free_y")+
  clean_background+
  theme(strip.background=element_rect(colour="white",
                                      fill="white"),
        strip.text = element_text(size = 12, color = "black"))+
  scale_x_discrete(labels=c("A2" = "M0", "B" = "M6",
                            "C2" = "M12"))
dev.off()
pd <- position_dodge(0.4) 
BMI %>% group_by(Group) %>% mutate(Group = as.factor(Group)) %>% 
  ggplot(aes(x=Month, y=mean, group=Group, color=Group)) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.1, position=pd, color="gray") +
  geom_point() +
  geom_line() +
  clean_background

