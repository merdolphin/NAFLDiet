source("NAFLDiet_sample.R")
library(latticeExtra)
library(ggridges)
library(viridis)
library(gridExtra)
library(compositions)

### vegdist(t(tmp1), binary = TRUE) is the Sørensen index of dissimilarity

relative_abundance_MGS_data[,rownames(sam127)] -> re_ab_127

re_ab_127  %>% dplyr::select(ends_with("M0")) -> re_ab_127.M0
re_ab_127  %>% dplyr::select(ends_with("M6")) -> re_ab_127.M6
re_ab_127  %>% dplyr::select(ends_with("M12")) -> re_ab_127.M12


##########################################
##############################
############################## Average beta diversity 
#############################
############################################

summary_mean_beta <- data.frame(matrix(nrow=0, ncol=4))
for (month in c("M0", "M6", "M12")){
  for (group in c(1,2,3)){
    sam127 %>% dplyr::filter(Group==group & Month==month) %>% rownames() -> tmp
    re_ab_127[,colnames(re_ab_127) %in% tmp] -> tmp1
    vegdist(t(tmp1), binary = TRUE) -> tmp2
    c(round(mean(tmp2),3), round(sd(tmp2),3), month, group) %>% rbind(summary_mean_beta) -> summary_mean_beta
  }
}


colnames(summary_mean_beta) <-c("mean", "sd", "Month","Group")
summary_mean_beta 

pdf("AverageBetaDiversityMatrix.pdf")
summary_mean_beta %>% 
  mutate(Month=factor(Month,levels=c("M0","M6","M12")),
         mean = as.numeric(mean)) %>% 
  ggplot(aes(x=Month, y=Group, color=mean, label=round(mean,2)))+
  geom_point(size=40)  +
  clean_background +
  scale_color_gradient2(midpoint=0.54, 
                        low="blue", 
                        mid="white",
                        high="red", 
                        space ="Lab" )

  
dev.off()

pdf("AverageBetaDiversityAtThreeTimePoints.pdf")
pd <- position_dodge(0.4) 
summary_mean_beta %>% mutate(mean=as.numeric(mean), 
                             sd=as.numeric(sd), 
                             Month=factor(Month, levels=c("M0","M6","M12")), 
                             Group = as.factor(Group)) %>% 
  ggplot(aes(x=Month, y=mean, fill=Group, color=Group, group=interaction(Month, Group), alpha=0.9)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.1, position=pd, color="gray") +
  geom_point(position=pd, size=3, shape=21, fill="white")+
  geom_line(aes(group=Group, color=Group), position=pd, size=1.2) +
  scale_color_manual(values = c("blue","palegreen4","red")) + 
  clean_background +
  ylim(0.2,0.9) +
  labs(x = "Month",
       y = "Mean Sørensen index of beta diversity",
       title = "Beta diversity at three time points") +
  theme(legend.position="bottom",
        axis.text = element_text(face="bold")) +
  scale_fill_discrete(name = "Groups", labels  = c("1: LC-PUFA", "2: HND", "3: NNR-recommend")) +
  guides(alpha = FALSE)


dev.off()

################ 
###############   beta diversity change ##########
##################



beta.M0 <- vegdist(t(re_ab_127.M0), binary = TRUE)
beta.M6 <- vegdist(t(re_ab_127.M6), binary = TRUE)
beta.M12 <- vegdist(t(re_ab_127.M12), binary = TRUE)


as.matrix(beta.M6)-as.matrix(beta.M0) 


i = 1;
inc <- function(x)
{
  eval.parent(substitute(x <- x + 1))
}

MATS <- list()
summary_value_beta <- data.frame(matrix(nrow=0, ncol=3))
for (month in c("M0", "M6", "M12")){
  for (group in c(1,2,3)){
    sam127 %>% dplyr::filter(Group==group & Month==month) %>% rownames() -> tmp
    re_ab_127[,colnames(re_ab_127) %in% tmp] -> tmp1
    vegdist(t(tmp1), binary = TRUE)  -> tmp2
    MATS[[i]] <- as.matrix(tmp2)
    inc(i)
    cbind(melt(as.matrix(tmp2)),rep(month, nrow(as.matrix(tmp2))), rep(group,nrow(as.matrix(tmp2)))) %>% rbind(summary_value_beta) -> summary_value_beta
  }
}
colnames(summary_value_beta) <-c("Sample1", "Sample2", "value", "Month","Group")
## M6-M0
beta_M60g1 <- MATS[[4]] -MATS[[1]]
beta_M60g2 <- MATS[[5]] -MATS[[2]]
beta_M60g3 <- MATS[[6]] -MATS[[3]]


rbind(melt(beta_M60g1) %>% mutate(Group = rep(1, nrow(melt(beta_M60g1)))),
      melt(beta_M60g2) %>% mutate(Group = rep(2, nrow(melt(beta_M60g2))))
) %>% rbind(
      melt(beta_M60g3) %>% mutate(Group = rep(3, nrow(melt(beta_M60g3))))
) -> beta_M60g123


beta_M60g123 %>% dplyr::select(Group, value)  %>% group_by(Group) %>%
  summarize(mean=mean(value)) 


pdf("betaChangeM60.pdf", width=15, height = 5)           
betaChangeM60 <- beta_M60g123 %>% dplyr::select(Group, value) %>%  
  ggplot(aes(x=value, y=as.factor(Group), fill = 0.5 - abs(0.5-after_stat(ecdf)), 
           group=Group)) +
  stat_density_ridges(geom = "density_ridges_gradient", 
                      quantile_lines = TRUE,
                      quantiles = 2,
                      calc_ecdf = TRUE) +
  xlim(-0.1,0.1) +
  scale_fill_viridis(option = "H")+
  geom_vline(aes(xintercept=0), linetype="dashed", color="black") +
  clean_background + 
  labs(x="Beta diversity change",
       y="The density of the inviduals") +  
  ggtitle(expression(~Delta[M6-M0]))+
  guides(linetype="none", color="none") +
  theme(legend.position="bottom",
        axis.text = element_text(face="bold")) 
plot(betaChangeM60)
dev.off()

#### M12-M6

beta_M126g1 <- MATS[[7]] -MATS[[4]]
beta_M126g2 <- MATS[[8]] -MATS[[5]]
beta_M126g3 <- MATS[[9]] -MATS[[6]]


rbind(melt(beta_M126g1) %>% mutate(Group = rep(1, nrow(melt(beta_M126g1)))),
      melt(beta_M126g2) %>% mutate(Group = rep(2, nrow(melt(beta_M126g2))))
) %>% rbind(
  melt(beta_M126g3) %>% mutate(Group = rep(3, nrow(melt(beta_M126g3))))
) -> beta_M126g123

beta_M126g123 %>% dplyr::select(Group, value)  %>% group_by(Group) %>%
  summarize(mean=mean(value)) 

beta_M126g123

pdf("betaChangeM126.pdf", width=15, height = 5)           
betaChangeM126 <- beta_M126g123 %>% dplyr::select(Group, value) %>%  
  ggplot(aes(x=value, y=as.factor(Group), fill = 0.5 - abs(0.5-after_stat(ecdf)), 
             group=Group)) +
  stat_density_ridges(geom = "density_ridges_gradient", 
                      quantile_lines = TRUE,
                      quantiles = 2,
                      calc_ecdf = TRUE) +
  xlim(-0.1,0.1) +
  scale_fill_viridis(option = "H")+
  geom_vline(aes(xintercept=0), linetype="dashed", color="black") +
  clean_background + 
  labs(x="Beta diversity change",
       y="The density of the inviduals") +  
  ggtitle(expression(~Delta[M12-M6]))+
  guides(linetype="none", color="none") +
  theme(legend.position="bottom",
        axis.text = element_text(face="bold")) 
plot(betaChangeM126)
dev.off()

##### M12 - M0
beta_M120g1 <- MATS[[7]] -MATS[[1]]
beta_M120g2 <- MATS[[8]] -MATS[[2]]
beta_M120g3 <- MATS[[9]] -MATS[[3]]


rbind(melt(beta_M120g1) %>% mutate(Group = rep(1, nrow(melt(beta_M120g1)))),
      melt(beta_M120g2) %>% mutate(Group = rep(2, nrow(melt(beta_M120g2))))
) %>% rbind(
  melt(beta_M120g3) %>% mutate(Group = rep(3, nrow(melt(beta_M120g3))))
) -> beta_M120g123

beta_M120g123 %>% dplyr::select(Group, value)  %>% group_by(Group) %>%
  summarize(mean=mean(value)) 

beta_M120g123


pdf("betaChangeM120.pdf", width=15, height = 5)           
betaChangeM120 <- beta_M120g123 %>% dplyr::select(Group, value) %>%  
  ggplot(aes(x=value, y=as.factor(Group), fill = 0.5 - abs(0.5-after_stat(ecdf)), 
             group=Group)) +
  stat_density_ridges(geom = "density_ridges_gradient", 
                      quantile_lines = TRUE,
                      quantiles = 2,
                      calc_ecdf = TRUE) +
  xlim(-0.1,0.1) +
  scale_fill_viridis(option = "H")+
  geom_vline(aes(xintercept=0), linetype="dashed", color="black") +
  clean_background + 
  labs(x="Beta diversity change",
       y="The density of the inviduals") +  
  ggtitle(expression(~Delta[M12-M0]))+
  guides(linetype="none", color="none") +
  theme(legend.position="bottom",
        axis.text = element_text(face="bold")) 
plot(betaChangeM120)
dev.off()
pdf("try.pdf")
ggarrange(betaChangeM60, betaChangeM126, betaChangeM120,
        ncol=1)
dev.off()


################################################
########################
######################## Beta diversity distribution
########################
#########################################


coul <- colorRampPalette(brewer.pal(9, "PiYG"))(25)
tmp2 %>% as.matrix()  
tmp2 %>% as.matrix() %>% heatmap.2(scale = "none", col = bluered(100), 
                                   trace = "none", density.info = "none")


sam127 %>% dplyr::filter(Group==1 & Month=="M12") %>% rownames() 



summary_value_beta %>% dplyr::filter(Month=="M0")-> summary_value_beta.M0

summary_value_beta.M0 %>% head()

pdf("try.pdf")
levelplot(value ~ Sample1 * Sample2, summary_value_beta.M0,
          panel = panel.levelplot.points, cex = 1.2) + 
  layer_(panel.2dsmoother(..., n = 20))

dev.off()

summary_value_beta %>% head()

summary_value_beta %>% dplyr::filter(value!=0 & value !=1) %>%
  dplyr::select(value) %>% 
  dudi.pca(scale=TRUE, scannf=FALSE,
           center = TRUE,
           nf=5) -> beta.pca


pdf("try.pdf")
fviz_pca_ind(beta.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping)
)

dev.off()


dplyr::group_by(Group, Month) %>% mutate(value=as.numeric(value),
                                         Month=factor(Month,levels=c("M0","M6","M12")),
                                         Group=as.factor(Group)) %>%
  pivot_wider(names_from = c(Month, Group),
              values_from = value) %>% head(n=7)

dudi.pca(nf=5,
         scannf=FALSE,) -> beta.pca

pdf("try.pdf")

  summary_value_beta %>% dplyr::filter(Month=="M0" & Group==1 & value !=0)%>% 
  ggplot(aes(x=Sample1, y=Sample2, color=value)) +
  geom_point()+
    scale_color_gradient2(midpoint=0.5, 
                         low="blue", 
                         mid="white",
                         high="red", 
                         space ="Lab" )+
  clean_background
dev.off()



pdf("beta_div_M0.pdf", width=5, height=5)
summary_value_beta %>% dplyr::filter(Month=="M0")-> summary_beta_M0
summary_beta_M0 %>%
  group_by(Group) %>% summarize(mean=mean(value)) -> mean.data

mean.data %>% dplyr::select(mean) %>% colMeans() -> all_mean

data=summary_beta_M0 %>% group_by(Group) %>%  
  summarise(value1=median(value))


summary_beta_M0 %>%
  ggplot(aes(x=value, y=Group, fill = 0.5 - abs(0.5-after_stat(ecdf)), group=Group
             )) +
   stat_density_ridges(geom = "density_ridges_gradient", 
                       calc_ecdf = TRUE) +
   xlim(0.25,0.9) +
  ylim(1,5)+
  scale_fill_viridis(option = "H")+
 geom_vline(aes(xintercept=all_mean), linetype="dashed", color="black") +
  clean_background + 
  labs(x="Beta diversity distribution",
       y="The density of the inviduals") +  
   guides(linetype="none", color="none") +
  theme(legend.position="bottom",
        axis.text = element_text(face="bold")) 


dev.off()

pdf("beta_div_M6.pdf", width=5, height=5)
summary_value_beta %>% dplyr::filter(Month=="M6")-> summary_beta_M6
summary_beta_M6 %>%
  group_by(Group) %>% summarize(mean=mean(value)) -> mean.data

mean.data %>% dplyr::select(mean) %>% colMeans() -> all_mean

data=summary_beta_M6 %>% group_by(Group) %>%  
  summarise(value1=median(value))


summary_beta_M6 %>%
  ggplot(aes(x=value, y=Group, fill = 0.5 - abs(0.5-after_stat(ecdf)), group=Group
  )) +
  stat_density_ridges(geom = "density_ridges_gradient", 
                      calc_ecdf = TRUE) +
  xlim(0.25,0.9) +
  ylim(1,5)+
  scale_fill_viridis(option = "H")+
  geom_vline(aes(xintercept=all_mean), linetype="dashed", color="black") +
  clean_background + 
  labs(x="Beta diversity distribution",
       y="The density of the inviduals") +  
  guides(linetype="none", color="none") +
  theme(legend.position="bottom",
        axis.text = element_text(face="bold")) 


dev.off()

pdf("beta_div_M12.pdf", width=5, height=5)

summary_value_beta %>% dplyr::filter(Month=="M12")-> summary_beta_M12
summary_beta_M12 %>%
  group_by(Group) %>% summarize(mean=mean(value)) -> mean.data

mean.data %>% dplyr::select(mean) %>% colMeans() -> all_mean

data=summary_beta_M12 %>% group_by(Group) %>%  
  summarise(value1=median(value))


summary_beta_M12 %>%
  ggplot(aes(x=value, y=Group, fill = 0.5 - abs(0.5-after_stat(ecdf)), group=Group
  )) +
  stat_density_ridges(geom = "density_ridges_gradient", 
                      calc_ecdf = TRUE) +
  xlim(0.25,0.9) +
  ylim(1,5)+
  scale_fill_viridis(option = "H")+
  geom_vline(aes(xintercept=all_mean), linetype="dashed", color="black") +
  clean_background + 
  labs(x="Beta diversity distribution",
       y="The density of the inviduals") +  
  guides(linetype="none", color="none") +
  theme(legend.position="bottom",
        axis.text = element_text(face="bold")) 


dev.off()


##########################################
##############################
############################## ordination
#############################
#############################################
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(microbiomeutil) # some utility tools
library(microbiomeutil)
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling
devtools::install_github("microsud/microbiomeutilities")
1
p1 <- plot_taxa_cv(ex1b, plot.type="scatter")
n


library(plotly)



raw_abundance_MGS_data[,rownames(sam127)] ->   raw_ab_127


raw_ab_127  %>% dplyr::select(ends_with("M0")) -> raw_ab_127.M0
raw_ab_127  %>% dplyr::select(ends_with("M6")) -> raw_ab_127.M6
raw_ab_127  %>% dplyr::select(ends_with("M12")) -> raw_ab_127.M12

##################
help("metaMDS")
re_ab_127_NMDS <- metaMDS(t(re_ab_127.M0), k=3)

re_ab_127_NMDS$points %>% as.data.frame() -> tmp_NMDS

cbind(tmp_NMDS,sam127) -> NMDS.2
NMDS.2
NMDS.2 %>% plot_ly(x = ~MDS1, y = ~MDS2, z = ~MDS3) %>%
  add_markers(color=~Group) 


ggplot(NMDS.2, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point()

##################
### ###############################
######################################

z <- betadiver(t(re_ab_127_M0),"z")
mod <- with(sam127, betadisper(z, Group))
mod$distances
mod

pdf("try.pdf")
plot(b_m0)
dev.off()

