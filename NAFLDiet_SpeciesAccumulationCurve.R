source("NAFLDiet_sample.R")

############# species accumulation curve ######
############# species accumulation curve ######
#############                             ######
relative_abundance_matrix %>% t() %>% as.data.frame.matrix() -> re_ab_matrix.df

cbind(re_ab_matrix.df[rownames(sam127),], sam127)-> abundance_with_group 


abundance_with_group  %>% dplyr::filter(Month == "M0") %>% 
  dplyr::select(c(1:6809)) %>% specaccum() -> sac_M0
abundance_with_group  %>% dplyr::filter(Month=="M6") %>% 
  dplyr::select(c(1:6809)) %>% specaccum() -> sac_M6
abundance_with_group  %>% dplyr::filter(Month=="M12") %>% 
  dplyr::select(c(1:6809)) %>% specaccum() -> sac_M12

(max(sac_M6$richness)-max(sac_M0$richness))/max(sac_M0$richness)*100
(max(sac_M12$richness)-max(sac_M6$richness))/max(sac_M6$richness)*100
(max(sac_M12$richness)-max(sac_M0$richness))/max(sac_M0$richness)*100


sac_df0 <- data.frame(Site = sac_M0$sites,
                      Richness = sac_M0$richness,
                      sd=sac_M0$sd,
                      col1="red")
sac_df6 <- data.frame(Site = sac_M6$sites, 
                      Richness = sac_M6$richness,
                      sd=sac_M6$sd,
                      col1="yellow")
sac_df12 <- data.frame(Site = sac_M12$sites, 
                       Richness = sac_M12$richness,
                       sd=sac_M12$sd,
                       col1="green")

pdf("Richness_2.pdf")
ggplot() +
  geom_line(data=sac_df0,aes(x=Site, y=Richness+1.5*sd)) +
  geom_line(data=sac_df0,aes(x=Site, y=Richness-1.5*sd)) +
  geom_line(data=sac_df6,aes(x=Site, y=Richness+1.5*sd)) +
  geom_line(data=sac_df6,aes(x=Site, y=Richness-1.5*sd)) +
  geom_line(data=sac_df12,aes(x=Site, y=Richness+1.5*sd)) +
  geom_ribbon(data=sac_df0, aes(x=Site, ymin=(Richness-1.5*sd), 
                                ymax=(Richness+1.5*sd)), fill='red', alpha=0.7)+
  geom_ribbon(data=sac_df6, aes(x=Site, ymin=(Richness-1.5*sd), 
                                ymax=(Richness+1.5*sd)), fill='green', alpha=0.7)+
  geom_ribbon(data=sac_df12, aes(x=Site, ymin=(Richness-1.5*sd), 
                                 ymax=(Richness+1.5*sd)), fill='yellow', alpha=0.7) +
  geom_line(data=sac_df0,aes(x=Site, y=Richness), color='red') +
  geom_line(data=sac_df6, aes(x=Site, y=Richness), color='green') +
  geom_line(data=sac_df12, aes(x=Site, y=Richness), color='yellow') +
  geom_line(data=sac_df12,aes(x=Site, y=Richness-1.5*sd)) +
  clean_background
dev.off()


pdf("Richness_1.pdf")
ggplot() +
  geom_line(data=sac_df0,aes(x=Site, y=Richness+1.5*sd)) +
  geom_line(data=sac_df0,aes(x=Site, y=Richness-1.5*sd)) +
  geom_line(data=sac_df6,aes(x=Site, y=Richness+1.5*sd)) +
  geom_line(data=sac_df6,aes(x=Site, y=Richness-1.5*sd)) +
  geom_line(data=sac_df12,aes(x=Site, y=Richness+1.5*sd)) +
  geom_ribbon(data=sac_df0, aes(x=Site, ymin=(Richness-1.5*sd), 
                                ymax=(Richness+1.5*sd)), fill='red', alpha=0.7)+
  geom_ribbon(data=sac_df6, aes(x=Site, ymin=(Richness-1.5*sd), 
                                ymax=(Richness+1.5*sd)), fill='green', alpha=0.7)+
  geom_ribbon(data=sac_df12, aes(x=Site, ymin=(Richness-1.5*sd), 
                                 ymax=(Richness+1.5*sd)), fill='yellow', alpha=0.7) +
  geom_line(data=sac_df0,aes(x=Site, y=Richness), color='red') +
  geom_line(data=sac_df6, aes(x=Site, y=Richness), color='green') +
  geom_line(data=sac_df12, aes(x=Site, y=Richness), color='yellow') +
  geom_line(data=sac_df12,aes(x=Site, y=Richness-1.5*sd)) +
  xlim(125,128)+
  ylim(2065,2130)+
  clean_background +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  labs(x = "Samples")
dev.off()

max(sac_df0$Richness)-max(sac_df6$Richness)
max(sac_df6$Richness)-max(sac_df12$Richness)

plot(sac_M0, ci.type="polygon",ci.col="yellow", ylim=c(0,3700))
plot(sac_M6, ci.type="polygon",ci.col="red",ylim=c(0,3700))


dev.off()
