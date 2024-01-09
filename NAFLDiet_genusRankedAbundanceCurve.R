source("NAFLDiet_sample.R")
################################
########                #######
####### Rank of genera  ######
######                  ######
################################

cbind(raw_tax_matrix_127,relative_abundance_MGS_data_127) %>%
  dplyr::select(c(genus,matches("lunliv"))) %>%
  group_by(genus) %>% summarise(across(everything(),sum))-> 
  re_ab_127_genus

re_ab_127_genus %>% column_to_rownames(var="genus") -> re_ab_127_genus.df
### use genus raw counts ###
genus_re_ab <- read_xlsx("all-tables_2023.xlsx",
                         sheet = "genus-downsized-counts")

genus_re_ab%>% dplyr::select(c(2,rownames(sam127))) %>%
  column_to_rownames(var="genus") -> genus_re_ab.127

genus_re_ab.127 %>%  as.data.frame() %>% 
  dplyr::select( ends_with("M0") ) %>% t() %>% as.data.frame() -> 
  genus_re_ab.127.M0

genus_re_ab.127 %>%  as.data.frame() %>% 
  dplyr::select( ends_with("M6") ) %>% t() %>% as.data.frame() -> 
  genus_re_ab.127.M6

genus_re_ab.127 %>%  as.data.frame() %>% 
  dplyr::select( ends_with("M12") ) %>% t() %>% as.data.frame() -> 
  genus_re_ab.127.M12

dim(genus_re_ab.127.M0)

### M0 
rad.lognormal(genus_re_ab.127.M0[1,])$y %>% enframe() %>% arrange(desc(value)) %>%
  mutate(no = row_number(),value1=as.numeric(value))   %>%
  dplyr::select(no, value1) %>%
  setNames(replace(names(.), names(.)=="value",paste0(names(.),1))) -> 
  rad127.M0
for (i in 2:127){
  mod <- rad.lognormal(genus_re_ab.127.M0[i,])
  valuei <- paste0("value",i)
  mod$y %>% enframe() %>% arrange(desc(value)) %>%
    mutate(no = row_number(), value=as.numeric(value)) %>%
    dplyr::select(no, value) -> tmp0
  colnames(tmp0) <-c("no",valuei)
  merge(tmp0, rad127.M0, by = "no", all=TRUE) -> rad127.M0
}

## M6
rad.lognormal(genus_re_ab.127.M6[1,])$y %>% enframe() %>% arrange(desc(value)) %>%
  mutate(no = row_number(),value1=as.numeric(value))   %>%
  dplyr::select(no, value1) %>%
  setNames(replace(names(.), names(.)=="value",paste0(names(.),1))) -> 
  rad127.M6

for (i in 2:127){
  mod <- rad.lognormal(genus_re_ab.127.M6[i,])
  valuei <- paste0("value",i)
  mod$y %>% enframe() %>% arrange(desc(value)) %>%
    mutate(no = row_number(), value=as.numeric(value)) %>%
    dplyr::select(no, value) -> tmp0
  colnames(tmp0) <-c("no",valuei)
  merge(tmp0, rad127.M6, by = "no", all=TRUE) -> rad127.M6
}

## M12
rad.lognormal(genus_re_ab.127.M12[1,])$y %>% enframe() %>% arrange(desc(value)) %>%
  mutate(no = row_number(),value1=as.numeric(value))   %>%
  dplyr::select(no, value1) %>%
  setNames(replace(names(.), names(.)=="value",paste0(names(.),1))) -> 
  rad127.M12

for (i in 2:127){
  mod <- rad.lognormal(genus_re_ab.127.M12[i,])
  valuei <- paste0("value",i)
  mod$y %>% enframe() %>% arrange(desc(value)) %>%
    mutate(no = row_number(), value=as.numeric(value)) %>%
    dplyr::select(no, value) -> tmp0
  colnames(tmp0) <-c("no",valuei)
  merge(tmp0, rad127.M12, by = "no", all=TRUE) -> rad127.M12
}


ltseed=c("blank", "solid", "dashed", "dotted", "dotdash", 
         "longdash", "twodash", "1F", "F1", "4C88C488", 
         "12345678")
shape127 <- rep( c(4,8,15,16,17,18,21,22,3,42), length.out=127)

lt127 <- rep(ltseed, length.out=127)


#col127 <- rep(RColorBrewer::brewer.pal(6, "Paste1"), length.out=127)
col127 <- rep(rainbow(12), length.out=127)


clean_background <- theme(plot.background = element_rect("white"),
                          panel.background = element_rect("white"),
                          panel.grid = element_line("white"),
                          axis.line = element_line("black"),
                          axis.text = element_text(size = 12),
                          axis.title = element_text(color = "gray25"),
                          legend.text = element_text(size = 12),
                          legend.key = element_rect("white"),
                          panel.border = element_rect(colour = "black", fill=NA, size=0.8))



months <- c("M0", "M6", "M12")
rad.tmp <- data.frame()

for (i in 1:length(months)){
  
  rad127 = eval(as.symbol(paste0("rad127.",months[i])))
  
  rad127[is.na(rad127)] <- 0
  
  apply(rad127[,-1], 2, function(x) 100*x/sum(x)) -> rad127.per
  
  cbind(rad127[,1], rad127.per) %>% as.data.frame() -> rad127.1
  
  rad127.1 %>%  
    pivot_longer(cols = c(-1)) %>%
    mutate(name = as.factor(name))%>%group_by(name) -> rad127.2
  
  cbind(rad127.2, monthNo = rep(i, nrow(rad127.2))) %>% 
    rbind(rad.tmp) -> rad.tmp
  
  plot <-  rad127.1 %>%  
    pivot_longer(cols = c(-1)) %>%
    mutate(name = as.factor(name))%>%group_by(name) %>%
    ggplot(aes(x=V1, y=value, color=name))+
    geom_point(aes(shape=name)) +
    geom_line(aes(linetype=name)) +
    scale_color_manual(values=col127) +
    scale_shape_manual(values =shape127)+
    scale_linetype_manual(values=lt127)+
    scale_x_continuous(trans="log2", 
                       breaks = c(1,2,5,10,20,50,100,300))+
    scale_y_continuous(limits = c(0,55),breaks=seq(0,55, by=10))+
    clean_background+
    labs(x="Rank at the geneus level", y="Ranked abundance (%)") +
    theme(legend.position = "none",
          axis.text = element_text(face="bold")) 
  
}
rad.tmp %>% head()

  rad.tmp %>% ggplot(aes(x=V1, y=value, color=name))+
  geom_point(aes(shape=name)) +
  geom_line(aes(linetype=name)) +
  scale_color_manual(values=col127) +
  scale_shape_manual(values =shape127)+
  scale_linetype_manual(values=lt127) +
  clean_background +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.text = element_text(face="bold")) +
    scale_x_continuous(trans="log2", 
                       breaks = c(1,2,5,10,20,50,100,300))+
    scale_y_continuous(limits = c(0,55),breaks=seq(0,55, by=10))+
  clean_background+
  labs(x="Rank at the geneus level", y="Ranked abundance") +
  theme(legend.position = "none" ) +
  facet_wrap(~monthNo, ncol = 1, strip.position = "right")
  

help("facet_wrap")
plot_gg(mtplot,  raytrace = FALSE, preview = TRUE)

plot_gg(mtplot, windowsize = c(800, 800), 
        zoom = 0.85, phi = 35, theta = 30, sunangle = 225, soliddepth = -100)



pdf("try.pdf", width=10,height=5)
rad127.1 %>%  
  pivot_longer(cols = c(-1)) %>%
  mutate(name = as.factor(name))%>%group_by(name) %>%
  ggplot(aes(x=V1, y=value, color=name))+
  geom_point(aes(shape=name)) +
  geom_line(aes(linetype=name)) +
  scale_color_manual(values=col127) +
  scale_shape_manual(values =shape127)+
  scale_linetype_manual(values=lt127)+
  scale_x_continuous(trans="log2", 
                     breaks = c(1,2,5,10,20,50,100,300))+
  scale_y_continuous(limits = c(0,55),breaks=seq(0,55, by=10))+
  clean_background+
  labs(x="Rank at the geneus level", y="Ranked abundance") +
  theme(legend.position = "none" )  

dev.off()


