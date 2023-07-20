setwd('C:/Users/Administrator/Desktop')
library(ggplot2)
library(reshape)
dat <- read.csv('phylum top 10.csv')
group <- read.csv('group 1.csv')
dat$Taxonomy <- factor(dat$Taxonomy, levels = rev(dat$Taxonomy))
dat <- reshape::melt(dat, id = 'Taxonomy')
names(group)[1] <- 'variable'
dat <- merge(dat, group, by = 'variable')
color <- c("black","orange","cyan","pink","deeppink3","navy", "darkorchid1","gold", "limegreen", "dodgerblue3", "firebrick1",'#993399','#CC99FF','skyblue','#CC0000','#FF9999','#FF3300','#FFCC33','#006600','#33CC00',"#00B76D")
p<-ggplot(dat,aes(x=times,y=100*value,fill=Taxonomy))+
  geom_col(position="stack",width=18)+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black'))+
  facet_wrap(~group, scales = 'free_x',  ncol = 4) +
  theme(strip.text=element_text(size=9), strip.background=element_rect(fill="white",color="black",linetype = 1))+
  scale_fill_manual(values = color) + 
  labs(x = '', y = 'Relative Abundance(%)')
p

