## Prolonged treatment-free remission in chronic myeloid leukemia patients with previous BCR-ABL1 kinase domain mutations. Claudiani et al, Haematologica (2020). 

The following code plots data for patients 4, 9 and 10 in Figure 1 of this paper. 

```
rm(list=ls())
library(myelo)
library(tidyverse)
head(d<-claudianiPon)
tc=function(sz) theme_classic(base_size=sz)
gy=ylab("BCR-ABL Percent")
sbb=theme(strip.background=element_blank())
library(ggrepel)
d%>%ggplot(aes(x=Month,y=Prct,col=Ponatinib,label=T315I))+facet_wrap(Pt~.,ncol=1,scale="free")+geom_line(size=1)+
  geom_text_repel(nudge_y=-10,nudge_x=-3) +
  gy+tc(14)+sbb+scale_y_log10()+theme(legend.position="top")
ggsave("../docs/claudianiPonLog.png",width=6,height=6)

```

![](../../docs/claudianiPonLog.png)

Percentages shown are of T315I mutations. Before ponatinib (light blue) a different TKI was used. 
After it, no TKI was used. 


