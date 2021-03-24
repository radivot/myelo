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

## Interpretation

Patient 4 (and perhaps 9) had CML (including a small T315I subclone) under immunological control 
that was lost by immune system exhaustion by an infection, similar to the Hong Kong 
flu of 1968 releasing latent A-bomb-induced CML clones in Hiroshima females in 1969-1974, 
see Radiat Environ Biophys. 2021;60(1):41-47. 
In the presence of TKI, loss of immuno-control allowed the pre-existing T315I clone to grow relative to WT.  Ponatinib
then drove the load back down to a size that could be controlled immunologically.  

Patient 10 had CML only under TKI control, which was lost with T315I, so the climb up was 100% T315I.  

