library(tidyverse)  
library(deSolve)
library(myelo)
(x0=craigIC[c(1,8)])
(parsQ=craigPars[c("Qss","Aqss","tauS","fQ","the2","s2")])
parsQ["kapDel"]=craigPars["kapss"]+craigPars["kapDel"]
parsQ
attach(as.list(parsQ))
Q=seq(0.01,2,0.01)
fbeta=function(Q) fQ/(1+(Q/the2)^s2)
beta=fbeta(Q)
betaSS=fbeta(Qss)
tc=function(sz) theme_classic(base_size=sz)
gx=xlab(quote(paste("Q in ",10^6," HSC/(kg body mass)")))
gy=ylab("Beta (fraction leaving Q per day)")
tibble(Q,beta)%>%ggplot(aes(x=Q,y=beta))+geom_line(size=1)+
  geom_vline(aes(xintercept=Qss))+geom_hline(aes(yintercept=betaSS))+gx+gy+tc(14)
ggsave("~/Results/myelo/Qprof.pdf",width=4, height=4)

thresh=0.1
fbetaDB=function(Q,thresh,betaSS,fQ) {
  ifelse(Q>thresh,betaSS,fQ-((fQ-betaSS)/thresh)*Q)
}

# (betaDB=fbetaDB(Q,thresh,betaSS,fQ))
# tibble(Q,betaDB)%>%ggplot(aes(x=Q,y=betaDB))+geom_line(size=1)+
#   geom_vline(aes(xintercept=Qss))+geom_hline(aes(yintercept=betaSS))+gx+gy+tc(14)
# ggsave("~/Results/myelo/QdbProf.pdf",width=5, height=5)

(betaDB=fbetaDB(Q,thresh,betaSS,fQ=2*betaSS))
tibble(Q,betaDB)%>%ggplot(aes(x=Q,y=betaDB))+geom_line(size=1)+
  geom_vline(aes(xintercept=Qss))+gx+gy+tc(14)+ylim(0,0.1)
ggsave("~/Results/myelo/QdbProf.pdf",width=4, height=4)

detach(as.list(parsQ))
