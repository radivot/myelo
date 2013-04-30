# function definitions to be sourced into SphaseAgentSchedules.R

simulateGrowth<-function(g,init=FALSE,times=seq(0, 672, 24),events=NULL) {
	K=as.list(t(g$params)[ifelse(init,1,2),])
	outs=ode(g$X0,times=times,
			func="derivsCyc8d2",	dllname = "myelo",initfunc = "parmsCyc8d2",
			parms=as.numeric(K),events=events,rtol=1e-5,atol=1e-5)
	g$ss=outs #********** steady state outs
	nout=dim(g$ss)[1]
	X0 = g$ss[1,2:9] # first row
	Xf = g$ss[nout,2:9] # last row
	ratio=sum(Xf[-3])/sum(X0[-3])
	totals=apply(g$ss[,2:9],1,sum)
	totL=apply(g$ss[,c(2:3,5:9)],1,sum)
	tmp=g$ss[,2:9]*(1/totals)
	g$ss=cbind(g$ss,fM=tmp[,1],fG1=tmp[,2],fMat=tmp[,3],fS=apply(tmp[,4:7],1,sum),fG2=tmp[,8],tots=totals,totL=totL)
	SSE=(10-ratio)^2
	if(init) g$SSE$initial=SSE
	if(init) g$X0=c(tmp[nout,],drug=0) else g$X0=c(Xf,drug=0)
	g$SSE$final=SSE
	g
}


fopt <- function(pars,model) {
	pars=exp(pars)
	model$params[names(pars),"final"]=pars
#  model=simulateData(model) 
	model=simulateGrowth(model) 
	model$SSE$final   
}


fitModel<-function(g) {
	g=simulateGrowth(g,init=TRUE)  # computes the initial SSE 
	p0=g$params[g$params[,"opt"],"initial"]
	names(p0)<-row.names(g$params)[g$params[,"opt"]]
	p0=sapply(p0,log)
#	opt<-optim(p0,fopt,hessian=TRUE,model=g,control=list(maxit=5000) ) 
	opt<-optim(p0,fopt,method="BFGS",hessian=TRUE,model=g)       
	opar<-exp(opt$par)
	g$params[names(opar),"final"]=opar
	g=simulateGrowth(g)
#  nout=dim(g$ss)[1]
#  g$X0 = g$ss[nout,2:4] # last row
	g
}


mkg<-function(pt,sched,Teff) {
  g=NULL
  g$data=pt
  g$X0=c(rep((1-g$data["LI"])/4,3),rep(g$data["LI"]/4,4),(1-g$data["LI"])/4,0)
  names(g$X0)<-c("M","G1","Mat","Sa","Sb","Sc","Sd","G2","drug")
  paramIC=c(
      kS=4/8,
      k2=1/3,
      kM=1,
      fSurv=1,   # hold fixed for now
      k1=1/12, 
      kMat=1/13,    # optimize
      kLeaveMarrow=1/2,
      kdG1mat=1/2,
      kdSmat=1/2
  )
  g$params=data.frame(initial=paramIC,final=paramIC,opt=c(F,F,F,F,F,T,F,F,F)) 
  g=simulateGrowth(g,init=TRUE)
  g=fitModel(g)
  g$params=cbind(g$params,taus=1/g$params[,2])
#  print(g$params)
  n=dim(g$ss)[1]
  ncol=dim(g$ss)[2]
  sim0=1e11*g$ss[(n-1):n,] # sim to 672 hours (28 days) by 24h chuncks is default in simulateGrowth
  sim0[,1]=c(-24,-.1) # show values going back one day. Start of what was 10 was set to 10^12. 
  g$X0=1e11*g$X0
  tg=g
  g=simulateGrowth(g,times = seq(0, 7*4*24),events=list(data=sched))
  outs=rbind(sim0,g$ss)
  list(g=g,days=outs[,1]/24,totals=outs[,ncol])
}


mkTmCrs<-function(pt,Teff) {
  g=NULL
  g$data=pt
  g$X0=c(rep((1-g$data["LI"])/4,3),rep(g$data["LI"]/4,4),(1-g$data["LI"])/4,0)
  names(g$X0)<-c("M","G1","Mat","Sa","Sb","Sc","Sd","G2","drug")
  paramIC=c(
      kS=4/8,
      k2=1/3,
      kM=1,
      fSurv=1,   # hold fixed for now
      k1=1/12, 
      kMat=1/13,    # optimize
      kLeaveMarrow=1/2,
      kdG1mat=1/2,
      kdSmat=1/2
  )
  g$params=data.frame(initial=paramIC,final=paramIC,opt=c(F,F,F,F,F,T,F,F,F)) 
  g=simulateGrowth(g,init=TRUE)
  g=fitModel(g)
  g$params=cbind(g$params,taus=1/g$params[,2])
  n=dim(g$ss)[1]
  ncol=dim(g$ss)[2]
  sim0=g$ss[(n-1):n,] # sim to 672 hours (28 days) by 24h chuncks is default in simulateGrowth
  sim0[,1]=c(-24,-.1) # show values going back one day. Start of what was 10 was set to 10^12. 
  sim0[,c(2:9,ncol)]=1e11*sim0[,c(2:9,ncol)]
  g$X0=sim0[2,2:10]
  g$X0[9]=1
  g=simulateGrowth(g,times=0:Teff)
  g$X0[9]=0
  (sim=rbind(sim0,g$ss))
  g=simulateGrowth(g,times=Teff:96)
  (sim=rbind(sim,g$ss))
  plot(sim[,1],sim[,14],xlab="hours",ylab="S fraction",type="l",log="y")
}

