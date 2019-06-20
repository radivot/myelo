library(shiny)

# Define UI
ui <- fluidPage(
  
  # App title ----
  titlePanel("Model of Morgan Craig et al Bull Marh Biol (2016)"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the LSC multiplier ----
      sliderInput(inputId = "mult",
                  label = "Chemo dose in mg/m2:",
                  min = .1,
                  max = 10,
                  step= 0.1,
                  value = 4) # initial value => Fig. 9
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      plotOutput(outputId = "distPlot")
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({
    
    library(tidyverse)
    library(deSolve)
    library(rodeoExt)
    
    craig16<-function(Time, State, Pars) {  # with Chemo and GCSF subQ
      with(as.list(c(State, Pars)), {
        fbeta=function(Q) fQ/(1+(Q/the2)^s2)
        fkapN=function(G1) kapss + (kapss-kapmin)*(G1^s1-G1ss^s1)/(G1^s1+G1ss^s1)
        fetaNP=function(G1) etaNPss  + (etaNPss -etaNPmin)*(bNP/G1ss)*(G1-G1ss)/(G1+bNP)
        # next func reduces to above when Cp=0 (no chemo)
        # fetaNPchemo=function(G1,Cp) etaNPInf  + (fetaNP(G1)-etaNPInf)/(1+(Cp/EC50)^Sc)
        fetaNPchemo=function(G1,Cp) etaNPInf  + (fetaNP(G1)-etaNPInf)/(1+(Cp/EM50)^Sc)
        Vn=function(G1) 1 + (Vmax-1)*(G1-G1ss)/(G1-G1ss+bV)
        GBF=G2/(V*(Nr+N))
        GBFss=G2ss/(V*(Nrss+Nss))
        phiNr=function(GBF) phiNrss  + (phiNrmax-phiNrss)*(GBF-GBFss)/(GBF-GBFss+bG)
        beta=fbeta(Q)
        kapG1=fkapN(G1)
        
        dIc=0 #Constant chemo infusion rate   
        dIg=0 # same for GCSF.  Changes only by events
        
        if (Time < 0) {
          dQ=-(beta+kapG1+kapDel)*Q + AQss*beta*Q
          dNr=-gamNr*Nr-phiNr(GBF)*Nr + (An/1e3)*kapG1*Q  # 1e3 maps units of Q to N
          dN= phiNr(GBF)*Nr-gamNss*N
          # dG1=0  # 4  unbound circulating GCSF
          # dG2=0   #5  #bound GCSF
          dG1=Gprod - kren*G1 - k12g*((Nr+N)*V-G2)*G1^Pow+k21g*G2
          dG2=k12g*((Nr+N)*V-G2)*G1^Pow -kint*G2 - k21g*G2
          dTn=0  #6
          dAn=0  #7
          dAq=0  #8
          dCp=0  #9
          dCf=0
          dCs1=0
          dCs2=0
          dGs=0 #skin pool that feeds into G1
        }	else {
          
          Qst=lagvalue(Time - tauS,1)
          betast=fbeta(Qst)
          dQ=-(beta+kapG1+kapDel)*Q + Aq*betast*Qst
          
          tauNM=Tn-tauNP  # Tn now time to reserves
          G1nmt=lagvalue(Time - tauNM,4)
          Vrat=Vn(G1)/Vn(G1nmt)
          dTn=1-Vrat
          
          G1nt=lagvalue(Time - Tn,4)
          Qnt=lagvalue(Time - Tn,1)
          kapG1nt=fkapN(G1nt)
          dNr=(An/1e3)*kapG1nt*Qnt*Vrat-gamNr*Nr-phiNr(GBF)*Nr # 1e3 maps units of S to N
          
          dN = phiNr(GBF)*Nr-gamNss*N
          
          dG1=Gprod + ka*Gs + Ig - kren*G1 - k12g*((Nr+N)*V-G2)*G1^Pow + k21g*G2
          dG2=k12g*((Nr+N)*V-G2)*G1^Pow - kint*G2 - k21g*G2
          
          Cpnt=lagvalue(Time - Tn,9)
          Cpnmt=lagvalue(Time - tauNM,9)
          
          etaNPnt=fetaNPchemo(G1nt,Cpnt)   # etaNPnt=fetaNP(G1nt)
          etaNPnmt=fetaNPchemo(G1nmt,Cpnmt)     # etaNPnmt=fetaNP(G1nmt)
          dAn=An*((1-dTn)*(etaNPnmt-etaNPnt)-gamNMss*dTn)
          
          Cpst=lagvalue(Time - tauS,9)
          dAq=Aq*hQ*(Cpst-Cp)
          
          dCp=k21*Cf+k31*Cs1-(k12+k13+kelC)*Cp  + Ic
          dCf=k12*Cp+k42*Cs2-(k21+k24)*Cf
          dCs1=k13*Cp-k31*Cs1
          dCs2=k24*Cf-k42*Cs2
          dGs=-ka*Gs + Ig # G in skin, add F*Dg/Vd to Gs at each injection
        }
        list(c(dQ,dNr,dN,dG1,dG2,dTn,dAn,dAq,dCp,dCf,dCs1,dCs2,dGs,dIc,dIg),c(ANC=N*8.19))
      })
    }
    
    craigPars16=c(gamSss = 0.1, tauS = 2.8, AQss = 1.5116, fQ = 8, s2 = 2, the2 = 0.0809, 
                  kapDel = 0.0146, kapss = 0.0073325, kapmin = 0.0052359, s1 = 1.5, 
                  etaNPss = 1.6647, bNP = 0.022868, etaNPmin = 1.406, tauNP = 7.3074, 
                  Vmax = 7.867, bV = 0.24611, aNM = 3.9, gamNMss = 0.1577, phiNrss = 0.364, 
                  phiNrmax = 4.1335, bG = 0.00018924, gamNr = 0.0063661, gamNss = 2.1875, 
                  G1ss = 0.025, GBFss = 1.5823e-05, Gprod = 0.014161, V = 0.525, 
                  kren = 0.16139, kint = 462.42, k12g = 2.2423, k21g = 184.87, 
                  Pow = 1.4608, Qss = 1.1, betaQss = 0.043, Nss = 0.3761, Ncircss = 0.22, 
                  Nrss = 2.26, Npss = 0.93, Nmss = 4.51, G2ss = 2.1899e-05, tauNr = 2.7, 
                  tauNcircss = 0.4571429, tauhalf = 7.6, ANss = 103780, bbarV = 0.031283, 
                  phiRatio = 11.356, phiMin = 0.020056, theta = 0.15096, Cko = 0.25, 
                  mu = 0.84458, Vd300 = 4754.7, F300 = 0.64466, ka300 = 8.0236, 
                  Vd375 = 2322.9, F375 = 0.49964, ka375 = 6.6133, Vd750 = 2178, 
                  F750 = 0.75, ka750 = 5.143, Vd = 2178, F = 0.75, ka = 5.143, 
                  k21 = 18.2222, k31 = 0.699, k12 = 90.2752, k13 = 8.2936, kelC = 132.0734, 
                  k42 = 62.5607, k24 = 9.2296, BSA = 1.723, hQ = 0.0079657, EC50 = 0.7539, 
                  EM50 = 24.65253, Sc = 0.89816, etaNPInf = 0)
    
    (x0=c(Q=craigPars16[["Qss"]],Nr=craigPars16[["Nrss"]],N=craigPars16[["Nss"]],
          G1=craigPars16[["G1ss"]],G2=craigPars16[["G2ss"]],
          Tn=craigPars16[["tauNP"]]+craigPars16[["aNM"]],  #here time Tn is to reserves
          An=craigPars16[["ANss"]],Aq=craigPars16[["AQss"]],Cp=0,Cf=0,Cs1=0,Cs2=0,Gs=0,Ic=0,Ig=0))
    
    (gtimes=as.numeric(t(outer(seq(0,80,14),4:13,"+"))))
    n=length(gtimes)
    (eventG=tibble(var=rep("Gs",n),
                   time=gtimes,
                   # value=rep(craigPars16[["F300"]]*300e3/craigPars16[["Vd300"]],n),
                   value=rep(craigPars16[["F750"]]*750e3/craigPars16[["Vd750"]],n),
                   method=rep("add",n)))
    
    (ctimes=seq(0,80,14))
    nc=length(ctimes)
    (delt=round(1/24,2)) # 1-hour chemo infusions
    dose=4*craigPars16[["BSA"]]*1e3# ug of chemo per injection (D=4 mg/m2)
    infusionRate=dose/delt # this feeds straight into dCp
    (eventCon=tibble(var=rep("Ic",nc),
                     time=ctimes,
                     value=rep(infusionRate,nc),
                     method=rep("rep",nc)))
    
    (eventCoff=tibble(var=rep("Ic",nc),
                      time=ctimes+delt,
                      value=rep(0,nc),
                      method=rep("rep",nc)))
    (eventdat=as.data.frame(bind_rows(eventG,eventCon,eventCoff)%>%arrange(time)))
    
    times <- seq(-15,85,by=.01)
    yout <- dede(x0,times = times, func = craig16,	parms = craigPars16,
                 events=list(data=eventdat),method="lsodar")
    
    myPlot=function(yout,cc) {
      D=data.frame(yout)
      D%>%filter(time>-.1,time<1)
      head(D)
      tail(D)
      d=D%>%select(time:Cp,Cs1,ANC)%>%gather(key="Lab",value="Value",-time)%>%
        mutate(Lab=factor(Lab,levels=c("Q","Nr","N","G1","G2","Tn","An","Aq","Cp","Cs1","ANC")))
      tc=function(sz) theme_classic(base_size=sz)
      gx=xlab("Days")
      sbb=theme(strip.background=element_blank())
      g=d%>%ggplot(aes(x=time,y=Value))+facet_grid(Lab~.,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
      # g=d%>%ggplot(aes(x=time,y=Value))+facet_wrap(Lab~.,ncol=2,scales = "free")+geom_line(size=1)+gx+tc(14)+sbb+cc
      print(g)
    }
    cc=coord_cartesian(xlim=c(-3,85))#clips high errorbars
    # At slider default 4, the following plot should be the same as the grey
    # bundle in Fig. 9, but ANC spike a bit high rises a little early (50% blasts at ~175 days)
    myPlot(yout,cc)
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
# shinyApp(ui = ui, server = server,options=list(display.mode="showcase"))
