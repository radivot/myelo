raue10<-function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
				v1=kon*Epo*EpoR;
				v2=kon*Kd*EpoEpoR;
				v3=kt*Bmax;
				v4=kt*EpoR;
				v5=ke*EpoEpoR;
				v6=kex*EpoEpoRi;
				v7=kdi*EpoEpoRi;
				v8=kde*EpoEpoRi;
				dEpo = -v1+v2+v6;
				dEpoR = -v1+v2+v3-v4+v6;
				dEpoEpoR = v1-v2-v5;
				dEpoEpoRi = v5-v6-v7-v8;
				ddEpoi = v7;
				ddEpoe = v8;
				y1 = scale*(Epo+dEpoe);
				y2 = scale*EpoEpoR;
				y3 = scale*(EpoEpoRi+dEpoi);
				return(list(c(dEpo,dEpoR,dEpoEpoR,dEpoEpoRi,ddEpoi,ddEpoe),
							c(y1=y1,y2=y2,y3=y3)))
			})
}
