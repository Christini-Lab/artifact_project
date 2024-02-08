/*
 *  SCN5A_function.cc
 *  Created by Gildas Loussouarn
 *  Copyright (c) 2019. All rights reserved.
 *
 */ 

#include "WT_HHSCN5A_global_function.h"  //header file for Na+ channel parameters


//**************************************************************************************
//	FUNCTION TO COMPUTE SCN5A NA CURRENT--HH O'Hara Rudy
//**************************************************************************************


	 double WT_INa_function (int n)		
{	

minf[n]=1/(1+exp(-(V[n]+39.57)/9.871)); 
mtau[n]=1/(6.765*exp(+(V[n]+11.64)/34.77)+8.552*exp(-(V[n]+77.42)/5.955));

jinf[n]=1/(1+exp(+(V[n]+82.9)/6.086));
jtau[n]=2.038+1/(0.02136*exp(-(V[n]+100.6)/8.281)+0.3052*exp((V[n]+0.9941)/38.45));

hinf[n]=1/(1+exp(+(V[n]+82.9)/6.086));
htau[n]=1/(1.432e-5*exp(-(V[n]+1.196)/6.285)+6.149*exp(+(V[n]+0.5096)/20.27));

if (t==0)
	{
  	m[n]=0.0;
	h[n]=1.0;
	j[n]=1.0;
	}
   
if (back == 0)
		{
		mtemp[n] = m[n];
		htemp[n] = h[n];
		jtemp[n] = j[n];
		}
	else
		{
		m[n] = mtemp[n];
		h[n] = htemp[n];
		j[n] = jtemp[n];
		}

m[n]=minf[n]-(minf[n]-m[n])*exp(-dt/mtau[n]);
h[n]=hinf[n]-(hinf[n]-h[n])*exp(-dt/htau[n]);
j[n]=jinf[n]-(jinf[n]-j[n])*exp(-dt/jtau[n]);

if (varmax<pow((m[n]- mtemp[n])*(m[n]- mtemp[n]),0.5)) varmax= pow((m[n]- mtemp[n])*(m[n]- mtemp[n]),0.5);
if (varmax<pow((h[n]- htemp[n])*(h[n]- htemp[n]),0.5)) varmax= pow((h[n]- htemp[n])*(h[n]- htemp[n]),0.5);
if (varmax<pow((j[n]- jtemp[n])*(j[n]- jtemp[n]),0.5)) varmax= pow((j[n]- jtemp[n])*(j[n]- jtemp[n]),0.5);

ENa=(R*T/(zNa*F))*log(Naout/Nain)*1000;	
return  GNa[n]*m[n]*m[n]*m[n]*j[n]*h[n]*(V[n]-ENa); // nA

}



