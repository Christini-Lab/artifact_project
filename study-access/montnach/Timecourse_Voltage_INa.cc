/*
 *  Current_Voltage.cc
 *  
 *
 *  Created by Gildas Loussouarn on June 06 2020.
 * 
 *
 */

//**************************************************************************************
//	FUNCTION TO COMPUTE A CURRENT in response to VOLTAGE 
//**************************************************************************************

double Timecourse_Voltage_INa ()					//start main program

{
long int goprint,new1,j,numerostep=0;
double vhold=-100, Rm[15], im[15], Vcmd, Vmem, is[15],cherchemin[15][30],nbstep=30;
double startstep=1, endstep=51, endprot=(1+nbstep)*2000; 
FILE *dataout[15],*dataout2;
accelerate=0;
back=0;
dt=1e-4;		
dtmax=10;      // time step size
t = 0.0;  
t1 = 0.0;        // Time (ms) 
  
for (n=0;n<11;n+=1) 
	{
	im[n]= 0.0;	
	is[n]=0.0;
	for (j=0;j<31;j+=1) cherchemin[n][j]=0.0;
	sprintf(filename[n],"INabad_%.3f_Rs%.1f_Cm%.1f",GNa[n],Rs,Cm);
	sprintf(filename2,"tableau des Imin");
	std::cout << "creation of: "<<filename[n]<<std::endl;
	dataout[n] = fopen(filename[n], "w");	// open data file to write 	
	}
	dataout2 = fopen(filename2, "w");	// open data file to write
	
do 		//time loop
	{
	t=t1;
	if (timeadapt==true) 
		{
		if (accelerate == 1) dt=dt*2;	
		for (j=1;j<=4;j+=1) dtlimit[j]=dtmax;
		if (floor(t+dt) != floor(t)) dtlimit[1]=(1+(int)t)-t;  //not to miss a 1 time
		if (floor((t+dt)*20) != floor(t*20)) dtlimit[2]=0.05*(1+(int)(t*20))-t; //not to miss a 0.05 time
		if (fmod(t,2000)<startstep && fmod(t+dt,2000)>startstep) dtlimit[3]= startstep - fmod(t,2000);
		if (fmod(t,2000)<endstep && fmod(t+dt,2000)>endstep) dtlimit[4]= endstep - fmod(t,2000);
		for (j=1;j<=4;j+=1) dt = std::min(dt, dtlimit[j]);
		}

	if (t+dt>2000 && fmod(t+dt,2000)>=startstep && fmod(t+dt,2000)<=endstep) Vcmd=(-80+5*floor((t+dt-2000)/2000)); 
		else 
		{
		Vcmd=vhold;
		new1=1;
		}

	do
		{
		varmax = 0;
		accelerate = 0;
		for (n=0;n<11;n+=1) 
			{
			if (t==0) V[n] = vhold;
			if (back == 0)  Vtemp[n] = V[n]; else V[n] = Vtemp[n];
			im[n]=WT_INa_function(n);
			if (im[n] < cherchemin[n][numerostep])  cherchemin[n][numerostep] = im[n];
			if (Rs != 0) 
				{
				V[n] = V[n] + (((Vcmd-V[n])*1e-3)/(Rs*1e6 * Cm*1e-12) - (im[n]*1e-9)/(Cm*1e-12))*1000 * (dt*1e-3); // V en mV - Ra en MOhm - im en nA
				is[n] = (Vcmd-V[n])*1e-3/(Rs*1e6)*1e9;// is, im  et dis en nA -  Cm en pF - Ra en MOhm - dt en ms
				if (varmax<pow((V[n]- Vtemp[n])/100*(V[n]- Vtemp[n])/100,0.5)) varmax= pow((V[n]- Vtemp[n])/100*(V[n]- Vtemp[n])/100,0.5);
				}
			if (Rs == 0) 
				{
				is[n] = im[n];				
				V[n] = Vcmd;
				}
			}
		back=0;
		if (varmax>varmaxmax && timeadapt==true) 
			{
			back=1;
			dt = dt/2;
			}
		else if (varmax<varmaxmax/2 && timeadapt==true)	accelerate=1;
		}
	while (back==1 && timeadapt==true);		
	
	if (new1 == 1 && Vcmd != vhold && Vmem != Vcmd) 
		{
		std::cout << "Depol to: "<<Vcmd<<" mV"<<std::endl;
		numerostep+=1;
		Vmem=Vcmd;
		new1=0;
		}		
 	goprint=0;	
	t1 = t+dt;
	if (floor(t1*20) != floor(t*20) && t1>=2000 && fmod(t1,2000)>=0 && fmod(t1,2000)<=6) goprint=1; 
	if (floor(t1) != floor(t) && t1>=2000 && fmod(t1,2000)>=6 && fmod(t1,2000)<=100.9) goprint=1; 	
	for (n=0;n<11;n+=1) if (goprint ==1) fprintf(dataout[n], "%.2f\t %g\t %g\t %g\t %g\t %g\t\r", t+dt, Vcmd, V[n], im[n], is[n], dt);		
	}
while (t1<=endprot);
for (n=0;n<11;n+=1) 
	{
	for (j=1;j<31;j+=1) fprintf(dataout2, "%g\t" , cherchemin[n][j]);
	fprintf(dataout2, "\r");
	fclose (dataout[n]);
	}
	fclose (dataout2);
return (0);
}

   
