int  Get_steady_state;

int n;
bool timeadapt=true;
double t,t1,dt,dtlimit[5];
double V[15],Vtemp[15],dtmax=10;
double back=0, accelerate=0, varmax=0, varmaxmax=1e-5;
double clock_start,clock_end;
double Rs = 5;// MOhm
double Cm = 20;//pF
double Nain=10.0, Naout= 145.0, ENa;  	//  Na+ concentration
