
/*
 *  Main program
 *  
 *
 *  Created by  Gildas Loussouarn on Oct 28 2019.
 *  
 *
 */ 


#include <string.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include "INa_STEP.h"	
#include "WT_HHSCN5A_function4RK.cc" 
#include "Timecourse_Voltage_INa.cc" 			
	
using namespace std;	
		 
int main ()
{



double WT_INa_function (n);			
double Timecourse_Voltage_INa ();	
        

clock_start=clock();
Timecourse_Voltage_INa ();		
clock_end=clock();
cout << "Time used: "<<(clock_end-clock_start)/(double)CLOCKS_PER_SEC<<"s"<<endl;
        
return (1);

}

  
