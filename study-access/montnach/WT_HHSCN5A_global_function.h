/*
 *  SCN5A_global_function.h
 *  
 *
 *  Created by Gildas Loussouarn on Wed Jan 29 2019.
 *  .
 *
 */

//#include <Carbon/Carbon.h>


double zNa =1, F = 96485, R = 8.314, T = 293; // PK = 17.41e-8 cm.s-1 dans oocytes PMC1189096
double m[11],h[11],j[11],mtemp[11],htemp[11],jtemp[11],minf[11],hinf[11],jinf[11],mtau[11],htau[11],jtau[11];
double GNa[11] = {0.002,0.004,0.006,0.02,0.04,0.06,0.2,0.4,0.6,2,4};
//double GNa[11] = {0.1,0.15,1,1.5,0.04,0.06,0.2,0.4,0.6,2,4};
char filename[11][25],filename2[11];
