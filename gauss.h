#ifndef FUNCIONESGAUSS_H
#define FUNCIONESGAUSS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double Gauss();

double Gauss(double sigma){
	
	double aleat,p,r;
	int prueba=0;
	
	while(prueba==0){
	aleat=rand()*1.0/RAND_MAX;
	p=exp(-aleat*aleat/(2*sigma*sigma))/(4*sigma*asin(1));
	r=rand()*1.0/RAND_MAX;
	if(r<p)
		prueba=1;
	}
	return aleat;
	}
	
#endif
