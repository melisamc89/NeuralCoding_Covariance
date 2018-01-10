#ifndef FUNCIONESCIR_H
#define FUNCIONESCIR_H

double circuito3(double (*I1)(double*,double*),double (*I2)(double*,double*),double (*I4)(double,double*),double *v,double *constantes,double I);
double circuito4(double (*I1)(double*,double*),double (*I2)(double*,double*),double (*I3)(double*,double*),double (*I4)(double,double*),double *v,double *constantes,double I);
double dinamica1(double (*estacionario)(double,double*),double (*tau)(double,double*),double v,double *constantes,double variable);
double dinamica2(double(*alfa)(double),double(*beta)(double),double v,double *constantes,double variable);

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

double circuito3(double (*I1)(double*,double*),double (*I2)(double*,double*),double (*I3)(double,double*),double *v,double *constantes,double I){
	double x;
	x = (-I1(v,constantes)-I2(v,constantes)-I3(v[0],constantes)+I);
	return x;
}

double circuito4(double (*I1)(double*,double*),double (*I2)(double*,double*),double (*I3)(double*,double*),double (*I4)(double,double*),double *v,double *constantes,double I){
	double x;
	x = (-I1(v,constantes)-I2(v,constantes)-I3(v,constantes)-I4(v[0],constantes)+I);
	return x;
}

double dinamica1(double (*estacionario)(double,double*),double (*tau)(double,double*),double v,double *constantes,double variable){
	double x;
	x = (estacionario(v,constantes)-variable)/tau(v,constantes);
	return x;
}

double dinamica2(double(*alfa)(double),double(*beta)(double),double v,double *constantes,double variable){
	double x;
	x = alfa(v)*(1-variable)-beta(v)*variable;
	return x;
}

#endif

