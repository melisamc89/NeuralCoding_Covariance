#ifndef FUNCIONESDIN_H
#define FUNCIONESDIN_H

//Modelo de Golomb

double golomb_F0(double *v,double *constantes,double I);
double golomb_F1(double *v,double *constantes);
double golomb_F2(double *v,double *constantes);
double golomb_F3(double *v,double *constantes);
double golomb_F4(double *v,double *constantes);

//Modelo de Hodgkin-Huxley

double hh_F0(double *v,double *constantes,double I);
double hh_F1(double *v,double *constantes);
double hh_F2(double *v,double *constantes);
double hh_F3(double *v,double *constantes);

//Modelo de Morris-Lecar

double ml_F0(double *v, double *constantes,double I);
double ml_F1(double *v,double *constantes);


#endif

#include <stdio.h>
#include "funciones_circuito.h"
#include "funciones_biologicas.h"
#include "constantes.h"

//Modelo de Golomb

double golomb_F0(double *v,double *constantes,double I){
	double x;
	x = -(golomb_cor_na(v,constantes)+golomb_cor_k(v,constantes)+golomb_cor_d(v,constantes)+golomb_cor_l(v[0],constantes)-I)/golomb_C;
	return x;
}

double golomb_F1(double *v,double *constantes){
	double x;
	x = (golomb_hinf(v[0],constantes)-v[1])/golomb_tauh(v[0],constantes);
	return x;
}

double golomb_F2(double *v,double *constantes){
	double x;
	x = (golomb_ninf(v[0],constantes)-v[2])/golomb_taun(v[0],constantes);
	return x;
}

double golomb_F3(double *v,double *constantes){
	double x;
	x = (golomb_ainf(v[0],constantes)-v[3])/golomb_taua(v[0],constantes);
	return x;
}

double golomb_F4(double *v,double *constantes){
	double x;
	x = (golomb_binf(v[0],constantes)-v[4])/golomb_taub(v[0],constantes);
	return x;
}

//Modelo de Hodgkin-Huxley

double hh_F0(double *v,double *constantes,double I){
    double x;
    x =circuito3(hh_cor_na,hh_cor_k,hh_cor_l,v,constantes,I)/hh_C;
    return  x;
}

double hh_F1(double *v,double *constantes){
	double x;
	x = dinamica1(hh_minf,hh_tm,v[0],constantes,v[1]);
    //x = dinamica2(hh_alfam,hh_betam,v[0],constantes,v[1]);
    return  x;
}

double hh_F2(double *v,double *constantes){
	double x;
	x = dinamica1(hh_hinf,hh_tm,v[0],constantes,v[2]);
    //x = dinamica2(hh_alfah,hh_betah,v[0],constantes,v[2]);
    return  x;
}

double hh_F3(double *v,double *constantes){
	double x;
	x = dinamica1(hh_ninf,hh_tn,v[0],constantes,v[3]);
    //x = dinamica2(hh_alfan,hh_betan,v[0],constantes,v[3]);
    return  x;
}

//Modelo de Morris-Lecar


double ml_F0(double *v,double *constantes,double corriente){
    double x;
    x=circuito3(ml_cor_k,ml_cor_ca,ml_cor_l,v,constantes,corriente)/ml_C;
    return x;
}

double ml_F1(double *v,double *constantes){
	 double x;
	 x=constantes[10]*dinamica1(ml_winf,ml_tw,v[0],constantes,v[1]);
	 return x;
}
