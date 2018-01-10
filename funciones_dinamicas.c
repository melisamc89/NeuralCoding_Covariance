#include <stdio.h>
#include "funciones_circuito.h"
#include "funciones_biologicas.h"
#include "constantes.h"

//Modelo de Golomb

double golomb_F0(double *v,double *constantes,double I){
	double x;
	x = circuito4(golomb_cor_na,golomb_cor_k,golomb_cor_d,golomb_cor_l,v,constantes,I)/golomb_C;
	return x;
}

double golomb_F1(double *v,double *constantes){
	double x;
	x = dinamica1(golomb_hinf,golomb_tauh,v[0],constantes,v[1]);
	return x;
}

double golomb_F2(double *v,double *constantes){
	double x;
	x = dinamica1(golomb_ninf,golomb_taun,v[0],constantes,v[2]);
	return x;
}

double golomb_F3(double *v,double *constantes){
	double x;
	x = dinamica1(golomb_ainf,golomb_taua,v[0],constantes,v[3]);
	return x;
}

double golomb_F4(double *v,double *constantes){
	double x;
	x = dinamica1(golomb_binf,golomb_taub,v[0],constantes,v[4]);
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
    x = dinamica2(hh_alfam,hh_betam,v[0],constantes,v[1]);
    return  x;
}

double hh_F2(double *v,double *constantes){
	double x;
    x = dinamica2(hh_alfah,hh_betah,v[0],constantes,v[2]);
    return  x;
}

double hh_F3(double *v,double *constantes){
	double x;
    x = dinamica2(hh_alfan,hh_betan,v[0],constantes,v[3]);
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

