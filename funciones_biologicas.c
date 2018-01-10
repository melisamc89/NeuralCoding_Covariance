#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "constantes.h"

//Funciones del modelo golomb

void golomb_init_array(double *array){

	array[0]=golomb_GL;
	array[1]=golomb_VL;
	array[2]=golomb_GNA;
	array[3]=golomb_VNA;
	array[4]=golomb_GK;
	array[5]=golomb_VK;
	array[6]=golomb_GD;
	array[7]=golomb_THETA_M;
	array[8]=golomb_SIGMA_M;
	array[9]=golomb_THETA_H;
	array[10]=golomb_SIGMA_H;
	array[11]=golomb_THETA_TH;
	array[12]=golomb_SIGMA_TH;
	array[13]=golomb_THETA_N;
	array[14]=golomb_SIGMA_N;
	array[15]=golomb_THETA_A;
	array[16]=golomb_SIGMA_A;
	array[17]=golomb_TAU_A;
	array[18]=golomb_THETA_B;
	array[19]=golomb_SIGMA_B;
	array[20]=golomb_TAU_B;

}

void golomb_variables(double *array,int *variables,double *valor){
	
	int i;
	int n = sizeof(variables)/sizeof(double);

	if(n > golomb_CONST){
		printf("array modificador fuera de rango");}

	for(i=0;i<n;i++){
		array[variables[i]] = valor[i];
	}
}

double golomb_minf(double v,double *constantes){
	double x;
	x = 1./(1. + exp(-(v-constantes[7])/constantes[8]));
	return x;
}

double golomb_cor_na(double *v,double *constantes){
	double x;
	x = constantes[2]*golomb_minf(v[0],constantes)*golomb_minf(v[0],constantes)*golomb_minf(v[0],constantes)*v[1]*(v[0]-constantes[3]);
	return x;
}

double golomb_cor_k(double *v,double *constantes){
	double x;
	x = constantes[4]*v[2]*v[2]*(v[0]-constantes[5]);
	return x;
}

double golomb_cor_d(double *v,double *constantes){
	double x;
	x = constantes[6]*v[3]*v[3]*v[3]*v[4]*(v[0]-constantes[5]);
	return x;
}

double golomb_cor_l(double v,double *constantes){
	double x;
	x = constantes[0]*(v-constantes[1]);
	return x;
}

double golomb_hinf(double v,double *constantes){
	 return 1./(1+exp(-(v-constantes[9])/constantes[10]));
}

double golomb_tauh(double v,double *constantes){
	return 0.5+14/(1+exp(-(v-constantes[11])/constantes[12]));
}

double golomb_ninf(double v,double *constantes){
	return 1./(1+exp(-(v-constantes[13])/constantes[14]));
}

double golomb_taun(double v,double *constantes){
	double a,b,x;
	a = 0.087 + 11.4 / (1+exp((v+14.6)/8.6));
	b = 0.087 + 11.4 / (1+exp(-(v-13)/18.7));
	x = a*b;
	return x;
}

double golomb_ainf(double v,double *constantes){
	return 1./(1+exp(-(v-constantes[15])/constantes[16]));
}

double golomb_taua(double v,double *constantes){
	return constantes[17];
}

double golomb_binf(double v,double *constantes){
	return 1./(1+exp(-(v-constantes[18])/constantes[19]));
}

double golomb_taub(double v,double *constantes){
	return constantes[20];
}

//Funciones del modelo Hodgkin-Huxley

void hh_init_array(double *array){

	array[0]=hh_GL;
	array[1]=hh_VL;
	array[2]=hh_GNA;
	array[3]=hh_VNA;
	array[4]=hh_GK;
	array[5]=hh_VK;

}


double hh_cor_na(double *v,double *constantes){
	return constantes[2]*v[1]*v[1]*v[1]*v[2]*(v[0]-constantes[3]);
}

double hh_cor_k(double *v,double *constantes){
	return constantes[4]*v[3]*v[3]*v[3]*v[3]*(v[0]-constantes[5]);
}

double hh_cor_l(double v,double *constantes){
	return constantes[0]*(v-constantes[1]);
}

double hh_minf(double v){
	return hh_alfam(v)/(alfam(v)+betam(v));
	}
	
double hh_hinf(double v){
	return hh_alfah(v)/(alfah(v)+betah(v));
	}
	
double hh_ninf(double v){
	return hh_alfan(v)/(alfan(v)+betan(v));
	}

double hh_tm(double v){
	return 1/(alfam(v)+betam(v));
	}
	
double hh_th(double v){
	return 1/(alfah(v)+betah(v));
	}
	
double hh_tn(double v){
	return 1/(alfan(v)+betan(v));
	}

/*double hh_alfam(double pot){
    return 0.1*(40+pot)/(1.0-exp(-(40+pot)/10));
}

double hh_betam(double pot){
    return 4.0 * exp(-(pot+65)/18.0);
}

double hh_alfah(double pot){
    return 0.07 * exp(-(pot+65) / 20.0);
}

double hh_betah(double pot){
    return 1.0 /(1.0+exp(-(35+pot)/10));
}

double hh_alfan(double pot){
    return  0.01*(pot+55)/ (1.0-exp(-(pot+55)/10));
}

double hh_betan(double pot){
    return 0.125 * exp(-(pot+35) / 80.0);
}

*/
double hh_alfam(double pot){
    return (2.5 - 0.1 * pot) / (exp(2.5 - 0.1 * pot) - 1.0);
}

double hh_betam(double pot){
    return 4.0 * exp(-pot/18.0);
}

double hh_alfah(double pot){
    return 0.07 * exp(-pot / 20.0);
}

double hh_betah(double pot){
    return 1.0 / (exp(3.0 - 0.1 * pot) + 1.0);
}

double hh_alfan(double pot){
    return  (0.1 - 0.01 * pot)/ (exp(1.0 - 0.1 * pot) - 1.0);
}

double hh_betan(double pot){
    return 0.125 * exp(-pot / 80.0);
}

//Funciones del modelo de Morris-Lecar

void ml_init_array(double *array){

	array[0]=ml_GCA;
	array[1]=ml_VCA;
	array[2]=ml_GK;
	array[3]=ml_VK;
	array[4]=ml_GL;
	array[5]=ml_VL;
	array[6]=ml_V1;
	array[7]=ml_V2;
	array[8]=ml_V3;
	array[9]=ml_V4;
	array[10]=ml_phi;
}

double ml_minf(double v,double *constantes);

double ml_cor_k(double *v,double *contantes){
	double x;
	x = contantes[2]*v[1]*(v[0]-contantes[3]);
	return x;
}

double ml_cor_ca(double *v,double *constantes){
	double x;
	x = constantes[0]*ml_minf(v[0],constantes)*(v[0]-constantes[1]);
	return x;
}

double ml_cor_l(double v,double *constantes){
	double x;
	x = constantes[4]*(v-constantes[5]);
	return x;
}

double ml_minf(double v,double *constantes){
	return 0.5*(1.0+tanh((v-constantes[6])/constantes[7]));
}

double ml_winf(double v,double *constantes){
	return 0.5*(1.0+tanh((v-constantes[8])/constantes[9]));
}

double ml_tw(double v,double *constantes){
	return 1.0/cosh((v-constantes[8])/(2.0*constantes[9]));
}

