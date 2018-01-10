#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "funciones_biologicas.h"
#include "funciones_circuito.h"
#include "funciones_dinamicas.h"
#include "constantes.h"

//include's de c++
#include <string>
#include <sstream>
using namespace std;

template<class T>
inline string to_string(const T& t){
    stringstream ss;
    ss << t;
    return ss.str();
}

//#define TMAX 2000000.0
#define TTRANS 500.0
#define dt 0.1
#define DIM 2
#define v0 -31.19
#define w0 0.1227
#define ml_CONST 11
//#define IMAX 7
#define di 1
#define UMBRAL 20
#define TISI 12
#define VECTOR 0.5
#define dD 0.1

void espera(){
    int d;
    printf("Apriete el 0 para continuar.");
    scanf("%i", &d);
}

int main(/*int argc, char *argv[]*/) {

	double iext;
	double D;
	double TMAX,DMAX,IMAX,IMIN;
	printf("Escriba el tiempo de simulacion en ms\n");
	scanf("%lf",&TMAX);
	printf("Escriba el valor de corriente inyectada minima.\n");
	scanf("%lf",&IMIN);
	printf("Escriba el valor de corriente inyectada maxima.\n");
	scanf("%lf",&IMAX);
	printf("Escriba el valor de la amplitud del ruido maximo\n");
	scanf("%lf",&DMAX);
	/*iext=atof(argv[1]);
	D=atoi(argv[2]);
	TMAX = atoi(argv[3]);*/

	double *v;
	double vant;
	double t;
	double aleatorio;
	double *constantes_biologicas;
	double estimulo,frec;
	int testimulo,NSpikes;

    testimulo = VECTOR/dt;
 
    v = (double*)malloc(DIM*sizeof(double));
	constantes_biologicas = (double*)malloc(ml_CONST * sizeof(double));

	/*
	v[0]: potencial de membrana
	v[1]: apertura de potasio
	*/

    void RKS(double *v,double *constantes,double I,double aleatorio,double D);

    v[0] = v0;
    v[1] = w0;
	ml_init_array(constantes_biologicas);

    const gsl_rng_type * T;
    gsl_rng * r;
     
    double sigma = 1.0;
     
    gsl_rng_env_setup();
     
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    FILE *tasadisparo;
    tasadisparo = fopen("/home/meli/Maestria/Modelos/Simulaciones/tasadedisparo.dat","w");

	/*iext=0;
    while(iext <= IMAX) { 
      printf("iext = %lf\n",iext);
      D=0;
		while(D<=DMAX){ 
			printf(" D = %lf\n",D);
			t = 0.0;
			NSpikes = 0;  
			while(t < TMAX){
				vant = v[0];  
				aleatorio=gsl_ran_gaussian(r,sigma);
				estimulo=iext+aleatorio;
				RKS(v,constantes_biologicas,iext,aleatorio,D);
				if(v[0]>UMBRAL && vant<UMBRAL){
					NSpikes++;
				}
				t=t+dt;
			}
			
			frec = 1000*NSpikes*1.0/TMAX;
			fprintf(tasadisparo,"%lf\t",frec);
			D+=dD;
		}
		fprintf(tasadisparo,"\n");
		iext+=di;
	}*/
	
	D=0;
    while(D <= DMAX) { 
      printf("D = %lf\n",D);
      iext=IMIN;
		while(iext<=IMAX){ 
			printf(" iext = %lf\n",iext);
			t = 0.0;
			NSpikes = 0;  
			while(t < TMAX){
				vant = v[0];  
				aleatorio=gsl_ran_gaussian(r,sigma);
				estimulo=iext+aleatorio;
				RKS(v,constantes_biologicas,iext,aleatorio,D);
				if(v[0]>UMBRAL && vant<UMBRAL){
					NSpikes++;
				}
				t=t+dt;
			}
			
			frec = 1000*NSpikes*1.0/TMAX;
			fprintf(tasadisparo,"%lf\t",frec);
			iext+=di;
		}
		fprintf(tasadisparo,"\n");
		D+=dD;
	}
	
	fclose(tasadisparo);
	free(v);
	free(constantes_biologicas);
	gsl_rng_free (r);

	return 0;
}
//---------------------Runge-Kutta Estocastico----------------------------------

void RKS(double *v,double *constantes,double corriente,double aleatorio,double D){

    double *vaux;
    double *K1,*K2;
    int i;

    vaux = (double*)malloc(2*sizeof(double));
    K1 = (double*)malloc(2*sizeof(double));
    K2 = (double*)malloc(2*sizeof(double));

    for(i=0;i<2;i++){
        vaux[i]=v[i];
    }

	K1[0] = dt*ml_F0(vaux,constantes,corriente);
	K1[1] = dt*ml_F1(vaux,constantes);

    vaux[0] = v[0]+K1[0]+sqrt(2*D*dt)*aleatorio;
    vaux[1] = v[1]+K1[1];

    K2[0] = dt*ml_F0(vaux,constantes,corriente);
	K2[1] = dt*ml_F1(vaux,constantes);
    
	v[0]= v[0] + 0.5*(K1[0]+K2[0])+sqrt(2*D*dt)*aleatorio;
	v[1]= v[1] + 0.5*(K1[1]+K2[1]);

    free(vaux);
	free(K1);
	free(K2);
}

