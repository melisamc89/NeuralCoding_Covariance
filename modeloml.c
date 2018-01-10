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

#define TMAX 500.0
#define dt 0.1
#define DIM 2
#define v0 -31.19
#define w0 0.1227
#define ml_CONST 11
//#define D 3.5

#define outputtrace "traza.dat"

void espera(){
    int d;
    printf("Apriete el 0 para continuar.");
    scanf("%i", &d);
}

int main(){

	FILE *ptrace,*tiempo;

	double *v;
    double t,intervalo;
    double aleatorio;
    double D=0;
	double *constantes_biologicas;
    
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
    
    if((ptrace = fopen(outputtrace, "w")) == NULL){
        printf("No puedo abrir el archivo %s.\n", outputtrace);
        espera();
        exit(1);
    }

    /*string outputtrace = to_string("tr_s_") + to_string(iext) + to_string(".dat");

        if((ptrace = fopen(outputtrace.c_str(), "w")) == NULL)
        {
            printf("No puedo abrir el archivo %s.\n", outputtrace.c_str());
            exit(1);
        }*/


    const gsl_rng_type * T;
    gsl_rng * r;
     
    double sigma = 1.0;
     
    gsl_rng_env_setup();
     
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

	t = 0.0;   
    while(t < TMAX){    
		aleatorio=0;
	    //aleatorio=gsl_ran_gaussian(r,sigma);
	    RKS(v,constantes_biologicas,IEXT,aleatorio,D);
		fprintf(ptrace, "%lf\t%lf\t%lf\n", t,v[0],v[1]);
	    t=t+dt;
	}
	 
	 free(v);
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

