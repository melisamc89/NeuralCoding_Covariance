#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include "funciones_biologicas.h"
#include "funciones_circuito.h"
#include "funciones_dinamicas.h"
#include "constantes.h"

#define TMAX 1000.0// en milisegundos
#define TTRANS 300.0
#define UMBRAL 20
#define dt 0.1
#define di 0.1
#define v0 5.354558
#define m0 0.097361
#define h0 0.406031
#define n0 0.401889
#define IMAX 15
#define DIM 4
#define hh_CONST 6

#define outputtrace "hh.dat"
#define outputfrec "frec.dat"

void espera(){
    int d;
    printf("Apriete el 0 para continuar.");
    scanf("%i", &d);
}

int main(){

    FILE *pfrec,*ptrace;

    double *V;
    double aleatorio=0;
    double D=0;
	double t;
	double *constantes_biologicas;
	int contador;
	double i;
	double vant;
	double frec;

    V = (double *) malloc (DIM * sizeof(double));
	constantes_biologicas = (double*)malloc(hh_CONST * sizeof(double));

	void RKS(double *v,double *constantes,double I,double aleatorio,double D);

    V[0] = v0;
    V[1] = m0;
    V[2] = h0;
    V[3] = n0;
	hh_init_array(constantes_biologicas);    

    if((ptrace = fopen(outputtrace, "w")) == NULL){
        printf("No puedo abrir el archivo %s.\n", outputtrace);
        espera();
        exit(1);
    }

    if((pfrec = fopen(outputfrec, "w")) == NULL){
        printf("No puedo abrir el archivo %s.\n", outputfrec);
        espera();
        exit(1);
    }
        
	i = IMAX;
	//while(i<=IMAX){
		contador=0;	
    	t = 0.0;
		while(t < TMAX ){
			vant = V[0];
			RKS(V,constantes_biologicas,i,aleatorio,D);
			if (t>TTRANS && V[0]>UMBRAL && vant<UMBRAL)
				contador++;
		    fprintf(ptrace, "%f\t%f\t%f\t%f\t%f\n", t, V[0], V[1], V[2], V[3]);
		    t += dt;
		//}
	frec =1000*contador/(TMAX-TTRANS);
	fprintf(pfrec,"%lf\t%lf\n",i,frec);
	i+=di;
	}

    free(V);
    //fclose(ptrace);
	fclose(pfrec);

    return 0;
}


void RKS(double *v,double *constantes,double corriente,double aleatorio,double D){

    double *vaux;
    double *K1,*K2;
    int i;

    vaux = (double*)malloc(4*sizeof(double));
    K1 = (double*)malloc(4*sizeof(double));
    K2 = (double*)malloc(4*sizeof(double));

    for(i=0;i<4;i++){
        vaux[i]=v[i];
    }

	K1[0] = dt*hh_F0(vaux,constantes,corriente);
	K1[1] = dt*hh_F1(vaux,constantes);
	K1[2] = dt*hh_F2(vaux,constantes);
	K1[3] = dt*hh_F3(vaux,constantes);

    vaux[0] = v[0]+K1[0]+sqrt(2*D*dt)*aleatorio;
    vaux[1] = v[1]+K1[1];
    vaux[2] = v[2]+K1[2];
    vaux[3] = v[3]+K1[3];

    K2[0] = dt*hh_F0(vaux,constantes,corriente);
	K2[1] = dt*hh_F1(vaux,constantes);
	K2[2] = dt*hh_F2(vaux,constantes);
	K2[3] = dt*hh_F3(vaux,constantes);
    
	v[0]= v[0] + 0.5*(K1[0]+K2[0])+sqrt(2*D*dt)*aleatorio;
	v[1]= v[1] + 0.5*(K1[1]+K2[1]);
	v[2]= v[2] + 0.5*(K1[2]+K2[2]);
	v[3]= v[3] + 0.5*(K1[3]+K2[3]);

    free(vaux);
	free(K1);
	free(K2);
}


