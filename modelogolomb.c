#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include "funciones_biologicas.h"
#include "funciones_circuito.h"
#include "funciones_dinamicas.h"
#include "constantes.h"

#define TMAX 3000
#define NSPIKES 100
#define IMAX 4.2
#define TTRANS 0
#define UMBRAL -20
#define dt 0.5
#define di 0.1
#define DIM 5

#define V0 0.02
#define H0 0.02
#define N0 0.02
#define A0 0.02
#define B0 0.02
#define golomb_CANTIDADES_VARIABLES 2
#define golomb_CONST 21

#define outputtrace "traza.dat"
#define outputfrec "frec.dat"

void RKS(double *v,double *constantes,double corriente,double aleatorio, double D);

int main(){

    FILE *ptrace,*pfrec;

    double *V,Vant,iext,frec,t;
	int spike;
	int *posiciones;
	double *constantes_biologicas,*valores;
	double random=0;
	double D=0;

    V = (double *) malloc (DIM * sizeof(double));
	posiciones = (int*) malloc (golomb_CANTIDADES_VARIABLES * sizeof(int));
	valores = (double*)malloc(golomb_CANTIDADES_VARIABLES * sizeof(double));
	constantes_biologicas = (double*)malloc(golomb_CONST*sizeof(double));

    V[0] = V0;
    V[1] = H0;
    V[2] = N0;
    V[3] = A0;
    V[4] = B0;
	golomb_init_array(constantes_biologicas);    

    if((ptrace = fopen(outputtrace, "w")) == NULL){
        printf("No puedo abrir el archivo %s.\n", outputtrace);
        exit(1);
    }
    
    if((pfrec = fopen(outputfrec, "w")) == NULL){
        printf("No puedo abrir el archivo %s.\n", outputfrec);
        exit(1);
    }
    
	iext=IMAX;
	int s=1;
	while(iext<=IMAX){

		t=0.;
		spike=0;

		while(t < TMAX ){
			Vant = V[0];
			RKS(V,constantes_biologicas,iext,random,D);
			if(V[0]>UMBRAL && Vant<UMBRAL)
				spike++;
			fprintf(ptrace,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, V[0], V[1], V[2], V[3],V[4]);
			t += dt;
			}
		
		/*frec = 1000.0*spike/(TMAX-TTRANS);
		fprintf(pfrec,"%lf\t%lf\n", iext,frec);
		if(iext==IMAX)
			s=0;
		if(s==1)
			iext+=di;
		else*/
			iext+=di;
		}

    free(V);
	free(posiciones);
	free(valores);
	free(constantes_biologicas);	
	fclose(pfrec);
    fclose(ptrace);

    return 0;
}

/*.............Algoritmo de integracion...Runge-Kutta estocastico..............*/


void RKS(double *v,double *constantes,double corriente,double aleatorio,double D){

    double *vaux;
    double *K1,*K2;
    int i;

   	vaux = (double*)malloc(DIM*sizeof(double));
   	K1 = (double*)malloc(DIM*sizeof(double));
    K2 = (double*)malloc(DIM*sizeof(double));

    for(i=0;i<DIM;i++){
       	vaux[i]=v[i];
    }

	K1[0] = dt*golomb_F0(vaux,constantes,corriente);
	K1[1] = dt*golomb_F1(vaux,constantes);
	K1[2] = dt*golomb_F2(vaux,constantes);
	K1[3] = dt*golomb_F3(vaux,constantes);
	K1[4] = dt*golomb_F4(vaux,constantes);

    vaux[0] = v[0]+K1[0]+sqrt(2*D*dt)*aleatorio;
    vaux[1] = v[1]+K1[1];
    vaux[2] = v[2]+K1[2];
   	vaux[3] = v[3]+K1[3];
   	vaux[4] = v[4]+K1[4];

   	K2[0] = dt*golomb_F0(vaux,constantes,corriente);
	K2[1] = dt*golomb_F1(vaux,constantes);
	K2[2] = dt*golomb_F2(vaux,constantes);
	K2[3] = dt*golomb_F3(vaux,constantes);
	K2[4] = dt*golomb_F4(vaux,constantes);
    
	v[0]= v[0] + 0.5*(K1[0]+K2[0])+sqrt(2*D*dt)*aleatorio;
	v[1]= v[1] + 0.5*(K1[1]+K2[1]);
	v[2]= v[2] + 0.5*(K1[2]+K2[2]);
	v[3]= v[3] + 0.5*(K1[3]+K2[3]);
	v[4]= v[4] + 0.5*(K1[4]+K2[4]);

    free(vaux);
	free(K1);
	free(K2);
}
