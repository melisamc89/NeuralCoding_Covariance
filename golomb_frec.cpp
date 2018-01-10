#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include "funciones_biologicas.h"
#include "funciones_circuito.h"
#include "funciones_dinamicas.h"
#include "constantes.h"
#include "gauss.h"

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

#define TMAX 50000.0
#define TTRANS 50
#define UMBRAL 20
#define IMAX 10.0
#define dt 0.1
#define di 0.1
#define DIM 5
#define TISI 12
#define VECTOR 0.5

#define V0 0.02
#define H0 0.02
#define N0 0.02
#define A0 0.02
#define B0 0.02
#define golomb_CANTIDADES_VARIABLES 2
#define golomb_CONST 21

void espera(){
        int d;
        printf("Apriete el 0 para continuar.");
        scanf("%d", &d);
}

void RK2(double *v,double *constantes,double corriente);

int main(){
	
	FILE *ptrace,*pspike,*info,*frecuencia;

	double iext=0,sigma;
	printf("Escriba el valor de la amplitud del ruido.\n");
	scanf("%lf",&sigma);

	double *v;
	double vant,t;
	double *constantes_biologicas,*valores;
	int *posiciones;
	double estimulo;
	int spike,testimulo,contador,NSpikes;
	double frec;
	string directorio1 = to_string("/home/meli/Documentos/Maestria en Fisica/Modelos/Simulaciones/");

    v = (double *) malloc (DIM * sizeof(double));
	posiciones = (int*) malloc (golomb_CANTIDADES_VARIABLES * sizeof(int));
	valores = (double*)malloc(golomb_CANTIDADES_VARIABLES * sizeof(double));
	constantes_biologicas = (double*)malloc(golomb_CONST*sizeof(double));
	
    v[0] = V0;
    v[1] = H0;
    v[2] = N0;
    v[3] = A0;
    v[4] = B0;
	golomb_init_array(constantes_biologicas);    

    string outputfrec=directorio1+ to_string("g_FI.dat"); 
    if((frecuencia = fopen(outputfrec.c_str(), "w")) == NULL){
            printf("No puedo abrir el archivo curva FI\n");
            exit(1);
        }

	for(int i=0;i<2*IMAX/di;i++){ 
		
		if(i<IMAX/di)
			iext=i*di;
		if(i>=IMAX/di)
			iext=2*IMAX-i*di;
		printf("%lf\n",iext);

		t = 0.0;
		NSpikes = 0;

		while(t < TMAX){
			vant = v[0];
			RK2(v,constantes_biologicas,iext);
			if(v[0]>UMBRAL && vant<UMBRAL && t>TTRANS){
				spike = 1;
				NSpikes++;
			}
			t=t+dt;
		}
		
        frec = 1000*NSpikes/(TMAX-TTRANS);
        fprintf(frecuencia, "%lf\t%lf\n", iext,frec);
	}
	
	free(v);
	fclose(frecuencia);
	free(posiciones);
	free(valores);
	free(constantes_biologicas);	

    return 0;
}

/*.............Algoritmo de integracion...Runge-Kutta estocastico..............*/

void RK2(double *v,double *constantes,double corriente){

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
    
    vaux[0] = v[0]+K1[0];
    vaux[1] = v[1]+K1[1];
    vaux[2] = v[2]+K1[2];
   	vaux[3] = v[3]+K1[3];
   	vaux[4] = v[4]+K1[4];

   	K2[0] = dt*golomb_F0(vaux,constantes,corriente);
	K2[1] = dt*golomb_F1(vaux,constantes);
	K2[2] = dt*golomb_F2(vaux,constantes);
	K2[3] = dt*golomb_F3(vaux,constantes);
	K2[4] = dt*golomb_F4(vaux,constantes);

	v[0]= v[0] + 0.5*(K1[0]+K2[0]);
	v[1]= v[1] + 0.5*(K1[1]+K2[1]);
	v[2]= v[2] + 0.5*(K1[2]+K2[2]);
	v[3]= v[3] + 0.5*(K1[3]+K2[3]);
	v[4]= v[4] + 0.5*(K1[4]+K2[4]);

    free(vaux);
	free(K1);
	free(K2);
}


