#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include "funciones_biologicas.h"
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

#define TMAX 500000.0
#define TTRANS 50
#define UMBRAL -20
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

#define TAU 10
#define I0 0

void espera(){
        int d;
        printf("Apriete el 0 para continuar.");
        scanf("%d", &d);
}

void RK2(double *v,double *constantes,double corriente);

int main(){
	
	FILE *ptrace,*pspike,*info,*frecuencia;

	double iext,sigma;
	/*printf("Escriba el tiempo de simulacion en ms\n");
	scanf("%lf",&TMAX);*/
	printf("Escriba el valor de corriente inyectada.\n");
	scanf("%lf",&iext);
	printf("Escriba el valor de la amplitud del ruido.\n");
	scanf("%lf",&sigma);


	double *v;
	double vant,t;
	double *constantes_biologicas,*valores;
	int *posiciones;
	double estimulo;
	int spike,testimulo,contador,NSpikes;
	double frec;
	string directorio1 = to_string("/home/meli/Documentos/Maestria en Fisica/Simulaciones/");

	testimulo=VECTOR/dt;

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

    /*string outputfrec=directorio1+ to_string("g_FI.dat"); 
    if((frecuencia = fopen(outputfrec.c_str(), "w")) == NULL){
            printf("No puedo abrir el archivo curva FI\n");
            exit(1);
        }*/

	/*for(int i=0;i<2*IMAX/di;i++){ 
		
		if(i<IMAX/di)
			iext=i*di;
		if(i>=IMAX/di)
			iext=2*IMAX-i*di;
		printf("%lf\n",iext);*/
		
		
		string output=directorio1+to_string("g_dinam_") + to_string(iext) + to_string("_")+to_string(sigma)+to_string(".dat");
			if((ptrace = fopen(output.c_str(), "w")) == NULL){
				printf("No puedo abrir el archivo %s.\n", output.c_str());
				exit(1);
			}
			
		string output2=directorio1+ to_string("g_stim_")+ to_string(iext) + to_string("_")+to_string(sigma)+to_string(".dat");
			if((pspike = fopen(output2.c_str(), "w")) == NULL){
				printf("No puedo abrir el archivo %s.\n", output.c_str());
				exit(1);
			}
			
		string output_info=directorio1+ to_string("g_inf_") + to_string(iext)+to_string("_")+to_string(sigma)+to_string(".dat");
			if((info = fopen(output_info.c_str(), "w")) == NULL){
				printf("No puedo abrir el archivo %s.\n", output_info.c_str());
				exit(1);
			}
        
		int iwrite=1;
		t = 0.0;
		contador = 0; 
		NSpikes = 0;

		while(t < TMAX){
			vant = v[0];
			if (iwrite == 1){
				spike = 0;
				iwrite=0;
				}    
			//iext=iext*(1-dt/TAU)+I0*dt/TAU+Gauss(sigma)*sqrt(dt)*sigma;
			RK2(v,constantes_biologicas,iext);
			fprintf(ptrace, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t,v[0],v[1],v[2],v[3],v[4]);
			if(v[0]>UMBRAL && vant<UMBRAL && t>TTRANS){
				spike = 1;
				NSpikes++;
			}
			if(contador == testimulo){
				estimulo = estimulo-iext;
				fprintf(pspike,"%lf\t%lf\t%lf\t%i\n",t,estimulo,v[0],spike);
				contador = 0;
				iwrite=1;
			}
			contador++;
			t=t+dt;
		}
		
        frec = 1000*NSpikes/(TMAX-TTRANS);
        fprintf(info,"%lf\t%lf\t%lf\t%i\t%lf\t%lf\n",iext,sigma,TMAX,NSpikes,frec,dt);
        //fprintf(frecuencia, "%lf\t%lf\n", iext,frec);
		fclose(info); 
	//}
	
	free(v);
	fclose(ptrace);
	fclose(pspike);
	//fclose(frecuencia);
	free(posiciones);
	free(valores);
	free(constantes_biologicas);	

    return 0;
}

/*.............Algoritmo de integracion...Runge-Kutta..............*/

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


