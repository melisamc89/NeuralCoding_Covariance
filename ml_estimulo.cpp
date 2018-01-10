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
#define IMAX 7
#define di 0.1
#define UMBRAL 20
#define TISI 12
#define VECTOR 0.5

void espera(){
    int d;
    printf("Apriete el 0 para continuar.");
    scanf("%i", &d);
}

int main(/*int argc, char *argv[]*/) {

	FILE *ptrace,*pspike,*info;

	double iext,D,TMAX;
	printf("Escriba el tiempo de simulacion en ms\n");
	scanf("%lf",&TMAX);
	printf("Escriba el valor de corriente inyectada.\n");
	scanf("%lf",&iext);
	printf("Escriba el valor de la amplitud del ruido\n");
	scanf("%lf",&D);
	/*iext=atof(argv[1]);
	D=atoi(argv[2]);
	TMAX = atoi(argv[3]);*/

	double *v;
	double vant,tant,interspike;
    double t,intervalo;
    double aleatorio;
	double *constantes_biologicas;
	double estimulo,frec;
	int spike,burst,testimulo,contador,NSpikes;
	string directorio1 = to_string("/home/meli/Maestria/Modelos/Simulaciones/");
	string modelo = to_string("ml2_");

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

	//while(iext<=IMAX){ 

		string output=directorio1+ modelo + to_string(iext) + to_string("_")+to_string(D)+to_string(".dat");
        if((pspike = fopen(output.c_str(), "w")) == NULL){
            printf("No puedo abrir el archivo %s.\n", output.c_str());
            exit(1);
        }

		printf("%lf\n",iext);
		
		/*string outputtrace=to_string("ml_") + to_string(iext) + to_string(".dat");
        if((ptrace = fopen(outputtrace.c_str(), "w")) == NULL){
            printf("No puedo abrir el archivo %s.\n", outputtrace.c_str());
            exit(1);
        }*/
		
		int iwrite=1;
		t = 0.0;
		tant=0.0;
		contador = 0; 
		NSpikes = 0;
		  
		while(t < TMAX){
			vant = v[0];
			if (iwrite == 1){
				spike = 0;
				iwrite=0;
				}
			burst=0;    
			aleatorio=gsl_ran_gaussian(r,sigma);
			estimulo=iext+aleatorio;
			RKS(v,constantes_biologicas,iext,aleatorio,D);
			//fprintf(ptrace, "%lf\t%lf\t%lf\n", t,v[0],v[1]);
			
			if(v[0]>UMBRAL && vant<UMBRAL){
				spike = 1;
				NSpikes++;
				interspike = tant-t;
				if(interspike>TISI){
					burst=1;				
				}
				tant=t;
			}

			if(contador == testimulo){
				estimulo = estimulo-iext;
				fprintf(pspike,"%lf\t%lf\t%lf\t%i\t%i\n",t,estimulo,v[0],spike,burst);
				contador = 0;
				iwrite=1;
			}
			contador++;
			t=t+dt;
		}
		
		string output_info=directorio1+ modelo + to_string(iext)+to_string("_")+to_string(D)+to_string("_info.dat");
        if((info = fopen(output_info.c_str(), "w")) == NULL){
            printf("No puedo abrir el archivo %s.\n", output_info.c_str());
            exit(1);
        }
        
        frec = 1000*NSpikes/TMAX;
        fprintf(info,"%lf\t%lf\t%lf\t%i\t%lf\t%lf\n",iext,D,TMAX,NSpikes,frec,VECTOR);
		iext=iext+di;
		fclose(info);
	//}
	
	free(v);
	//fclose(ptrace);
	fclose(pspike);
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

