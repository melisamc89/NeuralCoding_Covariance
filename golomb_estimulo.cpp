#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
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

//#define TMAX 500000.0
#define TTRANS 0
#define UMBRAL 20
#define IMAX 16
#define dt 0.1
#define di 1
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

int main(){
	
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
	double *constantes_biologicas,*valores;
	int *posiciones;
	double estimulo;
	int spike,burst,testimulo,contador,NSpikes,Nburst;
	double frec;
	string directorio1 = to_string("/home/meli/Documentos/Maestria en Fisica/Modelos/Simulaciones/");
	string modelo = to_string("g_");


    v = (double *) malloc (DIM * sizeof(double));
	posiciones = (int*) malloc (golomb_CANTIDADES_VARIABLES * sizeof(int));
	valores = (double*)malloc(golomb_CANTIDADES_VARIABLES * sizeof(double));
	constantes_biologicas = (double*)malloc(golomb_CONST*sizeof(double));

    /*
    V[0]: potencial de membrana
    V[1]: h
    V[2]: n
    V[3]: a
    V[4]: b */

	void RKS(double *v,double *constantes,double corriente,double aleatorio,double D);

    v[0] = V0;
    v[1] = H0;
    v[2] = N0;
    v[3] = A0;
    v[4] = B0;
	golomb_init_array(constantes_biologicas);    

    /*const gsl_rng_type * T;
    gsl_rng * r;
     
    double sigma = 1.0;
     
    gsl_rng_env_setup();
     
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);*/

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
			//aleatorio=0
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
	free(posiciones);
	free(valores);
	free(constantes_biologicas);
	gsl_rng_free (r);	

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


