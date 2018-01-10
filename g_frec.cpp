#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include "random.h"
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


//#define TMAX 25000000.0
#define NSPIKES 100
#define TTRANS 0
#define UMBRAL -20
#define dt 0.01
#define di 0.01
#define DIM 5
#define VECTOR 0.5

#define golomb_C 1.0
#define golomb_GL 0.25
#define golomb_VL -70
#define golomb_GNA 112.5
#define golomb_VNA 50
#define golomb_GK 225
#define golomb_VK -90
#define golomb_GD 0
#define golomb_THETA_M -24
#define golomb_SIGMA_M 11.5
#define golomb_THETA_H -58.3
#define golomb_SIGMA_H -6.7
#define golomb_THETA_TH -27
#define golomb_SIGMA_TH -15
#define golomb_THETA_N -12.4
#define golomb_SIGMA_N 6.8
#define golomb_THETA_A -50
#define golomb_SIGMA_A 20
#define golomb_TAU_A 2
#define golomb_THETA_B -70
#define golomb_SIGMA_B -6
#define golomb_TAU_B 150

#define V0 0.02
#define H0 0.02
#define N0 0.02
#define A0 0.02
#define B0 0.02


void RKS(double *v,double corriente,double aleatorio,double D);
void RK(double *v,double corriente);
double golomb_F0(double *v,double I);
double golomb_F1(double *v);
double golomb_F2(double *v);
double golomb_F3(double *v);
double golomb_F4(double *v);
double golomb_minf(double v);
double golomb_cor_na(double *v);
double golomb_cor_k(double *v);
double golomb_cor_d(double *v);
double golomb_cor_l(double v);
double golomb_hinf(double v);
double golomb_tauh(double v);
double golomb_ninf(double v);
double golomb_taun(double v);
double golomb_ainf(double v);
double golomb_taua(double v);
double golomb_binf(double v);
double golomb_taub(double v);

int main(){
	
	FILE *ptrace,*pspike,*info,*frecuencia;
	double iext,sigma,TMAX,IMAX;
	printf("Escriba tiempo maximo de integracion.\n");
	scanf("%lf",&TMAX);
	printf("Escriba el valor maximo de corriente inyectada.\n");
	scanf("%lf",&IMAX);
	printf("Escriba el valor de la amplitud del ruido.\n");
	scanf("%lf",&sigma);

	double *v;
	double vant,t;
	int NSpikes;
	double frec;
	int entero;
	Normaldev_BM aleat=Normaldev_BM(0,sigma,entero);
	string directorio1 = to_string("/home/meli/Documentos/Maestria en Fisica/Simulaciones/Sin corriente lenta/");
	
    v = (double *) malloc (DIM * sizeof(double));
	
    v[0] = V0;
    v[1] = H0;
    v[2] = N0;
    v[3] = A0;
    v[4] = B0; 

   string outputfrec=directorio1+ to_string("FI_gd=")+to_string(golomb_GD)+to_string("_thm=")+to_string(golomb_THETA_M)+to_string(".dat"); 
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
			RK(v,iext+aleat.dev());
			if(v[0]>UMBRAL && vant<UMBRAL && t>TTRANS){
				NSpikes++;
			}
		t=t+dt;
		}
		frec = 1000*NSpikes/(TMAX-TTRANS);
		fprintf(frecuencia, "%lf\t%lf\n", iext,frec); 
	}
	
	free(v);
	fclose(frecuencia);

    return 0;
}


/*.............Algoritmo de integracion...Runge-Kutta estocastico..............*/


void RK(double *v,double corriente){

    double *vaux;
    double *K1,*K2;
    int i;

   vaux = (double*)malloc(DIM*sizeof(double));
   K1 = (double*)malloc(DIM*sizeof(double));
   K2 = (double*)malloc(DIM*sizeof(double));

    for(i=0;i<DIM;i++)
       	vaux[i]=v[i];
    
    K1[0] = dt*golomb_F0(vaux,corriente);
    K1[1] = dt*golomb_F1(vaux);
    K1[2] = dt*golomb_F2(vaux);
    K1[3] = dt*golomb_F3(vaux);
    K1[4] = dt*golomb_F4(vaux);

   
    for(i=0;i<DIM;i++)
       	vaux[i]=v[i]+K1[i];
    
    K2[0] = dt*golomb_F0(vaux,corriente);
    K2[1] = dt*golomb_F1(vaux);
    K2[2] = dt*golomb_F2(vaux);
    K2[3] = dt*golomb_F3(vaux);
    K2[4] = dt*golomb_F4(vaux);
    
    for(i=0;i<DIM;i++)
       	v[i]=v[i]+0.5*(K1[i]+K2[i]);

    free(vaux);
    free(K1);
    free(K2);
}

void RKS(double *v,double corriente,double aleatorio,double D){

    double *vaux;
    double *K1,*K2;
    int i;

   	vaux = (double*)malloc(DIM*sizeof(double));
   	K1 = (double*)malloc(DIM*sizeof(double));
    K2 = (double*)malloc(DIM*sizeof(double));

    for(i=0;i<DIM;i++){
       	vaux[i]=v[i];
    }

    K1[0] = dt*golomb_F0(vaux,corriente);
    K1[1] = dt*golomb_F1(vaux);
    K1[2] = dt*golomb_F2(vaux);
    K1[3] = dt*golomb_F3(vaux);
    K1[4] = dt*golomb_F4(vaux);

    vaux[0] = v[0]+K1[0]+sqrt(2*D*dt)*aleatorio;
    vaux[1] = v[1]+K1[1];
    vaux[2] = v[2]+K1[2];
   	vaux[3] = v[3]+K1[3];
   	vaux[4] = v[4]+K1[4];

    K2[0] = dt*golomb_F0(vaux,corriente);
    K2[1] = dt*golomb_F1(vaux);
    K2[2] = dt*golomb_F2(vaux);
    K2[3] = dt*golomb_F3(vaux);
    K2[4] = dt*golomb_F4(vaux);
    
	v[0]= v[0] + 0.5*(K1[0]+K2[0])+sqrt(2*D*dt)*aleatorio;
	v[1]= v[1] + 0.5*(K1[1]+K2[1]);
	v[2]= v[2] + 0.5*(K1[2]+K2[2]);
	v[3]= v[3] + 0.5*(K1[3]+K2[3]);
	v[4]= v[4] + 0.5*(K1[4]+K2[4]);

    free(vaux);
	free(K1);
	free(K2);
}


double golomb_F0(double *v,double I){
	double x;
	x = -(golomb_cor_na(v)+golomb_cor_k(v)+golomb_cor_d(v)+golomb_cor_l(v[0])-I)/golomb_C;
	return x;
}

double golomb_F1(double *v){
	double x;
	x = (golomb_hinf(v[0])-v[1])/golomb_tauh(v[0]);
	return x;
}

double golomb_F2(double *v){
	double x;
	x = (golomb_ninf(v[0])-v[2])/golomb_taun(v[0]);
	return x;
}

double golomb_F3(double *v){
	double x;
	x = (golomb_ainf(v[0])-v[3])/golomb_taua(v[0]);
	return x;
}

double golomb_F4(double *v){
	double x;
	x = (golomb_binf(v[0])-v[4])/golomb_taub(v[0]);
	return x;
}

double golomb_minf(double v){
	double x;
	x = 1./(1. + exp(-(v-golomb_THETA_M)/golomb_SIGMA_M));
	return x;
}

double golomb_cor_na(double *v){
	double x;
	x = golomb_GNA*golomb_minf(v[0])*golomb_minf(v[0])*golomb_minf(v[0])*v[1]*(v[0]-golomb_VNA);
	return x;
}

double golomb_cor_k(double *v){
	double x;
	x = golomb_GK*v[2]*v[2]*(v[0]-golomb_VK);
	return x;
}

double golomb_cor_d(double *v){
	double x;
	x = golomb_GD*v[3]*v[3]*v[3]*v[4]*(v[0]-golomb_VK);
	return x;
}

double golomb_cor_l(double v){
	double x;
	x = golomb_GL*(v-golomb_VL);
	return x;
}

double golomb_hinf(double v){
	 return 1./(1+exp(-(v-golomb_THETA_H)/golomb_SIGMA_H));
}

double golomb_tauh(double v){
	return 0.5+14/(1+exp(-(v-golomb_THETA_TH)/golomb_SIGMA_TH));
}

double golomb_ninf(double v){
	return 1./(1+exp(-(v-golomb_THETA_N)/golomb_SIGMA_N));
}

double golomb_taun(double v){
	double a,b,x;
	a = 0.087 + 11.4 / (1+exp((v+14.6)/8.6));
	b = 0.087 + 11.4 / (1+exp(-(v-13)/18.7));
	x = a*b;
	return x;
}
double golomb_ainf(double v){
	return 1./(1+exp(-(v-golomb_THETA_A)/golomb_SIGMA_A));
}

double golomb_taua(double v){
	return golomb_TAU_A;
}

double golomb_binf(double v){
	return 1./(1+exp(-(v-golomb_THETA_B)/golomb_SIGMA_B));
}

double golomb_taub(double v){
	return golomb_TAU_B;
}

