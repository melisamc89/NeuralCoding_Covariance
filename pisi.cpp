#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <string>
#include <sstream>
using namespace std;

#define BINE 1

template<class T>
inline string to_string(const T& t){
    stringstream ss;
    ss << t;
    return ss.str();
}

typedef struct {

double valores[2];
int posiciones[2];

}ExtDouble;

typedef struct {

int valores[2];
int posiciones[2];

}ExtInt;

double *ISIs(double iext,double sigma);
void *InterSpikeInterval(FILE *datos,double *interval,int NSpikes,int lines);
int *Pisi(double iext,double sigma,double *data);
void FindExtDouble(double *vector,ExtDouble &ext);
void FindExtInt(int *vector,ExtInt &ext);
	
int main(){

	double iext=1;
	double sigma=1;
	double *ISpikeI;
	ISpikeI=ISIs(iext,sigma);
	ExtDouble extremos_isi;
	FindExtDouble(ISpikeI,extremos_isi);
	int *histograma=Pisi(iext,sigma,ISpikeI);
	
	return 0;
}

/*double SelectiveThreshold(int *datos){
	
	ExtInt extremos_hist;
	FindExtInt(histograma,extremos_hist);
	int maxpos=extremos_hist.posiciones[0];
	printf("%d\n",maxpos);
	
	int *auxiliar;
	auxiliar=(int*)malloc((histograma[0]-maxpos+1)*sizeof(int));
	auxiliar[0]=histograma[0]-maxpos+1;
	printf("%d\n",auxiliar[0]);
	for(int i=1;i<auxiliar[0];i++)
		auxiliar[i]=histograma[i+maxpos];
	ExtInt extremos_new;
	FindExtInt(auxiliar,extremos_new);
	printf("%lf\n",extremos_new.posiciones[1]*0.5+extremos_isi.valores[1]);
	
	}*/

double *ISIs(double iext,double sigma){
	
	FILE *datos,*dataisi;
	string directorio1 = to_string("/home/meli/Documentos/Maestria en Fisica/Simulaciones/");
	string directorio2 = to_string("/home/meli/Documentos/Maestria en Fisica/ISI/");

	FILE *inf;		
	string input2=directorio1+to_string("g_inf_") + to_string(iext) + to_string("_")+to_string(sigma)+to_string(".dat");
		if((inf = fopen(input2.c_str(), "r")) == NULL){
			printf("No puedo abrir el archivo %s.\n", input2.c_str());
			exit(1);
		}
	double TMAX,frec,dt;
	int NSpikes;
	fscanf(inf,"%lf\t%lf\t%lf\t%i\t%lf\t%lf\n",&iext,&sigma,&TMAX,&NSpikes,&frec,&dt);
	int lines = TMAX/dt;
	fclose(inf);
	
	string input=directorio1+to_string("g_stim_") + to_string(iext) + to_string("_")+to_string(sigma)+to_string(".dat");
			if((datos = fopen(input.c_str(), "r")) == NULL){
				printf("No puedo abrir el archivo %s.\n", input.c_str());
				exit(1);
			}
			
	string output=directorio2+to_string("g_isi_") + to_string(iext) + to_string("_")+to_string(sigma)+to_string(".dat");
			if((dataisi = fopen(output.c_str(), "w")) == NULL){
				printf("No puedo abrir el archivo %s.\n", output.c_str());
				exit(1);
			}

	double *isi;
	isi=(double*)malloc((NSpikes+1)*(sizeof(double)));
	isi[0]=NSpikes;
	InterSpikeInterval(datos,isi,NSpikes,lines);	
	fclose(datos);
	
	for(int i=0;i<NSpikes;i++)
		fprintf(dataisi,"%lf\n",isi[i]);
	fclose(dataisi);
	return isi;
	}
	
void *InterSpikeInterval(FILE *datos,double *interval,int NSpikes,int lines){
	
	int count=1;
	double t,tant=0,stimulus,v;
	int spike;
	for(int i=0;i<lines;i++){
		fscanf(datos,"%lf\t%lf\t%lf\t%i\n",&t,&stimulus,&v,&spike);
		if(spike==1){
			interval[count]=t-tant;
			tant=t;
			count++;
			}
		}
	}
	
int *Pisi(double iext,double sigma,double *data){
	
	FILE *histograma;
	string directorio2 = to_string("/home/meli/Documentos/Maestria en Fisica/ISI/");
	string output=directorio2+to_string("g_pisi_") + to_string(iext) + to_string("_")+to_string(sigma)+to_string(".dat");
			if((histograma = fopen(output.c_str(), "w")) == NULL){
				printf("No puedo abrir el archivo %s.\n", output.c_str());
				exit(1);
			}
	ExtDouble maxmin;
	FindExtDouble(data,maxmin);
	int tam=(maxmin.valores[0]-maxmin.valores[1])*2;
	
	int *distribucion;
	distribucion = (int*)malloc(tam*sizeof(int));
	distribucion[0]=tam;
	for(int i=1;i<tam;i++)
		distribucion[i]=0;
	for(int i=1;i<data[0];i++){
		int a=2*(data[i]-maxmin.valores[1]);
		if(a>0 && a<tam)
			distribucion[a]++;
		}
	for(int i=0;i<tam;i++)
		fprintf(histograma,"%lf\t%d\n", i*0.5+maxmin.valores[1],distribucion[i]);
	fclose(histograma);	
	
	return distribucion;
	}

void FindExtDouble(double *vector,ExtDouble &ext){

	double max=vector[1];
	double min=vector[1];
	int posmax=1;
	int posmin=1;
	int size=vector[0];
	for(int i=1;i<size;i++){
		if(vector[i]>max){
			max=vector[i];
			posmax=i;
			}
		if(vector[i]<min){
			min=vector[i];
			posmin=i;
			}
		}

	ext.valores[0]=max;
	ext.valores[1]=min;
	ext.posiciones[0]=posmax;
	ext.posiciones[1]=posmin;

	}
	
void FindExtInt(int *vector,ExtInt &ext){

	int max=vector[1];
	int min=vector[1];
	int posmax=1;
	int posmin=1;
	int size=vector[0];
	for(int i=1;i<size;i++){
		if(vector[i]>max){
			max=vector[i];
			posmax=i;
			}
		if(vector[i]<min){
			min=vector[i];
			posmin=i;
			}
		}

	ext.valores[0]=max;
	ext.valores[1]=min;
	ext.posiciones[0]=posmax;
	ext.posiciones[1]=posmin;

	}
	
