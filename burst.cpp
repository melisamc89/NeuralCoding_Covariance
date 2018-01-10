#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include <fstream>
#include <iostream>
#include "iomatrix.h"
#include "nrutilc.h"
#include "LUdcmp.h"
#include "jacobi.h"

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

void StimulusCorrelationMatrix(const VecDoub &vector, MatDoub &CorMat);
int SpikeVector(const VecDoub& stimulus,const VecInt &spike, VecDoub& VecSpike,int n);
int CantSpike(VecInt &spike,int n);
void SpikeTriggeredStatistics(const VecDoub& VecSpike, VecDoub& STA, VecDoub &SigmaSTA, MatDoub& CovarianceMatrix);
void SpikeTriggeredAverage(const VecDoub& VecSpike,VecDoub &STA,VecDoub &SigmaSTA);
void SpikeTriggeredCovariance(const VecDoub& VecSpike, MatDoub &CovarianceMatrix, const VecDoub& STA);

int BurstVector(const VecDoub& tiempo,const VecDoub& stimulus,const VecInt &spike,double Umbral,int cantidad,VecDoub& VecSpike);

void MultiplicaMatriz(const MatDoub &primera,const MatDoub &segunda,MatDoub &producto);
void ImprimeMat(FILE *archivo,const MatDoub &matrix,double dt);

int MaxBurst(VecInt spike,VecDoub t,double umbral);
void CantBurst(const VecInt &spike,const VecDoub &t,double umbral,VecInt &burst);

#define PRIOR 50000
#define NVECTOR 400

int main(int argc, char *argv[]) {
	
	double iext,sigma;
	printf("Escriba el valor de corriente inyectada.\n");
	scanf("%lf",&iext);
	printf("Escriba el valor de la amplitud del ruido\n");
	scanf("%lf",&sigma);
	
	double umbral=50;
	
	//directorios que se ultilizaran en el programa
	string directorio1 = to_string("/home/meli/Documentos/Maestria en Fisica/Simulaciones/Sin corriente lenta/");
	string directorio2 = to_string("/home/meli/Documentos/Maestria en Fisica/ISI/");
	string directorio3 = to_string("/home/meli/Documentos/Maestria en Fisica/Procesamiento/Sin corriente lenta/");
	
	//-----------------------------------------------------------------------------------------------------------------
	//lectura de la inforacion que se tiene del archivo que se desea analizar
	//-----------------------------------------------------------------------------------------------------------------
	FILE *inf;		
	string input2=directorio1+to_string("inf_") + to_string(iext) + to_string("_")+to_string(sigma)+to_string(".dat");
		if((inf = fopen(input2.c_str(), "r")) == NULL){
			printf("No puedo abrir el archivo %s.\n", input2.c_str());
			exit(1);
		}
	double TMAX,frec,dt;
	int NSpikes;
	fscanf(inf,"%lf\t%lf\t%lf\t%i\t%lf\t%lf\n",&iext,&sigma,&TMAX,&NSpikes,&frec,&dt);
	int lines = TMAX/dt;
	cout << "Cantidad de lineas para el analisis : " << lines << endl;
	fclose(inf);
	//-------------------------------------------------------------------------------------------------------------------
	//lectura de archivo de datos
	//-------------------------------------------------------------------------------------------------------------------
	FILE *datos;
	string input=directorio1+to_string("stim_") + to_string(iext) + to_string("_")+to_string(sigma)+to_string(".dat");
			if((datos = fopen(input.c_str(), "r")) == NULL){
				printf("No puedo abrir el archivo %s.\n", input.c_str());
				exit(1);
			}
	VecDoub stimulus(lines),t(lines);
	double v=0;
	VecInt spike(lines);
	int s=0;
	for(int i=0;i<lines;i++){
		fscanf(datos,"%lf\t%lf\t%lf\t%i\n",&t[i],&stimulus[i],&v,&spike[i]);
		if(spike[i]==1)
			s++;
		}
	cout << "Cantidad de Spikes: " << s << endl;
	fclose(datos);
	
	//------------------------------------------------------------------------------------------------------------------
	//Construccion de matriz de correlacion del estimulo
	//------------------------------------------------------------------------------------------------------------------
	VecDoub auxstimulus(PRIOR);
	for(int i=0;i<PRIOR;i++)
		auxstimulus[i]=stimulus[i];
	MatDoub cprior(NVECTOR,NVECTOR);
	StimulusCorrelationMatrix(auxstimulus,cprior);
	//------------------------------------------------------------------------------------------------------------------
	//Caracterizacion de Burst
	//------------------------------------------------------------------------------------------------------------------
	
	int maxspike=MaxBurst(spike,t,umbral);//busco cual es la maxima cantidad de burst en el archivo de datos
	VecInt cantburst(maxspike+1);
	CantBurst(spike,t,umbral,cantburst);//calculo la cantidad de burst que hay con n spikes (considerando hasta el maximo)
	FILE *archivoburst;
	string outputdata=directorio3+to_string("data_") + to_string(iext) + to_string("_")+to_string(sigma)+to_string(".dat");
			if((archivoburst= fopen(outputdata.c_str(), "w")) == NULL){
				printf("No puedo abrir el archivo %s.\n", outputdata.c_str());
				exit(1);
			}
	fprintf(archivoburst,"El numero de spikes es : %d\n", NSpikes);
	fprintf(archivoburst,"El burst de mayor tamaño tiene: %d Spikes\n", maxspike);
	for(int indice=0;indice<maxspike+1;indice++){
		fprintf(archivoburst,"Con %d spikes hay %d burst\n",indice,cantburst[indice]);
		printf("Con %d spikes hay %d burst\n",indice,cantburst[indice]);
		}
	fprintf(archivoburst,"Umbral utilizado %lf\n",umbral);
	fclose(archivoburst);
	
	//--------------------------------------------------------------------------------------------------------------------
	//Analisis de covarianza para cada conjunto de Burst con n Spikes
	//--------------------------------------------------------------------------------------------------------------------
	
	for(int indice=0;indice<=maxspike;indice++){ //solo me importan los disparos que van en conjuntos de al menos dos
		cout << indice << " " << cantburst[indice] << endl;
		if(cantburst[indice]!=0){ //si no hay burst, no puedo realizar el analisis estadistico
		
			VecDoub vecspike(cantburst[indice]*NVECTOR);
			for(int i=0;i<cantburst[indice]*NVECTOR;i++)
			    vecspike[i]=0;
			int NBurstprueba = BurstVector(t,stimulus,spike,umbral,indice,vecspike);//genero el vector con los estimulos que preceden al burst de indice spikes
			
			/*FILE *stim;
			string archstim = directorio3+"g_estimulo_"+to_string(iext)+to_string("_")+to_string(sigma)+to_string("_")+to_string(indice)+to_string(".dat");
			if((stim = fopen(archstim.c_str(), "w")) == NULL){
				printf("No puedo abrir el archivo STA.\n");
				exit(1);
				}
			for(int i=0;i<cantburst[indice]*NVECTOR;i++)
			  fprintf(stim,"%d \t %lf \t %d\n", i,vecspike[i],i/NVECTOR);
			fclose(stim);*/	
			
			//cout << "Calculando STA y matriz de covarianza" << endl;	
			VecDoub STA(NVECTOR),SigmaSTA(NVECTOR);
			MatDoub matcovariance(NVECTOR,NVECTOR);
			SpikeTriggeredStatistics(vecspike,STA,SigmaSTA,matcovariance);		
			FILE *STAverage;
			string archivo1 = directorio3+to_string("STA_")+to_string(iext)+to_string("_")+to_string(sigma)+to_string("_")+to_string(indice)+to_string(".dat");
			if((STAverage = fopen(archivo1.c_str(), "w")) == NULL){
				printf("No puedo abrir el archivo STA.\n");
				exit(1);
				}
			for(int i=0;i<NVECTOR;i++)
				fprintf(STAverage,"%lf\t%lf\t%lf\n",(NVECTOR/2-i)*dt,STA[i],SigmaSTA[i]);
			fclose(STAverage);	
			
			//cout << "Descomposicion LU" << endl;
			MatDoub diagonalizar(NVECTOR,NVECTOR);
			//LUdcmp auxiliar=LUdcmp(cprior);
			//auxiliar.solve(matcovariance,diagonalizar);
			
			/*FILE *arch;
			string archivon = directorio3+"g_mat_"+to_string(iext)+to_string("_")+to_string(sigma)+to_string("_")+to_string(indice)+to_string("_STA.dat");
			if((arch = fopen(archivon.c_str(), "w")) == NULL){
				printf("No puedo abrir el archivo STA.\n");
				exit(1);
				}
			for(int i=0;i<NVECTOR;i++)
			  for(int j=0;j<NVECTOR;j++)
				fprintf(arch,"%d\t%d\t%lf\t%lf\t%lf\n",i,j,cprior[i][j],matcovariance[i][j],diagonalizar[i][j]);
			fclose(arch);*/	
			
			//-------------------------------------------CALCULO AUTOVALORES
			
			//cout << "Jacobi" << endl;
			//Jacobi jac(diagonalizar);
			Jacobi jac(matcovariance);
			eigsrt(jac.d,&jac.v);
			
			//cout << "Escritura de autovalores" << endl;
			FILE *eigenvalues;
			string archivo2 = directorio3+"eigenval"+to_string(iext)+to_string("_")+to_string(sigma)+to_string("_")+to_string(indice)+to_string(".dat");
			if((eigenvalues = fopen(archivo2.c_str(), "w")) == NULL){
				printf("No puedo abrir el archivo eigenval.\n");
				exit(1);
				}
			for(int i=0;i<NVECTOR;i++){
				fprintf(eigenvalues,"%d\t%lf\n",i,jac.d[i]);
				}
			fclose(eigenvalues);
			
			//cout << "Imprime autovectores" << endl;
			FILE *eigenvectors;
			string archivo3 = directorio3+"eigenvec_"+to_string(iext)+to_string("_")+to_string(sigma)+to_string("_")+to_string(indice)+to_string(".dat");
			if((eigenvectors = fopen(archivo3.c_str(), "w")) == NULL){
				printf("No puedo abrir el archivo eigenvec.\n");
				exit(1);
				}
			ImprimeMat(eigenvectors,jac.v,dt);
			fclose(eigenvectors);
			
			}
		}

	return 0;
	}

int MaxBurst(VecInt spike,VecDoub t,double umbral){
	
	double tant=t[0];
	int contador=1;
	int max_spike_x_burst=0;
	for(int i=0;i<spike.size();i++){
		if(spike[i]==1){		
			if(fabs(t[i]-tant)<=umbral){
				contador++;
				}
			else{
				if(contador>max_spike_x_burst)
					max_spike_x_burst=contador;
				contador=1;
				}
			tant=t[i];		
			}
		}
	return max_spike_x_burst;	
	}
	
void CantBurst(const VecInt &spike,const VecDoub &t,double umbral,VecInt &burst){
	
	for(int i=0;i<burst.size();i++)
		burst[i]=0;
		
	double tant=t[0];
	int contador=1;
	for(int i=0;i<spike.size();i++){
		if(spike[i]==1){
			if(fabs(t[i]-tant)<=umbral)
				contador++;
			else{
				burst[contador]+=1;
				contador=1;
				}
			tant=t[i];	
			}
		}
	}


void StimulusCorrelationMatrix(const VecDoub &vector, MatDoub &CorMat){
	
	int lineas = vector.size();
	int n = CorMat.nrows();
	double aux1,aux2;
	int m=lineas/n;
	
	for(int i=0;i<n;i++){
		for(int j=0;j<=i;j++){
			aux1=0;
			aux2=0;
			for(int k=0;k<lineas-n;k++){
				aux1 += vector[k+i]*vector[k+j];
				aux2 += vector[k];
			}
			CorMat[i][j]=aux1/(lineas-n)-(aux2*aux2)/((lineas-n)*(lineas-n));
			}
		}
	for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
			CorMat[i][j]=CorMat[j][i];
			}
		}	
	}

int BurstVector(const VecDoub& tiempo,const VecDoub& stimulus,const VecInt &spike,double umbral,int cantidad,VecDoub& VecSpike){
	
	int j=0;
	double tant=tiempo[0];
	int contador=1;
	int count=0;
	int auxiliar=0;
	  
	for (int k=NVECTOR;k<spike.size();k++){
		auxiliar++;
		if(spike[k]==20){
			if(fabs(tiempo[k]-tant)<=umbral)
				contador++;
			else{
				if(contador==cantidad){
					for(int i=0;i<NVECTOR;i++)
						VecSpike[j*NVECTOR+i]=stimulus[k-auxiliar-i+NVECTOR/2];
					j++;
					}
				contador=1;
				auxiliar=0;
				}
			tant=tiempo[k];
			}
		}
	return j;
	}
	
int SpikeVector(const VecDoub& stimulus,const VecInt &spike, VecDoub& VecSpike, int n){
	
	int lineas = stimulus.size();

	int j=0;
	for (int k=n;k<lineas;k++){
		if(spike[k]==1){
			for(int i=0;i<n;i++)
				VecSpike[j*n+i]=stimulus[k-i+n/2+1];
			j++;
			}
		}
	return j;
	}
	
int CantSpike(VecInt &spike,int n){
	int cant=0;
	for(int i=n;i<spike.size();i++)
		if(spike[i]==1) cant++;
	return cant;
	}
	
void SpikeTriggeredStatistics(const VecDoub& VecSpike, VecDoub& STA, VecDoub &SigmaSTA, MatDoub& CovarianceMatrix){
	
	int n=CovarianceMatrix.nrows();
	
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++)
			CovarianceMatrix[i][j]=0;
		STA[i]=0;
		SigmaSTA[i]=0;
	}
	
	SpikeTriggeredAverage(VecSpike,STA,SigmaSTA);	
	SpikeTriggeredCovariance(VecSpike,CovarianceMatrix,STA);
	
	}
	
void SpikeTriggeredAverage(const VecDoub& VecSpike, VecDoub& STA, VecDoub &SigmaSTA){
	
	int n = STA.size();
	int m = VecSpike.size()/n;

	for(int i=0;i<n;i++){
		STA[i]=0;
		SigmaSTA[i]=0;
		}
	
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			STA[i]+=VecSpike[j*n+i];
			SigmaSTA[i]+=VecSpike[j*n+i]*VecSpike[j*n+i];
		}
		STA[i]=STA[i]/m;
		SigmaSTA[i]=sqrt(SigmaSTA[i]/m-STA[i]*STA[i])/sqrt(m-1);
		}
	}

void SpikeTriggeredCovariance(const VecDoub &VecSpike, MatDoub &CovarianceMatrix, const VecDoub &STA){
	
	double suma;
	int n=CovarianceMatrix.nrows();
	int m = VecSpike.size()/n;		
	for(int i=0;i<n;i++){
		for(int j=0;j<=i;j++){
			suma=0;
			for(int k=0;k<m;k++){
				suma+=VecSpike[k*n+i]*VecSpike[k*n+j];
				}
			CovarianceMatrix[i][j]=suma/m-STA[i]*STA[j];
			}
		}
		
	for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
			CovarianceMatrix[i][j]=CovarianceMatrix[j][i];
			}
		}	
	}

void MultiplicaMatriz(const MatDoub &primera,const MatDoub &segunda,MatDoub &producto){

	double elemento;
	if(primera.ncols()!=segunda.nrows())
		cout << "No se puede multiplicar matrices" << endl;
	
	int p=primera.ncols();
	int n=producto.nrows();
	int m=producto.ncols();

	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			elemento=0;
			for(int k=0;k<p;k++){
				elemento += primera[i][k]*segunda[k][j];			
			}
		producto[i][j]=elemento;
		}
	}

}

void ImprimeMat(FILE *archivo,const MatDoub &matrix,double dt){
	
	int n=matrix.nrows();
	int m=matrix.ncols();
	
	for(int i=0;i<n;i++){
		fprintf(archivo,"%lf\t",(NVECTOR/2-i)*dt);
		for(int j=0;j<m;j++){
			//fprintf(archivo,"%i\t%i\t%lf\n",i,j,matrix[i][j]);
			fprintf(archivo,"%lf\t",matrix[i][j]);
			}
			fprintf(archivo,"\n");
		}
	}
