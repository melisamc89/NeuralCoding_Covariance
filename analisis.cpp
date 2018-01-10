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
int SpikeMatrix(const VecDoub& stimulus,const VecInt &spike, MatDoub& MatSpike);
int SpikeVector(const VecDoub& stimulus,const VecInt &spike, VecDoub& VecSpike,int n);
int CantSpike(VecInt &spike,int n);
void SpikeTriggeredStatistics(const VecDoub& VecSpike, VecDoub& STA, VecDoub &SigmaSTA, MatDoub& CovarianceMatrix);
void SpikeTriggeredAverage(const VecDoub& VecSpike,VecDoub &STA,VecDoub &SigmaSTA,int n);
void SpikeTriggeredCovariance(const VecDoub& VecSpike, MatDoub &CovarianceMatrix, const VecDoub& STA);

void MultiplicaMatriz(const MatDoub &primera,const MatDoub &segunda,MatDoub &producto);
void ImprimeMat(FILE *archivo,const MatDoub &matrix,double dt);

#define NVECTOR 200
#define PRIOR 3000
//#define IMAX 7
		
int main(int argc, char *argv[]) {

	FILE *datos,*info;
	double iext,D,TMAX,frec,dt;
	int NSpike;
	char *model;
	/*printf("Escriba el nombre del modelo que desea analisis: ml1, ml2, hh, golomb\n");
	scanf("%c",model);*/
	printf("Escriba el valor de corriente inyectada.\n");
	scanf("%lf",&iext);
	printf("Escriba el valor de la amplitud del ruido\n");
	scanf("%lf",&D);
	
	model = argv[1];
	/*iext=atof(argv[2]);
	D=atoi(argv[3]);*/

	cout << model << "\t" << iext << "\t" << D << endl;

	string directorio1 = to_string("/home/meli/Maestria/Modelos/Simulaciones/");
	string directorio2 = to_string("/home/meli/Maestria/Correlacion Inversa/Resultados/");
	string modelo = to_string(model)+to_string("_");
	
	//while(iext<IMAX){
		
		//-----------------------------LECTURA DE ARCHIVO DE INFORMACION
		//cout << iext << endl;
	string archivo1= directorio1+modelo+to_string(iext)+to_string("_")+to_string(D)+to_string("_info.dat");
		if((info = fopen(archivo1.c_str(), "r")) == NULL){
			printf("No puedo abrir el archivo %s.\n", archivo1.c_str());
			exit(1);
			}
		fscanf(info,"%lf\t%lf\t%lf\t%i\t%lf\t%lf\n",&iext,&D,&TMAX,&NSpike,&frec,&dt);
		fclose(info);

		cout << "El numero de spikes de lectura es:" << NSpike << endl;
		
		//-----------------------------------LECTURA DE ARCHIVO DE DATOS

		//cout << "Lectura de archivo" << endl;
		string archivo2=directorio1+modelo+to_string(iext)+to_string("_")+to_string(D)+to_string(".dat");
		if((datos = fopen(archivo2.c_str(), "r")) == NULL){
			printf("No puedo abrir el archivo %s.\n", archivo2.c_str());
			exit(1);
			}	
		int NLINE = TMAX/dt;
		double v,t;
		VecDoub stimulus(NLINE);
		VecInt spike(NLINE),burst(NLINE);
		for(int i=0;i<NLINE;i++){
			fscanf(datos,"%lf\t%lf\t%lf\t%i\t%i\n",&t,&stimulus[i],&v,&spike[i],&burst[i]);
		}
		fclose(datos);
		//------------------------------------CONSTRUCCION MATRIZ CPRIOR
		
		cout << "matrix Cprior" << endl;
		VecDoub auxstimulus(PRIOR);
		for(int i=0;i<PRIOR;i++)
			auxstimulus[i]=stimulus[i];
		MatDoub cprior(NVECTOR,NVECTOR);
		StimulusCorrelationMatrix(auxstimulus,cprior);
		
		//-------------------------------------CONSTRUCCION MATRIZ SPIKE
		
		cout << "MatSpike" << endl;
		NSpike = CantSpike(spike,NVECTOR);
		VecDoub vecspike(NSpike*NVECTOR);
		int Nspikeprueba = SpikeVector(stimulus,spike,vecspike,NVECTOR);
		cout << "El numero de spikes en la matriz de spikes es : " << Nspikeprueba << endl; 
		if(Nspikeprueba!=NSpike){
			printf("Error en conteo de spikes\n");
			return -1;
			}

		//-------------------------CALCULO DE STA Y MATRIZ DE COVARIANZA
		
		cout << "Calculando STA y matriz de covarianza" << endl;	
		VecDoub STA(NVECTOR),SigmaSTA(NVECTOR);
		MatDoub matcovariance(NVECTOR,NVECTOR);
		SpikeTriggeredStatistics(vecspike,STA,SigmaSTA,matcovariance);		
		FILE *STAverage;
		string archivo3 = directorio2+modelo+to_string(iext)+to_string("_")+to_string(D)+to_string("_STA.dat");
		if((STAverage = fopen(archivo3.c_str(), "w")) == NULL){
			printf("No puedo abrir el archivo STA.\n");
			exit(1);
			}
		for(int i=0;i<NVECTOR;i++)
			fprintf(STAverage,"%lf\t%lf\t%lf\n",(i-NVECTOR/2+1)*dt,STA[i],SigmaSTA[i]);
		fclose(STAverage);
		
		//---------------------------------------------DESCOMPOSICION LU
		
		cout << "Descomposicion LU" << endl;
		MatDoub diagonalizar(NVECTOR,NVECTOR);
		/*LUdcmp auxiliar=LUdcmp(cprior);
		auxiliar.solve(matcovariance,diagonalizar);*/
		
		//-------------------------------------------CALCULO AUTOVALORES
		
		cout << "Jacobi" << endl;
		//Jacobi jac(diagonalizar);
		Jacobi jac(matcovariance);
		eigsrt(jac.d,&jac.v);
		
		//--------------------------ESCRITURA AUTOVALORES Y AUTOVECTORES
		
		cout << "Escritura de autovalores" << endl;
		FILE *eigenvalues;
		string archivo4=directorio2+modelo+to_string(iext)+to_string("_")+to_string(D)+to_string("_eigenval.dat");
		if((eigenvalues = fopen(archivo4.c_str(), "w")) == NULL){
			printf("No puedo abrir el archivo eigenval.\n");
			exit(1);
			}
		for(int i=0;i<NVECTOR;i++){
			fprintf(eigenvalues,"%lf\t%lf\n",(i-NVECTOR/2+1)*dt,jac.d[i]);
			}
		fclose(eigenvalues);
		
		cout << "Imprime autovectores" << endl;
		FILE *eigenvectors;
		string archivo5=directorio2+modelo+ to_string(iext)+to_string("_")+to_string(D)+to_string("_eigenvec.dat");
		if((eigenvectors = fopen(archivo5.c_str(), "w")) == NULL){
			printf("No puedo abrir el archivo eigenvec.\n");
			exit(1);
			}
		ImprimeMat(eigenvectors,jac.v,dt);
		fclose(eigenvectors);
	//}
	
		cout << "imprimiendo matrices" << endl;
		FILE *matrices;
		string archivo6 = directorio2+modelo+to_string(iext)+to_string("_")+to_string(D)+to_string("_matrices.dat");
		if((matrices = fopen(archivo6.c_str(), "w")) == NULL){
			printf("No puedo abrir el archivo matrices.\n");
			exit(1);
			}
		for(int i = 0; i< NVECTOR; i++){
			for(int j=0;j<NVECTOR;j++)
				fprintf(matrices,"%i\t%i\t%lf\t%lf\t%lf\n", i,j,cprior[i][j],matcovariance[i][j],diagonalizar[i][j]);
			fprintf(matrices,"\n");}
		fclose(matrices);
		
return 0;
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
	
	SpikeTriggeredAverage(VecSpike,STA,SigmaSTA,n);	
	SpikeTriggeredCovariance(VecSpike,CovarianceMatrix,STA);
	
	}
	
void SpikeTriggeredAverage(const VecDoub& VecSpike, VecDoub& STA, VecDoub &SigmaSTA, int n){
	
	int m = VecSpike.size();
	int NSpike = m/n;
	
	for(int i=0;i<n;i++){
		STA[i]=0;
		SigmaSTA[i]=0;
		}
		
	for(int i=0;i<n;i++){
		for(int j=0;j<NSpike;j++){
			STA[i]+=VecSpike[j*n+i];
			SigmaSTA[i]+=VecSpike[j*n+i]*VecSpike[j*n+i];
		}
		STA[i]=STA[i]/NSpike;
		SigmaSTA[i]=sqrt(SigmaSTA[i]/NSpike-STA[i]*STA[i])/sqrt(NSpike-1);
		}
	}

void SpikeTriggeredCovariance(const VecDoub &VecSpike, MatDoub &CovarianceMatrix, const VecDoub &STA){
	
	double suma;
	int n=CovarianceMatrix.nrows();
	int l = VecSpike.size()/10;
	//cout << " l = " << l << endl;
	int NSpike = VecSpike.size()/n;
			
	for(int i=0;i<n;i++){
		for(int j=0;j<=i;j++){
			suma=0;
			for(int k=0;k<NSpike-1;k++){
				suma+=VecSpike[k+i]*VecSpike[k+j];
				}
			CovarianceMatrix[i][j]=suma/NSpike-STA[i]*STA[j];
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
		fprintf(archivo,"%lf\t",(i-NVECTOR/2+1)*dt);
		for(int j=0;j<m;j++){
			//fprintf(archivo,"%i\t%i\t%lf\n",i,j,matrix[i][j]);
			fprintf(archivo,"%lf\t",matrix[i][j]);
			}
			fprintf(archivo,"\n");
		}
	}
