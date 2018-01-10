#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

void RKS(double *v,double *constantes,double c,double corriente,double aleatorio,double D){

    double *vaux;
    double *K1,*K2;
    int i;

   	vaux = (double*)malloc(DIM*sizeof(double));
   	K1 = (double*)malloc(DIM*sizeof(double));
    K2 = (double*)malloc(DIM*sizeof(double));

    for(i=0;i<DIM;i++){
       	vaux[i]=v[i];
    }

	K1[0] = dt*F0(vaux,constantes,c,corriente,aleatorio);
	K1[1] = dt*F1(vaux,constantes);
	K1[2] = dt*F2(vaux,constantes);
	K1[3] = dt*F3(vaux,constantes);
	K1[4] = dt*F4(vaux,constantes);

    vaux[0] = v[0]+K1[0]+sqrt(2*D*dt)*aleatorio;
    vaux[1] = v[1]+K1[1];
    vaux[2] = v[2]+K1[2];
   	vaux[3] = v[3]+K1[3];
   	vaux[4] = v[4]+K1[4];

   	K2[0] = dt*F0(vaux,constantes,c,corriente,aleatorio);
	K2[1] = dt*F1(vaux,constantes);
	K2[2] = dt*F2(vaux,constantes);
	K2[3] = dt*F3(vaux,constantes);
	K2[4] = dt*F4(vaux,constantes);
    
	v[0]= v[0] + 0.5*(K1[0]+K2[0])+sqrt(2*D*dt)*aleatorio;
	v[1]= v[1] + 0.5*(K1[1]+K2[1]);
	v[2]= v[2] + 0.5*(K1[2]+K2[2]);
	v[3]= v[3] + 0.5*(K1[3]+K2[3]);
	v[4]= v[4] + 0.5*(K1[4]+K2[4]);

    free(vaux);
	free(K1);
	free(K2);
}

