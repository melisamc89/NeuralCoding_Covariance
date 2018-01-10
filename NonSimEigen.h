#ifndef NONSIMEIGEN_H
#define NONSIMEIGEN_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <limits>
#include <complex>
#include "nrutilc.h"

typedef NRVec<int> VecInt, VecInt_O, VecInt_IO;
typedef const NRVec<int>VecInt_I;

typedef NRVec<double> VecDoub, VecDoub_O, VecDoub_IO;
typedef const NRVec<double> VecDoub_I;

typedef NRVec<complex<double> > VecComplex;

typedef NRMat<int> MatInt, MatInt_O, MatInt_IO;
typedef const NRMat<int> MatInt_I;

typedef NRMat<double> MatDoub, MatDoub_O, MatDoub_IO;
typedef const NRMat<double> MatDoub_I;

struct Unsymmeig {
	int n;
	MatDoub a,zz;
	VecComplex wri;
	VecDoub scale;
	VecInt perm;
	bool yesvecs,hessen;
	Unsymmeig (MatDoub_I &aa, bool yesvec=true, bool hessenb=false) :
	n(aa.nrows()), a(aa), zz(n,n,0.0), wri(n), scale(n,1.0), perm(n),
	yesvecs(yesvec), hessen(hessenb) {
		balance();
		if (!hessen) elmhes();
		if (yesvecs) {
			for (int i=0;i<n;i++)
				zz[i][i]=1.0;
			if (!hessen) eltran();
			hqr2();
			balbak();
			sortvecs();
		} else {
			hqr();
			sort();
		}
		}
	void balance();
	void elmhes();
	void eltran();
	void hqr();
	void hqr2();
	void balbak();
	void sort();
	void sortvecs();
};

void Unsymmeig::balance(){
	 
	const double RADIX = numeric_limits<double>::radix; 
	bool done=false; 
	double sqrdx=RADIX*RADIX; 
	while (!done){ 
		done=true; 
		for (int i=0;i<n;i++){ 
			double r=0.0,c=0.0; 
			for (int j=0;j<n;j++) 
			if (j != i){ 
				c += abs(a[j][i]); 
				r += abs(a[i][j]); 
			} 
			if (c != 0.0 && r != 0.0){ 
				double g=r/RADIX; 
				double f=1.0; 
				double s=c+r; 
				while (c<g){ 
				  f *= RADIX; 
				  c *= sqrdx;  
				  } 
				  g=r*RADIX; 
				  while (c>g){ 
					  f /= RADIX; 
					  c /= sqrdx; 
				  } 
				  if ((c+r)/f < 0.95*s) { 
					  done=false; 
					  g=1.0/f; 
					  scale[i] *= f; 
					  for (int j=0;j<n;j++) a[i][j] *= g;
					  for (int j=0;j<n;j++) a[j][i] *= f; 
				  } 
			  } 
		  } 
	  } 
  } 
  
 void Unsymmeig::balbak() { 
	for (int i=0;i<n;i++) 
		for (int j=0;j<n;j++) 
			zz[i][j] *= scale[i]; 
} 

void Unsymmeig::elmhes()  { 
	
	for (int m=1;m<n-1;m++) { 
		double x=0.0; 
		int i=m; 
		for (int j=m;j<n;j++){ 
			if (abs(a[j][m-1]) > abs(x)){ 
				x=a[j][m-1]; 
				i=j; 
			} 
		} 
		perm[m]=i; 
		if (i != m) {  
			for (int j=m-1;j<n;j++) SWAP(a[i][j],a[m][j]); 
			for (int j=0;j<n;j++) SWAP(a[j][i],a[j][m]); 
		} 
		if (x != 0.0){ 
			for (i=m+1;i<n;i++){ 
				double y=a[i][m-1]; 
				if (y != 0.0) { 
					y /= x; 
					a[i][m-1]=y; 
					for (int j=m;j<n;j++) a[i][j] -= y*a[m][j]; 
					for (int j=0;j<n;j++) a[j][m] += y*a[j][i]; 
				} 
			} 
		} 
	} 
} 

void Unsymmeig::eltran(){ 

	for (int mp=n-2;mp>0;mp--) { 
		for (int k=mp+1;k<n;k++) 
			zz[k][mp]=a[k][mp-1]; 
		int i=perm[mp]; 
		if (i != mp) { 
			for (int j=mp;j<n;j++) { 
				zz[mp][j]=zz[i][j]; 
				zz[i][j]=0.0; 
			} 
			zz[i][mp]=1.0; 
		} 
	} 
} 

void eigsrt(VecComplex &d, MatDoub_IO *v=NULL){
	
	int k;
	int n=d.size();
	for (int i=0;i<n-1;i++){
		complex p=d[k=i];
		for (int j=i;j<n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			if (v != NULL)
			for (int j=0;j<n;j++) {
				p=(*v)[j][i];
				(*v)[j][i]=(*v)[j][k];
				(*v)[j][k]=p;
			}
		}
	}
}

void Unsymmeig::sort() {
	if (yesvecs)
		eigsrt(wri,&zz);
	else
		eigsrt(wri);
}

#endif
