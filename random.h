#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "nrutilc.h"

typedef unsigned long long int Ullong;
typedef unsigned int Uint;

struct Ran {
	Ullong u,v,w;
	Ran(Ullong j) : v(4101842887655102017LL), w(1) {
	u = j ^ v; int64();
	v = u; int64();
	w = v; int64();
}
	inline Ullong int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
	inline Uint int32() { return (Uint)int64(); }
};

struct Normaldev_BM : Ran{
	
	double mu,sig;
	double storedval;
	Normaldev_BM(double mmu, double ssig, Ullong i)
	: Ran(i), mu(mmu), sig(ssig), storedval(0.){};
	double dev(){
		double v1,v2,rsq,fac;
		if (storedval == 0.){
			do {
			v1=2.0*doub()-1.0;
			v2=2.0*doub()-1.0;
			rsq=v1*v1+v2*v2;
			}while (rsq >= 1.0 || rsq == 0.0);
		
			fac=sqrt(-2.0*log(rsq)/rsq);
			storedval = v1*fac;
			return mu + sig*v2*fac;
		}
		else{
			fac = storedval;
			storedval = 0.;
			return mu + sig*fac;
		}
	}
};
