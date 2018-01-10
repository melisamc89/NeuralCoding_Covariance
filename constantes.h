#ifndef FUNCIONESDIN_H
#define FUNCIONESDIN_H

//Constantes del modelo golomb

#define golomb_CONST 21
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

//Constantes del modelo de Hodgkin-Huxel

#define hh_C 1.0
#define hh_CONST 6
/*
#define hh_GL 0.3
#define hh_VL 10.6
#define hh_GNA 90.0
#define hh_VNA 115.0
#define hh_GK 36.0
#define hh_VK -12.0*/

#define hh_GL 0.3
#define hh_VL -54.4
#define hh_GNA 120.0
#define hh_VNA 50.0
#define hh_GK 36.0
#define hh_VK -77.0

//Constantes del modelo de Morris-Lecar

/*parametros para tipo 2
#define V3 0.0
#define V4 30.0
#define phi 0.2

parametros para tipo 1
#define V3 10.0
#define V4 14.0
#define phi 0.3333*/

#define ml_C 1.0
#define ml_CONST 11

#define ml_GCA 1.1
#define ml_VCA 100.0
#define ml_GK 2.0
#define ml_VK -70.0
#define ml_GL 0.5
#define ml_VL -50.0
#define ml_V1 -1.0
#define ml_V2 15.0
#define ml_V3 10.0
#define ml_V4 14.0
#define ml_phi 0.3333

#endif
