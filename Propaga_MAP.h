/*********************************************************************************
Copyright 2010, 2011, 2012 Stefano Sinigardi, Graziano Servizi, Giorgio Turchetti
*********************************************************************************/


#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
#pragma warning(disable : 593)
#pragma warning(disable : 869)
#pragma warning(disable : 981)
#endif

#define DEBUG

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cstring>
#define C						2.99792458e+10		// cm / s
#define MP						1.6726231e-24		// g
#define MP_MEV					938.272013			// proton mass [MeV/c^2]
#define CHARGE					4.803262e-10		// statC    valore usato da Turchetti; nb: e' impreciso negli ultimi due decimali
#define FROM_TESLA_TO_GAUSS		1.0e+4
#define FROM_VOLT_TO_STATVOLT	3.335640951982e-3	// 1 statvolt = 299.792458 volts.
#define DA_ERG_A_MEV			6.241509744512e+5	// conversione mia come sotto descritta
#define FROM_M_TO_CM			1.0e+2

#define gamma_rel_inv(x)		(1.0 / sqrt(1.0 + x[3]*x[3]+x[4]*x[4]+x[5]*x[5]))
#define gamma_rel(x)			(sqrt(1.0 + x[3]*x[3]+x[4]*x[4]+x[5]*x[5]))		// gamma relativistico definito in funzione dei gamma*beta usati nei file

#define NUMERO_PARAMETRI_MINIMO_SIMULAZIONE  8

#define N_TIPI_MAGNETICI			 4
#define N_PARAMETRI_LATTICINO		 4
#define N_DIMENSIONI_SPAZIO_FASI	 6

#define DRIFT						"O"
#define SOLENOID					"S"
#define FOCUSING					"F"
#define DEFOCUSING					"D"
#define _DRIFT_						 0
#define _DRIFT_CHECK_				 0.5
#define _SOLENOID_					 1
#define _SOLENOID_CHECK_			 1.5
#define _FOCUSING_					 2
#define _FOCUSING_CHECK_			 2.5
#define _DEFOCUSING_				 3
#define _DEFOCUSING_CHECK_			 3.5

class Elemento_magnetico;
class Parametri;
class Propaga_maps;
void create_gnuplot_file();
