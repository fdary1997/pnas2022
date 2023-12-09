#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef TYPES_DECLARATION

typedef unsigned int  Int_Type;
typedef unsigned long Long_Type;

typedef struct {
	double x;
	double y;
	double z;
} Vector;

typedef struct {
	double psi;
	double tau;
} Ribbon;

typedef struct {

  Int_Type ar;

	Vector *tt;
	Vector *bb;
	Vector *nn;

	double *f_tt;
	double *f_bb;
	
	double zeta;
	double link; 
	double twist;
	double writhe;
	double energy;	
	double torsion;
	double torque; 
	double curvature;
	double end_to_end;
	
} Statistics;

typedef struct {

  Int_Type d;
  Int_Type h;
  Int_Type m;
  Int_Type s;

} CPU_Time;


#endif
#define TYPES_DECLARATION
#define PI 3.141592653589793238462643383279502884
#define SQRT3 1.732050807568877293527446341505872366
