#ifndef __BPNET_GEOM__
#define __BPNET_GEOM__
#include <stdio.h>
#include <math.h>
#define PI 3.141592654
double dist3d(double x1, double y1, double z1, double x2, double y2, double z2){
      return sqrt ( (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2) + (z1-z2) * (z1-z2) );
}

double distsqr3d(double x1, double y1, double z1, double x2, double y2, double z2){
      return ( (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2) + (z1-z2) * (z1-z2) );
}

void addpoint3d(double* x1, double* y1, double* z1, 
	    const double x2, const double y2, const double z2){
      *x1 += *x1 + x2;
      *y1 += *y1 + y2;
      *z1 += *z1 + z2;
}

void getplane3d(double* A, double* B, double* C, double* D,
	    double x1, double y1, double z1,
	    double x2, double y2, double z2,
	    double x3, double y3, double z3,
	    ){
      double yz = ((y2-y1)*(z3-z1)) - ((y3-y1)*(z2-z1));
      double xz = ((x2-x1)*(z3-z1)) - ((x3-x1)*(z2-z1));
      double xy = ((x2-x1)*(y3-y1)) - ((x3-x1)*(y2-y1));
      *A = yz;
      *B = -xz;
      *C = xy;
      *D = (-x1)*yz - (-y1)*xz + (-z1)*xy;
}

double plane_perpdist(const double A, const double B, const double C, const double D,
	    const double x1, const double y1, const double z1){
      double denominator = sqrt(A*A + B*B + C*C);
      double dist = (A * x1 + B * y1 + C * z1 + D)/ denominator;
      return fabs(dist);
}

double plane_plane_angledeg(double A1, double B1, double C1, double D1,
	    double A2, double B2, double C2, double D2){
      double deno1 =  sqrt(A1*A1 + B1*B1 + C1*C1);
      double deno2 =  sqrt(A2*A2 + B2*B2 + C2*C2);
      double cos_alpha = (A1 * A2 + B1 * B2 + C1 * C2) / (deno1 * deno2);
      double alpha     =(180.0 * acos(cos_alpha))/PI;
      return alpha;
}

#endif // __BPNET_GEOM__

