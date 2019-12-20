#ifndef UTILS_H
#define UTILS_H

typedef struct Points {double x, y, z;} Point;  


// Primitives
extern double calcDot_c(Point p1, Point p2);
extern Point calcAdd_c(Point p1, Point p2);
extern Point calcMult_c(double a,Point p1);
extern double checkBlockage_c(Point ctrPoint1, Point ctrPoint2, Point ctrPointB, 
    Point normalVectB, Point* vertsB, const double areaB);
extern double calcArea3DPoly_c( int n, Point* V, Point N );
extern double calcAngle_c( Point v1, Point v2 );
extern void calcRodriguesMtx_c(const double angle, const double* k_pt, double* R_pt);


// Local usage
extern void initPoint(Point *pt, double x, double y, double z);
// extern void initSimplePlane(SimplePlane *plane, Point normalVect, Point ctrPoint, double L, double W);

// #include "lapacke.h"
// extern void print_matrix( char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda );
// extern void print_int_vector( char* desc, lapack_int n, lapack_int* a );
// extern void testlapacke();

#endif