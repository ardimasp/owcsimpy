
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 

#include "utils.h"

// #include "lapacke.h"

/* Auxiliary routine: printing a matrix */
// void print_matrix( char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda ) 
// {
//         lapack_int i, j;
//         printf( "\n %s\n", desc );
//         for( i = 0; i < m; i++ ) {
//                 for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
//                 printf( "\n" );
//         }
// }
 
// /* Auxiliary routine: printing a vector of integers */
// void print_int_vector( char* desc, lapack_int n, lapack_int* a ) 
// {
//         lapack_int j;
//         printf( "\n %s\n", desc );
//         for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
//         printf( "\n" );
// }

// void testlapacke()
// {
//     // export LD_LIBRARY_PATH=/anaconda3/envs/owcsimpy-dev/lib/:$LD_LIBRARY_PATH
//     int N = 5;
//     int NRHS = 3;
//     int LDA = N;
//     int LDB = NRHS;
//     lapack_int n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
//     lapack_int ipiv[5];
//     double a[5*5] = {
//         6.80, -6.05, -0.45,  8.32, -9.67,
//        -2.11, -3.30,  2.58,  2.71, -5.14,
//         5.66, 5.36, -2.70,  4.35, -7.26,
//         5.97, -4.44,  0.27, -7.17, 6.08,
//         8.23, 1.08,  9.04,  2.14, -6.87
//     };
//     double b[3*5] = {
//         4.02, -1.56, 9.81,
//         6.19,  4.00, -4.09,
//        -8.22, -8.67, -4.57,
//        -7.57,  1.75, -8.61,
//        -3.03,  2.86, 8.99
//     };

//     double aNorm;
//        double rcond;
//        char ONE_NORM = '1';
//        lapack_int NROWS = n;
//        lapack_int NCOLS = n;
//        lapack_int LEADING_DIMENSION_A = n;
 
//               /* Print Entry Matrix */
//         print_matrix( "Entry Matrix A", n, n, a, lda );
//         /* Print Right Rand Side */
//         print_matrix( "Right Rand Side", n, nrhs, b, ldb );
//         printf( "\n" );
//         /* Executable statements */
//         printf( "LAPACKE_dgecon Example Program Results\n" );
//         aNorm = LAPACKE_dlange(LAPACK_ROW_MAJOR, ONE_NORM, NROWS, NCOLS, a, LEADING_DIMENSION_A);
//         info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, NROWS, NCOLS, a, LEADING_DIMENSION_A, ipiv);
//         info = LAPACKE_dgecon(LAPACK_ROW_MAJOR, ONE_NORM, n, a, LEADING_DIMENSION_A, aNorm, &rcond); // aNorm should be 35.019999999999996
//         double work[4*N];
//         int iwork[N];
//         //info = LAPACKE_dgecon_work(LAPACK_ROW_MAJOR, ONE_NORM, n, a, LEADING_DIMENSION_A, aNorm, &rcond, work, iwork); // aNorm should be 35.019999999999996
//         //dgecon_( &ONE_NORM, &n, a, &LEADING_DIMENSION_A, &aNorm, &rcond, work, iwork, &info );
//         /* Check for the exact singularity */
//               if (info == 0)
//               {
//                      printf("LAPACKE_dgecon completed SUCCESSFULLY...\n");
//               }
//               else if ( info < 0 )
//               {
//             printf( "Element %d of A had an illegal value\n", -info );
//             exit( 1 );
//         }
//               else
//               {
//             printf( "Unrecognized value of INFO = %d\n", info );
//             exit( 1 );
//               }
 
//         /* Print solution */
//        printf("LAPACKE_dlange / One-norm of A = %lf\n", aNorm);
//         printf("LAPACKE_dgecon / RCOND of A    = %f\n", rcond);

// }

double calcDot_c(Point p1, Point p2)
{
    return p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
}

Point calcAdd_c(Point p1, Point p2)
{
    Point res;
    res.x = p1.x+p2.x;
    res.y = p1.y+p2.y;
    res.z = p1.z+p2.z;
    return res;
}

Point calcMult_c(double a,Point p1)
{
    Point res;
    res.x = a*p1.x;
    res.y = a*p1.y;
    res.z = a*p1.z;
    return res;
}

// FIXME: change the type into bool or char
double checkBlockage_c(Point ctrPoint1, Point ctrPoint2, Point ctrPointB, 
    Point normalVectB, Point vertsB[4], const double areaB)
{
    
    const double verysmall = 1e-6;
    int i;

    // FIXME: static, only handle rect, update for future release
    Point verts[5]; 
    Point triVerts[4]; 
    // Point *verts = malloc(5*sizeof(Point)); // rectangle
    // Point *triVerts = malloc(4*sizeof(Point)); // triangle

    for (i=0;i<4;i++)
    {
        verts[i] = vertsB[i];
        // printf("%f,%f,%f \n",vertsB[i].x,vertsB[i].y,vertsB[i].z);
    }

    verts[4] = vertsB[0];

    // for (int i=0;i<5;i++)
    //     printf("%f,%f,%f \n",verts[i].x,verts[i].y,verts[i].z);
    // double area = calcArea3DPoly_c(4,verts,normalVectB);

    Point u,w;
    u = calcAdd_c(ctrPoint2,calcMult_c(-1,ctrPoint1));
    w = calcAdd_c(ctrPoint1,calcMult_c(-1,ctrPointB));
    // u.x = ctrPoint2.x-ctrPoint1.x;
    // u.y = ctrPoint2.y-ctrPoint1.y;
    // u.z = ctrPoint2.z-ctrPoint1.z;
    // w.x = ctrPoint1.x-ctrPointB.x;
    // w.y = ctrPoint1.y-ctrPointB.y;
    // w.z = ctrPoint1.z-ctrPointB.z;

    double dot = calcDot_c(normalVectB,u);

    // Sum of triangles' areas
    double totalarea = 0;
    // double area = 0; // Debug-purpose

    if (fabs(dot) < verysmall) // Line segment and the planeB is parallel
    {
        // printf("dot:%f\n",dot);
        // free(verts);
        // free(triVerts);
        return 0; // Not blocked
    }
    else
    {
        // Parametric value of the intersecting point
        double ti = -1*calcDot_c(normalVectB,w)/calcDot_c(normalVectB,u);

        // if(!(0 <= ti && ti <= 1) )
        if(!(0 < ti && ti < 1) )
        {
            /* 
            This is when the intersecting point is not between plane1 and plane2. 
            In addition, a strict inequality is needed to handle the case when 
            the intersencting point comes from the same plane as the tail or head 
            points. 
            */
            // free(verts);
            // free(triVerts);
            return 0;
        }
        else
        {
            // The intersecting point
            Point intersectingPoint = calcAdd_c(ctrPoint1,calcMult_c(ti,u));

            // Calculate triangles' areas
            for(i=0;i<4;i++) // loop for 4 triangles 
            {
                // Assign triangles' vertices
                triVerts[0] = intersectingPoint;
                triVerts[1] = verts[i];
                triVerts[2] = verts[i+1];
                triVerts[3] = intersectingPoint;
                totalarea += calcArea3DPoly_c(3,triVerts,normalVectB);
                
                // Debug
                // area = calcArea3DPoly_c(3,triVerts,normalVectB);
                // totalarea += area;
                // printf("i:%d,area:%f\n",i,area);
            }

            // printf("totalarea:%f\n",i,totalarea);

            // free(verts);
            // free(triVerts);

            if(fabs(totalarea-areaB) < verysmall)
                return 1; // Blocked
            else
                return 0; // Not blocked
            // return totalarea;
        }

    }

    
}

// see: http://geomalgorithms.com/a01-_area.html
double calcArea3DPoly_c( int n, Point* V, Point N )
{
    
    float area = 0;
    float an, ax, ay, az; // abs value of normal and its coords
    int  coord;           // coord to ignore: 1=x, 2=y, 3=z
    int  i, j, k;         // loop indices

    if (n < 3) return 0;  // a degenerate polygon

    // select largest abs coordinate to ignore for projection
    ax = (N.x>0 ? N.x : -N.x);    // abs x-coord
    ay = (N.y>0 ? N.y : -N.y);    // abs y-coord
    az = (N.z>0 ? N.z : -N.z);    // abs z-coord

    coord = 3;                    // ignore z-coord
    if (ax > ay) {
        if (ax > az) coord = 1;   // ignore x-coord
    }
    else if (ay > az) coord = 2;  // ignore y-coord

    // compute area of the 2D projection
    switch (coord) {
      case 1:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i].y * (V[j].z - V[k].z));
        break;
      case 2:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i].z * (V[j].x - V[k].x));
        break;
      case 3:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i].x * (V[j].y - V[k].y));
        break;
    }
    switch (coord) {    // wrap-around term
      case 1:
        area += (V[n].y * (V[1].z - V[n-1].z));
        break;
      case 2:
        area += (V[n].z * (V[1].x - V[n-1].x));
        break;
      case 3:
        area += (V[n].x * (V[1].y - V[n-1].y));
        break;
    }

    // scale to get area before projection
    an = sqrt( ax*ax + ay*ay + az*az); // length of normal vector
    switch (coord) {
      case 1:
        area *= (an / (2 * N.x));
        break;
      case 2:
        area *= (an / (2 * N.y));
        break;
      case 3:
        area *= (an / (2 * N.z));
    }
    return fabs(area); // I add absolute value here
}

double calcAngle_c( Point v1, Point v2 )
{
    
    double length_v1 = v1.x*v1.x+v1.y*v1.y+v1.z*v1.z;
    double length_v2 = v2.x*v2.x+v2.y*v2.y+v2.z*v2.z;

    double dot = (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z)/sqrt(length_v1*length_v2);

    // Clipping to -1 to 1 for numerical stability
    dot = ( dot < -1.0 ? -1.0 : ( dot > 1.0 ? 1.0 : dot ) );

    double angle = acos(dot);
       
    return angle;
}

// see: 
// https://github.com/chtran/computer_vision/blob/master/proj3/code/vlfeat/vl/rodrigues.c
// https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
void calcRodriguesMtx_c(const double angle,const double* k_pt, double* R_pt)
{
    #define k(i)   k_pt[(i)]
    #define R(i,j) R_pt[(i)+3*(j)]

    const double verysmall = 1e-6;

    double kx = k_pt[0];
    double ky = k_pt[1];
    double kz = k_pt[2];

    double length_k = sqrt(kx*kx +
                           ky*ky +
                           kz*kz  ) ;

    if (length_k < verysmall)
    {
        R(0,0) = 1.0 ; R(0,1) = 0.0 ; R(0,2) = 0.0 ;
        R(1,0) = 0.0 ; R(1,1) = 1.0 ; R(1,2) = 0.0 ;
        R(2,0) = 0.0 ; R(2,1) = 0.0 ; R(2,2) = 1.0 ;

    }
    else 
    {
        double x = kx / length_k ;
        double y = ky / length_k ;
        double z = kz / length_k ;

        double xx = x*x ;
        double xy = x*y ;
        double xz = x*z ;
        double yy = y*y ;
        double yz = y*z ;
        double zz = z*z ;

        const double yx = xy ;
        // const double zx = xz ;
        // const double zy = yz ;

        double sth  = sin(-angle) ;
        double cth  = cos(-angle) ;
        double mcth = 1.0 - cth ;

        R(0,0) = 1          - mcth * (yy+zz) ;
        R(1,0) =     sth*z  + mcth * xy ;
        R(2,0) =   - sth*y  + mcth * xz ;

        R(0,1) =   - sth*z  + mcth * yx ;
        R(1,1) = 1          - mcth * (zz+xx) ;
        R(2,1) =     sth*x  + mcth * yz ;

        R(0,2) =     sth*y  + mcth * xz ;
        R(1,2) =   - sth*x  + mcth * yz ;
        R(2,2) = 1 - mcth * (xx+yy) ;

    }
        
    return;
}


void initPoint(Point *pt, double x, double y, double z)
{
    pt->x = x;
    pt->y = y;
    pt->z = z;
}

// JUNK


// void getPowerDelayMtxAll_c(int N,
//     BareDetector collector,PointSource emitter,
//     SimplePlane* planes, double* powMtx,double* delayMtx,
//     double* powVect_t,double* delayVect_t,
//     double* powVect_r,double* delayVect_r)
// {
//     /*
//         Refer to eq. (18). 
//         Assumptions:
//             m =1 
//             FoV = pi/2

//         Currently assume no blockage to test the speed.

//     */

//     #define isClose(a,b)   fabs((a) - (b)) <= (1e-08 + 1e-05 * fabs(b))
//     #define P(i,j) powMtx[(j)+N*(i)] // column-wise (row-major)
//     #define D(i,j) delayMtx[(j)+N*(i)] // column-wise (row-major)

//     const int speed_of_light = 299792458;

//     int i;
//     int k;

//     #pragma omp parallel for
//     for(i=0;i<N;i++)
//     { // Iteration over the row (collector)

//         // This is very important to define it inside the loop if you 
//         // want to use openmp
//         char isVisible; // Need the smallest size. Hence using char
//         double vartheta,psi;

//         char isVisible_t; // Need the smallest size. Hence using char
//         double vartheta_t,psi_t;

//         char isVisible_r; // Need the smallest size. Hence using char
//         double vartheta_r,psi_r;

//         // For t_f
//         isVisible_t = 1;
//         vartheta_t = calcAngle_c(emitter.normalVect,calcAdd_c(planes[i].ctrPoint,calcMult_c(-1,emitter.ctrPoint)));
//         psi_t = calcAngle_c(planes[i].normalVect,calcAdd_c(emitter.ctrPoint,calcMult_c(-1,planes[i].ctrPoint)));

//         // Visible only when 0 <= vartheta <= pi/2 and 0<= psi <= pi/2
//         if (!(vartheta_t>=0 && vartheta_t < M_PI/2) || !(psi_t>=0 && psi_t < M_PI/2))
//         {
//             isVisible_t = 0;
//         }
//         if (isClose(vartheta_t,M_PI/2) || isClose(psi_t,M_PI/2))
//         {
//             isVisible_t = 0;
//         }

//         // distance
//         Point dist = calcAdd_c(planes[i].ctrPoint,calcMult_c(-1,emitter.ctrPoint));
//         double d_ki_square = calcDot_c(dist,dist);

//         if(!isVisible_t)
//         {
//             powVect_t[i] = 0;
//             delayVect_t[i] = 0;
//         }
//         else
//         {
//             powVect_t[i] = ((emitter.m+1)/(2*M_PI))*pow(cos(vartheta_t),emitter.m)*planes[i].area*cos(psi_t)/d_ki_square;
//             delayVect_t[i] = sqrt(d_ki_square)/speed_of_light;
//         }
//         // End for t_f

//         // For r_f
//         isVisible_r = 1;
//         vartheta_r = calcAngle_c(planes[i].normalVect,calcAdd_c(collector.ctrPoint,calcMult_c(-1,planes[i].ctrPoint)));
//         psi_r = calcAngle_c(collector.normalVect,calcAdd_c(planes[i].ctrPoint,calcMult_c(-1,collector.ctrPoint)));

//         // Visible only when 0 <= vartheta <= pi/2 and 0<= psi <= pi/2
//         if (!(vartheta_r>=0 && vartheta_r < M_PI/2) || !(psi_r>=0 && psi_r < collector.FoV))
//         {
//             isVisible_r = 0;
//         }
//         if (isClose(vartheta_r,M_PI/2) || isClose(psi_r,collector.FoV))
//         {
//             isVisible_r = 0;
//         }

//         // distance
//         dist = calcAdd_c(collector.ctrPoint,calcMult_c(-1,planes[i].ctrPoint));
//         d_ki_square = calcDot_c(dist,dist);

//         if(!isVisible_r)
//         {
//             powVect_r[i] = 0;
//             delayVect_r[i] = 0;
//         }
//         else
//         {
//             powVect_r[i] = (1/M_PI)*cos(vartheta_r)*collector.area*cos(psi_r)/d_ki_square;
//             delayVect_r[i] = sqrt(d_ki_square)/speed_of_light;
//         }
//         // End for r_f

//         for(k=0;k<N;k++)
//         { // Iteration over the column (emitter)
//             if(i!=k)
//             {
//                 // Calculate visibility based on angles
//                 isVisible = 1;
//                 vartheta = calcAngle_c(planes[k].normalVect,calcAdd_c(planes[i].ctrPoint,calcMult_c(-1,planes[k].ctrPoint)));
//                 psi = calcAngle_c(planes[i].normalVect,calcAdd_c(planes[k].ctrPoint,calcMult_c(-1,planes[i].ctrPoint)));

//                 // Visible only when 0 <= vartheta <= pi/2 and 0<= psi <= pi/2
//                 if (!(vartheta>=0 && vartheta < M_PI/2) || !(psi>=0 && psi < M_PI/2))
//                 {
//                     isVisible = 0;
//                 }
//                 if (isClose(vartheta,M_PI/2) || isClose(psi,M_PI/2))
//                 {
//                     isVisible = 0;
//                 }

//                 // distance
//                 dist = calcAdd_c(planes[i].ctrPoint,calcMult_c(-1,planes[k].ctrPoint));
//                 d_ki_square = calcDot_c(dist,dist);
                
//                 if(!isVisible)
//                 {
//                     P(i,k) = 0;
//                     D(i,k) = 0;
//                 }
//                 else
//                 {
//                     P(i,k) = (1/M_PI)*cos(vartheta)*planes[i].area*cos(psi)/d_ki_square;
//                     D(i,k) = sqrt(d_ki_square)/speed_of_light;
//                 }
//             }
//             else
//             {
//                 P(i,k) = 0;
//                 D(i,k) = 0;
//             }
//         }
//     } 
// }

