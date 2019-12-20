#ifndef CIRUTILS_C_H
#define CIRUTILS_C_H
	
#include "utils.h"

    typedef struct CompletePlanes {
       Point normalVect; // normal vector
       Point ctrPoint; // center point
       Point verts[4]; // vertices
       double m; // Lambertian mode
       double area; // area
       double FoV; // area
       double reflectivity; // area
    } CompletePlane;

   extern void getPowerDelayMtxAll_c(int N,int Nblocking,
        CompletePlane collector,CompletePlane emitter,
        CompletePlane* planes, CompletePlane* blockingplanes, double* powMtx,double* delayMtx,
        double* powVect_t,double* delayVect_t,
        double* powVect_r,double* delayVect_r);

    // extern void calcCIRFreqDom_c(
	   //  int Nplanes,
	   //  int Nblocking,
	   //  int Ntx,
	   //  int Nrx,
	   //  int Nf,
	   //  double f[],
	   //  double tx_normalVects[],
	   //  double tx_ctrPoints[],
	   //  double tx_ms[],
	   //  double rx_normalVects[],
	   //  double rx_ctrPoints[],
	   //  double rx_areas[],
	   //  double rx_FoVs[],
	   //  double planes_normalVects[], 
	   //  double planes_ctrPoints[], 
	   //  double planes_areas[],
	   //  double planes_reflectivities[],
	   //  double blocking_normalVects[], 
	   //  double blocking_ctrPoints[], 
	   //  double blocking_verts[], 
	   //  double blocking_areas[],
	   //  double *Hlos_real,
	   //  double *Hlos_imag,
	   //  double *Hdiff_real,
	   //  double *Hdiff_imag
    // );

    extern int get_array_len(Point room, double delta_t, int MAX_BOUNCE);
    extern int get_array_loc(double t, double del_t);
    extern void getPowerDelayLOSGen_c(int Nblocking, CompletePlane collector, CompletePlane emitter, 
        CompletePlane* blockingplanes, double* powVect,double* delayVect);


    extern void zero_bounce_power(int Nblocking, double del_t,
        CompletePlane emitter,CompletePlane collector,
        CompletePlane* planes, CompletePlane* blockingplanes,
        double h_t[]);
    extern void first_bounce_matrix(int Nplanes, int Nblocking, int ARRAYLENGTH, int c_M, double del_t,
        CompletePlane emitter,CompletePlane collector,
        CompletePlane* planes, CompletePlane* blockingplanes, 
        double* THE_MATRIX, int* M_start, int* M_end);
    extern void update_CIR_RX(int Nplanes, int Nblocking, int ARRAYLENGTH,  double del_t, CompletePlane emitter,CompletePlane collector,
        CompletePlane* planes, CompletePlane* blockingplanes,int c_M, 
        double* THE_MATRIX, int* M_start, int* M_end, double h_t[]);
    extern void update_matrix(int Nplanes, int Nblocking, int ARRAYLENGTH,  double del_t,
        CompletePlane* planes, CompletePlane* blockingplanes,int* c_M,int* n_M, 
        double* THE_MATRIX, int* M_start, int* M_end);
    extern void calcCIRTimeDom_c(
        int Nplanes,
        int Nblocking,
        int Ntx,
        int Nrx,
        double timeSampling,
        int MAX_BOUNCE,
        int ARRAYLENGTH,
        double tx_normalVects[],
        double tx_ctrPoints[],
        double tx_ms[],
        double rx_normalVects[],
        double rx_ctrPoints[],
        double rx_areas[],
        double rx_FoVs[],
        double planes_normalVects[], 
        double planes_ctrPoints[], 
        double planes_areas[],
        double planes_reflectivities[],
        double blocking_normalVects[], 
        double blocking_ctrPoints[], 
        double blocking_verts[], 
        double blocking_areas[],
        double ht_los[], // don't forget to zero initialization
        double ht_diff[] // don't forget to zero initialization
        );

#endif