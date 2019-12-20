#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <omp.h>

// #include "mkl_cblas.h"
// #include "mkl_lapacke.h"

#include "cirutils_c.h"
#include "utils.h"

/* Following is for the frequency domain approach.

*/


void getPowerDelayMtxAll_c(int N,int Nblocking,
    CompletePlane collector,CompletePlane emitter,
    CompletePlane* planes, CompletePlane* blockingplanes, double* powMtx,double* delayMtx,
    double* powVect_t,double* delayVect_t,
    double* powVect_r,double* delayVect_r)
{
    /* Calculate the power and delay matrices and vectors.
    
    Parameters
    ----------
    N: int
        Number of reflecting planes
    Nblocking: int
        Number of blocking planes
    collector: CompletePlane
        PD.
    emitter: CompletePlane
        LED.
    planes: CompletePlane*
        pointer to array of reflecting planes
    blockingplanes: CompletePlane*
        pointer to array of blocking planes

    Outputs
    -------
    powMtx, delayMtx: double*
        NxN matrices
    powVect_t, delayVect_t: double*
        Nx1 vector w.r.t. t_f
    powVect_r, delayVect_r: double*
        Nx1 vector w.r.t. r_f

    Notes
    -----
    OpenMP enabled
    Exploit symmetricity

    */
    

    #define isClose(a,b)   fabs((a) - (b)) <= (1e-08 + 1e-05 * fabs(b))
    #define P(i,j) powMtx[(j)+N*(i)] // column-wise (row-major)
    #define D(i,j) delayMtx[(j)+N*(i)] // column-wise (row-major)

    int i;
    int k;


    #pragma omp parallel for num_threads(omp_get_max_threads()) private(i,k) shared(powMtx,delayMtx,powVect_t,delayVect_t,powVect_r,delayVect_r,N,Nblocking,collector,emitter,planes,blockingplanes)
    for(i=0;i<N;i++)
    { // Iteration over the row (collector)

        // For t_f
        getPowerDelayLOSGen_c(Nblocking, planes[i], emitter, blockingplanes, &powVect_t[i],&delayVect_t[i]);

        // For r_f
        getPowerDelayLOSGen_c(Nblocking, collector, planes[i], blockingplanes, &powVect_r[i],&delayVect_r[i]);
        
        // For H
        for(k=i;k<N;k++) // Try symmetric, start from k=i
        { // Iteration over the column (emitter)
            if(i!=k)
            {
                // Calculate visibility based on angles
                getPowerDelayLOSGen_c(Nblocking, planes[i], planes[k], blockingplanes, &P(i,k),&D(i,k));
                
                // Adjustment for the symmetricity by means of
                // dividing by planes[i].area and multiplied by planes[k].area.
                // Recall that vartheta and psi are the same and m=1.
                P(k,i) = P(i,k)*planes[k].area/planes[i].area;
                D(k,i) = D(i,k);
                
            }
            else
            {
                P(i,k) = 0;
                D(i,k) = 0;
                P(k,i) = 0;
                D(k,i) = 0;
            }
        }
    } 
    // }
}


// void calcCIRFreqDom_c(
//     /* 
//     FIXME:
//     These crazy arguments should be shortened in the future.
//     Writing it like this makes it easier to interface with 
//     memviews in cython. 

//     */
//     int Nplanes,
//     int Nblocking,
//     int Ntx,
//     int Nrx,
//     int Nf,
//     double f[],
//     double tx_normalVects[],
//     double tx_ctrPoints[],
//     double tx_ms[],
//     double rx_normalVects[],
//     double rx_ctrPoints[],
//     double rx_areas[],
//     double rx_FoVs[],
//     double planes_normalVects[], 
//     double planes_ctrPoints[], 
//     double planes_areas[],
//     double planes_reflectivities[],
//     double blocking_normalVects[], 
//     double blocking_ctrPoints[], 
//     double blocking_verts[], 
//     double blocking_areas[],
//     double *Hlos_real,
//     double *Hlos_imag,
//     double *Hdiff_real,
//     double *Hdiff_imag
//     )
// {
//     /* Calculate CIR based on the frequency domain approach.
    
//     Notes
//     -----
//     This is a prototype for a full C implementation.
//     MKL is assumed.
//     This might be useful if we want to improve the code using CUDA.

//     */
    

//     #define H(i,j) H_f[(j)+Nplanes*(i)] // column-wise (row-major)
//     #define G(i,j) Grho[(j)+Nplanes*(i)] // column-wise (row-major)
//     #define Idt(i,j) IdtMtx[(j)+Nplanes*(i)] // column-wise (row-major)

//     // Mem alloc
//     // SimplePlane *planes = malloc(Nplanes*sizeof(SimplePlane));
//     // SimplePlane *blocking_planes = malloc(Nblocking*sizeof(SimplePlane));
//     CompletePlane *planes = malloc(Nplanes*sizeof(CompletePlane));
//     CompletePlane *blocking_planes = malloc(Nblocking*sizeof(CompletePlane));

//     double *powMtx = malloc(Nplanes*Nplanes*sizeof(double));
//     double *delayMtx = malloc(Nplanes*Nplanes*sizeof(double));
//     double *powVect_t = malloc(Nplanes*sizeof(double));
//     double *delayVect_t = malloc(Nplanes*sizeof(double));
//     double *powVect_r = malloc(Nplanes*sizeof(double));
//     double *delayVect_r = malloc(Nplanes*sizeof(double));
//     double powVect_los;
//     double delayVect_los;

//     MKL_Complex16 *H_f = malloc(Nplanes*Nplanes*sizeof(MKL_Complex16));
//     MKL_Complex16 *tmp_H_f = malloc(Nplanes*Nplanes*sizeof(MKL_Complex16));
//     MKL_Complex16 *Grho = calloc(Nplanes*Nplanes,sizeof(MKL_Complex16));
//     MKL_Complex16 *IdtMtx = calloc(Nplanes*Nplanes,sizeof(MKL_Complex16));
//     MKL_INT *ipiv = malloc(Nplanes*sizeof(MKL_INT));
//     MKL_Complex16 *t_f = malloc(Nplanes*sizeof(MKL_Complex16));
//     MKL_Complex16 *r_f = malloc(Nplanes*sizeof(MKL_Complex16));
//     MKL_Complex16 *s_f = malloc(Nplanes*sizeof(MKL_Complex16));
    

//     if (
//         planes == NULL || blocking_planes==NULL ||
//         powMtx == NULL || delayMtx==NULL ||
//         powVect_t == NULL || delayVect_t==NULL ||
//         powVect_r == NULL || delayVect_r==NULL ||
//         H_f == NULL || t_f==NULL || Grho == NULL ||
//         IdtMtx == NULL || s_f == NULL ||
//         r_f == NULL  
//         ) { 
//         printf("Memory allocation failed!\n"); 
//         exit(EXIT_FAILURE); 
//     }
//     else
//     {

//         int i,j,k;

//         #pragma omp parallel for num_threads(omp_get_max_threads()) private(i) shared(planes,IdtMtx,Grho,blocking_planes)
//         for(i=0;i<Nplanes;i++)
//         {
//             planes[i].normalVect.x = planes_normalVects[3*i+0];
//             planes[i].normalVect.y = planes_normalVects[3*i+1];
//             planes[i].normalVect.z = planes_normalVects[3*i+2];

//             planes[i].ctrPoint.x = planes_ctrPoints[3*i+0];
//             planes[i].ctrPoint.y = planes_ctrPoints[3*i+1];
//             planes[i].ctrPoint.z = planes_ctrPoints[3*i+2];

//             planes[i].area = planes_areas[i];
//             planes[i].m = 1;
//             planes[i].FoV = M_PI/2;
//             planes[i].reflectivity = planes_reflectivities[i];

//             Idt(i,i).real = 1;
//             G(i,i).real = planes_reflectivities[i];

//             if(i < Nblocking)
//             {
//                 blocking_planes[i].normalVect.x = blocking_normalVects[3*i+0];
//                 blocking_planes[i].normalVect.y = blocking_normalVects[3*i+1];
//                 blocking_planes[i].normalVect.z = blocking_normalVects[3*i+2];

//                 blocking_planes[i].ctrPoint.x = blocking_ctrPoints[3*i+0];
//                 blocking_planes[i].ctrPoint.y = blocking_ctrPoints[3*i+1];
//                 blocking_planes[i].ctrPoint.z = blocking_ctrPoints[3*i+2];

//                 blocking_planes[i].verts[0].x = blocking_verts[12*i+0];
//                 blocking_planes[i].verts[0].y = blocking_verts[12*i+1];
//                 blocking_planes[i].verts[0].z = blocking_verts[12*i+2];
//                 blocking_planes[i].verts[1].x = blocking_verts[12*i+0+3];
//                 blocking_planes[i].verts[1].y = blocking_verts[12*i+1+3];
//                 blocking_planes[i].verts[1].z = blocking_verts[12*i+2+3];
//                 blocking_planes[i].verts[2].x = blocking_verts[12*i+0+6];
//                 blocking_planes[i].verts[2].y = blocking_verts[12*i+1+6];
//                 blocking_planes[i].verts[2].z = blocking_verts[12*i+2+6];
//                 blocking_planes[i].verts[3].x = blocking_verts[12*i+0+9];
//                 blocking_planes[i].verts[3].y = blocking_verts[12*i+1+9];
//                 blocking_planes[i].verts[3].z = blocking_verts[12*i+2+9];

//                 blocking_planes[i].area = blocking_areas[i];

//             }
//         }

        
//         // PointSource emitter;
//         CompletePlane emitter;
//         emitter.normalVect.x = tx_normalVects[0];
//         emitter.normalVect.y = tx_normalVects[1];
//         emitter.normalVect.z = tx_normalVects[2];
//         emitter.ctrPoint.x = tx_ctrPoints[0];
//         emitter.ctrPoint.y = tx_ctrPoints[1];
//         emitter.ctrPoint.z = tx_ctrPoints[2];
//         emitter.m = tx_ms[0];

//         // BareDetector collector;
//         CompletePlane collector;
//         collector.normalVect.x = rx_normalVects[0];
//         collector.normalVect.y = rx_normalVects[1];
//         collector.normalVect.z = rx_normalVects[2];
//         collector.ctrPoint.x = rx_ctrPoints[0];
//         collector.ctrPoint.y = rx_ctrPoints[1];
//         collector.ctrPoint.z = rx_ctrPoints[2];
//         collector.area = rx_areas[0];
//         collector.FoV = rx_FoVs[0];

//         getPowerDelayMtxAll_c(Nplanes,Nblocking,collector,emitter,
//             planes,blocking_planes,powMtx,delayMtx,
//             powVect_t,delayVect_t,
//             powVect_r,delayVect_r);

//         // getPowerDelayLOS_c(Nblocking,collector,emitter, 
//         //     blocking_planes, &powVect_los,&delayVect_los);
//         getPowerDelayLOSGen_c(Nblocking,collector,emitter, 
//             blocking_planes, &powVect_los,&delayVect_los);

//         // print_matrix_double(powMtx,Nplanes,Nplanes);

//         MKL_Complex16 alpha,beta;
//         alpha.real = 1;
//         alpha.imag = 0;
//         beta.real = 0;
//         beta.imag = 0;

//         double fi;
//         MKL_INT info;
//         MKL_INT n = Nplanes;
//         MKL_INT nrhs = 1;
//         MKL_Complex16 Hdiff_tmp;

//         for(j=0; j<Nf;j++)
//         {
//             fi = f[j];
//             #pragma omp parallel for num_threads(omp_get_max_threads()) private(i,k) shared(t_f,r_f,H_f)
//             for(i=0; i<Nplanes; i++)
//             {
//                 t_f[i].real = powVect_t[i]*cos(-1*2*M_PI*delayVect_t[i]*fi);
//                 t_f[i].imag = powVect_t[i]*sin(-1*2*M_PI*delayVect_t[i]*fi);
//                 r_f[i].real = powVect_r[i]*cos(-1*2*M_PI*delayVect_r[i]*fi);
//                 r_f[i].imag = powVect_r[i]*sin(-1*2*M_PI*delayVect_r[i]*fi);
//                 for(k=0; k<Nplanes; k++)
//                 {
//                     H(i,k).real = powMtx[(k)+Nplanes*(i)]*cos(-1*2*M_PI*delayMtx[(k)+Nplanes*(i)]*fi);
//                     H(i,k).imag = powMtx[(k)+Nplanes*(i)]*sin(-1*2*M_PI*delayMtx[(k)+Nplanes*(i)]*fi);
//                 }
//             }


//             // https://software.intel.com/en-us/mkl-developer-reference-c-cblas-gemm
//             // H_f@Grho
//             cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
//                 Nplanes, Nplanes, Nplanes, &alpha, H_f, Nplanes, Grho, Nplanes,
//                       &beta, tmp_H_f, Nplanes);
//             // print_matrix_double_complex2(tmp_H_f,200,250,201,280,Nplanes);

//             // I-H_f@Grhp
//             beta.real = -1;
//             cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
//                 Nplanes, Nplanes, Nplanes, &alpha, IdtMtx, Nplanes, IdtMtx, Nplanes,
//                       &beta, tmp_H_f, Nplanes);

//             // https://software.intel.com/en-us/node/520973
//             // solve s_f with A=(I-H_f@Grhp), b= t_f and A@r_f=t_f
//             info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, n, nrhs,tmp_H_f, n,ipiv, t_f,nrhs);
//             // LAPACKE_zgesv(LAPACK_ROW_MAJOR, n, nrhs,tmp_H_f, n,ipiv, t_f,nrhs);

//             if(info != 0)
//                 exit(EXIT_FAILURE);

//             // Grho@s_f 
//             beta.real = 0;
//             cblas_zgemv(CblasRowMajor,CblasNoTrans,Nplanes,Nplanes,
//                  &alpha,Grho, Nplanes,t_f,1,&beta,s_f,1);

//             // Dot product, r_f.s_f
//             cblas_zdotu_sub(n, r_f, nrhs,s_f, nrhs, &Hdiff_tmp);

//             // Assign the outputs
//             Hdiff_real[j] = Hdiff_tmp.real;
//             Hdiff_imag[j] = Hdiff_tmp.imag;
//             Hlos_real[j] = powVect_los*cos(-1*2*M_PI*delayVect_los*fi);
//             Hlos_imag[j] = powVect_los*sin(-1*2*M_PI*delayVect_los*fi);
//         }

//         // Free memories
//         free(planes);
//         free(blocking_planes);

//         free(powMtx);
//         free(delayMtx);
//         free(powVect_t);
//         free(delayVect_t);
//         free(powVect_r);
//         free(delayVect_r);

//         free(H_f);
//         free(tmp_H_f);
//         free(IdtMtx);
//         free(Grho);
//         free(ipiv);
//         free(t_f);
//         free(r_f);
//         free(s_f);
        
//     } 


// }

/* Following is for the time domain approach.

*/

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) > (Y) ? (X) : (Y))
#define TM(i,j,k) THE_MATRIX[(j)+ARRAYLENGTH*(i)+ARRAYLENGTH*Nplanes*(k)] // column-wise (row-major)
#define MS(i,j) M_start[(j)+Nplanes*(i)] // column-wise (row-major)
#define ME(i,j) M_end[(j)+Nplanes*(i)] // column-wise (row-major)


int get_array_len(Point room, double delta_t, int MAX_BOUNCE)
{
    const int speed_of_light = 299792458;
    int ARRAYLENGTH = (MAX_BOUNCE+1)*sqrt(
        room.x*room.x + room.y *room.y  + room.z*room.z
        )/(delta_t*speed_of_light);

    return ARRAYLENGTH;
}

int get_array_loc(double t, double del_t)
{
    return t/del_t;
}

void getPowerDelayLOSGen_c(int Nblocking, CompletePlane collector, CompletePlane emitter, 
    CompletePlane* blockingplanes, double* powVect,double* delayVect)
{
    #define isClose(a,b)   fabs((a) - (b)) <= (1e-08 + 1e-05 * fabs(b))

    const int speed_of_light = 299792458;

    char isVisible; // Need the smallest size. Hence using char
    double vartheta,psi;

    int j;
    double isblocked;
    
    isVisible = 1;
    isblocked = 0;

    // Check whether emitter and collector are too close
    Point segment = calcAdd_c(collector.ctrPoint,calcMult_c(-1,emitter.ctrPoint));
    
    if(isClose(segment.x,0) && isClose(segment.y,0) && isClose(segment.z,0))
    {
        // They are too close. Assign 0 (invisible).
        isVisible = 0;
    }
    else
    { 
        vartheta = calcAngle_c(emitter.normalVect,calcAdd_c(collector.ctrPoint,calcMult_c(-1,emitter.ctrPoint)));
        psi = calcAngle_c(collector.normalVect,calcAdd_c(emitter.ctrPoint,calcMult_c(-1,collector.ctrPoint)));

        // Visible only when 0 <= vartheta <= pi/2 and 0<= psi <= pi/2
        if (!(vartheta>=0 && vartheta < M_PI/2) || !(psi>=0 && psi < collector.FoV))
        {
            isVisible = 0;
        }
        if (isClose(vartheta,M_PI/2) || isClose(psi,collector.FoV))
        {
            isVisible = 0;
        }

        // Check blockage 
        if(isVisible)
        {
            for(j=0;j<Nblocking;j++)
            {
                isblocked = checkBlockage_c(emitter.ctrPoint,collector.ctrPoint,blockingplanes[j].ctrPoint, 
                    blockingplanes[j].normalVect, blockingplanes[j].verts,blockingplanes[j].area);
                if(isblocked) break;
            }
        }
    }

    isVisible = isVisible && !isblocked;

    // distance
    Point dist = calcAdd_c(collector.ctrPoint,calcMult_c(-1,emitter.ctrPoint));
    double d_ki_square = calcDot_c(dist,dist);

    if(!isVisible)
    {
        powVect[0] = 0;
        delayVect[0] = 0;
    }
    else
    {
        powVect[0] = ((emitter.m+1)/(2*M_PI))*pow(cos(vartheta),emitter.m)*collector.area*cos(psi)/d_ki_square;
        delayVect[0] = sqrt(d_ki_square)/speed_of_light;
    }  

}


void zero_bounce_power(int Nblocking, double del_t,
    CompletePlane emitter,CompletePlane collector,
    CompletePlane* planes, CompletePlane* blockingplanes,
    double h_t[]) 
{

    double attenuation;
    double delay;
    getPowerDelayLOSGen_c(Nblocking,collector,emitter, 
            blockingplanes, &attenuation,&delay);

    if(attenuation!=0 && delay !=0)
        h_t[get_array_loc(delay,del_t)] += attenuation;

}

void first_bounce_matrix(int Nplanes, int Nblocking, int ARRAYLENGTH, int c_M, double del_t, 
    CompletePlane emitter,CompletePlane collector,
    CompletePlane* planes, CompletePlane* blockingplanes, 
    double* THE_MATRIX, int* M_start, int* M_end)
{
   
    int i;
    double attenuation;
    double delay;
    int t_i;
    CompletePlane plane;
    for(i=0;i<Nplanes;i++)
    {
        plane = planes[i];
        getPowerDelayLOSGen_c(Nblocking,plane,emitter, 
            blockingplanes, &attenuation,&delay);

        // printf("%.2e\n", delay);
        if(attenuation!=0 && delay !=0)
        {
            t_i = get_array_loc(delay,del_t);

            TM(i,t_i,c_M) += attenuation;
            MS(0,i) = min(MS(0,i),t_i);
            ME(0,i) = max(ME(0,i),t_i);
        }

    }
}

void update_CIR_RX(int Nplanes, int Nblocking, int ARRAYLENGTH, double del_t, CompletePlane emitter,CompletePlane collector,
    CompletePlane* planes, CompletePlane* blockingplanes,int c_M, 
    double* THE_MATRIX, int* M_start, int* M_end, double h_t[])
{
    int i;
    int count;
    double attenuation;
    double delay;
    int t_i;
    CompletePlane plane;

    for(i=0;i<Nplanes;i++)
    {
        plane = planes[i];
        getPowerDelayLOSGen_c(Nblocking,collector,plane, 
            blockingplanes, &attenuation,&delay);

        if(attenuation!=0 && delay !=0)
        {
            attenuation *= plane.reflectivity;
        
            t_i = get_array_loc(delay,del_t);

            for(count= MS(c_M,i);count < (ME(c_M,i) + 1);count++)
            {
                h_t[count+t_i] += attenuation*TM(i,count,c_M);
            }
        }

    }

}

void update_matrix(int Nplanes, int Nblocking, int ARRAYLENGTH, double del_t, 
    CompletePlane* planes, CompletePlane* blockingplanes,int* c_M,int* n_M, 
    double* THE_MATRIX, int* M_start, int* M_end)
{

    int i;
    int k;
    int count;
    int start_i;
    int end_i;
    double attenuation;
    double delay;
    int t_i;

    for(k=0;k<Nplanes;k++) 
    {
        /*
        Following the notation in the frequency domain, I tend to think of 
        k is the index for the transmitter and i is for the receiver.
        */
        for(i=0;i<Nplanes;i++)
        {
            // Calculate power and delay between k-th plane and i-th plane
            // TODO: Use getPowerDelayMtxAll_c instead and call this func
            // once. That should be more efficient.
            getPowerDelayLOSGen_c(Nblocking,planes[i],planes[k], 
                blockingplanes, &attenuation,&delay);

            
            if(attenuation!=0 && delay !=0)
            {
                attenuation *= planes[i].reflectivity;

                t_i = get_array_loc(delay,del_t);
                start_i =  MS(*c_M,k);
                end_i =  ME(*c_M,k);

                for(count=start_i;count<(end_i+1);count++)
                {
                    TM(i,count+t_i,*n_M) += attenuation*TM(k,count,*c_M);
                }

                MS(*n_M,i) = min(MS(*n_M,i),start_i+t_i);
                ME(*n_M,i) = max(ME(*n_M,i),end_i+t_i);
            }
        }


    }

    // Swap the active index of THE_MATRIX to be processed next
    int tmp = *n_M;
    *n_M = *c_M;
    *c_M = tmp;

}

void calcCIRTimeDom_c(
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
    )
{
    /*
    
    Returns
    -------
    ht_los: double[]
        CIR from LOS component
    ht_diff: double[]
        CIR from diffuse channel

    Notes
    -----
        Notations refer to https://github.com/UCaNLabUMB/CandLES whose comments,
        conventions, names, etc. are readable. So, I'll keep it the same. This 
        is also beneficial for those who are already familiar with the libraries 
        and want to compare it with the python version.

    */

    int c_M = 0; // index of currently active matrix
    int n_M = 1; // index of next or to-be-active matrix
    int bounce_no;
    int i;
    int j;

    // Mem alloc
    CompletePlane *planes = malloc(Nplanes*sizeof(CompletePlane));
    CompletePlane *blocking_planes = malloc(Nblocking*sizeof(CompletePlane));

    double* THE_MATRIX; // Store information about attenuation power
    // Store information of which index of CIR that would be assigned
    int* M_start; 
    int* M_end;
    if(MAX_BOUNCE == 1)
    {
        THE_MATRIX = calloc(Nplanes*ARRAYLENGTH*1,sizeof(double));
        M_start = calloc(1*Nplanes,sizeof(int));
        M_end = calloc(1*Nplanes,sizeof(int));
    }
    else // MAX_BOUNCE > 1
    {
        THE_MATRIX = calloc(Nplanes*ARRAYLENGTH*2,sizeof(double));
        M_start = calloc(2*Nplanes,sizeof(int));
        M_end = calloc(2*Nplanes,sizeof(int));
    }

    if(THE_MATRIX==NULL || M_start == NULL || M_end == NULL ||
        planes == NULL || blocking_planes == NULL)
    {
        printf("Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }
    else{

        if(MAX_BOUNCE == 1)
        {  
            for(i = 0; i< Nplanes; i++)
                // M_start[i] = ARRAYLENGTH-1;
                M_start[i] = ARRAYLENGTH;
        }
        else{
            for(i = 0; i< 2*Nplanes; i++)
                // M_start[i] = ARRAYLENGTH-1;
                M_start[i] = ARRAYLENGTH;
        }

        for(i=0;i<Nplanes;i++)
        {
            planes[i].normalVect.x = planes_normalVects[3*i+0];
            planes[i].normalVect.y = planes_normalVects[3*i+1];
            planes[i].normalVect.z = planes_normalVects[3*i+2];

            planes[i].ctrPoint.x = planes_ctrPoints[3*i+0];
            planes[i].ctrPoint.y = planes_ctrPoints[3*i+1];
            planes[i].ctrPoint.z = planes_ctrPoints[3*i+2];

            planes[i].area = planes_areas[i];
            planes[i].m = 1;
            planes[i].FoV = M_PI/2;
            planes[i].reflectivity = planes_reflectivities[i];

            if(i < Nblocking)
            {
                blocking_planes[i].normalVect.x = blocking_normalVects[3*i+0];
                blocking_planes[i].normalVect.y = blocking_normalVects[3*i+1];
                blocking_planes[i].normalVect.z = blocking_normalVects[3*i+2];

                blocking_planes[i].ctrPoint.x = blocking_ctrPoints[3*i+0];
                blocking_planes[i].ctrPoint.y = blocking_ctrPoints[3*i+1];
                blocking_planes[i].ctrPoint.z = blocking_ctrPoints[3*i+2];

                blocking_planes[i].verts[0].x = blocking_verts[12*i+0];
                blocking_planes[i].verts[0].y = blocking_verts[12*i+1];
                blocking_planes[i].verts[0].z = blocking_verts[12*i+2];
                blocking_planes[i].verts[1].x = blocking_verts[12*i+0+3];
                blocking_planes[i].verts[1].y = blocking_verts[12*i+1+3];
                blocking_planes[i].verts[1].z = blocking_verts[12*i+2+3];
                blocking_planes[i].verts[2].x = blocking_verts[12*i+0+6];
                blocking_planes[i].verts[2].y = blocking_verts[12*i+1+6];
                blocking_planes[i].verts[2].z = blocking_verts[12*i+2+6];
                blocking_planes[i].verts[3].x = blocking_verts[12*i+0+9];
                blocking_planes[i].verts[3].y = blocking_verts[12*i+1+9];
                blocking_planes[i].verts[3].z = blocking_verts[12*i+2+9];

                blocking_planes[i].area = blocking_areas[i];

            }
        }

        CompletePlane emitter;
        emitter.normalVect.x = tx_normalVects[0];
        emitter.normalVect.y = tx_normalVects[1];
        emitter.normalVect.z = tx_normalVects[2];
        emitter.ctrPoint.x = tx_ctrPoints[0];
        emitter.ctrPoint.y = tx_ctrPoints[1];
        emitter.ctrPoint.z = tx_ctrPoints[2];
        emitter.m = tx_ms[0];

        CompletePlane collector;
        collector.normalVect.x = rx_normalVects[0];
        collector.normalVect.y = rx_normalVects[1];
        collector.normalVect.z = rx_normalVects[2];
        collector.ctrPoint.x = rx_ctrPoints[0];
        collector.ctrPoint.y = rx_ctrPoints[1];
        collector.ctrPoint.z = rx_ctrPoints[2];
        collector.area = rx_areas[0];
        collector.FoV = rx_FoVs[0];

        // LOS channel
        zero_bounce_power(Nblocking,timeSampling,emitter,collector,
            planes,blocking_planes,ht_los); 

        // First order reflection
        first_bounce_matrix(Nplanes,Nblocking,ARRAYLENGTH,c_M,timeSampling, 
            emitter,collector,planes,blocking_planes, 
            THE_MATRIX, M_start,M_end);

        update_CIR_RX(Nplanes,Nblocking,ARRAYLENGTH,timeSampling,
            emitter,collector,
            planes,blocking_planes,c_M, 
            THE_MATRIX,M_start,M_end,ht_diff);
        
        // For higher order of reflections
        for(bounce_no=2;bounce_no<(MAX_BOUNCE+1);bounce_no++)
        {

            // Reset the index matrices
            for(i=0;i<Nplanes;i++)
            {
                // MS(n_M,i) = ARRAYLENGTH-1;
                MS(n_M,i) = ARRAYLENGTH;
                ME(n_M,i) = 0;
                for(j=0;j<ARRAYLENGTH;j++)
                    TM(i,j,n_M) = 0;
            }

            update_matrix(Nplanes,Nblocking,ARRAYLENGTH,timeSampling, 
                planes,blocking_planes,&c_M,&n_M, 
                THE_MATRIX,M_start,M_end);
            
            if(bounce_no >= 0)
            {
                update_CIR_RX(Nplanes,Nblocking,ARRAYLENGTH,timeSampling,
                    emitter,collector,
                    planes,blocking_planes,c_M, 
                    THE_MATRIX,M_start,M_end,ht_diff);

            }
        }
        
        free(THE_MATRIX);
        free(M_start);
        free(M_end);
    }

}

