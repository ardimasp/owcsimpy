cdef extern from "cirutils_c.h":
    
    ctypedef struct Point:
        double x 
        double y
        double z

    ctypedef struct SimplePlane:
        Point normalVect # normal vector
        Point ctrPoint # center point
        Point verts[4] # vertices
        double area # area

    ctypedef struct PointSource:
        Point normalVect # normal vector
        Point ctrPoint # center point
        double m # Lambertian mode

    ctypedef struct BareDetector:
        Point normalVect # normal vector
        Point ctrPoint # center point
        double area # area
        double FoV # area

    ctypedef struct CompletePlane:
        Point normalVect # normal vector
        Point ctrPoint # center point
        Point verts[4] # vertices
        double area # area
        double m # area
        double FoV # area
        double reflectivity # area


    # void getPowerDelayLOS_c(int Nblocking,BareDetector collector, 
    #     PointSource emitter, SimplePlane* blockingplanes, 
    #     double* powVect,double* delayVect) nogil;
    # void getPowerDelayMtxAll_c(int N,int Nblocking,
    #     BareDetector collector,PointSource emitter,
    #     SimplePlane* planes, SimplePlane* blockingplanes, double* powMtx,double* delayMtx,
    #     double* powVect_t,double* delayVect_t,
    #     double* powVect_r,double* delayVect_r) nogil;

    void getPowerDelayMtxAll_c(int N,int Nblocking,
        CompletePlane collector,CompletePlane emitter,
        CompletePlane* planes, CompletePlane* blockingplanes, double* powMtx,double* delayMtx,
        double* powVect_t,double* delayVect_t,
        double* powVect_r,double* delayVect_r) nogil;

    void calcCIRFreqDom_c(
        int Nplanes,
        int Nblocking,
        int Ntx,
        int Nrx,
        int Nf,
        double f[],
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
        double *Hlos_real,
        double *Hlos_imag,
        double *Hdiff_real,
        double *Hdiff_imag
    ) nogil;

    void getPowerDelayLOSGen_c(int Nblocking, CompletePlane collector, CompletePlane emitter, 
        CompletePlane* blockingplanes, double* powVect,double* delayVect) nogil;

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
        double ht_los[], # don't forget to zero initialization
        double ht_diff[] # don't forget to zero initialization
    ) nogil;
