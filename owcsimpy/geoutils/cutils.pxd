cdef extern from "utils.h":
    ctypedef struct Point:
        double x 
        double y
        double z
    # ctypedef struct SimplePlane:
    #     Point normalVect # normal vector
    #     Point ctrPoint # center point
    #     Point verts[4] # vertices
    #     double L # length
    #     double W # width
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


    double calcDot_c(Point p1, Point p2) nogil;
    Point calcAdd_c(Point p1, Point p2) nogil;
    Point calcMult_c(double a,Point p1) nogil;

    double checkBlockage_c(Point ctrPoint1, Point ctrPoint2, Point ctrPointB, 
        Point normalVectB, Point vertsB[4], const double areaB) nogil; 
    double calcArea3DPoly_c( int n, Point* V, Point N ) nogil;
    double calcAngle_c( Point v1, Point v2 ) nogil;
    void calcRodriguesMtx_c(const double angle, const double* k_pt, double* R_pt) nogil;

    # void getPowerDelayMtx_c(int N,SimplePlane* planes, double* powMtx,double* delayMtx) nogil;
    # void getPowerDelayVect_t_c(int N,PointSource emitter, SimplePlane* planes, double* powVect,double* delayVect) nogil;
    # void getPowerDelayVect_r_c(int N,BareDetector collector, SimplePlane* planes, double* powVect,double* delayVect) nogil;
    # void getPowerDelayLOS_c(int Nblocking,BareDetector collector, PointSource emitter, SimplePlane* blockingplanes, double* powVect,double* delayVect) nogil;
    # void getPowerDelayMtxAll_c(int N,int Nblocking,
    #     BareDetector collector,PointSource emitter,
    #     SimplePlane* planes, SimplePlane* blockingplanes, double* powMtx,double* delayMtx,
    #     double* powVect_t,double* delayVect_t,
    #     double* powVect_r,double* delayVect_r) nogil;

    # void testlapacke();
    
    # JUNK
    # void partitionWalls_c(int N, double Lx, double Wy, double Hz, int Px, int Py, int Pz, SimplePlane* planes)
    # void partitionWallsFaster_c(int N, double Lx, double Wy, double Hz, int Px, int Py, int Pz, SimplePlane* planes)
    # void visibilityMtx_c(int N, double Lx, double Wy, double Hz, int Px, int Py, int Pz, double* Mtx)
