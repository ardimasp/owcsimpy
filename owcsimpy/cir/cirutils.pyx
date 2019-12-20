# distutils: language = c
# cython: cdivision = True
# cython: boundscheck = False
# cython: wraparound = False
# cython: profile = False

import cython
import numpy as np
import scipy as sp
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.math cimport M_PI,cos,sin
from libc.math cimport sqrt,pow,fabs
from scipy.linalg.cython_lapack cimport dsysv
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

from cython.parallel import prange, threadid

from scipy.linalg.lapack import get_lapack_funcs
from scipy.linalg import get_blas_funcs

from . cimport cirutils

ctypedef np.int_t DTYPE_int_t
ctypedef np.float64_t DTYPE_float64_t

cimport openmp

from owcsimpy.misc import flatten

def calcCIRFreqDom(f,emitter,collector,planes,blockingplanes=[],
    numReflections=np.inf,BLASLAPACK=False):
    """
    Parameters
    ----------
    f: list
        list of freq. bins.
    emitter: tuple
        a simple PointSource.
    collector: tuple
        a simple BareDetector.
    planes: list
        list of SimplePlane.
    blockingplanes: list, optional
        list of SimplePlane (default value is empty list).
    numReflections: int or inf
        Number of reflections.
    BLASLAPACK: bool
        Indicator to use BLAS or LAPACK libs from scipy.

    Returns
    -------
    Hlos: ndarray
        Freq. response for LOS channel.
    Hdiff: ndarray
        Freq. response for diffuse channel.

    Notes
    -----
    Refs: 
    .. [1] H. Schulze, "Frequency-Domain Simulation of the Indoor Wireless 
       Optical Communication Channel," in IEEE Transactions on Communications, 
       vol. 64, no. 6, pp. 2551-2562, June 2016.

    The implementation of this function is carried out such that we can store power
    and delay matrices between N planes. This is more suitable if we want to generate 
    a massive dataset due to we can later calculate at frequency bins we want to 
    calculate. 

    If you want to further optimize it, combine the loop in getPowerDelayMtxAll_c 
    and O(N^2) under the 'for idxf in range(Nf):' loop. Don't forget to vectorize 
    the operation and take advantage of map-lambda function from Python.

    """

    # Casting into a list
    if isinstance(f,np.ndarray):
        f = f.tolist()

    if not isinstance(f,list):
        f = [f]

    # Number of reflecting and blocking planes, respectively
    cdef int Nplanes = len(planes)
    cdef int Nblocking = len(blockingplanes)

    # Interface to C modules
    # Input
    cdef CompletePlane *planes_c = <CompletePlane *> PyMem_Malloc(Nplanes * sizeof(CompletePlane))
    cdef CompletePlane *blockingplanes_c = <CompletePlane *> PyMem_Malloc(Nblocking * sizeof(CompletePlane))
    # Output
    cdef double *powMtx = <double *> PyMem_Malloc(Nplanes*Nplanes * sizeof(double))
    cdef double *delayMtx = <double *> PyMem_Malloc(Nplanes*Nplanes * sizeof(double))
    cdef double *powVect_t = <double *> PyMem_Malloc(Nplanes * sizeof(double))
    cdef double *delayVect_t = <double *> PyMem_Malloc(Nplanes * sizeof(double))
    cdef double *powVect_r = <double *> PyMem_Malloc(Nplanes * sizeof(double))
    cdef double *delayVect_r = <double *> PyMem_Malloc(Nplanes * sizeof(double))
    
    cdef double powVect_LOS
    cdef double delayVect_LOS
    if (not powMtx or not delayMtx or not powVect_t or not delayVect_t 
        or not powVect_r or not delayVect_r):  
        raise MemoryError()

    cdef CompletePlane emitter_c
    cdef CompletePlane collector_c

    # Temporary variable on cython
    # The intrinsic transfer matrix Eq. (18)
    cdef complex [:,:] H_f = np.array(Nplanes*Nplanes*[1j*1.]).reshape(Nplanes,Nplanes) 
    # Transmitter transfer vector Eq. (17)
    cdef complex [:,:] t_f = np.array(Nplanes*[1j*1.]).reshape(Nplanes,1)
    # Receiver transfer vector Eq. (20)
    cdef complex [:,:] r_f = np.array(Nplanes*[1j*1.]).reshape(1,Nplanes)
    # Recursive vector Fig. 2
    cdef complex [:,:] s_f = np.array(Nplanes*[1j*1.]).reshape(Nplanes,1)


    # for the iterative method (refer to (53) in the paper)
    cdef complex [:,:] sM = np.array(Nplanes*[1j*1.]).reshape(Nplanes,1)
    cdef int M 
    cdef int RUN_INF_REFLECTIONS = 1
    if numReflections != np.inf :
        RUN_INF_REFLECTIONS = 0
        M = numReflections-1

    cdef int i,k
    cdef double[:] f_cy = np.array(f).astype("float64");
    cdef double fi
    cdef int Nf = len(f)
    cdef complex [:] Hdiff = np.array(Nf*[1j*1.])
    cdef complex [:] Hlos = np.array(Nf*[1j*1.])

    # Initialization
    for idx in range(Nplanes):
        """
        Recall that we define the simpler plane tuple as following.

        (normalVect: ndarray(3,) ---> i=0
            ctrPoint: ndarray(3,) ---> i=1 
            verts: ndarray(4,3) ---> i=2
            area: double    ---> i=3
            m: double   ---> i=4
            FoV: double ---> i=5
            reflectivity: double)   ---> i=6
        

        see :mod:`~owcsimpy.geoobjects.bases.vector_py.Vector_py.getSimplePlane`
        
        I will think how to tidy this part later. 

        """
        planes_c[idx].normalVect = {'x':planes[idx][0][0],'y':planes[idx][0][1],'z':planes[idx][0][2]}
        planes_c[idx].ctrPoint = {'x':planes[idx][1][0],'y':planes[idx][1][1],'z':planes[idx][1][2]}
        # Recall that there are 4 vertices
        planes_c[idx].verts[0] = {'x':planes[idx][2][0][0],'y':planes[idx][2][0][1],'z':planes[idx][2][0][2]}
        planes_c[idx].verts[1] = {'x':planes[idx][2][1][0],'y':planes[idx][2][1][1],'z':planes[idx][2][1][2]}
        planes_c[idx].verts[2] = {'x':planes[idx][2][2][0],'y':planes[idx][2][2][1],'z':planes[idx][2][2][2]}
        planes_c[idx].verts[3] = {'x':planes[idx][2][3][0],'y':planes[idx][2][3][1],'z':planes[idx][2][3][2]}
        planes_c[idx].area = planes[idx][3]
        
        # planes_c[idx].m = 1
        # planes_c[idx].FoV = np.pi/2
        planes_c[idx].m = planes[idx][4]
        planes_c[idx].FoV = planes[idx][5]
        planes_c[idx].reflectivity = planes[idx][6]

        if idx < Nblocking:
            blockingplanes_c[idx].normalVect = {'x':blockingplanes[idx][0][0],'y':blockingplanes[idx][0][1],'z':blockingplanes[idx][0][2]}
            blockingplanes_c[idx].ctrPoint = {'x':blockingplanes[idx][1][0],'y':blockingplanes[idx][1][1],'z':blockingplanes[idx][1][2]}
            # Recall that there are 4 vertices
            blockingplanes_c[idx].verts[0] = {'x':blockingplanes[idx][2][0][0],'y':blockingplanes[idx][2][0][1],'z':blockingplanes[idx][2][0][2]}
            blockingplanes_c[idx].verts[1] = {'x':blockingplanes[idx][2][1][0],'y':blockingplanes[idx][2][1][1],'z':blockingplanes[idx][2][1][2]}
            blockingplanes_c[idx].verts[2] = {'x':blockingplanes[idx][2][2][0],'y':blockingplanes[idx][2][2][1],'z':blockingplanes[idx][2][2][2]}
            blockingplanes_c[idx].verts[3] = {'x':blockingplanes[idx][2][3][0],'y':blockingplanes[idx][2][3][1],'z':blockingplanes[idx][2][3][2]}
            blockingplanes_c[idx].area = blockingplanes[idx][3]

    
    emitter_c.normalVect = {'x':emitter[0][0],'y':emitter[0][1],'z':emitter[0][2]}
    emitter_c.ctrPoint = {'x':emitter[1][0],'y':emitter[1][1],'z':emitter[1][2]}
    emitter_c.m = emitter[2]

    collector_c.normalVect = {'x':collector[0][0],'y':collector[0][1],'z':collector[0][2]}
    collector_c.ctrPoint = {'x':collector[1][0],'y':collector[1][1],'z':collector[1][2]}
    collector_c.area = collector[2]
    collector_c.FoV = collector[3]

    # Get power and delay values for H_f, t_f and r_f
    getPowerDelayMtxAll_c(Nplanes,Nblocking,collector_c,emitter_c,planes_c,blockingplanes_c,&powMtx[0],&delayMtx[0],
        &powVect_t[0],&delayVect_t[0],&powVect_r[0],&delayVect_r[0])
    
    getPowerDelayLOSGen_c(Nblocking,collector_c,emitter_c,blockingplanes_c,&powVect_LOS,&delayVect_LOS)

    # Get identity matrices
    cdef double [:] reflectivities = np.array([planes[idx][6] for idx in range(Nplanes)])
    Grho = np.diag(np.array(reflectivities))
    I = np.eye(Nplanes)
    makearray = np.array
    
    # Choose either numpy or scipy solution to solve linear equation system
    solve = np.linalg.solve
    # solve = sp.linalg.solve

    # Invoke necessary BLAS and LAPACK modules
    xgesvx = get_lapack_funcs(('gesvx'),dtype=np.complex128)
    xgemm = get_blas_funcs('gemm', dtype=np.complex128) # matrix-matrix mult
    xgemv = get_blas_funcs('gemv', dtype=np.complex128) # matrix-vector mult
    xdotu = get_blas_funcs('dotu', dtype=np.complex128) # dot prod

    
    for idxf in range(Nf):
        fi = f_cy[idxf]

        # LOS
        Hlos[idxf] = (powVect_LOS*
            cos(-1*2*M_PI*delayVect_LOS*fi)
            +1j*(powVect_LOS*sin(-1*2*M_PI*delayVect_LOS*fi)))

        # Diffuse
        # This can be 
        # for i in range(Nplanes):
        for i in prange(Nplanes, nogil=True,num_threads=openmp.omp_get_max_threads()):
            t_f[i,0] = (powVect_t[i]*
                cos(-1*2*M_PI*delayVect_t[i]*fi)
                +1j*(powVect_t[i]*sin(-1*2*M_PI*delayVect_t[i]*fi)))
            r_f[0,i] = (powVect_r[i]*
                cos(-1*2*M_PI*delayVect_r[i]*fi)
                +1j*(powVect_r[i]*sin(-1*2*M_PI*delayVect_r[i]*fi)))
            for k in range(Nplanes):
                H_f[i,k] = (powMtx[(k)+Nplanes*(i)]*
                    cos(-1*2*M_PI*delayMtx[(k)+Nplanes*(i)]*fi)
                    +1j*(powMtx[(k)+Nplanes*(i)]*sin(-1*2*M_PI*delayMtx[(k)+Nplanes*(i)]*fi)))
        

        if RUN_INF_REFLECTIONS:
            if not BLASLAPACK:
                # Using numpy
                s_f = solve(I-makearray(H_f)@Grho,makearray(t_f)) # --> this is the bottleneck
                """
                Basically, it uses LU decomposition approach. If we can exploit 
                positive definiteness of the matrix we can use another approach. 

                see: http://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html

                """
                Hdiff[idxf] = makearray(r_f)@Grho@makearray(s_f)
            else:
                # Using LAPACK and BLAS
                s_f = xgesvx(I-xgemm(1,makearray(H_f),Grho), makearray(t_f))[7]
                Hdiff[idxf] = xdotu(makearray(r_f),xgemv(1,Grho,makearray(s_f)))
        else: # if finite reflections

            # Initialization
            sM[:] = t_f

            if not BLASLAPACK:
                # Using numpy
                if M > 0:
                    for m in range(1,M+1):
                        sM = makearray(t_f) + (makearray(H_f)@Grho)@makearray(sM)

                Hdiff[idxf] = makearray(r_f)@Grho@makearray(sM)
            else:
                # Using LAPACK and BLAS
                if M > 0:
                    for m in range(1,M+1):
                        sM = makearray(t_f) + makearray(xgemv(1,xgemm(1,makearray(H_f),Grho),makearray(sM))).reshape(Nplanes,1)

                Hdiff[idxf] = xdotu(makearray(r_f),xgemv(1,Grho,makearray(sM)))


    PyMem_Free(planes_c)
    PyMem_Free(blockingplanes_c)
    PyMem_Free(powMtx)
    PyMem_Free(delayMtx)
    PyMem_Free(powVect_t)
    PyMem_Free(delayVect_t)
    PyMem_Free(powVect_r)
    PyMem_Free(delayVect_r)

    return np.array(Hlos),np.array(Hdiff)

# def calcCIRFreqDom_withc(f,emitter,collector,planes,blockingplanes):
#     """
#     This is a prototype if we want to use full C implementation. 
#     Might be useful if we want to extend it using CUDA.

#     """

#     if isinstance(f,np.ndarray):
#         f = f.tolist()

#     if not isinstance(f,list):
#         f = [f]

#     cdef int Nplanes = len(planes)
#     cdef int Nblocking = len(blockingplanes)
#     cdef int Ntx = 1
#     cdef int Nrx = 1
#     cdef int Nf = len(f)

#     planes_normalVects_i = list(flatten([planes[idx][0] for idx in range(len(planes))]))
#     planes_ctrPoints_i = list(flatten([planes[idx][1] for idx in range(len(planes))]))
#     planes_verts_i = list(flatten([planes[idx][2] for idx in range(len(planes))]))
#     planes_areas_i = list(flatten([planes[idx][3] for idx in range(len(planes))]))
#     planes_reflectivites_i = list(flatten([planes[idx][6] for idx in range(len(planes))]))

#     blocking_normalVects_i = list(flatten([blockingplanes[idx][0] for idx in range(len(blockingplanes))]))
#     blocking_ctrPoints_i = list(flatten([blockingplanes[idx][1] for idx in range(len(blockingplanes))]))
#     blocking_verts_i = list(flatten([blockingplanes[idx][2] for idx in range(len(blockingplanes))]))
#     blocking_areas_i = list(flatten([blockingplanes[idx][3] for idx in range(len(blockingplanes))]))

#     tx_normalVects_i = list(flatten([emitter[0]]))
#     tx_ctrPoints_i = list(flatten([emitter[1]]))
#     tx_ms_i = list(flatten([emitter[2]]))

#     rx_normalVects_i = list(flatten([collector[0]]))
#     rx_ctrPoints_i = list(flatten([collector[1]]))
#     rx_areas_i = list(flatten([collector[2]]))
#     rx_FoVs_i = list(flatten([collector[3]]))

#     cdef double[:] f_memview = np.array(f).astype("float64");
#     cdef double[:] tx_normalVects = np.array(tx_normalVects_i)
#     cdef double[:] tx_ctrPoints = np.array(tx_ctrPoints_i)
#     cdef double[:] tx_ms = np.array(tx_ms_i).astype("float64")
#     cdef double[:] rx_normalVects = np.array(rx_normalVects_i)
#     cdef double[:] rx_ctrPoints = np.array(rx_ctrPoints_i)
#     cdef double[:] rx_areas = np.array(rx_areas_i)
#     cdef double[:] rx_FoVs = np.array(rx_FoVs_i)
#     cdef double[:] planes_normalVects = np.array(planes_normalVects_i)
#     cdef double[:] planes_ctrPoints = np.array(planes_ctrPoints_i)
#     cdef double[:] planes_areas = np.array(planes_areas_i)
#     # cdef double[:] planes_reflectivities = np.array(reflectivities)
#     cdef double[:] planes_reflectivities = np.array(planes_reflectivites_i)
#     cdef double[:] blocking_normalVects = np.array(blocking_normalVects_i)
#     cdef double[:] blocking_ctrPoints = np.array(blocking_ctrPoints_i)
#     cdef double[:] blocking_verts = np.array(blocking_verts_i)
#     cdef double[:] blocking_areas = np.array(blocking_areas_i)

#     cdef double[:] Hlos_real = np.array(Nf*[1.])
#     cdef double[:] Hlos_imag = np.array(Nf*[1.])
#     cdef double[:] Hdiff_real = np.array(Nf*[1.])
#     cdef double[:] Hdiff_imag = np.array(Nf*[1.])

#     calcCIRFreqDom_c(Nplanes,Nblocking,Ntx,Nrx,Nf,&f_memview[0],
#         &tx_normalVects[0],&tx_ctrPoints[0],&tx_ms[0],
#         &rx_normalVects[0],&rx_ctrPoints[0],&rx_areas[0],&rx_FoVs[0],
#         &planes_normalVects[0],&planes_ctrPoints[0],&planes_areas[0],&planes_reflectivities[0],
#         &blocking_normalVects[0],&blocking_ctrPoints[0],&blocking_verts[0],&blocking_areas[0],
#         &Hlos_real[0],&Hlos_imag[0],&Hdiff_real[0],&Hdiff_imag[0]);

#     return (np.array(Hlos_real)+1j*np.array(Hlos_imag),
#         np.array(Hdiff_real)+1j*np.array(Hdiff_imag))

def getMtxFreqDom(emitter,collector,planes,blockingplanes=[]):
    """
    Parameters
    ----------
    emitter: tuple
        a simple PointSource.
    collector: tuple
        a simple BareDetector.
    planes: list
        list of SimplePlane.
    blockingplanes: list, optional
        list of SimplePlane (default value is empty list).

    Returns
    -------
    H_f: ndarray
        NxN transfer matrix with f = 1.
    t_f: ndarray
        Nx1 transmit vector with f = 1.
    r_f: ndarray
        Nx1 received vector with f = 1.

    """

    # Number of reflecting and blocking planes, respectively
    cdef int Nplanes = len(planes)
    cdef int Nblocking = len(blockingplanes)

    # Interface to C modules
    # Input
    cdef CompletePlane *planes_c = <CompletePlane *> PyMem_Malloc(Nplanes * sizeof(CompletePlane))
    cdef CompletePlane *blockingplanes_c = <CompletePlane *> PyMem_Malloc(Nblocking * sizeof(CompletePlane))
    # Output
    cdef double *powMtx = <double *> PyMem_Malloc(Nplanes*Nplanes * sizeof(double))
    cdef double *delayMtx = <double *> PyMem_Malloc(Nplanes*Nplanes * sizeof(double))
    cdef double *powVect_t = <double *> PyMem_Malloc(Nplanes * sizeof(double))
    cdef double *delayVect_t = <double *> PyMem_Malloc(Nplanes * sizeof(double))
    cdef double *powVect_r = <double *> PyMem_Malloc(Nplanes * sizeof(double))
    cdef double *delayVect_r = <double *> PyMem_Malloc(Nplanes * sizeof(double))
    
    cdef double powVect_LOS
    cdef double delayVect_LOS
    if (not powMtx or not delayMtx or not powVect_t or not delayVect_t 
        or not powVect_r or not delayVect_r):  
        raise MemoryError()

    cdef CompletePlane emitter_c
    cdef CompletePlane collector_c

    # Temporary variable on cython
    # The intrinsic transfer matrix Eq. (18)
    cdef complex [:,:] H_f = np.array(Nplanes*Nplanes*[1j*1.]).reshape(Nplanes,Nplanes) 
    # Transmitter transfer vector Eq. (17)
    cdef complex [:,:] t_f = np.array(Nplanes*[1j*1.]).reshape(Nplanes,1)
    # Receiver transfer vector Eq. (20)
    cdef complex [:,:] r_f = np.array(Nplanes*[1j*1.]).reshape(1,Nplanes)
    
    cdef int i,k
    
    # Initialization
    for idx in range(Nplanes):
        """
        Recall that we define the simpler plane tuple as following.

        (normalVect: ndarray(3,) ---> i=0
            ctrPoint: ndarray(3,) ---> i=1 
            verts: ndarray(4,3) ---> i=2
            area: double    ---> i=3
            m: double   ---> i=4
            FoV: double ---> i=5
            reflectivity: double)   ---> i=6
        

        see :mod:`~owcsimpy.geoobjects.bases.vector_py.Vector_py.getSimplePlane`
        
        I will think how to tidy this part later. 

        """
        planes_c[idx].normalVect = {'x':planes[idx][0][0],'y':planes[idx][0][1],'z':planes[idx][0][2]}
        planes_c[idx].ctrPoint = {'x':planes[idx][1][0],'y':planes[idx][1][1],'z':planes[idx][1][2]}
        # Recall that there are 4 vertices
        planes_c[idx].verts[0] = {'x':planes[idx][2][0][0],'y':planes[idx][2][0][1],'z':planes[idx][2][0][2]}
        planes_c[idx].verts[1] = {'x':planes[idx][2][1][0],'y':planes[idx][2][1][1],'z':planes[idx][2][1][2]}
        planes_c[idx].verts[2] = {'x':planes[idx][2][2][0],'y':planes[idx][2][2][1],'z':planes[idx][2][2][2]}
        planes_c[idx].verts[3] = {'x':planes[idx][2][3][0],'y':planes[idx][2][3][1],'z':planes[idx][2][3][2]}
        planes_c[idx].area = planes[idx][3]
        
        # planes_c[idx].m = 1
        # planes_c[idx].FoV = np.pi/2
        planes_c[idx].m = planes[idx][4]
        planes_c[idx].FoV = planes[idx][5]
        planes_c[idx].reflectivity = planes[idx][6]

        if idx < Nblocking:
            blockingplanes_c[idx].normalVect = {'x':blockingplanes[idx][0][0],'y':blockingplanes[idx][0][1],'z':blockingplanes[idx][0][2]}
            blockingplanes_c[idx].ctrPoint = {'x':blockingplanes[idx][1][0],'y':blockingplanes[idx][1][1],'z':blockingplanes[idx][1][2]}
            # Recall that there are 4 vertices
            blockingplanes_c[idx].verts[0] = {'x':blockingplanes[idx][2][0][0],'y':blockingplanes[idx][2][0][1],'z':blockingplanes[idx][2][0][2]}
            blockingplanes_c[idx].verts[1] = {'x':blockingplanes[idx][2][1][0],'y':blockingplanes[idx][2][1][1],'z':blockingplanes[idx][2][1][2]}
            blockingplanes_c[idx].verts[2] = {'x':blockingplanes[idx][2][2][0],'y':blockingplanes[idx][2][2][1],'z':blockingplanes[idx][2][2][2]}
            blockingplanes_c[idx].verts[3] = {'x':blockingplanes[idx][2][3][0],'y':blockingplanes[idx][2][3][1],'z':blockingplanes[idx][2][3][2]}
            blockingplanes_c[idx].area = blockingplanes[idx][3]

    
    emitter_c.normalVect = {'x':emitter[0][0],'y':emitter[0][1],'z':emitter[0][2]}
    emitter_c.ctrPoint = {'x':emitter[1][0],'y':emitter[1][1],'z':emitter[1][2]}
    emitter_c.m = emitter[2]

    collector_c.normalVect = {'x':collector[0][0],'y':collector[0][1],'z':collector[0][2]}
    collector_c.ctrPoint = {'x':collector[1][0],'y':collector[1][1],'z':collector[1][2]}
    collector_c.area = collector[2]
    collector_c.FoV = collector[3]

    # Get power and delay values for H_f, t_f and r_f
    getPowerDelayMtxAll_c(Nplanes,Nblocking,collector_c,emitter_c,planes_c,blockingplanes_c,&powMtx[0],&delayMtx[0],
        &powVect_t[0],&delayVect_t[0],&powVect_r[0],&delayVect_r[0])
    
    getPowerDelayLOSGen_c(Nblocking,collector_c,emitter_c,blockingplanes_c,&powVect_LOS,&delayVect_LOS)

    # Get identity matrices
    cdef double [:] reflectivities = np.array([planes[idx][6] for idx in range(Nplanes)])
    Grho = np.diag(np.array(reflectivities))
    I = np.eye(Nplanes)
    makearray = np.array
    
    # Choose either numpy or scipy solution to solve linear equation system
    solve = np.linalg.solve
    # solve = sp.linalg.solve

    # Invoke necessary BLAS and LAPACK modules
    xgesvx = get_lapack_funcs(('gesvx'),dtype=np.complex128)
    xgemm = get_blas_funcs('gemm', dtype=np.complex128) # matrix-matrix mult
    xgemv = get_blas_funcs('gemv', dtype=np.complex128) # matrix-vector mult
    xdotu = get_blas_funcs('dotu', dtype=np.complex128) # dot prod

   
    Hlos = (powVect_LOS*
            cos(delayVect_LOS)
            +1j*(powVect_LOS*sin(delayVect_LOS)))

    # This can be 
    # for i in range(Nplanes):
    for i in prange(Nplanes, nogil=True,num_threads=openmp.omp_get_max_threads()):
        t_f[i,0] = (powVect_t[i]*
            cos(delayVect_t[i])
            +1j*(powVect_t[i]*sin(delayVect_t[i])))
        r_f[0,i] = (powVect_r[i]*
            cos(delayVect_r[i])
            +1j*(powVect_r[i]*sin(delayVect_r[i])))
        for k in range(Nplanes):
            H_f[i,k] = (powMtx[(k)+Nplanes*(i)]*
                cos(delayMtx[(k)+Nplanes*(i)])
                +1j*(powMtx[(k)+Nplanes*(i)]*sin(delayMtx[(k)+Nplanes*(i)])))
    

    PyMem_Free(planes_c)
    PyMem_Free(blockingplanes_c)
    PyMem_Free(powMtx)
    PyMem_Free(delayMtx)
    PyMem_Free(powVect_t)
    PyMem_Free(delayVect_t)
    PyMem_Free(powVect_r)
    PyMem_Free(delayVect_r)

    return np.array(H_f),np.array(t_f),np.array(r_f),Hlos

def calcCIRFreqDom_fromMtx(f,Nplanes_saved,reflectivities_saved,
    H_f_saved,t_f_saved,r_f_saved,H_los_saved,numReflections=np.inf,BLASLAPACK=False):
    """
    Parameters
    ----------
    f: list
        list of freq. bins.
    emitter: tuple
        a simple PointSource.
    collector: tuple
        a simple BareDetector.
    planes: list
        list of SimplePlane.
    blockingplanes: list, optional
        list of SimplePlane (default value is empty list).
    numReflections: int or inf
        Number of reflections.
    BLASLAPACK: bool
        Indicator to use BLAS or LAPACK libs from scipy.

    Returns
    -------
    Hlos: ndarray
        Freq. response for LOS channel.
    Hdiff: ndarray
        Freq. response for diffuse channel.

    Notes
    -----
    Refs: 
    .. [1] H. Schulze, "Frequency-Domain Simulation of the Indoor Wireless 
       Optical Communication Channel," in IEEE Transactions on Communications, 
       vol. 64, no. 6, pp. 2551-2562, June 2016.

    The implementation of this function is carried out such that we can store power
    and delay matrices between N planes. This is more suitable if we want to generate 
    a massive dataset due to we can later calculate at frequency bins we want to 
    calculate. 

    If you want to further optimize it, combine the loop in getPowerDelayMtxAll_c 
    and O(N^2) under the 'for idxf in range(Nf):' loop. Don't forget to vectorize 
    the operation and take advantage of map-lambda function from Python.

    """

    # Casting into a list
    if isinstance(f,np.ndarray):
        f = f.tolist()

    if not isinstance(f,list):
        f = [f]

    # Number of reflecting and blocking planes, respectively
    cdef int Nplanes = Nplanes_saved

    #

    cdef double [:,:] powMtx = np.array(np.abs(H_f_saved)).reshape(Nplanes,Nplanes) 
    cdef double [:,:] delayMtx = np.array(np.angle(H_f_saved)).reshape(Nplanes,Nplanes) 
    cdef double [:,:] powVect_t = np.array(np.abs(t_f_saved)).reshape(Nplanes,1) 
    cdef double [:,:] delayVect_t = np.array(np.angle(t_f_saved)).reshape(Nplanes,1) 
    cdef double [:,:] powVect_r = np.array(np.abs(r_f_saved)).reshape(1,Nplanes) 
    cdef double [:,:] delayVect_r = np.array(np.angle(r_f_saved)).reshape(1,Nplanes) 
    cdef double powVect_LOS = np.abs(H_los_saved)
    cdef double delayVect_LOS = np.angle(H_los_saved)


    # Temporary variable on cython
    # The intrinsic transfer matrix Eq. (18)
    cdef complex [:,:] H_f = np.array(Nplanes*Nplanes*[1j*1.]).reshape(Nplanes,Nplanes) 
    # Transmitter transfer vector Eq. (17)
    cdef complex [:,:] t_f = np.array(Nplanes*[1j*1.]).reshape(Nplanes,1)
    # Receiver transfer vector Eq. (20)
    cdef complex [:,:] r_f = np.array(Nplanes*[1j*1.]).reshape(1,Nplanes)
    # Recursive vector Fig. 2
    cdef complex [:,:] s_f = np.array(Nplanes*[1j*1.]).reshape(Nplanes,1)


    # for the iterative method (refer to (53) in the paper)
    cdef complex [:,:] sM = np.array(Nplanes*[1j*1.]).reshape(Nplanes,1)
    cdef int M 
    cdef int RUN_INF_REFLECTIONS = 1
    if numReflections != np.inf :
        RUN_INF_REFLECTIONS = 0
        M = numReflections-1

    cdef int i,k
    cdef double[:] f_cy = np.array(f).astype("float64");
    cdef double fi
    cdef int Nf = len(f)
    cdef complex [:] Hdiff = np.array(Nf*[1j*1.])
    cdef complex [:] Hlos = np.array(Nf*[1j*1.])

    
    # Get identity matrices
    cdef double [:] reflectivities = reflectivities_saved
    Grho = np.diag(np.array(reflectivities))
    I = np.eye(Nplanes)
    makearray = np.array
    
    # Choose either numpy or scipy solution to solve linear equation system
    solve = np.linalg.solve
    # solve = sp.linalg.solve

    # Invoke necessary BLAS and LAPACK modules
    xgesvx = get_lapack_funcs(('gesvx'),dtype=np.complex128)
    xgemm = get_blas_funcs('gemm', dtype=np.complex128) # matrix-matrix mult
    xgemv = get_blas_funcs('gemv', dtype=np.complex128) # matrix-vector mult
    xdotu = get_blas_funcs('dotu', dtype=np.complex128) # dot prod

    
    for idxf in range(Nf):
        fi = f_cy[idxf]

        # LOS
        Hlos[idxf] = (np.abs(H_los_saved)*
            cos(-1*2*M_PI*np.angle(H_los_saved)*fi)
            +1j*(np.abs(H_los_saved)*sin(-1*2*M_PI*np.angle(H_los_saved)*fi)))

        # Diffuse
        # This can be 
        for i in range(Nplanes):
            t_f[i,0] = (powVect_t[i,0]*
                cos(-1*2*M_PI*delayVect_t[i,0]*fi)
                +1j*(powVect_t[i,0]*sin(-1*2*M_PI*delayVect_t[i,0]*fi)))
            r_f[0,i] = (powVect_r[0,i]*
                cos(-1*2*M_PI*delayVect_r[0,i]*fi)
                +1j*(powVect_r[0,i]*sin(-1*2*M_PI*delayVect_r[0,i]*fi)))
            for k in range(Nplanes):
                H_f[i,k] = (powMtx[i,k]*
                    cos(-1*2*M_PI*delayMtx[i,k]*fi)
                    +1j*(powMtx[i,k]*sin(-1*2*M_PI*delayMtx[i,k]*fi)))
        
       
        if RUN_INF_REFLECTIONS:
            if not BLASLAPACK:
                # Using numpy
                s_f = solve(I-makearray(H_f)@Grho,makearray(t_f)) # --> this is the bottleneck
                """
                Basically, it uses LU decomposition approach. If we can exploit 
                positive definiteness of the matrix we can use another approach. 

                see: http://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html

                """
                Hdiff[idxf] = makearray(r_f)@Grho@makearray(s_f)
            else:
                # Using LAPACK and BLAS
                s_f = xgesvx(I-xgemm(1,makearray(H_f),Grho), makearray(t_f))[7]
                Hdiff[idxf] = xdotu(makearray(r_f),xgemv(1,Grho,makearray(s_f)))
        else: # if finite reflections

            # Initialization
            sM[:] = t_f

            if not BLASLAPACK:
                # Using numpy
                if M > 0:
                    for m in range(1,M+1):
                        sM = makearray(t_f) + (makearray(H_f)@Grho)@makearray(sM)

                Hdiff[idxf] = makearray(r_f)@Grho@makearray(sM)
            else:
                # Using LAPACK and BLAS
                if M > 0:
                    for m in range(1,M+1):
                        sM = makearray(t_f) + makearray(xgemv(1,xgemm(1,makearray(H_f),Grho),makearray(sM))).reshape(Nplanes,1)

                Hdiff[idxf] = xdotu(makearray(r_f),xgemv(1,Grho,makearray(sM)))


    return np.array(Hlos),np.array(Hdiff)

def calcCIRTimeDom(timesampling,numReflections,furthestVert,
    emitter,collector,planes,blockingplanes):
    """
    Parameters
    ----------
    timesampling: double
        Time sampling in s
    numReflections: int
        Number of reflections.
    furthestVert: array-like(3,)
        The furthest vertice. This is used to calculate 
        the required array length of the CIR.
    emitter: tuple
        a simple PointSource.
    collector: tuple
        a simple BareDetector.
    planes: list
        list of SimplePlane.
    blockingplanes: list, optional
        list of SimplePlane (default value is empty list).

    Returns
    -------
    ht_los: ndarray
        CIR for LOS channel.
    ht_diff: ndarray
        CIR for diffuse channel.

    Notes
    -----
    Refs:
    .. [1] J. B. Carruthers and P. Kannan, "Iterative site-based modeling for 
       wireless infrared channels," in IEEE Transactions on Antennas and Propagation, 
       vol. 50, no. 5, pp. 759-765, May 2002. 
    .. [2] https://github.com/UCaNLabUMB/CandLES
    .. [3] http://iss.bu.edu/bwc/irsimit/

    """

    # Number of reflecting and blocking planes
    cdef int Nplanes = len(planes)
    cdef int Nblocking = len(blockingplanes)
    # Number of Tx and Rx
    # Can only handle 1 now
    cdef int Ntx = 1
    cdef int Nrx = 1
    
    # Convert into list
    planes_normalVects_i = list(flatten([planes[idx][0] for idx in range(len(planes))]))
    planes_ctrPoints_i = list(flatten([planes[idx][1] for idx in range(len(planes))]))
    planes_verts_i = list(flatten([planes[idx][2] for idx in range(len(planes))]))
    planes_areas_i = list(flatten([planes[idx][3] for idx in range(len(planes))]))
    planes_reflectivites_i = list(flatten([planes[idx][6] for idx in range(len(planes))]))

    blocking_normalVects_i = list(flatten([blockingplanes[idx][0] for idx in range(len(blockingplanes))]))
    blocking_ctrPoints_i = list(flatten([blockingplanes[idx][1] for idx in range(len(blockingplanes))]))
    blocking_verts_i = list(flatten([blockingplanes[idx][2] for idx in range(len(blockingplanes))]))
    blocking_areas_i = list(flatten([blockingplanes[idx][3] for idx in range(len(blockingplanes))]))

    tx_normalVects_i = list(flatten([emitter[0]]))
    tx_ctrPoints_i = list(flatten([emitter[1]]))
    tx_ms_i = list(flatten([emitter[2]]))

    rx_normalVects_i = list(flatten([collector[0]]))
    rx_ctrPoints_i = list(flatten([collector[1]]))
    rx_areas_i = list(flatten([collector[2]]))
    rx_FoVs_i = list(flatten([collector[3]]))

    # Get a cython objects from the list above
    cdef double del_t = np.array(timesampling).astype("float64")
    cdef int MAX_BOUNCE = numReflections
    cdef double[:] tx_normalVects = np.array(tx_normalVects_i).astype("float64")
    cdef double[:] tx_ctrPoints = np.array(tx_ctrPoints_i).astype("float64")
    cdef double[:] tx_ms = np.array(tx_ms_i).astype("float64")
    cdef double[:] rx_normalVects = np.array(rx_normalVects_i).astype("float64")
    cdef double[:] rx_ctrPoints = np.array(rx_ctrPoints_i).astype("float64")
    cdef double[:] rx_areas = np.array(rx_areas_i).astype("float64")
    cdef double[:] rx_FoVs = np.array(rx_FoVs_i).astype("float64")
    cdef double[:] planes_normalVects = np.array(planes_normalVects_i).astype("float64")
    cdef double[:] planes_ctrPoints = np.array(planes_ctrPoints_i).astype("float64")
    cdef double[:] planes_areas = np.array(planes_areas_i).astype("float64")
    cdef double[:] planes_reflectivities = np.array(planes_reflectivites_i).astype("float64")
    cdef double[:] blocking_normalVects = np.array(blocking_normalVects_i).astype("float64")
    cdef double[:] blocking_ctrPoints = np.array(blocking_ctrPoints_i).astype("float64")
    cdef double[:] blocking_verts = np.array(blocking_verts_i).astype("float64")
    cdef double[:] blocking_areas = np.array(blocking_areas_i).astype("float64")

    # Calculate the required array length for the CIR
    cdef int ARRAYLENGTH = np.int(
            numReflections + np.ceil(
                (numReflections+1)*np.linalg.norm(
                    np.array([furthestVert[0],furthestVert[1],furthestVert[2]]))/(timesampling*299792458) ))

    # Zero initialization (IMPORTANT!)
    cdef double[:] ht_los = np.zeros(ARRAYLENGTH)
    cdef double[:] ht_diff = np.zeros(ARRAYLENGTH)

    calcCIRTimeDom_c(
        Nplanes,
        Nblocking,
        Ntx,
        Nrx,
        del_t,
        MAX_BOUNCE,
        ARRAYLENGTH,
        &tx_normalVects[0],
        &tx_ctrPoints[0],
        &tx_ms[0],
        &rx_normalVects[0],
        &rx_ctrPoints[0],
        &rx_areas[0],
        &rx_FoVs[0],
        &planes_normalVects[0], 
        &planes_ctrPoints[0], 
        &planes_areas[0],
        &planes_reflectivities[0],
        &blocking_normalVects[0], 
        &blocking_ctrPoints[0], 
        &blocking_verts[0], 
        &blocking_areas[0],
        &ht_los[0], # don't forget to zero initialization
        &ht_diff[0] # don't forget to zero initialization
    )

    return (np.array(ht_los),np.array(ht_diff))

