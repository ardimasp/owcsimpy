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
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

from . cimport cutils

ctypedef np.int_t DTYPE_int_t
ctypedef np.float64_t DTYPE_float64_t


def checkBlockage(plane1,plane2,planeB,nargout=1):
    """ To check whether two planes are blocked by another plane.

    Parameters
    ----------
    plane1,plane2,planeB: SimplePlane
        plane1 and plane2 are the tested planes that will be checked whether they 
        are blocked by planeB. SimplePlane refers to a tuple of normalVect, 
        ctrPoint, verts and area of a rectangular plane 
        (see :class:`~owcsimpy.geoobjects.bases.rectplane_py.RectPlane_py`).
    nargout: int (optional)
        If nargout

    Returns
    -------
    isBlocked: bool
        `True` denotes the blockage event occurs; otherwise, it returns `False`.
    intersectingPoint: ndarray(3,)
        The intersecting point on the planeB. 

    Notes
    -----
    The blocakge event is defined when the intersecting point of the line 
    segment from/to plane1 to/from plane2 intersects planeB and the intersecting 
    point is inside the convex hull of the planeB, and it lies between 
    the line segment of plane1 and plane2.

    Instead of using a convex hull, comparing the area of 4 triangles and the square
    is faster. The triangles are formed from the intersecting point and two vertices 
    of the plane.

    Examples
    --------

    .. plot::
        :format: doctest
        :include-source: True

        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> 
        >>> from owcsimpy.geoobjects.bases.rectplane_py import RectPlane_py as RectPlane
        >>> from owcsimpy.geoutils.draw import draw
        >>> from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
        >>> from owcsimpy.geoutils.cutils import checkBlockage
        >>> 
        >>> # Planes
        >>> planes = []
        >>> 
        >>> normalVect = Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)]),which='spherical')
        >>> ctrPoint = np.array([0,2,1])
        >>> planes.append(RectPlane(normalVect,ctrPoint=ctrPoint,
        >>>     RodriguesAngle=np.deg2rad(0),dimensions=[1.7,1]))
        >>> 
        >>> normalVect = Vector(np.array([1,np.deg2rad(90),np.deg2rad(180)]),which='spherical')
        >>> ctrPoint = np.array([4,2,2])
        >>> planes.append(RectPlane(normalVect,ctrPoint=ctrPoint,
        >>>     RodriguesAngle=np.deg2rad(0),dimensions=[1.7,1]))
        >>> 
        >>> # This should not block
        >>> normalVect = Vector(np.array([1,np.deg2rad(115),np.deg2rad(222)]),which='spherical')
        >>> ctrPoint = np.array([2,4,1])
        >>> planes.append(RectPlane(normalVect,ctrPoint=ctrPoint,RodriguesAngle=np.deg2rad(0),dimensions=[1.7,1]))
        >>> 
        >>> # This should block
        >>> normalVect = Vector(np.array([1,np.deg2rad(115),np.deg2rad(222)]),which='spherical')
        >>> ctrPoint = np.array([2,2,1])
        >>> planes.append(RectPlane(normalVect,ctrPoint=ctrPoint,RodriguesAngle=np.deg2rad(0),dimensions=[1.7,1]))
        >>> 
        >>> # This should not block
        >>> normalVect = Vector(np.array([1,np.deg2rad(90),np.deg2rad(180)]),which='spherical')
        >>> ctrPoint = np.array([5,2,3])
        >>> planes.append(RectPlane(normalVect,ctrPoint=ctrPoint,RodriguesAngle=np.deg2rad(0),dimensions=[1.7,1]))
        >>> 
        >>> # Prepare canvases
        >>> fig,axs = draw(subplots=True,figsize=(14,6),nrows=1,ncols=3,xlim=[0,5],ylim=[0,5],zlim=[0,3],
        >>>     azim=-32,elev=25)
        >>> 
        >>> for idx in range(2,5):
        >>>     colors = ['red','red','black']
        >>>     fig,axs[idx-2] = draw(figure=fig,axes=axs[idx-2],planes=[planes[0],planes[1],planes[idx]],
        >>>         colors=colors,facecolors=colors)
        >>> 
        >>>     ray = np.append(planes[0].ctrPoint,planes[1].ctrPoint).reshape(2,-1)
        >>> 
        >>>     # Draw ray
        >>>     axs[idx-2].plot(ray[:,0],ray[:,1],ray[:,2],color='k');
        >>> 
        >>>     # Get simple planes from each plane
        >>>     plane1 = planes[0].getSimplePlane()
        >>>     plane2 = planes[1].getSimplePlane()
        >>>     planeB = planes[idx].getSimplePlane()
        >>> 
        >>>     # Check whether the blockage occurs
        >>>     isBlocked,intersectingPoint = checkBlockage(plane1,plane2,planeB,nargout=2)
        >>> 
        >>>     if intersectingPoint is not None:
        >>>         xi,yi,zi = intersectingPoint
        >>>         axs[idx-2].scatter(xi,yi,zi)
        >>>         
        >>>     axs[idx-2].set_title("isBlocked: {}".format(isBlocked))
        >>> 
        >>> plt.show()


    """


    # Unpack the tuples
    # _,ctrPoint1,_,_ = plane1
    # _,ctrPoint2,_,_ = plane2
    # normalVectB,ctrPointB,vertsB,areaB = planeB
    _,ctrPoint1,*junk = plane1
    _,ctrPoint2,*junk = plane2
    normalVectB,ctrPointB,vertsB,areaB,*junk = planeB

    cdef Point ctrPoint1_c, ctrPoint2_c 
    cdef Point normalVectB_c, ctrPointB_c
    cdef Point vertsB_c[4]
    cdef double areaB_c 

    ctrPoint1_c = {'x':ctrPoint1[0],'y':ctrPoint1[1],'z':ctrPoint1[2]}
    ctrPoint2_c = {'x':ctrPoint2[0],'y':ctrPoint2[1],'z':ctrPoint2[2]}
    ctrPointB_c = {'x':ctrPointB[0],'y':ctrPointB[1],'z':ctrPointB[2]}
    normalVectB_c = {'x':normalVectB[0],'y':normalVectB[1],'z':normalVectB[2]}
    
    for idx in range(4):
        vertsB_c[idx] = {'x':vertsB[idx,0],'y':vertsB[idx,1],'z':vertsB[idx,2]}
    
    areaB_c = areaB

    isBlocked = checkBlockage_c(ctrPoint1_c, ctrPoint2_c, ctrPointB_c, normalVectB_c, vertsB_c, areaB_c)
    # totalarea = checkBlockage_c(ctrPoint1_c, ctrPoint2_c, ctrPointB_c, normalVectB_c, vertsB_c, areaB_c)

    # if np.allclose(areaB,totalarea):
    #     isBlocked = True
    # else:
    #     isBlocked = False

    if nargout == 1:
        return bool(isBlocked)
    elif nargout == 2:
    # elif nargout == 2 or nargout == 3:
        # FIXME: Calculation of intersectingPoint is pure python
        # Unpack the tuples
        # _,ctrPoint1,_,_ = plane1
        # _,ctrPoint2,_,_ = plane2
        # normalVectB,ctrPointB,vertsB,areaB = planeB
        _,ctrPoint1,*junk = plane1
        _,ctrPoint2,*junk = plane2
        normalVectB,ctrPointB,vertsB,areaB,*junk = planeB

        # Check intersection
        # see http://geomalgorithms.com/a05-_intersect-1.html
        u = ctrPoint2-ctrPoint1
        w = ctrPoint1-ctrPointB

        if np.allclose(np.dot(normalVectB,u),0):
            # Line segment and the planeB is parallel
            intersectingPoint = None
        else:

            # Parametric value of the intersecting point
            ti = -np.dot(normalVectB,w)/np.dot(normalVectB,u)
            
            # The intersecting point
            intersectingPoint = ctrPoint1+ti*u

        return bool(isBlocked),intersectingPoint
        # if nargout == 2:
        #     return bool(isBlocked),intersectingPoint
        # elif nargout == 3:
        #     return bool(isBlocked),intersectingPoint,totalarea

def calcArea3DPoly(int n,
    np.ndarray[dtype=double, ndim=1] Vx,
    np.ndarray[dtype=double, ndim=1] Vy,
    np.ndarray[dtype=double, ndim=1] Vz,
    np.ndarray[dtype=double, ndim=1] Nnp
    ):
    
    cdef Point *ps = <Point *> malloc((n+1) * sizeof(Point))
    cdef Point N = {'x':Nnp[0],'y':Nnp[1],'z':Nnp[2]}

    for idx in range(n+1):
        ps[idx] = {'x':Vx[idx],'y':Vy[idx],'z':Vz[idx]}
    
    area = abs(calcArea3DPoly_c(n,ps,N))

    free(ps)

    return area

def calcAngle(
    np.ndarray v1,
    np.ndarray v2):

    # Casting to float
    v1 = v1.astype(float)
    v2 = v2.astype(float)
  
    cdef Point v1_c = {'x':v1[0],'y':v1[1],'z':v1[2]}
    cdef Point v2_c = {'x':v2[0],'y':v2[1],'z':v2[2]}

    cdef double angle = calcAngle_c(v1_c,v2_c);

    return angle

def calcRodriguesMtx(double angle, np.ndarray k_in):

    cdef double *R_pt = <double *> malloc(9 * sizeof(double))
    
    cdef double k[3]
    k[:] = k_in
    cdef double *k_pt = k

    calcRodriguesMtx_c(angle, k_pt, R_pt);

    R = np.array([[R_pt[idx] for idx in range(9)]]).reshape(3,3)
    free(R_pt)
    
    return R
