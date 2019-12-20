import numpy as np

def calcAngle_py(v1,v2):
    
    v1 = v1/np.linalg.norm(v1)
    v2 = v2/np.linalg.norm(v2)
    
    return np.arccos(np.dot(v1,v2).clip(-1,1))

def calcRodriguesMtx_py(theta,k):
    """
    Get Rodrigues' rotation matrix for 
    theta rad rotation (counter clockwise) w.r.t. a unit vector k
    
    see:
    https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    
    """
    
#     assert k.size == 3 and not np.allclose(k,np.zeros(3))
    
    # Ensure that k is a unit vector
    k = k/np.linalg.norm(k)
    
    # Cross-product matrix
    K = np.array([
            [    0,-k[2], k[1]],
            [ k[2],    0,-k[0]],
            [-k[1], k[0],   0]])
    
    return np.eye(3)+ K*np.sin(theta) + (K@K*(1-np.cos(theta)))

def calcAreaPoly_py(verts):
    """
    ref: http://geomalgorithms.com/a01-_area.html
    """

    L = verts.shape[0]
    assert L > 2, 'length of verts must be at least 3'

    # Get sequence of index
    idx1 = np.arange(L)
    idx2 = np.arange(1,L+1)
    idx2[-1] = 0

    # Get the unit normar vector
    n = np.cross(verts[1]-verts[0],verts[2]-verts[0])
    n /= np.linalg.norm(n)

    return abs(np.dot(np.sum(np.cross(verts[idx1],verts[idx2]),axis=0),n)/2)

    # assert L == 4, "Currently handle 4 vertices"

    # # Get the unit normar vector
    # n = np.cross(verts[1]-verts[0],verts[2]-verts[0])
    # n /= np.linalg.norm(n)

    # return abs(np.dot(np.cross(verts[2]-verts[0],verts[3]-verts[1]),n)/2)

def checkBlockage(plane1,plane2,planeB,nargout=1):
    """
    Check wheter plane1 and plane2 are blocked by planeB.

    plane1, plane2, planeB are simple planes defined in the 
    rectplane class.


    """

    from owcsimpy.geoutils.cutils import calcArea3DPoly

    # Unpack the tuples
    _,ctrPoint1,_,_ = plane1
    _,ctrPoint2,_,_ = plane2
    normalVectB,ctrPointB,vertsB,areaB = planeB

    # Check intersection
    # see http://geomalgorithms.com/a05-_intersect-1.html
    u = ctrPoint2-ctrPoint1
    w = ctrPoint1-ctrPointB

    # print(np.dot(normalVectB,u))

    if np.allclose(np.dot(normalVectB,u),0):
        # Line segment and the planeB is parallel
        isBlocked = False 
        intersectingPoint = None
    else:

        # Parametric value of the intersecting point
        ti = -np.dot(normalVectB,w)/np.dot(normalVectB,u)
        
        # The intersecting point
        intersectingPoint = ctrPoint1+ti*u

        # insert V[0] as V[n], see
        # http://geomalgorithms.com/a01-_area.html
        verts = np.append(vertsB,vertsB[0]).reshape(-1,3) # num of col is 3

        listTriangleVerts = [np.array([intersectingPoint.tolist(),verts[idx].tolist(),
                      verts[idx+1].tolist(),intersectingPoint.tolist()]) 
                     for idx in range(4)] 

        totalarea = sum([calcArea3DPoly(3,triVerts[:,0],
            triVerts[:,1],triVerts[:,2],normalVectB) 
        for triVerts in listTriangleVerts])

        if np.allclose(areaB,totalarea):
            isBlocked = True
        else:
            isBlocked = False

    if nargout == 1:
        return isBlocked
    elif nargout == 2:
        return isBlocked,intersectingPoint
    elif nargout == 3:
        return isBlocked,intersectingPoint,totalarea

def checkBlockageQHull(plane1,plane2,planeB,nargout=1):
    """
    Using ConvexHull.

    This is slower, but worth having as it's more intuitive.
    We can use this to have a third opinion.

    Check wheter plane1 and plane2 are blocked by planeB.

    plane1, plane2, planeB are simple planes defined in the 
    rectplane class.


    """

    from scipy.spatial import Delaunay

    # Unpack the tuples
    _,ctrPoint1,_,_ = plane1
    _,ctrPoint2,_,_ = plane2
    normalVectB,ctrPointB,vertsB,areaB = planeB

    # Check intersection
    # see http://geomalgorithms.com/a05-_intersect-1.html
    u = ctrPoint2-ctrPoint1
    w = ctrPoint1-ctrPointB

    if np.allclose(np.dot(normalVectB,u),0):
        # Line segment and the planeB is parallel
        isBlocked = False 
        intersectingPoint = None
    else:

        # Parametric value of the intersecting point
        ti = -np.dot(normalVectB,w)/np.dot(normalVectB,u)
        
        # The intersecting point
        intersectingPoint = ctrPoint1+ti*u

        verts = vertsB
        epsilon = 1e-3
        # for vert in vertsB:
        #     verts = np.append(verts,vert+normalVectB*epsilon)
        verts = np.append(verts,verts+normalVectB*epsilon)
        
        verts = verts.reshape(-1,3)

        delaunay = Delaunay(verts)
        isBlocked = delaunay.find_simplex(intersectingPoint)>=0 

    if nargout == 1:
        return isBlocked
    elif nargout == 2:
        return isBlocked,intersectingPoint



    