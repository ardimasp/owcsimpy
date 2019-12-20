import copy
import math
import itertools
import numpy as np

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.rectplane_py import RectPlane_py as RectPlane
from owcsimpy.geoobjects.bases.paramline_py import ParamLine_py as Line
from owcsimpy.geoutils.cutils import calcRodriguesMtx as getRodriguesMtx
from owcsimpy.geoutils.cutils import calcAngle

class Cube_py(object):
    """A 3D cube.

    Parameters
    ----------
    normalVect: Vector_py
        Its normal vector.
    ctrPoint: ndarray(3,)
        Defines the center point of the cube.
    dimensions: array-like with 3 elements
        The dimensions are defined as, for example, [L,W,H], where L is the length,
        W is the width and H is the height.
    RodriguesAngle: float
        The angle of Rodrigues rotation w.r.t. the normal vector. This follows 
        the right-hand rule.
    invert: bool
        Invert is used to reverse the direction of the normal vectors of each faces.
        True means the direction is inward.
    identity: list
        Identity of the object. It is used to track down which parent object this object 
        is from. Mainly used for debugging purposes.
    ms: dict
        A dictionary of the Lambertian mode of each face
    FoVs: dict
        A dictionary of the field of view in rads.
    reflectivities: dict
        A dictionary of reflectivities of the surface.

    Attributes
    ----------
    normalVect: ndarray(3,)
    ctrPoint: ndarray(3,)
    listPlanes: list of :class:`~owcsimpy.geoobjects.bases.vector_py.RectPlane_py`
        List of planes on each direction.
    L: float
        Length
    W: float
        Width
    H: float
        Height
    identity: list
    ms: dict
    FoVs: dict
    reflectivities: dict

    Notes
    -----
    Unlike :class:`~owcsimpy.geoobjects.bases.vector_py.Vector_py` which 
    has two different ways of instantiating, it is sufficient in our 
    application to instantiate a cube with the angles information and 
    its dimensions.
    
    Following is the original position of a cube, i.e., the cube faces 
    upwards. 

    The directions of faces are defined based on 6 directions, i.e.,:
        *   B: Bottom, p0, downwards, (polar,az) = (180,0) deg
        *   T: Top, p1, upwards, (polar,az) = (0,0) deg
        *   S: South, p2, -, (polar,az) = (90,180) deg
        *   N: North, p3, -, (polar,az) = (90,0) deg
        *   E: East, p4, -, (polar,az) = (90,-90) deg
        *   W: West, p5, -, (polar,az) = (90,90) deg
    
    .. code-block:: text

                                       z                                            
                                       I                                            
                                       I                                            
                                       I                                            
                                       I                                            
                                       I                                            
                                       I                                            
                                       I                                            
                                p1 (T) I                                            
                                       I                                            
                                       I                                            
                                       I                                            
                                       I                                            
                                       I                                            
          p2 (S)               ;*IIVII*V                           y                
          ***;           ;**IVFVVVVVVVVMFVI*;                    *;                 
             ;****  **IVVFVVVVVVVVVVVVVMVVVVVFVII*;         ;****  ;*; p5 (W)       
               **I$FVVVVVVVVVVVVVVVVVVVMVVVVVVVVVVFVII*;;****  ;****                
               I*IVFVVVVVVVVVVVVVVVVVVVNVVVVVVVVVVVVVVVV$VI*I**;                    
               *;;;;*IIVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVFVVI                   
               *;;;;;;;;*IVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVII**                   
               *;;;;;;;;;;;;*IIVVVVVVVVVVVVVVVVVVVVVVVVVVII*;;;;*                   
               *;;;;;;;;;;;;;;;;**IVVVVVVVVVVVVVVVVVVII*;;;;;;;;*                   
               *;;;;;;;;;;;;;;;;;;;;*IIVVVVVVVVVFVI*;;;;;;;;;;;;*                   
                ;**;;;;;;;II*;;;;;;;;;;;;*IVFVI*;;;;;;;;;;;;;;;;*                   
                   ;*IIII*;;;;;;;;;;;;;;;;;;*;;;;;;;;;;II;;;;;;;*                   
                 ;**** ***;;;;;;;;;;;;;;;;;;*;;;;;;;;;;;;IIII***                    
        p4 (E)  ****;        ;***;;;;;;;;;;;;;;*;;;;;;;;;;;;;********               
                              ;**;;;;;;;;;;;*;;;;;;;;;;**;  p3 (N) ****             
                                 ;***;;;;;;;*;;;;;;***;               ****          
                                     ;*I;;;;*;;;***                       ;****     
                                       I;***I**;                               ; x  
                                       I    ;                                       
                                p0 (B) I                                            
                                                                                                                                                                        
    
    Examples
    --------
    .. plot::
        :format: doctest
        :include-source: True

        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> 
        >>> from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
        >>> from owcsimpy.geoobjects.bases.cube_py import Cube_py as Cube
        >>> from owcsimpy.geoutils.draw import draw
        >>> 
        >>> def genCube(polar,az,Rod,ctrPoint):
        >>>     
        >>>     cube = Cube(
        >>>         Vector(np.array([1,np.deg2rad(polar),np.deg2rad(az)])),
        >>>         ctrPoint = ctrPoint*np.ones(3),
        >>>         dimensions = [2,1,1],
        >>>         RodriguesAngle = np.deg2rad(Rod)
        >>>     )
        >>> 
        >>>     return cube
        >>> 
        >>> 
        >>> # Will draw 4 different canvases
        >>> fig,axs = draw(subplots=True,figsize=(14,6),nrows=1,ncols=4,xlim=[-2,2],ylim=[-2,2],zlim=[-2,2])
        >>> 
        >>> # Original position
        >>> polar,az,Rod=0,0,0 # polar, azimuth, Rodrigues
        >>> ctrPoint=0
        >>> cube = genCube(polar,az,Rod,ctrPoint)
        >>> 
        >>> fig,axs[0] = draw(figure=fig,axes=axs[0],cubes=cube,colors='blue',facecolors='blue')
        >>> axs[0].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,ctrPoint))
        >>> 
        >>> # Polar and azimuth
        >>> # Copy previous object as a reference (black)
        >>> fig,axs[1] = draw(figure=fig,axes=axs[1],cubes=cube,colors='black',facecolors='black')
        >>> 
        >>> polar,az,Rod=90,45,0 # polar, azimuth, Rodrigues
        >>> ctrPoint=0
        >>> cube = genCube(polar,az,Rod,ctrPoint)
        >>> 
        >>> fig,axs[1] = draw(figure=fig,axes=axs[1],cubes=cube,colors='blue',facecolors='blue')
        >>> axs[1].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,ctrPoint))
        >>> 
        >>> # Rodrigues
        >>> # Copy previous object as a reference (black)
        >>> fig,axs[2] = draw(figure=fig,axes=axs[2],cubes=cube,colors='black',facecolors='black')
        >>> 
        >>> polar,az,Rod=90,45,30 # polar, azimuth, Rodrigues
        >>> ctrPoint=0
        >>> cube = genCube(polar,az,Rod,ctrPoint)
        >>> 
        >>> fig,axs[2] = draw(figure=fig,axes=axs[2],cubes=cube,colors='blue',facecolors='blue')
        >>> axs[2].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,ctrPoint))
        >>> 
        >>> # Translation
        >>> # Copy previous object as a reference (black)
        >>> fig,axs[3] = draw(figure=fig,axes=axs[3],cubes=cube,colors='black',facecolors='black')
        >>> 
        >>> polar,az,Rod=90,45,30 # polar, azimuth, Rodrigues
        >>> ctrPoint=0.5
        >>> cube = genCube(polar,az,Rod,ctrPoint)
        >>> 
        >>> fig,axs[3] = draw(figure=fig,axes=axs[3],cubes=cube,colors='blue',facecolors='blue')
        >>> axs[3].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,ctrPoint))
        >>> 
        >>> 
        >>> plt.show()


    The left figure shows a cube in its original position (no rotation and translation). 
    Then, transformations are carried out w.r.t. the spherical coordinates rotation, 
    Rodrigues' rotation and translation in order. 

    """

    def __init__(self,normalVect,ctrPoint,dimensions,RodriguesAngle=0,invert=False,
        identity=[np.inf],ms={'p0':1,'p1':1,'p2':1,'p3':1,'p4':1,'p5':1},
        FoVs={'p0':np.pi/2,'p1':np.pi/2,'p2':np.pi/2,'p3':np.pi/2,'p4':np.pi/2,'p5':np.pi/2},
        reflectivities={'p0':1,'p1':1,'p2':1,'p3':1,'p4':1,'p5':1},):

        assert ctrPoint.size == 3 

        assert isinstance(ms,dict) and isinstance(FoVs,dict) and isinstance(reflectivities,dict), "must be a dictionary"
        assert len(ms) == 6 and len(FoVs) == 6 and len(reflectivities) == 6, (
           "the length must be 6")
        assert sorted( 
           [string.lower() for string in list(ms.keys())]
           ) == sorted(['p0','p1','p2','p3','p4','p5']), "keys name are wrong"

        assert sorted( 
           [string.lower() for string in list(FoVs.keys())]
           ) == sorted(['p0','p1','p2','p3','p4','p5']), "keys name are wrong"

        assert sorted( 
           [string.lower() for string in list(reflectivities.keys())]
           ) == sorted(['p0','p1','p2','p3','p4','p5']), "keys name are wrong"


        # Make sure that normalVect is the instance of Vector
        if not isinstance(normalVect,Vector):
            # If it is not the instance of Vector, assume it contains the cartesian 
            # coordinates
            normalVect = Vector(normalVect,which='cartesian')

        self.L,self.W,self.H = dimensions # must be size of 3

        # Normal vector and center point of the cube
        self.normalVect = normalVect()
        self.ctrPoint = ctrPoint

        # Directions of the normal vectors of each face
        if invert:
            polars   = np.deg2rad([  0,180, 90, 90, 90, 90]).tolist()
            azimuths = np.deg2rad([  0,  0,  0,180, 90,-90]).tolist()
        else:
            polars   = np.deg2rad([180,  0, 90, 90, 90, 90]).tolist()
            azimuths = np.deg2rad([  0,  0,180,  0,-90, 90]).tolist()

        # Cast dict of reflectivities into list.
        # Assume they are in order.
        rhos = list(reflectivities.values())

        # Generate planes for each face
        self.listPlanes=[]
        for idx in range(6):
            self.listPlanes.append(
                RectPlane(
                    Vector(np.array([1,polars[idx],azimuths[idx]])),
                    verts=getVertsFaces(idx,dimensions),
                    identity=[identity,idx],
                    m=list(ms.values())[idx],
                    FoV=list(FoVs.values())[idx],
                    reflectivity=rhos[idx]
                    )
                )

        # Rotation
        for idx in range(6):
            # Polar angle
            self.listPlanes[idx] = self.listPlanes[idx].rotate(normalVect.spherical[1],np.array([0,1,0]))
            # Azimuth angle
            self.listPlanes[idx] = self.listPlanes[idx].rotate(normalVect.spherical[2],np.array([0,0,1]))
            # Rodrigues angle
            if isinstance(RodriguesAngle,tuple):
                self.listPlanes[idx] = self.listPlanes[idx].rotate(RodriguesAngle[0],RodriguesAngle[1])
            else:
                self.listPlanes[idx] = self.listPlanes[idx].rotate(RodriguesAngle,self.normalVect)

        # Translation
        for idx in range(6):
            self.listPlanes[idx] = self.listPlanes[idx].translate(self.ctrPoint)

        self.identity = identity if isinstance(identity,list) else [identity]
        self.ms,self.FoVs,self.reflectivities = ms,FoVs,reflectivities

    def rotate(self,theta,refVector):
        """ Rotate the cube w.r.t. to a reference vector.

        Parameters
        ----------
        theta: float
        refVector: ndarray(3,) or Vector_py

        Returns
        -------
        Cube_py

        """
        if not isinstance(refVector,Vector):
            # If it is not the instance of Vector, assume it contains the cartesian 
            # coordinates
            refVector = Vector(refVector,which='cartesian')

        # Rodrigues' rotation
        Rd = getRodriguesMtx(theta,refVector())
        normalVect = (Rd@self.normalVect.reshape(3,1)).reshape(-1)

        newCube =  Cube_py(Vector(normalVect,which='cartesian'), ctrPoint=self.ctrPoint, 
            dimensions=[self.L,self.W,self.H],RodriguesAngle=(theta,refVector),
            identity=self.identity,ms=self.ms,FoVs=self.FoVs,
            reflectivities=self.reflectivities)

        newCube.identity = self.identity

        return newCube


    def translate(self,newCtrPoint):
        """ Translate the cube to a new center point.

        Parameters
        ----------
        newCtrPoint: ndarray(3,)

        Returns
        -------
        Cube_py

        """

        newCube = Cube_py(self.normalVect, ctrPoint=newCtrPoint, 
            dimensions=[self.L,self.W,self.H],
            identity=self.identity,ms=self.ms,FoVs=self.FoVs,
            reflectivities=self.reflectivities)

        newCube.identity = self.identity
        return newCube

    def getPartition(self,Ps=[1,1,1],delta=None):
        """Get a list of partitioned planes of each face.
        
        Parameters
        ----------
        Ps: list
            List of number of partition of each face. 
        delta: list
            Define the partition based on partition lengths

        Returns
        -------
        list:
            A list of partitioned planes. Each plane is an instant of 
            RectPlane_py.
        
        See Also
        --------
        :mod:`owcsimpy.geoobjects.bases.rectplane_py.RectPlane_py.getPartition`
        
        Examples
        --------

        .. plot::
            :format: doctest
            :include-source: True

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> 
            >>> from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
            >>> from owcsimpy.geoobjects.bases.cube_py import Cube_py as Cube
            >>> from owcsimpy.geoutils.draw import draw
            >>> cube = Cube(
            >>>     Vector(np.array([1,np.deg2rad(90),np.deg2rad(90)])),
            >>>     ctrPoint = np.array([0.5,0.5,0.5]),
            >>>     dimensions = [2,1,1],
            >>>     RodriguesAngle = np.deg2rad(30)
            >>> )
            >>> planes = cube.getPartition(delta=0.5)
            >>> fig,ax = draw(planes=planes,alphas=0.2,xlim=[-2,2],ylim=[-2,2],zlim=[-2,2])
            >>> plt.show()



        """

        if delta == None:
            # Casting to list if it is a constant
            Ps = [Ps] if not isinstance(Ps,list) else Ps

            assert isinstance(Ps,list)
            assert 0 < len(Ps) <= 3

            if len(Ps) == 1:
                #  If the length is one then replicate
                Ps = 3*Ps 
            else: 
                raise ValueError('Either assign one for all or all of them at once')

            # Unpack
            Px,Py,Pz = Ps
        else:
            # Casting to a list
            delta = [delta] if not isinstance(delta,list) else delta
            
            assert isinstance(delta,list)
            assert 0 < len(delta) <= 3

            if len(delta) == 1:
                #  If the length is one then replicate
                delta = 3*delta 
            else: 
                raise ValueError('Either assign one for all or all of them at once')

            # Unpack
            delta_x,delta_y,delta_z = delta

            assert 0 < delta_x <= self.L and 0 < delta_y <= self.W and 0 < delta_z <= self.H

            Px = int(self.L/delta_x)
            Py = int(self.W/delta_y)
            Pz = int(self.H/delta_z)

        partitionedPlanes = []

        Ps = [[Px,Py],[Px,Py],[Py,Pz],[Py,Pz],[Px,Pz],[Px,Pz]]

        for idx in range(6): # iteration over each face
            partitionedPlanes.append(self.listPlanes[idx].getPartition(Ps=Ps[idx]))

        return list(itertools.chain.from_iterable(partitionedPlanes))



def getVertsFaces(idx,dimensions):
    """ Get verticese of faces. 

    Implemented as a lookup table.
    
    Parameters
    ----------
    idx: int {0,1,2,3,4,5}
        w.r.t. the initial position
        idx:0 -> bottom face
        idx:1 -> top face
        idx:2 -> south face
        idx:3 -> north face
        idx:4 -> east face
        idx:5 -> west face

    """

    L,W,H = dimensions

    if idx == 0: # bottom face
        return np.array([
            [L/2,W/2,-H/2],
            [-L/2,W/2,-H/2],
            [-L/2,-W/2,-H/2],
            [L/2,-W/2,-H/2],
            ])
    elif idx == 1: # top face
        return np.array([
            [L/2,W/2,H/2],
            [-L/2,W/2,H/2],
            [-L/2,-W/2,H/2],
            [L/2,-W/2,H/2],
            ])
    elif idx == 2: # south face
        return np.array([
            [-L/2,W/2,H/2],
            [-L/2,-W/2,H/2],
            [-L/2,-W/2,-H/2],
            [-L/2,W/2,-H/2],
            ])
    elif idx == 3: # north face
        return np.array([
            [L/2,W/2,H/2],
            [L/2,-W/2,H/2],
            [L/2,-W/2,-H/2],
            [L/2,W/2,-H/2],
            ])
    elif idx == 4: # east face
        return np.array([
            [L/2,-W/2,H/2],
            [-L/2,-W/2,H/2],
            [-L/2,-W/2,-H/2],
            [L/2,-W/2,-H/2],
            ])
    elif idx == 5: # west face
        return np.array([
            [L/2,W/2,H/2],
            [-L/2,W/2,H/2],
            [-L/2,W/2,-H/2],
            [L/2,W/2,-H/2],
            ])
