import math
import itertools
import numpy as np

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.cube_py import Cube_py as Cube

class HumanCube_py(Cube):
    """ A simple human model with a 3D cube.

    The human is assumed to have a surface touching xy-plane at z=0. 
    Therefore, the location will be defined as 2D array.
    
    Dimensions are now defined such that the cube's original position has 
    the polar angle of 90deg and azimuth of 0 deg (see the notes below).

    HumanCube_py is inherited from 
    :class:`~owcsimpy.geoobjects.bases.cube_py.Cube_py`

    See Also
    --------
    :class:`~owcsimpy.geoobjects.bases.cube_py.Cube_py`

    Parameters
    ----------
    loc: ndarray(2,)
        Location in xy-plane.
    dimensions: ndarray(3,)
    direction: float
        The direction of the human is modeled as the azimuth angle.
    identity: int
    reflectivities: dict
        The keys are 'shirt' and 'hair'.

    Notes
    -----
    .. code-block:: text


                                     z                                                  
                                     ;                                                  
                                     V                                                  
                                     V                                                  
                          W          V ;**I**;                                          
                              ;;**IIVMFVVVVVVFFVIII**;;                                 
                     ;;*IIIVVFVVVVVVVMVVVVVVVVVVVVVVVV$$VI;                             
                *IVVVFVVVVVVVVVVVVVVVVVVVVVVVVVVVIII**;;;;*                             
                *;**IIIVFVVVVVVVVVVVVVVVFVIII**;;;;;;;;;;;;                             
                I;;;;;;;;**IIVVFVFVVII**;;;;;;;;;;;;;;;;;*                              
                *;;;;;;;;;;;;;;;I;;;;;;;;;;;;;;;;;;;;;;;;I                              
                *;;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;;;;;;;I                              
                *;;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;;;;;;;I                              
                *;;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;;;;;;;I                              
                *;;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;;;;;;;*                              
                ;;;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;;;;;;;I                              
                ;*;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;;;;;;;*                              
            H    *;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;;;;;;;*                              
                 I;;;;;;;;;;;;;;I;;;;;;;;;;;;;III;;;;;;;;I            y                  
                 I;;;;;;;;;;;;;;I;;;;;;;;;;;;;;;*IIIII;;;*        ***;                  
                 *;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;;;;*IIV**;*****                      
                 *;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;;;;;;;I*IVI***;                      
                 *;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;;;;;;;I       ******                 
                 *;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;;;;;;;*             ******           
                 *;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;;;;;;;*                  ;*****;     
                 I;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;;;;;**;                       ;** n   
                 *;;;;;;;;;;;;;;*;;;;;;;;;;;;;;;;;;****                                 
                 ***;;;;;;;;;;;;*;;;;;;;;;;;;;****;                                     
                   ;***;;;;;;;;;*;;;;;;;;;;***V**                                       
                       ****;;;;;I;;;;;;***;     ;****                                   
                   L       ****;I;****;              *****                              
                               *I*;                      ;****                          
                                                              ****;                     
                                                                  ****;                 
                                                                      ;****;            
                                                                           ****;        
                                                                               *****    
                                                                                    * x  
                                                                                

    Examples
    --------
    .. plot:: 
            :format: doctest
            :include-source: True

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> 
            >>> from owcsimpy.geoobjects.models.humancube_py import HumanCube_py as Human
            >>> from owcsimpy.geoutils.draw import draw
            >>> 
            >>> rho_keys = ['shirt','hair']
            >>> rho_vals = [0.3,0.6]
            >>> reflectivities = {rho_keys[i]:rho_vals[i] for i in range(2)}
            >>> 
            >>> persons = []
            >>> persons.append(Human(
            >>>     direction=np.deg2rad(45),
            >>>     loc=np.array([2,2]),
            >>>     dimensions=[0.3,0.5,1.7],
            >>>     identity=2,
            >>>     reflectivities=reflectivities))
            >>> persons.append(Human(
            >>>     direction=np.deg2rad(180),
            >>>     loc=np.array([4,1]),
            >>>     dimensions=[0.3,0.5,1.7],
            >>>     identity=2,
            >>>     reflectivities=reflectivities))
            >>> draw(models3d=persons,xlim=[0,5],ylim=[0,4],zlim=[0,3],
            >>>     facecolors=['red','blue'],colors=['red','blue']);
            >>>
            >>> plt.show() 

    """

    def __init__(self,direction,loc,dimensions,identity,
        reflectivities={'shirt':1,'hair':1}):

        assert isinstance(reflectivities,dict), "must be a dictionary"
        assert len(reflectivities) == 2, (
           "the length must be 2")
        assert sorted( 
           [string.lower() for string in list(reflectivities.keys())]
           ) == sorted(['shirt','hair']), "keys name are wrong"

        L,W,H = dimensions

        super().__init__(
            normalVect=Vector(np.array([1,np.pi/2,direction])),
            ctrPoint=np.array([loc[0],loc[1],H/2]),
            dimensions=np.array([H,W,L]),
            identity=identity,
            reflectivities={
            'p0':reflectivities['shirt'],
            'p1':reflectivities['shirt'],
            'p2':reflectivities['hair'],
            'p3':reflectivities['shirt'],
            'p4':reflectivities['shirt'],
            'p5':reflectivities['shirt']}
            )

    def getOuterVerts(self):
        """ Get outer vertices of a 3D object.
        
        This is used to check whether a 3D object intersects with another 3D object by 
        means of convex hull.

        """

        listVerts = [self.listPlanes[0].verts.tolist(),self.listPlanes[1].verts.tolist()]
        listVerts = list(itertools.chain.from_iterable(listVerts))

        return np.array(listVerts)

    def getPartition(self,Ps=[1,1,1],delta=None):
        """Get a list of partitioned planes of each face.

        Overriding the method from the parent class as we don't need the bottom face.

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


        # The bottom face is the fourth face
        for idx in [0,1,2,4,5]: # iteration over each face
            partitionedPlanes.append(self.listPlanes[idx].getPartition(Ps=Ps[idx]))

        return list(itertools.chain.from_iterable(partitionedPlanes))







