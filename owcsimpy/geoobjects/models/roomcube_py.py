import math
import itertools
import numpy as np

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.cube_py import Cube_py as Cube

class RoomCube_py(Cube):
    """ A room model using a cube.

    The room is defined on the postive axis.

    RoomCube_py is inherited from 
    :class:`~owcsimpy.geoobjects.bases.cube_py.Cube_py`
    
    Parameters
    ----------
    dimensions: ndarray(3,)
    identity: int
    reflectivities: dict
        The keys are 'b','t','s','n','e' and 'w'. 
        See the notes below for more description.


    See Also
    --------
    :class:`~owcsimpy.geoobjects.bases.cube_py.Cube_py`

    Notes
    -----
    Room configuration:

        *   B: Bottom
        *   T: Top
        *   S: South
        *   N: North
        *   E: East
        *   W: West

    .. code-block:: text

            Z
               $.I                                                  
           I:F$                                                 
            :              NNNNNNNNNNNNNNNNNNNNNNNNN            
            :           NNNNN                    NN$            
            :         NNN   N                  NNN $            
            :      NNN      N                 NN   $            
            :   NNN         N               NN     $            
            :NNN            N             NN       $            
            INNNNNNNNNNNNNNN$NNNNNNNNNNNNN$        $            
            I               N       Y     $        $            
            I               N   N$F       $        $            
            I               N  N:::       $        $            
            I               N $IFN$       $        $            
            I               VIV           $        $            
        H   I             NVVNNNNNNNNNNNNN$NNNNNNNN$            
            I           NVV               $       $             
            I          FVN                $     N$              
            I        $VN                  $    NN               
            I      $VN                    $   NN   W             
            I    NV$                      $  NN                 
            I  NVF                        $ $N                  
            INVV                          $$                 N  
           NIIV$$$$$$$$$$$$$$$$$$$$$$$$$$FIVVVVVVVVVVVVVVVVVI:.V   X
          o  $N            L                                 NIF 
                                                                
                                                            
                                                                
                                T                                
                                                                
                                IN                              
                              N:.:N            W                 
                             N:...:$   NFI*V                    
                            NI**.:*IF:.....N                    
                               $.F   $:...I                     
                               $.F  V:INF*$                     
                     NV$       $.FF:*N      FVN                 
                 NFI:..INNNNNNNF.::VNNNNNNNNI..:I$              
           S    NI:....:IIIIII*..*IIIIIIIIII:....:IN  N          
                   NFI:I    NI:I.F          I:IF                
                           V:IN$.F                              
                      V:IF:*N  $.F                              
                      :...:$NII*.*IIN                           
                     F...:*I$$.....F                            
                     VV$N     N:..$                             
                  E            N*N                              
                                 
                                B 

    Examples
    --------

    .. plot::
            :format: doctest
            :include-source: True

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> 
            >>> from owcsimpy.geoobjects.models.roomcube_py import RoomCube_py as Room
            >>> from owcsimpy.geoutils.draw import draw
            >>> 
            >>> rho_keys = ['b','t','s','n','e','w']
            >>> rho_vals = [0.1,0.5,0.5,0.5,0.5,0.5]
            >>> reflectivities = {rho_keys[i]:rho_vals[i] for i in range(6)}
            >>> 
            >>> fig,axs = draw(subplots=True,nrows=1,ncols=2,figsize=(14,8),
            >>>     xlim=[0,5],ylim=[0,4],zlim=[0,3]);
            >>> 
            >>> room = Room(dimensions=[5,4,3],identity=1,reflectivities=reflectivities)
            >>> 
            >>> # Partition with the delta of 1 m
            >>> subplanes = room.getPartition(delta=1)
            >>> 
            >>> fig,axs[0]=draw(figure=fig,axes=axs[0],models3d=room);
            >>> fig,axs[1]=draw(figure=fig,axes=axs[1],planes=subplanes,xlim=[0,5],ylim=[0,4],zlim=[0,3],
            >>>      alphas=[0.5*plane.reflectivity for plane in subplanes],
            >>>     lengths=0.3);
            >>>
            >>> plt.show()



    """
    def __init__(self,dimensions,
        identity,reflectivities={'b':1,'t':1,'s':1,'n':1,'e':1,'w':1}):

        assert isinstance(reflectivities,dict), "must be a dictionary"
        assert len(reflectivities) == 6, (
           "the length must be 6")
        assert sorted( 
           [string.lower() for string in list(reflectivities.keys())]
           ) == sorted(['b','t','s','n','e','w']), "keys name are wrong"

        L,W,H = dimensions
        super().__init__(
            normalVect=Vector(np.array([1,0,0])),
            ctrPoint=np.array([L/2,W/2,H/2]),
            dimensions=dimensions,
            invert=True,
            identity=identity,
            reflectivities={
            'p0':reflectivities['b'],
            'p1':reflectivities['t'],
            'p2':reflectivities['s'],
            'p3':reflectivities['n'],
            'p4':reflectivities['e'],
            'p5':reflectivities['w']},
            )

    def getOuterVerts(self):
        """ Get outer vertices of a 3D object.
        
        This is used to check whether a 3D object intersects with another 3D object by 
        means of convex hull.

        """

        listVerts = [self.listPlanes[0].verts.tolist(),self.listPlanes[1].verts.tolist()]
        listVerts = list(itertools.chain.from_iterable(listVerts))

        return np.array(listVerts)
