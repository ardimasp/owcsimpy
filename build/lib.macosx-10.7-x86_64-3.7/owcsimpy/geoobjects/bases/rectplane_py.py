import copy
import math
import numpy as np

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.paramline_py import ParamLine_py as Line
from owcsimpy.geoutils.cutils import calcRodriguesMtx as getRodriguesMtx
from owcsimpy.geoutils.cutils import calcAngle

class RectPlane_py(object):
    """
    A 3D rectangular plane.

    Parameters
    ----------
    normalVect: Vector_py
        Its normal vector.
    ctrPoint: ndarray(3,)
        Defines the center point of the rectangular plane.
    RodriguesAngle: float
        The angle of Rodrigues rotation w.r.t. the normal vector. This follows 
        the right-hand rule.
    dimensions: array-like with 2 elements
        The dimensions are defined as, for example, [L,W], where L is the length
        and W is the width.
    verts: ndarray(4,3)
        The four vertices of the plane. 
    identity: list
        Identity of the object. It is used to track down which parent object this object 
        is from. Mainly used for debugging purposes.
    m: float
        The Lambertian mode.
    FoV: float
        Field of view in rads.
    reflectivity: float
        Reflectivity of the surface.

    Attributes
    ----------
    normalVect: ndarray(3,)
        Normal vector as a **ndarray** type
    ctrPoint: ndarray(3,)
    verts: ndarray(4,3)
    area: float
        Area of the square plane
    L: float
        Length
    W: float
        Width
    identity: list
    m: float
    FoV: float
    reflectivity: float
        

    Raises
    ------
    NotImplementedError:   
        When the parameters don't follow neither of two methods.


    Notes
    -----
    There are two different ways of instantiating a rectangular plane, i.e.,:
    
    1.   defining its normal vector, its center point, its Rodrigues angle 
         (w.r.t. its normal vector), and its dimensions (length and width), 
    2.   defining its four vertices.
    
    Following is the original position of a rectangular plane, i.e., 
    the plane faces upwards.

    .. code-block:: text

                                    z                                                        
                                    :                                               
                                    *                                               
                                    *                                               
                                    *                                               
                                    *                                               
                                    *                                               
                                    ^                                               
                               n   +++                                               
                                  +++++                                               
                                    *                                               
                                    *                              y                 
                                    *                          :::                  
                 W            :***: *                    :::::::                    
                      :*IIVFNN$$$$$NNVII*:         ::::::                           
               :**IV$$N$$$$$$$$$$$$$M$$$$NN$VVII*:::                                
          *IVF$N$$$$$$$$$$$$$$$$$$$$M$$NNNNNNN$$N$FVI*::                            
           :*IV$NN$$$$$$$$$$$$$$$$$$MNNN$$$$$$$$$$$$$$$N$VVI*:                      
                 *II$NN$$$$$$$$$$$$$$$$NNNNNNN$$$$$$$$$$$$$$NN$I                    
                      :*IF$N$$$$$$$$$$$$$$$$$NNNNNNN$$$NN$VI:                       
                           :*IV$N$$$$$$$$$$$$$$$$$$N$VVI::                          
                L             :*IV$NN$$$$$$$N$VI*:     :::::::                    
                                      *IV$$VI*:                :::::::              
                                                                     :::::::        
                                                                           ::::::  x 
                                                                                    

    See Also
    --------
    :mod:`~owcsimpy.geoutils.cutils.checkBlockage`
    
    Examples
    --------

    .. plot::
       :format: doctest
       :include-source: True

        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> 
        >>> from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
        >>> from owcsimpy.geoobjects.bases.rectplane_py import RectPlane_py as RectPlane
        >>> from owcsimpy.geoutils.draw import draw
        >>> 
        >>> def genPlane(polar,az,Rod,translation):
        >>>     
        >>>     v = Vector(np.array([1,np.deg2rad(polar),np.deg2rad(az)]))
        >>>     ctrPoint=np.array(3*[translation])
        >>>     plane = RectPlane(normalVect=v,ctrPoint=ctrPoint,
        >>>         RodriguesAngle=np.deg2rad(Rod),dimensions=[2,1])
        >>> 
        >>>     return plane
        >>> 
        >>> 
        >>> # Will draw 4 different canvases
        >>> fig,axs = draw(subplots=True,figsize=(14,6),nrows=1,ncols=4,xlim=[-2,2],ylim=[-2,2],zlim=[-2,2])
        >>> 
        >>> # Original position
        >>> polar,az,Rod=0,0,0 # polar, azimuth, Rodrigues
        >>> translation=0
        >>> plane = genPlane(polar,az,Rod,translation)
        >>> 
        >>> fig,axs[0] = draw(figure=fig,axes=axs[0],planes=plane,colors='blue',facecolors='blue')
        >>> axs[0].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,translation))
        >>> 
        >>> # Polar and azimuth
        >>> # Copy previous object as a reference (black)
        >>> fig,axs[1] = draw(figure=fig,axes=axs[1],planes=plane,colors='black',facecolors='black')
        >>> 
        >>> polar,az,Rod=90,45,0 # polar, azimuth, Rodrigues
        >>> translation=0
        >>> plane = genPlane(polar,az,Rod,translation)
        >>> 
        >>> fig,axs[1] = draw(figure=fig,axes=axs[1],planes=plane,colors='blue',facecolors='blue')
        >>> axs[1].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,translation))
        >>> 
        >>> # Rodrigues
        >>> # Copy previous object as a reference (black)
        >>> fig,axs[2] = draw(figure=fig,axes=axs[2],planes=plane,colors='black',facecolors='black')
        >>> 
        >>> polar,az,Rod=90,45,30 # polar, azimuth, Rodrigues
        >>> translation=0
        >>> plane = genPlane(polar,az,Rod,translation)
        >>> 
        >>> fig,axs[2] = draw(figure=fig,axes=axs[2],planes=plane,colors='blue',facecolors='blue')
        >>> axs[2].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,translation))
        >>> 
        >>> # Translation
        >>> # Copy previous object as a reference (black)
        >>> fig,axs[3] = draw(figure=fig,axes=axs[3],planes=plane,colors='black',facecolors='black')
        >>> 
        >>> polar,az,Rod=90,45,30 # polar, azimuth, Rodrigues
        >>> translation=0.5
        >>> plane = genPlane(polar,az,Rod,translation)
        >>> 
        >>> fig,axs[3] = draw(figure=fig,axes=axs[3],planes=plane,colors='blue',facecolors='blue')
        >>> axs[3].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,translation))
        >>> 
        >>> 
        >>> plt.show()


    Angles are defined as polar, azimuth and Rodrigues' angles in order. The left figure shows 
    a rectangular plane in its original position (no rotation and translation). Then, transformations 
    are carried out w.r.t. the spherical coordinates rotation, Rodrigues' rotation and translation 
    in order. 
    
    """


    def __init__(self,normalVect,
                 ctrPoint=None,RodriguesAngle=0,dimensions=None,
                 verts=None,
                 identity=[np.inf],m=1,FoV=np.pi/2,reflectivity=1):
        
        # assert isinstance(normalVect,Vector)
        # assert (RodriguesAngle != None and dimensions != None) or (verts != None), (
        # assert (RodriguesAngle != None and dimensions != None) or (verts != None), (
            # "Please specify following the first or the second method")
        
        assert m >= 0 and 0 <= FoV <= np.pi/2

        # Make sure that normalVect is the instance of Vector
        if not isinstance(normalVect,Vector):
            # If it is not the instance of Vector, assume it contains the cartesian 
            # coordinates
            normalVect = Vector(normalVect,which='cartesian')

        if (dimensions != None):
            # The first method
            
            # Initial
            L,W = dimensions # must be size of 2
            
            x = np.array([ L/2,-L/2,-L/2, L/2])
            y = np.array([ W/2, W/2,-W/2,-W/2])
            z = np.array([   0,   0,   0,   0])
            
            verts = np.array([[xi,yi,zi] for xi,yi,zi in zip(x,y,z)])
            
            # Get polar and azimuth angles of the normalVect
            polar,azimuth = normalVect.spherical[1:]
            
            # Rotate based on polar and azimuth
            # Get rotation matrix
            R = getRodriguesMtx(azimuth,np.array([0,0,1]))@getRodriguesMtx(polar,np.array([0,1,0]))
            verts = np.array([(R@vert.reshape(3,1)).reshape(-1) for vert in verts])
            
            # Rodrigues' rotation
            Rd = getRodriguesMtx(RodriguesAngle,normalVect())
            verts = np.array([(Rd@vert.reshape(3,1)).reshape(-1) for vert in verts])
            
            # Translate each verts
            verts = verts+np.array(ctrPoint)

            self.normalVect = normalVect()
            self.ctrPoint = np.array(ctrPoint) # make sure that it is ndarray
            self.verts = verts
            self.area = L*W
            self.L,self.W = L,W
            
        # elif verts != None:
        elif isinstance(verts,np.ndarray):
            # The second method
            self.normalVect = normalVect()
            self.verts = verts
            self.area = np.linalg.norm(verts[0]-verts[1])*np.linalg.norm(
                verts[1]-verts[2]) # assuming that it is a rectangular
            
            # Adjust the center point
            self.ctrPoint = (verts[0]+verts[2])/2

            self.L = np.linalg.norm(verts[0]-verts[1])
            self.W = np.linalg.norm(verts[1]-verts[2])
        
        else:
            raise NotImplementedError('Undefined method!')
            
        self.identity = identity if isinstance(identity,list) else [identity]
        self.m,self.FoV,self.reflectivity = m,FoV,reflectivity

        # store the normal vector as an instance of the Vector class
        self.__normalVect = normalVect

        # Check whether normalVect and verts form a plane
        v1 = self.normalVect
        
        checkPlane = [math.isclose(np.pi/2,calcAngle(v1,vert-self.ctrPoint)) 
                      for vert in self.verts]
        
        """
        A valid rectangular plane is the normal vector and 
        the vectors formed by the vertices and the center point are 
        orthogonal. Also, the lenghts of the opposite sides must be 
        the same.

        """
        if not np.array(checkPlane).all() and (
            np.allclose(verts[0]-verts[1],verts[2]-verts[3])) and (
            np.allclose(verts[1]-verts[2],verts[3]-verts[0])):
            print("The vertices are {}".format(self.verts))
            raise ValueError('Not a valid plane!')
    
    def rotate(self,theta,refVector):
        """ Rotate the plane w.r.t. to a reference vector.

        Parameters
        ----------
        theta: float
        refVector: ndarray(3,) or Vector_py

        Returns
        -------
        RectPlane_py

        Examples
        --------
        
        .. plot::
            :format: doctest
            :include-source: True

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> 
            >>> from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
            >>> from owcsimpy.geoobjects.bases.rectplane_py import RectPlane_py as RectPlane
            >>> from owcsimpy.geoutils.draw import draw
            >>> 
            >>> # Prepare a canvas
            >>> fig,ax = draw(xlim=[-2,2],ylim=[-2,2],zlim=[-2,2])
            >>> v = Vector(np.array([1,np.deg2rad(90),np.deg2rad(90)]))
            >>> plane = RectPlane(normalVect=v,ctrPoint=np.zeros(3),dimensions=[2,1])
            >>> fig,ax = draw(figure=fig,axes=ax,planes=plane,facecolors='black',colors='black')
            >>> plane = plane.rotate(np.deg2rad(45),plane.normalVect)
            >>> fig,ax = draw(figure=fig,axes=ax,planes=plane,facecolors='red',colors='red')

        """
        if not isinstance(refVector,Vector):
            # If it is not the instance of Vector, assume it contains the cartesian 
            # coordinates
            refVector = Vector(refVector,which='cartesian')

        # Rodrigues' rotation
        Rd = getRodriguesMtx(theta,refVector())
        normalVect = (Rd@self.normalVect.reshape(3,1)).reshape(-1)
        verts = np.array([(Rd@vert.reshape(3,1)).reshape(-1) for vert in self.verts])

        newPlane = RectPlane_py(normalVect,verts=verts,
            identity=self.identity,m=self.m,FoV=self.FoV,reflectivity=self.reflectivity)

        newPlane.identity = self.identity
        return newPlane

    def translate(self,newCtrPoint):
        """ Translate the plane to a new center point.

        Parameters
        ----------
        newCtrPoint: ndarray(3,)

        Returns
        -------
        RectPlane_py

        Examples
        --------
        .. plot::
            :format: doctest
            :include-source: True

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> 
            >>> from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
            >>> from owcsimpy.geoobjects.bases.rectplane_py import RectPlane_py as RectPlane
            >>> from owcsimpy.geoutils.draw import draw
            >>> 
            >>> # Prepare a canvas
            >>> fig,ax = draw(xlim=[-2,2],ylim=[-2,2],zlim=[-2,2])
            >>> v = Vector(np.array([1,np.deg2rad(90),np.deg2rad(90)]))
            >>> plane = RectPlane(normalVect=v,ctrPoint=np.zeros(3),dimensions=[2,1])
            >>> fig,ax = draw(figure=fig,axes=ax,planes=plane,facecolors='black',colors='black')
            >>> plane = plane.translate(np.ones(3))
            >>> fig,ax = draw(figure=fig,axes=ax,planes=plane,facecolors='red',colors='red')

        """

        verts = self.verts+newCtrPoint
        newPlane =  RectPlane_py(self.normalVect,verts=verts,
            identity=self.identity,m=self.m,FoV=self.FoV,reflectivity=self.reflectivity)

        newPlane.identity = self.identity
        return newPlane

    def getSimplePlane(self):
        """
        Extract necessary information as ndarray without methods.

        A simpleplane will be simply as a tuple of 
        normalVect, ctrPoint, verts and area. 
        
        Returns
        -------
        (revised)
        tuple:
            (normalVect: ndarray(3,)
                ctrPoint: ndarray(3,) 
                verts: ndarray(4,3)
                area: double
                m: double
                FoV: double
                reflectivity: double)
        
        (previous)
        tuple:
            (normalVect: ndarray(3,)
            ctrPoint: ndarray(3,) 
            verts: ndarray(4,3)
            area: double)

        Notes
        -----
        The order of the output matters.

        """

        # return (self.normalVect, 
        #     self.ctrPoint, 
        #     self.verts,
        #     self.area)

        return (self.normalVect, 
            self.ctrPoint, 
            self.verts,
            self.area,
            self.m,
            self.FoV,
            self.reflectivity)

    def getPartition(self,Ps=[1,1],delta=None):
        """Get a list of partitioned planes.
        
        Parameters
        ----------
        Ps: list
            List of number of partition of each side. 
        delta: list
            Define the partition based on partition lengths

        Returns
        -------
        list:
            A list of partitioned planes. Each plane is an instant of 
            RectPlane_py.
        
        Notes
        -----
        .. code-block:: text

                  0                             1
                   _________,_________,_________,
                   |        |         |         |
                   |        |         |         |
                   |________|_________|_________|
                   |        |         |         |
            l_03   |        |         |         | l_12
                   |________|_________|_________|
                   |        |         |         |
                   |        |         |         |
                  \./_______|_________|________\./
                 3                              2
        


        Steps:
            1. Get parametric lines l_03 and l_12
            2. Get row-major parametric lines whose initial 
               and end points are obtained from l_03 and l_12
            3. Get matrix of vertices
            4. Create a subplane from the matrix

        For future works, it would be more efficient if we have private attributes 
        on the spherical and Rodrigues' rotations, then we can partition the original 
        position rectangular plane and do rotation and translation afterwards.

        Examples
        --------

        .. plot::
            :format: doctest
            :include-source: True

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> 
            >>> from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
            >>> from owcsimpy.geoobjects.bases.rectplane_py import RectPlane_py as RectPlane
            >>> from owcsimpy.geoutils.draw import draw
            >>> 
            >>> # Original position
            >>> polar,az,Rod=90,-90,30 # polar, azimuth, Rodrigues
            >>> translation=0.5
            >>> v = Vector(np.array([1,np.deg2rad(polar),np.deg2rad(az)]))
            >>> ctrPoint=np.array(3*[translation])
            >>> plane = RectPlane(normalVect=v,ctrPoint=ctrPoint,
            >>>     RodriguesAngle=np.deg2rad(Rod),dimensions=[2,1])
            >>> 
            >>> subplanes = plane.getPartition(2)
            >>> fig,ax = draw(planes=subplanes,xlim=[-2,2],ylim=[-2,2],zlim=[-2,2])
            >>> # For reference
            >>> draw(figure=fig,axes=ax,planes=plane,facecolors='black',colors='black')
            >>> 
            >>> plt.show()
        
        """

        if delta == None:
            # Casting to list if it is a constant
            Ps = [Ps] if not isinstance(Ps,list) else Ps

            assert isinstance(Ps,list)
            assert 0 < len(Ps) <= 2

            if len(Ps) == 1:
                #  If the length is one then replicate
                Ps = 2*Ps 

            # Unpack
            Px,Py = Ps
        else:
            # Casting to a list
            delta = [delta] if not isinstance(delta,list) else delta
            
            assert isinstance(delta,list)
            assert 0 < len(delta) <= 2

            delta = 2*delta if len(delta) == 1 else delta

            # Unpack
            delta_x,delta_y = delta

            assert 0 < delta_x <= self.L and 0 < delta_y <= self.W

            Px = int(self.L/delta_x)
            Py = int(self.W/delta_y)

        # Get the parametric lines
        l_12 = Line(self.verts[1],self.verts[2])
        l_03 = Line(self.verts[0],self.verts[3])

        tx = np.linspace(0,1,Px+1)
        ty = np.linspace(0,1,Py+1)

        # Get row-major lines
        lrow = [Line(l_03.getPoint(typ),l_12.getPoint(typ)) for typ in ty]

        # Get matrix of vertices
        tmp = [list(map(lambda t: l.getPoint(t), tx)) for l in lrow]
        M = np.array([[tmp[idr][idc].tolist() for idc in range(Px+1)] for idr in range(Py+1)])

        # Get list of indexex of vertices of each partitioned plane
        listIdxVerts = [[(idr,idc),(idr,idc+1),(idr+1,idc+1),(idr+1,idc)] for idr in range(Py) for idc in range(Px)]
        

        # FIXME: transform this into a list comprehension would 
        # make it faster in exchange for readibility
        subplanes = []
        for idxNum,idxVerts in enumerate(listIdxVerts):
            verts = np.array([M[idxVerts[idv][0],idxVerts[idv][1],:].tolist() for idv in range(4)])
            normalVect = Vector(self.normalVect,which='cartesian')
            subplanes.append(RectPlane_py(normalVect,verts=verts,
                identity=self.identity+[idxNum],m=self.m,FoV=self.FoV,reflectivity=self.reflectivity))
        
        return subplanes












        
       