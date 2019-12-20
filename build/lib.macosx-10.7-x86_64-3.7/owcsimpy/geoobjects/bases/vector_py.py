import numpy as np
import math 

from owcsimpy.geoutils.cutils import calcRodriguesMtx

class Vector_py(object):
    """ A 3D vector.

    Parameters
    ----------

    coord: ndarray(3,)
        coord can be either represents spherical coordinates, i.e., 
        np.array([radius,polar,azimuth]) or Cartesian coordinates, 
        i.e., np.array([radius,polar,azimuth]). The default method is
        the spherical one. Change the parameter **which** to change 
        the instantiation method.
    refPoint: ndarray(3,)
        refPoint is the reference point of the vector.
    which: {'spherical','cartesian'}
        **which** defines which representation that is used. 
        **which** is either 'spherical' (default) or 'cartesian'. 
    
    Attributes
    ----------
    spherical: ndarray(3,)
        Spherical coordinates of the vector.
    cartesian: ndarray(3,)
        Cartesian coordinates of the vector.
    refPoint: ndarray(3,)
        The reference point of the vector.

    Raises
    ------

    NotImplementedError: 
        When **which** is neither 'spherical' nor 'cartesian'.

    See Also
    --------

    Notes
    -----
    This class supports a callable method (with empty argument) that returns its 
    Cartesian coordinates. 
    
    Examples
    --------

    .. plot::
       :format: doctest
       :include-source: True

       >>> import numpy as np
       >>> import matplotlib.pyplot as plt
       >>> from owcsimpy.geoutils.draw import draw
       >>> from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
       >>> r = 0.5; polar = np.deg2rad(45); azimuth = np.deg2rad(25)
       >>> v1 = Vector(np.array([r,polar,azimuth]),refPoint=np.array([0.25,0.25,0]),which='spherical')
       >>> # Normalize length and rotate
       >>> v2 = v1.normalize().rotate(theta=np.deg2rad(30),refVector=np.array([0,0,1]))
       >>> # Translate
       >>> v3 = v1.translate(v1.refPoint+np.array([0,0,0.5]))
       >>> # Draw
       >>> fig,ax = draw(vectors=[v1,v2,v3],colors=['black','red','blue'],figsize=(5,6))
       >>> plt.show()
       >>> import matplotlib.pyplot as plt

    """
    
    def __init__(self,coord,refPoint=np.zeros(3),which='spherical'):
        
        # assert coord.size == 3 
        try: 
            assert coord.size == 3
        except:
            assert len(coord) == 3 
        
        if which.lower() == 'spherical':
            radius,polar,azimuth = coord
            x = radius*np.sin(polar)*np.cos(azimuth)
            y = radius*np.sin(polar)*np.sin(azimuth)
            z = radius*np.cos(polar)
        elif which.lower() == 'cartesian':
            x,y,z = coord 
            radius = np.linalg.norm(coord)
            if math.isclose(radius,0):
                polar,azimuth = 0,0
            else:
                polar = np.arccos((z)/radius)
                azimuth = np.arctan2(y,x)
        else:
            raise NotImplementedError("This class can only accept 'spherical' or 'cartesian'")
            
        self.spherical = np.array([radius,polar,azimuth])
        self.cartesian = np.array([x,y,z]) 
        self.refPoint = refPoint

    def __call__(self):
        """ 
        Callable class.

        Returns
        =======

        Cartesian coordinates.

        """
        return self.cartesian
            
    def normalize(self):
        """ 
        Normalize the vector into a unit length vector.
        
        Returns
        =======

        A new vector object.

        """
        newVector = Vector_py(np.array([1,self.spherical[1],self.spherical[2]]),
            refPoint=self.refPoint,which='spherical')
        return newVector

    def rotate(self,theta,refVector):
        """ 
        Rotate a matrix with theta rads w.r.t. refVector. 
        
        Parameters
        ==========
        theta: float
            The rotatton angle theta (in rads) follows the right-hand rule.
        refVector: ndarray(3,)
            An axis of rotation about which the vector rotates.

        Returns
        =======

        A new vector object.

        """
        # Rotation matrix
        R = calcRodriguesMtx(theta,refVector)
        newCoord = ((R@(self.cartesian).reshape(3,1)).reshape(-1))
 
        newVector = Vector_py(newCoord,
            refPoint=self.refPoint,which='cartesian')
        return newVector

    def translate(self,newRefPoint):
        """ 
        Shift the refPoint of the vector to newRefPoint. 
        
        Parameters
        ==========
        newRefPoint: ndarray(3,)
            The new reference point to which the vector will be translated.

        Returns
        =======

        A new vector object.

        """

        return Vector_py(self.cartesian,refPoint=newRefPoint,
            which='cartesian')
        

        