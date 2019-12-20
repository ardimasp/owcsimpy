import math
import itertools
import numpy as np

from owcsimpy.geoobjects.bases.circle_py import Circle_py as Circle
from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector

class BareDetector_py(Circle):
    """ A bare detector model.

    A bare detector is modeled by a tranparent circle in 3D space. 
    It means that the detector does block any path. A circle is 
    used only for illustration. 

    Parameters
    ----------
    polar: float
        Polar angle of the normal vector in rads
    azimuth: float
        Azimuth angle in rads
    loc: ndarray(3,)
        Location of the point source.
    area: float
        Detector area in m^2.

    Attributes
    ----------
    normalVect: ndarray(3,)
    loc: ndarray(3,)
    radius: float
    polar: float
        Polar angle of the normal vector
    azimuth: float
        Azimuth angle of the normal vector

    Examples
    --------
    .. plot::
            :format: doctest
            :include-source: True

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> 
            >>> from owcsimpy.geoobjects.models.baredetector_py import BareDetector_py as BareDetector
            >>> from owcsimpy.geoutils.draw import draw
            >>>
            >>> pd = BareDetector(np.deg2rad(45),np.deg2rad(30),np.array([1,2,0]),area=1e-4)
            >>>
            >>> draw(circles=pd,scales=5e3,xlim=[0,5],ylim=[0,4],zlim=[0,3])
            >>>
            >>> plt.show()

    """

    def __init__(self,polar,azimuth,loc,area,FoV=np.pi/2):

        normalVect = Vector(np.array([1,polar,azimuth]),which='spherical')
        radius = np.sqrt(area/np.pi)
        super().__init__(normalVect,loc,radius)

        self.loc = loc
        self.area = area
        self.FoV = FoV

    def getSimpleBareDetector(self):
        """ Get a simple bare detector.

        Returns
        -------
        tuple: 
            (normalVect: ndarray(3,),ctrPoint: ndarray(3,), area: float, FoV: float)

        Notes
        -----
        The order of the output matters.
        
        """

        return (self.normalVect,self.ctrPoint,self.area,self.FoV)


