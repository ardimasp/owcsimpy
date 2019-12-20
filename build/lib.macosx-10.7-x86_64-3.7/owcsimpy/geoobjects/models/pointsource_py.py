import math
import itertools
import numpy as np

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector

class PointSource_py(Vector):
    """ A point source models.

    Mainly used for modeling an LED.

    HumanCube_py is inherited from 
    :class:`~owcsimpy.geoobjects.bases.vector_py.Vector_py`

    See Also
    --------
    :class:`~owcsimpy.geoobjects.bases.vector_py.Vector_py`

    Parameters
    ----------

    polar: float
        Polar angle of the normal vector in rads
    azimuth: float
        Azimuth angle in rads
    loc: ndarray(3,)
        Location of the point source.

    Attributes
    ----------
    loc: ndarray(3,)
        Location of the point source

    Examples
    --------

    .. plot:: 
            :format: doctest
            :include-source: True

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> 
            >>> from owcsimpy.geoobjects.models.pointsource_py import PointSource_py as PointSource
            >>> from owcsimpy.geoutils.draw import draw
            >>>
            >>> led = PointSource(np.pi,0,np.array([2.5,2,3]))
            >>>
            >>> draw(vectors=led,xlim=[0,5],ylim=[0,4],zlim=[0,3])
            >>>
            >>> plt.show()


    """

    def __init__(self,polar,azimuth,loc,m=1):

        coord = np.array([1,polar,azimuth])
        super().__init__(coord=coord,refPoint=loc,which='spherical')

        self.loc = loc
        self.m = m

    def getSimplePointSource(self):
        """ Get a simple point source.

        Returns
        -------
        tuple: 
            (normalVect: ndarray(3,),ctrPoint: ndarray(3,), m: float)

        Notes
        -----
        The order of the output matters.
        
        """

        return (self.cartesian,self.loc,self.m)

