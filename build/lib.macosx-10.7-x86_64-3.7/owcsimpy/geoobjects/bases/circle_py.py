import math
import numpy as np

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoutils.cutils import calcRodriguesMtx as getRodriguesMtx
from owcsimpy.geoutils.cutils import calcAngle

class Circle_py(object):
    """A circle in 3D plane.

    Parameters
    ----------
    normalVect: Vector_py
        Normal vector.
    ctrPoint: ndarray(3,)
        Center point.
    radius: float

    Attributes
    ----------
    normalVect: ndarray(3,)
    ctrPoint: ndarray(3,)
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
        >>> from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
        >>> from owcsimpy.geoobjects.bases.circle_py import Circle_py as Circle
        >>> from owcsimpy.geoutils.draw import draw
        >>> 
        >>> normalVect = Vector(np.array([1,np.deg2rad(30),np.deg2rad(45)]))
        >>> ctrPoint = 0.5*np.ones(3)
        >>> 
        >>> circle = Circle(normalVect,ctrPoint,radius=0.25)
        >>> 
        >>> fig,ax = draw(circles=circle,figsize=(5,6))
        >>> plt.show()

    """

    def __init__(self,normalVect,ctrPoint,radius):

        assert isinstance(normalVect,Vector)
        assert radius > 0
        assert ctrPoint.size == 3

        self.normalVect = normalVect()
        self.ctrPoint = ctrPoint
        self.radius = radius
        self.polar,self.azimuth = normalVect.spherical[1:]



