import numpy as np

class ParamLine_py(object):
    """ A 3D parametric line.

    Parameters
    ----------

    P0: ndarray(3,)
        A tail point.
    P1: ndarray(3,)
        A head point.

    Attributes
    ----------
    P0
    P1
    u: ndarray(3,)
        Show the direction of the parametric line

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from owcsimpy.geoutils.draw import draw
    >>> from owcsimpy.geoobjects.bases.paramline_py import ParamLine_py as Line
    >>> # Generate a line l
    >>> l = Line(np.array([0.5,0.5,0.5]),np.ones(3))
    >>> # Draw
    >>> fig,ax = draw(lines=l,figsize=(5,6))
    >>> # Get a point at t = 0.25
    >>> P = l.getPoint(0.25)
    >>> print("Point at t=0.25 is {}".format(P))
    Point at t=0.25 is [0.625 0.625 0.625]
    >>> # Draw
    >>> x,y,z = P
    >>> ax.scatter(x,y,z)
    >>> plt.show()  

    .. plot::
       :format: doctest

       >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> from owcsimpy.geoutils.draw import draw
        >>> from owcsimpy.geoobjects.bases.paramline_py import ParamLine_py as Line
        >>> # Generate a line l
        >>> l = Line(np.array([0.5,0.5,0.5]),np.ones(3))
        >>> # Draw
        >>> fig,ax = draw(lines=l,figsize=(5,6))
        >>> # Get a point at t = 0.25
        >>> P = l.getPoint(0.25)
        >>> print("Point at t=0.25 is {}".format(P))
        >>> # Draw
        >>> x,y,z = P
        >>> ax.scatter(x,y,z)
        >>> plt.show()  


    Notes
    -----
    A parametric line, :math:`l(t)`, is defined as:


    .. math:: l(t) = P_0 + \mathbf{u} t, \mathrm{where}\ \mathbf{u} = P_1-P_0.


    """


    def __init__(self,P0,P1):
        
        assert P0.size == 3 and P1.size == 3

        self.P0,self.P1 = P0,P1
        self.u = P1-P0

    def getPoint(self,t0):
        """ Get a point at t=t0.
        
        Returns
        -------
        ndarray(3,)
            Return P0+u t0

        """

        return self.P0+self.u*t0

    def isValid(self,P):
        """ Check wheter the point P is in line.

        Returns
        -------
        bool

        """ 

        up = P-self.P0 # u prime
        return True if np.allclose(up,np.zeros(3)) else np.allclose(
            up/np.linalg.norm(up),self.u/np.linalg.norm(self.u))

    def getParam(self,P):
        """ Get t of P. Or, find t s.t. P0+u t = P
        
        Returns
        -------
        float or None
            Return None if P is not in the line.

        """

        # Check validity of P
        if self.isValid(P):
            idxnotzero = np.where(self.u!=0)[0]
            return (P-self.P0)[idxnotzero]/self.u[idxnotzero]
        else:
            return None


