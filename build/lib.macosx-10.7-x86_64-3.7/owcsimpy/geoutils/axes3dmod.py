# Axes3DMod.py
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.collections as mcoll
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
plt.style.use('seaborn-whitegrid')
# We can't expect that the users have latex by default
# plt.rc('text', usetex=True)
plt.rc('font', family='serif')


class Axes3DMod(Axes3D):
    """
    Inherited class from Axes3D
    This class is defined to enable animations on our shapes class.
    
    
    Modifications
    -------------
        1. Overwrite the `add_collection3d` method
            add the return values

    """
    def __init__(
            self, fig, rect=None, *args,
            azim=-60, elev=30, zscale=None, sharez=None, proj_type='persp',**kwargs):
        super().__init__(fig, rect=rect, *args, azim=azim, elev=elev, 
            zscale=zscale, sharez=sharez, proj_type=proj_type, **kwargs)


    def add_collection3d(self, col, zs=0, zdir='z'):
        '''
        Add a 3D collection object to the plot.

        2D collection types are converted to a 3D version by
        modifying the object and adding z coordinate information.

        Supported are:
            - PolyCollection
            - LineCollection
            - PatchCollection
        '''
        zvals = np.atleast_1d(zs)
        if len(zvals) > 0 :
            zsortval = min(zvals)
        else :
            zsortval = 0   # FIXME: Fairly arbitrary. Is there a better value?

        # FIXME: use issubclass() (although, then a 3D collection
        #       object would also pass.)  Maybe have a collection3d
        #       abstract class to test for and exclude?
        if type(col) is mcoll.PolyCollection:
            art3d.poly_collection_2d_to_3d(col, zs=zs, zdir=zdir)
            col.set_sort_zpos(zsortval)
        elif type(col) is mcoll.LineCollection:
            art3d.line_collection_2d_to_3d(col, zs=zs, zdir=zdir)
            col.set_sort_zpos(zsortval)
        elif type(col) is mcoll.PatchCollection:
            art3d.patch_collection_2d_to_3d(col, zs=zs, zdir=zdir)
            col.set_sort_zpos(zsortval)

        return super().add_collection(col)
    
    def get_arrow(self, tail, head,
               length=1, arrow_length_ratio=.3, normalize=False,nargout=1,
               **kwargs):
        
        def calc_arrow(uvw, angle=15):
            """
            To calculate the arrow head. uvw should be a unit vector.
            We normalize it here:
            """
            # get unit direction vector perpendicular to (u,v,w)
            norm = np.linalg.norm(uvw[:2])
            if norm > 0:
                x = uvw[1] / norm
                y = -uvw[0] / norm
            else:
                x, y = 0, 1

            # compute the two arrowhead direction unit vectors
            ra = np.radians(angle)
            c = np.cos(ra)
            s = np.sin(ra)

            # construct the rotation matrices
            Rpos = np.array([[c+(x**2)*(1-c),      x*y*(1-c),  y*s],
                             [     y*x*(1-c), c+(y**2)*(1-c), -x*s],
                             [          -y*s,            x*s,    c]])
            # opposite rotation negates all the sin terms
            Rneg = Rpos.copy()
            Rneg[[0,1,2,2],[2,2,0,1]] = -Rneg[[0,1,2,2],[2,2,0,1]]

            # multiply them to get the rotated vector
            return Rpos.dot(uvw), Rneg.dot(uvw)
        
        assert tail.size == 3 and head.size == 3
        
        shaft_dt = np.linspace(0, length, num=2)
        arrow_dt = shaft_dt * arrow_length_ratio
        shaft_dt -= length

        XYZ = tail.reshape(1,3)
        UVW = head.reshape(1,3).astype(float)
        
        # Normalize rows of UVW
        norm = np.linalg.norm(UVW, axis=1)

        # If any row of UVW is all zeros, don't make a quiver for it
        mask = norm > 0
        XYZ = XYZ[mask]
        if normalize:
            UVW = UVW[mask] / norm[mask].reshape((-1, 1))
        else:
            UVW = UVW[mask]
        
        if len(XYZ) > 0:
            # compute the shaft lines all at once with an outer product
            shafts = (XYZ - np.multiply.outer(shaft_dt, UVW)).swapaxes(0, 1)
            # compute head direction vectors, n heads by 2 sides by 3 dimensions
            head_dirs = np.array([calc_arrow(d) for d in UVW])
            # compute all head lines at once, starting from where the shaft ends
            heads = shafts[:, :1] - np.multiply.outer(arrow_dt, head_dirs)
            # stack left and right head lines together
            heads.shape = (len(arrow_dt), -1, 3)
            # transpose to get a list of lines
            heads = heads.swapaxes(0, 1)
            
            """
            Rodrigues' rotation so that 
            we can project the heads to x-y plane properly
            """
            # Get the vector representation of the shaft
            # notice that the convention is reversed
            k = shafts[0,0,:]-shafts[0,1,:] 
            if not np.allclose(np.linalg.norm(k),0):
                k = k/np.linalg.norm(k)
        
                # Cross-product matrix
                K = np.array([
                        [    0,-k[2], k[1]],
                        [ k[2],    0,-k[0]],
                        [-k[1], k[0],   0]])
                
                theta = np.pi/2
                R = np.eye(3)+ K*np.sin(theta) + (K@K*(1-np.cos(theta)))
            else:
                R = np.eye(3)
            
            
            for idx in range(2):
                # Shift to origin and after rotation shift it back to the original
                # reference point.
                heads[idx,:,:] = ((heads[idx,:,:]-shafts[0,1,:])@R.T)+shafts[0,1,:] 
            
            lines = [*shafts, *heads]
        else:
            lines = []
            
        arrow = art3d.Line3DCollection(lines, **kwargs)
        
        if nargout == 1:
            return arrow
        elif nargout == 2:
            return arrow,XYZ
           
    def draw_arrow(self, tail, head,
               length=1, arrow_length_ratio=.3, normalize=False,
               **kwargs):
        """
        Taken from super().quiver() and simplified.
        
        Plot a 3D field of arrows.

        call signatures::

            arrow(tail, head, **kwargs)

        Arguments:

            tail:
                ndarray(3,) which contains position of the tail point in 
                x, y, z coordinates 

            head:
                ndarray(3,) which contains position of the tail point in 
                x, y, z coordinates

        Keyword arguments:

            *length*: [1.0 | float]
                The length of each quiver, default to 1.0, the unit is
                the same with the axes

            *arrow_length_ratio*: [0.3 | float]
                The ratio of the arrow head with respect to the quiver,
                default to 0.3

            *normalize*: bool
                When True, all of the arrows will be the same length. This
                defaults to False, where the arrows will be different lengths
                depending on the values of u,v,w.

        Any additional keyword arguments are delegated to
        :class:`~matplotlib.collections.LineCollection`

        """
        
        arrow,XYZ = self.get_arrow(tail, head,
               length=length, arrow_length_ratio=.3, normalize=False, nargout=2,
               **kwargs)

        # had_data = self.has_data()
        # self.auto_scale_xyz(XYZ[:, 0], XYZ[:, 1], XYZ[:, 2], had_data)
        
        self.add_collection(arrow)

        return arrow




