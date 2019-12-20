import functools
import warnings

from matplotlib import docstring
import matplotlib.artist as martist
from matplotlib.axes._axes import Axes
from matplotlib.gridspec import GridSpec, SubplotSpec
import matplotlib._layoutbox as layoutbox

import inspect

from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.collections as mcoll
import numpy as np

class SubplotBase(object):
    """
    Base class for subplots, which are :class:`Axes` instances with
    additional methods to facilitate generating and manipulating a set
    of :class:`Axes` within a figure.
    """

    def __init__(self, fig, *args, **kwargs):
        """
        *fig* is a :class:`matplotlib.figure.Figure` instance.

        *args* is the tuple (*numRows*, *numCols*, *plotNum*), where
        the array of subplots in the figure has dimensions *numRows*,
        *numCols*, and where *plotNum* is the number of the subplot
        being created.  *plotNum* starts at 1 in the upper left
        corner and increases to the right.

        If *numRows* <= *numCols* <= *plotNum* < 10, *args* can be the
        decimal integer *numRows* * 100 + *numCols* * 10 + *plotNum*.
        """

        # print(inspect.getmodule(fig))
        self.figure = fig

        if len(args) == 1:
            if isinstance(args[0], SubplotSpec):
                self._subplotspec = args[0]
            else:
                try:
                    s = str(int(args[0]))
                    rows, cols, num = map(int, s)
                except ValueError:
                    raise ValueError('Single argument to subplot must be '
                        'a 3-digit integer')
                self._subplotspec = GridSpec(rows, cols,
                                             figure=self.figure)[num - 1]
                # num - 1 for converting from MATLAB to python indexing
        elif len(args) == 3:
            rows, cols, num = args
            rows = int(rows)
            cols = int(cols)
            if isinstance(num, tuple) and len(num) == 2:
                num = [int(n) for n in num]
                self._subplotspec = GridSpec(
                        rows, cols,
                        figure=self.figure)[(num[0] - 1):num[1]]
            else:
                if num < 1 or num > rows*cols:
                    raise ValueError(
                        ("num must be 1 <= num <= {maxn}, not {num}"
                        ).format(maxn=rows*cols, num=num))
                self._subplotspec = GridSpec(
                        rows, cols, figure=self.figure)[int(num) - 1]
                # num - 1 for converting from MATLAB to python indexing
        else:
            raise ValueError('Illegal argument(s) to subplot: %s' % (args,))

        self.update_params()

        # _axes_class is set in the subplot_class_factory
        self._axes_class.__init__(self, fig, self.figbox, **kwargs)
        # add a layout box to this, for both the full axis, and the poss
        # of the axis.  We need both because the axes may become smaller
        # due to parasitic axes and hence no longer fill the subplotspec.
        if self._subplotspec._layoutbox is None:
            self._layoutbox = None
            self._poslayoutbox = None
        else:
            name = self._subplotspec._layoutbox.name + '.ax'
            name = name + layoutbox.seq_id()
            self._layoutbox = layoutbox.LayoutBox(
                    parent=self._subplotspec._layoutbox,
                    name=name,
                    artist=self)
            self._poslayoutbox = layoutbox.LayoutBox(
                    parent=self._layoutbox,
                    name=self._layoutbox.name+'.pos',
                    pos=True, subplot=True, artist=self)

    def __reduce__(self):
        # get the first axes class which does not inherit from a subplotbase
        axes_class = next(
            c for c in type(self).__mro__
            if issubclass(c, Axes) and not issubclass(c, SubplotBase))
        return (_picklable_subplot_class_constructor,
                (axes_class,),
                self.__getstate__())

    def get_geometry(self):
        """get the subplot geometry, e.g., 2,2,3"""
        rows, cols, num1, num2 = self.get_subplotspec().get_geometry()
        return rows, cols, num1 + 1  # for compatibility

    # COVERAGE NOTE: Never used internally or from examples
    def change_geometry(self, numrows, numcols, num):
        """change subplot geometry, e.g., from 1,1,1 to 2,2,3"""
        self._subplotspec = GridSpec(numrows, numcols,
                                     figure=self.figure)[num - 1]
        self.update_params()
        self.set_position(self.figbox)

    def get_subplotspec(self):
        """get the SubplotSpec instance associated with the subplot"""
        return self._subplotspec

    def set_subplotspec(self, subplotspec):
        """set the SubplotSpec instance associated with the subplot"""
        self._subplotspec = subplotspec

    def get_gridspec(self):
        """get the GridSpec instance associated with the subplot"""
        return self._subplotspec.get_gridspec()

    def update_params(self):
        """update the subplot position from fig.subplotpars"""

        self.figbox, self.rowNum, self.colNum, self.numRows, self.numCols = \
            self.get_subplotspec().get_position(self.figure,
                                                return_all=True)

    def is_first_col(self):
        return self.colNum == 0

    def is_first_row(self):
        return self.rowNum == 0

    def is_last_row(self):
        return self.rowNum == self.numRows - 1

    def is_last_col(self):
        return self.colNum == self.numCols - 1

    # COVERAGE NOTE: Never used internally.
    def label_outer(self):
        """Only show "outer" labels and tick labels.

        x-labels are only kept for subplots on the last row; y-labels only for
        subplots on the first column.
        """
        lastrow = self.is_last_row()
        firstcol = self.is_first_col()
        if not lastrow:
            for label in self.get_xticklabels(which="both"):
                label.set_visible(False)
            self.get_xaxis().get_offset_text().set_visible(False)
            self.set_xlabel("")
        if not firstcol:
            for label in self.get_yticklabels(which="both"):
                label.set_visible(False)
            self.get_yaxis().get_offset_text().set_visible(False)
            self.set_ylabel("")

    def _make_twin_axes(self, *kl, **kwargs):
        """
        Make a twinx axes of self. This is used for twinx and twiny.
        """
        from matplotlib.projections import process_projection_requirements
        if 'sharex' in kwargs and 'sharey' in kwargs:
            # The following line is added in v2.2 to avoid breaking Seaborn,
            # which currently uses this internal API.
            if kwargs["sharex"] is not self and kwargs["sharey"] is not self:
                raise ValueError("Twinned Axes may share only one axis.")
        kl = (self.get_subplotspec(),) + kl
        projection_class, kwargs, key = process_projection_requirements(
            self.figure, *kl, **kwargs)

        ax2 = subplot_class_factory(projection_class)(self.figure,
                                                      *kl, **kwargs)
        self.figure.add_subplot(ax2)
        self.set_adjustable('datalim')
        ax2.set_adjustable('datalim')

        if self._layoutbox is not None and ax2._layoutbox is not None:
            # make the layout boxes be explicitly the same
            ax2._layoutbox.constrain_same(self._layoutbox)
            ax2._poslayoutbox.constrain_same(self._poslayoutbox)

        self._twinned_axes.join(self, ax2)
        return ax2
    
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
            
            # print('heads')
            # print(heads)
            # print('shafts')
            # print(shafts)
            """
            Rodrigues' rotation so that 
            we can project the heads to x-y plane properly
            """
            # Get the vector representation of the shaft
            # notice that the convention is reversed
            k = shafts[0,0,:]-shafts[0,1,:] 
            k = k/np.linalg.norm(k)
    
            # Cross-product matrix
            K = np.array([
                    [    0,-k[2], k[1]],
                    [ k[2],    0,-k[0]],
                    [-k[1], k[0],   0]])
            
            theta = np.pi/2
            R = np.eye(3)+ K*np.sin(theta) + (K@K*(1-np.cos(theta)))
            
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



# this here to support cartopy which was using a private part of the
# API to register their Axes subclasses.

# In 3.1 this should be changed to a dict subclass that warns on use
# In 3.3 to a dict subclass that raises a useful exception on use
# In 3.4 should be removed

# The slow timeline is to give cartopy enough time to get several
# release out before we break them.
_subplot_classes = {}


@functools.lru_cache(None)
def subplot_class_factory(axes_class=None,fig=None):
    """
    This makes a new class that inherits from `.SubplotBase` and the
    given axes_class (which is assumed to be a subclass of `.axes.Axes`).
    This is perhaps a little bit roundabout to make a new class on
    the fly like this, but it means that a new Subplot class does
    not have to be created for every type of Axes.
    """
    if axes_class is None:
        axes_class = Axes
    try:
        # Avoid creating two different instances of GeoAxesSubplot...
        # Only a temporary backcompat fix.  This should be removed in
        # 3.4
        return next(cls for cls in SubplotBase.__subclasses__()
                    if cls.__bases__ == (SubplotBase, axes_class))
    except StopIteration:
        return type("%sSubplot" % axes_class.__name__,
                    (SubplotBase, axes_class),
                    {'_axes_class': axes_class})


# This is provided for backward compatibility
Subplot = subplot_class_factory()


def _picklable_subplot_class_constructor(axes_class):
    """
    This stub class exists to return the appropriate subplot class when called
    with an axes class. This is purely to allow pickling of Axes and Subplots.
    """
    subplot_class = subplot_class_factory(axes_class)
    return subplot_class.__new__(subplot_class)


docstring.interpd.update(Axes=martist.kwdoc(Axes))
docstring.dedent_interpd(Axes.__init__)

docstring.interpd.update(Subplot=martist.kwdoc(Axes))
