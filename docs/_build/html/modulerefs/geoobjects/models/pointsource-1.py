import matplotlib.pyplot as plt
import numpy as np

from owcsimpy.geoobjects.models.pointsource_py import PointSource_py as PointSource
from owcsimpy.geoutils.draw import draw
# >>>
led = PointSource(np.pi,0,np.array([2.5,2,3]))
# >>>
draw(vectors=led,xlim=[0,5],ylim=[0,4],zlim=[0,3])
# >>>
plt.show()
