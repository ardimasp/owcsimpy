import matplotlib.pyplot as plt
import numpy as np

from owcsimpy.geoobjects.models.baredetector_py import BareDetector_py as BareDetector
from owcsimpy.geoutils.draw import draw
# >>>
pd = BareDetector(np.deg2rad(45),np.deg2rad(30),np.array([1,2,0]),area=1e-4)
# >>>
draw(circles=pd,scales=5e3,xlim=[0,5],ylim=[0,4],zlim=[0,3])
# >>>
plt.show()
