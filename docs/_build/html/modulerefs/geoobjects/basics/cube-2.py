import matplotlib.pyplot as plt
import numpy as np

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.cube_py import Cube_py as Cube
from owcsimpy.geoutils.draw import draw
cube = Cube(
    Vector(np.array([1,np.deg2rad(90),np.deg2rad(90)])),
    ctrPoint = np.array([0.5,0.5,0.5]),
    dimensions = [2,1,1],
    RodriguesAngle = np.deg2rad(30)
)
planes = cube.getPartition(delta=0.5)
fig,ax = draw(planes=planes,alphas=0.2,xlim=[-2,2],ylim=[-2,2],zlim=[-2,2])
plt.show()
