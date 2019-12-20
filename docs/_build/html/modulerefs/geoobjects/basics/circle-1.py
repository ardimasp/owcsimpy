import matplotlib.pyplot as plt
import numpy as np

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.circle_py import Circle_py as Circle
from owcsimpy.geoutils.draw import draw

normalVect = Vector(np.array([1,np.deg2rad(30),np.deg2rad(45)]))
ctrPoint = 0.5*np.ones(3)

circle = Circle(normalVect,ctrPoint,radius=0.25)

fig,ax = draw(circles=circle,figsize=(5,6))
plt.show()
