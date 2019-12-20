import numpy as np
import matplotlib.pyplot as plt
from owcsimpy.geoutils.draw import draw
from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
r = 0.5; polar = np.deg2rad(45); azimuth = np.deg2rad(25)
v1 = Vector(np.array([r,polar,azimuth]),refPoint=np.array([0.25,0.25,0]),which='spherical')
# Normalize length and rotate
v2 = v1.normalize().rotate(theta=np.deg2rad(30),refVector=np.array([0,0,1]))
# Translate
v3 = v1.translate(v1.refPoint+np.array([0,0,0.5]))
# Draw
fig,ax = draw(vectors=[v1,v2,v3],colors=['black','red','blue'],figsize=(5,6))
plt.show()
import matplotlib.pyplot as plt
