import numpy as np
import matplotlib.pyplot as plt
from owcsimpy.geoutils.draw import draw
from owcsimpy.geoobjects.bases.paramline_py import ParamLine_py as Line
# Generate a line l
l = Line(np.array([0.5,0.5,0.5]),np.ones(3))
# Draw
fig,ax = draw(lines=l,figsize=(5,6))
# Get a point at t = 0.25
P = l.getPoint(0.25)
print("Point at t=0.25 is {}".format(P))
# Draw
x,y,z = P
ax.scatter(x,y,z)
plt.show()  
