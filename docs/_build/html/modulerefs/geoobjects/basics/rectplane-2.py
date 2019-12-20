import matplotlib.pyplot as plt
import numpy as np

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.rectplane_py import RectPlane_py as RectPlane
from owcsimpy.geoutils.draw import draw

# Original position
polar,az,Rod=90,-90,30 # polar, azimuth, Rodrigues
translation=0.5
v = Vector(np.array([1,np.deg2rad(polar),np.deg2rad(az)]))
ctrPoint=np.array(3*[translation])
plane = RectPlane(normalVect=v,ctrPoint=ctrPoint,
    RodriguesAngle=np.deg2rad(Rod),dimensions=[2,1])

subplanes = plane.getPartition(2)
fig,ax = draw(planes=subplanes,xlim=[-2,2],ylim=[-2,2],zlim=[-2,2])
# For reference
draw(figure=fig,axes=ax,planes=plane,facecolors='black',colors='black')

plt.show()
