import matplotlib.pyplot as plt
import numpy as np

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.rectplane_py import RectPlane_py as RectPlane
from owcsimpy.geoutils.draw import draw

# Prepare a canvas
fig,ax = draw(xlim=[-2,2],ylim=[-2,2],zlim=[-2,2])
v = Vector(np.array([1,np.deg2rad(90),np.deg2rad(90)]))
plane = RectPlane(normalVect=v,ctrPoint=np.zeros(3),dimensions=[2,1])
fig,ax = draw(figure=fig,axes=ax,planes=plane,facecolors='black',colors='black')
plane = plane.translate(np.ones(3))
fig,ax = draw(figure=fig,axes=ax,planes=plane,facecolors='red',colors='red')
