import numpy as np
import matplotlib.pyplot as plt

from owcsimpy.geoobjects.bases.rectplane_py import RectPlane_py as RectPlane
from owcsimpy.geoutils.draw import draw
from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoutils.cutils import checkBlockage

# Planes
planes = []

normalVect = Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)]),which='spherical')
ctrPoint = np.array([0,2,1])
planes.append(RectPlane(normalVect,ctrPoint=ctrPoint,
    RodriguesAngle=np.deg2rad(0),dimensions=[1.7,1]))

normalVect = Vector(np.array([1,np.deg2rad(90),np.deg2rad(180)]),which='spherical')
ctrPoint = np.array([4,2,2])
planes.append(RectPlane(normalVect,ctrPoint=ctrPoint,
    RodriguesAngle=np.deg2rad(0),dimensions=[1.7,1]))

# This should not block
normalVect = Vector(np.array([1,np.deg2rad(115),np.deg2rad(222)]),which='spherical')
ctrPoint = np.array([2,4,1])
planes.append(RectPlane(normalVect,ctrPoint=ctrPoint,RodriguesAngle=np.deg2rad(0),dimensions=[1.7,1]))

# This should block
normalVect = Vector(np.array([1,np.deg2rad(115),np.deg2rad(222)]),which='spherical')
ctrPoint = np.array([2,2,1])
planes.append(RectPlane(normalVect,ctrPoint=ctrPoint,RodriguesAngle=np.deg2rad(0),dimensions=[1.7,1]))

# This should not block
normalVect = Vector(np.array([1,np.deg2rad(90),np.deg2rad(180)]),which='spherical')
ctrPoint = np.array([5,2,3])
planes.append(RectPlane(normalVect,ctrPoint=ctrPoint,RodriguesAngle=np.deg2rad(0),dimensions=[1.7,1]))

# Prepare canvases
fig,axs = draw(subplots=True,figsize=(14,6),nrows=1,ncols=3,xlim=[0,5],ylim=[0,5],zlim=[0,3],
    azim=-32,elev=25)

for idx in range(2,5):
    colors = ['red','red','black']
    fig,axs[idx-2] = draw(figure=fig,axes=axs[idx-2],planes=[planes[0],planes[1],planes[idx]],
        colors=colors,facecolors=colors)

    ray = np.append(planes[0].ctrPoint,planes[1].ctrPoint).reshape(2,-1)

    # Draw ray
    axs[idx-2].plot(ray[:,0],ray[:,1],ray[:,2],color='k');

    # Get simple planes from each plane
    plane1 = planes[0].getSimplePlane()
    plane2 = planes[1].getSimplePlane()
    planeB = planes[idx].getSimplePlane()

    # Check whether the blockage occurs
    isBlocked,intersectingPoint = checkBlockage(plane1,plane2,planeB,nargout=2)

    if intersectingPoint is not None:
        xi,yi,zi = intersectingPoint
        axs[idx-2].scatter(xi,yi,zi)
        
    axs[idx-2].set_title("isBlocked: {}".format(isBlocked))

plt.show()
