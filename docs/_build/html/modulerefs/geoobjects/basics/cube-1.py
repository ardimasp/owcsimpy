import matplotlib.pyplot as plt
import numpy as np

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.cube_py import Cube_py as Cube
from owcsimpy.geoutils.draw import draw

def genCube(polar,az,Rod,ctrPoint):
    
    cube = Cube(
        Vector(np.array([1,np.deg2rad(polar),np.deg2rad(az)])),
        ctrPoint = ctrPoint*np.ones(3),
        dimensions = [2,1,1],
        RodriguesAngle = np.deg2rad(Rod)
    )

    return cube


# Will draw 4 different canvases
fig,axs = draw(subplots=True,figsize=(14,6),nrows=1,ncols=4,xlim=[-2,2],ylim=[-2,2],zlim=[-2,2])

# Original position
polar,az,Rod=0,0,0 # polar, azimuth, Rodrigues
ctrPoint=0
cube = genCube(polar,az,Rod,ctrPoint)

fig,axs[0] = draw(figure=fig,axes=axs[0],cubes=cube,colors='blue',facecolors='blue')
axs[0].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,ctrPoint))

# Polar and azimuth
# Copy previous object as a reference (black)
fig,axs[1] = draw(figure=fig,axes=axs[1],cubes=cube,colors='black',facecolors='black')

polar,az,Rod=90,45,0 # polar, azimuth, Rodrigues
ctrPoint=0
cube = genCube(polar,az,Rod,ctrPoint)

fig,axs[1] = draw(figure=fig,axes=axs[1],cubes=cube,colors='blue',facecolors='blue')
axs[1].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,ctrPoint))

# Rodrigues
# Copy previous object as a reference (black)
fig,axs[2] = draw(figure=fig,axes=axs[2],cubes=cube,colors='black',facecolors='black')

polar,az,Rod=90,45,30 # polar, azimuth, Rodrigues
ctrPoint=0
cube = genCube(polar,az,Rod,ctrPoint)

fig,axs[2] = draw(figure=fig,axes=axs[2],cubes=cube,colors='blue',facecolors='blue')
axs[2].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,ctrPoint))

# Translation
# Copy previous object as a reference (black)
fig,axs[3] = draw(figure=fig,axes=axs[3],cubes=cube,colors='black',facecolors='black')

polar,az,Rod=90,45,30 # polar, azimuth, Rodrigues
ctrPoint=0.5
cube = genCube(polar,az,Rod,ctrPoint)

fig,axs[3] = draw(figure=fig,axes=axs[3],cubes=cube,colors='blue',facecolors='blue')
axs[3].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,ctrPoint))


plt.show()
