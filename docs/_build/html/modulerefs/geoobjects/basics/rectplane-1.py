import matplotlib.pyplot as plt
import numpy as np

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.rectplane_py import RectPlane_py as RectPlane
from owcsimpy.geoutils.draw import draw

def genPlane(polar,az,Rod,translation):
    
    v = Vector(np.array([1,np.deg2rad(polar),np.deg2rad(az)]))
    ctrPoint=np.array(3*[translation])
    plane = RectPlane(normalVect=v,ctrPoint=ctrPoint,
        RodriguesAngle=np.deg2rad(Rod),dimensions=[2,1])

    return plane


# Will draw 4 different canvases
fig,axs = draw(subplots=True,figsize=(14,6),nrows=1,ncols=4,xlim=[-2,2],ylim=[-2,2],zlim=[-2,2])

# Original position
polar,az,Rod=0,0,0 # polar, azimuth, Rodrigues
translation=0
plane = genPlane(polar,az,Rod,translation)

fig,axs[0] = draw(figure=fig,axes=axs[0],planes=plane,colors='blue',facecolors='blue')
axs[0].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,translation))

# Polar and azimuth
# Copy previous object as a reference (black)
fig,axs[1] = draw(figure=fig,axes=axs[1],planes=plane,colors='black',facecolors='black')

polar,az,Rod=90,45,0 # polar, azimuth, Rodrigues
translation=0
plane = genPlane(polar,az,Rod,translation)

fig,axs[1] = draw(figure=fig,axes=axs[1],planes=plane,colors='blue',facecolors='blue')
axs[1].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,translation))

# Rodrigues
# Copy previous object as a reference (black)
fig,axs[2] = draw(figure=fig,axes=axs[2],planes=plane,colors='black',facecolors='black')

polar,az,Rod=90,45,30 # polar, azimuth, Rodrigues
translation=0
plane = genPlane(polar,az,Rod,translation)

fig,axs[2] = draw(figure=fig,axes=axs[2],planes=plane,colors='blue',facecolors='blue')
axs[2].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,translation))

# Translation
# Copy previous object as a reference (black)
fig,axs[3] = draw(figure=fig,axes=axs[3],planes=plane,colors='black',facecolors='black')

polar,az,Rod=90,45,30 # polar, azimuth, Rodrigues
translation=0.5
plane = genPlane(polar,az,Rod,translation)

fig,axs[3] = draw(figure=fig,axes=axs[3],planes=plane,colors='blue',facecolors='blue')
axs[3].set_title("angles=({},{},{}), x=y=z={}".format(polar,az,Rod,translation))


plt.show()
