# import
from owcsimpy.geoobjects.models.typicalsimpleofficeenv_py import TypicalSimpleOfficeEnv_py as TypicalOfficeEnv
from owcsimpy.cir.freqdomaincir import FreqDomainCIR
from owcsimpy.cir.timedomaincir import TimeDomainCIR
from owcsimpy.cir.spheremodelcir import SphereModelCIR
from scipy.constants import speed_of_light
from owcsimpy.cir.cirutils import calcCIRFreqDom_fromMtx
import numpy as np

import matplotlib.pyplot as plt

import mpl_toolkits.mplot3d.art3d as art3d

from JSAnimation.IPython_display import display_animation
from matplotlib import animation

from owcsimpy.geoutils.mobility import getArrayMobility 

roomLength,roomWidth,roomHeight=4,3,3

position, furnitureConfig, humanLoc, humanDirection, activity = \
    'sitting','cornerYFacingWall',[roomLength/2,roomWidth/2],np.deg2rad(45),'reading'

office = TypicalOfficeEnv(
    roomDim=[roomLength,roomWidth,roomHeight],
    humanLoc=humanLoc,humanDirection=humanDirection,
    mode='rx',
    activity=activity,
    position=position,
    furnitureConfig=furnitureConfig
)

# activity = 'standing'
# minSamples = 2000; # 1 sample = 1 ms for activity = 'standing'
# traj = getArrayMobility(activity = 'standing',minSamples=minSamples)

# arrayLoc = traj.get('arrayLoc')
# arrayVel = traj.get('arrayVel')
# arrayTime = traj.get('arrayTime')
# arrayDir = traj.get('arrayDir')
# arrayTheta0 = traj.get('arrayTheta0')
# arrayOmega0 = traj.get('arrayOmega0')

# # farrayLoc, arrayVel, and arrayDir are empty lists if activity = 'standing'
# arrayLoc = np.ones((minSamples,2))*office.human.ctrPoint[0:2] # static
# arrayVel = np.zeros((minSamples,)) # zero velocity
# arrayDir = np.ones((minSamples,))*humanDirection # static direction

fig,ax = office.draw();
fig.set_size_inches(5, 4)

plt.show()

# # Since the time resolution is 1 ms, we need to slice it so that the animation won't take too long
# numSlice = 10 # in the order of (1 -> 1 ms, 10 -> 10 ms, etc.) 

# arrayLocSliced = arrayLoc[0::numSlice]
# arrayVelSliced = arrayVel[0::numSlice]
# arrayTimeSliced = arrayTime[0::numSlice]
# arrayDirSliced = arrayDir[0::numSlice]
# arrayTheta0Sliced = arrayTheta0[0::numSlice]
# arrayOmega0Sliced = arrayOmega0[0::numSlice]


# fig,ax,anim = office.getAnimation(
#     arrayLocSliced,arrayDirSliced,arrayTheta0Sliced,arrayOmega0Sliced,arrayTimeSliced,interval=1)
# fig.set_size_inches(5, 4)

# plt.show()

# # Case 1
# position, furnitureConfig, humanLoc, humanDirection, activity = \
#     'standing','cornerYFacingWall',[0,0],np.deg2rad(0),'calling'

# office1 = TypicalOfficeEnv(
#     roomDim=[roomLength,roomWidth,roomHeight],
#     humanLoc=humanLoc,humanDirection=humanDirection,
#     mode='rx',
#     activity=activity,
#     position=position,
#     furnitureConfig=furnitureConfig
# )

# Nframes = 60

# arrayLoc = np.transpose([np.linspace(0,2,Nframes),np.zeros((Nframes,))])
# arrayDir = np.zeros((Nframes,))
# arrayPolar = np.linspace(0,np.pi/4,Nframes)

# fig,ax,anim = office1.getAnimation(arrayLoc,arrayDir,arrayPolar)

# # display_animation(anim, default_mode='once')
# plt.show()
