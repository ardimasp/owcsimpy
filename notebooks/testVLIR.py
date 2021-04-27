from owcsimpy.geoobjects.models.typicalsimpleofficeenv_py import TypicalSimpleOfficeEnv_py as TypicalOfficeEnv
from owcsimpy.cir.freqdomaincir import FreqDomainCIR
from owcsimpy.cir.timedomaincir import TimeDomainCIR
from owcsimpy.cir.spheremodelcir import SphereModelCIR
from scipy.constants import speed_of_light
from owcsimpy.cir.cirutils import calcCIRFreqDom_fromMtx
import numpy as np

import matplotlib.pyplot as plt

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
FONT_SIZE = 24

plt.rc('font', size=FONT_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=FONT_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=FONT_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=FONT_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=FONT_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=FONT_SIZE)    # legend fontsize
plt.rc('figure', titlesize=FONT_SIZE)  # fontsize of the figure title

roomLength,roomWidth,roomHeight=4,3,3

# Case 1
position, furnitureConfig, humanLoc, humanDirection, activity = \
    'standing','cornerYFacingWall',[2.5,1.5],np.deg2rad(-45),'calling'

# position, furnitureConfig, humanLoc, humanDirection, activity = \
#     'sitting','cornerYFacingWall',[3.5,2.5],np.deg2rad(-135),'usbdongle'


# position, furnitureConfig, humanLoc, humanDirection, activity = \
#     'standing','cornerYFacingWall',[3.5,2.5],np.deg2rad(-135),'calling'


office1 = TypicalOfficeEnv(
    roomDim=[roomLength,roomWidth,roomHeight],
    humanLoc=humanLoc,humanDirection=humanDirection,
    mode='rx',
    activity=activity,
    position=position,
    furnitureConfig=furnitureConfig
)

fig,ax = office1.draw();
fig.set_size_inches(10, 8)
# ax.set_aspect(0.16)

ax.set_xticks([0,4])
ax.set_yticks([0,3])
ax.set_zticks([0,3])

plt.show()