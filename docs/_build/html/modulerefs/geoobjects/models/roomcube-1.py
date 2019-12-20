import matplotlib.pyplot as plt
import numpy as np

from owcsimpy.geoobjects.models.roomcube_py import RoomCube_py as Room
from owcsimpy.geoutils.draw import draw

rho_keys = ['b','t','s','n','e','w']
rho_vals = [0.1,0.5,0.5,0.5,0.5,0.5]
reflectivities = {rho_keys[i]:rho_vals[i] for i in range(6)}

fig,axs = draw(subplots=True,nrows=1,ncols=2,figsize=(14,8),
    xlim=[0,5],ylim=[0,4],zlim=[0,3]);

room = Room(dimensions=[5,4,3],identity=1,reflectivities=reflectivities)

# Partition with the delta of 1 m
subplanes = room.getPartition(delta=1)

fig,axs[0]=draw(figure=fig,axes=axs[0],models3d=room);
fig,axs[1]=draw(figure=fig,axes=axs[1],planes=subplanes,xlim=[0,5],ylim=[0,4],zlim=[0,3],
     alphas=[0.5*plane.reflectivity for plane in subplanes],
    lengths=0.3);
# >>>
plt.show()
