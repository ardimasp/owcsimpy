import matplotlib.pyplot as plt
import numpy as np

from owcsimpy.geoobjects.models.humancube_py import HumanCube_py as Human
from owcsimpy.geoutils.draw import draw

rho_keys = ['shirt','hair']
rho_vals = [0.3,0.6]
reflectivities = {rho_keys[i]:rho_vals[i] for i in range(2)}

persons = []
persons.append(Human(
    direction=np.deg2rad(45),
    loc=np.array([2,2]),
    dimensions=[0.3,0.5,1.7],
    identity=2,
    reflectivities=reflectivities))
persons.append(Human(
    direction=np.deg2rad(180),
    loc=np.array([4,1]),
    dimensions=[0.3,0.5,1.7],
    identity=2,
    reflectivities=reflectivities))
draw(models3d=persons,xlim=[0,5],ylim=[0,4],zlim=[0,3],
    facecolors=['red','blue'],colors=['red','blue']);
# >>>
plt.show() 
