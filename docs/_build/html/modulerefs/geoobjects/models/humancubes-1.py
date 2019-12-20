import matplotlib.pyplot as plt
import numpy as np

from owcsimpy.geoobjects.models.humancubes_py import HumanCubes_py as Human
from owcsimpy.geoutils.draw import draw
person = Human(np.array([1,1]),np.deg2rad(30),reflectivities={'hair':0.7,'face':0.5,'shirt':0.3})

draw(models3d=person,figsize=(6,6),azim=-100,elev=25,xlim=[0,3],ylim=[0,3],zlim=[0,3]);

plt.show()
