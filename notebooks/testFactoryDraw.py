import numpy as np
import matplotlib.pyplot as plt

from owcsimpy.geoutils.draw import draw
from owcsimpy.geoobjects.models.humanwithactivity_py import HumanWithActivity_py as Human

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.paramline_py import ParamLine_py as Line
from owcsimpy.geoobjects.bases.circle_py import Circle_py as Circle
from owcsimpy.geoobjects.bases.rectplane_py import RectPlane_py as RectPlane
from owcsimpy.geoobjects.bases.cube_py import Cube_py as Cube

from owcsimpy.geoutils.cutils import calcRodriguesMtx as getRodriguesMtx
from owcsimpy.misc import flatten

import itertools

if __name__ == "__main__":

    room_L,room_W,room_H = [100,50,10]

    x = np.arange(0,room_L,5)
    y = np.arange(0,room_W,5)
    # xv,yv = np.meshgrid(np.arange(0,room_L,20),np.arange(0,room_W,10))

    x = (x[1]-x[0])/2+x;
    y = (y[1]-y[0])/2+y;

    locs = np.array(list(itertools.product(x,y,room_H*np.ones(np.size(x)))))

    vectors=[]

    for loc in locs:
        vectors.append(Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array(loc)))


    cube0 = Cube(
        Vector(np.array([1,np.deg2rad(0),np.deg2rad(0)])),
        ctrPoint = np.array([10,5,0.5]),
        dimensions = [14,6,1]
        )

    cube1 = Cube(
        Vector(np.array([1,np.deg2rad(0),np.deg2rad(0)])),
        ctrPoint = np.array([10,45,0.5]),
        dimensions = [14,6,1]
        )

    cube2 = Cube(
        Vector(np.array([1,np.deg2rad(0),np.deg2rad(0)])),
        ctrPoint = np.array([90,45,0.5]),
        dimensions = [14,6,1]
        )

    cube3 = Cube(
        Vector(np.array([1,np.deg2rad(0),np.deg2rad(0)])),
        ctrPoint = np.array([90,5,0.5]),
        dimensions = [14,6,1]
        )

    agv = Cube(
        Vector(np.array([1,np.deg2rad(0),np.deg2rad(0)])),
        ctrPoint = np.array([60,10,0.5]),
        dimensions = [1,0.5,1.5]
        )

    fig,ax = draw(xlim=[0,room_L],ylim=[0,room_W],zlim=[0,room_H])

    fig, ax = draw(figure=fig,axes=ax,
        vectors=vectors,
        colors='blue')

    fig, ax = draw(figure=fig,axes=ax,
        cubes=agv,
        enablevect=False,
        alphas=1)

    planesAGV = agv.getPartition(Ps=1)
    fig, ax = draw(figure=fig,axes=ax,
        colors='black',
        lengths=5,
        vectors=Vector(coord=planesAGV[3].normalVect,refPoint=planesAGV[3].ctrPoint,which='cartesian'))

    fig, ax = draw(figure=fig,axes=ax,
        planes=cube0.getPartition(Ps=1),
        alphas=1,
        enablevect=False,
        facecolors='gray')

    fig, ax = draw(figure=fig,axes=ax,
        planes=cube1.getPartition(Ps=1),
        alphas=1,
        enablevect=False,
        facecolors='gray')

    fig, ax = draw(figure=fig,axes=ax,
        planes=cube2.getPartition(Ps=1),
        alphas=1,
        enablevect=False,
        facecolors='gray')

    fig, ax = draw(figure=fig,axes=ax,
        planes=cube3.getPartition(Ps=1),
        alphas=1,
        enablevect=False,
        facecolors='gray')

    ax.plot([10,90,90,10,10],[10,10,40,40,10])
    ax.set_aspect(0.75)
    fig.set_size_inches(8, 6)
    ax.set_zticks([0,10])
    plt.show() 