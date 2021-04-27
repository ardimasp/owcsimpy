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

from owcsimpy.geoobjects.bases.cube_py import Cube_py as Cube
from owcsimpy.geoobjects.bases.rectplane_py import RectPlane_py as RectPlane
from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoutils.draw import draw
from owcsimpy.geoobjects.models.humanwithactivity_py import HumanWithActivity_py as HumanWithActivity

from owcsimpy.geoobjects.models.baredetector_py import BareDetector_py as BareDetector

from owcsimpy.geoutils.cutils import calcRodriguesMtx as getRodriguesMtx
from owcsimpy.misc import flatten

from owcsimpy.geoobjects.models.roomcube_py import RoomCube_py as Room

class HumanSitting(object):
    def __init__(self,loc,direction,activity='calling',isHuman=True):

        offsetHeight = 0;

        if isHuman:
            self.head = Cube(
                Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
                ctrPoint = np.array([0,0,offsetHeight+1.3]),
                dimensions = [0.2,0.15,0.15]
            )
            self.head = self.head.rotate(direction,np.array([0,0,1]));
            self.head = self.head.translate(np.array([*loc,self.head.ctrPoint[2]]));

            self.body = Cube(
                Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
                ctrPoint = np.array([0,0,offsetHeight+0.8]),
                dimensions = [0.8,0.4,0.15],
            )
            self.body = self.body.rotate(direction,np.array([0,0,1]));
            self.body = self.body.translate(np.array([*loc,self.body.ctrPoint[2]]));

        else:
            self.head = None
            self.body = None

        self.chair = Cube(
            Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
            ctrPoint = np.array([0,0,0.2]),
            dimensions = [0.4,0.4,0.4]
        )
        self.chair = self.chair.rotate(direction,np.array([0,0,1]));
        self.chair = self.chair.translate(np.array([*loc,0])+self.chair.ctrPoint);

        Rd = getRodriguesMtx(direction,np.array([0,0,1]))
        shiftedVect = (Rd@np.array([0.2-0.075,0,0]).reshape(3,1)).reshape(-1)
        self.chair = self.chair.translate(shiftedVect+self.chair.ctrPoint)        

        self.table = Cube(
            Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
            ctrPoint = np.array([0,0,0.3]),
            dimensions = [0.6,1.2,0.9],
        )
        self.table = self.table.rotate(direction,np.array([0,0,1]));
        self.table = self.table.translate(np.array([*loc,0])+self.table.ctrPoint);

        Rd = getRodriguesMtx(direction,np.array([0,0,1]))
        shiftedVect = (Rd@np.array([1,0,0]).reshape(3,1)).reshape(-1)
        self.table = self.table.translate(shiftedVect+self.table.ctrPoint)
        
        # The table is modeled a simple plane
        self.table = self.table.listPlanes[2]

        self.ctrPoint = np.array([0,0,1.2])+np.array([*loc,0])
        self.listPlanes = self.getPartition(Ps=1)
        self.activity = activity

        self.mode = None
        self.pd = None
        self.led = None

        self.isHuman = isHuman;
        
        if isHuman:
            if activity.lower() == 'calling':
                polar = np.deg2rad(45)
                azimuth = self.listPlanes[1]._RectPlane_py__normalVect.spherical[2]+np.pi
                loc = self.head.ctrPoint+0.15*self.head.listPlanes[4].normalVect
            elif activity.lower() == 'reading':
                polar = np.deg2rad(50)
                azimuth = self.listPlanes[1]._RectPlane_py__normalVect.spherical[2]+np.pi
                loc = (self.ctrPoint
                       +0.5*self.listPlanes[1].normalVect
                       +np.array([0,0,-0.2+offsetHeight]))
            elif activity.lower() == 'usbdongle':
                polar = np.deg2rad(0)
                azimuth = self.listPlanes[1]._RectPlane_py__normalVect.spherical[2]+np.pi
                loc = (self.ctrPoint
                       +1*self.listPlanes[1].normalVect
                       +np.array([0,0,-0.6]))
            else:
                raise ValueError("'activity' is not recognized! \
                Should be either 'calling' or 'reading'")

            self.pd = BareDetector(polar,azimuth,loc,area=1e-4,FoV=np.deg2rad(90))


    def getPartition(self,Ps=1):

        planes=[]
        if self.head:
            planes.append(self.head.getPartition(Ps=Ps))
        if self.body:
            planes.append(self.body.getPartition(Ps=Ps))
        planes.append(self.chair.getPartition(Ps=Ps))
        planes.append(self.table.getPartition(Ps=Ps))
        
        planes = list(flatten(planes))
        return planes

class Cluster(object):
    def __init__(self,location):

        activities = ['calling','reading','usbdongle'];

        self.human0 = HumanSitting(location+[-1.45,0.6],np.deg2rad(0),activity=activities[np.random.randint(3,size=1)[0]],isHuman=np.random.binomial(1,0.5));
        self.human1 = HumanSitting(location+[1.45,0.6],np.deg2rad(180),activity=activities[np.random.randint(3,size=1)[0]],isHuman=np.random.binomial(1,0.5));
        self.human2 = HumanSitting(location+[1.45,-0.6],np.deg2rad(180),activity=activities[np.random.randint(3,size=1)[0]],isHuman=np.random.binomial(1,0.5));
        self.human3 = HumanSitting(location+[-1.45,-0.6],np.deg2rad(0),activity=activities[np.random.randint(3,size=1)[0]],isHuman=np.random.binomial(1,0.5));

    def getPlanes(self):
        planes = []
        planes.append(self.human0.getPartition());
        planes.append(self.human1.getPartition());
        planes.append(self.human2.getPartition());
        planes.append(self.human3.getPartition());

        planes = list(flatten(planes))
        return planes

if __name__ == "__main__":

    np.random.seed(1412);

    room_L,room_W,room_H = [20,20,3]

    rho_keys = ['b','t','s','n','e','w']
    rho_vals = [0.,0.,0.,0.,0.,0.]
    reflectivities = {rho_keys[i]:rho_vals[i] for i in range(6)}

    room = Room(dimensions=[20,20,3],identity=1,reflectivities=reflectivities)

    clusters = [
        Cluster(np.array([ 2.5,2.5])),
        Cluster(np.array([ 7.5,2.5])),
        Cluster(np.array([12.5,2.5])),
        Cluster(np.array([17.5,2.5])),
        Cluster(np.array([ 2.5,7.5])),
        Cluster(np.array([ 7.5,7.5])),
        Cluster(np.array([12.5,7.5])),
        Cluster(np.array([17.5,7.5])),
        Cluster(np.array([ 2.5,12.5])),
        Cluster(np.array([ 7.5,12.5])),
        Cluster(np.array([12.5,12.5])),
        Cluster(np.array([17.5,12.5])),
        Cluster(np.array([ 2.5,17.5])),
        Cluster(np.array([ 7.5,17.5])),
        Cluster(np.array([12.5,17.5])),
        Cluster(np.array([17.5,17.5]))
    ]

    planes = []
    pds = []
    for cluster in clusters:
        planes.append(cluster.getPlanes())

        if cluster.human0.isHuman:
            pds.append(cluster.human0.pd)
        if cluster.human1.isHuman:
            pds.append(cluster.human1.pd)
        if cluster.human2.isHuman:
            pds.append(cluster.human2.pd)
        if cluster.human3.isHuman:
            pds.append(cluster.human3.pd)

    planes = list(flatten(planes))
    pds = list(flatten(pds))

    vectors = [
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 2.5,2.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([   5,2.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 7.5,2.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  10,2.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([12.5,2.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  15,2.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([17.5,2.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 2.5,5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([   5,5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 7.5,5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  10,5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([12.5,5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  15,5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([17.5,5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 2.5,7.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([   5,7.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 7.5,7.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  10,7.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([12.5,7.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  15,7.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([17.5,7.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 2.5,10,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([   5,10,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 7.5,10,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  10,10,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([12.5,10,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  15,10,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([17.5,10,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 2.5,12.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([   5,12.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 7.5,12.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  10,12.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([12.5,12.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  15,12.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([17.5,12.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 2.5,15,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([   5,15,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 7.5,15,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  10,15,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([12.5,15,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  15,15,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([17.5,15,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 2.5,17.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([   5,17.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([ 7.5,17.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  10,17.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([12.5,17.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([  15,17.5,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([17.5,17.5,3]))
    ]


    fig,ax = draw(xlim=[0,room_L],ylim=[0,room_W],zlim=[0,room_H])

    fig, ax = draw(figure=fig,axes=ax,
        planes=planes,
        facecolors='grey',
        alphas=1,
        enablevect=False)

    # fig, ax = draw(figure=fig,axes=ax,
    #     circles=pds,scales=5e3)

    # fig, ax = draw(figure=fig,axes=ax,
    #     vectors=vectors,
    #     colors='blue')

    human = HumanWithActivity(loc=np.array([15,3]),direction=np.deg2rad(-30),
                mode='rx',activity='calling',position='standing')

    fig, ax = draw(figure=fig,axes=ax,
        planes=human.listPlanes,
        facecolors='grey',
        alphas=1,
        enablevect=False)

    # fig, ax = draw(figure=fig,axes=ax,
    #     colors='black',
    #     vectors=Vector(coord=human.normalVect,refPoint=human.ctrPoint,which='cartesian'))

    # fig, ax = draw(figure=fig,axes=ax,
    #     circles=human.pd,scales=5e3)

    # fig, ax = draw(figure=fig,axes=ax,
    #     models3d=room,
    #     enablevect=False)

    # r = 1;
    # u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    # x = r*np.cos(u)*np.sin(v) + 5
    # y = r*np.sin(u)*np.sin(v) + 5
    # z = r*np.cos(v) + 3
    # ax.plot_wireframe(x, y, z, color="r", alpha=0.5)

    # r = 1;
    # u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    # x = r*np.cos(u)*np.sin(v) + 15
    # y = r*np.sin(u)*np.sin(v) + 5
    # z = r*np.cos(v) + 3
    # ax.plot_wireframe(x, y, z, color="r", alpha=0.5)

    # r = 1;
    # u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    # x = r*np.cos(u)*np.sin(v) + 5
    # y = r*np.sin(u)*np.sin(v) + 15
    # z = r*np.cos(v) + 3
    # ax.plot_wireframe(x, y, z, color="r", alpha=0.5)

    # r = 1;
    # u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    # x = r*np.cos(u)*np.sin(v) + 15
    # y = r*np.sin(u)*np.sin(v) + 15
    # z = r*np.cos(v) + 3
    # ax.plot_wireframe(x, y, z, color="r", alpha=0.5)



    fig.set_size_inches(20, 16)
    ax.set_aspect(0.16)
    # ax.set_aspect('auto')
    # ax.set_xticks(np.arange(0,25,5))
    # ax.set_yticks(np.arange(0,25,5))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    # plt.xticks([],[])
    # ax.axis('off')

    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    # ax.get_zaxis().set_visible(False)

    # frame1 = plt.gca()
    # frame1.axes.get_xaxis().set_visible(False)
    # frame1.axes.get_yaxis().set_visible(False)

    # ax.set_position([20,20,3])

    plt.show() 