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

class HumanSitting(object):
    def __init__(self,loc,direction,
        activity='calling',
        isHuman=True):

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

        # self.table = Cube(
        #     Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
        #     ctrPoint = np.array([0,0,0.3]),
        #     dimensions = [0.6,1.2,0.9],
        # )
        # self.table = self.table.rotate(direction,np.array([0,0,1]));
        # self.table = self.table.translate(np.array([*loc,0])+self.table.ctrPoint);

        # Rd = getRodriguesMtx(direction,np.array([0,0,1]))
        # shiftedVect = (Rd@np.array([1,0,0]).reshape(3,1)).reshape(-1)
        # self.table = self.table.translate(shiftedVect+self.table.ctrPoint)
        
        # # The table is modeled a simple plane
        # self.table = self.table.listPlanes[2]

        self.ctrPoint = np.array([0,0,1.2])+np.array([*loc,0])
        self.listPlanes = self.getPartition(Ps=1)
        self.activity = activity

        self.mode = None
        self.pd = None
        self.led = None
        
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
        # planes.append(self.table.getPartition(Ps=Ps))
        
        planes = list(flatten(planes))
        return planes

if __name__ == "__main__":

    np.random.seed(1412);


    room_L,room_W,room_H = [10,10,3] 

    v = Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)]))
    plane1 = RectPlane(normalVect=v,ctrPoint=np.array([5,5,1.5]),
        RodriguesAngle=np.deg2rad(0),dimensions=[3,10])

    v = Vector(np.array([1,np.deg2rad(90),np.deg2rad(90)]))
    plane2 = RectPlane(normalVect=v,ctrPoint=np.array([2.5,7,1.5]),
        RodriguesAngle=np.deg2rad(0),dimensions=[3,5])

    v = Vector(np.array([1,np.deg2rad(90),np.deg2rad(90)]))
    plane3 = RectPlane(normalVect=v,ctrPoint=np.array([2.5,4,1.5]),
        RodriguesAngle=np.deg2rad(0),dimensions=[3,5])

    v = Vector(np.array([1,np.deg2rad(0),np.deg2rad(0)]))
    plane_desk = RectPlane(normalVect=v,ctrPoint=np.array([7.5,5,0.5]),
        RodriguesAngle=np.deg2rad(0),dimensions=[1,2])

    cube1 = Cube(
        Vector(np.array([1,np.deg2rad(0),np.deg2rad(0)])),
        ctrPoint = np.array([8,6,0.25]),
        dimensions = [0.01,0.01,0.5],
        RodriguesAngle = np.deg2rad(0)
    )
    planesCube1 = cube1.getPartition(1)

    cube2 = Cube(
        Vector(np.array([1,np.deg2rad(0),np.deg2rad(0)])),
        ctrPoint = np.array([7,6,0.25]),
        dimensions = [0.01,0.01,0.5],
        RodriguesAngle = np.deg2rad(0)
    )
    planesCube2 = cube2.getPartition(1)

    cube3 = Cube(
        Vector(np.array([1,np.deg2rad(0),np.deg2rad(0)])),
        ctrPoint = np.array([7,4,0.25]),
        dimensions = [0.01,0.01,0.5],
        RodriguesAngle = np.deg2rad(0)
    )
    planesCube3 = cube3.getPartition(1)

    cube4 = Cube(
        Vector(np.array([1,np.deg2rad(0),np.deg2rad(0)])),
        ctrPoint = np.array([8,4,0.25]),
        dimensions = [0.01,0.01,0.5],
        RodriguesAngle = np.deg2rad(0)
    )
    planesCube4 = cube4.getPartition(1)

    human1 = HumanSitting(loc=np.array([9.8,9.8]),direction=np.deg2rad(-90))
    human2 = HumanSitting(loc=np.array([9.4,9.8]),direction=np.deg2rad(-90),isHuman=False)    
    human2a = HumanSitting(loc=np.array([9.0,9.8]),direction=np.deg2rad(-90),isHuman=False)    
    human3 = HumanSitting(loc=np.array([8.6,9.8]),direction=np.deg2rad(-90),activity='reading')    

    human_sitting = HumanWithActivity(loc=[1.2,0.6],direction=np.deg2rad(180),
             activity='usbdongle',position='sitting',mode='rx')

    human_standing = HumanWithActivity(loc=[7,2],direction=np.deg2rad(30),
             activity='reading',position='standing',mode='rx')

    fig,ax = draw(xlim=[0,room_L],ylim=[0,room_W],zlim=[0,room_H])

    fig, ax = draw(figure=fig,axes=ax,
        planes=[plane1, plane2, plane3],
        facecolors='gray',
        alphas=0.7,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        planes=human1.getPartition(),
        facecolors='gray',
        alphas=1.,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        planes=human2.getPartition(),
        facecolors='gray',
        alphas=1.,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        planes=human2a.getPartition(),
        facecolors='gray',
        alphas=1.,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        planes=human3.getPartition(),
        facecolors='gray',
        alphas=1.,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        planes=human_sitting.listPlanes,
        facecolors='gray',
        alphas=1.,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        planes=human_standing.listPlanes,
        facecolors='gray',
        alphas=1.,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        planes=plane_desk,
        facecolors='gray',
        alphas=1.,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        circles=human1.pd,scales=5e3)

    fig, ax = draw(figure=fig,axes=ax,
        circles=human3.pd,scales=5e3)

    fig, ax = draw(figure=fig,axes=ax,
        circles=human_sitting.pd,scales=5e3)

    fig, ax = draw(figure=fig,axes=ax,
        circles=human_standing.pd,scales=5e3)

    fig, ax = draw(figure=fig,axes=ax,
        planes=planesCube1,
        facecolors='black',
        alphas=1.,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        planes=planesCube2,
        facecolors='black',
        alphas=1.,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        planes=planesCube2,
        facecolors='black',
        alphas=1.,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        planes=planesCube3,
        facecolors='black',
        alphas=1.,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        planes=planesCube4,
        facecolors='black',
        alphas=1.,
        enablevect=False)

    ap1 = Vector(np.array([1,np.deg2rad(180),np.deg2rad(0)]), 
        refPoint=np.array([7.5,7.5,3]))

    ap2 = Vector(np.array([1,np.deg2rad(180),np.deg2rad(0)]), 
        refPoint=np.array([7.5,5,3]))

    ap3 = Vector(np.array([1,np.deg2rad(180),np.deg2rad(0)]), 
        refPoint=np.array([7.5,2.5,3]))

    ap4 = Vector(np.array([1,np.deg2rad(180),np.deg2rad(0)]), 
        refPoint=np.array([2.5,8,3]))

    ap5 = Vector(np.array([1,np.deg2rad(180),np.deg2rad(0)]), 
        refPoint=np.array([2.5,5.5,3]))

    ap6 = Vector(np.array([1,np.deg2rad(180),np.deg2rad(0)]), 
        refPoint=np.array([2.5,2,3]))

    fig, ax = draw(figure=fig,axes=ax,
        vectors=[ap1,ap2,ap3,ap4,ap5,ap6],
        colors='blue')

    r = 1;
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = r*np.cos(u)*np.sin(v) + 5
    y = r*np.sin(u)*np.sin(v) + 5
    z = r*np.cos(v) + 1.5
    ax.plot_wireframe(x, y, z, color="r", alpha=0.5)

    # N=50
    # stride=2
    # # ax = axes[0,0]
    # u = np.linspace(0, 2 * np.pi, N)
    # v = np.linspace(0, np.pi, N)
    # x = 1*np.outer(np.cos(u), np.sin(v)) + 10
    # y = 1*np.outer(np.sin(u), np.sin(v)) + 5
    # z = 1*np.outer(np.ones(np.size(u)), np.cos(v)) + 1.5
    # ax.plot_surface(x, y, z, linewidth=0.0, cstride=stride, rstride=stride, alpha=1)
    # ax.set_title('{0}x{0} data points, stride={1}'.format(N,stride))


    fig.set_size_inches(10, 8)
    ax.set_aspect(0.3)
    # ax.set_aspect('auto')
    ax.set_zticks([0,3])
    plt.show()