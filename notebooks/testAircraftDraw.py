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

        self.backrest = Cube(
            Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
            ctrPoint = np.array([-0.1,0,0.65]),
            dimensions = [1.3,0.4,0.05]
            )
        # self.backrest = self.backrest.rotate(direction,np.array([0,0,1]));
        # self.backrest = self.backrest.translate(np.array([*loc,0]));
        self.backrest = self.backrest.rotate(direction,np.array([0,0,1]));
        self.backrest = self.backrest.translate(np.array([*loc,0])+self.backrest.ctrPoint);

        # Rd = getRodriguesMtx(direction,np.array([0,0,1]))
        # shiftedVect = (Rd@np.array([0.2-0.075,0,0]).reshape(3,1)).reshape(-1)
        # self.backrest = self.backrest.translate(shiftedVect+self.backrest.ctrPoint)  

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
            else:
                raise ValueError("'activity' is not recognized! \
                Should be either 'calling' or 'reading'")

            self.pd = BareDetector(polar,azimuth,loc,area=1e-4,FoV=np.deg2rad(90))


    def getPartition(self,Ps=1):

        planes=[]
        planes.append(self.backrest.getPartition(Ps=Ps))
        if self.head:
            planes.append(self.head.getPartition(Ps=Ps))
        if self.body:
            planes.append(self.body.getPartition(Ps=Ps))
        planes.append(self.chair.getPartition(Ps=Ps))
        
        planes = list(flatten(planes))
        return planes


if __name__ == "__main__":

    np.random.seed(1412);

    room_L,room_W,room_H = [5,5,3]

    human1 = HumanSitting(loc=np.array([2.5,2.5]),direction=np.deg2rad(0),activity='reading')
    human2 = HumanSitting(loc=np.array([2.5,3]),direction=np.deg2rad(0),isHuman=False)
    human2a = HumanSitting(loc=np.array([2.5,1]),direction=np.deg2rad(0),isHuman=False)
    human3 = HumanSitting(loc=np.array([1.9,2.5]),direction=np.deg2rad(0))
    human4 = HumanSitting(loc=np.array([1.9,3]),direction=np.deg2rad(0),isHuman=False)
    human4a = HumanSitting(loc=np.array([1.9,1]),direction=np.deg2rad(0),isHuman=False)


    plane1 = RectPlane(
        normalVect=Vector(np.array([1,np.deg2rad(0),np.deg2rad(0)])),
        ctrPoint=np.array([2.5,2.5,2]),
        RodriguesAngle=np.deg2rad(0),dimensions=[2,1])

    vector1 = Vector(np.array([0.4,np.deg2rad(180),np.deg2rad(0)]),
        refPoint=np.array([2.8,2.5,2]));

    vector2 = Vector(np.array([0.4,np.deg2rad(180),np.deg2rad(0)]),
        refPoint=np.array([2.0,2.5,2]));

    plane2 = RectPlane(
        normalVect=Vector(np.array([1,np.deg2rad(45),np.deg2rad(90)])),
        ctrPoint=np.array([2.5,1.65,2.35]),
        RodriguesAngle=np.deg2rad(0),dimensions=[1,2])

    plane3 = RectPlane(
        normalVect=Vector(np.array([1,np.deg2rad(45),np.deg2rad(90)])),
        ctrPoint=np.array([2.5,2.65,2.35]),
        RodriguesAngle=np.deg2rad(0),dimensions=[1,2])

    plane4 = RectPlane(
        normalVect=Vector(np.array([1,np.deg2rad(60),np.deg2rad(90)])),
        ctrPoint=np.array([2.5,3.25,1.55]),
        RodriguesAngle=np.deg2rad(0),dimensions=[1,2])

    plane5 = RectPlane(
        normalVect=Vector(np.array([1,np.deg2rad(90),np.deg2rad(90)])),
        ctrPoint=np.array([2.5,3.5,0.55]),
        RodriguesAngle=np.deg2rad(0),dimensions=[1.15,2])

    fig,ax = draw(xlim=[0,room_L],ylim=[0,room_W],zlim=[0,room_H],azim=-30,elev=-1)

    fig, ax = draw(figure=fig,axes=ax,
        planes=list(flatten([
            plane1,
            plane2,
            plane3,
            plane4,
            plane5,
            ])),
        facecolors='grey',
        alphas=0.3,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        planes=list(flatten([
            human1.getPartition(),
            human2.getPartition(),
            human2a.getPartition(),
            human3.getPartition(),
            human4.getPartition(),
            human4a.getPartition(),
            ])),
        facecolors='grey',
        alphas=1,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        vectors=list(flatten([
            vector1,
            vector2,
            ])),
        colors='blue')

    fig, ax = draw(figure=fig,axes=ax,
        circles=list(flatten([
            human1.pd,
            human3.pd,
            ])),
        scales=5e2)


    fig.set_size_inches(10, 8)
    ax.set_aspect(0.6)
    ax.axis('off')
    # ax.set_aspect('auto')
    # ax.set_xticks(np.arange(0,25,5))
    # ax.set_yticks(np.arange(0,25,5))
    # ax.set_zticks([0,3])
    plt.show() 

