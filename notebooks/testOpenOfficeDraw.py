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

class HumanSitting(object):
    def __init__(self,loc,direction,isHuman=True):

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
        self.human0 = HumanSitting(location+[-1.45,0.6],np.deg2rad(0),isHuman=np.random.binomial(1,0.5));
        self.human1 = HumanSitting(location+[1.45,0.6],np.deg2rad(180),isHuman=np.random.binomial(1,0.5));
        self.human2 = HumanSitting(location+[1.45,-0.6],np.deg2rad(180),isHuman=np.random.binomial(1,0.5));
        self.human3 = HumanSitting(location+[-1.45,-0.6],np.deg2rad(0),isHuman=np.random.binomial(1,0.5));

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
    # room = Room(dimensions=[30,10,3],identity=1)

    # room_L,room_W,room_H = [30,10,3] 
    room_L,room_W,room_H = [25,10,3] 

    # plane = RectPlane(
    #     normalVect=Vector(np.array([1,np.deg2rad(0),np.deg2rad(0)])),
    #     ctrPoint=np.array([0,0,0.9]),
    #     dimensions=[0.6,1.2])

    # human = HumanSitting(np.array([2,2]),np.deg2rad(0),isHuman=False);
    # human = HumanSitting(np.array([2,2]),np.deg2rad(30));

    cluster0 = Cluster(np.array([3,2]));
    cluster1 = Cluster(np.array([3,8]));
    cluster2 = Cluster(np.array([8,2]));
    cluster3 = Cluster(np.array([8,8]));

    cluster4 = Cluster(np.array([18,2]));
    cluster5 = Cluster(np.array([18,8]));
    cluster6 = Cluster(np.array([23,2]));
    cluster7 = Cluster(np.array([23,8]));

    planes = [cluster0.getPlanes(),
              cluster1.getPlanes(),
              cluster2.getPlanes(),
              cluster3.getPlanes(),
              cluster4.getPlanes(),
              cluster5.getPlanes(),
              cluster6.getPlanes(),
              cluster7.getPlanes(),
        ]
    planes = list(flatten(planes))

    # vectors = [
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([3,2,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([3,5,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([3,8,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([8,2,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([8,5,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([8,8,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([13,2,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([13,5,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([13,8,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([18,2,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([18,5,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([18,8,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([23,2,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([23,5,3])),
    #     Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([23,8,3])),
    # ]

    vectors = [
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([3,3,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([3,7,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([8,3,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([8,7,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([13,3,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([13,7,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([18,3,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([18,7,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([23,3,3])),
        Vector(np.array([1,np.deg2rad(180),0]),refPoint=np.array([23,7,3])),
    ]

    human = Human(loc=np.array([12,8]),direction=np.deg2rad(-30),
                          mode='rx',activity='calling',position='standing')


    humanBlocking = Human(loc=np.array([11,4]),direction=np.deg2rad(20),
                          mode='rx',activity='calling',position='standing')

    # planes = [planes,human.listPlanes]
    # planes = list(flatten(planes))

    # print(human.body.ctrPoint)
    # print(human.table.ctrPoint)

    fig,ax = draw(xlim=[0,room_L],ylim=[0,room_W],zlim=[0,room_H])

    # fig, ax = draw(figure=fig,axes=ax,
    #     planes=plane,
    #     enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        planes=planes,
        facecolors='gray',
        alphas=1,
        enablevect=False)

        # vectors=Vector(coord=human.normalVect,refPoint=human.ctrPoint,which='cartesian'),
        # planes=human.getPartition(),

    fig, ax = draw(figure=fig,axes=ax,
        planes=human.listPlanes,
        alphas=1,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        colors='black',
        vectors=Vector(coord=human.normalVect,refPoint=human.ctrPoint,which='cartesian'))

    fig, ax = draw(figure=fig,axes=ax,
        vectors=vectors,
        colors='blue')

    fig, ax = draw(figure=fig,axes=ax,
             circles=[human.pd],scales=5e3);

    fig, ax = draw(figure=fig,axes=ax,
        planes=humanBlocking.listPlanes,
        facecolors='gray',
        alphas=1,
        enablevect=False)

    fig, ax = draw(figure=fig,axes=ax,
        colors='black',
        vectors=Vector(coord=humanBlocking.normalVect,refPoint=humanBlocking.ctrPoint,which='cartesian'))

    fig.set_size_inches(10, 8)
    ax.set_aspect(0.75)
    # ax.set_aspect('auto')
    ax.set_zticks([0,3])
    plt.show()
