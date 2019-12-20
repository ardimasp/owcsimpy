import numpy as np
import matplotlib.pyplot as plt

from owcsimpy.geoutils.draw import draw
from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.cube_py import Cube_py as Cube
from owcsimpy.misc import flatten

class HumanCubes_py(object):
    """ This is a quick implementation of a more realistic human body model.

    Parameters
    ----------
    loc: ndarray(2,)
        Location of the human in xy-plane.
    direction: double
        Azimuth angle.
    reflectivities: dict
        The default value is {'hair':1,'face':1,'shirt':1}

    Notes
    -----
    The dimensions are fixed, which is 1.8 m height and 0.4 m width.

    Examples
    --------
    .. plot:: 
            :format: doctest
            :include-source: True

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> 
            >>> from owcsimpy.geoobjects.models.humancubes_py import HumanCubes_py as Human
            >>> from owcsimpy.geoutils.draw import draw
            >>> person = Human(np.array([1,1]),np.deg2rad(30),reflectivities={'hair':0.7,'face':0.5,'shirt':0.3})
            >>> 
            >>> draw(models3d=person,figsize=(6,6),azim=-100,elev=25,xlim=[0,3],ylim=[0,3],zlim=[0,3]);
            >>> 
            >>> plt.show()
    
    """
    
    def __init__(self,loc,direction,
                reflectivities={'hair':1,'face':1,'shirt':1}):
        
        assert len(loc)==2
        
        self.head = Cube(
            Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
            ctrPoint = np.array([0,0,1.7]),
            dimensions = [0.2,0.15,0.15],
            reflectivities={'p0':reflectivities['hair'],
                            'p1':reflectivities['face'],
                            'p2':reflectivities['hair'],
                            'p3':reflectivities['face'],
                            'p4':reflectivities['face'],
                            'p5':reflectivities['face']}
        )
        self.head = self.head.rotate(direction,np.array([0,0,1]));
        self.head = self.head.translate(np.array([*loc,self.head.ctrPoint[2]]));
        

        self.body = Cube(
            Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
            ctrPoint = np.array([0,0,1.2]),
            dimensions = [0.8,0.4,0.15],
            reflectivities={'p0':reflectivities['shirt'],
                            'p1':reflectivities['shirt'],
                            'p2':reflectivities['shirt'],
                            'p3':reflectivities['shirt'],
                            'p4':reflectivities['shirt'],
                            'p5':reflectivities['shirt']}
        )
        self.body = self.body.rotate(direction,np.array([0,0,1]));
        self.body = self.body.translate(np.array([*loc,self.body.ctrPoint[2]]));
        

        self.leg = Cube(
            Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
            ctrPoint = np.array([0,0,0.4]),
            dimensions = [0.8,0.3,0.15],
            reflectivities={'p0':reflectivities['shirt'],
                            'p1':reflectivities['shirt'],
                            'p2':reflectivities['shirt'],
                            'p3':reflectivities['shirt'],
                            'p4':reflectivities['shirt'],
                            'p5':reflectivities['shirt']}
        )
        self.leg = self.leg.rotate(direction,np.array([0,0,1]));
        self.leg = self.leg.translate(np.array([*loc,self.leg.ctrPoint[2]]));
        
        self.ctrPoint = np.array([0,0,1.2])+np.array([*loc,0])
        self.normalVect = self.head.normalVect;
        
        # This is for draw method
        self.listPlanes = self.getPartition(Ps=1)
        
    def getPartition(self,Ps=1,delta=None):
        
        planes=[]
        if delta == None:
            planes.append(self.head.getPartition(Ps=Ps))
            planes.append(self.body.getPartition(Ps=Ps))
            planes.append(self.leg.getPartition(Ps=Ps))
        else:
            if delta > 0.15:
                delta = 0.15
            planes.append(self.head.getPartition(delta=delta))
            planes.append(self.body.getPartition(delta=delta))
            planes.append(self.leg.getPartition(delta=delta))
        
        planes = list(flatten(planes))
        return planes
        
    def draw(self):
        fig,ax = draw(cubes=[self.head,self.body,self.leg],
             xlim=[0,3],ylim=[0,3],zlim=[0,3],
             lengths=0,enablevect='False'
                     );
        fig,ax = draw(figure=fig,axes=ax,
                      vectors=Vector(
                          coord=self.normalVect,refPoint=self.ctrPoint,which='cartesian'))
        
        