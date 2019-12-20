# SimpleOfficeEnv_py
import numpy as np
import matplotlib.pyplot as plt

from owcsimpy.geoobjects.models.humanwithactivity_py import HumanWithActivity_py as Human
from owcsimpy.geoutils.draw import draw

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.cube_py import Cube_py as Cube
from owcsimpy.geoobjects.models.pointsource_py import PointSource_py as PointSource
from owcsimpy.geoobjects.models.baredetector_py import BareDetector_py as BareDetector
from owcsimpy.geoobjects.models.roomcube_py import RoomCube_py as Room
from owcsimpy.geoutils.cutils import calcRodriguesMtx as getRodriguesMtx
from owcsimpy.misc import flatten
from scipy.spatial import Delaunay

class SimpleOfficeEnv_py(object):
    """ A model of simple office environments.
    
    A simple office environment comprises only 
        - a LED or a PD 
        - a human model with different activities and positions
        - a furniture (a table and a chair)
        
    Parameters
    ----------
    roomDim: list
        The dimensions of the room in [Length,Width,Height].
    humanLoc
        The human's location.
    humanDirection
        The humand's direction.
    chairLoc
        The chair's location. Note that the location of the table is
        defined relative to the chair's location.
    chairDirection
        The chair's direction.
    mode: {'tx','rx'}
        Describe whether transmitting (tx) or receiving (rx)
    activity: {'calling','reading','usbdongle'}
        'calling' denotes the LED or the PD is near the right ear.
        Meanwhile, 'reading' denotes the device is in front of the 
        human's chest. 'usbdongle' denotes the use of LiFi usb dongles
        attached to a laptop on a desk.
    position: {'standing','sitting'}
        The human's positions. 
    
    Attributes
    ----------
    human: HumanWithActivity_py
        A human model with different positions and activities
    room: RoomCube_py
        A simple room model with a cube.
    chair: Cube_py
        A chair model with a cube.
    table: RectPlane_py
        A table model with a rectangular plane.
    tableAsCube: Cube_py
        A dummy variable to model the table as a cube. This is mostly used
        for the validation.
    led: PointSource_py
        A LED model.
    pd: BareDetector_py
        A PD model.
    listPlanes: list
        List of RectPlane_py
    isValid: bool
        To show that all elements are properly placed inside the room.
        It means that there is no element that is out of bounds and 
        no element that intersects to one another.
        
    Examples
    --------
    .. plot::
            :format: doctest
            :include-source: True
            
            >>> import numpy as np
            >>> from owcsimpy.geoobjects.models.simpleofficeenv_py import SimpleOfficeEnv_py as OfficeEnv
            >>> 
            >>> office = OfficeEnv(roomDim=[4,3,3],
            >>>                    humanLoc=[1,1.5],humanDirection=np.deg2rad(90),
            >>>                    activity='usbdongle',position='sitting')
            >>> 
            >>> 
            >>> office.draw(displayRoom=False)
            >>> print('is valid:', office.isValid)
    
    """
    
    def __init__(self,
                 roomDim,
                 humanLoc,humanDirection,
                 chairLoc=None,chairDirection=0,
                 mode='rx',
                 activity='calling',
                 position='standing'):
        
        # Assertion checks
        assert len(roomDim)==3 and np.alltrue(np.array(roomDim)>0) 
            
        assert mode.lower() == 'tx' or mode.lower() == 'rx',\
            ("'mode' is not recognized! Should be either 'tx' or 'rx'")
        
        assert activity.lower() == 'calling' or activity.lower() == 'reading' or activity.lower() == 'usbdongle',\
            ("'activity' is not recognized! Should be either 'calling', 'reading' or 'usbdongle'")
        
        assert position.lower() == 'standing' or position.lower() == 'sitting',\
            ("'mode' is not recognized! Should be either 'standing' or 'sitting'") 
        
        assert not(activity.lower() == 'usbdongle' and position.lower()=='standing'),\
            ("cannot choose the `usbdongle` aciticity and the `standing` position.")
        
        
        # Human and other elements
        self.human = Human(loc=humanLoc,direction=humanDirection,
                          mode=mode,activity=activity,position=position)
        
        
        # Room
        # Reflectivities
        # Except the floor, which is modeled as a pinewood, other 
        # parts are modeled with plasters
        rho_keys = ['b','t','s','n','e','w']
        
        if mode.lower() == 'tx':
            # IR-LED (uplink)
            rho_vals = [0.92,0.83,0.83,0.83,0.83,0.83]
        elif mode.lower() == 'rx':
            # VL-LED (downlink)
            rho_vals = [0.54,0.76,0.76,0.76,0.76,0.76]
            
        reflectivities = {rho_keys[i]:rho_vals[i] for i in range(len(rho_keys))}

        self.room = Room(dimensions=roomDim,identity=1,reflectivities=reflectivities)
        
        # A chair and a table
        self.chair, self.table, self.tableAsCube = None,None,None
        if chairLoc:
            self.chair = Cube(
                Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
                ctrPoint = np.array([0,0,0.2]),
                dimensions = [0.4,0.4,0.4],
                reflectivities={'p0':0.04,
                                'p1':0.04,
                                'p2':0.04,
                                'p3':0.04,
                                'p4':0.04,
                                'p5':0.04}
            )
            self.chair = self.chair.rotate(chairDirection,np.array([0,0,1]));
            self.chair = self.chair.translate(np.array([*chairLoc,0])+self.chair.ctrPoint);

            Rd = getRodriguesMtx(chairDirection,np.array([0,0,1]))
            shiftedVect = (Rd@np.array([0.2,0,0]).reshape(3,1)).reshape(-1)
            self.chair = self.chair.translate(shiftedVect+self.chair.ctrPoint)

            # The table's location is defined based on 
            self.table = Cube(
                Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
                ctrPoint = np.array([0,0,0.3]),
                dimensions = [0.6,1.2,0.9],
                reflectivities={'p0':0.04,
                                'p1':0.04,
                                'p2':0.04,
                                'p3':0.04,
                                'p4':0.04,
                                'p5':0.04}
            )
            self.table = self.table.rotate(chairDirection,np.array([0,0,1]));
            self.table = self.table.translate(np.array([*chairLoc,0])+self.table.ctrPoint);

            Rd = getRodriguesMtx(chairDirection,np.array([0,0,1]))
            shiftedVect = (Rd@np.array([1,0,0]).reshape(3,1)).reshape(-1)
            self.table = self.table.translate(shiftedVect+self.table.ctrPoint)

            # The table is modeled a simple plane
            self.tableAsCube = self.table
            self.table = self.table.listPlanes[2]
        
        
        # LED or PD
        self.led,self.pd = None,None
    
        if mode.lower() == 'rx':
            self.led = PointSource(np.pi,0,
                                   np.array([self.room.L/2,self.room.W/2,self.room.H]),
                                   m=-np.log(2)/np.log(np.cos(np.deg2rad(40.))))
        elif mode.lower() == 'tx':
            self.pd = BareDetector(np.pi,0,
                                   np.array([self.room.L/2,self.room.W/2,self.room.H]),
                                   area=1e-4,FoV=np.deg2rad(85))
        
        # Get all planes
        self.listPlanes = []
        self.listPlanes.append(self.human.listPlanes)
        self.listPlanes.append(self.room.listPlanes)
        if self.chair and self.table:
            self.listPlanes.append(self.chair.getPartition(Ps=1))
            self.listPlanes.append(self.table.getPartition(Ps=1))
        self.listPlanes = list(flatten(self.listPlanes))

        # Get blocking object
        self.blockingObj = []
        if self.chair:
            self.blockingObj.append(self.chair)
        if self.table:
            self.blockingObj.append(self.table)
        if self.human:
            self.blockingObj.append(self.human)
        
        # Validate whether a vertice is NOT out of bounds
        roomVerts = np.append(self.room.listPlanes[0].verts,self.room.listPlanes[1].verts).reshape([-1,3])
        
        delaunay = Delaunay(roomVerts)

        verts = np.array(list(flatten([plane.verts.tolist() for plane in self.listPlanes]))).reshape([-1,3])
        
        if self.human.led:
            verts = np.append(verts,self.human.led.loc).reshape([-1,3])
        elif self.human.pd:
            verts = np.append(verts,self.human.pd.loc).reshape([-1,3])
        
#         # Check feasibility of human
#         verts = np.array(list(flatten([plane.verts.tolist() for plane in self.human.listPlanes]))).reshape([-1,3])
        
#         # Check feasibility of chair
#         verts = np.array(list(flatten([plane.verts.tolist() for plane in self.chair.getPartition(Ps=1)]))).reshape([-1,3])
        
#         # Check feasibility of table
#         verts = np.array(list(flatten([plane.verts.tolist() for plane in self.table.getPartition(Ps=1)]))).reshape([-1,3])
        
#         # Check feasibility of LED
#         verts = self.human.led.loc
        
#         # Check feasibility of PD
#         verts = self.human.pd.loc
    
        self.isValid = np.alltrue(delaunay.find_simplex(verts)>=0)
        
        
        
        # Check whether the human intersects with the furniture
        if self.chair and self.table:    
            chairVerts = np.append(self.chair.listPlanes[0].verts,self.chair.listPlanes[1].verts).reshape([-1,3])
            tableVerts = np.append(self.tableAsCube.listPlanes[0].verts,self.tableAsCube.listPlanes[1].verts).reshape([-1,3])
            
            verts = np.array(list(flatten([plane.verts.tolist() for plane in self.human.listPlanes]))).reshape([-1,3])
            
            if self.human.led:
                verts = np.append(verts,self.human.led.loc).reshape([-1,3])
            elif self.human.pd:
                verts = np.append(verts,self.human.pd.loc).reshape([-1,3])
            
            self.isValid = np.alltrue(Delaunay(chairVerts).find_simplex(verts)<0) and np.alltrue(Delaunay(tableVerts).find_simplex(verts)<0)
            
    def draw(self,
             displayRoom=False,
             displayReflectivity=False):
        
        """ Draw.
        
        Parameters
        ----------
        displayRoom: bool
        displayReflectivity: bool
        
        
        """
        
        if displayRoom:
            if displayReflectivity:
                alphas=[0.5*plane.reflectivity for plane in self.listPlanes]
            else:
                alphas=[1 for plane in self.listPlanes]
                
            fig,ax = draw(xlim=[0,self.room.L],ylim=[0,self.room.W],zlim=[0,self.room.H],
                         planes=self.listPlanes,
                         vectors=Vector(coord=self.human.normalVect,refPoint=self.human.ctrPoint,which='cartesian'),
                         enablevect=False,
                         alphas=alphas);
        else:
            listPlanes = []
            listPlanes.append(self.human.listPlanes)
            if self.chair and self.table:
                listPlanes.append(self.chair.getPartition(Ps=1))
                listPlanes.append(self.table.getPartition(Ps=1))
            listPlanes = list(flatten(listPlanes))
            
            if displayReflectivity:
                alphas=[0.5*plane.reflectivity for plane in listPlanes]
            else:
                alphas=[1 for plane in listPlanes]
            
            fig,ax = draw(xlim=[0,self.room.L],ylim=[0,self.room.W],zlim=[0,self.room.H],
                         planes=listPlanes,
                         vectors=Vector(coord=self.human.normalVect,refPoint=self.human.ctrPoint,which='cartesian'),
                         enablevect=False,
                         alphas=alphas);
        
        fig, ax = draw(figure=fig,axes=ax,
             circles=[self.human.pd,self.pd],scales=5e3,
             vectors=[self.human.led,self.led]);

        return fig,ax

        