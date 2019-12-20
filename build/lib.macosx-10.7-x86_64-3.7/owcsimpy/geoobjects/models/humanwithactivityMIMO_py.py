# HumanWithActivityMIMO_py
import numpy as np
import matplotlib.pyplot as plt

from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
from owcsimpy.geoobjects.bases.cube_py import Cube_py as Cube
from owcsimpy.geoobjects.models.pointsource_py import PointSource_py as PointSource
from owcsimpy.geoobjects.models.baredetector_py import BareDetector_py as BareDetector
from owcsimpy.geoutils.cutils import calcRodriguesMtx as getRodriguesMtx
from owcsimpy.misc import flatten

class HumanWithActivityMIMO_py(object):
    """ This is a quick implementation of a more realistic human body model and 
    its activity and its position.

    Parameters
    ----------
    loc: ndarray(2,)
        Location of the human in xy-plane.
    direction: double
        Azimuth angle.
    reflectivities: dict
        The default value is 
            {'skin':0.66,'shirt':0.64,'furniture':0.04} for tx mode (uplink IR)
            {'skin':0.42,'shirt':0.61,'furniture':0.04} for rx mode (downlink VL)
    mode: {'tx','rx'}
        Describe whether transmitting (tx) or receiving (rx)
    activity: {'calling','reading','usbdongle'}
        'calling' denotes the LED or the PD is near the right ear.
        Meanwhile, 'reading' denotes the device is in front of the 
        human's chest. 'usbdongle' denotes the use of LiFi usb dongles
        attached to a laptop on a desk.
        
    Attributes
    ----------
    head: Cube_py
        Modeling the head
    body: Cube_py
        Modeling the body
    leg: Cube_py
        Modeling the leg
    chair: Cube_py
        Modeling the chair
    table: RectPlane_py
        Modeling the table
    ctrPoint: ndarray(3,)
        The center point of the human model
    normalVect: ndarray(3,)
        The human's direction
    listPlanes: list
        List of RectPlane_py
    mode: {'tx','rx'}
        It is either tx mode or rx mode
    pd: BareDetector_py
        Modeling the PD if the mode is 'rx'
    led: PointSource_py
        Modeling the LED if the mode is 'tx'

    Notes
    -----
    The dimensions are fixed.
        a. head (15cm x 15cm x 20cm)
        b. body (40cm x 15cm x 80cm)
        c. legs (30cm x 15cm x 80cm)
        d. chair (40cm x 40cm x 40cm)
        e. desk (220cm x 90cm) with the height of 60cm

    The geometries and locations of the LED and PD are also fixed.
        

    Examples
    --------
    .. plot:: 
            :format: doctest
            :include-source: True

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> 
            >>> from owcsimpy.geoobjects.models.humanwithactivity_py import HumanWithActivity_py as Human
            >>> from owcsimpy.geoutils.draw import draw
            >>> # Calling-sitting-rx
            >>> human_sitting = Human(loc=[2,2],direction=np.deg2rad(45),
            >>>              activity='calling',position='sitting',mode='rx')
            >>>  
            >>> # Calling-standing-tx
            >>> human_standing = Human(loc=[2,2],direction=np.deg2rad(45),
            >>>              activity='calling',position='standing',mode='tx')
            >>>  
            >>> # Plot
            >>> fig,axs = draw(subplots=True,figsize=(6,4),nrows=1,ncols=2,xlim=[0,3],ylim=[0,3],zlim=[0,3])
            >>>  
            >>> draw(figure=fig,axes=axs[0],planes=human_sitting.listPlanes,
            >>>      vectors=Vector(coord=human_sitting.normalVect,refPoint=human_sitting.ctrPoint,which='cartesian'),
            >>>      enablevect='False',
            >>>      alphas=[plane.reflectivity for plane in human_sitting.listPlanes]);
            >>>  
            >>> draw(figure=fig,axes=axs[0],
            >>>      circles=human_sitting.pd,scales=5e3,
            >>>      vectors=human_sitting.led);
            >>>  
            >>> draw(figure=fig,axes=axs[1],planes=human_standing.listPlanes,
            >>>      vectors=Vector(coord=human_standing.normalVect,refPoint=human_standing.ctrPoint,which='cartesian'),
            >>>      enablevect='False',
            >>>      alphas=[plane.reflectivity for plane in human_standing.listPlanes]);
            >>>  
            >>> draw(figure=fig,axes=axs[1],
            >>>      circles=human_standing.pd,scales=5e3,
            >>>      vectors=human_standing.led);
            >>> 
            >>> plt.show()
    
    """
    
    def __init__(self,loc,direction,
                reflectivities={'skin':1,'shirt':1,'furniture':1},
                mode='rx',
                activity='calling',
                position='standing',
                m=1,
                area=1e-4,FoV=np.deg2rad(90)):
        
        # Assertion checks
        assert len(loc)==2 and m > 0 and area > 0 and 0 <= FoV <= np.pi/2
        
#         assert sorted( 
#            [string.lower() for string in list(reflectivities.keys())]
#            ) == sorted(['chair','face','hair','shirt','table']), "keys name are wrong"
        
        assert mode.lower() == 'tx' or mode.lower() == 'rx',\
            ("'mode' is not recognized! Should be either 'tx' or 'rx'")
        
        assert activity.lower() == 'calling' or activity.lower() == 'reading' or activity.lower() == 'usbdongle',\
            ("'activity' is not recognized! Should be either 'calling', 'reading' or 'usbdongle'")
        
        assert position.lower() == 'standing' or position.lower() == 'sitting',\
            ("'mode' is not recognized! Should be either 'standing' or 'sitting'") 
        
        assert not(activity.lower() == 'usbdongle' and position.lower()=='standing'),\
            ("cannot choose the `usbdongle` aciticity and the `standing` position.")
        
        # 
        offsetHeight = 0.4 if position.lower() == 'standing' else 0
        
        if mode.lower() == 'tx':
            # IR-LED (uplink)
            reflectivities={'skin':0.66,'shirt':0.64,'furniture':0.04}
        elif mode.lower() == 'rx':
            # VL-LED (downling)
            reflectivities={'skin':0.42,'shirt':0.61,'furniture':0.04}
        
        # Head part
        self.head = Cube(
            Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
            ctrPoint = np.array([0,0,offsetHeight+1.3]),
            dimensions = [0.2,0.15,0.15],
            reflectivities={'p0':reflectivities['skin'],
                            'p1':reflectivities['skin'],
                            'p2':reflectivities['skin'],
                            'p3':reflectivities['skin'],
                            'p4':reflectivities['skin'],
                            'p5':reflectivities['skin']}
        )
        self.head = self.head.rotate(direction,np.array([0,0,1]));
        self.head = self.head.translate(np.array([*loc,self.head.ctrPoint[2]]));
        
        # Body part
        self.body = Cube(
            Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
            ctrPoint = np.array([0,0,offsetHeight+0.8]),
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
        
        if position.lower() == 'sitting':
            self.chair = Cube(
                Vector(np.array([1,np.deg2rad(90),np.deg2rad(0)])),
                ctrPoint = np.array([0,0,0.2]),
                dimensions = [0.4,0.4,0.4],
                reflectivities={'p0':reflectivities['furniture'],
                                'p1':reflectivities['shirt'],
                                'p2':reflectivities['shirt'],
                                'p3':reflectivities['furniture'],
                                'p4':reflectivities['furniture'],
                                'p5':reflectivities['furniture']}
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
                reflectivities={'p0':reflectivities['furniture'],
                                'p1':reflectivities['furniture'],
                                'p2':reflectivities['furniture'],
                                'p3':reflectivities['furniture'],
                                'p4':reflectivities['furniture'],
                                'p5':reflectivities['furniture']}
            )
            self.table = self.table.rotate(direction,np.array([0,0,1]));
            self.table = self.table.translate(np.array([*loc,0])+self.table.ctrPoint);

            Rd = getRodriguesMtx(direction,np.array([0,0,1]))
            shiftedVect = (Rd@np.array([1,0,0]).reshape(3,1)).reshape(-1)
            self.table = self.table.translate(shiftedVect+self.table.ctrPoint)
            
            # The table is modeled a simple plane
            self.table = self.table.listPlanes[2]
            
            self.leg = None
            
        elif position.lower() == 'standing':
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
        
            self.chair, self.table = None,None
            
        
        self.ctrPoint = np.array([0,0,1.2])+np.array([*loc,0])
        self.normalVect = self.head.normalVect;
        
        self.listPlanes = self.getPartition(Ps=1)
        
        self.mode = None
        self.pd0 = None
        self.pd1 = None
        self.led0 = None
        self.led1 = None
        
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
            
        if mode.lower() == 'rx':
            self.mode = mode.lower()
            # PD
            self.pd0 = BareDetector(polar,azimuth,loc,area=area,FoV=FoV)
            self.pd1 = BareDetector(polar,azimuth,loc+np.array([0.05,0,0]),area=area,FoV=FoV)
        elif mode.lower() == 'tx':
            self.mode = mode.lower()
            # LED
            self.led0 = PointSource(polar,azimuth,loc,m=m)
            self.led1 = PointSource(polar,azimuth,loc+np.array([0.05,0,0]),m=m)
        else:
            raise ValueError("'mode' is not recognized! \
            Should be either 'tx' or 'rx'")
        
    def getPartition(self,Ps=1,delta=None):
        """ Partition all objects.
        
        Parameters
        ----------
        Ps: list
            List of number of partition of each side. 
        delta: list
            Define the partition based on partition lengths

        Returns
        -------
        list:
            A list of partitioned planes. Each plane is an instant of 
            RectPlane_py.
        
        """
        
        planes=[]
        if delta == None:
            planes.append(self.head.getPartition(Ps=Ps))
            planes.append(self.body.getPartition(Ps=Ps))
            if self.leg:
                planes.append(self.leg.getPartition(Ps=Ps))
            if self.chair:
                planes.append(self.chair.getPartition(Ps=Ps))
            if self.table:
                planes.append(self.table.getPartition(Ps=Ps))
        else:
            if delta > 0.15:
                delta = 0.15
            planes.append(self.head.getPartition(delta=delta))
            planes.append(self.body.getPartition(delta=delta))
            if self.leg:
                planes.append(self.leg.getPartition(delta=delta))
            if self.chair:
                planes.append(self.chair.getPartition(delta=delta))
            if self.table:
                planes.append(self.table.getPartition(delta=delta))
        
        planes = list(flatten(planes))
        return planes