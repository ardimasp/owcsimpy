import numpy as np

from owcsimpy.geoobjects.models.humanwithactivity_py import HumanWithActivity_py as Human
from owcsimpy.geoobjects.models.simpleofficeenv_py import SimpleOfficeEnv_py as SimpleOfficeEnv
from owcsimpy.misc import flatten
from owcsimpy.geoutils.draw import draw
from owcsimpy.geoobjects.models.pointsource_py import PointSource_py as PointSource
from owcsimpy.geoobjects.models.baredetector_py import BareDetector_py as BareDetector
import mpl_toolkits.mplot3d.art3d as art3d

from JSAnimation.IPython_display import display_animation
from matplotlib import animation

class TypicalSimpleOfficeEnv_py(SimpleOfficeEnv):
    """ A typical realization of a simple office. 
    
    The chair and the table are located in the corner of the room.
    
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
    furnitureConfig: {'cornerXFacingWall','cornerXNotFacingWall','cornerYFacingWall','cornerYNotFacingWall'}
        The direction of the chair and the table.
        
    Attributes
    ----------
    See
    :class:`~owcsimpy.geoobjects.model.simpleofficeenv_py.SimpleOfficeEnv_py`
    
    """
    
    def __init__(self,
                 roomDim,
                 humanLoc=None,humanDirection=0,
                 mode='rx',
                 activity='calling',
                 position='standing',
                 furnitureConfig='cornerXFacingWall'
                ):

        chairLoc,chairDirection = None,0

        if position.lower()=='standing':
            if furnitureConfig == 'cornerXFacingWall':
                chairLoc=[0.6,1.45]
                chairDirection=np.deg2rad(270)
            elif furnitureConfig == 'cornerXNotFacingWall':
                chairLoc=[0.6,0]
                chairDirection=np.deg2rad(90)
            elif furnitureConfig == 'cornerYFacingWall':
                chairLoc=[1.45,0.6]
                chairDirection=np.deg2rad(180)
            elif furnitureConfig == 'cornerYNotFacingWall':
                chairLoc=[0,0.6]
                chairDirection=np.deg2rad(0)
        elif position.lower()=='sitting':
            # When sitting, the configuration of human will be overridden
            if furnitureConfig == 'cornerXFacingWall':
                humanLoc=[0.6,1.45+0.075]
                humanDirection=np.deg2rad(270)
            elif furnitureConfig == 'cornerXNotFacingWall':
                humanLoc=[0.6,0+0.075]
                humanDirection=np.deg2rad(90)
            elif furnitureConfig == 'cornerYFacingWall':
                humanLoc=[1.45+0.075,0.6]
                humanDirection=np.deg2rad(180)
            elif furnitureConfig == 'cornerYNotFacingWall':
                humanLoc=[0+0.075,0.6]
                humanDirection=np.deg2rad(0)
        
        super().__init__(
            roomDim=roomDim,
            humanLoc=humanLoc,humanDirection=humanDirection,
            chairLoc=chairLoc,chairDirection=chairDirection,
            mode=mode,
            activity=activity,
            position=position
        )

        if mode.lower() == 'rx':
            self.pdPolarInit = self.human.pd.polar
            self.pdAzimuthInit = self.human.pd.azimuth
        elif mode.lower() == 'tx':
            self.ledPolarInit = self.human.led.polar
            self.ledAzimuthInit = self.human.led.azimuth

    def updateHuman(self,newLoc,newDir,theta0,omega0):

        """ Update the human's configuration for animation purposes. 
        
        Parameters
        ----------
        newLoc
            New location
        newDir
            New direction
        theta0
            New delta of the polar angle of UE
        omega0
            New delta of the azimuth angle of UE
            
        
        """

        # Human 
        self.human = Human(loc=newLoc,direction=newDir,
                          mode=self.human.mode,activity=self.human.activity,position=self.human.position)

        if self.human.mode.lower() == 'rx':
            loc = self.human.pd.loc
            area = self.human.pd.area
            FoV = self.human.pd.FoV
            az = self.human.head.listPlanes[1]._RectPlane_py__normalVect.spherical[2] + np.pi
            self.human.pd = BareDetector(self.pdPolarInit+theta0,az+omega0,loc,area=area,FoV=FoV)
        elif self.human.mode.lower() == 'tx':
            loc = self.human.led.loc
            m = self.human.led.m
            az = self.human.head.listPlanes[1]._RectPlane_py__normalVect.spherical[2] + np.pi
            self.human.led = PointSource(self.ledPolarInit+theta0,az+omega0,loc,m=m)

        super().updateListPlanesAndCheckValidity()


    def drawToAnimation(self):
        
        """ Draw. Assuming that the human is standing.
        
        Parameters
        ----------
        
        
        """
        
    
        listPlanes = []
        if self.chair and self.table:
            listPlanes.append(self.chair.getPartition(Ps=1))
            listPlanes.append(self.table.getPartition(Ps=1))
        listPlanes = list(flatten(listPlanes))
        
        alphas=[1 for plane in listPlanes]
        
        fig,ax = draw(xlim=[0,self.room.L],ylim=[0,self.room.W],zlim=[0,self.room.H],
                     planes=listPlanes,
                     enablevect=False,
                     alphas=alphas);

        fig, ax = draw(figure=fig,axes=ax,
             vectors=[self.led],colors='blue');
        

        return fig,ax

    def getAnimation(self,arrayLoc,arrayDir,arrayTheta0,arrayOmega0,arrayTime,interval=50):
        """ Create animation based on arrays of locations, direactions, and polar angles
        
        Parameters
        ----------
        arrayLoc: ndarray(N,2)
            Locations of the human over time
        arrayDir: list
            Directions of the human over time
        arrayTheta0: list
            Delta of polar angles of the over time
        arrayOmega0: list
            Delta of azimuth angles of the over time
        interval:
            Interval of animation in ms
            
        
        """


        fig,ax = self.drawToAnimation();

        label = ax.text(0,self.room.L,self.room.H,'t= 0s')

        colls = []

        for plane in self.human.head.listPlanes:
            verts = [list(map(lambda vert: (vert[0],vert[1],vert[2]), 
                    plane.verts))]
            col = art3d.Poly3DCollection(verts,edgecolor='k',facecolor='r')
            colls.append(ax.add_collection3d(col))

        for plane in self.human.body.listPlanes:
            verts = [list(map(lambda vert: (vert[0],vert[1],vert[2]), 
                    plane.verts))]
            col = art3d.Poly3DCollection(verts,edgecolor='k',facecolor='r')
            colls.append(ax.add_collection3d(col))
            
        for plane in self.human.leg.listPlanes:
            verts = [list(map(lambda vert: (vert[0],vert[1],vert[2]), 
                    plane.verts))]
            col = art3d.Poly3DCollection(verts,edgecolor='k',facecolor='r')
            colls.append(ax.add_collection3d(col))
            
        arrow_dir = ax.get_arrow(self.human.ctrPoint,self.human.normalVect,color='black')
        colls.append(ax.add_collection(arrow_dir))

        arrow_pd = ax.get_arrow(self.human.pd.ctrPoint,self.human.pd.normalVect,color='red')
        colls.append(ax.add_collection(arrow_pd))


        def animate(i):
    
            self.updateHuman(arrayLoc[i],arrayDir[i],arrayTheta0[i],arrayOmega0[i])
            
            text = 't = %.3f s' % (arrayTime[i])
            label.set_text(text)

            offset = 0;
            for j,plane in enumerate(self.human.head.listPlanes):
                verts = [list(map(lambda vert: (vert[0],vert[1],vert[2]), 
                        plane.verts))]
                colls[offset+j].set_verts(verts)
            
            offset = 6;
            for j,plane in enumerate(self.human.body.listPlanes):
                verts = [list(map(lambda vert: (vert[0],vert[1],vert[2]), 
                        plane.verts))]
                colls[offset+j].set_verts(verts)
            
            offset = 12;
            for j,plane in enumerate(self.human.leg.listPlanes):
                verts = [list(map(lambda vert: (vert[0],vert[1],vert[2]), 
                        plane.verts))]
                colls[offset+j].set_verts(verts)
            
            arrow = ax.get_arrow(self.human.ctrPoint,self.human.normalVect)
            colls[18].set_segments([arrow._segments3d[k,:,:] for k in range(3)])
            
            arrow = ax.get_arrow(self.human.pd.ctrPoint,self.human.pd.normalVect)
            colls[19].set_segments([arrow._segments3d[k,:,:] for k in range(3)])
            
            return colls

        anim = animation.FuncAnimation(fig, animate, frames=arrayDir.shape[0], interval=interval)

        return fig,ax,anim
        # display_animation(anim, default_mode='once')
