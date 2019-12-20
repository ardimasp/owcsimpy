import numpy as np

from owcsimpy.geoobjects.models.simpleofficeenv_py import SimpleOfficeEnv_py as SimpleOfficeEnv

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


