from owcsimpy.geoutils.mobility.randommobility_py import *
from owcsimpy.geoutils.mobility.randomorientation_py import *

#############################################################
import numpy as np

from owcsimpy.geoutils.mobility.randommobility_py import RandomWaypoint,TruncatedLevyWalk,RandomDirection
from owcsimpy.geoutils.mobility.randomorientation_py import RandOrientationArdimas as RandomOrientation
from owcsimpy.misc import flatten

def getArrayMobility(randomMobModel = 'rwp', activity = 'walking', dimensions = (0,1), velocityMinMax = (0.1,1), waitTimeMax = 1., excludedArea = 1., minSamples=30):
    """ Obtain arrays of locations, velocities, time instances, directions, and 
    changes of polar and azimuth angles of orientation of the UE.

    Parameters
    ----------
    randomMobModel: {'rwp','rd','lw'}
        'rwp': random waypoint
        'rd' : random direction
        'lw' : Levy-walk
    activity: {'sitting','standing','walking'}
        self-explanotory
    dimensions: tuple(2,)
        Areas of mobility of the human. 
        Ex: dimensions = (roomLength,roomWidth)
    velocityMinMax: tuple(2,)
        The minimum and maximum of velocities.
    waitTimeMax: double
        The maximum waiting time in seconds.
    excludedArea: double
        Area that will be excluded during simulation.
        Currently, it is as simple as a circle centered at the origin.
        Therefore, it is now only suitable for the simple office environment.
    minSamples: double
        The minimum number of generated samples. Default is 30.


    """

    assert randomMobModel.lower() in {'rwp','rd','lw'}
    assert activity.lower() in {'sitting','standing','walking'}
    assert len(dimensions) == len(velocityMinMax) == 2
    assert waitTimeMax >= 0
    assert excludedArea >= 0

    if activity.lower() == 'walking':
        if randomMobModel.lower() == 'rwp':
            rmm = RandomWaypoint(1, dimensions=dimensions, velocity=velocityMinMax, wt_max=waitTimeMax)
        elif randomMobModel.lower() == 'rd':
            rmm = RandomDirection(1, dimensions=dimensions, velocity=velocityMinMax, wt_max=waitTimeMax)
        elif randomMobModel.lower() == 'lw': 
            rmm = TruncatedLevyWalk(1, dimensions=dimensions, WT_MAX=waitTimeMax)

        # Make an iterable object from rmm
        rmmIter = iter(rmm)

        while True:
            arrayLoc = []
            arrayVel = []
            arrayTime = [0]

            while True:
                loc = next(rmmIter)
                
                if np.linalg.norm(loc) < excludedArea:
                    break
                else:
                    arrayLoc.append(np.array(loc))
                    arrayVel.append(np.array(rmm.velocity))
                    
            
            for i in range(len(arrayLoc)-1):
                if arrayVel[i] == 0:
                    deltaTime = 0
                else:
                    deltaTime = np.linalg.norm(arrayLoc[i+1]-arrayLoc[i])/arrayVel[i]
                    
                arrayTime.append(arrayTime[i]+deltaTime)

            arrayVel = np.array(list(flatten(arrayVel)))
            arrayTime = np.array(list(flatten(arrayTime)))
            arrayLoc = np.array(list(flatten(arrayLoc))).reshape((-1,2))

            if arrayLoc.shape[0] >= minSamples:
                break

        arrayDir = []
        for i in range(arrayLoc.shape[0]-1):
            arrayDir.append(getDirection(arrayLoc[i],arrayLoc[i+1]))

        arrayDir.append(arrayDir[-1])
        arrayDir = np.array(arrayDir)

        ro = RandomOrientation('walking')
        arrayTheta0,arrayOmega0 = ro.genSamples(arrayTime)  
    else: # {'sitting','standing'}
        ro = RandomOrientation(activity)
        arrayTime = np.arange(minSamples)*0.001; # 1 ms resolution
        arrayTheta0,arrayOmega0 = ro.genSamples(arrayTime) 

        arrayLoc, arrayVel, arrayDir = [],[],[]


    return {'arrayLoc':arrayLoc,'arrayVel':arrayVel,'arrayTime':arrayTime,'arrayDir':arrayDir,'arrayTheta0':arrayTheta0,'arrayOmega0':arrayOmega0}


def getAngle(vector):
    vector = list(vector)
    
    return np.arctan2(vector[1],vector[0])

def getDirection(vector0,vector1):
    vector = vector1-vector0
    
    return getAngle(vector)

