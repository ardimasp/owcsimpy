# examples:
# python datasetSimpleOffice.py --config-name standing_downlink_01 --mode cir --min-index 0 --max-index 0


__version__ = "3.7.3"

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import speed_of_light
from tqdm import tqdm
import pickle
import logging
logging.getLogger().setLevel(logging.INFO)

import pprint
import json
import sys
import os
import string
from time import sleep

from owcsimpy.geoobjects.models.typicalsimpleofficeenv_py import TypicalSimpleOfficeEnv_py as TypicalOfficeEnv
from owcsimpy.cir.timedomaincir import TimeDomainCIR

from owcsimpy.geoutils.mobility import getArrayMobility 

roomLength,roomWidth,roomHeight=4,3,3

__usage__ = """Usage: datasetSimpleOffice [options] ...
Collections of scripts related to the dataset for the simple office scenario.

Options:
  -n                            the name of the working, default = SimpleOffice
  --mode                        generate config, cir files, of animation, mode = ('config','cir','animation','matlab')
  --config-name                 the name of configuration
  --config-list                 list of the available configuration name
  --config-generate             the number of config samples generated
  --config-filename             the filename of generated config
  --number-realizations         the number of realizations
  --min-index                   set the minimum index, default = 0
  --max-index                   set the maximum index, default = depends on the config file
  -h, --help                    show this help message


"""

class ArgumentError(Exception):
    pass

def searchListDicts(listIn,nameConfig):

    for idx,item in enumerate(listIn):
        if item.get('name').lower()==nameConfig:
            return idx

    return -1

def console():
    try:
        outputParsing = parse_args()

    except(ArgumentError):
        sys.stderr.write(__usage__)

    configDirName = outputParsing.get('configDirName')
    datasetDirName = outputParsing.get('datasetDirName')
    configName = outputParsing.get('configName')

    with open(configDirName + 'config.json') as json_file:
        json_data = json.load(json_file)

    listConfigNames = []
    for data in json_data:
        listConfigNames.append(data.get('name'))

    if not( configName in listConfigNames):
        logging.error('List of config names: %s', listConfigNames)
        raise SystemExit

    if not os.path.exists(configDirName + configName):
        logging.info("Create a folder named '%s' ", configDirName + configName)
        os.makedirs(configDirName + configName + '/',exist_ok=True)

    if not os.path.exists(datasetDirName + configName):
        logging.info("Create a folder named '%s' ", datasetDirName + configName)
        os.makedirs(datasetDirName + configName + '/',exist_ok=True)

    idx = searchListDicts(json_data,configName)

    configContent = json_data[idx]

    filename = configDirName + configName + '/' + outputParsing.get('configFilename')
    datasetname = datasetDirName + configName + '/' + outputParsing.get('datasetFilename')
    
    #################################################

    if outputParsing.get('mode').lower() == 'config':

        if os.path.exists(filename):
            answer = None
            while answer not in ("y", "n"):
                answer = input("# '%s' exists. Do you want to replace it? [y/n]: " % filename)
                if answer[0].lower() == "y":
                     answer = 'y'
                elif answer[0].lower() == "n":
                     raise SystemExit
                else:
                    print("> Please enter yes or no.")

        ################################################
        if configContent.get('mode').lower() == 'walking':

            randomMobModel = configContent.get('mobility_model')
            activity = 'walking'
            dimensions = (roomLength-0.2,roomWidth-0.2)
            velocityMinMax = (configContent.get('velocity_min'),configContent.get('velocity_max'))
            waitTimeMax = configContent.get('wait_time_max')
            excludedArea = configContent.get('excluded_area')

            listTraj = []
            for _ in range(outputParsing.get('numRealizations')):
                traj = getArrayMobility(
                    randomMobModel=randomMobModel,
                    activity=activity,
                    dimensions=dimensions,
                    velocityMinMax=velocityMinMax,
                    waitTimeMax=waitTimeMax,
                    excludedArea=excludedArea) 

                listTraj.append(traj)         

        elif configContent.get('mode').lower() == 'standing':

            minSamples = 5000; # 1 sample = 1 ms 
            
            listTraj = []
            for _ in range(outputParsing.get('numRealizations')):
                listTraj.append(getArrayMobility(activity = 'standing',minSamples=minSamples))

        elif configContent.get('mode').lower() == 'sitting':
            
            minSamples = 5000; # 1 sample = 1 ms 
            
            listTraj = []
            for _ in range(outputParsing.get('numRealizations')):
                listTraj.append(getArrayMobility(activity = 'sitting',minSamples=minSamples))


        else:
            info.error('mode: %s is not registered',configContent.get('mode'))
            raise SystemExit

        with open(filename,'wb') as pickle_file:
            pickle.dump(listTraj,pickle_file)

        logging.info("'%s' is created",filename)
        
    #################################################
    elif outputParsing.get('mode').lower() in {'cir','animation'}:

        # if os.path.exists(datasetname):
        #     answer = None
        #     while answer not in ("y", "n"):
        #         answer = input("# '%s' exists. Do you want to replace it? [y/n]: " % filename)
        #         if answer[0].lower() == "y":
        #              answer = 'y'
        #         elif answer[0].lower() == "n":
        #              raise SystemExit
        #         else:
        #             print("> Please enter yes or no.")

        if not os.path.exists(filename):
            logging.error("'%s' does not exist. Run mode='config' first.",filename)
            raise SystemExit

        with open(filename,'rb') as pickle_file:
            listTraj = pickle.load(pickle_file)

        if configContent.get('mode').lower() == 'walking':
            
            position = 'standing'
            furnitureConfig = 'cornerYFacingWall'
            humanLoc = [roomLength/2,roomWidth/2]
            humanDirection = np.deg2rad(45)
            activity = configContent.get('activity')
            if configContent.get('downlink_uplink').lower() == 'downlink':
                modeRxOrTx = 'rx' 
            elif configContent.get('downlink_uplink').lower() == 'uplink':
                modeRxOrTx = 'tx'
            else:
                info.error('downlink_uplink: %s is not registered',configContent.get('downlink_uplink'))
                raise SystemExit

            office = TypicalOfficeEnv(
                roomDim=[roomLength,roomWidth,roomHeight],
                humanLoc=humanLoc,humanDirection=humanDirection,
                mode=modeRxOrTx,
                activity=activity,
                position=position,
                furnitureConfig=furnitureConfig
            )

            if outputParsing.get('mode').lower() == 'animation':
                print('Press Ctrl+C to exit. Close the figure window for the next animation.')
                print("There are %d animations" % len(listTraj))
                for i,traj in enumerate(listTraj):

                    logging.info("Animation: %d",i)

                    arrayLoc = traj.get('arrayLoc')
                    arrayVel = traj.get('arrayVel')
                    arrayTime = traj.get('arrayTime')
                    arrayDir = traj.get('arrayDir')
                    arrayTheta0 = traj.get('arrayTheta0')
                    arrayOmega0 = traj.get('arrayOmega0')

                    # Since the time resolution is 1 ms, we need to slice it so that the animation won't take too long
                    numSlice = 100 # in the order of (1 -> 1 ms, 10 -> 10 ms, etc.) 

                    arrayLocSliced = arrayLoc[0::numSlice]
                    arrayVelSliced = arrayVel[0::numSlice]
                    arrayTimeSliced = arrayTime[0::numSlice]
                    arrayDirSliced = arrayDir[0::numSlice]
                    arrayTheta0Sliced = arrayTheta0[0::numSlice]
                    arrayOmega0Sliced = arrayOmega0[0::numSlice]

                    fig,ax,anim = office.getAnimation(
                        arrayLocSliced,arrayDirSliced,arrayTheta0Sliced,arrayOmega0Sliced,arrayTimeSliced,interval=1)
                    fig.set_size_inches(5, 4)

                    plt.show()

        elif configContent.get('mode').lower() == 'standing':
            
            position = 'standing'
            furnitureConfig = 'cornerYFacingWall'
            humanLoc = [configContent.get('loc_x'),configContent.get('loc_y')]
            humanDirection = np.deg2rad(configContent.get('direction_degree'))
            activity = configContent.get('activity')
            if configContent.get('downlink_uplink').lower() == 'downlink':
                modeRxOrTx = 'rx' 
            elif configContent.get('downlink_uplink').lower() == 'uplink':
                modeRxOrTx = 'tx'
            else:
                info.error('downlink_uplink: %s is not registered',configContent.get('downlink_uplink'))
                raise SystemExit

            office = TypicalOfficeEnv(
                roomDim=[roomLength,roomWidth,roomHeight],
                humanLoc=humanLoc,humanDirection=humanDirection,
                mode=modeRxOrTx,
                activity=activity,
                position=position,
                furnitureConfig=furnitureConfig
            )

            if outputParsing.get('mode').lower() == 'animation':
                print('Press Ctrl+C to exit. Close the figure window for the next animation.')
                print("There are %d animations" % len(listTraj))
                for i,traj in enumerate(listTraj):

                    logging.info("Animation: %d",i)

                    arrayLoc = traj.get('arrayLoc')
                    arrayVel = traj.get('arrayVel')
                    arrayTime = traj.get('arrayTime')
                    arrayDir = traj.get('arrayDir')
                    arrayTheta0 = traj.get('arrayTheta0')
                    arrayOmega0 = traj.get('arrayOmega0')

                    # farrayLoc, arrayVel, and arrayDir are empty lists if activity = 'standing'
                    arrayLoc = np.ones((len(arrayTime),2))*office.human.ctrPoint[0:2] # static
                    arrayVel = np.zeros((len(arrayTime),)) # zero velocity
                    arrayDir = np.ones((len(arrayTime),))*humanDirection # static direction


                    # Since the time resolution is 1 ms, we need to slice it so that the animation won't take too long
                    numSlice = 10 # in the order of (1 -> 1 ms, 10 -> 10 ms, etc.) 

                    arrayLocSliced = arrayLoc[0::numSlice]
                    arrayVelSliced = arrayVel[0::numSlice]
                    arrayTimeSliced = arrayTime[0::numSlice]
                    arrayDirSliced = arrayDir[0::numSlice]
                    arrayTheta0Sliced = arrayTheta0[0::numSlice]
                    arrayOmega0Sliced = arrayOmega0[0::numSlice]


                    fig,ax,anim = office.getAnimation(
                        arrayLocSliced,arrayDirSliced,arrayTheta0Sliced,arrayOmega0Sliced,arrayTimeSliced,interval=1)
                    fig.set_size_inches(5, 4)

                    plt.show()

        elif configContent.get('mode').lower() == 'sitting':
            
            position = 'sitting'
            furnitureConfig = 'cornerYFacingWall'
            humanLoc = [configContent.get('loc_x'),configContent.get('loc_y')]
            humanDirection = np.deg2rad(configContent.get('direction_degree'))
            activity = configContent.get('activity')
            if configContent.get('downlink_uplink').lower() == 'downlink':
                modeRxOrTx = 'rx' 
            elif configContent.get('downlink_uplink').lower() == 'uplink':
                modeRxOrTx = 'tx'
            else:
                info.error('downlink_uplink: %s is not registered',configContent.get('downlink_uplink'))
                raise SystemExit

            office = TypicalOfficeEnv(
                roomDim=[roomLength,roomWidth,roomHeight],
                humanLoc=humanLoc,humanDirection=humanDirection,
                mode=modeRxOrTx,
                activity=activity,
                position=position,
                furnitureConfig=furnitureConfig
            )

            if outputParsing.get('mode').lower() == 'animation':
                print('Press Ctrl+C to exit. Close the figure window for the next animation.')
                print("There are %d animations" % len(listTraj))
                for i,traj in enumerate(listTraj):

                    logging.info("Animation: %d",i)

                    arrayLoc = traj.get('arrayLoc')
                    arrayVel = traj.get('arrayVel')
                    arrayTime = traj.get('arrayTime')
                    arrayDir = traj.get('arrayDir')
                    arrayTheta0 = traj.get('arrayTheta0')
                    arrayOmega0 = traj.get('arrayOmega0')

                    # farrayLoc, arrayVel, and arrayDir are empty lists if activity = 'standing'
                    arrayLoc = np.ones((len(arrayTime),2))*office.human.ctrPoint[0:2] # static
                    arrayVel = np.zeros((len(arrayTime),)) # zero velocity
                    arrayDir = np.ones((len(arrayTime),))*humanDirection # static direction


                    # Since the time resolution is 1 ms, we need to slice it so that the animation won't take too long
                    numSlice = 10 # in the order of (1 -> 1 ms, 10 -> 10 ms, etc.) 

                    arrayLocSliced = arrayLoc[0::numSlice]
                    arrayVelSliced = arrayVel[0::numSlice]
                    arrayTimeSliced = arrayTime[0::numSlice]
                    arrayDirSliced = arrayDir[0::numSlice]
                    arrayTheta0Sliced = arrayTheta0[0::numSlice]
                    arrayOmega0Sliced = arrayOmega0[0::numSlice]


                    fig,ax,anim = office.getAnimation(
                        arrayLocSliced,arrayDirSliced,arrayTheta0Sliced,arrayOmega0Sliced,arrayTimeSliced,interval=1)
                    fig.set_size_inches(5, 4)

                    plt.show()

        else:
            info.error('mode: %s is not registered',configContent.get('mode'))
            raise SystemExit

        if outputParsing.get('mode').lower() == 'cir':

            minIdx = outputParsing.get('minIdx')
            maxIdx = outputParsing.get('maxIdx')

            if maxIdx > len(listTraj):
                maxIdx = len(listTraj)-1

            if minIdx > maxIdx:
                print(minIdx)
                print(maxIdx)
                raise ArgumentError

            if not (0 <= minIdx <= 999999) and not (0 <= maxIdx <= 999999):
                raise ArgumentError

            for idx in tqdm(range(minIdx,maxIdx+1)):
                
                traj = listTraj[idx]

                arrayLoc = traj.get('arrayLoc')
                arrayVel = traj.get('arrayVel')
                arrayTime = traj.get('arrayTime')
                arrayDir = traj.get('arrayDir')
                arrayTheta0 = traj.get('arrayTheta0')
                arrayOmega0 = traj.get('arrayOmega0')

                if configContent.get('mode').lower() in {'sitting','standing'}:
                    # farrayLoc, arrayVel, and arrayDir are empty lists if activity = 'standing'
                    arrayLoc = np.ones((len(arrayTime),2))*office.human.ctrPoint[0:2] # static
                    arrayVel = np.zeros((len(arrayTime),)) # zero velocity
                    arrayDir = np.ones((len(arrayTime),))*humanDirection # static direction

                listResult = []

                for idy in tqdm(range(0,len(arrayTime))):

                    office.updateHuman(arrayLoc[idy],arrayDir[idy],arrayTheta0[idy],arrayOmega0[idy])

                    timeSampling = 1e-9
                    tdcir = TimeDomainCIR(timeSampling=timeSampling)

                    if configContent.get('downlink_uplink').lower() == 'downlink':
                        ht_los, ht_diff = tdcir.calc(office.led,office.human.pd,office.blockingObj,office.room,
                                                           numReflections=3,partitionDist=speed_of_light*timeSampling)
                    elif configContent.get('downlink_uplink').lower() == 'uplink':
                        ht_los, ht_diff = tdcir.calc(office.human.led,office.pd,office.blockingObj,office.room,
                                                           numReflections=3,partitionDist=speed_of_light*timeSampling)
                    else:
                        info.error('downlink_uplink: %s is not registered',configContent.get('downlink_uplink'))
                        raise SystemExit

                    listResult.append({
                        'ht_los':ht_los,
                        'ht_diff':ht_diff,
                        't':tdcir.t})

                with open(datasetname+str(idx)+'.pickle','wb') as pickle_file:
                    pickle.dump(listResult,pickle_file)

                logging.info("'%s' is created",datasetname+str(idx)+'.pickle')

    elif outputParsing.get('mode').lower() == 'matlab':
        
    else:
        logging.error("The entry of 'mode' is wrong")
        sys.stderr.write(__usage__)
        raise SystemExit




def parse_args():

    args = sys.argv[1:]

    # initialization
    rootDir = './../dataset/'
    workingDir = 'SimpleOffice'
    configName = ''
    numConfigGen = 0
    mode = 'config'
    configFilename = 'samples.pickle'
    datasetFilename = 'samples'
    numRealizations = 1000

    minIdx = 0;
    maxIdx = int(1e20-1); # I don't think we will have this many samples

    showListConfig = False

    # check wheter the dataset folder exists
    # if not create one
    if not os.path.exists(os.getcwd() + '/../dataset'):
        os.makedirs(os.getcwd() + '/../dataset')

    i = 0
    while i < len(args):
        op = args[i]
        if op in ("-h","--help"):
            raise ArgumentError
        elif op in ("-n"):
            workingDir = str(args[i+1])
            i += 1
        elif op in ("--config-list"):
            showListConfig = True
        elif op in ("--config-name"):
            configName = str(args[i+1])
            i += 1
        elif op in ("--config-generate"):
            numConfigGen = int(args[i+1])
            i += 1
        elif op in ("--config-filename"):
            configFilename = str(args[i+1])
            i += 1
        elif op in ("--mode"):
            mode = str(args[i+1])
            i += 1
        elif op in ("--number-realizations"):
            numRealizations = int(args[i+1])
            i += 1
        elif op in ("--min-index"):
            minIdx = int(args[i+1])
            i += 1
        elif op in ("--max-index"):
            maxIdx = int(args[i+1])
            i += 1
        i += 1

    datasetDirName = rootDir + workingDir + '/dataset/'
    configDirName = rootDir + workingDir + '/config/'

    # create datasetDirName folder if not exist
    if not os.path.exists(datasetDirName):
        logging.info("Creating folder '%s'" % datasetDirName)
        os.makedirs(datasetDirName,exist_ok=True)

    if not os.path.exists(configDirName):
        logging.info("Creating folder '%s'" % configDirName)
        os.makedirs(configDirName,exist_ok=True)

    if showListConfig:
        try:
            with open(configDirName+'config.json') as json_file:
                json_data = json.load(json_file)

            # pprint.pprint(json_data)

            print('# List of config names:')
            for data in json_data:
                print('## '+ data.get('name'))

        except:
            logging.error("There is something wrong with '%s'",configDirName+'config.json')
            raise SystemExit

        raise SystemExit

    return {'datasetDirName':datasetDirName,
            'configDirName':configDirName,
            'configName':configName,
            'configFilename':configFilename,
            'datasetFilename':datasetFilename,
            'numConfigGen':numConfigGen,
            'mode':mode,
            'minIdx':minIdx,
            'maxIdx':maxIdx,
            'numRealizations':numRealizations}

if __name__ == "__main__":
    try:
        
        console()

    except (KeyboardInterrupt, SystemExit):
        pass
