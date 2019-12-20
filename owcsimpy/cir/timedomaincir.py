import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.constants import speed_of_light

from warnings import warn

# from owcsimpy.geoobjects.models.roomcube_py import RoomCube_py as Room
# from owcsimpy.geoobjects.models.humancube_py import HumanCube_py as Human
# from owcsimpy.geoobjects.models.pointsource_py import PointSource_py as PointSource
# from owcsimpy.geoobjects.models.baredetector_py import BareDetector_py as BareDetector
from owcsimpy.cir.cirutils import calcCIRTimeDom
from owcsimpy.misc import flatten

class TimeDomainCIR(object):
    """ CIR calculation by means of the time-domain approach [1].
    
    Parameters
    ----------
    timeSampling: double (optional=0.1 ns)

    Attributes
    ----------
    timeSampling: double
    f: list
        List of frequency bins
    ht_los: array-like
        Time-domain CIR for LOS
    ht_diff: array-like
        Time-domain CIR for diffuse channel
    Hf_los: array-like
        Freq response for LOS channel
    Hf_diff: array-like
        Freq response for diffuse channel


    Notes
    -----
    The implementation and the variables' name follow [2]. IMHO, [2] is more 
    readable compared to [3]. 

    Refs:
        [1] J. B. Carruthers and P. Kannan, "Iterative site-based modeling for wireless infrared channels," in IEEE Transactions on Antennas and Propagation, vol. 50, no. 5, pp. 759-765, May 2002. 
        
        [2] https://github.com/UCaNLabUMB/CandLES
        
        [3] http://iss.bu.edu/bwc/irsimit/

    """

    def __init__(self,timeSampling=1e-10):

        self.timeSampling = timeSampling
        self.t = []
        self.f = []
        # time domain
        self.ht_los = []
        self.ht_diff = []
        # frequency domain
        self.Hf_los = []
        self.Hf_diff = []

    def calc(self,LEDs,PDs,blockingObj,reflectingObj,
        partitionDist=speed_of_light*1e-9,numReflections=3,verbose=False):
        """ Calcute the CIR.
        
        Parameters
        ----------
        LEDs: list
            List of LEDs. But, currently supports 1 LED.
        PDs: list
            List of PDs. Currently supports 1 PD.
        blockingObj: list
            List of blocking objects. 
        reflectingObj: list
            List of reflecting objects, for example a room.
            The naming is due to we might need a genera case when
            we can do the assumptions of infinite room such that 
            what matters most are the ceiling and the floor.
            Note that the blockingObj will also be treated as a 
            relfecting object.
        partitionDist: double, optional
            Delta distance based on which we partition a plane 
            (the default value is 1ns * c, c: speed of light).
            The default value is c*timesampling, timesampling is 
            assumed to 1ns by default. See Schulze's paper. 
        numReflections: inf or integer, optional
            Denoting the number of reflections (the default is 3).


        """


        # Transform into list
        if not isinstance(LEDs,list):
            LEDs = [LEDs]
        
        if not isinstance(PDs,list):
            PDs = [PDs]

        if not isinstance(blockingObj,list):
            blockingObj = [blockingObj]

        if not isinstance(reflectingObj,list):
            reflectingObj = [reflectingObj]

        # Currently support 1 LED and 1 PD
        assert len(LEDs)==1 and len(PDs)==1

        # Partition assuming that all obj in reflectingObj
        # has the getPartition method
        reflectingPlanes=[]
        for obj in reflectingObj:
            reflectingPlanes.append(obj.getPartition(delta=partitionDist))

        # don't forget to include blockingObj as reflecting objects
        for obj in blockingObj:
            reflectingPlanes.append(obj.getPartition(delta=partitionDist))

        # flatten the list
        reflectingPlanes = list(flatten(reflectingPlanes))

        # Partition blockingObj
        blockingPlanes=[]
        for obj in blockingObj:
            # The number of partition on each side is one
            blockingPlanes.append(obj.getPartition(Ps=1))

        blockingPlanes = list(flatten(blockingPlanes))

        # Get simple planes
        simpleReflectingPlanes = [plane.getSimplePlane() for plane in reflectingPlanes]
        simpleBlockingPlanes = [plane.getSimplePlane() for plane in blockingPlanes]
        emitters = [led.getSimplePointSource() for led in LEDs]
        collectors = [pd.getSimpleBareDetector() for pd in PDs]

        # Get the farthest vertix from origin
        dummyPlanes = []
        for obj in reflectingObj:
            dummyPlanes.append(obj.getPartition(Ps=1))

        dummyPlanes = list(flatten(dummyPlanes))

        # Assuming each dummyPlane has the vertixes attributes
        # Collect all vertices
        verts = np.array(list(flatten([plane.verts for plane in dummyPlanes]))).reshape(-1,3)

        # Calculate norm relative to origin
        listOfNorm = list(map(lambda v: np.linalg.norm(v), verts))

        outerVert = verts[listOfNorm.index(max(listOfNorm))]

        if verbose:
            print("> Calculating CIR......(wait)")
        
        self.ht_los,self.ht_diff = calcCIRTimeDom(self.timeSampling,numReflections,outerVert,
            emitters[0],collectors[0],simpleReflectingPlanes,simpleBlockingPlanes)
        
        if verbose:
            print("> Finish calculating CIR :)")

        self.t = np.arange(self.ht_los.size)*self.timeSampling

        return self.ht_los,self.ht_diff

    def transform(self,NFFT=[]):
        """ Transform to freq. domain.

        Returns
        -------
        f: array-like
            Frequency bins.
        Hf_los: array-like
            Freq. response of LOS channel.
        Hf_diff: array-like
            Freq. response of diffuse channel 
        
        """
        
        freqSampling = 1/self.timeSampling
        L = self.ht_los.size
        if NFFT:
            N = NFFT
        else:
            N = L+1 # DFT point

        f = np.arange(N)*freqSampling/N # frequency bins
        
        Hlos = np.fft.fft(self.ht_los.reshape(-1),n=N)
        Hdiff = np.fft.fft(self.ht_diff.reshape(-1),n=N)

        # self.f, self.Hf_los, self.Hf_diff = f,Hlos,Hdiff
        self.f = f
        self.Hf_los, self.Hf_diff = Hlos,Hdiff

        return self.Hf_los, self.Hf_diff

    def plot(self,domain='time'):
        """ Plot the CIR.

        Parameters
        ----------
        domain: {'frequency','time'}, optional
            Domain to plot (the default value is 'time').
        
        Raises
        ------
        Exception:
            If the the method transform hasn't been invoked yet.

        """

        fig, ax = plt.subplots()

        if domain == 'frequency':
            if self.Hf_los == [] or self.Hf_diff == []:
                raise Exception(
                    """The freq-domain transformation hasn't been 
                    carried out yet!""")
            else:

                ax.plot(self.f/1e6,10*np.log10(np.abs(self.Hf_los+self.Hf_diff)**2))
                ax.set_xlabel(r"$f$ [MHz]"); 
                ax.set_xlim([self.f[0]/1e6,self.f[int(len(self.f)/2)]/1e6])
                ax.set_ylabel(r"$\vert H(f) \vert^2$ [dB]");
                
        
        elif domain == 'time':

            # t = np.arange(self.ht_los.size)*self.timeSampling
            ax.plot(self.t,(self.ht_los+self.ht_diff)/self.timeSampling)
            ax.set_xlabel(r"$t$ [s]"); ax.set_ylabel(r"$h(t)$ [1/s]");
            



