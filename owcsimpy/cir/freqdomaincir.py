import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.constants import speed_of_light

from warnings import warn

# from owcsimpy.geoobjects.models.roomcube_py import RoomCube_py as Room
# from owcsimpy.geoobjects.models.humancube_py import HumanCube_py as Human
# from owcsimpy.geoobjects.models.pointsource_py import PointSource_py as PointSource
# from owcsimpy.geoobjects.models.baredetector_py import BareDetector_py as BareDetector
from owcsimpy.misc import flatten
from owcsimpy.cir.cirutils import calcCIRFreqDom
from owcsimpy.cir.cirutils import getMtxFreqDom
from owcsimpy.cir.cirutils import calcCIRFreqDom_fromMtx


class FreqDomainCIR(object):
    """ CIR calculation by means of the frequency-domain approach [1].
    
    Parameters
    ----------
    N: int 
        Number of frequency bins
    freqSampling: double
        Frequency sampling
    
    Attributes
    ----------
    N: int
        Number of frequency bins
    freqSampling: double
        Frequency sampling
    timeSampling: double
        Time sampling
    f: array-like
        Frequency bins
    t: array-like
        Time bins
    Hf_los: array-like
        Freq response for LOS channel
    Hf_diff: array-like
        Freq response for diffuse channel
    ht_los: array-like
        Time-domain CIR for LOS
    ht_diff: array-like
        Time-domain CIR for diffuse channel

    Notes
    -----
    Refs: 
    [1] H. Schulze, "Frequency-Domain Simulation of the Indoor Wireless Optical Communication Channel," in IEEE Transactions on Communications, vol. 64, no. 6, pp. 2551-2562, June 2016.

    """

    def __init__(self,N,freqSampling):

        assert N > 0 and freqSampling > 0
        if not isinstance(N,int):
            warn("N is not int. N will be casted into int")

        self.N = np.int(N)
        self.freqSampling = freqSampling
        self.timeSampling = 1/freqSampling
        self.f = np.arange(N)*freqSampling/N
        self.t = self.timeSampling*np.arange(N)
        self.fi = []
        # Freq. response
        self.Hf_los = []
        self.Hf_diff = []
        # time-domain CIR
        self.ht_los = []
        self.ht_diff = []

        # Saved matrices and vectors
        self.H_f = []
        self.t_f = []
        self.r_f = []
        self.H_los = []

    def calc(self,LEDs,PDs,blockingObj,reflectingObj,
        partitionDist=speed_of_light*1e-9,numReflections=math.inf,verbose=False):
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
        partitionDist: double (optional)
            Delta distance based on which we partition a plane.
            The default value is c*timesampling, timesampling is 
            assumed to 1ns by default. See Schulze's paper.
        numReflections: inf or integer (optional)
            Denoting the number of reflections.


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

        if verbose:
            print("> Calculating CIR......(wait)")
        
        if len(self.H_f) != 0:
            # self.Hf_los,self.Hf_diff = calcCIRFreqDom_fromMtx(self.f,
            #     emitters[0],collectors[0],simpleReflectingPlanes,
            #     self.H_f,self.t_f,self.r_f,self.H_los,blockingplanes=simpleBlockingPlanes)
            self.Hf_los,self.Hf_diff = calcCIRFreqDom_fromMtx(self.f,
                self.Nplanes,self.reflectivities,
                self.H_f,self.t_f,self.r_f,self.H_los)
        else:
            self.Hf_los,self.Hf_diff = calcCIRFreqDom(self.f,
                emitters[0],collectors[0],simpleReflectingPlanes,simpleBlockingPlanes)
        
        if verbose:
            print("> Finish calculating CIR :)")

        return self.Hf_los,self.Hf_diff

    def getMtx(self,LEDs,PDs,blockingObj,reflectingObj,
        partitionDist=speed_of_light*1e-9,verbose=False):
        """ Get matrices for CIR calc.
        
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
        partitionDist: double (optional)
            Delta distance based on which we partition a plane.
            The default value is c*timesampling, timesampling is 
            assumed to 1ns by default. See Schulze's paper.
        
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

        if verbose:
            print("> Calculating matrices......(wait)")
        
        H_f,t_f,r_f,H_los = getMtxFreqDom(
            emitters[0],collectors[0],simpleReflectingPlanes,simpleBlockingPlanes)
        
        if verbose:
            print("> Finish calculating matrices :)")

        self.Nplanes = len(reflectingPlanes)
        self.H_f = H_f
        self.t_f = t_f
        self.r_f = r_f
        self.H_los = H_los
        self.reflectivities = np.array([plane[6] for plane in simpleReflectingPlanes])

        return H_f,t_f,r_f,H_los,self.reflectivities
    
    def transform(self,window="hann",L=10):
        """ Transform to time domain.
        
        Parameters
        ----------
        window: tuple or string
            see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.get_window.html
        L: int
            Number of samples in window

        Returns
        -------
        ht_los: array-like
            CIR of LOS channel.
        ht_diff: array-like
            CIR of diffuse channel 
        
        """
        from scipy.signal import get_window

        assert type(L)==int and L >=1

        if window == 'None':
            W = 1
        else:
            w = get_window(window,L)
            # Normalization
            w = w/((np.abs(np.fft.rfft(w)))[0])
            
            # the factor of 2 is important to align it with the freq. response
            W = np.fft.rfft(w,n=2*self.N)
            # f_W = np.fft.rfftfreq(2*self.N, d=self.timeSampling/2) # frequency bins for W

            # Aligning with 
            W = np.delete(W, -1) # the already-calculated frequency response
            # f_W = np.delete(f_W, -1)

            # _, gd = group_delay((np.array(w),np.array(1)))
            # meangd = gd.mean()

        # Filtering
        Hfiltered_los = self.Hf_los*W
        Hfiltered_diff = self.Hf_diff*W

        # Transform to time domain and clip the negative values
        ht_los = np.fft.irfft(Hfiltered_los,n=self.N).clip(0,np.inf)
        ht_diff = np.fft.irfft(Hfiltered_diff,n=self.N).clip(0,np.inf)

        # Further normalization

        # Junk
        # ht_los = np.fft.irfft(self.Hlos,n=self.N)
        # ht_diff = np.fft.irfft(self.Hdiff,n=self.N)

        self.ht_los = ht_los
        self.ht_diff = ht_diff

        return self.ht_los,self.ht_diff

    def plot(self,domain='frequency'):
        """ Plot the CIR.

        Parameters
        ----------
        domain: {'frequency','time'}, optional
            Domain to plot (the default value is 'frequency').

        Raises
        ------
        Exception:
            If the the method transform hasn't been invoked yet.

        """

        fig, ax = plt.subplots()
        
        if domain == 'time':
            if self.ht_los == [] or self.ht_diff == []:
                raise Exception(
                    """The time-domain transformation hasn't been 
                    carried out yet!""")
            else:

                ax.plot(self.t,(self.ht_los+self.ht_diff)/self.timeSampling)
                ax.set_xlabel(r"$t$ [s]"); ax.set_ylabel(r"$h(t)$ [1/s]");
        
        elif domain == 'frequency':

            ax.plot(self.f/1e6,10*np.log10(np.abs(self.Hf_los+self.Hf_diff)**2))
            ax.set_xlabel(r"$f$ [MHz]"); 
            ax.set_ylabel(r"$\vert H(f) \vert^2$ [dB]");










