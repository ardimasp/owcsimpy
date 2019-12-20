import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.constants import speed_of_light

from warnings import warn

from owcsimpy.misc import flatten

from owcsimpy.cir.freqdomaincir import FreqDomainCIR


class SphereModelCIR(object):
    """ CIR calculation by means of the sphere model [1].
    
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
    [1] V. Jungnickel, V. Pohl, S. Nönnig, C. von Helmolt, "A physical model of the wireless infrared communication channel", IEEE J. Sel. Areas Commun., vol. 20, no. 3, pp. 631-640, Apr. 2002.
    [2] H. Schulze, "Frequency-Domain Simulation of the Indoor Wireless Optical Communication Channel," in IEEE Transactions on Communications, vol. 64, no. 6, pp. 2551-2562, June 2016.
    
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

    def calc(self,LEDs,PDs,blockingObj,reflectingObj):
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

        Notes
        -----
        This is a quick implementation of [1] based on the approx. given by (12)-(16) in [2]

        """

        # Quickly get the propertiest of the elements
        L,W,H = reflectingObj.L, reflectingObj.W, reflectingObj.H
        fdcir = FreqDomainCIR(N=self.N,freqSampling=self.freqSampling)
        *tmp,H_los,reflectivities = fdcir.getMtx(LEDs,PDs,blockingObj,reflectingObj)

        Aroom = 2*L*W+2*L*H+2*H*W
        Vroom = L*W*H

        Arx = PDs.area

        # Calculate the average of the reflectivities
        rho = np.average(reflectivities)

        # Calculate (12) in [2]
        eta = Arx/Aroom*rho/(1-rho)

        # Calculate (15) in [2]
        tau = -1/np.log(rho)*4*27/(Aroom*speed_of_light)

        # Calculate (16) in [2]
        Hf_diff = eta/(1+1j*2*np.pi*fdcir.f*tau)

        Hf_los = np.array([
            np.abs(H_los)*np.cos(-1*2*np.pi*np.angle(H_los)*fi)
            +1j*(np.abs(H_los)*np.sin(-1*2*np.pi*np.angle(H_los)*fi)) 
            for fi in self.f])

        return Hf_los, Hf_diff






