import numpy as np
from numpy.random import rand

# define a Uniform Distribution
U = lambda MIN, MAX, SAMPLES: rand(*SAMPLES.shape) * (MAX - MIN) + MIN

class RandOrientationArdimas(object):
    """ A base class to generate a random orientation based on [1].

    Parameters
    ----------
    position: {'sitting','walking','standing'} 
        Activity of the user
    seed: int
        Seed number of random generator (1412 default)


    Attributes
    ----------
    position
    paramTheta: tuple(3,)
        Parameters of theta0 based on Table I in [1].
        paramTheta = (A,f,sigma)
    paramOmega: tuple(3,)
        Parameters of omega0 based on Table I in [1].
        paramOmega = (A,f,sigma)

    Notes
    -----
    Refs: 
    [1] A. A. Purwita, M. D. Soltani, M. Safari and H. Haas, "Terminal Orientation in OFDM-Based LiFi Systems," in IEEE Transactions on Wireless Communications, vol. 18, no. 8, pp. 4003-4016, Aug. 2019.

    """

    def __init__(self,position,seed=1412):

        # Assertion check
        # position
        assert position.lower() in {'sitting','walking','standing'}

        # assign default values of paramTheta and paramOmega based on position
        if position.lower() == 'sitting':
            self.paramTheta = (1.88,0.67,5.91)
            self.paramOmega = (1.31,1.46,3.22)
        elif position.lower() == 'standing':
            self.paramTheta = (3.22,1.86,7.59)
            self.paramOmega = (1.31,1.46,3.22)
        elif position.lower() == 'walking':
            self.paramTheta = (3.22,1.86,7.59)
            self.paramOmega = (3.15,1.71,9.48)

        self.seed = seed;

        # np.random.seed(self.seed)

    def genSamples(self,time):
        """ Generate samples based on time instance.
        
        Parameters:
        time: array
            Time instance.

        """

        time = np.array(time)

        A,f,sigma = self.paramTheta
        theta0 = A*np.sin(2*np.pi*f*time+U(-np.pi,np.pi,np.array([1]))) + sigma*np.random.randn()

        A,f,sigma = self.paramOmega
        omega0 = A*np.sin(2*np.pi*f*time+U(-np.pi,np.pi,np.array([1]))) + sigma*np.random.randn()

        return np.deg2rad(theta0),np.deg2rad(omega0)





