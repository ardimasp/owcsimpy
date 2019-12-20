Channel Impulse Response (CIR)
==============================

Introduction
------------

It is obvious that the calculation of CIR is essential, yet it is easily overlooked due to its complexity. First, I will focus on the implementation of the time-domain [Carruthers2002]_ and frequency-domain approaches [Schulze2016]_. 

The main advantage of the frequency-domain approach is that it considers infinite number of reflections. However, when it comes to time-domain CIR, freq. domain approach needs windowing filter to transform it to time domain. To give you a concrete example of what the disadvantage of the freq. domain approach is that a flat frequency response will not give us an impulse in time-domain if the time bin is not correct. In fact, it gives us a sinc func (can be a negative value). Therefore, a careful filtering is necessary. On the other hand, the time-domain approach is more resilient in terms of transforming to other domain.

.. tip ::
	
	In case you don't know, there's a rule of thumb on how small we can sufficiently partition a plane, i.e., :math:`\Delta x \leq c \Delta t`, where :math:`c` is the speed of light and :math:`\Delta t` is time sampling. More detailed discussions can be found in [Schulze2016], [Barry1993]_ and [Schulze2018]_ (especially the latter one).


CIR for IM/DD
-------------

* Time domain (:class:`~owcsimpy.cir.timedomaincir.TimeDomainCIR`)
* Freq. domain (:class:`~owcsimpy.cir.freqdomaincir.FreqDomainCIR`)

Check following notebook for examples.

.. raw:: html

    <iframe src="../../../../../notebooks/UseCase_1_CIR.html" height="345px" width="100%"></iframe>


References
----------

.. [Carruthers2002] `J. B. Carruthers and P. Kannan, "Iterative site-based modeling for wireless infrared channels," in IEEE Transactions on Antennas and Propagation, vol. 50, no. 5, pp. 759-765, May 2002. <https://ieeexplore.ieee.org/abstract/document/1011244>`_
.. [Schulze2016] `H. Schulze, "Frequency-Domain Simulation of the Indoor Wireless Optical Communication Channel," in IEEE Transactions on Communications, vol. 64, no. 6, pp. 2551-2562, June 2016. <https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7457357>`_
.. [Barry1993] `J. R. Barry, J. M. Kahn, W. J. Krause, E. A. Lee and D. G. Messerschmitt, "Simulation of multipath impulse response for indoor wireless optical channels," in IEEE Journal on Selected Areas in Communications, vol. 11, no. 3, pp. 367-379, April 1993. <https://ieeexplore.ieee.org/document/219552>`_
.. [Schulze2018] `Schulze, Henrik: 'FEM simulations for the wireless optical indoor communication channel', IET Optoelectronics, 2018, 12, (2), p. 94-105. <https://digital-library.theiet.org/content/journals/10.1049/iet-opt.2017.0089>`_


.. Submodules
.. ~~~~~~~~~~

.. toctree::
    :hidden:

    timedomaincir.rst
    freqdomaincir.rst
