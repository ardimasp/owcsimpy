Geometry Objects: Basics
========================

Introduction
------------

This basic geometry modules provide base classes with necessary methods for 
higher-level models. We envision that the higher-level models will be something
like Minecraft objects which are formed by stacking cubes. In addition, we will 
focus on point source models to represent optical devices such as LED and PD.  


Available Entities
------------------

The following entities are currently available in the **basic** geometry module:

* Vector (:class:`~owcsimpy.geoobjects.bases.vector_py.Vector_py`)
	A vector is mainly used to represent an LED that is modeled as a point source. 
* Parametric line (:class:`~owcsimpy.geoobjects.bases.paramline_py.ParamLine_py`)
	Mainly used for partitioning a rectangular plane. This class will be potentially
	used to represent mobility of users.
* Circle (:class:`~owcsimpy.geoobjects.bases.circle_py.Circle_py`)
	A bare detector will be modeled using a circle.
* Rectangular plane (:class:`~owcsimpy.geoobjects.bases.rectplane_py.RectPlane_py`)
	Smallest element of higher-level models for the calculation of CIR.
* Cube (:class:`~owcsimpy.geoobjects.bases.cube_py.Cube_py`)
	This will be a base class for higher-level models, such as, blocking objects, 
	rooms, or furnitures. The main method used from this class is `getPartition`.

Related Modules
---------------

Miscellaneous Notes
-------------------

The suffix `_py` is used to prepare in case we want to migrate the classes to C/C++.

Submodules
~~~~~~~~~~

.. toctree::
    :maxdepth: 2

    vector.rst
    paramline.rst
    circle.rst
    rectplane.rst
    cube.rst
