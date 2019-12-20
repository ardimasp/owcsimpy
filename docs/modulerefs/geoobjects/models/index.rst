Geometry Objects: Models
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

* Point source (:class:`~owcsimpy.geoobjects.models.pointsource_py.PointSource_py`)
* Bare detector (:class:`~owcsimpy.geoobjects.models.baredetector_py.BareDetector_py`)
* Room using cube (:class:`~owcsimpy.geoobjects.models.roomcube_py.RoomCube_py`)
* Human using cube (:class:`~owcsimpy.geoobjects.models.humancube_py.HumanCube_py`)
* Human (more realistic) (:class:`~owcsimpy.geoobjects.models.humancubes_py.HumanCubes_py`)
	 

Related Modules
---------------

Miscellaneous Notes
-------------------

The suffix `_py` is used to prepare in case we want to migrate the classes to C/C++.

Submodules
~~~~~~~~~~

.. toctree::
    :maxdepth: 1

    pointsource.rst
    baredetector.rst
    roomcube.rst
    humancube.rst
    humancubes.rst


