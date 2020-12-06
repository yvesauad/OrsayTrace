.. -*- mode: rst -*-

|python_version|_ |pypi_version|_ |rtd|_

.. |pypi_version| image:: http://img.shields.io/pypi/v/orsaytrace.svg?style=flat
.. _pypi_version: https://pypi.python.org/pypi/orsaytrace

.. |python_version| image:: https://img.shields.io/pypi/pyversions/orsaytrace.svg?style=flat
.. _python_version: https://pypi.python.org/pypi/orsaytrace

.. |rtd| image:: https://readthedocs.org/projects/orsaytrace/badge/?version=latest
.. _rtd: https://readthedocs.org/projects/orsaytrace/?badge=latest


OrsayTrace
----------

Orsaytrace is an open source python library for flexible optical simulations in :math:`3D` dimensions.

Orsaytrace makes it easy to add elements or transformations sequentially. Multiple shapes can be created using
building blocks such as spheres, prisms, parabolas and spatial transformations like rotations or translations.

Data analysis are done by the means of plans, which save photons that cross them during simulation run. Plans
are powerful objects and can be instantiated as conditional plans, which guarantees that appended photons
possesses a given attribute in a certain range. A plethora of pre build functions are already implemented
for photons subset, but users can build their own
using all available photon objects.

Contributing 
************

Contributing is encouraged for everyone. Please be in touch if you wish to help.
