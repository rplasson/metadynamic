Metadynamic I/O modules
==========================

Modules
-------

Several modules are available for dealing with exchanges between the simulations and input or output
files:

  - :mod:`.inputs` for reading parameters from json files.
  - :mod:`.outputs` for the organization of output files and folders
  - :mod:`.logger` for logging facilities to standard output, text files, or hdf5 files.
  - :mod:`.hdf5` for writing simulation results and parameters in a structured hdf5 file.
  - :mod:`.network` for writing and reading chemical reaction networks in graphviz dot format.


inputs module
-------------

Provides
~~~~~~~~

 - :exc:`.LockedError`: Exception raised when attempting to modify a locked Readerclass
 - :class:`.Castreader`: Extend a :class:`.caster.Caster` for dealing with :class:`.Readerclass`, converting them as dictionary
 - :class:`.Readerclass`: dataclass with interface for reading its data from json files
 - :class:`.RuleParam`: Parameters for a reaction rule
 - :class:`.RulesetParam`: Parameters for a reaction ruleset
 - :class:`.StatParam`: Parameters for statistics
 - :class:`.MapParam`: Parameters for map statistics
 - :class:`.Param`: Run parameters
 - :class:`.DotParam`: Parameters for graphviz CRN representation

Module documentation
~~~~~~~~~~~~~~~~~~~~

.. automodule:: metadynamic.inputs
   :members:


outputs module
--------------

Provides
~~~~~~~~

 - :class:`.Output`


Module documentation
~~~~~~~~~~~~~~~~~~~~

.. automodule:: metadynamic.outputs
   :members:


logger module
-------------

Provides
~~~~~~~~

 - :class:`.Timer`: Simple class for tracking passed processing time from a checkpoint
 - :class:`.Log`: High level class for MPI-aware logging, with log save in hdf5.
 - :obj:`.LOGGER`: Global Log object

Module documentation
~~~~~~~~~~~~~~~~~~~~


.. automodule:: metadynamic.logger
   :members:


hdf5 module
-----------

Provides
~~~~~~~~

 - :class:`.ResultWriter`: storage for all the simulation data and results

Module documentation
~~~~~~~~~~~~~~~~~~~~

.. automodule:: metadynamic.hdf5
   :members:


network module
--------------

Provides
~~~~~~~~

 - :class:`.Data2dot`: generate a graphviz Digraph from dict description of compounds/reactions
 - :class:`Scaler`: for scaling data to graph dimensions.

Module documentation
~~~~~~~~~~~~~~~~~~~~


.. automodule:: metadynamic.network
   :members:



