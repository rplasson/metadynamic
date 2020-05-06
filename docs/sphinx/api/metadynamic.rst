######################
Metadynamic module API
######################


********************************
Metadynamic stochastic simulator
********************************

Provides
--------

    - :func:`.launch`: function for launching a simulation from a .json parameter file, a .hdf5 result
      file from a previous run, with eventual additional parameters, and storing the results in a
      new .hdf5 file.
    - :class:`.ResultReader`: class for reading and extracting data from a .hdf5 result file.
    - :class:`.System`: class for creating, running, and directly controlling a simmulation.
    - :class:`.Crn`: class describing a Chemical Reaction Network.
    - :obj:`.LOGGER`: global object for logging messages
    - :obj:`.MPI_STATUS`: global object for getting information about the MPI status


Module documentation
--------------------

.. automodule:: metadynamic
   :members:


***
CLI
***

.. autoprogram:: metarun:get_parser()
   :prog: metarun

      
********************************
Metadynamic interface submodules
********************************


.. automodule:: metadynamic.launcher
   :members:

.. automodule:: metadynamic.system
   :members:

.. automodule:: metadynamic.results
   :members:


*******************************
Metadynamic modeling submodules
*******************************

.. automodule:: metadynamic.proba
   :members:

.. automodule:: metadynamic.ruleset
   :members:


*******************************
Metadynamic chemical submodules
*******************************


.. automodule:: metadynamic.collector
   :members:

.. automodule:: metadynamic.chemical
   :members:


**************************      
Metadynamic I/O submodules
**************************
      
.. automodule:: metadynamic.inputs
   :members:

.. automodule:: metadynamic.outputs
   :members:

.. automodule:: metadynamic.logger
   :members:

.. automodule:: metadynamic.hdf5
   :members:

.. automodule:: metadynamic.network
   :members:


********************************      
Metadynamic utilities submodules
********************************
      
.. automodule:: metadynamic.mpi
   :members:

.. automodule:: metadynamic.ends
   :members:

.. automodule:: metadynamic.caster
   :members:

.. automodule:: metadynamic.inval
   :members:


