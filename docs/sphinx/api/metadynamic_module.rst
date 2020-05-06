Metadynamic module API
======================

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

