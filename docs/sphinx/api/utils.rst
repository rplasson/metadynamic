Utilities modules
=================

Modules
-------

Several modules are available for generic operations:

  - :mod:`.mpi` for handling of MPI runs
  - :mod:`.ends` for handling exceptions raised for signaling simulation termination (correct or
    incorrect)
  - :mod:`.caster` for automatic type casting
  - :mod:`.inval` for defining typed invalid values


mpi module
----------

Provides
~~~~~~~~

 - :class:`.Cont`: switchable "continue" flag
 - :class:`.MpiGate`: MPI barrier that can be opened/closed by any thread, for on-demand requests of
   synchronization
 - :class:`.MpiStatus`: Interface to general MPI operations.
 - :obj:`.MPI_STATUS`: global MpiStatus object
 - :obj:`.MPI_GATE`: global MpiGate object


Module documentation
~~~~~~~~~~~~~~~~~~~~

.. automodule:: metadynamic.mpi
   :members:

      
ends module
-----------

Provides
~~~~~~~~

 - :exc:`.Finished` is the parent ending exception
 - :exc:`.HappyEnding` is intended for normal end of the computation
 - :exc:`.Aborted` is intended to signal that the computation ended earlier than expected,
   but lead to a nonetheless correct result
 - :exc:`.BadEnding` is intended to signal that something when wrong during the computation
 - :exc:`.InputError` is intended to signal problems with file read or write, unrelated to
   the computation itself.

Inheritance diagram
~~~~~~~~~~~~~~~~~~~~~~~

.. inheritance-diagram:: metadynamic.ends
   :parts: 1

Module documentation
~~~~~~~~~~~~~~~~~~~~

.. automodule:: metadynamic.ends
   :members:


caster module
-------------

Provides
~~~~~~~~

 - :class:`.Caster`


Module documentation
~~~~~~~~~~~~~~~~~~~~
      
.. automodule:: metadynamic.caster
   :members:


inval module
------------

Provides
~~~~~~~~

 - :exc:`.InvalidError`: error to be raised when an invalid value is found where unexpected
 - :class:`.Invalid`: generic Invalid class
 - :func:`.isvalid`: function testing the validity of an object
 - :data:`.invalidint`: invalid int object.  isvalid(isvalidint) returns False;
   isinstance(invalidint, int) returns True
 - :data:`.invalidfloat`: invalid float object.  isvalid(isvalidfloat) returns False;
   isinstance(invalidfloat, float) returns True
 - :data:`.invalidstr`: invalid str object.  isvalid(isvalidstr) returns False;
   isinstance(invalidstr, str) returns True

Inheritance diagram
~~~~~~~~~~~~~~~~~~~~~~~

.. inheritance-diagram:: metadynamic.inval
   :parts: 1			    

Module documentation
~~~~~~~~~~~~~~~~~~~~
      
.. automodule:: metadynamic.inval
   :members:


