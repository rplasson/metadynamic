Modeling modules
================


Modules
-------

Several modules are available for performing the simulations:

  - :mod:`.proba` implements the Gillespie's algorithm.
  - :mod:`.ruleset` provides the interface for defining sets of rules for creating chemical systems.


    
proba module
------------
      
Provides
~~~~~~~~

 - :class:`.Probalist`

Module documentation
~~~~~~~~~~~~~~~~~~~~

.. automodule:: metadynamic.proba
   :members:


      
ruleset module
--------------
      
Provides
~~~~~~~~

 - classes for the description of a rule model:
    - :class:`.Parameters`: Maintain a set of parameters {name:values}
      that can be linked with each other.
    - :class:`.Descriptor`: Tool for describing a compound from its name
      (categories and properties)
    - :class:`.Rule`: Describes a given rule for building a reaction
    - :class:`.Model`: Full description of a rule model
 - generators to be used for building models:
    - :obj:`.Paramrel` generators: :func:`.parmul` , :func:`.arrhenius`, :func:`.linproc`.
    - :obj:`.ConstBuilder` generators: :func:`.kinvar` , :func:`.kalternate`, :func:`.kdualchoice`.
    - :obj:`.VariantBuilder` generators: :func:`.novariant_gen`, :func:`.singlevariant`,
      :func:`.rangevariant`.
    - :obj:`.ProdBuilder` generators: :func:`.joiner`, :func:`.splitter`.

Module documentation
~~~~~~~~~~~~~~~~~~~~

.. automodule:: metadynamic.ruleset
   :members:

