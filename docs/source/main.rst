main
=====


User input
----------
This file contains the ``ExploreHalo`` class which carries our halo merger tree operations.

Upon initialisation the user is prompted to enter information about the halo we want to analyse. These are the instance variables ``self.quantity`` and ``self.how_large``. The former is the quantity to search for in halo catalogue and the latter is the size of the quantity possessed by our halo, relative to the maximum. For instance, if we set ``self.quantity = 'mass' `` and ``self.how_large = 3``, we would be looking at the dark matter halo with the third to largest mass. 

A list of the methods of this class is provided below:

* ``ExploreHalo.load_halo()``
  Subsequently, the ``ExploreHalo.load_halo()`` function loads the ``Arbor`` object of the halo of interest.
  .. autofunction:: ExploreHalo.load_halo



To retrieve a list of random ingredients,
you can use the ``lumache.get_random_ingredients()`` function:

.. autofunction:: lumache.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

