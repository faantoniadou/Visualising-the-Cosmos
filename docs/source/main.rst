main
===========

This file contains the ``ExploreHalo`` class which carries our halo merger tree operations.

Upon initialisation the user is prompted to enter information about the halo we want to analyse. These are the instance variables ``self.quantity`` and ``self.how_large``. The former is the quantity to search for in halo catalogue and the latter is the size of the quantity possessed by our halo, relative to the maximum. For instance, if we set ``self.quantity = 'mass' `` and ``self.how_large = 3``, we would be looking at the dark matter halo with the third to largest mass. 

A list of the methods of this class is provided below:

* ``ExploreHalo.load_halo()`` :

  Subsequently, the ``ExploreHalo.load_halo()`` function loads the ``Arbor`` object of the halo of interest and obtains the ``TreeNode`` (halo) of interest (target halo).
  
.. autofunction:: ExploreHalo.load_halo
  
* ``ExploreHalo.get_properties()`` :

  Gets target halo and its progenitors' properties. This includes halo position, halo radius, halo mass, progenitors' snapshot numbers and progenitor     positions. 
  
.. autofunction:: ExploreHalo.get_properties

* ``ExploreHalo.get_neighbour()`` :

  Gets neighbouring halo based on criteria. The user is prompted to input how far from the target halo they would like to search for to obtain another halo's position. This is done in units of the target halo's diameter (e.g. inputting the number 6 would mean that we would like to search for neighbouring halos at a distance equal to 6 times the diameter of the target halo). Then the neighbouring halo to select is picked based on its mass. However, this method can be modified to use a sphere object (see https://ytree.readthedocs.io/en/latest/Arbor.html?) instead. This is used in the flyby method described below.
  
.. autofunction:: ExploreHalo.get_neighbour
  
* ``ExploreHalo.neighbour_path()``:

  Creates path between the target halo and neighbouring halo. This is used later for the camera to navigate through. The user may change the number of frames to be used in the movie. This defines the number of points to be used in the array of camera positions constructed in this method.
  
.. autofunction:: ExploreHalo.neighbour_path

* ``ExploreHalo.make_prog_path()``:

  Makes a smooth camera path for a time series movie of the target halo. This uses the savgol filter (https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html) to produce a smooth camera path, in case any sharp jumps in camera focus position occur. A possible cause of this could be major mergers.

.. autofunction:: ExploreHalo.make_prog_path

* ``ExploreHalo.make_halo_region()``:

  Create data object to define the spherical region around the target halo (see more at https://yt-project.org/doc/developing/creating_datatypes.html?highlight=data%20object). This can be used in `ExploreHalo.get_neighbour()`` as well.

.. autofunction:: ExploreHalo.get_flybys_energy

* ``ExploreHalo.get_flybys_energy()``:

  Loops through halos in the region created above and identifies the type of interaction of the target halo with these. In particular, it identifies flybys based on certain criteria as described in my report. The main criterion is that the total energy of the two-body system is negative, such that the two halos are not energetically bound and thus will not merge in the future. All parameters of interest are then stored in a Pandas ``DataFrame`` saved in a .csv file with the node's number. This takes quite a while to run, so it could benefit from parallelism strategies.
  
.. autofunction:: ExploreHalo.get_flybys_energy

* ``ExploreHalo.flyby_posns()``:

  Plots graphs to visualise the position of the flyby halo identified (requires user input of the .csv file name). It plots both a 2D plot and a 3D diagram of the flyby's position.
 
.. autofunction:: ExploreHalo.flyby_posns

* ``ExploreHalo.flyby_vis``

  Plots projection plots to visualise the flyby's path. The position of the flyby is annotated.
  
.. autofunction:: ExploreHalo.flyby_vis

* ``ExploreHalo.plot_camera_path``

  Used to plot the path of the camera in progenitor time series movies.
  
.. autofunction:: ExploreHalo.plot_camera_path









