Visualising the Cosmic Web and Dark Matter Halos
================================================

**This project was developed as part of my Senior Honours Project at the University of Edinburgh. If you use this code please cite it as in my GitHub repository.** 
Our software tool is written in Python using the yt (http://yt-project.org) software toolkit and ytree (http://ytree.readthedocs.io) as well as data from GADGET.

The code is broken down into four classes, explained in each section of our documentation.:
The :ref:`main` class involves all operations that extract halo properties from our simulation, accepting various parameters as user input. 
The :ref:`movie_maker` class sets the scene using yt and renders it to produe images from our simulation, according to user input and the operation we want to carry out. The :ref:`dataset_basics` class allows us to view and understand the structure of the covering grids in our dataset, while the :ref:`simple_plots` class plots a few slices and projections.
Check out the :doc:`usage` section for further information.

.. note::

   This project is under active development.

Contents
--------

.. toctree::

   usage
      main
      movie maker
      dataset_basics
      simple_plot
