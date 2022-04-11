import yt
#yt.enable_parallelism()
# parallel = False
import ytree
import time
import glob
import numpy as np
import math
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.animation as animation
import matplotlib.image as mpimg

from multiprocessing import Pool

from yt import YTArray
from yt.visualization.volume_rendering.api import Scene, Camera
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.render_source import VolumeSource
import datetime
import unyt
from unyt import Mpc, kpc
from matplotlib import rc_context
import warnings
import yt
# from yt.utilities.parallel_tools.parallel_analysis_interface \
#     import communication_system

today = datetime.datetime.now()

# The number 4, below, is the number of processes to parallelize over, which
# is generally equal to the number of MPI tasks the job is launched with.
# If num_procs is set to zero or a negative number, the for loop below
# will be run such that each iteration of the loop is done by a single MPI
# task. Put another way, setting it to zero means that no matter how many
# MPI tasks the job is run with, num_procs will default to the number of
# MPI tasks automatically.
#num_procs = 2

# fns is a list of all the simulation data files in the directory.
fns = glob.glob("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/*.h5")
fns.sort()
count =0
# ts = yt.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/*.h5")      # load all files

#fig = plt.figure()

# This dict will store information collected in the loop, below.
# Inside the loop each task will have a local copy of the dict, but
# the dict will be combined once the loop finishes.
#my_storage = {}

# In this example, because the storage option is used in the
# parallel_objects function, the loop yields a tuple, which gets used
# as (sto, fn) inside the loop.
# In the loop, sto is essentially my_storage, but a local copy of it.
# If data does not need to be combined after the loop is done, the line
# would look like:
#       for fn in parallel_objects(fns, num_procs):
for fn in fns:#yt.parallel_objects(fns, num_procs):#, storage = my_storage):

    # Open a data file, remembering that fn is different on each task.
    ds = yt.load(fn)
    dd = ds.all_data()

    # This copies fn and the min/max of density to the local copy of
    # my_storage
    # sto.result_id = fn
    # sto.result = dd.quantities.extrema("nbody_mass")

    # Makes and saves a plot of the gas density.
    p = yt.ProjectionPlot(ds, "z",("grid", "nbody_mass"))#, cmap='arbre')
    
    p.set_zlim(("grid", "nbody_mass"), 10e67, 10e70)
    # nicen up the plot by setting the background color to the minimum of the colorbar
    p.set_background_color(("grid", "nbody_mass"))
    # hide the colorbar:
    p.hide_colorbar()
    p.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True)
    p.annotate_scale(corner="upper_right")
    
    p.save(f"{count}_proj.png")
    count+=1

# At this point, as the loop exits, the local copies of my_storage are
# combined such that all tasks now have an identical and full version of
# my_storage. Until this point, each task is unaware of what the other
# tasks have produced.
# Below, the values in my_storage are printed by only one task. The other
# tasks do nothing.
# if yt.is_root():
#     for fn, vals in sorted(my_storage.items()):
#         print (fn, vals)
        
        


# for ds in fns.piter():
#     # Load the data and create a single plot
#     ds = yt.load(fn)  # load data
#     p = yt.ProjectionPlot(ds, "z", ("grid", "nbody_mass"), width=(80000, "kpc"))
#     p.set_cmap(field=("grid", "nbody_mass"), cmap="kamae")
#     # Ensure the colorbar limits match for all plots
#     p.set_zlim(("grid", "nbody_mass"), 10e67, 10e70)
#     # nicen up the plot by setting the background color to the minimum of the colorbar
#     p.set_background_color(("grid", "nbody_mass"))
#     # hide the colorbar:
#     p.hide_colorbar()

#     # hide the axes, while still keeping the background color correct:
#     p.hide_axes(draw_frame=True)
#     p.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True)
#     p.annotate_scale(corner="upper_right")

#     # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
#     plot = p.plots[("grid", "nbody_mass")]
#     plot.figure = fig
#     plot.axes = grid[i].axes
#     plot.cax = grid.cbar_axes[i]

#     # Finally, this actually redraws the plot.
#     p._setup_plots()
# plt.savefig("{:%Y_%m_%d_%H_%M_%S}".format(today)+"_multi.png")