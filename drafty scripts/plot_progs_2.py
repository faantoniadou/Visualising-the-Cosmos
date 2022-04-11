import yt
yt.enable_parallelism()
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
#from tqdm import tqdm
import warnings
import yt
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import communication_system
    
    

today = datetime.datetime.now()

# volume rendering only uses a power of two number of processors
num_procs = 16

# quantity maximum halo to find:
quantity = str(input('Maximum quantity to search for in halo catalogue: '))

def get_halo_centre():
    print('where am I?')
    # Finds halo with the maximum of the quantity (field) the user provides.
    a = ytree.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5")

    
    # find most massive halo
    i = np.argmax(a[quantity])             #['mass'])

    my_tree = a[i]

    halo_position = my_tree["position"]
    radius = my_tree["virial_radius"].to("unitary")
    
    # get progenitors 
    prog_positions = my_tree["prog", "position"]     #give position of most massive progenitor
    
    prog_redshifts = my_tree["prog", "redshift"]
    prog_snaps = my_tree["prog", "Snap_idx"]

    width = 3*radius

    return halo_position, width, prog_positions, prog_snaps

t1 = time.time()



# empty list of files to append 
files = []

# we focus on the halo we want from the beginning
final_position, width, prog_positions, prog_snaps = get_halo_centre()
print(len(prog_snaps))

# load files in which halo is found
for fnumber, in prog_positions:
    files.append(f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{fnumber:03}_covering_grid.h5")

# fns is a list of all the simulation data files in the directory.
fns = glob.glob(files)
fns.sort()

for fn in yt.parallel_objects(fns, num_procs):
    print(fn)
    # Open a data file, remembering that fn is different on each task.
    ds = yt.load(fn)
    dd = ds.all_data()
    
    # Create a volume rendering
    sc = yt.create_scene(ds, field=("grid", "nbody_mass"))
    sc.grey_opacity = True

    source = sc[0]
    cam = sc.camera
    
    # for resolution purposes
    source.set_use_ghost_zones(True)

    # set initial camera setting
    cam.width = width * 3
    cam.focus = prog_positions[i]
    
    # Now increase the resolution
    sc.camera.resolution = (1024, 1024)
    
    # Set up a custom transfer function using the TransferFunctionHelper.
    # We use 10 Gaussians evenly spaced logarithmically between the min and max
    # field values.
    tfh = TransferFunctionHelper(ds)
    tfh.set_field("nbody_mass")
    tfh.set_log(True)
    tfh.set_bounds()
    tfh.build_transfer_function()
    tfh.tf.add_layers(10, colormap="cmyt.arbre")

    # Grab the first render source and set it to use the new transfer function
    render_source = sc.get_source()
    render_source.transfer_function = tfh.tf
    
    '''
    Moving process begins here
    '''
    
    # save an image at the starting position
    
    sc.save(f"progenitor_{fn}_{quantity}.png", sigma_clip=4.)    
    
    
    
    
    
    
    
    
    
    
    

#num_procs = 4       # number of processors to use

# for i in range(len(prog_snaps)-8):#yt.parallel_objects(prog_snaps, num_procs):        # loop through snapshots and plot the halo's progenitors
#     # Load the dataset
#     print(i)
#     ds = yt.load(make_grid_name(prog_snaps[i]))
    
#     # Create a volume rendering
#     sc = yt.create_scene(ds, field=("grid", "nbody_mass"))
#     sc.grey_opacity = True

#     source = sc[0]
#     cam = sc.camera
    
#     # for resolution purposes
#     source.set_use_ghost_zones(True)

#     # set initial camera setting
#     cam.width = width * 3
#     cam.focus = prog_positions[i]
    
#     # Now increase the resolution
#     sc.camera.resolution = (1024, 1024)
    
#     # Set up a custom transfer function using the TransferFunctionHelper.
#     # We use 10 Gaussians evenly spaced logarithmically between the min and max
#     # field values.
#     tfh = TransferFunctionHelper(ds)
#     tfh.set_field("nbody_mass")
#     tfh.set_log(True)
#     tfh.set_bounds()
#     tfh.build_transfer_function()
#     tfh.tf.add_layers(10, colormap="cmyt.arbre")

#     # Grab the first render source and set it to use the new transfer function
#     render_source = sc.get_source()
#     render_source.transfer_function = tfh.tf
    
#     '''
#     Moving process begins here
#     '''
    
#     # save an image at the starting position
    
#     sc.save(f"progenitor_{i}_{quantity}.png", sigma_clip=4.)    

# t2 = time.time()

# # if yt.is_root():
# #     print ("BigStuff took %.5e sec" % (t2 - t1))

# #animation()