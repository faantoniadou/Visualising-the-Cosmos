import yt
#from yt.mods import *           # deprecated
import ytree

import numpy as np
import math
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.animation as animation
import matplotlib.image as mpimg

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

#yt.enable_parallelism()

nframes = int(input('Number of frames: '))
factor = float(input('Zoom factor: '))
# quantity maximum halo to find:
quantity = str(input('Maximum quantity to search for in halo catalogue: '))
how_large = int(input(f'Order to largest {quantity} (e.g. 1 is maximum and 2 is second to largest): '))
#comm = communication_system.communicators[-1]
print('where am I?')

# this makes the file name for the covering grids directory
def make_grid_name(fnumber):
    print('where am I?')
    grid_name = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{fnumber:03}_covering_grid.h5"
    return grid_name

# this makes file name for halo trees 
# def get_halos():
#     print('where am I?')
#     halo_dir = "/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5"
#     return halo_dir

def get_halo_centre():
    print('where am I?')
    # Finds halo with the maximum of the quantity (field) the user provides.
    a = ytree.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5")
    
    values = (a[quantity]).tolist()
    
    values.sort()
    i = np.where(a[quantity] == (values[-how_large]))[0][0]           #['mass'])
    my_tree = a[i]

    halo_position = my_tree["position"]
    radius = my_tree["virial_radius"].to("unitary")
    width = 3*radius
    print(i)
    return halo_position, width

def zoom_halo():
    print('Zooming now')
    
    # Load the dataset
    ds = yt.load(make_grid_name(101))

    # Create a volume rendering
    sc = yt.create_scene(ds, field=("grid", "nbody_mass"))
    sc.grey_opacity = True

    source = sc[0]
    cam = sc.camera
    # Now increase the resolution
    sc.camera.resolution = (1024, 1024)
    
    # for resolution purposes
    source.set_use_ghost_zones(True)
    
    # we focus on the halo we want from the beginning
    final_position, final_width = get_halo_centre()

    # set initial camera setting
    #cam.width = final_width * 10
    #cam.focus = final_position
    
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
    frame = 0
    sc.save("camera_movement_%04i.png" % frame, sigma_clip=4.)
    frame += 1
    
    for _ in cam.iter_move(final_position, nframes):
        print('saved')
        sc.render()
        sc.save("camera_movement_%04i.png" % frame)
        frame += 1
        
        cam.width[:2] = cam.width[:2] / factor



def rotate_halo():
    print("Let's dance")
    
    # we focus on the halo we want from the beginning
    final_position, final_width = get_halo_centre()
    print(final_position, final_width)
    # Load the dataset
    ds = yt.load(make_grid_name(101))

    # Create a volume rendering
    sc = yt.create_scene(ds, field=("grid", "nbody_mass"))
    sc.grey_opacity = True

    source = sc[0]
    cam = sc.camera
    
    # for resolution purposes
    source.set_use_ghost_zones(True)

    # set initial camera setting
    cam.width = final_width * 2
    cam.focus = final_position
    
    # Set up a custom transfer function using the TransferFunctionHelper.
    # We use 10 Gaussians evenly spaced logarithmically between the min and max
    # field values.
    tfh = TransferFunctionHelper(ds)
    tfh.set_field(("grid", "nbody_mass"))
    tfh.set_log(True)
    tfh.set_bounds()
    tfh.build_transfer_function()
    tfh.tf.add_layers(10, colormap="cmyt.arbre")

    # Grab the first render source and set it to use the new transfer function
    render_source = sc.get_source()
    render_source.transfer_function = tfh.tf
    
    '''
    rotating process begins here
    '''
    # save an image at the starting position
    frame = 0
    sc.save(f"rotation_{frame}_{quantity}.png", sigma_clip=4.)
    frame += 1
    
    for i in cam.iter_rotate(2*np.pi, 100):
        #cam.width = final_width * 4
        # im = sc.render()
        sc.camera.resolution = (1024, 1024)
        sc.save(f"rotation_{i}_{quantity}_order={how_large}.png", sigma_clip=4.)
        
rotate_halo()
    