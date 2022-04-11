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
from yt.visualization.volume_rendering.api import Scene, Camera, create_volume_source
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.render_source import VolumeSource
import datetime
import unyt
from unyt import Mpc, kpc
from matplotlib import rc_context
import warnings
import yt
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import communication_system
    
# to create smoother camera path:
from scipy.interpolate import CubicSpline           # method 1
from scipy.signal import savgol_filter              # method 2
from matplotlib.pyplot import figure

    
    

today = datetime.datetime.now()


# quantity maximum halo to find:
quantity = 'mass' #str(input('Maximum quantity to search for in halo catalogue: '))

# if we want to start from a specific redshift form the snap catalogue onwards
starting_snap = 0 #int(input("Starting snap: "))



# Load the dataset
# this makes the file name for the covering grids directory
def make_grid_name(fnumber):
    grid_name = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{fnumber:03}_covering_grid.h5"
    return grid_name



def get_halo_centre():
    print('where am I?')
    # Finds halo with the maximum of the quantity (field) the user provides.
    a = ytree.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5")
    
    # find halo with the maximum of the quantity we want
    i = np.argmax(a[quantity])             #['mass'])
    my_tree = a[i]

    halo_position = my_tree["position"]
    radius = my_tree["virial_radius"].to("unitary")
    
    # get progenitors 
    all_positions = my_tree["prog", "position"]     #give position of most massive progenitor
    
    prog_redshifts = my_tree["prog", "redshift"]
    #print(prog_redshifts)
    prog_snaps = my_tree["prog", "Snap_idx"]

    width = 3*radius

    return halo_position, width, all_positions, prog_snaps

get_halo_centre()

def make_path(data_points):
    '''
    Makes a smooth camera path
    '''

    # we focus on the halo we want from the beginning
    final_position, width, prog_positions, prog_snaps = get_halo_centre()

    camera_path = np.empty(prog_positions.shape)

    window_size = 2*int((len(prog_positions)/4)) + 1            # window size must always be odd
    pol_degree = 6

    smoother_pos_x = savgol_filter(prog_positions[:,0], window_size, pol_degree)
    smoother_pos_y = savgol_filter(prog_positions[:,1], window_size, pol_degree)
    smoother_pos_z = savgol_filter(prog_positions[:,2], window_size, pol_degree)

    # set new camera paths
    camera_path[:,0], camera_path[:,1], camera_path[:,2] = smoother_pos_x, smoother_pos_y, smoother_pos_z
    #camera_path = np.concatenate(camera_path, axis=1)
    
    # visualise how data points were smoothed
    x = np.arange(len(prog_positions))
    
    plt.plot(x, prog_positions[:,0], label='x_actual',linewidth=1)
    plt.plot(x, smoother_pos_x, label='x_smooth',linewidth=1)

    plt.plot(x, prog_positions[:,1], label='y_actual',linewidth=1)
    plt.plot(x, smoother_pos_y, label='y_smooth',linewidth=1)

    plt.plot(x, prog_positions[:,2], label='z_actual',linewidth=1)
    plt.plot(x, smoother_pos_z, label='z_smooth',linewidth=1)

    plt.legend()

    plt.savefig(f'smoother_points_{quantity}_.png')

    return camera_path



#num_procs = 4       # number of processors to use
def time_haloes():
    final_position, width, all_positions, prog_snaps = get_halo_centre()
    prog_positions = make_path(all_positions)
    
    for i in range(len(prog_snaps)-starting_snap):#yt.parallel_objects(prog_snaps, num_procs):        # loop through snapshots and plot the halo's progenitors
        # Load the dataset
        
        ds = yt.load(make_grid_name(prog_snaps[i+starting_snap]))
        dd = ds.sphere(prog_positions[i+starting_snap], width)
        
        # Create a volume rendering
        #sc = yt.create_scene(ds, field=("grid", "nbody_mass"))
        sc = Scene()
        vol = create_volume_source(ds, field=("grid", "nbody_mass") )
        
        #source = sc[0]
        #cam = sc.camera
        cam = sc.add_camera(dd)#, lens_type='perspective')
        # for resolution purposes
        #source.set_use_ghost_zones(True)

        # set initial camera setting
        cam.width = width * 3
        cam.focus = prog_positions[i+starting_snap]
        
        # increase the resolution
        sc.camera.resolution = (1024, 1024)
        
        # Set up a custom transfer function using the TransferFunctionHelper.
        # We use 10 Gaussians evenly spaced logarithmically between the min and max
        # field values.
        
        ad = ds.all_data()
        mi, ma = ad.quantities.extrema('nbody_mass')
        print(mi,ma)
        #we bump up our minimum to cut out some of the background
        bounds = (2.0, np.log10(ma))
        
        tfh = TransferFunctionHelper(ds)
        tfh.set_field(("grid", "nbody_mass"))
        tfh.set_log(True)
        tfh.set_bounds(bounds)
        tfh.build_transfer_function()
        tfh.tf.add_layers(10, colormap="cmyt.arbre")

            
        render_source = sc.add_source(vol)
        render_source.transfer_function = tfh.tf #transfer_function = tfh.tf
        
        #tf = render_source.transfer_function
        #tf.clear()
        #tf.add_layers(10, colormap='cmyt.arbre')
        
        # save an image at the starting position
        
        sc.save(f"progenitor_{i+starting_snap}_{quantity}.png", sigma_clip=4.)    
       
 
time_haloes()