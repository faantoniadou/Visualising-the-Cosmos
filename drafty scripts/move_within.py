'''
This script is used to move between haloes. We will try to place the camera within the halo 
'''
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
quantity = 'mass'#str(input('Maximum quantity to search for in halo catalogue: '))
frames = int(input('Number of frames:'))
# if we want to start from a specific redshift form the snap catalogue onwards
#starting_snap = int(input("Starting snap: "))
how_large = int(input(f'Order to largest {quantity} (e.g. 1 is maximum and 2 is second to largest): '))

zoom_factor = int(input("Zoom out factor: "))
how_far = float(input("How far should I search (in units of the halo's radius)? "))


# Load the dataset
# this makes the file name for the covering grids directory
def make_grid_name(fnumber):
    grid_name = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{fnumber:03}_covering_grid.h5"
    return grid_name



def get_halo_centre():
    '''
    This finds the centre of the halo pairs we want to navigate between
    
    find haloes below a mass threshold
    if linalg norm is less than 5* radius or most massive
    make path
    navigate
    if not
    display error message
    '''
    
    # Finds halo with the maximum of the quantity (field) the user provides.
    a = ytree.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5")
    
   # find halo with the order of maximum of the quantity we want
    values = (a[quantity]).tolist()
    
    values.sort()
    i = np.where(a[quantity] == (values[-how_large]))[0][0]           #['mass'])
    #i = np.argmax(a[quantity])
    print(i)
    my_tree = a[i]
    
    initial_position = my_tree["position"].to("unitary")
    initial_radius = my_tree["virial_radius"].to("unitary")

    # get progenitors 
    #all_positions = my_tree["prog", "position"]     # give position of most massive progenitor
    
    prog_redshifts = my_tree["prog", "redshift"]
    prog_snaps = my_tree["prog", "Snap_idx"]

    width = 2*initial_radius
    
    
    #empty list for halo distances
    # distances = []
    # positions = []

    # for halo in a.select_halos("(tree['forest', 'mass'].to('Msun') > 6e14) & (tree['forest', 'mass'].to('Msun') > 690072000000000.0)"):        # exclude the most massive halo
    #     pair_pos = halo["position"].to("unitary")
    #     distance = np.linalg.norm(np.array(pair_pos) - np.array(initial_position))
    #     positions.append(np.array(pair_pos))
    #     distances.append(distance)
    
    # if len(distances) == 0:
    #     print('No haloes found.')
    # # now we want to choose a halo that is neither too close nor too far from the target
    # idx = (np.abs(distances - how_far * np.array(width))).argmin()
    # final_position = positions[idx]
    initial_position, final_position = [0.7368401, 0.19819638, 0.0234113 ], [0.0433627, 0.0624088, 0.95199627]
    
    print(np.array(initial_position), np.array(final_position))
    #return halo_position, width, all_positions, prog_snaps, haloes
    
    return np.array(initial_position), width, prog_snaps, np.array(final_position)







def make_path():#no_steps):
    '''
    Makes a smooth camera path
    '''
    no_steps = frames
    initial_position, width, prog_snaps, final_position = get_halo_centre()
    
    x_initial, y_initial, z_initial = np.array(initial_position[0]), np.array(initial_position[1]), np.array(initial_position[2])
    x_final, y_final, z_final = np.array(final_position[0]), np.array(final_position[1]), np.array(final_position[2])
    
    # step_x = (np.array(x_final) - np.array(x_initial))/no_steps
    # step_y = (np.array(y_final) - np.array(y_initial))/no_steps
    # step_z = (np.array(z_final) - np.array(z_initial))/no_steps
    
    
    # x_path, y_path, z_path = np.arange(x_initial, x_final+step_x, step_x), np.arange(y_initial, y_final+step_y, step_y), \
    #     np.arange(z_initial, z_final+step_z, step_z)
    
    x_path, y_path, z_path = np.linspace(x_initial, x_final, no_steps), np.linspace(y_initial, y_final, no_steps), \
        np.linspace(z_initial, z_final, no_steps)
        

    
    camera_path = np.vstack((x_path, y_path, z_path)).T
    #camera_path = np.empty((no_steps, 3)
    
    #camera_path[:,0], camera_path[:,1], camera_path[:,2] = x_path, y_path, z_path
    
    return initial_position, width, prog_snaps, final_position, camera_path


#num_procs = 4       # number of processors to use
def move():
    '''
    This function executes the camera move between the haloes
    '''
    
    final_position, width, all_positions, prog_snaps, prog_positions = make_path()

    ds = yt.load(make_grid_name(101))
    
    # Create a volume rendering
    sc = yt.create_scene(ds, field=("grid", "nbody_mass"))
    sc.grey_opacity = True

    source = sc[0]
    cam = sc.camera
    # for resolution purposes
    source.set_use_ghost_zones(True)
    
    # set initial camera setting
    cam.width = width * zoom_factor
    
    ad = ds.all_data()
    mi, ma = ad.quantities.extrema("nbody_mass")
    
    tfh = TransferFunctionHelper(ds)
    tfh.set_field(("grid", "nbody_mass"))
    tfh.set_log(True)
    tfh.set_bounds()
    tfh.build_transfer_function()
    
    tf = yt.ColorTransferFunction((np.log10(mi), np.log10(ma)))
    tfh.tf.add_layers(16, colormap="cmyt.arbre")

    render_source = sc.get_source()
    render_source.transfer_function = tfh.tf
    
    counter=0
    for pos in range(len(prog_positions)-2):          #yt.parallel_objects(prog_snaps, num_procs):        # loop through snapshots and plot the halo's progenitors
        counter+=1
        cam.position = prog_positions[pos]
        cam.focus = prog_positions[pos+2]
        print(cam.focus)
        # increase the resolution
        sc.camera.resolution = (1024, 1024)
        
        # save an image at the starting position
        sc.save(f"position_zoom={zoom_factor}_frames={frames}_counter={counter}_far={how_far}_order={how_large}.png", sigma_clip=4., render=False)    
    print(np.array(final_position),prog_positions[-1])

move()