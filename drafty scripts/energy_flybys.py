'''
This script is used to move between haloes. We will do this by finding all haloes above a specific threshold (e.g. above a specific mass)
and then looping through to identify pairs of haloes with similar properties.
'''
import yt
#yt.enable_parallelism()
# parallel = False
import ytree
import time
import glob
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.animation as animation
import matplotlib.image as mpimg

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
from scipy.signal import savgol_filter
from matplotlib.pyplot import figure
from yt.utilities.cosmology import Cosmology

co = Cosmology()
H = co.hubble_parameter(1.0).in_units("unitary/s/unitary")
    
    

today = datetime.datetime.now()

quantity = 'mass' #str(input('Maximum quantity to search for in halo catalogue: '))

how_large = 40 #int(input(f'Order to largest {quantity} (e.g. 1 is maximum and 2 is second to largest): '))

G_value = (6.67e-11)#.in_units('m**3 * kg ** (-1) * s ** (-2)').to('unitary**3 * kg * s ** (-2)')

# Load the dataset
# this makes the file name for the covering grids directory
def make_grid_name(fnumber):
    grid_name = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{fnumber:03}_covering_grid.h5"
    return grid_name



def get_flybys():
    '''
    This finds flybys based on criteria
    Primary halo is chosen based on mass
    '''
    
    # Finds halo with the maximum of the quantity (field) the user provides.
    a = ytree.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5")
    
    values = (a[quantity]).tolist()
    values.sort()
    i = np.where(a[quantity] == (values[-how_large]))[0][0]           #['mass'])
    
    print("This is different")
    # i = np.abs(np.array(a['mass'].to('Msun')) - 1.5e12).argmin()

    #i = np.argmax(a[quantity])
    my_tree = a[i]
    
    initial_position = my_tree["position"].to("unitary").value
    initial_radius = my_tree["virial_radius"].value
    primary_mass = my_tree["mass"].to('Msun').value
    print("{:e}".format(my_tree["mass"].to('Msun')))
    ds = a.ytds
    
    
    G = ds.quan(G_value, 'm**3 * kg ** (-1) * s ** (-2)').to('unitary**3 * Msun ** (-1) * s ** (-2)')
    
    
    sphere = ds.sphere(initial_position, (initial_radius*10, "unitary"))
    #sphere.force_periodicty() 
   
    prog_positions = np.array(my_tree["prog", "position"].to("unitary"))
    main_vel = my_tree['velocity'].to('unitary/s').value
    
    for node in a.get_nodes_from_selection(sphere):
        #print (node["position"])
        node_redshifts = node["prog", "redshift"]
        node_prog_positions = np.array(node["prog", "position"].to("unitary"))
        node_radius = node["virial_radius"].value
        node_snaps = node["prog", "Snap_idx"]
        node_mass = node["mass"].to('Msun').value
        node_vel = node['velocity'].to('unitary/s').value
        #print(node_prog_positions)
        distances = []
        secondary_x, secondary_y, secondary_z = [], [], []           # probably redundant
        snaps = []
        primary_x, primary_y, primary_z = [], [], []
        R12 = np.array(node['position'].to('unitary') - my_tree["position"].to('unitary'))
        
        E = node_mass * primary_mass * (1/2 * (np.linalg.norm((node['velocity'].to('unitary/s').value - 
                                                              my_tree['velocity'].to('unitary/s').value + 
                                                              H.value * R12)))**2/(primary_mass + node_mass))
        
        GPE = node_mass * primary_mass * G.value/(np.linalg.norm(R12))
        
        # GET BINDING ENERGY BETWEEN HALOS
        E_binding = E - GPE
        print(E, GPE)
        print(E_binding)
        #node_mass * primary_mass * (1/2 * np.linalg.norm((node['velocity'].to('m/s').value - my_tree['velocity'].to('m/s').value + H.value * R12.value))**2/(primary_mass + node_mass) - 6.67e-11/(np.linalg.norm(R12.value)))
        if E_binding < 0:
            if node['mass'].value >= 0.1 * primary_mass:        # do this for haloes with a mass comparable to that of the primary one
                print('Suspected flyby found!')
                pair_dist = np.linalg.norm(np.array(node['position'].to("unitary")) - np.array(initial_position))
                            
                # find distance between progenitors of node and primary halo for as long as the node existed
                for i, (main_pos, node_pos) in enumerate(zip(prog_positions, node_prog_positions)):
                    # let us account for periodic boundary conditions:
                    vec = main_pos - node_pos
                    box_size = 1
                    abs_dist = np.linalg.norm(np.subtract(np.mod(np.add(vec, np.multiply(box_size, box_size/2)), box_size),
                                                        np.multiply(box_size, box_size/2)))
                    distances.append(abs_dist)
                    
                    secondary_x.append(node_pos[0])
                    secondary_y.append(node_pos[1])
                    secondary_z.append(node_pos[2])
                    
                    snaps.append(node_snaps[i])
                    
                    primary_x.append(main_pos[0])
                    primary_y.append(main_pos[1])
                    primary_z.append(main_pos[2])
        
                
                df = pd.DataFrame({"x": primary_x, "y": primary_y, "z": primary_z, "x'": secondary_x, "y'": secondary_y, "z'": secondary_z, 
                                "Distances" : distances, "Snapshot" : snaps})
                df.to_csv(f'flyby_dists_{node}_{quantity}={how_large}_energy.csv', index=False)
                
        

get_flybys()




def make_path():
    '''
    Makes a smooth camera path
    '''
    no_steps = frames
    initial_position, width, prog_snaps, final_position = get_halo_centre()
    
    x_initial, y_initial, z_initial = np.array(initial_position[0]), np.array(initial_position[1]), np.array(initial_position[2])
    x_final, y_final, z_final = np.array(final_position[0]), np.array(final_position[1]), np.array(final_position[2])
    
    x_path, y_path, z_path = np.linspace(x_initial, x_final, no_steps), np.linspace(y_initial, y_final, no_steps), \
        np.linspace(z_initial, z_final, no_steps)
        

    
    camera_path = np.vstack((x_path, y_path, z_path)).T

    
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
    
    tfh = TransferFunctionHelper(ds)
    tfh.set_field(("grid", "nbody_mass"))
    tfh.set_log(True)
    tfh.set_bounds()
    tfh.build_transfer_function()
    ad = ds.all_data()
    mi, ma = ad.quantities.extrema("nbody_mass")
    tf = yt.ColorTransferFunction((np.log10(mi), np.log10(ma)))
    tfh.tf.add_layers(10, colormap="cmyt.arbre")
    
    render_source = sc.get_source()
    render_source.transfer_function = tfh.tf
    
    counter=0
    for pos in prog_positions:          #yt.parallel_objects(prog_snaps, num_procs):        # loop through snapshots and plot the halo's progenitors
        counter+=1
        cam.focus = pos
        print(cam.focus)
        # increase the resolution
        sc.camera.resolution = (1024, 1024)
        
        # save an image at the starting position
        sc.save(f"position_zoom={zoom_factor}_frames={frames}_counter={counter}_far={how_far}_order={how_large}.png", sigma_clip=4.)    
    print(np.array(final_position),prog_positions[-1])

# move()