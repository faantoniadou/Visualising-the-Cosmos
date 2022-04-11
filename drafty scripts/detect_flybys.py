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

    
    

today = datetime.datetime.now()

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
    quantity = str(input('Maximum quantity to search for in halo catalogue: '))

    how_large = int(input(f'Order to largest {quantity} (e.g. 1 is maximum and 2 is second to largest): '))
    
    # Finds halo with the maximum of the quantity (field) the user provides.
    a = ytree.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5")
    
    values = (a[quantity]).tolist()
    values.sort()
    i = np.where(a[quantity] == (values[-how_large]))[0][0]           #['mass'])
    
    
    # i = np.abs(np.array(a['mass'].to('Msun')) - 1.5e12).argmin()

    #i = np.argmax(a[quantity])
    my_tree = a[i]
    
    initial_position = my_tree["position"]
    initial_radius = (my_tree["virial_radius"].to("unitary")).value
    primary_radii = (my_tree["prog", "virial_radius"].to("unitary")).value

    primary_mass = my_tree["mass"].value
    print("{:e}".format(my_tree["mass"].to('Msun')))
    ds = a.ytds
    sphere = ds.sphere(initial_position, (initial_radius*10, "unitary"))
        
    prog_positions = np.array(my_tree["prog", "position"].to("unitary"))
    # found = 0 
    for node in a.get_nodes_from_selection(sphere):
        #print (node["position"])
        node_redshifts = node["prog", "redshift"]
        node_prog_positions = np.array(node["prog", "position"].to("unitary"))
        node_radii = (node["prog", "virial_radius"].to("unitary")).value
        node_snaps = node["prog", "Snap_idx"]
        #print(node_prog_positions)
        
        distances, secondary_x, secondary_y, secondary_z = [], [], [], []           # can probably be simplified
        main_radii, radii, redshifts, positions, snaps = [], [], [], [], []
        primary_x, primary_y, primary_z = [], [], []
        
        
        
        if node['mass'].value >= 0.1 * primary_mass:        # do this for haloes with a mass comparable to that of the primary one
            pair_dist = np.linalg.norm(np.array(node['position'].to("unitary")) - np.array(initial_position))

            # find distance between progenitors of node and primary halo for as long as the node existed
            for i, (main_pos, node_pos) in enumerate(zip(prog_positions, node_prog_positions)):
                # let us account for periodic boundary conditions:
                vec = main_pos - node_pos
                box_size = 1
                abs_dist = np.linalg.norm(np.subtract(np.mod(np.add(vec, np.multiply(box_size, box_size/2)), box_size),
                                                      np.multiply(box_size, box_size/2)))
                distances.append(abs_dist)
                positions.append(node_pos)
                snaps.append(node_snaps[i])
                radii.append(node_radii[i])
                main_radii.append(primary_radii[i])
                redshifts.append(node_redshifts[i])
                
                # for pandas dataframe purposes:
                
                secondary_x.append(node_pos[0])
                secondary_y.append(node_pos[1])
                secondary_z.append(node_pos[2])               
                
                primary_x.append(main_pos[0])
                primary_y.append(main_pos[1])
                primary_z.append(main_pos[2])
            
            print("finished")
            print(np.any(np.r_[True, np.array(distances)[1:] < np.array(distances)[:-1]] & np.r_[np.array(distances)[:-1] < np.array(distances)[1:], True]==True))
            # for i in range(len(distances)):
            for j in range(3, len(distances)-2):
                # this checks if there are local minima in the distance array
                if np.any(np.r_[True, np.array(distances)[1:] < np.array(distances)[:-1]] & 
                          np.r_[np.array(distances)[:-1] < np.array(distances)[1:], True]==True):
                # if distances[i] <= (node_radii[i] + main_radii[i]):     # for internal flybys
                # if (distances[j-1] > distances[j] and distances[j-2] > distances[j]
                #     and distances[j] < distances[j+1] and distances[j-1] < distances[j+2]):
            # if np.any(distances < initial_radius + node_radius):
                    # if found <= 3:
                    print('Suspected flyby found!')
                    # found += 1
                    # print(found)
                    df = pd.DataFrame({"x": primary_x, "y": primary_y, "z": primary_z, "x'": secondary_x, "y'": secondary_y, "z'": secondary_z, 
                                    "Distances" : distances, "Snapshot" : snaps, "Flyby Radius": radii, "Primary Radius": main_radii, 
                                    "Redshift" : redshifts})
                    name = f'{node}_{quantity}={how_large}_dist.csv'
                    print(f"Saving as {name}")
                    
                    df.to_csv(name, index=False)
                    break
                # break
                

        

def get_flyby_positions(filename):
    '''
    This function reads the csv files produced in flyby detection and produces smooth paths
    '''
    df = pd.read_csv(filename)
    
    radius = np.array(df["Primary Radius"])[0]          # to be used in width
    redshifts = np.array(df['Redshift'])
    snaps = np.arange(101) #np.array(df['Snapshot'])
    
    x_posns = np.array(df["x'"])
    y_posns = np.array(df["y'"])
    z_posns = np.array(df["z'"])
    
    x_main = np.array(df["x"])
    y_main = np.array(df["y"])
    z_main = np.array(df["z"])
    
    
    
    window_size = 2*int((len(x_posns)/4)) + 1            # window size must always be odd
    pol_degree = 8
    
    camera_path = np.empty((len(x_posns),3))

    rel_x = np.subtract(np.mod(np.add(x_main - x_posns, np.multiply(1, 1/2)), 1), np.multiply(1, 1/2))
    rel_y = np.subtract(np.mod(np.add(y_main - y_posns, np.multiply(1, 1/2)), 1), np.multiply(1, 1/2))
    rel_z = np.subtract(np.mod(np.add(z_main - z_posns, np.multiply(1, 1/2)), 1), np.multiply(1, 1/2))
    
    window_size = window_size = 2*int((len(rel_x)/4)) + 1
    pol_degree = 8
    
    resultant_smooth_x = savgol_filter(rel_x, window_size, pol_degree)
    resultant_smooth_y = savgol_filter(rel_y, window_size, pol_degree)
    resultant_smooth_z = savgol_filter(rel_z, window_size, pol_degree)
    
    smoother_pos_x = savgol_filter(x_posns, window_size, pol_degree)
    smoother_pos_y = savgol_filter(y_posns, window_size, pol_degree)
    smoother_pos_z = savgol_filter(z_posns, window_size, pol_degree)
    
    camera_path[:,0], camera_path[:,1], camera_path[:,2] = smoother_pos_x, smoother_pos_y, smoother_pos_z
    
    width = radius * 2
    
    plt.plot(redshifts, resultant_smooth_x, label='smoothed x-coordinate',linewidth=1, marker='.')
    plt.xlabel('Redshift')
    plt.title('X-coordinate of the flyby relative to the primary halo')
    plt.ylabel('Relative Position')
    plt.savefig(f'X_{filename}_path_boundary.png')
    plt.clf()
    
    plt.plot(redshifts, resultant_smooth_y, label='smoothed y-coordinate',linewidth=1, marker='.')
    plt.xlabel('Redshift')
    plt.title('Y-coordinate of the flyby relative to the primary halo')
    plt.ylabel('Relative Position')
    plt.savefig(f'Y_{filename}_path_boundary.png')
    plt.clf()
    
    plt.plot(redshifts, resultant_smooth_z, label='smoothed z-coordinate',linewidth=1, marker='.')
    plt.xlabel('Redshift')
    plt.title('Z-coordinate of the flyby relative to the primary halo')
    plt.ylabel('Relative Position')
    plt.savefig(f'Z_{filename}_path_boundary.png')
    plt.clf()
    
    
    
    # all together:
    plt.plot(redshifts, resultant_smooth_x, label='x-coordinate',linewidth=1, marker='.')
    plt.plot(redshifts, resultant_smooth_y, label='y-coordinate',linewidth=1, marker='.')
    plt.plot(redshifts, resultant_smooth_z, label='z-coordinate',linewidth=1, marker='.')
    
    plt.xlabel('Redshift')
    plt.ylabel('Relative Position')
    plt.legend()
    # plt.ylim(np.min((np.min(resultant_smooth_x), np.min(resultant_smooth_y), np.min(resultant_smooth_z))),
            #  np.max((np.max(resultant_smooth_x), np.max(resultant_smooth_y), np.max(resultant_smooth_z))))
    plt.title('Flyby position relative to the primary halo')
    plt.savefig(f'{filename}_path_boundary.png')
    
    plt.clf()
    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111, projection='3d')

    # Data for a three-dimensional line
    zline = resultant_smooth_z
    xline = resultant_smooth_x
    yline = resultant_smooth_y
    img = ax.scatter3D(xline, yline, zline, c=redshifts, cmap='jet', marker='o', label='smoothed path')
    ax.plot3D(xline, yline, zline)
    
    # Data for three-dimensional scattered points
    zdata = z_posns - z_main
    xdata = x_posns - x_main
    ydata = y_posns - y_main
    
    ax.set_xlabel('X Path')
    ax.set_ylabel('Y Path')
    ax.set_zlabel('Z Path')

    img = ax.scatter3D(rel_x, rel_y, rel_z, c=redshifts, cmap='jet', marker='x', label='actual path')
    # ax.plot3D(rel_x, rel_y, rel_z)
    
    cbar = fig.colorbar(img, fraction=0.03, pad=0.04)
    cbar.set_label('Redshift')
    plt.legend(loc=2, prop={'size': 9}, bbox_to_anchor=(0.1,0.9))
    fig.savefig(f'{name}_3Dpath_boundary.png')

    return snaps, width, camera_path






#num_procs = 4       # number of processors to use
def time_haloes(filename):
    nlayers = int(input('Number of colour layers: '))
    zoom_factor = float(input('Zoom factor: '))

    # if we want to start from a specific redshift form the snap catalogue onwards
    starting_snap = int(input("Starting snap: "))
    
    prog_snaps, width, prog_positions = get_flyby_positions(filename)
    
    for i in range(len(prog_snaps)-starting_snap):#yt.parallel_objects(prog_snaps, num_procs):        # loop through snapshots and plot the halo's progenitors
        # Load the dataset
        ds = yt.load(make_grid_name(prog_snaps[100-i+starting_snap]))
        
        # Create a volume rendering
        sc = yt.create_scene(ds, field=("grid", "nbody_mass"))
        sc.grey_opacity = True

        source = sc[0]
        cam = sc.camera
        #cam = sc.add_camera(ds, lens_type='perspective')
        # for resolution purposes
        source.set_use_ghost_zones(True)

        # set initial camera setting
        cam.width = width * zoom_factor
        cam.focus = prog_positions[i+starting_snap]
        
        # increase the resolution
        sc.camera.resolution = (1024, 1024)
        
        # Set up a custom transfer function using the TransferFunctionHelper.
        # We use 10 Gaussians evenly spaced logarithmically between the min and max
        # field values.
        
        #we bump up our minimum to cut out some of the background
        
        ad = ds.all_data()
        mi, ma = ad.quantities.extrema(("grid", "nbody_mass"))
        
        tfh = TransferFunctionHelper(ds)
        tfh.set_field()
        tfh.set_log(True)
        tfh.set_bounds()
        tfh.build_transfer_function()
        tf = yt.ColorTransferFunction((np.log10(mi), np.log10(ma)))
        tfh.tf.add_layers(nlayers, colormap="cmyt.arbre")

        # Grab the first render source and set it to use the new transfer function
        render_source = sc.get_source()
        render_source.transfer_function = tfh.tf
        
        # save an image at the starting position
        #sc.annotate_axes()
        sc.save(f"{filename}_{quantity}={how_large}.png", sigma_clip=4.)    
       
 
 
get_flybys()
# name =  str(input('Flybys file name: ')) #'flyby_dists_TreeNode[621479483]_mass=6000_dist.csv'
# get_flyby_positions(name)
# time_haloes(name)