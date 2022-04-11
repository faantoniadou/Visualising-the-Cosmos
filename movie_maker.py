import yt
#yt.enable_parallelism()
import ytree
from yt import YTArray
from yt.visualization.volume_rendering.api import Scene, Camera, create_volume_source
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.render_source import VolumeSource

import unyt
from unyt import Mpc, kpc
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import communication_system
from yt.utilities.cosmology import Cosmology

import time
import pandas as pd
import numpy as np
from numpy.polynomial import polynomial
import datetime
import math
import warnings

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.animation as animation
import matplotlib.image as mpimg
from matplotlib.pyplot import figure
from matplotlib import rc_context, cm
from matplotlib.colors import Normalize 
    
from tqdm import tqdm
# to create smoother camera path:
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline, interpn
from scipy.stats import gaussian_kde

from IPython.display import Image

from main import ExploreHalo


co = Cosmology()


today = datetime.datetime.now()


def make_grid_name(fnumber=int(101)):
    # this makes the file name for the covering grids directory
    return yt.load(f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{fnumber:03}_covering_grid.h5")


class Setting:
    '''
    We need input to be path positions
    for rotarion we do it in here at a centre (render once)
    for time evolution we need to render snaps too 
    for moving between points we change focus and need path input
    '''
    def __init__(self, zoom_factor=5., width=1., units="unitary", fnumber=101):           ### check if units are right
        self.field = ("grid", "nbody_mass")
        self.fnumber = fnumber
        self.ds = make_grid_name(self.fnumber) #
        
        self.left_corner = self.ds.domain_left_edge
        self.right_corner = self.ds.domain_right_edge
        self.ad = self.ds.all_data()
        self.mi, self.ma = self.ad.quantities.extrema("nbody_mass")
        self.zoom_factor = zoom_factor

        self.width = self.ds.quan(width, self.units)
        self.units = units
        self.path = None
        
        self.sc = yt.create_scene(self.ds, field=self.field)
        self.source = self.sc[0]
        self.cam = self.sc.add_camera(self.ds)

        # make it pretty
        self.scene.grey_opacity = True
        # for resolution purposes
        self.source.set_use_ghost_zones(True)
        self.cam.width = self.width * self.zoom_factor


    def transfer_save(self, filename, s_clip=4., nlayers=10, cmap="cmyt.arbre"):
        tfh = TransferFunctionHelper(self.ds)
        tfh.set_field(self.field)

        tfh.set_log(True)
        tfh.set_bounds()
        tfh.build_transfer_function()
        
        # Set up transfer function
        tf = yt.ColorTransferFunction((np.log10(self.mi), np.log10(self.ma)))
        tf.add_layers(nlayers, colormap=cmap)
        
        render_source = self.sc.get_source()
        render_source.transfer_function = tfh.tf
        
        self.scene.render()
        self.sc.camera.resolution = (1024, 1024)
        self.scene.save(filename, sigma_clip=s_clip)



class Animations(Setting, ExploreHalo):
    # class with Functions for making animations

    def __init__(self, nlayers=10, fnumber=101):
        Setting.__init__(self, zoom_factor=self.zoom_factor, fnumber=self.fnumber)
        ExploreHalo.__init__(self)
        
        self.fnumber = fnumber
        self.zoom_factor = 5.
        self.nlayers = nlayers
        

    def single_halo(self, nlayers=10):
        halo_position, halo_radius, halo_mass, prog_snaps, prog_positions = self.get_properties()
        self.cam.focus = halo_position
        self.zoom_factor = float(input("Zoom factor: "))
        self.transfer_save(f"position_zoom={zoom_factor}_{self.quantity}={self.how_large}_colours={no_layers}.png",
                           s_clip=4., nlayers=self.nlayers)
        
        
    def linear_path(self, how_far=10, no_steps=200):
        # makes a smooth camera path between 2 points

        initial_position, width, prog_snaps, final_position = self.get_neighbour(how_far)
        
        
        x_initial, y_initial, z_initial = initial_position[0], initial_position[1], initial_position[2]
        x_final, y_final, z_final = final_position[0], final_position[1], final_position[2]
        
        x_path, y_path, z_path = np.linspace(x_initial, x_final, no_steps), np.linspace(y_initial, y_final, no_steps),\
            np.linspace(z_initial, z_final, no_steps)
        
        camera_path = np.column_stack((x_path, y_path, z_path))

        return initial_position, width, prog_snaps, final_position, camera_path


    def nav_movie(self, fnumber=101):
        # implements nav_progs.py file, navigating between neighbouring progenitors
        
        no_steps = int(input("Number of frames to move through"))
        this_far = float(input("How far should I search (in units of the halo's radius)? "))
        final_position, width, all_positions, prog_snaps, prog_positions = self.linear_path(this_far, no_steps)
        
        counter=0
        for pos in prog_positions:          #yt.parallel_objects(prog_snaps, num_procs):        # loop through snapshots and plot the halo's progenitors
            counter+=1
            self.cam.focus = pos
            # increase the resolution
            
            self.transfer_save((f"position_zoom={zoom_factor}_{self.quantity}={self.how_large}_frames={no_steps}_counter={counter}_far={this_far}.png"),
                               s_clip=4., nlayers=self.nlayers)


    def smooth_path(self, positions):
        # this returns a smooth camera path provided all positions we want the camera to pass through
        window_size = 2*int((len(positions)/4)) + 1
        pol_degree = 8

        smoother_pos_x = savgol_filter(positions[:,0], window_size, pol_degree)
        smoother_pos_y = savgol_filter(positions[:,1], window_size, pol_degree)
        smoother_pos_z = savgol_filter(positions[:,2], window_size, pol_degree)
        
        camera_path = np.column_stack((smoother_pos_x, smoother_pos_y, smoother_pos_z))
        return camera_path


    def time_movie(self, starting_snap=101):
        # movie of progenitors
        halo_position, halo_radius, halo_mass, prog_snaps, prog_positions = self.get_properties()
        prog_positions = self.smooth_path(prog_positions)
        self.zoom_factor = float(input("Zoom factor: "))
        
        for i in range(len(prog_snaps)-starting_snap):#yt.parallel_objects(prog_snaps, num_procs):        # loop through snapshots and plot the halo's progenitors
            # Load the dataset
            self.fnumber = prog_snaps[i+starting_snap]
            self.cam.focus = prog_positions[i]
            # increase the resolution
            
            self.transfer_save(f"progenitor_{i+starting_snap}_{self.quantity}_order={self.how_large}_zoomout={zoom_factor}.png",
                               s_clip=4., nlayers=self.nlayers)

    
    def rotate_halo(snap_no=101, degrees=2*np.pi, nframes=100):
        halo_position, halo_radius, halo_mass, prog_snaps, prog_positions = self.get_properties()
        
        # we focus on the halo we want from the beginning
        final_position, final_width = get_halo_centre()

        # set initial camera setting
        self.cam.width = halo_radius * 2
        self.cam.focus = halo_position
        
        '''
        rotating process begins here
        '''
        # save an image at the starting position
        frame = 0

        self.transfer_save(f"rotation_{frame}_{self.quantity}={self.how_large}.png", s_clip=4., nlayers=self.nlayers)
        
        for _ in cam.iter_rotate(degrees, nframes):
            frame += 1
            self.transfer_save(f"rotation_{frame}_{self.quantity}={self.how_large}.png", s_clip=4., nlayers=self.nlayers)


