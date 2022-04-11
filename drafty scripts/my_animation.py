import numpy as np
import math
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.animation as animation
import matplotlib.image as mpimg
import yt
from yt import YTArray
from yt.visualization.volume_rendering.api import Scene, Camera
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.render_source import VolumeSource
import datetime
import unyt
from unyt import Mpc, kpc
from matplotlib import rc_context

#yt.enable_parallelism()

# this makes the file name for the covering grids directory
def make_grid_name(frame):
    grid_name = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{frame:03}_covering_grid.h5"
    return grid_name

# this makes file name for halo trees 
def make_halo_name():
    halo_name = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5"
    return halo_name



class Setting(object):

    def __init__ (self, fname):

        self.field = ("grid", "nbody_mass")
        self.ds = yt.load(fname)
        self.left_corner = self.ds.domain_left_edge
        self.right_corner = self.ds.domain_right_edge
        self.ad = self.ds.all_data()
        self.mi, self.ma = self.ad.quantities.extrema("nbody_mass")

        # Create the region
        self.region = self.ds.box(self.left_corner, self.right_corner)
        
        self.scene = yt.create_scene(self.ds, lens_type="perspective", field=self.field)
        #self.scene.annotate_domain(self.grid, [1,1,1, 0.1])
        self.cam = self.scene.add_camera(self.ds, lens_type="perspective")

        #Image array to save images for the Animation
        self.ims = []

        # make it pretty
        self.scene.grey_opacity = True


    def set_camera (self, pos_vector, focus_vector, zoom_factor):
        # Set up the camera parameters: focus, width, resolution, and image orientation

        self.cam.position = self.ds.arr(pos_vector, "unitary")
        self.cam.focus = self.ds.arr(focus_vector, "unitary")
        # if we want focus to be at the centre of the domain we must do:
        # self.cam.focus = self.ds.domain_center
        #self.cam.zoom(zoom_factor)


    # rotate around a fixed point, by a given angle, in a given number of frames
    def rotate(self, angle, frames, centre, s_clip):

        self.fig = plt.figure()
        #It saves the image each time, before loading them into the self.ims list
        for _ in self.cam.iter_rotate(angle, frames, [0,0,1], centre):          # note: angle is in radians
            self.scene.save("Temporary", sigma_clip = s_clip)
            im = plt.imshow(mpimg.imread("Temporary.png"), animated = False)
            self.ims.append([im])


    #Function to animate a zoom in operation onto the focus
    def zoom(self, final_zoom, frames, s_clip):

        self.fig = plt.figure()
        #Same process as in Rotate
        for _ in self.cam.iter_zoom(final_zoom, frames):
            self.scene.save("Temporary", sigma_clip = s_clip)
            im = plt.imshow(mpimg.imread("Temporary.png"), animated = False)
            self.ims.append([im])


    #Moves the camera and animates in the process
    def move(self, final_pos, frames, s_clip):

        self.fig = plt.figure()
        final_pos = self.grid.arr(final_pos, 'code_length')
        #Same process as in Rotate
        for _ in self.cam.iter_move(final_pos, frames, exponential = True):
            self.scene.save("Temporary", sigma_clip = s_clip)
            im = plt.imshow(mpimg.imread("Temporary.png"), animated = False)
            self.ims.append([im])
            
    #Method that will create and show animation of what the moving camera sees up to that point
    #Frozen because there is no time evolution
    def frozen_animation(self, fps, fname):

        self.cam.switch_orientation([0,1,0], [0,0,1])
        ani = animation.ArtistAnimation(self.fig, self.ims, interval = 63, blit = False, repeat_delay = 0)
        mywriter = animation.writers['ffmpeg']
        writer = mywriter(fps = 16, metadata = dict(artist='Me'), bitrate = 1000)
        ani.save(fname, writer = writer)
        
        #User can mess around with various transfer functions
    def render_save (self, s_clip, fname):
        tfh = TransferFunctionHelper(self.ds)
        tfh.set_field(self.field)
        tfh.set_log(True)
        tfh.set_bounds()
        tfh.build_transfer_function()
        
        # Set up transfer function
        tf = yt.ColorTransferFunction((np.log10(self.mi), np.log10(self.ma)))
        tf.add_layers(6, w=0.05)
        
        self.scene.render()
        self.scene.save(fname, sigma_clip = s_clip)
    
    

#Class with Functions for making time-evolving animations
class Animation(object):

    def __init__(self, snapshot):
        self.fig = plt.figure()
        self.ims = []
        self.grid_nos = snapshot

    # simply animates the box (focus is by default at the centre of the box)
    def simple_animation(self, s_clip, cam_position, focus):

        for frame in range (self.grid_nos[0], self.grid_nos[1]+1):

            fname = make_grid_name(frame)
            
            sc = Setting(fname)

            sc.set_camera(cam_position, 1, focus)

            sc.scene.render_save(s_clip, self.grid_no)


    def save_animation(self, fps, grid_no):
        anim = animation.ArtistAnimation(self.fig, self.ims, interval = 63, blit = False, repeat_delay = 0)

        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)
        anim.save(f'{grid_no}.mp4', writer = writer)



def main():
    first = int(input('First grid number: '))
    last = int(input('Last grid number: '))
    
    grid_nos = np.arange(first,last+1,1).astype(int)
    fps = int(input('Frames per second: '))
    A1 = Animation(grid_nos)
    A1.simple_animation(4, [0,0,0], [0,0,1])
    A1.save_animation(fps, grid_no)


main()











































    # #Render only the members of the tree that are the progenitors of the given Halo_id. i.e. Animate the evolution
    # #of a particular Tree
    # def render_tree(self, domain, haloID, res, dist, s_clip, haloDir):
    #     #Return a list of lists containing the tree
    #     Pos_list = Most_Massive_Progenitor(haloID, domain[1], domain[0], haloDir)
    #     Pos_list.reverse()
    #     sources = Make_Tree(haloID, domain[1], haloDir)

    #     for frame in range (domain[0], domain[1]+1):
    #         index = frame - domain[0]

    #         fname = make_fname(frame, self.dirname)
    #         #Create the Box and all that
    #         sc = Setting(fname, res, sources[frame-1], self.is_big)

    #         if sources[frame-1] != []:
    #             #Set the camera postion and focus
    #             #But first, transform the Camera-postion and camera-focus vectors in unitary dimensions
    #             #dist must have the same dimensions as the box (This is essential as the positions have to be tranformed in Unitary dimaensions)
    #             Cam_Pos = np.array(Pos_list[index].to(dist.units))-np.array([0,-dist, 0])
    #             Focus_Pos = np.array(Pos_list[index].to(dist.units))
    #             Cam_Pos = Cam_Pos/sc.domain_length.value
    #             Focus_Pos = Focus_Pos/sc.domain_length.value

    #         else:
    #             Cam_Pos = np.array([0.5, 0, 0.5])
    #             Focus_Pos = np.array([0.5, 0.5, 0.5])

    #         #Set the camera postion and focus
    #         sc.set_camera(Cam_Pos, [0, 0, 1], Focus_Pos, 1)

    #         #If running on your computer, disable this line and enable the next 3 commands
    #         sc.scene.save('snapshots2/'+str(frame), sigma_clip=s_clip)

    #         #Save the image and store it in the list with the others
    #         #Sc.scene.save_annotated("Temporary", sigma_clip = 4.0, text_annotate = [[(.1, .1), text_str]])
    #         #im = plt.imshow(mpimg.imread("Temporary.png"), animated = False)
    #         #self.ims.append([im])


    # #with domain = (45, 50), it will save frames 45 to 50, included
    # #haloID must be that at the last frame of the domain
    # #dist = distance from the Halo. The camera is moved away from the Halo along the y direction. This can be changed
    # #Halo_directory = directory where hlists and scales are.
    # def Follow_Most_Massive(self, domain, Halo_id, res, dist, s_clip, Halo_directory):

    #     Pos_list = Most_Massive_Progenitor(haloID, domain[1], domain[0], Halo_directory)
    #     Pos_list.reverse()

    #     for frame in range (domain[0], domain[1]+1):
    #         index = frame - domain[0]
    #         fname = make_fname(frame, self.dirname)

    #         #Create_Box and all that
    #         sc = Setting(fname)

    #         #Set the camera postion and focus
    #         #But first, transform the Camera-postion and camera-focus vectors in unitary dimensions
    #         #Assuming that dist has the same dimensions as the box
    #         Cam_Pos = np.array(Pos_list[index].to(dist.units))-np.array([0,-dist, 0])
    #         Focus_Pos = np.array(Pos_list[index].to(dist.units))
    #         Cam_Pos = Cam_Pos/sc.domain_length.value
    #         Focus_Pos = Focus_Pos/sc.domain_length.value

    #         #Set the camera postion and focus
    #         sc.Set_Cam(Cam_Pos, [0, 0, 1], Focus_Pos, 1)

    #         #If running on your computer, disable this line and enable the next 3 commands
    #         sc.scene.save('snapshots2/'+str(frame), sigma_clip = s_clip)

    #         #Save the image and store it in the list with the others
    #         #Sc.scene.save("Temporary", sigma_clip = s_clip)
    #         #im = plt.imshow(mpimg.imread("Temporary.png"), animated = False)
    #         #self.ims.append([im])