import yt
#yt.enable_parallelism()
import ytree
import time
import glob
import pandas as pd
import numpy as np
import datetime
import math
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.animation as animation
import matplotlib.image as mpimg
from matplotlib.pyplot import figure
from matplotlib import rc_context
import warnings

from yt import YTArray
from yt.visualization.volume_rendering.api import Scene, Camera, create_volume_source
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.render_source import VolumeSource
import unyt
from unyt import Mpc, kpc
import yt
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import communication_system
    
from tqdm import tqdm
# to create smoother camera path:
from scipy.signal import savgol_filter
from numpy.polynomial import polynomial
from scipy.interpolate import UnivariateSpline

from matplotlib.pyplot import figure
from yt.utilities.cosmology import Cosmology
from IPython.display import Image
from scipy.stats import gaussian_kde
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn

plt.style.use("seaborn-whitegrid")

co = Cosmology()

# some constants
H = co.hubble_parameter(1.0).in_units("unitary/s/unitary")
G_value = 6.67e-11


today = datetime.datetime.now()


def make_grid_name(fnumber=int(101)):
    # this makes the file name for the covering grids directory
    return yt.load(f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{fnumber:03}_covering_grid.h5")
    

class ExploreHalo:
    def __init__(self, quantity='mass'):
        # some user inputs to be used in subsequent methods
        # quantity maximum halo to find:
        self.quantity = quantity #str(input('Maximum quantity to search for in halo catalogue: '))
        self.how_large = int(input(f'Order to largest {quantity} (e.g. 1 is maximum and 2 is second to largest): '))
        #self.grid_name = make_grid_name(101) #
        self.a = ytree.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5")        # load tree
        self.width = None
        self.filename = str(input("File name: ")) #'TreeNode[615163324]_mass=42_dist.csv' #None
        self.field = ("grid", "nbody_mass")
        self.ds = make_grid_name(101)
        

    def basics(self):

        print (f'Box size: {self.a.box_size}')
        print (f'Cosmological parameters: {self.a.hubble_constant, self.a.omega_matter, self.a.omega_lambda}')
        print (f'Field list: {self.a.field_list}')
        print (f'Derived field list: {self.a.derived_field_list}')
        # print (f'Index: {self.a["Snap_idx"]}')
        print (self.a.field_info['Mvir_all'])

        
        #print (f'Node field list: {my_tree.field_list}')

        # how many trees are there?
        # print (f'Number of trees: {a.size}')

        # # Individual trees can be accessed by indexing the Arbor object.
        # print(f'Individual tree: {self.a[0]}')
        # print(f'Position: {self.a["position"]}')

        # # A TreeNode is one halo in a merger tree. The number is the universal identifier associated with halo. 
        # # It is unique to the whole arbor. Fields can be accessed for any given TreeNode in the same dictionary-like fashion.
        # print(f'Mass of the tree: {my_tree["mass"]}')

        # Array slicing can also be used to select multiple TreeNode objects. This will return a generator that can be iterated over or cast to a list.
        # every_second_tree = list(a[::2])
        # print (every_second_tree[0]["mass"])

        # print (f'Access nodes: {my_tree["tree"]}')
        # # loop over nodes
        # for my_node in my_tree["tree"]:
        #     print (my_node, my_node["mass"])
            
        #print (f'All nodes in the tree: {list(my_tree["tree"])}')
        #print (f'Some nodes in the tree: {list(my_tree["tree"])[0:5]}')

        #A list of defined vector fields can be seen by doing:
        # print(f'Vector fields: {a.field_info.vector_fields}')
        # print(f'Position: {a["position"]}')
        # print(f"Mass: {a['mass']}, \n radius: {a['virial_radius']}, \n velocity: {a['velocity']}")
        # print(f"IDs: {a['Tree_root_ID'], a['halo_id']}")


    def load_halo(self):
        '''
        This loads the target halo 
        '''
        values = (self.a[self.quantity]).tolist()
        values.sort()
        
        i = np.where(self.a[self.quantity] == (values[-self.how_large]))[0][0]    # gets index 
        my_tree = self.a[i]
        
        return my_tree


    def get_properties(self):
        '''
        Gets target halo and its progenitors' properties
        '''
        my_tree = self.load_halo()
        
        halo_position = np.array(my_tree["position"].to("unitary"))
        halo_radius = my_tree["virial_radius"].to("unitary").value
        prog_redshifts = my_tree["prog", "redshift"]
        prog_snaps = my_tree["prog", "Snap_idx"]
        # get progenitor potisions
        prog_positions = np.array(my_tree["prog", "position"].to("unitary"))
        halo_mass = my_tree["mass"].to('Msun').value
        
        self.width = 2 * halo_radius            # minimum camera width to use in movies etc.

        return halo_position, halo_radius, halo_mass, prog_snaps, prog_positions


    def get_neighbour(self, how_far=5.):
        '''
        Gets neighbouring halo based on criteria
        This can be modified to use a sphere object instead since it is specific to one halo
        '''
        initial_position, _, _, prog_snaps, _ = self.get_properties()
        self.how_far = float(input("How far should I search (in units of the halo's diameter)? "))
        
        # empty list for halo distances
        distances = []
        positions = []

        for halo in a.select_halos("(tree['forest', 'mass'].to('Msun') > 6e14) & (tree['forest', 'mass'].to('Msun') > 690072000000000.0)"):        # exclude the most massive halo
            pair_pos = halo["position"].to("unitary")
            
            # calculate distance to target halo
            distance = np.linalg.norm(np.array(pair_pos) - np.array(initial_position))
            positions.append(np.array(pair_pos))
            distances.append(distance)
        
        if len(distances) == 0:
            print('No haloes found.')
            
        # now we want to choose a halo that is neither too close nor too far from the target
        idx = (np.abs(distances - self.how_far * np.array(self.width))).argmin()
        final_position = positions[idx]
        
        return initial_position, prog_snaps, final_position


    def neighbour_path(self, frames=100):
        '''
        Create path between neighbours
        '''
        self.frames = frames
        
        initial_position, prog_snaps, final_position, _ = self.get_neighbour()
        
        # create arrays for the path
        x_initial, y_initial, z_initial = np.array(initial_position[0]), np.array(initial_position[1]), np.array(initial_position[2])
        x_final, y_final, z_final = np.array(final_position[0]), np.array(final_position[1]), np.array(final_position[2])
        
        x_path, y_path, z_path = np.linspace(x_initial, x_final, self.frames), np.linspace(y_initial, y_final, self.frames), \
        np.linspace(z_initial, z_final, self.frames)
        
        camera_path = np.vstack((x_path, y_path, z_path)).T
        
        return initial_position, prog_snaps, final_position, camera_path


    def make_prog_path(self, pol_degree=10):
        '''
        Makes a smooth camera path for time series. Uses the savgol filter.
        Input: positions to navigate through
        '''
        
        position_points = self.get_properties()[-1]      # can make this an argument
        
        camera_path = np.empty(position_points.shape)
        window_size = 2*int((len(position_points)/4)) + 1            # window size must always be odd

        smoother_pos_x = savgol_filter(position_points[:,0], window_size, pol_degree)
        smoother_pos_y = savgol_filter(position_points[:,1], window_size, pol_degree)
        smoother_pos_z = savgol_filter(position_points[:,2], window_size, pol_degree)
        
        camera_path[:,0], camera_path[:,1], camera_path[:,2] = smoother_pos_x, smoother_pos_y, smoother_pos_z

        return camera_path


    def make_halo_region(self, how_far=10.):
        '''
        Create data object
        '''
        initial_position, initial_radius, _, _, _ = self.get_properties()
        ds = self.a.ytds
        sphere = ds.sphere(initial_position, (initial_radius*how_far, "unitary"))
        center = initial_position

        return sphere


    def get_flybys(self):
        '''
        This finds flybys based on criteria provided above
        '''

        my_tree = self.load_halo()
        
        sphere = self.make_halo_region()
        initial_position, initial_radius, primary_mass, _, prog_positions= self.get_properties()
        primary_radii = (my_tree["prog", "virial_radius"].to("unitary")).value
        primary_mass = my_tree["mass"].value
        
        found = 0           # to keep track of identified flybys
        my_momenta = my_tree["prog", "angular_momentum_magnitude"]
        
        for node in self.a.get_nodes_from_selection(sphere):
            node_redshifts = node["prog", "redshift"]
            node_prog_positions = np.array(node["prog", "position"].to("unitary"))
            node_radii = (node["prog", "virial_radius"].to("unitary")).value
            node_snaps = node["prog", "Snap_idx"]
            node_momenta = node["prog", "angular_momentum_magnitude"]

            distances, secondary_momenta, secondary_x, secondary_y, secondary_z = [], [], [], [], []           # can probably be simplified
            main_radii, main_momenta, radii, redshifts, positions, snaps = [], [], [], [], [], []
            primary_x, primary_y, primary_z = [], [], []            
            

            if (node['mass'].value >= 0.1 * primary_mass and node['mass'].value <= primary_mass):        # do this for haloes with a mass comparable to that of the primary one
                pair_dist = np.linalg.norm(np.array(node['position'].to("unitary")) - np.array(initial_position))
                
                # start timing to not waste too much time looking for flybys
                start = time.time()
                
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
                    secondary_momenta.append(node_momenta[i])
                    main_momenta.append(my_momenta[i])
                    
                    main_radii.append(primary_radii[i])
                    redshifts.append(node_redshifts[i])
                    
                    # for pandas dataframe purposes:
                    
                    secondary_x.append(node_pos[0])
                    secondary_y.append(node_pos[1])
                    secondary_z.append(node_pos[2])               
                    
                    primary_x.append(main_pos[0])
                    primary_y.append(main_pos[1])
                    primary_z.append(main_pos[2])
                    
                    now = time.time()
                    
                # remember time is in seconds
                if now-start >= 3600:
                    print("Time is up...")
                    break
                
                else:
                    # this checks for internal flybys and for minima in distance array
                    if (np.any((np.array(main_radii) + np.array(radii) - np.array(distances)) <= 0) and \
                        np.any(np.r_[True, np.array(distances)[1:] < np.array(distances)[:-1]] & \
                        np.r_[np.array(distances)[:-1] < np.array(distances)[1:], True]==True)):
        
                        print('Suspected flyby found!')

                        df = pd.DataFrame({"x": primary_x, "y": primary_y, "z": primary_z, "x'": secondary_x, "y'": secondary_y, "z'": secondary_z, 
                                        "Distances" : distances, "Snapshot" : snaps, "Flyby Radius": radii, "Primary Radius": main_radii, 
                                        "Redshift" : redshifts, "Primary AngMom": main_momenta, "Secondary AngMom": secondary_momenta})
                        
                        name = f'{node}_{self.quantity}={self.how_large}_dist.csv'
                        print(f"Saving as {name}")
                        
                        df.to_csv(name, index=False)
                        break 
                    continue

                break


    def flyby_posns(self):
        # self.filename = str(input("File name: ")) 
        # read the dataframe
        df = pd.read_csv(self.filename)
        
        fly_radius = np.array(df["Flyby Radius"])
        main_radius = np.array(df["Primary Radius"])         # to be used in width of camera
        redshifts = np.array(df['Redshift'])
        snaps = np.array(df['Snapshot'])
        dists = np.array(df['Distances'])
        
        x_posns = np.array(df["x'"])
        y_posns = np.array(df["y'"])
        z_posns = np.array(df["z'"])
        
        x_main = np.array(df["x"])
        y_main = np.array(df["y"])
        z_main = np.array(df["z"])
        
        
        window_size = 2*int((len(x_posns)/4)) + 1            # window size must always be odd
        pol_degree = 8
        
        camera_path = np.empty((len(x_posns),3))

        # relative positions of flyby accounting for boundary conditions
        rel_x = np.subtract(np.mod(np.add(x_main - x_posns, np.multiply(1, 1/2)), 1), np.multiply(1, 1/2))
        rel_y = np.subtract(np.mod(np.add(y_main - y_posns, np.multiply(1, 1/2)), 1), np.multiply(1, 1/2))
        rel_z = np.subtract(np.mod(np.add(z_main - z_posns, np.multiply(1, 1/2)), 1), np.multiply(1, 1/2))
        
        window_size = window_size = 2*int((len(rel_x)/4)) + 1
        pol_degree = 8
        
        # smooth path of flyby relative to primary halo
        resultant_smooth_x = savgol_filter(rel_x, window_size, pol_degree)
        resultant_smooth_y = savgol_filter(rel_y, window_size, pol_degree)
        resultant_smooth_z = savgol_filter(rel_z, window_size, pol_degree)
        
        # smooth path of the flyby
        smoother_pos_x = savgol_filter(x_posns, window_size, pol_degree)
        smoother_pos_y = savgol_filter(y_posns, window_size, pol_degree)
        smoother_pos_z = savgol_filter(z_posns, window_size, pol_degree)
        
        # smooth path of the main halo
        main_smoother_pos_x = savgol_filter(x_main, window_size, pol_degree)
        main_smoother_pos_y = savgol_filter(y_main, window_size, pol_degree)
        main_smoother_pos_z = savgol_filter(z_main, window_size, pol_degree)
        
        # camera path for rendering purposes
        camera_path[:,0], camera_path[:,1], camera_path[:,2] = smoother_pos_x, smoother_pos_y, smoother_pos_z
        
        # all together:
        plt.plot(redshifts, resultant_smooth_x, label='x-coordinate',linewidth=1, marker='.')
        plt.plot(redshifts, resultant_smooth_y, label='y-coordinate',linewidth=1, marker='.')
        plt.plot(redshifts, resultant_smooth_z, label='z-coordinate',linewidth=1, marker='.')
        
        plt.xlabel('Redshift')
        plt.ylabel('Relative Position')
        plt.xlim(0)
        plt.legend()
        
        plt.title('Flyby position relative to the primary halo')
        plt.savefig(f'{self.filename}_path_boundary.png')
        
        plt.clf()
        fig = plt.figure(figsize=(16, 16))
        ax = fig.add_subplot(111, projection='3d')

        # Data for a three-dimensional line
        zline = resultant_smooth_z
        xline = resultant_smooth_x
        yline = resultant_smooth_y
        # img = ax.scatter3D(xline, yline, zline, c=redshifts, cmap='jet', marker='o', label='smoothed path')
        ax.plot3D(xline, yline, zline, label='smoothed path', linewidth=0.5, c='k')
        
        # Data for three-dimensional scattered points
        zdata = z_posns - z_main
        xdata = x_posns - x_main
        ydata = y_posns - y_main
        
        ax.set_xlabel('X Path')
        ax.set_ylabel('Y Path')
        ax.set_zlabel('Z Path')

        img = ax.scatter3D(rel_x, rel_y, rel_z, c=redshifts, cmap='jet', marker='x', label='actual path')
        
        cbar = fig.colorbar(img, fraction=0.03, pad=0.04)
        cbar.set_label('Redshift')
        
        
        plt.legend(loc=2, prop={'size': 9}, bbox_to_anchor=(0.1,0.9))
        fig.savefig(f'{self.filename}_3Dpath_boundary.png')
    
        # make useful arrays
        res_smooth = np.column_stack((resultant_smooth_x, resultant_smooth_y, resultant_smooth_z))
        fly_actual = np.column_stack((x_posns, y_posns, z_posns))
        rel_boundary = np.column_stack((rel_x, rel_y, rel_z))
        main_actual = np.column_stack((x_main, y_main, z_main))
        actual_smooth = np.column_stack((main_smoother_pos_x, main_smoother_pos_y, main_smoother_pos_z))

        return dists, camera_path, res_smooth, actual_smooth, fly_actual, main_actual, rel_boundary, redshifts, snaps, main_radius, fly_radius


    def flyby_vis(self):
        my_tree = self.load_halo()
        df = pd.read_csv(self.filename)

        initial_position, initial_radius, primary_mass, _, prog_positions = self.get_properties()
        _, _, _, _, fly_actual, main_actual, _, _, snaps, _, fly_radius = self.flyby_posns()

        
        # for annotations in the projection 
        axis = "y"  # str(input("Projection Axis (x, y or z): "))
        slice1 = None
        slice2 = None
        
        if axis == "x":
            slice1 = int(1)
            slice2 = int(2)
        
        elif axis == "y":
            slice1 = int(2)
            slice2 = int(0)
            
        elif axis == "z":
            slice1 = int(0)
            slice2 = int(1)
            
            
        for i, snap in enumerate(snaps):
            # create coordinates in kpc
            # main halo coordinate
            center =  ((self.ds.quan(main_actual[i+3][0], 'unitary').to('kpc'), 
                        self.ds.quan(main_actual[i+3][1], 'unitary').to('kpc'),
                        self.ds.quan(main_actual[i+3][2], 'unitary').to('kpc')))
            
            # flyby coordinate
            fly_coord  = np.array((self.ds.quan(fly_actual[i][0], 'unitary').to('kpc'),
                                   self.ds.quan(fly_actual[i][1], 'unitary').to('kpc'),
                                   self.ds.quan(fly_actual[i][2], 'unitary').to('kpc'))) - np.array(center)

            # make projection plot
            ds_fly = yt.load(f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{snap:03}_covering_grid.h5")
        
            # can Try with another data container, like a sphere or disk.  (didn't work)
            p = yt.ProjectionPlot(
                ds_fly, axis, self.field,  center=center,  width=(10000, 'kpc'), weight_field=self.field#, data_source=region, #, origin='native'   
            )
            p.set_cmap(field=self.field, cmap="cmyt.arbre")
            p.set_zlim(self.field, 1e41, 1e46)
            p.set_axes_unit("kpc")
            p.hide_axes(draw_frame=True)
            buff_size = 4000                    # not sure if this makes a difference
            p.set_buff_size(buff_size)
            p.set_background_color(self.field)
            p.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True)
            p.annotate_scale(corner="upper_right")

            # we can to convert units back to kpc to make axes values meaningful (see documentation): 
            # https://yt-project.org/doc/reference/api/yt.visualization.plot_window.html?highlight=projectionplot#yt.visualization.plot_window.ProjectionPlot)

            # useful annotation
            p.annotate_sphere([fly_coord[slice1], fly_coord[slice2]], coord_system="plot", radius=(fly_radius[i]/2, "unitary"), 
                              circle_args={"color": "red"})

            # Save the image with the keyword.
            p.save(f"ole_{snap}_10000_{self.quantity}={self.how_large}")
            p.clear_annotations()                                           # clear annotations to avoid confusion


    def plot_camera_path(self, final_position, position_points, pol_degree=10):
        
        # visualise how data points were smoothed
        smoother_pos_x, smoother_pos_y, smoother_pos_z = self.make_prog_path()
        
        plt.plot(x, position_points[:,0], label='x_actual',linewidth=1, marker='x')
        plt.plot(x, smoother_pos_x, label='x_smooth',linewidth=1, marker='.')
        plt.xlabel('Time Step Index')
        plt.ylabel('x-coordinate')
        plt.legend()
        plt.savefig(f'smoother_points_x_{self.quantity}_.png')
        plt.clf()
        
        plt.plot(x, position_points[:,1], label='y_actual',linewidth=1, marker='x')
        plt.plot(x, smoother_pos_y, label='y_smooth',linewidth=1, marker='.')
        plt.xlabel('Time Step Index')
        plt.ylabel('y-coordinate')
        plt.legend()
        plt.savefig(f'smoother_points_y_{self.quantity}_.png')
        plt.clf()
        
        plt.plot(x, position_points[:,2], label='z_actual',linewidth=1, marker='x')
        plt.plot(x, smoother_pos_z, label='z_smooth',linewidth=1, marker='.')
        plt.xlabel('Time Step Index')
        plt.ylabel('z-coordinate')
        plt.legend()
        plt.savefig(f'smoother_points_z_{self.quantity}_order={self.how_large}_.png')
        plt.clf()


    def plot_arbor(self):
        print("Plotting arbor...")
        my_tree = self.load_halo()
        
        p = ytree.TreePlot(my_tree, dot_kwargs={'rankdir': 'LR', 'size': '"12,4"'})
        p.save(f'{self.quantity}_order-{self.how_large}_tree.png')


    def plot_flyby_changes(self, with_polynomial=False, pol_deg=8.):
        # a consequence of flybys might be angular momentum change of the main halo so we plot this
        my_tree = self.load_halo()
        prog_redshifts = my_tree["prog", "redshift"]
        prog_momenta = my_tree["prog", "angular_momentum_magnitude"]
        
        # best fit polynomial 
        # spline was giving me a strange error and stackoverflow says it's because of the big difference 
        # in magnitude between y and x axes. this might also not be useful
        if with_polynomial==True:
            z = np.polyfit(prog_redshifts[:10], prog_momenta[:10], pol_deg)
            p = np.poly1d(z)
            xp = np.linspace(np.min(prog_redshifts[:10]), np.max(prog_redshifts[:10]), 1000)
            plt.plot(xp, p(xp),lw=0.5, linestyle=':', c='g')
            
        plt.scatter(prog_redshifts[:10], prog_momenta[:10], marker='+', c='black')
        plt.plot(prog_redshifts[:10], prog_momenta[:10], lw=0.5, linestyle='--', c='k', alpha=0.6)
        plt.suptitle('Halo angular momentum magnitude as a function of redshift')
        plt.xlabel('Redshift')
        plt.ylabel(r'Anglular momentum magnitude ($Mpc~M_{\odot} km/(h^2 s)$)')
        plt.savefig(f'AngMomfor={self.how_large}')
        plt.clf()


    def plot_single_redshift(self):
        my_tree = self.load_halo()
        prog_redshifts = my_tree["prog", "redshift"]
        prog_masses = my_tree["prog", "mass"].to('Msun')
        plt.plot(prog_redshifts, prog_masses)
        plt.ylabel(r'Halo Mass ($M_{\odot}$)')
        plt.xlabel('Redshift')
        plt.yscale("log")
        plt.title('Halo mass as a function of redshift')
        plt.savefig(f'LOGRe_Redshift_VS_Mass_order={self.how_large}')
        plt.clf()


    def plot_many_halos(self):
        values = (self.a[self.quantity]).tolist()
        values.sort()
        redshifts = []
        masses = []
        snapshots = []
        
        for k in tqdm(range(1,self.how_large+1)):
            values = (self.a[self.quantity]).tolist()
            values.sort()
        
            idx = np.where(self.a[self.quantity] == (values[-k]))[0][0]    # gets index 
            one_tree = self.a[idx]
            
            prog_redshifts = np.array(one_tree["prog", "redshift"])
            prog_masses = np.array((one_tree["prog", "mass"].to('Msun')).value)
            prog_snapshots = np.array(one_tree["prog", "Snap_idx"])
            
            masses.extend(prog_masses)
            redshifts.extend(prog_redshifts)
            snapshots.extend(prog_snapshots)
            # plt.plot(prog_redshifts, prog_masses, label=f'Halo mass order from largest = {k}')
            # plt.yscale("log")
                    
        # Calculate the point density
        xy = np.vstack([redshifts, masses])
        z = gaussian_kde(xy)(xy)
        plt.scatter(prog_snapshots, masses, c=z, s=2, cmap='jet')
    
        # plt.hexbin(redshifts, masses, gridsize=(150,150), cmap ='YlOrRd') 
        # plt.hist2d(redshifts, masses, bins=(50, 50), cmap=plt.cm.jet)
        # plt.ylim(min(redshifts), max(redshifts))
        # plt.xlim(min(masses), min(mass.es))
        plt.colorbar()
        plt.yscale("log")
        plt.xlabel('Redshift')
        plt.ylabel(r'Halo Mass ($M_{\odot}$)')
        plt.xlim(0)
        
        plt.title('Halo mass as a function of redshift')
            
        plt.savefig(f'hexbin Many redshifts until {self.how_large}')
        plt.clf()
        
        # same figure but with snapshots instead of redshifts
        # xy = np.vstack([snapshots, masses])
        # z = gaussian_kde(xy)(xy)
        # plt.scatter(redshifts, masses, c=z, s=20, cmap='jet')
        # plt.yscale("log")
        # plt.colorbar()
    
        # # # plt.hexbin(redshifts, masses, gridsize = 50, cmap ='YlOrRd') 
        # # plt.hist2d(redshifts, masses, bins=(50, 50), cmap=plt.cm.jet)
        # # # plt.ylim(min(redshifts), max(redshifts))
        # # # plt.xlim(min(masses), min(mass.es))
        # plt.xlabel('Redshift')
        # plt.ylabel(r'Halo Mass ($M_{\odot}$)')
        # plt.xlim(0)
        
        # plt.title('Halo mass as a function of redshift')
            
        # plt.savefig(f'Many redshifts until {self.how_large}')



class DatasetBasics:
    def __init__(self, snap_number=101):
        self.ds = make_grid_name(snap_number)
        self.ad = self.ds.all_data()                                # this is a region describing the entire box
        self.redshift = self.ds.current_redshift
        self.time = self.ds.current_time.in_units("Gyr")
        
    def plot_quant(self):
        snaps = np.arange(102)
        redshifts = np.empty(len(snaps))
        times = np.empty(len(snaps))
        
        for i, snap in enumerate(snaps):
            now_ds = make_grid_name(snap)
            times[i] = now_ds.current_time.in_units("Gyr")
            redshifts[i] =  now_ds.current_redshift
        
        plt.title('Redshift corresponding to each Snapshot Number in the dataset')
        plt.xlabel('Snapshot Number')
        plt.ylabel('Redshift')
        # plt.xlim(np.min(snaps))
        # plt.ylim(np.min(redshifts))
        plt.scatter(snaps, redshifts, marker='.', c='black', s = 0.6)
        plt.plot(snaps, redshifts, linewidth=0.5, linestyle='--', c='k', alpha=0.6)
        plt.grid()
        plt.savefig('SnapRed')

        plt.clf()

        plt.title('Cosmic Time corresponding to each Snapshot Number in the dataset')
        plt.xlabel('Snapshot Number')
        plt.ylabel('Cosmic Time (Gyr)')
        # plt.xlim(np.min(snaps))
        # plt.ylim(np.min(times))
        plt.scatter(snaps, times, marker='.', c='black', s = 0.6)
        plt.plot(snaps, times, linewidth=0.5, linestyle='--', c='k', alpha=0.6)
        plt.grid()
        plt.savefig('SnapTime')
        plt.clf()
        
        
def test():    
    halo = ExploreHalo()
    # halo.plot_arbor()
    # halo.plot_redshift()
    # halo.basics()
    # halo.make_halo_region()
    halo.get_flybys()
    # halo.flyby_posns()

    # halo.plot_flyby_path()
    # halo.flyby_vis()
    
    # data = DatasetBasics()
    # data.plot_quant()
    # halo.plot_flyby_changes()
    
    # halo.plot_many_halos()
    
test()
            
            
# # class Setting:
#     '''
#     We need input to be path positions
#     for rotarion we do it in here at a centre (render once)
#     for time evolution we need to render snaps too 
#     for moving between points we change focus and need path input
#     '''
#     def __init__(self, zoom_factor=5., width=1., units="unitary", nlayers=12):           ### check if units are right
#         self.field = ("grid", "nbody_mass")
#         self.fnumber = fnumber
#         self.ds = make_grid_name(self.fnumber) #
        
#         self.left_corner = self.ds.domain_left_edge
#         self.right_corner = self.ds.domain_right_edge
#         self.ad = self.ds.all_data()
#         self.mi, self.ma = self.ad.quantities.extrema("nbody_mass")
#         self.nlayers = nlayers
#         self.zoom_factor = zoom_factor

#         self.width = self.ds.quan(width, self.units)
#         self.units = units
#         self.path = None
        
#         self.sc = yt.create_scene(self.ds, field=self.field)
#         self.source = self.sc[0]
#         self.cam = self.sc.add_camera(self.ds)

#         # make it pretty
#         self.scene.grey_opacity = True
#         # for resolution purposes
#         self.source.set_use_ghost_zones(True)
        
#         self.cam.width = self.width * self.zoom_factor


#     def transfer_func (self, s_clip=4., fname=101, cmap="cmyt.arbre"):
#         tfh = TransferFunctionHelper(self.ds)
#         tfh.set_field(self.field)
#         tfh.set_log(True)
#         tfh.set_bounds()
#         tfh.build_transfer_function()
        
#         # Set up transfer function
#         tf = yt.ColorTransferFunction((np.log10(self.mi), np.log10(self.ma)))
#         tf.add_layers(self.nlayers, colormap=cmap)
        
#         self.scene.render()
#         self.scene.save(fname, sigma_clip = s_clip)

        
        
        

# #Class with Functions for making time-evolving animations
# class Animation(object):
#     '''
#     can have list of prog snaps or just one prog snap and depending on its dimension we can loop through both positions
#     and snaps.
#     '''

#     def __init__(self, snapshot, nframes):
#         self.snaps
#         self.starting_snap
#         self.grid_nos = snapshot
#         self.nframes = nframes

#     # simply animates the box (focus is by default at the centre of the box)
#     def cosmos_animation(self, s_clip, cam_position, focus=[0,0,0]):

#         for i in range(len(self.snaps)-self.starting_snap):

#             fname = make_grid_name(frame)
            
#             sc = Setting(fname)

#             sc.set_camera(cam_position, 1, focus)
#             sc.camera.resolution = (1024, 1024)
#             sc.scene.render_save(s_clip, self.grid_no)
    
    
#     # Zoom in by a factor of 2 over 5 frames
#     for _ in cam.iter_zoom(2., 5):
#         print('zoom now')
#         sc.save("camera_movement`_%04i.png" % frame, sigma_clip = 4, render=False)
#         frame += 1

#     # Move to the position [-10.0, 10.0, -10.0] over 5 frames
#     pos = ds.arr([-10.0, 10.0, -10.0], "code_length")
#     for _ in cam.iter_move(pos, 5):
#         print('move now')
#         sc.save("camera_movement2_%04i.png" % frame, sigma_clip = 4, render=False)
#         frame += 1

#     def zoom():
#         pass
        
        
        
#     def rotate(angle, nframes=self.nframes):
#         #angle is input in radians
#         for _ in cam.iter_rotate(angle, nframes):#np.pi, 5):
#             print('rotate now')
#             sc.save("camera_movement3_%04i.png" % frame, sigma_clip = 4, render=False)
#             frame += 1

        
        
        
        
    
    


# class SimplePlots:
#     '''
#     Here we include projections and slices 
#     that can be either a movie or a multiplot or a single plot
#     we can also play around with focus but by default that is the origin if not we need input
#     '''
#     def __init__(self):           ### check if units are right
#         self.field = ("grid", "nbody_mass")
#         self.grid_name = make_grid_name(101)
    
    




