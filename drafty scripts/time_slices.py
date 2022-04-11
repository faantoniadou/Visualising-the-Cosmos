import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import datetime
from yt.units import kpc

import yt

today = datetime.datetime.now()

fns = [
    "/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_000_covering_grid.h5",
    "/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_033_covering_grid.h5",
    "/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_067_covering_grid.h5",
    "/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_101_covering_grid.h5",
]

fig = plt.figure()

# See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
# These choices of keyword arguments produce a four panel plot with a single
# shared narrow colorbar on the right hand side of the multipanel plot. Axes
# labels are drawn for all plots since we're slicing along different directions
# for each plot.
grid = AxesGrid(
    fig,
    (0.075, 0.075, 0.85, 0.85),
    nrows_ncols=(2, 2),
    axes_pad=0.05,
    label_mode="L",
    share_all=True,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="3%",
    cbar_pad="0%",
)


for i, fn in enumerate(fns):
    # Load the data and create a single plot
    ds = yt.load(fn)  # load data
    p = yt.SlicePlot(ds, "z", ("grid", "nbody_mass"), width=(80000, "kpc"))
    p.set_cmap(field=("grid", "nbody_mass"), cmap="octar")
    # Ensure the colorbar limits match for all plots
    p.set_zlim(("grid", "nbody_mass"), 1e36, 1e46)
    # nicen up the plot by setting the background color to the minimum of the colorbar
    p.set_background_color(("grid", "nbody_mass"))
    # hide the colorbar:
    p.hide_colorbar()

    # hide the axes, while still keeping the background color correct:
    p.hide_axes(draw_frame=True)
    p.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True)
    p.annotate_scale(corner="upper_right")

    # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots[("grid", "nbody_mass")]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

    # Finally, this actually redraws the plot.
    p._setup_plots()
plt.savefig("{:%Y_%m_%d_%H_%M_%S}".format(today)+"_multi.png")
