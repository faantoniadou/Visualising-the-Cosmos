from matplotlib import rc_context
from matplotlib.animation import FuncAnimation
import datetime
import os
import os.path
import yt
from yt.units import kpc

today = datetime.datetime.now()

ts = yt.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_*")

plot = yt.SlicePlot(ts[0], "z", ("grid", "nbody_mass"))
plot.set_zlim(("grid", "nbody_mass"), 10e40, 11e45)
plot.set_background_color(("grid", "nbody_mass"))
plot.set_cmap(field=("grid", "nbody_mass"), cmap="inferno")
plot.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True)
plot.annotate_scale(corner="upper_right")
fig = plot.plots[("grid", "nbody_mass")].figure

# animate must accept an integer frame number. We use the frame number
# to identify which dataset in the time series we want to load
def animate(i):
    ds = ts[i]
    plot._switch_ds(ds)


animation = FuncAnimation(fig, animate, frames=len(ts))

# Override matplotlib's defaults to get a nicer looking font
with rc_context({"mathtext.fontset": "stix"}):
    animation.save("{:%Y_%m_%d_%H_%M_%S}".format(today)+"animation.mp4")
