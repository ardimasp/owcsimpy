

# # JS Animation import is available at http://github.com/jakevdp/JSAnimation
# from JSAnimation.IPython_display import display_animation
# from matplotlib import animation

# import matplotlib.pyplot as plt
# import numpy as np

# # Set up the axes, making sure the axis ratio is equal
# fig = plt.figure(figsize=(6.5, 2.5))
# ax = fig.add_axes([0, 0, 1, 1], xlim=(-0.02, 13.02), ylim=(-0.02, 5.02),
#                   xticks=range(14), yticks=range(6), aspect='equal', frameon=False)
# ax.grid(True)

# # Define the shapes of the polygons
# P1 = np.array([[0, 0], [5, 0], [5, 2], [0, 0]])
# P2 = np.array([[0, 0], [8, 0], [8, 3], [0, 0]])
# P3 = np.array([[0, 0], [5, 0], [5, 1], [3, 1], [3, 2], [0, 2], [0, 0]])
# P4 = np.array([[0, 1], [3, 1], [3, 0], [5, 0], [5, 2], [0, 2], [0, 1]])

# # Draw the empty polygons for the animation
# kwds = dict(ec='k', alpha=0.5)
# patches = [ax.add_patch(plt.Polygon(0 * P1, fc='g', **kwds)),
#            ax.add_patch(plt.Polygon(0 * P2, fc='b', **kwds)),
#            ax.add_patch(plt.Polygon(0 * P3, fc='y', **kwds)),
#            ax.add_patch(plt.Polygon(0 * P4, fc='r', **kwds))]

# # This function moves the polygons as a function of the frame i
# Nframes = 30
# def animate(nframe):
#     f = nframe / (Nframes - 1.0)
#     patches[0].set_xy(P1 + (8 - 8 * f, 3 - 3 * f + 0.5 * np.sin(f * np.pi)))
#     patches[1].set_xy(P2 + (5 * f, 2 * f - 0.5 * np.sin(f * np.pi)))
#     patches[2].set_xy(P3 + (8 - 3 * f, 0))
#     patches[3].set_xy(P4 + (8, 1 - f))
#     return patches
    
# anim = animation.FuncAnimation(fig, animate, frames=Nframes, interval=50)
# display_animation(anim, default_mode='once')

"""
=========================
Simple animation examples
=========================

Two animations where the first is a random walk plot and
the second is an image animation.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def update_line(num, data, line):
    line.set_data(data[..., :num])
    return line,

###############################################################################

fig1 = plt.figure()

# Fixing random state for reproducibility
np.random.seed(19680801)

data = np.random.rand(2, 25)
l, = plt.plot([], [], 'r-')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel('x')
plt.title('test')
line_ani = animation.FuncAnimation(fig1, update_line, 25, fargs=(data, l),
                                   interval=50, blit=True)

# To save the animation, use the command: line_ani.save('lines.mp4')

###############################################################################

fig2 = plt.figure()

x = np.arange(-9, 10)
y = np.arange(-9, 10).reshape(-1, 1)
base = np.hypot(x, y)
ims = []
for add in np.arange(15):
    ims.append((plt.pcolor(x, y, base + add, norm=plt.Normalize(0, 30)),))

im_ani = animation.ArtistAnimation(fig2, ims, interval=50, repeat_delay=3000,
                                   blit=True)
# To save this second animation with some metadata, use the following command:
# im_ani.save('im.mp4', metadata={'artist':'Guido'})

plt.show()
