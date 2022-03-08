import numpy as np
import matplotlib.pyplot as plt
import sys

import matplotlib as mpl  
mpl.rc('font',family='Times New Roman')

font = {'family': 'Times New Roman',
    'color':  'k',
    'weight': 'normal',
    'size': 17,
    }
legend_tick_size = 17
font2 = {'family': 'Times New Roman',
    'weight': 'normal',
    'size': 17,
    }

from mpl_toolkits import axes_grid1
def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs):
    """Add a vertical color bar to an image plot."""
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)
    return im.axes.figure.colorbar(im, cax=cax, **kwargs)


xpixNum, ypixNum = 60, 60
resolution = 1.5/60

# image plane pixel number and resolution
M, N = 100, 100
resol_imgplane = 4/100.0 #xyregion/M
xyregion = M*resol_imgplane
img_center = (-0.17, 0.2)
ext = [-resolution*(xpixNum//2) +img_center[0], resolution*(xpixNum//2)+img_center[0], -resolution*(ypixNum//2)+img_center[1], resolution*(ypixNum//2)+img_center[1]]
extimg = [-xyregion/2+img_center[0], xyregion/2+img_center[0], -xyregion/2+img_center[1], xyregion/2+img_center[1]]

def load_petsc_vec(filename):
    with open(filename, 'r') as f:
        res = f.readlines()
    return np.array(res[2:]).astype(float)
    # return np.array(res[3:-1]).astype(float)


src1d = load_petsc_vec('./src1d.dat')
img1d = load_petsc_vec('./img1d.dat')
img1d_scat = load_petsc_vec('./img1d_scat.dat')

src = src1d.reshape(xpixNum, -1).T
img2 = img1d.reshape(M, -1).T
img = img1d_scat.reshape(M, -1).T

npmax = np.max(img2)


plt.figure(figsize = (16, 8))
plt.subplots_adjust(wspace=0.2, hspace=0, left = 0.1, right = 0.95, top = 0.95, bottom = 0.05)
ax2 = plt.subplot(132)
im = plt.imshow(img2, extent = extimg, cmap = "jet", vmax = npmax)
plt.title("Lensed image w/o scattering", fontdict = font)
# plt.xlabel("Arcsec")
cbar = add_colorbar(im)
cbar.ax.tick_params(labelsize=legend_tick_size)

ax3 = plt.subplot(133, sharex = ax2, sharey = ax2)
im = plt.imshow(img, extent = extimg, cmap = "jet", vmax = npmax)
plt.title("Lensed image w/ scattering", fontdict = font)


cbar = add_colorbar(im)
cbar.ax.tick_params(labelsize=legend_tick_size)

ax1 = plt.subplot(131)
im = plt.imshow(src, extent = ext, cmap = "jet")
plt.title("Source", fontdict = font)

cbar = add_colorbar(im)
cbar.ax.tick_params(labelsize=legend_tick_size)

ax1.tick_params(axis='both', labelsize = legend_tick_size, direction = "in")
ax2.tick_params(axis='both', labelsize = legend_tick_size, direction = "in")
ax3.tick_params(axis='both', labelsize = legend_tick_size, direction = "in")

ax1.set_ylabel(r'$y$ (arcsec)',fontdict=font, color='k')
ax1.set_xlabel(r'$x$ (arcsec)',fontdict=font, color='k')
ax2.set_xlabel(r'$x$ (arcsec)',fontdict=font, color='k')
ax3.set_xlabel(r'$x$ (arcsec)',fontdict=font, color='k')

plt.show()