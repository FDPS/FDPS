# $ python py/tipsy.py filename [pmax]
import sys
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt


# import re
import struct
def read_position_tipsy(filename):
    # open the file
    fin = open(filename, "rb")
    size = len(fin.read())
    fin.seek(0)

    # read header (28 = 8 + 4 * 5, additional padding of 4 bytes is introduced)
    little_endian = True
    time, nbodies, ndim, nsph, ndark, nstar = struct.unpack("<diiiii", fin.read(28))
    if (ndim < 1 or ndim > 3):
        little_endian = False
        fin.seek(0)
        time, nbodies, ndim, nsph, ndark, nstar = struct.unpack(">diiiii", fin.read(28))
    # remove padding of 4 bytes
    fin.read(4)

    # this function reads dark matter particles only
    # skip gas particles (12 floats = 48 bytes)
    fin.read(48 * nsph)

    # arrays for dark matter particles
    px = [0] * ndark
    py = [0] * ndark
    pz = [0] * ndark

    # read dark matter particles (9 floats = 36 bytes)
    if ndark > 0:
        for ii in range(ndark):
            if little_endian:
                mass, x, y, z, vx, vy, vz, eps, phi = struct.unpack("<fffffffff", fin.read(36))
            else:
                mass, x, y, z, vx, vy, vz, eps, phi = struct.unpack(">fffffffff", fin.read(36))

            px[ii] = x
            py[ii] = y
            pz[ii] = z

    return (px, py, pz)


def locate_panels(ax, nx, ny, share_xaxis, share_yaxis):
    margin = 0.12
    if (share_xaxis == False) or (share_yaxis == False):
        margin = 0.15

    xmin, xmax = margin, 1.0 - margin
    ymin, ymax = margin, 1.0 - margin
    xbin = (xmax - xmin) / nx
    ybin = (ymax - ymin) / ny
    xmargin, ymargin = 0, 0

    if share_yaxis == False:
        xmin = 0.0
        xbin = 1.0 / nx
        xmargin = xbin * margin

    if share_xaxis == False:
        ymin = 0.0
        ybin = 1.0 / ny
        ymargin = ybin * margin

    for ii in range(nx):
        xl = xmin + ii * xbin + xmargin

        for jj in range(ny):
            yl = ymin + jj * ybin + ymargin
            kk = ii * ny + jj
            ax[kk] = fig.add_axes((xl, yl, xbin - 2 * xmargin, ybin - 2 * ymargin))

            if share_xaxis == True:
                ax[kk].tick_params(labelbottom = "off")
                if jj == 0:
                    ax[kk].tick_params(labelbottom = "on")

            if share_yaxis == True:
                ax[kk].tick_params(labelleft = "off")
                if ii == 0:
                    ax[kk].tick_params(labelleft = "on")


# obtain input argument(s)
argv = sys.argv
argc = len(argv)
filename = argv[1]
set_range = True
if argc == 3:
    set_range = False
    pmax = float(argv[2])

# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
plt.rcParams['font.size'] = 14

# set number of panels
nxpanel, nypanel = 2, 2
ax = [0] * nxpanel * nypanel

# set figure size and its aspect ratio
Lx = 6
if nxpanel > 2 * nypanel:
    Lx *= 2
Ly = (Lx / nxpanel) * nypanel
fig = plt.figure(figsize = (Lx, Ly))

# set location of panels
locate_panels(ax, nxpanel, nypanel, True, True)

# get number of components and particles
input = "doc/" + filename + ".summary.txt"
fin = open(input, "r")
lines = fin.readlines()
fin.close()
idx = lines[1].find("\t")
kind = int(lines[1][0:idx])# number of components

# get partition
num  = [0] * kind
head = [0] * kind
head[0] = 0
for ii in range(kind):
    idx = lines[2 + ii].find("\n")
    num[ii] = int(lines[2 + ii][0:idx])
    if ii != 0:
        head[ii] = head[ii - 1] + num[ii - 1]

# read particle data
px, py, pz = read_position_tipsy("dat/" + filename + ".tipsy")

# set plot range
if set_range == True:
    pmax = max(np.abs(px))
    ymax = max(np.abs(py))
    zmax = max(np.abs(pz))
    if pmax < ymax:
        pmax = ymax
    if pmax < zmax:
        pmax = zmax

# sparse sampling if necessary
skip = 1
nmax = 1048576
if len(px) > nmax:
    skip = int(np.ceil(len(px) / nmax))

# set colors of dots
col = [0] * kind
for ii in range(kind):
    if ii % 6 == 0:
        col[ii] = "black"
    if ii % 6 == 1:
        col[ii] = "red"
    if ii % 6 == 2:
        col[ii] = "blue"
    if ii % 6 == 3:
        col[ii] = "magenta"
    if ii % 6 == 4:
        col[ii] = "green"
    if ii % 6 == 5:
        col[ii] = "brown"

# plot particle distribution
for ii in range(nxpanel):
    for jj in range(nypanel):
        idx = ii * nypanel + jj

        if (idx != (nxpanel * nypanel - 1)):
            # ii = 0, jj = 0: xy-plot
            # ii = 0, jj = 1: xz-plot
            # ii = 1, jj = 0: zy-plot
            # ii = 1, jj = 1: blank
            if ii == 0:
                xx = px
            if ii == 1:
                xx = pz
                ax[idx].tick_params(labelleft = False)
            if jj == 0:
                yy = py
            if jj == 1:
                yy = pz
                ax[idx].tick_params(labelbottom = False)

            # plot the data
            for kk in range(kind):
                ax[idx].plot(xx[head[kk]:(head[kk]+num[kk]):skip], yy[head[kk]:(head[kk]+num[kk]):skip], ",", color = col[kk])

            # set plot range
            ax[idx].set_xlim([-pmax, pmax])
            ax[idx].set_ylim([-pmax, pmax])

            # ax[idx].grid()
            ax[idx].tick_params(axis = "both", direction = "in", color = "black", bottom = "on", top = "on", left = "on", right = "on")

            # set label
            if (ii == 0) and (jj == 0):
                ax[idx].set_xlabel(r"$x$")
                ax[idx].set_ylabel(r"$y$")
            if (ii == 0) and (jj == 1):
                ax[idx].set_ylabel(r"$z$")
            if (ii == 1) and (jj == 0):
                ax[idx].set_xlabel(r"$z$")

        else:
            # remove box at the upper right corner
            ax[idx].spines["right"].set_color("none")
            ax[idx].spines["top"].set_color("none")
            ax[idx].tick_params(labelbottom = False, labelleft = False, bottom = False, left = False, right = False, top = False)


# output the figure
plt.savefig("dot.png", format = "png", dpi = 300, bbox_inches = "tight")
