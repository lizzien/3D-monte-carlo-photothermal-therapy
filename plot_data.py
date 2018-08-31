# Module to plot data

import numpy as np
import colorcet as cc
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from skimage import io
from mpl_toolkits.axes_grid1 import ImageGrid

def plot_spectra(data, labels, title, xlabel, ylabel,
                 xmin, xmax, ymin, ymax, colours, side_labels):

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    for id, d in enumerate(data):
        ax.plot(d[0], d[1], label=labels[id], color=colours[id])
        # find y position for GNP concentration labels
        iy = int(np.where(np.isclose(d[0], xmax))[0])
        y_pos = 0.95 * d[1][iy]
        # add concentration label
        ax.text(xmax + 5, y_pos, side_labels[id], fontsize=14)

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    fig.tight_layout()
    plt.show()
    
def plot_input_image(bv_file, gnp_file, xmax, zmax, title):
    bv_pic = io.imread(bv_file)
    gnp_pic = io.imread(gnp_file)

    bv = bv_pic[:, :, 0] < 255
    gnp = gnp_pic[:, :, 0] < 255
    mix = bv * gnp
    pic = np.zeros((bv_pic.shape[0], bv_pic.shape[1], 3), np.uint8)  # height, width
    pic.fill(255)

    print(pic.shape)      # gives #rows #cols #rgb
    # Identify GNP and blood vessels
    pic[bv] = [255, 0, 0]
    pic[gnp] = [0, 0, 255]
    pic[mix] = [200, 0, 200]
    # Show tissue structure
    fig, ax = plt.subplots()
    plt.imshow(pic)
    my_xticks = [np.arange(-2, 2, 0.25)]
    x_locs, x_labels = plt.xticks()           # Get locations and labels
    y_locs, y_labels = plt.yticks()
    plt.xticks(np.linspace(0, 200, 9), np.linspace(0, 2 * xmax, 9))
    plt.yticks(np.linspace(0, 200, 9), np.linspace(0, 2 * zmax, 9))
    # set legend colors and labels
    colors=['red','blue', (200/255, 0/255, 200/255)]
    labels=['Blood Vessels', 'GNPs', 'Blood Vessels & GNPs']
    plt.xlabel('x ' + r'$[cm]$')
    plt.ylabel('z ' + r'$[cm]$')
    # create a patch (proxy artist) for every color 
    patches = [mpatches.Patch(color=colors[i], label=labels[i])
               for i in range(len(labels))]
    # put those patched as legend-handles into the legend
    plt.legend(handles=patches, loc=1, fontsize=14)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    plt.show()
    
    return pic
    
# plot absorption map    
def plot_map(data, title, xlabel, ylabel):
    im = plt.imshow(data, origin='lower')
    clb = plt.colorbar(im)
    clb.set_label(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

# plot subplots of absorption maps
def plot_grid_maps(data, title, clb_label, xlabel, ylabel, xmax, zmax,
                   separate, x_min, x_max, y_min, y_max,
                   conc, power, wavelengths, min_val, max_val, file_label):

    n_conc = len(conc)
    n_power = len(power)
    n_lambda = len(wavelengths)
    
    # Replace zero values with 0.1 * the next smallest values
    # Allows a logarithmic colourbar to be used
    #flat_list = [item for m in data for row in m for item in row]
    #flat_list = filter(lambda i: i != 0, flat_list)
    #second_min = min(flat_list)
    #print('Replace zeros with: ' + str(second_min * 0.1))
    #for i,m in enumerate(data):
    #    for j,row in enumerate(m):
    #        for k,item in enumerate(row):
    #            if (item == 0):
    #                data[i][j][k] = second_min * 0.1
    #if (min_val == 0):
    #    min_val = 0.0001#second_min * 0.1
    
    plt.rc(['xtick','ytick'], labelsize=8)     # tick labels size
    
    # Find plots to put x and y labels on
    if (n_lambda % 2 == 0):
        xi = (n_lambda / 2) - 1
    else:
        xi = (n_lambda - 1) / 2
    if (n_power % 2 == 0):
        yi = n_power - (n_power / 2)
    else:
        yi = (n_power - 1) / 2

    map_count = 0
    if (separate):
        # Adapted from: https://stackoverflow.com/questions/13784201/
        #               matplotlib-2-subplots-1-colorbar/13784887
        # Date accessed: 15/08/2018

        for c in range(n_conc):
            fig = plt.figure(figsize=(10, 12))
            grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                             nrows_ncols=(n_power, n_lambda),
                             axes_pad=0.3,
                             share_all=True,
                             cbar_location="right",
                             cbar_mode="single",
                             cbar_size="5%",
                             cbar_pad=0.3)

            # Add data to image grid
            for j,ax in enumerate(grid):
                im = ax.imshow(data[map_count],
                               vmin=min_val, vmax=max_val,
                               origin='lower', #norm=colors.LogNorm(),
                               cmap='viridis',
                               extent=[x_min, x_max, y_min, y_max])
                # Calc. concentration, power, & wavelength indices of map
                p = int(j / n_lambda)
                lmda = j % n_lambda
                ax.set_title(r'$P = {}W, \lambda = {:d}nm$'
                             .format(power[p],
                             int(wavelengths[lmda])), fontsize=14)
                if (lmda == xi):
                    ax.set_xlabel(xlabel)
                if (p == yi):
                    ax.set_ylabel(ylabel)

                ax.set_xticklabels(np.linspace(0, 2 * xmax, 5), fontsize=14)
                ax.set_yticklabels(np.linspace(2 * zmax, 0, 5), fontsize=14)

                map_count = map_count + 1

            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, y_max)
            ax.set_xticks(np.linspace(x_min, x_max, 5))
            ax.set_yticks(np.linspace(y_min, y_max, 5))

            # Colorbar
            #formatter = ticker.LogFormatter(10, labelOnlyBase=False) 
            clb = ax.cax.colorbar(im)
                                  #, ticks=[min_val, 0.001, 0.01, 0.1, 1, 10,
                                  #100, 1000, max_val])
                                  #, format=ticker.FuncFormatter(fmt))
            clb.ax.tick_params(labelsize=14)
            ax.cax.toggle_label(True)
            clb.set_label_text(clb_label, fontsize=20)
            fig.subplots_adjust(top=0.95, bottom=0.05, right=0.95, left=0.05)
            plt.savefig('data/plots/absorbed_energy_den_{}_c{}.png'.format(
                        file_label, c),
                        bbox_inches='tight')
            fig.show()
        plt.show()
    else:
        fig = plt.figure(figsize=(10, 10))
        outer = gridspec.GridSpec(n_conc, 1, wspace=0.2, hspace=0.15)

        for c in range(n_conc):
            inner = gridspec.GridSpecFromSubplotSpec(n_power, n_lambda,
                            subplot_spec=outer[c],
                            wspace=0.1, hspace=0.2)

            for j in range(n_power * n_lambda):
                ax = plt.Subplot(fig, inner[j])
                #im = plt.imshow(Slice.T,origin='lower')
                # Use a perceptually linear colour map,
                # made by Peter Kovesi
                # https://peterkovesi.com/projects/colourmaps/
                # 15/08/2018
                im = ax.imshow(data[map_count],
                               vmin=min_val, vmax=max_val,
                               origin='lower',
                               #norm=colors.LogNorm(),
                               cmap='viridis',
                               extent=[x_min, x_max, y_min, y_max])
                               
                # Calculate concentration, power, & wavelength index of map
                p = int(j / n_lambda)
                lmda = j % n_lambda
                ax.set_title(r'$c = {}, P = {}W, \lambda = {}nm$'
                             .format(conc[c], power[p],
                             wavelengths[lmda]), fontsize=10)

                if (p == n_power -1 and c == n_conc - 1 and lmda == 0):
                    ax.set_xlabel('x ' + r'$[cm]$')
                    ax.set_ylabel('z ' + r'$[cm]$')
                ax.set_xticks(np.linspace(x_min, x_max, 5))
                ax.set_yticks(np.linspace(y_min, y_max, 5))
                ax.set_xticklabels(np.linspace(0, 2 * xmax, 5))
                ax.set_yticklabels(np.linspace(0, 2 * zmax, 5))
                if (j < (n_power - 1) * n_lambda or c != n_conc - 1):
                    ax.set_xticklabels([])
                if (lmda != 0):
                    ax.set_yticklabels([])
                fig.add_subplot(ax)
                map_count = map_count + 1

        cax = fig.add_axes([0.88, 0.05, 0.02, 0.9])
        clb = fig.colorbar(im, cax)#, format=ticker.FuncFormatter(fmt))
        clb.ax.tick_params(labelsize=12)
        clb.set_label(clb_label, fontsize=18)
        #fig.tight_layout()
        plt.subplots_adjust(left=0.05, right=0.88, top=0.95, bottom=0.05)
        plt.savefig('data/energy_map_Tue.png')
        plt.show()

# https://stackoverflow.com/questions/25983218
# /scientific-notation-colorbar-in-matplotlib
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
    
def plot_images():
    filename = ['data/fimage.dat','data/qimage.dat','data/uimage.dat']
    plt.figure()
    for i,file in enumerate(filename):
        f = open(file,'r')
        lines = f.readlines()
        f.close()
        images = []
        for j,line in enumerate(lines):
            data = line.split()
            images.append([])
            for d in data:
                images[j].append(float(d))

        plt.subplot(1,3,i+1)
        plt.imshow(images, cmap='gray')
        plt.colorbar()
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    #cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    #plt.colorbar(cax=cax)
    plt.show()

def format_plot():
    plt.rc(['xtick','ytick'], labelsize=16)     # tick labels size
    plt.rc('axes', labelsize=18)                # axes label size
    plt.rc('axes', titlesize=20)                # figure title size
    plt.rc('legend', fontsize=16)
    #plt.rc('figure', figsize=(9, 6))

# Following code is adapted from:
# https://www.datacamp.com/community/tutorials/matplotlib-3d-volumetric-data
# Written by Juan Nunez-Iglesias
# Date accessed: 24/07/2018
def remove_keymap_conflicts(new_keys_set):
    for prop in plt.rcParams:
        if prop.startswith('keymap.'):
            keys = plt.rcParams[prop]
            remove_list = set(keys) & new_keys_set
            for key in remove_list:
                keys.remove(key)

def multi_slice_viewer(volume):
    remove_keymap_conflicts({'j', 'k'})
    fig, ax = plt.subplots()
    ax.volume = volume
    ax.index = volume.shape[0] // 2
    im = ax.imshow(volume[ax.index],origin='lower')
    fig.colorbar(im)
    fig.canvas.mpl_connect('key_press_event', process_key)

def process_key(event):
    fig = event.canvas.figure
    ax = fig.axes[0]
    if event.key == 'j':
        previous_slice(ax)
    elif event.key == 'k':
        next_slice(ax)
    fig.canvas.draw()

def previous_slice(ax):
    volume = ax.volume
    ax.index = (ax.index - 1) % volume.shape[0]  # wrap around using %
    ax.images[0].set_array(volume[ax.index])

def next_slice(ax):
    volume = ax.volume
    ax.index = (ax.index + 1) % volume.shape[0]
    ax.images[0].set_array(volume[ax.index])
