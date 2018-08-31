# Program to read input and plot output for Monte Carlo photon
# transport simulation through a simple tissue model.

# Reads and plots spectra
# Reads tissue image
# Generates 2D array of optical parameters for structure
# Executes external Monte Carlo 3D grid light transport program
# Plots results of simulation as maps of absorbed energy density

import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from scipy.io import FortranFile
from skimage import io
from time import clock

from get_data import get_data, get_array, get_maps
from set_params import set_parameters
from write_data import set_density_grid, save_results
from plot_data import plot_spectra, plot_map, multi_slice_viewer
from plot_data import format_plot, plot_grid_maps
from plot_data import plot_images, plot_input_image
from process_data import takeClosest, get_energy_density, choose_maps

# Set parameter values
nphotons = 10**6      # number of photons
xmax = 1                # extent of grid
ymax = 1
zmax = 1
g = 0                   # scattering asymmetry
nxg = 201               # number of grid cells
nyg = 201
nzg = 201

# get the oxygenated and deoxygenated blood molar extinction coefficient
# e [cm-1/(moles/litre)]
# https://omlc.org/spectra/hemoglobin/summary.html
lambda_b = []
e_hbO2 = []
e_hb = []
get_data('data/h_spectra.txt', 2, [lambda_b, e_hbO2, e_hb])

# Convert extinction coefficient to absorption coefficient, mua [cm-1]:
# mua = ln(10) * ec * molar concentration
#     = ln(10) * ec * x / 64500
# ln(10) ~ 2.303
# molar concentration [mole/litre]: (x g/litre)/(64,500 g Hb/mole)
# x [g/litre]: number of grams per liter ~ 150
# 64,458: gram molecular weight of hemoglobin.
# hgb mass concentration = 150 g/litre
hgbMW = 64458
hgbMC = 150
mua_hbO2 = [ec * 2.303 * hgbMC / hgbMW for ec in e_hbO2]
mua_hb = [ec * 2.303 * hgbMC / hgbMW for ec in e_hb]

# Absorbance of Gold Nanoparticles
lambda_gnp = []
a_gnp = []
get_data('data/GNRs_ctab_160ul_orig.Sample.asc.txt', 88, [lambda_gnp, a_gnp])
mua_gnp = [a * 2.303 for a in a_gnp]
# A = ε l c
# l = 1cm
# c = 2.7 x 10−8 M

# Absorption coefficient of water [cm^-1]
# https://omlc.org/spectra/water/data/hale73.dat
lambda_water = []
mua_water = []
get_data('data/water_hale73.dat', 5, [lambda_water, mua_water])

# Absorption coefficient of fat [cm^-1]
# https://omlc.org/spectra/fat/
lambda_fat = []
mua_fat = []
get_data('data/fat.txt', 6, [lambda_fat, mua_fat])

# Format plots (change font size, etc.)
format_plot()

# Plot absorption coefficient [cm^-1]
data = [[lambda_b, mua_hbO2], [lambda_b, mua_hb],
        [lambda_gnp, mua_gnp], [lambda_water[:60], mua_water[:60]],
        [lambda_fat, mua_fat]]
labels = ['HbO2', 'Hb', 'GNP', 'Water', 'Fat']
plot_spectra(data, labels, 'Absorption Coefficients',
             'Wavelength ' + r'$[nm]$',
             'Absorption Coefficient ' + r'$[cm^{-1}]$',
             500, 1000, 10**(-4), 10**3,
             ['#d62728', '#2ca02c', '#ff7f0e', '#1f77b4', '#9467bd'],
             ['', '', '', '', ''])

# Generic tissue properties
# Values shown are for skin dermis layer
fvb = 0.002         # blood volume fraction
s = 0.39            # oxygen saturation of hemoglobin
fvw = 0.65          # water volume fraction
fvm = 0             # melanosome volume fraction
fvf = 0             # fat volume fraction
mus_reduced = 45.3  # µs.500nm' = reduced scattering coeff. at 500 nm [cm-1]
b_scatt = 1.292     # scattering power

# Choose wavelength range to test
# Select wavelengths from set of wavelengths for water
# because water data set has the least precision
#lambdas = [w for w in lambda_water[12:49]]
lambdas = [w for w in lambda_water[12:49] if w % 50 == 0]
i_blood, i_gnp, i_water, i_fat, reduced_scatter = ([] for i in range(5))
for w in lambdas:
    # Haemaglobin spectra only has even wavelength values
    # so choose closest wavelength
    i_blood.append(takeClosest(lambda_b, w))
    i_gnp.append(next(i for i, _ in enumerate(lambda_gnp) if
                        np.isclose(_, w)))
    i_water.append(next(i for i, _ in enumerate(lambda_water) if
                        np.isclose(_, w)))
    i_fat.append(next(i for i, _ in enumerate(lambda_fat) if
                        np.isclose(_, w)))
    reduced_scatter.append(mus_reduced * (w / 500)**-b_scatt)

# Plot generic tissue, blood, and gold nanoparticle absorption coefficients
mua_blood = [(s * hbO2) + ((1 - s) * hb) for hbO2,hb in zip(mua_hbO2,mua_hb)]
# mua_dry = (fvw * mua_water) + (fvf * mua_fat)
# + (fvm * mua_mel)     to add melanin 
mua_tissue = [(fvw * mua_water[wi]) + (fvw * mua_fat[fi])
              + (fvb * mua_blood[bi])
              for wi,fi,bi in zip(i_water, i_fat, i_blood)]

# Plot input image and get its shape
pic = plot_input_image('data/blood_vessels.png', 'data/gnps.png', xmax, zmax,
                       'Tissue Structure with Gold Nanoparticles')

# Initialise array for holding absorption coefficient values in grid
# The original Fortran90 MC code uses density*opacity to get absorption coefficient
# mu_a = rho * kappa
# Here, we know mu_a, so set kappa = 1 and input density as
# rho = mu_a / kappa
# rho = mu_a
kappa = 1               # opacity
density_grid = np.array([[0. for x in range(pic.shape[1])]
                             for y in range(pic.shape[0])])
print('density slice shape ' + str(density_grid.shape))

# Set parameter values to run simulation for
conc_factor = [1, 5, 10, 15, 20]               # Factor to multiply concentration by
data = []
for c in conc_factor:
    data.append([lambda_gnp, [c * m for m in mua_gnp]])
data.append([lambdas, mua_tissue])
data.append([lambda_b, mua_blood])
labels = ['GNPs', None, None, None, None, 'Tissue', 'Blood']
plot_spectra(data, labels, 'Absorption Coefficients',
             'Wavelength ' + r'$[nm]$',
             'Absorption Coefficient ' + r'$[cm^{-1}]$',
             500, 1000, 10**(-1), 5*10**2,
             ['#ff7f0e', '#ff7f0e', '#ff7f0e', '#ff7f0e', '#ff7f0e', '#1f77b4', '#d62728'],
             [r'$c_0$', r'$5c_0$', r'$10c_0$', r'$15c_0$', r'$20c_0$', '', ''])

power = [1, 5, 10, 15, 20]          # Power [W]
maps_xz, maps_yz, maps_xy = [],[],[]
run_mc = False                      # Run Monte Carlo code, or just plot data
if (run_mc):
    # loop over different parameter values
    for c in conc_factor:
        for lum in power:
            for i in range(len(lambdas)):
                # Set parameter values
                mus = reduced_scatter[i] / (1 - g)
                albedo = mus / (mus + mua_tissue[i])
                print("albedo = " + str(albedo)
                      + " mua = " + str(mua_tissue[i])
                      + " mus = " + str(mus))
                b = 0.25   # radius where beam amplitude is 1/e of axial value
                set_parameters(nphotons, kappa, albedo, g,
                               xmax, ymax, zmax, lum, b)
                # Multiply concentration of GNPs by some factor
                mua_gnp_f = [m * c for m in mua_gnp]
                set_density_grid(density_grid, pic, mua_tissue[i],
                                 mua_blood[i_blood[i]], mua_gnp_f[i_gnp[i]])
    
                # Run Monte Carlo 3D grid code
                print('c: ' + str(c) + ' p: ' + str(lum) + ' lambda: '
                      + str(lambdas[i]))
                start_time = clock()
                subprocess.call("make")
                subprocess.call("./mcgrid")
                print ("MC took", clock() - start_time, "s to run")
    
                # Plot projected image planes
                #plot_images()
    
                # Plot individual absorption map
                #absorb = np.empty([201,201])
                #get_array('data/absorb.dat', absorb)
                #plot_map(absorb, 'Absorption Fraction Map',
                         #'x [cm]', 'y [cm]')

                # Get absorption data
                absorb_xz = np.empty([201,201])
                absorb_yz = np.empty([201,201])
                absorb_xy = np.empty([201,201])
                get_array('data/absorption_count_xz.dat', absorb_xz)
                maps_xz.append(get_energy_density(absorb_xz, lum,
                                                  xmax, ymax, zmax,
                                                  nxg, nyg, nzg, nphotons))
                get_array('data/absorption_count_yz.dat', absorb_yz)
                maps_yz.append(get_energy_density(absorb_yz, lum,
                                                  xmax, ymax, zmax,
                                                  nxg, nyg, nzg, nphotons))
                get_array('data/absorption_count_xy.dat', absorb_xy)
                maps_xy.append(get_energy_density(absorb_xy, lum,
                                                  xmax, ymax, zmax,
                                                  nxg, nyg, nzg, nphotons))
    
    # Save maps to file
    save_results('data/map_data_xz.dat', maps_xz)
    save_results('data/map_data_yz.dat', maps_yz)
    save_results('data/map_data_xy.dat', maps_xy)

# Read map data from file
maps_xz = get_maps('data/map_data_xz.dat')
maps_yz = get_maps('data/map_data_yz.dat')
maps_xy = get_maps('data/map_data_xy.dat')
# Select only the wavelengths and powers you want to plot
maps, min_val, max_val = choose_maps(maps_xz, len(conc_factor),
                                     len(power), len(lambdas),
                                     [3, 4], [6, 7])
plot_grid_maps(maps, 'Absoprtion Maps',
               'Absorbed Energy Density ' +  r'$[W cm^{-3}]$',
               'x ' + r'$[cm]$', 'z ' + r'$[cm]$',
               xmax, zmax, True, 0, 200, 0, 200,
               conc_factor, power[3:5], lambdas[6:8],
               min_val, max_val, 'xz_67')

maps, min_val, max_val = choose_maps(maps_xz, len(conc_factor),
                                     len(power), len(lambdas),
                                     [1, 2, 3], [0, 1, 2, 3])
plot_grid_maps(maps, 'Absoprtion Maps',
               'Absorbed Energy Density ' +  r'$[W cm^{-3}]$',
               'x ' + r'$[cm]$', 'z ' + r'$[cm]$',
               xmax, zmax, True, 0, 200, 0, 200,
               conc_factor, power[1:4], lambdas[0:4],
               min_val, max_val, 'xz_0123')

maps, min_val, max_val = choose_maps(maps_xz, len(conc_factor),
                                     len(power), len(lambdas),
                                     [1, 2, 3], [4, 5, 6, 7])
plot_grid_maps(maps, 'Absoprtion Maps',
               'Absorbed Energy Density ' +  r'$[W cm^{-3}]$',
               'x ' + r'$[cm]$', 'z ' + r'$[cm]$',
               xmax, zmax, True, 0, 200, 0, 200,
               conc_factor, power[1:4], lambdas[4:8],
               min_val, max_val, 'xz_4567')

maps, min_val, max_val = choose_maps(maps_xz, len(conc_factor),
                                     len(power), len(lambdas),
                                     [1, 2, 3], [8, 9, 10])
plot_grid_maps(maps, 'Absoprtion Maps',
               'Absorbed Energy Density ' +  r'$[W cm^{-3}]$',
               'x ' + r'$[cm]$', 'z ' + r'$[cm]$',
               xmax, zmax, True, 0, 200, 0, 200,
               conc_factor, power[1:4], lambdas[8:11],
               min_val, max_val, 'xz_8910')

maps, min_val, max_val = choose_maps(maps_yz, len(conc_factor),
                                     len(power), len(lambdas),
                                     [1, 2, 3], [4, 5, 6, 7])
plot_grid_maps(maps, 'Absoprtion Maps',
               'Absorbed Energy Density ' +  r'$[W cm^{-3}]$',
               'y ' + r'$[cm]$', 'z ' + r'$[cm]$',
               xmax, zmax, True, 0, 200, 0, 200,
               conc_factor, power[1:4], lambdas[4:8],
               min_val, max_val, 'yz_4567')

maps, min_val, max_val = choose_maps(maps_xy, len(conc_factor),
                                     len(power), len(lambdas),
                                     [1, 2, 3], [4, 5, 6, 7])
plot_grid_maps(maps, 'Absoprtion Maps',
               'Absorbed Energy Density ' +  r'$[W cm^{-3}]$',
               'x ' + r'$[cm]$', 'y ' + r'$[cm]$',
               xmax, zmax, True, 0, 200, 0, 200,
               conc_factor, power[1:4], lambdas[4:8],
               min_val, max_val, 'xy_4567')

