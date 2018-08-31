# Module for writing data to file

import numpy as np

def set_density_grid(density_grid, pic, mua_tissue, mua_b, mua_gnp):
    f = open('data/abs_coeff_structure.dat', 'w')
    print("mu_a tissue: " + str(mua_tissue)
          + " mu_a blood: " + str(mua_b)
          + " mu_a gnp: " + str(mua_gnp))
    try:
        #f.write(str(density_grid.shape[0]) + '\t' + str(density_grid.shape[1]) + '\n')
        for i in range(pic.shape[0]):
            for j in range(pic.shape[1]):
                if np.array_equal(pic[i, j, :],[255, 0, 0]):
                    density_grid[i,j] = mua_b     #0.5
                elif np.array_equal(pic[i, j, :], [0, 0, 255]):
                    density_grid[i,j] = mua_gnp #0.5
                elif np.array_equal(pic[i, j, :], [200, 0, 200]):
                    density_grid[i,j] = 0.5 * mua_gnp + 0.5 * mua_b
                else:
                    density_grid[i,j] = mua_tissue #1.0
            for n in density_grid[i]:
                f.write(formatF(n) + '\t')
            #f.write(str(n) for num in line.split())
            if (i != pic.shape[0] - 1):
                f.write('\n')
    finally:
        f.close()
    
    return density_grid
    
def formatF(value):
    return "{0:.3f}".format(value)

def save_results(file, data):
    f = open(file, 'w')
    n_maps = len(data)
    n_rows = len(data[0])
    n_cols = len(data[0][0])
    try:
        f.write(str(n_maps) + '\t' + str(n_rows) + '\t'
                + str(n_cols) + '\n')
        for i,m in enumerate(data):
            for j,row in enumerate(m):
                for item in row:
                    f.write(formatF(item) + '\t')
                if ((i != n_maps - 1) or (j != n_rows - 1)):
                    f.write('\n')
    finally:
        f.close()
