# Module to read in data

import numpy as np

# read data from file
def get_data(file, data_start, params):
    f = open(file)
    lines = f.readlines()[data_start:]

    try:
        for line in lines:
            data = line.split()
            for i in range(len(params)):
                params[i].append(float(data[i]))

    finally:
        f.close()

    return params

def get_array(file, array):
    f = open(file)
    lines = f.readlines()

    try:
        for i,line in enumerate(lines):
          data = line.split()
          for j,d in enumerate(data):
            array[i,j] = float(d)

    finally:
        f.close()

    return array

def get_maps(file):
    f = open(file)
    lines = f.readlines()

    try:
        for i,line in enumerate(lines):
            data = line.split()
            if (i == 0):
                print(data)
                n_maps = int(data[0])
                n_rows = int(data[1])
                n_cols = int(data[2])
                maps = np.array([[[0. for x in range(n_cols)]
                                for y in range(n_rows)]
                                for z in range(n_maps)])
            else:
                # take away one because first line is the dimensions
                m = int((i - 1) / n_rows)
                row = (i - 1) - (m * n_rows)
                for col,d in enumerate(data):
                    maps[m][row][col] = float(d)
    finally:
        f.close()

    return maps
    

