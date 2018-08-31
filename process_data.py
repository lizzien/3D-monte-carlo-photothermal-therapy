# Module for processing data

from bisect import bisect_left

# Adapted from: https://stackoverflow.com/questions/12141150/
# from-list-of-integers-get-number-closest-to-a-given-value
# Accessed: 20/08/2018
def takeClosest(myList, myNumber):
    #Assumes myList is sorted. Returns index of closest value to myNumber.
    #If two numbers are equally close, return the index of smallest number.
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return pos
    else:
       return pos - 1

# absorption events
# normalise to get fractional density matrix of incident light absorbed
# in response to one unit of delivered power (or energy)
def get_energy_density(absorption_map, lum, xmax, ymax, zmax,
                       nxg, nyg, nzg, nphotons):
    dV = (xmax / nxg) * (ymax / nyg) * (zmax / nzg)
    for i,row in enumerate(absorption_map):
        for j,item in enumerate(row):
            absorption_map[i][j] = item * lum / (nphotons * dV)
    return absorption_map

def choose_maps(all_maps, n_conc, n_power, n_lambda, powers, lambdas):
    maps = []
    map_count = 0
    for c in range(n_conc):
        for j in range(n_power*n_lambda):
            p = int(j / n_lambda)
            lmda = j % n_lambda
            if (p in powers and lmda in lambdas):
                maps.append(all_maps[map_count])
            map_count = map_count + 1

    max_vals = [max(j) for i in maps for j in i]
    max_val = max(max_vals)
    min_vals = [min(j) for i in maps for j in i]
    min_val = min(min_vals)
    return maps, min_val, max_val
    
