# Module to write the parameter values to file

def set_parameters(nphotons, kappa, albedo, g, xmax, ymax, zmax, lum, b):

    # set parameter values
    f = open("params.par")
    newlines = ""
    try:
        for line in f:
            # number of photons
            if line[0:4] == "  np":
                newlines = newlines + "  nphotons={},\n".format(nphotons)
            # kappa (opacity, mass absorption coefficient)
            elif line[0:4] == "  ka":
                newlines = newlines + "  kappa={:.3f}d0,\n".format(kappa)
            elif line[0:4] == "  al":
                newlines = newlines + "  albedo={:.4f}d0,\n".format(albedo)
            elif line[0:4] == "  hg":
                newlines = newlines + "  hgg={:.3f}d0,\n".format(g)
            elif line[0:4] == "  xm":
                newlines = newlines + "  xmax={:.3f}d0,\n".format(xmax)
            elif line[0:4] == "  ym":
                newlines = newlines + "  ymax={:.3f}d0,\n".format(ymax)
            elif line[0:4] == "  zm":
                newlines = newlines + "  zmax={:.3f}d0,\n".format(zmax)
            elif line[0:4] == "  lu":
                newlines = newlines + "  lum={:.3f}d0,\n".format(lum)
            elif line[0:4] == "  b=":
                newlines = newlines + "  b={:.3f}d0/\n".format(b)
            else:
                newlines = newlines + line
    finally:
        f.close()

    f = open("params.par", 'w')
    try:
        f.writelines(newlines)
    finally:
        f.close()