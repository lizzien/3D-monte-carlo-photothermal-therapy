# 3D-monte-carlo-photothermal-therapy
Simulation of photothermal therapy with gold nanoparticles using Python and a 3D Monte Carlo radiation transport code in FORTRAN 77.

This directory contains the source codes for the 3D Cartesian grid code illuminated by a Gaussian beam.  The current setup models light transport in a tissue simple sample.  The output is maps of the absorbed energy density.

The Python files are:

        get_data.py
        plot_data.py
        process_data.py
        set_params.py
        write_data.py
        
    The program is run by running the Python file:

        mc3d.py

The FORTRAN files are:

    struc.f

    Code from Kenneth Wood, accessed 04.07.2018, source:
    http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html
    
        Most files were made double precision

        Adapted:

            density.f   -   read absorption coefficient from grid
            gridset.f   -   removed AU units
            iarray.f    -   added additional arrays to set up
            mcpolar.f   -   removed forced first scattering and peeling off
                            added weighted absorption and scattering
                            removed image plane
            stokes.f    -   uncommented isotropic scattering
            sourceph.f  -   Gaussian beam source
            sources.f   -   Gaussian beam centre at centre of top surface
            tauint2.f   -   added flux count out of top surface

        As original (except additional comments):

            peeloff.f
            ran2.f
            scatt1.f
            scattp.f
            taufind1.f
            tauint.f

The following code was from Kenneth Wood, accessed 04.07.2018, source:
http://www-star.st-and.ac.uk/~kw25/research/montecarlo/points/points.html

Include files are:

    Adapted:

        grid.txt        -   added additional arrays
        sources.txt     -   number of sources decreased to one

    As original:

        images.txt
        photon.txt

Input parameters are in:

    Adapted:

        params.par      -   added parameters

The file that compiles the FORTRAN 77 code and creates the executable file 'mcgrid' is:

    Adapted:

        Makefile
