# input file for transpor.cpp

#x and y positions of the gaussian beam, or of the infinite narrow beam
1.5 1.5
# standard desviation of the gaussian beam (cm)
#0.01

# number of photon packets
1000000
# number of overall layers ( only the ones of the phantom )
3
# Boundaries in x and y (cm)
# xmax xmin ymax ymin
3 0 3 0
# Maximum and minimum values of the voxels with finite dimensions for each layer
# ixmax ixmin iymax iymin
2 0 2 0
2 0 2 0
2 0 2 0
# Number of grid elements for the output structure
# ( for x, y, z, radius, angle with the normal of the surface )
100 100 100 
# refraction index of the first layer
1.0
# properties ( zlpos, dx, dy, dz, n, g, mua, mus) of every layer, starting and ending with an air layer
0.00  1    1 0.1 1.37 0.9 1 100 
0.1   1    1 0.1 1.37 0.0 1 10 
0.2   1    1 0.2 1.37 0.7 2 10
0.4   1    1   1 1    1   1
#This is the second part of the input file.
#Only needed if there are cubes with different 
#optical properties than those of the rest of the layer.
#
# number of layers where there are this different cubes.
1
# number of cubes with differen opt.prop. in each one of the number of layers defined above.
1
# special layer in one line, then the labels of the cube in another line( ix, iy, iz ), and in the
# next one the opt.prop( n, g, mua, mus).
2
2 2 1
1.40 0.995 266.0 473.0
