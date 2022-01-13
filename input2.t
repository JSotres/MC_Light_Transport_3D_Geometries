# input file for transpor.cpp


#x and y positions of the gaussian beam ( if desired )
0.726 0.726
# standard desviation of the gaussian beam (cm)
#0.1

# number of photon packets
10000000
# number of overall layers ( only the ones of the phantom )
2
# Width in x and y (cm), this is the range where the final 
# magnitudes will be represented
1.0 0.5 1.0 0.5
# Maximum and minimum values of the voxels with finite dimensions for each layer
# ixmax ixmin iymax iymin
8 6 8 6
121 119 10 3
# Number of grid elements for the output structure
# ( for x, y, z, radius, angle with the normal of the surface )
150 150 150 
# refraction index of the first layer
1.0
# properties ( zlpos, dx, dy, dz, n, g, mua, mus) of every layer, starting and ending with an air layer
0.000   0.1    0.1 0.006 1.4 0.787 19 480 
0.006   0.006  0.1 0.006 1.4 0.787 2.2 210
0.042   1.0
#This is the second part of the input file.
#Only needed if there are cubes with different 
#optical properties than those of the rest of the layer.
#
# number of layers where there are this different cubes.
1
# number of cubes with differen opt.prop. in each one of the number of layers defined above.
6
# special layer in one line, then the labels of the cube in another line( ix, iy, iz ), and in the
# next one the opt.prop( n, g, mua, mus).
2
120 6 1
1.4 0.995 266 473
120 7 1
1.4 0.995 266 473
120 8 1
1.4 0.995 266 473
120 9 1
1.4 0.995 266 473
120 5 1
1.4 0.995 266 473
120 4 1
1.4 0.995 266 473

