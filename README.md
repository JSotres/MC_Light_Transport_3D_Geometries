# A MONTE CARLO CODE FOR LIGHT TRANSPORT SUPPORTING 3D GEOMETRIES

Monte Carlo simulations of photon propagation offer a flexible yet rigorous approach toward photon transport in scattering fields. 
However, scattering is strongly dependent on the geometry of the sample. MCML [1] and most other photon transport codes only work 
in homogeneous or layered geometries. This project aimed at developing a new code that would address photon transport through 3D geometries. 

For this, I used as a starting point a previous model coded in C (MCML, [1]). The main difference, besides coding in C++, was that I 
represented the scattering field as a cuboid divided into layers consisting of cubic voxels. Each layer is infinitely wide, and 
is described by the following parameters: thickness, number and position of the voxels in the layer, dimensions of the voxels, 
and a general index of refraction, anisotropy factor, absorption coefficient and scattering coefficient. 
The refractive indices of the top and bottom ambient mediums are also needed. Then, in order to model specific geometries, 
labels and optical properties of voxels different from those of the rest of the layer can be specified. 
The simulated quantities are photon absorption, specular and diffuse reflectance and transmittance.

[1]  L.H. Wang and S.L. Jacques," Monte Carlo Modeling of Light Transport in Multi-layered Tissues in Standard C", University of Texas / M.D. Anderson Cancer Center, 1992 ( reprinted with corrections 1998). The software package can be downloaded from http://omlc.ogi.edu/software/mc.