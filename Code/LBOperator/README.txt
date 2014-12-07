Anisotropic Laplac-Beltrami operators source code
http://www.di.ens.fr/willow/research/anisotropy/
===========

Here you will find a Matlab implementation of the algorithm described
in the following paper:

Anisotropic Laplace-Beltrami Operators for Shape Analysis
Mathieu Andreux , Emanuele Rodolà, Mathieu Aubry, Daniel Cremers
NORDIA 2014

Note that this implementation has minor differences with the one used to generate the results shown in the paper.

For any questions or feedback regarding the source code please contact Mathieu Andreux mathieu.andreux@polytechnique.org. 



### RUNNING THE CODE:

1. Start by running the init script or add the 'code' folder to your matlab path


2. Compute the spectral decomposition of the Laplace-Beltrami operator by providing a set of vertices and triangular faces. Parameters, in particular the anisotropy strength may have to be tuned depending on the application intended. 
