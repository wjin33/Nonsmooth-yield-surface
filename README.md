# Nonsmooth-yield-surface
This project includes the implementation of the discrete damage model with multiple non-smooth yield surfaces and its application for engineering problems based on our paper: 

Jin, Wencheng, and Chloé Arson. "Micromechanics based discrete damage model with multiple non-smooth yield surfaces: theoretical formulation, numerical implementation and engineering applications." International Journal of Damage Mechanics 27.5 (2018): 611-639.

Please kindly cite above paper if you used any of the functions or algorithms listed in this Github repository, thank you


Matlab code:

DDM_CCP.m -> Matlab implementation of micromechanics based discrete damage model at Gauss Point, Closest Point Projection (return mapping) is used for iteration;


UMAT

UMAT_DDM_3D_CPP.for/UMAT_DDM_2D_CPP.for -> Abaqus UMAT implementation of discrete damage model using Closest Point Projection algorithm for 3 dimentiona/plane strain cases

Abaqus Input

input files used for the above mentioned paper