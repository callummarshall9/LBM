# LBM
Implementation of Lattice Boltzmann in 3D using D3Q15

# Lattice Boltzmann Implementation

1. Implementation, initialises using the equilibrium distributionb.

2. Time Step Algorithm
- Macroscopic moments are computed
- Equilibrium distribution is calculated.
- The macroscopic fields density and fluid speed is written to hard disk for visualisation or post processing.
- Colllision relaxation is performed.
- Streaming is done.
- The time step is incremented and repeat the time step algorithm.

Implemeted using The Lattice Boltzmann Method Principles and Practice. ISBN: 978-3-319-44647-9
