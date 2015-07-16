# Zeropoint-Test
Python program to calculate the Cosmic Variance expected in the 6dFGSv Zeropoint, using mock catalogues of 6dFGSv

### ZeropointTest.py

Performs a test of the cosmic variance in the 6dFGS zeropoint, and the resulting effect on the bulk flow, using the 6dFGSv LambdaCDM mock catalogues (described in Scrimgeour et al. 2015, in prep). 

The program contains 2 tests that can be called by main():

1) linearVelocityTest(): Uses the true linear peculiar velocities from the simulation, with no observational uncertainty. Normalies the true velocities to zero in the Great Circle.

2) observableEtaTest(): Converts the velocities into observable logarithmic quantity x = log_10(Dr/Dz), and normalises these to zero in the Great Circle.

### astrofunctions.py

Contains useful Python functions for astrophysics, e.g. for calculating RA and Dec from Cartesian coordinates.