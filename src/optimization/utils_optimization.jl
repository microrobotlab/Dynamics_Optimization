# UTILITY FUNCTIONS FOR OPTIMIZATION


"""
    packing_fraction_to_Np(packing_fraction, R, L)

Give the number of particles from the packing fraction knowing the particle radius and the dimensions of the space (assumed to have area `L` x `L`).

Packing fraction is used in optimization as a continuous variable instead of the number of particles which is discrete, it is defined as the 
area covered by all particles divided by the total area (cannot use elliptical shape in the current version, the area would not be `L` x `L`).
"""
packing_fraction_to_Np(packing_fraction, R, L) = round(Int, (packing_fraction * L^2) / (π * R^2))