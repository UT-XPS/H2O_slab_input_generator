# H2O_slab_input_generator

Generates a nonsense-ice structure as an input for MD simulations of bulk water. Has the option of adding a gap along the z-dimension to generate a slab. It is essential that energy-minimization must be performed on the structure before an MD run is initiated.

Note: for running simulations of bulk water, it is NECESSARY to start with a nonsensical ice structure and then minimize it. If the starting structure is a regular (reasonable) ice, then at room temperature the simulation will end up being stuck in the solid phase.
