# spin_tools
A library that contains functions to calculate energies of atomic states in a magnetic field.
TODO: Add antenna simulator.

# Install instructions
Install anaconda (miniforge for M1 Mac). Activate the conda environment in your terminal.

To install the required packages, use the provided conda environment file with
`conda env create --file environment.yaml`.
For M1 Mac, qutip (which is required) has to be compiled and installed from souce.

Install this library by navigating to its root directory and executing
`conda develop .`
If this does not work, try
`python setup.py develop`.

# Usage
Please look in the `tutorial` folder of this repository. It contains jupyter notebook that should explain the basics.