# NMR chemical shift predictions using graph network
This repositary contains all codes and data required to reproduce the work in Real-time NMR chemical shift precitions using graph convolutional network.
## Requirements
### Automatic database2database workflow for QM chemical shift calculations
The code can be found in code/NMR8K_to_DFT8K.

This workflow is essentially a driver for running Gaussian on a high-performance computing (HPC) cluster. The workflow grab 2D structures from a ase database and generate 3D conformers using RDkit. Gaussian 16 is then used to optimize embedded 3D structures and calculate chemical shifts. Gaussian calculations will be sent to the computing node using a queueing software. 

As such, it requires the following:
1. Gaussian16
2. Slurm queueing software
3. Atomic Simulation Environment (ASE)
4. Rdkit and Pybel
5. pandas, numpy
6. tqdm

The workflow should always be run on the head/login node or other node from which you can submit jobs to computing node.

### Graph convolutional neural network
The code and trained model can be found in code/predicting_model

To train the neural network model, it requires the following:
1. Tensorflow 1.12.0 (CPU)/Tensorflow-gpu 1.12.0 (GPU)
2. Keras 2.2.4
3. Pandas
4. Rdkit
5. scikit-learn

To train on a GPU, we recommend cuda 9.0.

## Organization of the code
