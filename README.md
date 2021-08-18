![CASCADE](cascadebanner.png)
===
[![DOI](https://zenodo.org/badge/180881021.svg)](https://zenodo.org/badge/latestdoi/180881021) 

# ChemicAl Shift CAlculation with DEep learning 

This repositary contains all codes and data required to reproduce the work: Real-time NMR chemical shift prediction using a graph convolutional network. The production webserver can be found here: http://nova.chem.colostate.edu/cascade/

For all work using CASCADE please cite the primary publication: Guan, Y.; Sowndarya, S. V. S.; Gallegos, L. C.; St. John P. C.; Paton, R. S. *Chem. Sci.* **2021**, [**DOI:** 10.1039/D1SC03343C](https://doi.org/10.1039/D1SC03343C)

# CASCADE: What is it?

**CASCADE** stands for **C**hemic**A**l **S**hift **CA**lculation with **DE**ep learning. It is a stereochemically-aware online calculator for NMR chemical shifts using a graph network approach developed at Colorado State University. Molecular input can be specified as SMILES or through the graphical interface. An automated workflow executes 3D structure embedding and MMFF conformer searching. The full ensemble of optimized conformations are passed to a trained graph neural network to predict the NMR chemical shift (in ppm) for each C atom.

The graph Network considers each atom as a node in a molecular graph that permits communication between atoms - this is through edges connecting atom pairs. Each molecule is represented by a 3D molecular graph, in which each pair of atoms are connected by a distance (i.e., as opposed to solely representing bonds in a 2D molecular graph). The architecture of the network used is shown below.

<p align="center">
  <img src="graph_network.png">
</p>
 
We evaluated the performance of three neural network models. The first model, DFTNN, takes DFT optimized structures as input and was trained against DFT calculated shielding tensors. The MAE for the 500 testing molecules is 1.26 ppm compared to  mPW1PW91/6-311+G(d,p) DFT calculations. The second model, ExpNN-dft, which was trained against experimental chemical shifts using DFT optimized structures as input through transfer learning achieves a MAE of 1.25 ppm for the 500 testing molecules compared to experimental observed chemical shifts. The eventual model, Expnn-ff, which serves as the backend predictor for this web app takes MMFF optimized structures as input and predict experimental chemical shift. The performance of the eventual model on the testing molecules is 1.43 ppm compared to the experimental 13C chemical shifts. The Scatter plot and histogram show the correlation between predicted and experimental chemical shift are shown below.

<p align="center">
  <img src="performance.png">
</p>

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
2. Keras 2.2.5
3. Pandas
4. Rdkit
5. scikit-learn

To train the model on a GPU, we recommend using cuda 9.0.

## Organization of the code
The NMR8K_to_DFT8K folder contains all codes required to run the automatic workflow that optimizes structures and calculates NMR chemical shifts. Details of how to run the workflow can be found in the directory.

The predicting_model folder contains all codes for the neural network. The predicting_model/C contains three models discussed in the paper (DFTNN, ExpNN-dft, and ExpNN-ff). Each directory includes a script of preprocessing data, a script of training the neural network, and a trained model. The predicting_model/C only contains one model for proton prediction, the DFTNN, the organization is same as that of DFTNN for C. In each subdirectory, the preprocess_pubchem.py read data from the data folder and covert it into the form that the neural network can take in. The organization of the data is explains as follows.

## Organization of the data
First, the data folder contains three dataset used in the paper, NMR8K, DFT8K, and Exp5K. The NMR8K folder contains a zip file of the sdf file for about 8K 2D structures sampled from NMRShiftDB. The DFT8K folder contains a zip file of the sdf file for the 8K DFT optimized 3D structures (DFT.sdf.gz), a zip file of the sdf file for the 8K FF optimized 3D structures (FF.sdf.gz), and a zip file of the csv file for the DFT calculated chemical shifts for all C and H. The organization of the Exp5K is the same as DFT8K, but the csv files contains the experimental observed C13 chemical shifts filtered by the DFT calculations. 

The CHESIRE folder contains optimized structures and experimental C13 chemical shifts for 24 molecules of the CHESHIRE probing set, which can be found on the CHESHIRE website (http://cheshirenmr.info/MoleculeSets.htm).

The Fig7_nprod folder contains data for the large molecules with MW > 500 found in the NMRShiftDB, which is used for figure 7 in the paper. 3D structures, experimental C13 chemical shifts, and neural network predicted chemical shifts can be found in the folder. 

## Web app
The ExpNN-ff model which predict experimental chemical shifts using a MMFF 3D structure was developed into a free access web app (CASCADE: http://nova.chem.colostate.edu/cascade/) to fullfill the daily usage of real-time C13 chemical shift predictions that is based on 3D structures and thus can systemitically consider stereochemistry and conformers. 



If you have any questions please open an issue here or contact us at patonlab@colostate.edu.
