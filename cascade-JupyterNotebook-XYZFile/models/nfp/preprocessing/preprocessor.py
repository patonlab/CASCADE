import logging, sys

import numpy as np
from tqdm import tqdm
from scipy.linalg import eigh

from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolToSmiles, AddHs

from . import features
from .features import Tokenizer
import time


class SmilesPreprocessor(object):
    """ Given a list of SMILES strings, encode these molecules as atom and
    connectivity feature matricies.

    Example:
    >>> preprocessor = SmilesPreprocessor(explicit_hs=False)
    >>> inputs = preprocessor.fit(data.smiles)
    """

    def __init__(self, explicit_hs=True, atom_features=None, bond_features=None):
        """

        explicit_hs : bool
            whether to tell RDkit to add H's to a molecule.
        atom_features : function
            A function applied to an rdkit.Atom that returns some
            representation (i.e., string, integer) for the Tokenizer class.
        bond_features : function
            A function applied to an rdkit Bond to return some description.

        """

        self.atom_tokenizer = Tokenizer()
        self.bond_tokenizer = Tokenizer()
        self.explicit_hs = explicit_hs

        if atom_features is None:
            atom_features = features.atom_features

        if bond_features is None:
            bond_features = features.bond_features_v1

        self.atom_features = atom_features
        self.bond_features = bond_features


    def fit(self, smiles_iterator):
        """ Fit an iterator of SMILES strings, creating new atom and bond
        tokens for unseen molecules. Returns a dictionary with 'atom' and
        'connectivity' entries """
        return list(self.preprocess(smiles_iterator, train=True))


    def predict(self, smiles_iterator):
        """ Uses previously determined atom and bond tokens to convert a SMILES
        iterator into 'atom' and 'connectivity' matrices. Ensures that atom and
        bond classes commute with previously determined results. """
        return list(self.preprocess(smiles_iterator, train=False))


    def preprocess(self, smiles_iterator, train=True):

        self.atom_tokenizer.train = train
        self.bond_tokenizer.train = train

        for smiles in tqdm(smiles_iterator):
            yield self.construct_feature_matrices(smiles)


    @property
    def atom_classes(self):
        """ The number of atom types found (includes the 0 null-atom type) """
        return self.atom_tokenizer.num_classes + 1


    @property
    def bond_classes(self):
        """ The number of bond types found (includes the 0 null-bond type) """
        return self.bond_tokenizer.num_classes + 1


    def construct_feature_matrices(self, smiles):
        """ construct a molecule from the given smiles string and return atom
        and bond classes.

        Returns
        dict with entries
        'n_atom' : number of atoms in the molecule
        'n_bond' : number of bonds in the molecule
        'atom' : (n_atom,) length list of atom classes
        'bond' : (n_bond,) list of bond classes
        'connectivity' : (n_bond, 2) array of source atom, target atom pairs.

        """

        mol = MolFromSmiles(smiles)
        if self.explicit_hs:
            mol = AddHs(mol)

        n_atom = len(mol.GetAtoms())
        n_bond = 2 * len(mol.GetBonds())

        # If its an isolated atom, add a self-link
        if n_bond == 0:
            n_bond = 1

        atom_feature_matrix = np.zeros(n_atom, dtype='int')
        bond_feature_matrix = np.zeros(n_bond, dtype='int')
        connectivity = np.zeros((n_bond, 2), dtype='int')

        bond_index = 0

        atom_seq = mol.GetAtoms()
        atoms = [atom_seq[i] for i in range(n_atom)]

        for n, atom in enumerate(atoms):

            # Atom Classes
            atom_feature_matrix[n] = self.atom_tokenizer(
                self.atom_features(atom))

            start_index = atom.GetIdx()

            for bond in atom.GetBonds():
                # Is the bond pointing at the target atom
                rev = bond.GetBeginAtomIdx() != start_index

                # Bond Classes
                bond_feature_matrix[n] = self.bond_tokenizer(
                    self.bond_features(bond, flipped=rev))

                # Connectivity
                if not rev:  # Original direction
                    connectivity[bond_index, 0] = bond.GetBeginAtomIdx()
                    connectivity[bond_index, 1] = bond.GetEndAtomIdx()

                else:  # Reversed
                    connectivity[bond_index, 0] = bond.GetEndAtomIdx()
                    connectivity[bond_index, 1] = bond.GetBeginAtomIdx()

                bond_index += 1


        return {
            'n_atom': n_atom,
            'n_bond': n_bond,
            'atom': atom_feature_matrix,
            'bond': bond_feature_matrix,
            'connectivity': connectivity,
        }


class ConnectivityAPreprocessor(object):
    """ Given a list of SMILES strings, encode these molecules as atom and
    connectivity feature matricies.

    Example:
    >>> preprocessor = SmilesPreprocessor(explicit_hs=False)
    >>> inputs = preprocessor.fit(data.smiles)
    """

    def __init__(self, explicit_hs=True, atom_features=None, bond_features=None):
        """

        explicit_hs : bool
            whether to tell RDkit to add H's to a molecule.
        atom_features : function
            A function applied to an rdkit.Atom that returns some
            representation (i.e., string, integer) for the Tokenizer class.
        bond_features : function
            A function applied to an rdkit Bond to return some description.

        """

        self.atom_tokenizer = Tokenizer()
        self.bond_tokenizer = Tokenizer()
        self.explicit_hs = explicit_hs

        if atom_features is None:
            atom_features = features.atom_features_v1

        if bond_features is None:
            bond_features = features.bond_features_v1

        self.atom_features = atom_features
        self.bond_features = bond_features


    def fit(self, smiles_iterator):
        """ Fit an iterator of SMILES strings, creating new atom and bond
        tokens for unseen molecules. Returns a dictionary with 'atom' and
        'connectivity' entries """
        return list(self.preprocess(smiles_iterator, train=True))


    def predict(self, smiles_iterator):
        """ Uses previously determined atom and bond tokens to convert a SMILES
        iterator into 'atom' and 'connectivity' matrices. Ensures that atom and
        bond classes commute with previously determined results. """
        return list(self.preprocess(smiles_iterator, train=False))


    def preprocess(self, smiles_iterator, train=True):

        self.atom_tokenizer.train = train
        self.bond_tokenizer.train = train

        for smiles in tqdm(smiles_iterator):
            yield self.construct_feature_matrices(smiles)


    @property
    def atom_classes(self):
        """ The number of atom types found (includes the 0 null-atom type) """
        return self.atom_tokenizer.num_classes + 1


    @property
    def bond_classes(self):
        """ The number of bond types found (includes the 0 null-bond type) """
        return self.bond_tokenizer.num_classes + 1


    def construct_feature_matrices(self, smiles):
        """ construct a molecule from the given smiles string and return atom
        and bond classes.

        Returns
        dict with entries
        'n_atom' : number of atoms in the molecule
        'n_bond' : number of bonds in the molecule
        'atom' : (n_atom,) length list of atom classes
        'bond' : (n_bond,) list of bond classes
        'connectivity' : (n_bond, 2) array of source atom, target atom pairs.

        """

        mol = MolFromSmiles(smiles)
        if self.explicit_hs:
            mol = AddHs(mol)

        n_atom = len(mol.GetAtoms())
        n_bond = 2 * len(mol.GetBonds())

        # If its an isolated atom, add a self-link
        if n_bond == 0:
            n_bond = 1

        atom_feature_matrix = np.zeros(n_atom, dtype='int')
        bond_feature_matrix = np.zeros(n_bond, dtype='int')
        connectivity = np.zeros((n_bond, 2), dtype='int')

        bond_index = 0

        atom_seq = mol.GetAtoms()
        atoms = [atom_seq[i] for i in range(n_atom)]

        for n, atom in enumerate(atoms):

            # Atom Classes
            atom_feature_matrix[n] = self.atom_tokenizer(
                self.atom_features(atom))

            start_index = atom.GetIdx()

            for bond in atom.GetBonds():
                # Is the bond pointing at the target atom
                rev = bond.GetBeginAtomIdx() != start_index

                # Bond Classes
                bond_feature_matrix[n] = self.bond_tokenizer(
                    self.bond_features(bond, flipped=rev))

                # Connectivity
                if not rev:  # Original direction
                    connectivity[bond_index, 0] = bond.GetBeginAtomIdx()
                    connectivity[bond_index, 1] = bond.GetEndAtomIdx()

                else:  # Reversed
                    connectivity[bond_index, 0] = bond.GetEndAtomIdx()
                    connectivity[bond_index, 1] = bond.GetBeginAtomIdx()

                bond_index += 1

        return {
            'n_atom': n_atom,
            'n_bond': n_bond,
            'atom': atom_feature_matrix,
            'bond': bond_feature_matrix,
            'connectivity': connectivity,
        }


class MolPreprocessor(SmilesPreprocessor):
    """ I should refactor this into a base class and separate
    SmilesPreprocessor classes. But the idea is that we only need to redefine
    the `construct_feature_matrices` method to have a working preprocessor that
    handles 3D structures.

    We'll pass an iterator of mol objects instead of SMILES strings this time,
    though.

    """

    def __init__(self, n_neighbors, cutoff, **kwargs):
        """ A preprocessor class that also returns distances between
        neighboring atoms. Adds edges for non-bonded atoms to include a maximum
        of n_neighbors around each atom """

        self.n_neighbors = n_neighbors
        self.cutoff = cutoff
        super(MolPreprocessor, self).__init__(**kwargs)


    def construct_feature_matrices(self, mol):
        """ Given an rdkit mol, return atom feature matrices, bond feature
        matrices, and connectivity matrices.

        Returns
        dict with entries
        'n_atom' : number of atoms in the molecule
        'n_bond' : number of edges (likely n_atom * n_neighbors)
        'atom' : (n_atom,) length list of atom classes
        'bond' : (n_bond,) list of bond classes. 0 for no bond
        'distance' : (n_bond,) list of bond distances
        'connectivity' : (n_bond, 2) array of source atom, target atom pairs.

        """

        n_atom = len(mol.GetAtoms())

        # n_bond is actually the number of atom-atom pairs, so this is defined
        # by the number of neighbors for each atom.
        #if there is cutoff,
        distance_matrix = Chem.Get3DDistanceMatrix(mol)

        if self.n_neighbors <= (n_atom - 1):
            n_bond = self.n_neighbors * n_atom
        else:
            # If there are fewer atoms than n_neighbors, all atoms will be
            # connected
            n_bond = distance_matrix[(distance_matrix < self.cutoff) & (distance_matrix != 0)].size

        if n_bond == 0: n_bond = 1

        # Initialize the matrices to be filled in during the following loop.
        atom_feature_matrix = np.zeros(n_atom, dtype='int')
        bond_feature_matrix = np.zeros(n_bond, dtype='int')
        bond_distance_matrix = np.zeros(n_bond, dtype=np.float32)
        connectivity = np.zeros((n_bond, 2), dtype='int')

        # Hopefully we've filtered out all problem mols by now.
        if mol is None:
            raise RuntimeError("Issue in loading mol")

        # Get a list of the atoms in the molecule.
        atom_seq = mol.GetAtoms()
        atoms = [atom_seq[i] for i in range(n_atom)]

        # Here we loop over each atom, and the inner loop iterates over each
        # neighbor of the current atom.
        bond_index = 0  # keep track of our current bond.
        for n, atom in enumerate(atoms):

            # update atom feature matrix
            atom_feature_matrix[n] = self.atom_tokenizer(
                self.atom_features(atom))

            # if n_neighbors is greater than total atoms, then each atom is a
            # neighbor.
            if (self.n_neighbors + 1) > len(mol.GetAtoms()):
                neighbor_end_index = len(mol.GetAtoms())
            else:
                neighbor_end_index = (self.n_neighbors + 1)

            distance_atom = distance_matrix[n, :]
            cutoff_end_index = distance_atom[distance_atom < self.cutoff].size

            end_index = min(neighbor_end_index, cutoff_end_index)

            # Loop over each of the nearest neighbors

            neighbor_inds = distance_matrix[n, :].argsort()[1:end_index]
            if len(neighbor_inds)==0: neighbor_inds = [n]
            for neighbor in neighbor_inds:

                # update bond feature matrix
                bond = mol.GetBondBetweenAtoms(n, int(neighbor))
                if bond is None:
                    bond_feature_matrix[bond_index] = 0
                else:
                    rev = False if bond.GetBeginAtomIdx() == n else True
                    bond_feature_matrix[bond_index] = self.bond_tokenizer(
                        self.bond_features(bond, flipped=rev))

                distance = distance_matrix[n, neighbor]
                bond_distance_matrix[bond_index] = distance

                # update connectivity matrix
                connectivity[bond_index, 0] = n
                connectivity[bond_index, 1] = neighbor

                bond_index += 1
        print(connectivity)

        return {
            'n_atom': n_atom,
            'n_bond': n_bond,
            'atom': atom_feature_matrix,
            'bond': bond_feature_matrix,
            'distance': bond_distance_matrix,
            'connectivity': connectivity,
        }


class MolBPreprocessor(MolPreprocessor):
    """
    This is a subclass of Molpreprocessor that preprocessor molecule with
    bond property target
    """
    def __init__(self, **kwargs):
        """
        A preprocessor class that also returns bond_target_matrix, besides the bond matrix
        returned by MolPreprocessor. The bond_target_matrix is then used as ref to reduce molecule
        to bond property
        """
        super(MolBPreprocessor, self).__init__(**kwargs)

    def construct_feature_matrices(self, entry):
        """
        Given an entry contining rdkit molecule, bond_index and for the target property,
        return atom
        feature matrices, bond feature matrices, distance matrices, connectivity matrices and bond
        ref matrices.

        returns
        dict with entries
        see MolPreproccessor
        'bond_index' : ref array to the bond index
        """
        mol, bond_index_array = entry

        n_atom = len(mol.GetAtoms())
        n_pro = len(bond_index_array)

        # n_bond is actually the number of atom-atom pairs, so this is defined
        # by the number of neighbors for each atom.
        #if there is cutoff,
        distance_matrix = Chem.Get3DDistanceMatrix(mol)

        if self.n_neighbors <= (n_atom - 1):
            n_bond = self.n_neighbors * n_atom
        else:
            # If there are fewer atoms than n_neighbors, all atoms will be
            # connected
            n_bond = distance_matrix[(distance_matrix < self.cutoff) & (distance_matrix != 0)].size

        if n_bond == 0: n_bond = 1

        # Initialize the matrices to be filled in during the following loop.
        atom_feature_matrix = np.zeros(n_atom, dtype='int')
        bond_feature_matrix = np.zeros(n_bond, dtype='int')
        bond_distance_matrix = np.zeros(n_bond, dtype=np.float32)
        bond_index_matrix = np.full(n_bond, -1, dtype='int')
        connectivity = np.zeros((n_bond, 2), dtype='int')

        # Hopefully we've filtered out all problem mols by now.
        if mol is None:
            raise RuntimeError("Issue in loading mol")

        # Get a list of the atoms in the molecule.
        atom_seq = mol.GetAtoms()
        atoms = [atom_seq[i] for i in range(n_atom)]

        # Here we loop over each atom, and the inner loop iterates over each
        # neighbor of the current atom.
        bond_index = 0  # keep track of our current bond.
        for n, atom in enumerate(atoms):
            # update atom feature matrix
            atom_feature_matrix[n] = self.atom_tokenizer(
                self.atom_features(atom))

            # if n_neighbors is greater than total atoms, then each atom is a
            # neighbor.
            if (self.n_neighbors + 1) > len(mol.GetAtoms()):
                neighbor_end_index = len(mol.GetAtoms())
            else:
                neighbor_end_index = (self.n_neighbors + 1)

            distance_atom = distance_matrix[n, :]
            cutoff_end_index = distance_atom[distance_atom < self.cutoff].size

            end_index = min(neighbor_end_index, cutoff_end_index)

            # Loop over each of the nearest neighbors

            neighbor_inds = distance_matrix[n, :].argsort()[1:end_index]
            if len(neighbor_inds)==0: neighbor_inds = [n]
            for neighbor in neighbor_inds:

                # update bond feature matrix
                bond = mol.GetBondBetweenAtoms(n, int(neighbor))
                if bond is None:
                    bond_feature_matrix[bond_index] = 0
                else:
                    rev = False if bond.GetBeginAtomIdx() == n else True
                    bond_feature_matrix[bond_index] = self.bond_tokenizer(
                        self.bond_features(bond, flipped=rev))
                    try:
                        bond_index_matrix[bond_index] = bond_index_array.tolist().index(bond.GetIdx())
                    except:
                        pass

                distance = distance_matrix[n, neighbor]
                bond_distance_matrix[bond_index] = distance

                # update connectivity matrix
                connectivity[bond_index, 0] = n
                connectivity[bond_index, 1] = neighbor

                bond_index += 1
        return {
            'n_atom': n_atom,
            'n_bond': n_bond,
            'n_pro': n_pro,
            'atom': atom_feature_matrix,
            'bond': bond_feature_matrix,
            'distance': bond_distance_matrix,
            'connectivity': connectivity,
            'bond_index': bond_index_matrix,
        }

class MolAPreprocessor(MolPreprocessor):
    """
    This is a subclass of Molpreprocessor that preprocessor molecule with
    bond property target
    """
    def __init__(self, **kwargs):
        """
        A preprocessor class that also returns bond_target_matrix, besides the bond matrix
        returned by MolPreprocessor. The bond_target_matrix is then used as ref to reduce molecule
        to bond property
        """
        super(MolAPreprocessor, self).__init__(**kwargs)

    def construct_feature_matrices(self, entry):
        """
        Given an entry contining rdkit molecule, bond_index and for the target property,
        return atom
        feature matrices, bond feature matrices, distance matrices, connectivity matrices and bond
        ref matrices.

        returns
        dict with entries
        see MolPreproccessor
        'bond_index' : ref array to the bond index
        """
        mol, atom_index_array = entry

        n_atom = len(mol.GetAtoms())
        n_pro = len(atom_index_array)

        # n_bond is actually the number of atom-atom pairs, so this is defined
        # by the number of neighbors for each atom.
        #if there is cutoff,
        distance_matrix = Chem.Get3DDistanceMatrix(mol)

        #if self.n_neighbors <= (n_atom - 1):
        #    n_bond = self.n_neighbors * n_atom
        #else:
            # If there are fewer atoms than n_neighbors, all atoms will be
            # connected
        n_bond = distance_matrix[(distance_matrix < self.cutoff) & (distance_matrix != 0)].size

        if n_bond == 0: n_bond = 1

        # Initialize the matrices to be filled in during the following loop.
        atom_feature_matrix = np.zeros(n_atom, dtype='int')
        bond_feature_matrix = np.zeros(n_bond, dtype='int')
        bond_distance_matrix = np.zeros(n_bond, dtype=np.float32)
        atom_index_matrix = np.full(n_atom, -1, dtype='int')
        connectivity = np.zeros((n_bond, 2), dtype='int')

        # Hopefully we've filtered out all problem mols by now.
        if mol is None:
            raise RuntimeError("Issue in loading mol")

        # Get a list of the atoms in the molecule.
        atom_seq = mol.GetAtoms()
        atoms = [atom_seq[i] for i in range(n_atom)]

        # Here we loop over each atom, and the inner loop iterates over each
        # neighbor of the current atom.
        bond_index = 0  # keep track of our current bond.
        for n, atom in enumerate(atoms):
            # update atom feature matrix
            atom_feature_matrix[n] = self.atom_tokenizer(
                self.atom_features(atom))
            try:
                atom_index_matrix[n] = atom_index_array.tolist().index(atom.GetIdx())
            except:
                pass
            # if n_neighbors is greater than total atoms, then each atom is a
            # neighbor.
            if (self.n_neighbors + 1) > len(mol.GetAtoms()):
                neighbor_end_index = len(mol.GetAtoms())
            else:
                neighbor_end_index = (self.n_neighbors + 1)

            distance_atom = distance_matrix[n, :]
            cutoff_end_index = distance_atom[distance_atom < self.cutoff].size

            end_index = min(neighbor_end_index, cutoff_end_index)

            # Loop over each of the nearest neighbors

            neighbor_inds = distance_matrix[n, :].argsort()[1:end_index]
            if len(neighbor_inds)==0: neighbor_inds = [n]
            for neighbor in neighbor_inds:

                # update bond feature matrix
                bond = mol.GetBondBetweenAtoms(n, int(neighbor))
                if bond is None:
                    bond_feature_matrix[bond_index] = 0
                else:
                    rev = False if bond.GetBeginAtomIdx() == n else True
                    bond_feature_matrix[bond_index] = self.bond_tokenizer(
                        self.bond_features(bond, flipped=rev))

                distance = distance_matrix[n, neighbor]
                bond_distance_matrix[bond_index] = distance

                # update connectivity matrix
                connectivity[bond_index, 0] = n
                connectivity[bond_index, 1] = neighbor

                bond_index += 1
        return {
            'n_atom': n_atom,
            'n_bond': n_bond,
            'n_pro': n_pro,
            'atom': atom_feature_matrix,
            'bond': bond_feature_matrix,
            'distance': bond_distance_matrix,
            'connectivity': connectivity,
            'atom_index': atom_index_matrix,
        }


# TODO: rewrite this
# class LaplacianSmilesPreprocessor(SmilesPreprocessor):
#     """ Extends the SmilesPreprocessor class to also return eigenvalues and
#     eigenvectors of the graph laplacian matrix.
#
#     Example:
#     >>> preprocessor = SmilesPreprocessor(
#     >>>     max_atoms=55, max_bonds=62, max_degree=4, explicit_hs=False)
#     >>> atom, connectivity, eigenvalues, eigenvectors = preprocessor.fit(
#             data.smiles)
#     """
#
#     def preprocess(self, smiles_iterator, train=True):
#
#         self.atom_tokenizer.train = train
#         self.bond_tokenizer.train = train
#
#         for smiles in tqdm(smiles_iterator):
#             G = self._mol_to_nx(smiles)
#             A = self._get_atom_feature_matrix(G)
#             C = self._get_connectivity_matrix(G)
#             W, V = self._get_laplacian_spectral_decomp(G)
#             yield A, C, W, V
#
#
#     def _get_laplacian_spectral_decomp(self, G):
#         """ Return the eigenvalues and eigenvectors of the graph G, padded to
#         `self.max_atoms`.
#         """
#
#         w0 = np.zeros((self.max_atoms, 1))
#         v0 = np.zeros((self.max_atoms, self.max_atoms))
#
#         w, v = eigh(nx.laplacian_matrix(G).todense())
#
#         num_atoms = len(v)
#
#         w0[:num_atoms, 0] = w
#         v0[:num_atoms, :num_atoms] = v
#
#         return w0, v0
#
#
#     def fit(self, smiles_iterator):
#         results = self._fit(smiles_iterator)
#         return {'atom': results[0],
#                 'connectivity': results[1],
#                 'w': results[2],
#                 'v': results[3]}
#
#
#     def predict(self, smiles_iterator):
#         results = self._predict(smiles_iterator)
#         return {'atom': results[0],
#                 'connectivity': results[1],
#                 'w': results[2],
#                 'v': results[3]}


def get_max_atom_bond_size(smiles_iterator, explicit_hs=True):
    """ Convienence function to get max_atoms, max_bonds for a set of input
    SMILES """

    max_atoms = 0
    max_bonds = 0
    for smiles in tqdm(smiles_iterator):
        mol = MolFromSmiles(smiles)
        if explicit_hs:
            mol = AddHs(mol)
        max_atoms = max([max_atoms, len(mol.GetAtoms())])
        max_bonds = max([max_bonds, len(mol.GetBonds())])

    return dict(max_atoms=max_atoms, max_bonds=max_bonds*2)


def canonicalize_smiles(smiles, isomeric=True, sanitize=True):
    try:
        mol = MolFromSmiles(smiles, sanitize=sanitize)
        return MolToSmiles(mol, isomericSmiles=isomeric)
    except Exception:
        pass
