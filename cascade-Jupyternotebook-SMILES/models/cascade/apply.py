import sys,io,os
import pandas as pd
import numpy as np
import gzip, pickle, argparse, warnings
import pickle
import math

from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ForwardSDMolSupplier

from itertools import islice

from nfp.preprocessing import MolAPreprocessor, GraphSequence
from .genConf import genConf

import keras
import keras.backend as K

from keras.callbacks import ModelCheckpoint, CSVLogger, LearningRateScheduler

from keras.layers import (Input, Embedding, Dense, BatchNormalization,
                                 Concatenate, Multiply, Add)

from keras.models import Model, load_model

from nfp.layers import (MessageLayer, GRUStep, Squeeze, EdgeNetwork,
                               ReduceBondToPro, ReduceBondToAtom,
                               GatherAtomToBond, ReduceAtomToPro)
from nfp.models import GraphModel

def to_C(atom):
    neighbors = [x.GetAtomicNum() for x in atom.GetNeighbors()]
    if 6 in neighbors:
        return True
    else:
        return False

def Mol_iter(df):
    for index,r in df.iterrows():
        yield(r['Mol'], r['atom_index'])

def preprocess_C(mols, preprocessor, keep_all_cf=False):
    mols_id = []
    confs_id = []
    mols_conf = []
    for i,m in enumerate(mols):
        try:
            mol,ids,nr = genConf(m, rms=-1, nc=200, efilter=10.0, rmspost=0.5)
        except Exception as e:
            msg = "cannot generate FF conformers"
            print(msg)
            mols_conf.append(None)
            continue

        mols_i = [Chem.RWMol(mol) for x in ids]
        if keep_all_cf:
            for m_i,id in zip(mols_i, ids):
                m_i.SetProp('E', str(id[0]))
                m_i.SetProp('_Name', '{}_{}'.format(i, id[1]))
                for ia in range(mol.GetNumAtoms()):
                    m_i.GetConformer().SetAtomPosition(ia, mol.GetConformer(id[1]).GetAtomPosition(ia))
            mols_conf.extend(mols_i)

        else:
            id = ids[0]
            m.SetProp('E', str(id[0]))
            m.SetProp('_Name', '{}_{}'.format(i, id[1]))

    df = []
    for m in mols_conf:
        if m:
            E = float(m.GetProp('E'))
            m_id, cf_id = [ int(x) for x in m.GetProp('_Name').split('_') ]
            Cs = [x for x in m.GetAtoms() if x.GetAtomicNum()==6]
            C_index = np.array([x.GetIdx() for x in Cs]).astype(int)
            df.append([m_id, m, m.GetNumAtoms(), C_index, E, cf_id])

    df = pd.DataFrame(df, columns=['mol_id', 'Mol', 'n_atoms', 'atom_index', 'relative_E', 'cf_id'])

    inputs = preprocessor.predict(Mol_iter(df))

    return inputs, df, mols_conf

class RBFSequence(GraphSequence):
    def process_data(self, batch_data):
        batch_data['distance_rbf'] = self.rbf_expansion(batch_data['distance'])

        offset = self._compute_stacked_offsets(
            batch_data['n_pro'], batch_data['n_atom'])

        offset = np.where(batch_data['atom_index']>=0, offset, 0)
        batch_data['atom_index'] += offset

        del batch_data['n_atom']
        del batch_data['n_bond']
        del batch_data['distance']
        return batch_data

    def _compute_stacked_offsets(self, sizes, repeats):
        return np.repeat(np.cumsum(np.hstack([0, sizes[:-1]])), repeats)

    def rbf_expansion(self, distances, mu=0, delta=0.1, kmax=256):
        k = np.arange(0, kmax)
        logits = -(np.atleast_2d(distances).T - (-mu + delta * k))**2 / delta
        return np.exp(logits)

def evaluate_C(inputs, preprocessor, model):
    batch_size = 32
    evaluate_sequence = RBFSequence(inputs, batch_size=batch_size)

    predicted = []

    for x in evaluate_sequence:
        out = model.predict_on_batch(x)
        out = np.concatenate(out)
        predicted.extend(out)

    return predicted

def predict_NMR_C(smiles, model):

    mols = [Chem.MolFromSmiles(smiles)]
    for m in mols: AllChem.EmbedMolecule(m, useRandomCoords=True)
    mols = [Chem.AddHs(m, addCoords=True) for m in mols]


    with open(os.path.join('cascade', 'preprocessor.p'), 'rb') as ft:
        preprocessor = pickle.load(ft)['preprocessor']

    inputs, df, mols = preprocess_C(mols, preprocessor, True)

    if len(inputs) == 0:
        raise RuntimeError('Failed to find any conformer for the given molecule')

    predicted = evaluate_C(inputs, preprocessor, model)

    spread_df = pd.DataFrame([], columns=['mol_id', 'atom_index', 'relative_E', 'cf_id'])
    for _,r in df.iterrows():
        mol_id=[r.mol_id]*len(r.atom_index)
        cf_id = [r.cf_id]*len(r.atom_index)
        E = [r.relative_E]*len(r.atom_index)
        df_mol = pd.DataFrame(data={'mol_id': mol_id, 'atom_index': r.atom_index, 'relative_E':E,  'cf_id': cf_id})
        spread_df = pd.concat([spread_df, df_mol], sort=True)
    spread_df['predicted'] = predicted
    spread_df['b_weight'] = spread_df.relative_E.apply(lambda x: math.exp(-x/(0.001987*298.15)))

    df_group = spread_df.set_index(['mol_id', 'atom_index', 'cf_id']).groupby(level=[0,1])
    final = []
    for (m_id, a_id),df in df_group:
        weighted_shift = df.apply(lambda x: x['b_weight']*x['predicted'], axis=1).sum()/df.b_weight.sum()
        final.append([m_id, a_id, weighted_shift])

    final = pd.DataFrame(final, columns=['mol_id', 'atom_index', 'Shift'])
    final['atom_index'] = final.atom_index.apply(lambda x: x+1)
    final = final.round(2).astype(dtype={'atom_index':'int'})
    spread_df['atom_index'] = spread_df.atom_index.apply(lambda x: x+1)
    spread_df = spread_df.round(2)

    return(mols, final, spread_df)


def preprocess_H(mols, preprocessor, keep_all_cf=False):
    mols_id = []
    confs_id = []
    mols_conf = []
    for i,m in enumerate(mols):
        try:
            mol,ids,nr = genConf(m, rms=-1, nc=200, efilter=10.0, rmspost=0.5)
        except Exception as e:
            msg = "cannot generate FF conformers"
            print(msg)
            mols_conf.append(None)
            continue

        mols_i = [Chem.RWMol(mol) for x in ids]
        if keep_all_cf:
            for m_i,id in zip(mols_i, ids):
                m_i.SetProp('E', str(id[0]))
                m_i.SetProp('_Name', '{}_{}'.format(i, id[1]))
                for ia in range(mol.GetNumAtoms()):
                    m_i.GetConformer().SetAtomPosition(ia, mol.GetConformer(id[1]).GetAtomPosition(ia))
            mols_conf.extend(mols_i)

        else:
            id = ids[0]
            m.SetProp('E', str(id[0]))
            m.SetProp('_Name', '{}_{}'.format(i, id[1]))

    df = []
    for m in mols_conf:
        if m:
            E = float(m.GetProp('E'))
            m_id, cf_id = [ int(x) for x in m.GetProp('_Name').split('_') ]
            Hs = [x for x in m.GetAtoms() if x.GetAtomicNum()==1]
            H_index = np.array([x.GetIdx() for x in Hs]).astype(int)
            df.append([m_id, m, m.GetNumAtoms(), H_index, E, cf_id])

    df = pd.DataFrame(df, columns=['mol_id', 'Mol', 'n_atoms', 'atom_index', 'relative_E', 'cf_id'])

    inputs = preprocessor.predict(Mol_iter(df))

    return inputs, df, mols_conf

def evaluate_H(inputs, preprocessor, model):
    batch_size = 32
    evaluate_sequence = RBFSequence(inputs, batch_size=batch_size)

    predicted = []

    for x in evaluate_sequence:
        out = model.predict_on_batch(x)
        out = np.concatenate(out)
        predicted.extend(out)

    return predicted

def predict_NMR_H(smiles, model):

    mols = [Chem.MolFromSmiles(smiles)]
    for m in mols: AllChem.EmbedMolecule(m, useRandomCoords=True)
    mols = [Chem.AddHs(m, addCoords=True) for m in mols]


    with open(os.path.join('cascade', 'preprocessor.p'), 'rb') as ft:
        preprocessor = pickle.load(ft)['preprocessor']

    inputs, df, mols = preprocess_H(mols, preprocessor, True)
    predicted = evaluate_H(inputs, preprocessor, model)

    spread_df = pd.DataFrame([], columns=['mol_id', 'atom_index', 'relative_E', 'cf_id'])
    for _,r in df.iterrows():
        mol_id=[r.mol_id]*len(r.atom_index)
        cf_id = [r.cf_id]*len(r.atom_index)
        E = [r.relative_E]*len(r.atom_index)
        df_mol = pd.DataFrame(data={'mol_id': mol_id, 'atom_index': r.atom_index, 'relative_E':E,  'cf_id': cf_id})
        spread_df = pd.concat([spread_df, df_mol], sort=True)
    spread_df['predicted'] = predicted
    spread_df['b_weight'] = spread_df.relative_E.apply(lambda x: math.exp(-x/(0.001987*298.15)))

    df_group = spread_df.set_index(['mol_id', 'atom_index', 'cf_id']).groupby(level=[0,1])
    final = []
    for (m_id, a_id),df in df_group:
        weighted_shift = df.apply(lambda x: x['b_weight']*x['predicted'], axis=1).sum()/df.b_weight.sum()
        final.append([m_id, a_id, weighted_shift])

    final = pd.DataFrame(final, columns=['mol_id', 'atom_index', 'Shift'])
    final['atom_index'] = final.atom_index.apply(lambda x: x+1)
    final = final.round(2).astype(dtype={'atom_index':'int'})
    spread_df['atom_index'] = spread_df.atom_index.apply(lambda x: x+1)
    spread_df = spread_df.round(2)

    return(mols, final, spread_df)

