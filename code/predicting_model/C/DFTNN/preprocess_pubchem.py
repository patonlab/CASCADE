import sys

#change the path into where the nfp folder is
sys.path.append('../..')

import pandas as pd
import numpy as np
import gzip

import warnings
from tqdm import tqdm

import gzip
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ForwardSDMolSupplier

from itertools import islice

from nfp.preprocessing import MolAPreprocessor, GraphSequence

mols = []
with gzip.open('../../../../data/DFT8K/DFT.sdf.gz', 'r') as sdfile:
    mol_supplier = ForwardSDMolSupplier(sdfile, removeHs=False, sanitize=False)
    for mol in tqdm(mol_supplier):
        if mol:
            mols += [(int(mol.GetProp('_Name')), mol, mol.GetNumAtoms())]

mols = pd.DataFrame(mols, columns=['mol_id', 'Mol', 'n_atoms'])
mols = mols.set_index('mol_id', drop=True)

df = pd.read_csv('../../../../data/DFT8K/DFT8K.csv.gz', index_col=0)
#only choose C and H
df = df.loc[df.atom_type == 6]

df['Mol'] = mols.reindex(df.mol_id).Mol.values

grouped_df = df.groupby(['mol_id'])
df_Shift = []
for mol_id,df in grouped_df:
    df_Shift.append([mol_id, df.atom_index.values.astype('int'), df.Shift.values.astype(np.float32)])
    if len(df.atom_index.values) != len(set(df.atom_index.values)):
        print(mol_id)

df_Shift = pd.DataFrame(df_Shift, columns=['mol_id', 'atom_index', 'Shift']) 

test = df_Shift.sample(n=500, random_state=666)
valid = df_Shift[~df_Shift.mol_id.isin(test.mol_id)].sample(n=500, random_state=666)
train = df_Shift[
    (~df_Shift.mol_id.isin(test.mol_id) & ~df_Shift.mol_id.isin(valid.mol_id))
              ]
test = test.set_index('mol_id')
valid = valid.set_index('mol_id')
train = train.set_index('mol_id')

test = mols.reindex(test.index).join(test[['atom_index', 'Shift']]).dropna()
valid = mols.reindex(valid.index).join(valid[['atom_index', 'Shift']]).dropna()
train = mols.reindex(train.index).join(train[['atom_index', 'Shift']]).dropna()

test.to_pickle('test.pkl.gz', compression='gzip')
valid.to_pickle('valid.pkl.gz', compression='gzip')
train.to_pickle('train.pkl.gz', compression='gzip')

# Preprocess molecules
def atomic_number_tokenizer(atom):
    return atom.GetAtomicNum()
def Mol_iter(df):
    for index,r in df.iterrows():
        yield(r['Mol'], r['atom_index'])

preprocessor = MolAPreprocessor(
    n_neighbors=100, cutoff=5, atom_features=atomic_number_tokenizer)

inputs_train = preprocessor.fit(Mol_iter(train))
inputs_valid = preprocessor.predict(Mol_iter(valid))
inputs_test = preprocessor.predict(Mol_iter(test))

import pickle
with open('processed_inputs.p', 'wb') as file:        
    pickle.dump({
        'inputs_train': inputs_train,
        'inputs_valid': inputs_valid,
        'inputs_test': inputs_test,
        'preprocessor': preprocessor,
    }, file)
