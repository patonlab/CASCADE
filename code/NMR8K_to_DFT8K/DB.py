from ase.db.sqlite import SQLite3Database
from ase import Atoms
import rdkit
from rdkit import Chem
import os,sys,glob
import random
import numpy as np
from tqdm import tqdm

"""
Classes to interact with database
"""

HARTREE_TO_KCAL = 627.509

class ASEDB(SQLite3Database):

    def __init__(self, db_file, create_indices=True, use_lock_file=True,
                 append=True, serial=False):
        if not append and os.path.isfile(db_file):
            os.remove(db_file)

        self.dbpath = db_file
        super(ASEDB, self).__init__(db_file, create_indices, use_lock_file, serial)

    def str_Ids(self):
        try:
            str_ids = []
            for row in self.select():
                str_ids.append(row.str_id)
        except:
            raise AttributeError('Current data base %s does not have column as str_id' %
                                 self.dbpath)
        else:
            return str_ids

    def ids(self):
        return list(range(1, len(self)+1))

    #FIXME cannot update atom and data at the same time
    def update_by_str(self, str_id, **key_values):
        try:
            row = self.get(str_id=str_id)
        except:
            e = sys.exc_info()[0]
            raise AttributeError("cannot get id from str_id: %s" % e)
        else:
            id = row.id

        if 'data' in key_values.keys():
            atoms = row.toatoms()
            
            kvp = row.key_value_pairs

            new_id = self.write(atoms, **kvp, data=key_values['data'])

            del key_values['data']
            del self[id]
        if 'atoms' in key_values.keys():
            atoms = key_values['atoms']
            if not isinstance(atoms, Atoms):
                raise TypeError('atoms must be an ase.Atoms object')
           
            kvp = row.key_value_pairs  
            
            if hasattr(row, 'data'):        
                new_id = self.write(atoms, **kvp, data=row.data)
            else:
                new_id = self.write(atoms, **kvp)

            del key_values['atoms']
            del self[id] 

        else: new_id = id
        
        if key_values: self.update(new_id, **key_values)
        return new_id

    def update_from_sdf(self, str_id, sdf):
        '''
        '''
        rdkit_m = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0]
        
        try:
            enthalpy = float(rdkit_m.GetProp('enthalpy'))
        except:
            enthalpy = None

        ase_m = Rdkit_2_ASE(rdkit_m)
        if enthalpy:
            new_id = self.update_by_str(str_id, atoms=ase_m, enthalpy=enthalpy)
        else:
            new_id = self.update_by_str(str_id, atoms=ase_m)

        return new_id


class NMRDB:
    """
    """
    def __init__(self, source=None):
        """
        """
        if not (isinstance(source, str) and '.db' in source):
            raise TypeError('The source database must be a .db file')
        
        self.db = ASEDB(source)


    def assign_computer(self, n_points, n_batches, calculator):
        """
        """
        
        rows = [x for x in self.db.select(filter=lambda x: x.opt is False)]
        
        if len(rows) > n_points:
            rows = random.sample(rows, n_points)

        if len(rows) == 0: return []

        for r in tqdm(rows, total=len(rows)): self.db.update(r.id, opt='running')

        m_list = [r.str_id for r in rows]

        batch_size = int(n_points/n_batches)
        m_batches = [m_list[batch_size*i:batch_size*(i+1)] for i in range(n_batches)]
        m_batches[-1].extend(m_list[batch_size*n_batches:])

        computers = []

        i = len(glob.glob('mol_batch*'))  
        for batch in tqdm(m_batches, total=n_batches):
            i += 1
            calculator_tmp = calculator.copy()
            calculator_tmp.load_batch(batch)
            calculator_tmp.s_type = 'mol'
            calculator_tmp.path = ('mol_batch_'+calculator_tmp.operation+str(i))
            computers.append(calculator_tmp)

        return computers

def ASE_2_Rdkit (ase_mol, id, sanitize=False):
    """
    """
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats('xyz', 'sdf')
    
    tmp_xyz = tempfile.NamedTemporaryFile(delete=True)
    ase.io.write(tmp_xyz.name, ase_mol, format='xyz')
    try:
        with open(os.devnull, 'w') as dev:
            with redirect_stdout(dev) as h:
                babel_mol = openbabel.OBMol()
                obConversion.ReadFile(babel_mol, tmp_xyz.name)
                tmp_xyz.close()

                tmp_sdf = tempfile.NamedTemporaryFile(delete=True)
                obConversion.WriteFile(babel_mol, tmp_sdf.name)

        rdkit_mol = Chem.SDMolSupplier(tmp_sdf.name, removeHs=False, sanitize=sanitize)[0]
        rdkit_mol.SetProp('_Name', str(id))
    except Exception as e:
        raise RuntimeError('Cannot convert ase mol to rdkit mol: ' + e)

    tmp_sdf.close()

    return rdkit_mol

def Rdkit_2_ASE (rdkit_mol):
    """
    """
    tmp_sdf = tempfile.NamedTemporaryFile(delete=True)
    writer = Chem.rdmolfiles.SDWriter(tmp_sdf.name)
    try:
        writer.write(rdkit_mol)
        writer.close()
        ase_mol = ase.io.read(tmp_sdf.name, index=0, format='sdf')
    except Exception as e:
        raise RuntimeError('Cannot convert rdkit_mol to ase_mol: ' + e)

    tmp_sdf.close()
    return ase_mol 

def canonicalize_smiles(smiles):
    """ Return a consisten SMILES representation for the given molecule """
    mol = rdkit.Chem.MolFromSmiles(smiles)
    return rdkit.Chem.MolToSmiles(mol)

    


        
        

        
                
            
