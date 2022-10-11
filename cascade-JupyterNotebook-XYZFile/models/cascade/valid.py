from rdkit import Chem

allowed_atoms = [1, 6, 7, 8, 9, 15, 16, 17]

def validate_smiles(smiles):

    #check if the smile string is Valid
    if not smiles: return False

    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        return False

    for a in mol.GetAtoms():
        if a.GetAtomicNum() not in allowed_atoms:
            return False
        if a.GetFormalCharge() != 0:
            return False
    return True
