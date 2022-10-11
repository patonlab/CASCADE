class Tokenizer(object):
    """ A class to turn arbitrary inputs into integer classes. """

    def __init__(self):
        # the default class for an unseen entry during test-time
        self._data = {'unk': 1}
        self.num_classes = 1
        self.train = True
        self.unknown = []

    def __call__(self, item):
        """ Check to see if the Tokenizer has seen `item` before, and if so,
        return the integer class associated with it. Otherwise, if we're
        training, create a new integer class, otherwise return the 'unknown'
        class.

        """
        try:
            return self._data[item]

        except KeyError:
            if self.train:
                self._add_token(item)
                return self(item)

            else:
                # Record the unknown item, then return the unknown label
                self.unknown += [item]
                return self._data['unk']

    def _add_token(self, item):
        self.num_classes += 1
        self._data[item] = self.num_classes

# The rest of the methods in this module are specific functions for computing
# atom and bond features. New ones can be easily added though, and these are
# passed directly to the Preprocessor class.

def get_ring_size(obj, max_size=12):
    if not obj.IsInRing():
        return 0
    else:
        for i in range(max_size):
            if obj.IsInRingSize(i):
                return i
        else:
            return 'max'

def atom_features(atom):
    return atom.GetAtomicNum()


def atom_features_v1(atom):
    """ Return an integer hash representing the atom type
    """

    return str((
        atom.GetSymbol(),
        atom.GetDegree(),
        atom.GetTotalNumHs(),
        atom.GetImplicitValence(),
        atom.GetIsAromatic(),
    ))


def atom_features_v2(atom):

    props = ['GetChiralTag', 'GetDegree', 'GetExplicitValence',
             'GetFormalCharge', 'GetHybridization', 'GetImplicitValence',
             'GetIsAromatic', 'GetNoImplicit', 'GetNumExplicitHs',
             'GetNumImplicitHs', 'GetNumRadicalElectrons', 'GetSymbol',
             'GetTotalDegree', 'GetTotalNumHs', 'GetTotalValence']

    atom_type = [getattr(atom, prop)() for prop in props]
    atom_type += [get_ring_size(atom)]

    return str(tuple(atom_type))


def bond_features_v1(bond, **kwargs):
    """ Return an integer hash representing the bond type.

    flipped : bool
        Only valid for 'v3' version, whether to swap the begin and end atom types

    """

    return str((
        bond.GetBondType(),
        bond.GetIsConjugated(),
        bond.IsInRing(),
        sorted([
            bond.GetBeginAtom().GetSymbol(),
            bond.GetEndAtom().GetSymbol()]),
        ))


def bond_features_v2(bond, **kwargs):

    return str((
        bond.GetBondType(),
        bond.GetIsConjugated(),
        bond.GetStereo(),
        get_ring_size(bond),
        sorted([
            bond.GetBeginAtom().GetSymbol(),
            bond.GetEndAtom().GetSymbol()]),
        ))


def bond_features_v3(bond, flipped=False):

    if not flipped:
        start_atom = atom_features(bond.GetBeginAtom())
        end_atom = atom_features(bond.GetEndAtom())

    else:
        start_atom = atom_features(bond.GetEndAtom())
        end_atom = atom_features(bond.GetBeginAtom())

    return str((
        bond.GetBondType(),
        bond.GetIsConjugated(),
        bond.GetStereo(),
        get_ring_size(bond),
        bond.GetEndAtom().GetSymbol(),
        start_atom,
        end_atom))
