# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 14:14:08 2018

This is a module handling file format converting

@author: Yanfei Guan
"""

import pyQRC as pq
from goodvibes import GoodVibes as GV
from pybel import readfile
from rdkit import Chem
import pybel
import os,abc,re

#Some useful arrays
periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
    "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

def elementID(massno):
        if massno < len(periodictable): return periodictable[massno]
        else: return "XX"

class file_parser(object):
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractproperty
    def rdkit_mol(self):
        pass
    
    @abc.abstractproperty
    def pybel_mol(self):
        pass
    
    def write_gau_com (self,
                       name='file_parser_auto_com',
                       route='',
                       footer='',
                       charge=0,
                       mult=1,
                       memory='96GB',
                       nprocs=24,
                       ):
        if not isinstance(route, basestring): raise TypeError("route must be a string")
        if not isinstance(footer, basestring): raise TypeError("footer must be a string")
        
        with open(name+'.com', 'w') as fileout:
            fileout.write("%mem="+memory+"\n")
            fileout.write("nprocshared="+str(nprocs)+"\n")
    
            fileout.write(route+"\n\n")
            fileout.write(name+"\n\n")
            fileout.write("%d %d\n" % (charge, mult))
    
            for j in range(0, self.rdkit_mol.GetNumAtoms()):
                atom, pos = [self.rdkit_mol.GetAtoms()[j], 
                             self.rdkit_mol.GetConformer().GetAtomPosition(j)]
                fileout.write('{:<4} {:10.6f} {:10.6f} {:10.6f}'.format(atom.GetSymbol(), pos.x, pos.y, pos.z))
                fileout.write("\n")
            
            fileout.write("\n")
            fileout.write(footer+"\n\n\n")
        
        return name+'.com'
    
    def write_sdf(self, output=None):
        if output is None: output = self.name+'.sdf'
        self._pybel_mol.write(format='sdf', filename=output, overwrite=True)
        
    def write_orca_inp(self,
                      name='file_parser_auto_in',
                      route='DLPNO-CCSD(T) Extrapolate(2/3,cc)',
                      charge=0,
                      mult=1,
                      memory='4GB',
                      nprocs=24,
                      ):
        if not isinstance(route, basestring): raise TypeError("route must be a string")
        with open(name+'.inp', 'w') as fileout:
            fileout.write("%maxcore "+str(int(memory.split('GB')[0])*1000)+"\n")
            fileout.write("% pal nprocs "+str(nprocs)+" end\n")
            fileout.write("! "+route+"\n\n")
            fileout.write("%scf maxiter 500\n   end\n% mdci\n")
            fileout.write("    TCutPairs 1e-4\n    TCutPNO 3.33e-7\n    TCutMKN 1e-3\n")
            fileout.write("  Density None\nend\n")
            fileout.write("% output\n  printlevel mini\n  print[ P_SCFInfo ] 1\n")
            fileout.write("  print[ P_SCFIterInfo ] 1\n  print[ P_OrbEn ] 0\n  print[ P_Cartesian ] 0\nend\n")
            fileout.write("% elprop\n  Dipole False\nend\n")

            fileout.write('*xyz '+str(charge)+" "+str(mult)+"\n")
            mol = self._rdkit_mol
            
            for j in range(0, mol.GetNumAtoms()):
                atom, pos = mol.GetAtoms()[j], mol.GetConformer().GetAtomPosition(j)
                fileout.write('{:<4} {:10.6f} {:10.6f} {:10.6f}'.format(atom.GetSymbol(), pos.x, pos.y, pos.z))
                fileout.write("\n")
            fileout.write("*\n")
            
    def get_rdkit_mol(self):
        '''
        If the object has a pybel molecule, use this to get rdkit_mol
        '''
        temp_sdf = self.name+'_temp.sdf'
        self.write_sdf(output=temp_sdf)
        self._rdkit_mol = Chem.SDMolSupplier(temp_sdf, removeHs=False, sanitize=False)[0]
        os.remove(temp_sdf)
        
    def get_pybel_mol(self):
        '''
        If the object has a rdkit molecule, use this to get a pybel molecule
        '''
        temp_sdf = self.name+'_temp.sdf'
        output = pybel.Outputfile('sdf', temp_sdf)
        output.write(self._rdkit_mol)
        self._pybel_mol = readfile('sdf', temp_sdf).next()
        os.remove(temp_sdf)


class g16_log(file_parser):
    def __init__(self, file):
        # default values for thermochemical calculations
        if '.log' not in file:
            raise TypeError('A g16 .log file must be provided')
            
        self.file=file
        self.name=os.path.basename(file)
        
        QH, freq_cutoff, temperature, freq_scale_factor, solv, spc = 'grimme', 100.0, 298.15, 1.0, 'none', False
        conc = GV.atmos / GV.GAS_CONSTANT / temperature

        self.get_termination() and self.getCPU()
        self.get_error() and self.getCPU()
        
        try:
            #pybel molecule object
            self._pybel_mol = next(readfile('g09', file))
        except:
            raise Exception('You log file is too bad to generate a molecule')
        else:
            self.get_charge()
            self.get_route()
            self.get_coords()
            self.freq = pq.getoutData(file)
         
        if self.termination:
            bbe = GV.calc_bbe(file, QH, freq_cutoff, temperature, conc, freq_scale_factor, solv, spc)
            if hasattr(bbe, 'scf_energy'): self._pybel_mol.data['scf_energy'] = bbe.scf_energy
            if hasattr(bbe, 'enthalpy'): self._pybel_mol.data['enthalpy'] = bbe.enthalpy 
            if hasattr(bbe, 'im_freq'): self._pybel_mol.data['im_freq'] = len(bbe.im_freq)
            if hasattr(self, 'CPU'):
                cpu = self.CPU[0]*3600*24 + self.CPU[1]*3600 + self.CPU[2]*60 + self.CPU[3]
                self._pybel_mol.data['CPU'] = '{}s'.format(cpu)
            
        self.get_rdkit_mol()

    def freq_done(self):
        if 'enthalpy' in self._pybel_mol.data:
            return True
        else:
            return False
        
    def get_termination(self):
        with open(self.file) as fh:
            for line in (fh):
                if line.find("Normal termination") > -1: 
                    self.termination = True
                    return True
            self.termination = False

    def get_error(self):
        with open(self.file) as fh:
            for line in (fh):
                if line.find("Error termination") > -1: 
                    self.error = line
                    return True
            self.error = None

    def get_route(self):
        with open(self.file) as fh:
            line = fh.readline()
            start = False
            route = ''
            while line:
                if '*********************' in line and (not start):
                    if 'Gaussian' in fh.readline(): 
                        start = True

                if start:
                    if '%mem' in line:
                        self.mem = line.rstrip().split('=')[1]
                    if '%nprocshared' in line:
                        self.nprocs = line.rstrip().split('=')[1]
                    if '-------------' in line:
                        line = fh.readline()
                        while '---------' not in line:
                            route += line.rstrip().lstrip()
                            line = fh.readline()
                        break

                line = fh.readline()

            self.route = route

    def getCPU(self):
        with open(self.file) as fh:
            for line in (fh):
                if line.find("Job cpu time") > -1:
                    days = int(line.split()[3])
                    hours = int(line.split()[5])
                    mins = int(line.split()[7])
                    secs = float(line.split()[9])
                    
                    self.CPU=[days, hours, mins, secs]
                    break

    def get_charge(self):
        with open(self.file) as fh:
            for line in (fh):
                if line.find("Charge = ") > -1:
                    self._pybel_mol.data['charge'] = int(line.split()[2])
                    self._pybel_mol.data['mult'] = int(line.split()[5].rstrip("\n")) 
                    break

    def get_coords(self):
        with open(self.file) as fh:
            starting = False
            found_coord = False
            for line in (fh):
                if line.find('orientation') > -1: 
                    starting = True
                    self.AtomsNum = []
                    self.AtomsType = []
                    self.Coords = []
                    sep = 0
                    found_coord = False
                if starting:
                    m = re.search('(\d+)\s+(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)',
                                  line)
                    if not m: continue
                    self.AtomsNum.append(int(m.group(2)))
                    self.AtomsType.append(elementID(int(m.group(2))))
                    self.Coords.append([float(m.group(4)), float(m.group(5)), float(m.group(6))])
                    found_coord = True
                if found_coord and line.find('-----------') > -1: 
                    starting = False
                    found_coord = False

    def fix_imag(self, new_com=None, amp=0.5, route=None):
        if not new_com:
            name = self.file.split('.log')[0]
            name += '_QRC'
            new_com = name + '.com'
        else:
            name = new_com.split('.com')[0]

        shift = []

        # Save the original Cartesian coordinates before they are altered
        orig_carts = []
        for atom in range(0,self.freq.NATOMS):
            orig_carts.append([self.freq.CARTESIANS[atom][0], self.freq.CARTESIANS[atom][1], self.freq.CARTESIANS[atom][2]])

        # could get rid of atomic units here, if zpe_rat definition is changed
        for mode in range(0,3*self.freq.NATOMS-6):
            if self.freq.FREQS[mode] < 0.0:
                shift.append(amp)
            else: shift.append(0.0)

            # The starting geometry is displaced along the each normal mode according to the random shift
            for atom in range(0,self.freq.NATOMS):
                for coord in range(0,3):
                    self.freq.CARTESIANS[atom][coord] = self.freq.CARTESIANS[atom][coord] + self.freq.NORMALMODE[mode][atom][coord] * shift[mode]

        if not route:
            route = self.route
            
            if ' opt' not in route:
                route += ' opt'

        self.write_com(new_com=new_com, route=route)

    def write_com(self, new_com=None, route=None, footer=None):
        if not new_com:
            name = self.file.split('.log')[0]
            new_com = name + '.com'
        else:
            name = new_com.split('.com')[0]

        if not route:
            route = self.route

        with open(new_com, 'w') as nc:
            nc.write('%nprocshared='+str(self.nprocs)+'\n%mem='+self.mem+'\n')
            nc.write(route+'\n\n'+name+'\n\n'+str(self.rdkit_mol.GetProp('charge'))+" "+str(self.rdkit_mol.GetProp('mult'))+'\n')
            for atom in range(0,self.freq.NATOMS):
                nc.write('{:>2} {:12.8f} {:12.8f} {:12.8f}\n'.format(self.AtomsType[atom], self.Coords[atom][0], self.Coords[atom][1], self.Coords[atom][2]))

            if footer:
                nc.write(footer)

            nc.write("\n\n")

    def get_NMR(self):
        self.NMR = []
        with open(self.file, 'r') as fh:
            for line in fh:
                if len(self.NMR) == len(self.AtomsNum): break
                m = re.search('Isotropic\s*=\s*(-?\d+\.\d+)', line)
                if not m: continue
                self.NMR.append(float(m.group(1)))
    
    @property
    def rdkit_mol(self):
        return self._rdkit_mol
    
    @property
    def pybel_mol(self):
        return self._pybel_mol
    
        
class sdf(file_parser):
    def __init__(self, file):
        if '.sdf' not in file:
            raise TypeError('A .sdf file must be provided')

        self.name = os.path.splitext(os.path.basename(sdf))[0]
        self._rdkit_mol = Chem.SDMolSupplier(file, removeHs=False, sanitize=False)
        self._pybel_mol = readfile('sdf', file).next()

#FIXME more functions required
class orca_out(file_parser):
    def __init__(self,file):
        if '.out' not in file:
            raise TypeError('A orca .out file must be provided')
    
        self.file = file
        self.name = os.path.basename(file)

        self.get_termination()
        self.get_nprocs()

        if self.termination: self.getCPU()

    def get_termination(self):
        with open(self.file) as fh:
            for line in (fh):
                if line.find("ORCA TERMINATED NORMALLY") > -1: 
                    self.termination = True
                    return True
            self.termination = False

    def get_nprocs(self):
        with open(self.file) as fh:
            for line in (fh):
                if line.find("pal nprocs") > -1:
                    match = re.search(r'nprocs (\d+)', line, re.I|re.M)
                    self.nprocs = int(match.group(1))
                    break

    def getCPU(self):
        with open(self.file) as fh:
            for line in (fh):
                if line.find("TOTAL RUN TIME") > -1:
                    days = int(line.split()[3])
                    hours = int(line.split()[5])
                    mins = int(line.split()[7])
                    secs = float(line.split()[9])
                    
                    self.CPU=[days, hours, mins, secs]
                    break
        

    #FIXME 
    def sp_energy(self):
        sp = GV.sp_energy(self.file)
        return sp 


    #FIXME
    @property
    def rdkit_mol(self):
        return self._rdkit_mol
    
    #FIXME
    @property
    def pybel_mol(self):
        pass

    
    #FIXME
    @rdkit_mol.setter
    def rdkit_mol(self, mol):
        self._rdkit_mol = mol
    

             








