# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 16:41:03 2018

@author: Yanfei-PC
"""

import os, sys, subprocess, warnings, re, shutil
from rdkit import Chem
import pybel
from file_parsers import g16_log
import pandas as pd

def write_com(mol,
              name='new',
              route='',
              footer='',
              calcdir='.', 
              charge=0,
              mult=1,
              memory='96GB',
              nprocs=24,
              confId=None
              ):
    if not isinstance(route, basestring): raise TypeError("route must be a string")
    
    with open(os.path.join(calcdir, name+'.com'), 'w') as fileout:
        fileout.write("%mem="+memory+"\n")
        fileout.write("%nprocshared="+str(nprocs)+"\n")
    
        fileout.write(route+"\n\n")
        fileout.write(name+"\n\n")
        fileout.write("%d %d\n" % (charge, mult))
    
        for j in range(0, mol.GetNumAtoms()):
            if confId:
                atom, pos = mol.GetAtoms()[j], mol.GetConformer(confId).GetAtomPosition(j)
            else:
                atom, pos = mol.GetAtoms()[j], mol.GetConformer().GetAtomPosition(j)
            fileout.write('{:<4} {:10.6f} {:10.6f} {:10.6f}'.format(atom.GetSymbol(), pos.x, pos.y, pos.z))
            fileout.write("\n")
        fileout.write("\n")

def write_orca_inp(mol,
                  name='file_parser_auto_in',
                  route='DLPNO-CCSD(T) Extrapolate(2/3,cc)',
                  charge=0,
                  mult=1,
                  memory='4GB',
                  nprocs=12,
                  ):
    if not isinstance(route, basestring): raise TypeError("route must be a string")
    with open(name+'.inp', 'w') as fileout:
        fileout.write("%maxcore "+str(int(memory.split('GB')[0])*1000)+"\n")
        fileout.write("% pal nprocs "+str(nprocs)+" end\n")
        fileout.write(route+"\n")
        fileout.write("%scf maxiter 500\n   end\n% mdci\n")
        fileout.write("    TCutPairs 1e-4\n    TCutPNO 3.33e-7\n    TCutMKN 1e-3\n")
        fileout.write("  Density None\nend\n")
        fileout.write("% output\n  printlevel mini\n  print[ P_SCFInfo ] 1\n")
        fileout.write("  print[ P_SCFIterInfo ] 1\n  print[ P_OrbEn ] 0\n  print[ P_Cartesian ] 0\nend\n")
        fileout.write("% elprop\n  Dipole False\nend\n")

        fileout.write('*xyz '+str(charge)+" "+str(mult)+"\n")
        
        for j in range(0, mol.GetNumAtoms()):
            atom, pos = mol.GetAtoms()[j], mol.GetConformer().GetAtomPosition(j)
            fileout.write('{:<4} {:10.6f} {:10.6f} {:10.6f}'.format(atom.GetSymbol(), pos.x, pos.y, pos.z))
            fileout.write("\n")
        fileout.write("*\n")
        
def sdf_to_com(sdf,calcdir,route,charge=0,mult=1):
    mol = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0]
    name = os.path.splitext(os.path.basename(sdf))[0]
    write_com(mol, name=name, route=route, calcdir=calcdir, charge=charge, mult=mult)
    return name

def run_gau_job(name,
            scratch,
            g16_root,
            calcdir='.',
            overwrite=False,
            ):

    if not isinstance(name, basestring): raise TypeError("filename must be a string")

    cwd = os.getcwd()
    com = os.path.join(cwd, calcdir, name+'.com')
    log = os.path.join(cwd, calcdir, name+'.log')    
    if os.path.isfile(log):
        if not overwrite:
            return        
    try:
        os.chdir(scratch)
    except:
        raise Exception('Cannot change into %s' % (scratch))

    grun = g16_root+"/g16/g16 < %s > %s" % (com, log)
    #grun = "/projects/yanfei@colostate.edu/g16/g16 < %s > %s" % (com, log)
    p = subprocess.Popen(grun, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p.wait()

    try:
        os.chdir(cwd)
    except:
        raise Exception('Cannot change back into working directory')

def run_orca_job(name,
                scratch,
                orca_root,
                calcdir='.',
                overwrite=False,):
    if not isinstance(name, basestring): raise TypeError("filename must be a string")

    cwd = os.getcwd()

    bare_inp = name+'.inp'
    bare_out = name+'.out'

    inp = os.path.join(calcdir, name+'.inp')
    out = os.path.join(calcdir, name+'.out')

    if os.path.isfile(out):
        if not overwrite:
            return

    shutil.copyfile(inp, os.path.join(scratch, bare_inp)) 

    try:
        os.chdir(scratch)
    except:
        raise Exception('Cannot change into %s' % (scratch))

    orun = orca_root+'/orca '+bare_inp+' > '+bare_out
    #orun = "/home/yanfei@colostate.edu/rpaton/orca_4_0_1_2_linux_x86-64/orca "+inp+' > '+out
    p = subprocess.Popen(orun, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p.wait()

    shutil.copyfile(bare_out, os.path.join(cwd, calcdir, bare_out))

    try:
        os.chdir(cwd)
    except:
        raise Exception('Cannot change back into working directory')    

def get_mult_from_DB (Db, Id, m_type):
    mult = 1 
    with open (Db, 'r') as csv:
        df = pd.read_csv(csv, dtype={'ID':'str', 'C_Struct':'str', 'R1_Struct':'str', 'R2_Struct':'str'})
        
        if m_type == 'radical':
            match_df = df.loc[(df['R1_Struct']==Id) | (df['R2_Struct']==Id)].iloc[0]

            if match_df['B_Type'] == 'C#C':
                if pd.isnull(match_df['Radical1']) or pd.isnull(match_df['Radical2']):
                    mult = 3
                else: mult = 2
            elif match_df['B_Type'] == 'C=C':
                if pd.isnull(match_df['Radical1']) or pd.isnull(match_df['Radical2']):
                    mult = 5
                else: mult = 3
            elif match_df['B_Type'] == 'C-C':
                if pd.isnull(match_df['Radical1']) or pd.isnull(match_df['Radical2']):
                    mult = 3
                else: mult = 2
            else: mult = 2
    return mult    

if __name__ == "__main__":
    g16_root = sys.argv[1]
    scratch = sys.argv[2]
        
    compounds = [f for f in os.listdir('compound') if os.path.isfile(os.path.join('compound', f))]
    radicals = [f for f in os.listdir('radical') if os.path.isfile(os.path.join('radical', f))]
    
    opt_compounds = [f for f in os.listdir('opt_compound') if os.path.isfile(os.path.join('opt_compound', f))]
    opt_radicals = [f for f in os.listdir('opt_radical') if os.path.isfile(os.path.join('opt_radical', f))]
    
    u_compounds = [x for x in compounds if x not in opt_compounds]
    u_radicals = [x for x in radicals if x not in opt_radicals]
    
    route = '#UB3LYP/6-31G(d) opt'
    
    
    for compound in u_compounds:
        sdf = os.path.join('compound', compound)
        com = sdf_to_com(sdf, 'run_compound', route)
        run_gau_job(com, scratch, g16_root, 'run_compound')
        name = os.path.basename(com)
        log = g16_log(os.path.join('run_compound', name+'.log'))
        
        if log.termination:
            output = pybel.Outputfile('sdf', os.path.join('opt_compound', name+'.sdf'))
            output.write(log.mol)
        else:
            msg = "compounds %s did not finish normally" % (name)
            warnings.warn(msg)
            
    
    for radical in u_radicals:
        sdf = os.path.join('radical', radical)
        name = sdf_to_com(sdf, 'run_radical', route, mult=2)
        run_gau_job(name, scratch, g16_root, 'run_radical')
        log = g16_log(os.path.join('run_radical', name+'.log'))
        
        if log.termination:
            output = pybel.Outputfile('sdf', os.path.join('opt_radical', name+'.sdf'))
            output.write(log.mol)
        else:
            msg = "radical %s did not finish normally" % (name)
            warnings.warn(msg)

   
