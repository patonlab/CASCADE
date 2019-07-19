"""
Created on Tue Jul 24

@author Yanfei Guan
"""
from genConf import genConf, postrmsd, energy_filter
from optimize_DB import write_com, run_gau_job, get_mult_from_DB
from file_parsers import g16_log
import os, sys, argparse, warnings
import time
from rdkit import Chem

parser = argparse.ArgumentParser()
parser.add_argument('-i','--isdf', nargs='+', required=True)
parser.add_argument('-nconf', type=int, required=False, 
                    help='number of conformers')
parser.add_argument('-nDFTconf', type=int, required=False, default=1,
                    help='Number of DFT conformers')
parser.add_argument('-rmspreFF', type=float, required=False, 
                    help='rms threshold pre FF optimization')
parser.add_argument('-rmspostFF', type=float, required=False, default=0.2, 
                    help='rms threshold post FF minimization')
parser.add_argument('-rmspostDFT', type=float, required=False, default=0.2,
                    help='rms threshold post DFT optimization')
parser.add_argument('-cutoffFF', type=float, required=False, default=10.0, 
                    help='energy window for FF kcal/mol')
parser.add_argument('-cutoffDFT', type=float, required=False, default=0.015,
                    help='energy window for DFT hartree')    
parser.add_argument('-confdir', type=str, required=False, default='Conformers',
                    help='conformation directory')
parser.add_argument('-DFTcharge', type=int, required=False, default=0,
                    help='charge for DFT calculation')
parser.add_argument('-DFTmult', type=int, required=False, default=1,
                    help='mult for DFT calculation')
parser.add_argument('-DFTlevel', type=str, required=False, default='M062x/Def2TZVP',
                    help='DFT optimization level of theory')
parser.add_argument('-g16root', type=str, required=True,
                    help='Gaussian 16 root')
parser.add_argument('-scratch', type=str, required=True,
                    help='Scratch place')
parser.add_argument('-log', type=str, required=False, default='genConf.log',
                    help='log file to record')
parser.add_argument('-NMR', action='store_true',
                    help='if used, do NMR calculation after optimization')
parser.add_argument('--skip_FF', action='store_true',
                    help='if used, skipFF')

args = parser.parse_args()
n_conf_dict = {}

route = "#"+args.DFTlevel+' opt=(maxcycle=1000) freq=noraman'
if args.NMR:
    NMR_route = '# NMR=GIAO mpw1pw91/6-311+g(d,p) scrf=(smd, solvent=Chloroform)'

if not os.path.isdir(args.confdir): os.mkdir(args.confdir)

DFT_run_conf_dir = args.confdir+'_DFTrun'
if not os.path.isdir(DFT_run_conf_dir): os.mkdir(DFT_run_conf_dir)

conf_log = open(args.log, 'w')
for sdf in args.isdf:
    name = os.path.splitext(os.path.basename(sdf))[0]
    out = os.path.join(args.confdir, name+'_FFConf.sdf')
    n_conf_dict[name] = {}
    
    if os.path.isfile(out):
        conf_log.write("FF conformers for %s were found in sdf file, now DFT optimization\n" % (name))
    elif args.skip_FF:
        mol = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0]
        writer = Chem.SDWriter(out)
        mol.SetProp('_Name', name+'_conf1')
        mol.SetProp('FFEnergies', "0 kcal/mol")
        mol.SetProp('ConfId', '1')
        writer.write(mol)

    else:
        m = Chem.SDMolSupplier(sdf, removeHs=False)[0]
        writer = Chem.SDWriter(out)

        mol,ids,nr = None,None,None
        try:
            mol,ids,nr = genConf(m, args.nconf, args.rmspreFF, args.cutoffFF, args.rmspostFF)        
        except Exception as e:
            msg = "cannot generate FF conformers for %s: %s\n" % (name, str(e))
            conf_log.write(msg)
            continue
        
            
        for en,id in ids:
            mol.SetProp('_Name', name+'_conf'+str(id))
            mol.SetProp('FFEnergies', str(en)+ " kcal/mol")
            mol.SetProp('ConfId', str(id))
            writer.write(mol, confId=id)
            
        writer.close()
        
        conf_log.write("FF conformer generating for %s is done, now DFT optimization\n" % (name))

    try:
        mol_confs = Chem.SDMolSupplier(out, removeHs=False, sanitize=False)
    except:
        continue

    n_conf_dict[name]['FF'] = len(mol_confs)    
    DFT_out = os.path.join(args.confdir, name+'_DFTConf.sdf')
    
    if os.path.isfile(DFT_out):
        mol_DFT_confs = Chem.SDMolSupplier(DFT_out, removeHs=False, sanitize=False)
        n_conf_dict[name]['DFT'] = len(mol_DFT_confs)
        conf_log.write("DFT conformers were found in the sdf file for %s " % (name))
        conf_log.write("%s: %s FF conformers %s DFT conformers\n" % 
                       (name, n_conf_dict[name]['FF'], n_conf_dict[name]['DFT']))
        continue

    mol_confs = [x for x in mol_confs][:args.nDFTconf]
    try:
        mol_DFT = Chem.Mol(mol_confs[0])
    except:
        continue

    DFT_diz = []
    DFT_FF = {}
    DFT_confs = []
    for conf in mol_confs:
        conf_name = conf.GetProp('_Name')
        log_name = os.path.join(DFT_run_conf_dir, conf_name+'.log')

        DFT_log = None
        if os.path.isfile(log_name):
            try:
                DFT_log = g16_log(log_name)
            except Exception as e:
                msg = "Log file %s too bad: %s\n" % (log_name, e)
                conf_log.write(msg)

        skip_opt = False
        if DFT_log and DFT_log.termination and DFT_log.freq_done(): 
            skip_opt = True
            log = DFT_log
        elif DFT_log and DFT_log.error:
            continue
        elif DFT_log:
            try:
                log_conf = DFT_log.rdkit_mol
                for prop in conf.GetPropNames():
                    log_conf.SetProp(prop, conf.GetProp(prop))
                log_conf.SetProp('_Name', conf_name)
                conf = log_conf
            except:
                continue

        if not skip_opt:
            try:
                if args.DFTmult == 999:
                    mult = get_mult_from_DB('BDE_DB.csv', name, 'radical')
                else:
                    mult = args.DFTmult
                write_com(conf, 
                          name=conf_name, 
                          route=route,
                          calcdir=DFT_run_conf_dir,
                          charge=args.DFTcharge,
                          mult=mult)

                run_gau_job(conf.GetProp('_Name'),
                        args.scratch,
                        args.g16root,
                        calcdir=DFT_run_conf_dir,
                        overwrite=True)
            except Exception as e:
                msg = "Cannot optimize %s by DFT method: %s\n" % (conf_name, e)
                conf_log.write(msg)
                continue

            times = 1
            keep_going = True
            break_outer = False

            while times <= 3 and keep_going:        
                try:
                    log = g16_log(os.path.join(DFT_run_conf_dir, conf_name+'.log'))
                except:
                    break_outer = True
                    break

                if log.termination:
                    if 'im_freq' in log.pybel_mol.data and (int(log.pybel_mol.data['im_freq']) > 0):
                        msg = "%s has %s imaginary frequency, fix and restart it\n" %(conf_name, log.pybel_mol.data['im_freq'])
                        conf_log.write(msg)
                        try:
                            log.fix_imag(new_com=os.path.join(DFT_run_conf_dir, conf_name+'.com'))
                            run_gau_job(conf.GetProp('_Name'),
                                    args.scratch,
                                    args.g16root,
                                    calcdir=DFT_run_conf_dir,
                                    overwrite=True)
                            msg = "frequency try times %s\n" % (times)
                        except Exception as e:
                            msg = "Cannot fix the imaginary frequency, skip it\n"
                            conf_log.write(msg)
                            keep_going = False
                    else: keep_going = False
                else: break
                    
                times += 1

            if break_outer:
                continue
    
        if log.termination:
            conf_DFT = log.rdkit_mol
            for prop in conf_DFT.GetPropNames():
                conf.SetProp(prop, conf_DFT.GetProp(prop))

            DFT_confs.append(conf)
            DFT_id = mol_DFT.AddConformer(conf_DFT.GetConformer(), assignId=True)
            DFT_e = conf_DFT.GetProp('scf_energy')
            DFT_diz.append([float(DFT_e), DFT_id])
            DFT_FF[DFT_id] = len(DFT_confs) - 1
        else:
            msg = "DFT Optimization for %s failed\n" % (conf_name)
            conf_log.write(msg)

    mol_confs = DFT_confs

    if len(DFT_diz) > 0:
        mol_DFT_ef, DFT_diz_ef = energy_filter(mol_DFT, DFT_diz, args.cutoffDFT)
        m_DFT, DFT_diz_postrms = postrmsd(mol_DFT_ef, DFT_diz_ef, args.rmspostDFT)
        n_conf_dict[name]['DFT'] = len(DFT_diz_postrms)
        out = os.path.join(args.confdir, name+'_DFTConf.sdf')
        writer = Chem.SDWriter(out)

        for en,id in DFT_diz_postrms:
            if args.NMR:
                name_tmp = mol_confs[DFT_FF[id]].GetProp('_Name') + '_NMR'
                log_name = name_tmp + '.log'
                skip_NMR = False
                nmr_log = None
                if os.path.isfile(log_name):
                    try:
                        nmr_log = g16_log(nmr_log_name)
                    except Exception as e:
                        msg = "Log file %s too bad: %s\n" % (nmr_log_name, e)
                        conf_log.write(msg)

                if nmr_log and nmr_log.termination: 
                    skip_NMR = True
                elif nmr_log and nmr_log.error:
                    continue

                if not skip_NMR:
                    try:
                        if args.DFTmult == 999:
                            mult = get_mult_from_DB('BDE_DB.csv', name, 'radical')
                        else:
                            mult = args.DFTmult

                        write_com(m_DFT, 
                                  name=name_tmp, 
                                  route=NMR_route,
                                  calcdir=DFT_run_conf_dir,
                                  charge=args.DFTcharge,
                                  mult=mult,
                                  confId=id)

                        run_gau_job(name_tmp,
                                args.scratch,
                                args.g16root,
                                calcdir=DFT_run_conf_dir,
                                overwrite=True)
                    except Exception as e:
                        msg = "Cannot get NMR %s by DFT method: %s\n" % (conf_name, e)
                        conf_log.write(msg)
                try:
                    log = g16_log(os.path.join(DFT_run_conf_dir, name_tmp+'.log'))
                except: continue

                if log.termination:
                    log.get_NMR()
                    nmr = ','.join([str(x) for x in log.NMR])

                    for prop in mol_confs[DFT_FF[id]].GetPropNames():
                        m_DFT.SetProp(prop, mol_confs[DFT_FF[id]].GetProp(prop))
                    m_DFT.SetProp('_Name', mol_confs[DFT_FF[id]].GetProp('_Name'))
                    m_DFT.SetProp('NMR', nmr)
                    writer.write(m_DFT, confId=id)
            else:
                for prop in mol_confs[DFT_FF[id]].GetPropNames():
                    m_DFT.SetProp(prop, mol_confs[DFT_FF[id]].GetProp(prop))
                m_DFT.SetProp('_Name', mol_confs[DFT_FF[id]].GetProp('_Name'))
                writer.write(m_DFT, confId=id)
                
        writer.close()    
        
    else:
        msg = "No conformer left for %s after DFT optimization, check that!" % (name)
        conf_log.write(msg)
        n_conf_dict[name]['DFT'] = 0

    conf_log.write("%s: %s FF conformers %s DFT conformers\n" % 
                   (name, n_conf_dict[name]['FF'], n_conf_dict[name]['DFT']))
        
conf_log.close()            
                        
                        

                
            




