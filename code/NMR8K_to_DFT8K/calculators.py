import rdkit
from rdkit import Chem
from DB import Rdkit_2_ASE, ASE_2_Rdkit
import os, copy, subprocess, re

class Calculator:
    def __init__(self):
        pass

    def copy(self):
        return copy.deepcopy(self)

    @property
    def path(self):
        return self._path
    
    @path.setter
    def path(self, path):
        self._path = path
        self.sh = self.sh.replace('__JOB__NAME__', self.path)

    @property
    def s_type(self):
        return self._s_type

    @s_type.setter
    def s_type(self, s_type):
        self._s_type = s_type
        if s_type == 'mol':
            self.command = self.command.replace('__MULT__', '1')
        if s_type == 'rad':
            self.command = self.command.replace('__MULT__', '2')

    @property
    def operation(self):
        return self._operation

    def load_batch(self):
        raise NotImplementedError('load_batch not implemented')

    def connect_db(self):
        raise NotImplementedError('connect_db not implemented')

    def disconnect_db(self):
        raise NotImplementedError('disconnect_db not implelemented')

    def launch(self):
        raise NotImplementedError('lauch not implelemented')


class SummitOpt(Calculator):
    """
    """
    def __init__(self):
        super(SummitOpt, self).__init__()

#set your own shell script for queueing software here.

        self.sh = """#!/bin/bash
#SBATCH -J __JOB__NAME__
#SBATCH -p shas
#SBATCH --qos normal
#SBATCH --account=csu8_summit1
#SBATCH -t 23:59:59
#SBATCH -N 1
#SBATCH --ntasks-per-node 24

g16root="/projects/$USER"
export g16root
. $g16root/g16/bsd/g16.profile
Sctchpath="/scratch/summit/$USER/$SLURM_JOB_ID"
mkdir $Sctchpath

export PYTHONPATH=$PYTHONPATH:/home/yanfei@colostate.edu/projects/anaconda3/envs/rdkit/lib/python2.7/site-packages/
"""
        self.command = '~/projects/anaconda3/envs/rdkit/bin/python ~/bin/genDFTConf.py -i __INPUT__ -g16root $g16root -scratch $Sctchpath -nconf 1000 -nDFTconf 1 -DFTmult __MULT__ -DFTlevel M062x/Def2TZVP -NMR'
        self.launch_command = 'sbatch __JOB__NAME__'
        self._operation='opt'

    def load_batch(self, batch):
        """
        """
        self.batch = batch
        batch_sdf = [str(i)+'.sdf' for i in batch]
        self.command = self.command.replace('__INPUT__', ' '.join(batch_sdf))

    def connect_db(self, db):
        """
        """
        self.db = db
        ase_rows = self.db.select(filter=lambda x: x.str_id in self.batch)
        if not os.path.isdir(self.path): os.mkdir(self.path)

        for row in ase_rows:
            ase_m = row.toatoms()
            rdkit_m = ASE_2_Rdkit(ase_m, row.id)

            name = str(row.str_id) + '.sdf'
            writer = Chem.SDWriter(os.path.join(self.path, name))
            writer.write(rdkit_m)
            writer.close()

    def launch(self):
        """
        """
        cwd = os.getcwd()
        os.chdir(self.path)
        sh_file = self.path + '.sh'
        self.launch_command = self.launch_command.replace('__JOB__NAME__', sh_file)
        with open(sh_file, 'w') as sh:
            sh.write(self.sh)
            sh.write(self.command)
        launch = self.launch_command.split()
        jobId = subprocess.check_output(launch)
        jobId = jobId.decode('ascii')

        self.jobId =[int(x) for x in jobId.rstrip().split() if x.isdigit()][0]
        os.chdir(cwd)

    def check_status(self):
        """
        """
        check_command = "sacct --format=\"JobID,Elapsed,time,state\" | grep -v \"bat+\" | grep -v \"ext+\" | grep %i" % (self.jobId)
        status_line = os.popen(check_command).read()
        self.status = status_line.split()[3]

    def on_queue(self):
        """
        """
        check_command = "sacct --format=\"JobID,JobName%60,Elapsed,time,state\" | grep -v \"bat+\" | grep -v \"ext+\" | grep -E \"RUNNING|PENDING\" | grep " + self.path
        status_line = os.popen(check_command).read()

        p = re.search(self.path+'\s', status_line, re.M|re.I)  
        
        if not p: return False
        
        return int(status_line.split()[0])

    def restart(self):
        cwd = os.getcwd()
        os.chdir(self.path)
        sh_file = self.path + '.sh'
        self.launch_command = self.launch_command.replace('__JOB__NAME__', sh_file)
        launch = self.launch_command.split()
        jobId = subprocess.check_output(launch)
        jobId = jobId.decode('ascii')

        self.jobId =[int(x) for x in jobId.rstrip().split() if x.isdigit()][0]
        os.chdir(cwd)

    def gather_data(self):

        conformer_path = os.path.join(self.path, 'Conformers')
        for str_id in self.batch:
            sdf_file = os.path.join(conformer_path, str(str_id) + '_DFTConf.sdf')
            if os.path.isfile(sdf_file):
                try:
                    row = self.db.get(str_id=str_id)
                except:
                    continue
                if row.opt is True: continue
                try:
                    mol = Chem.SDMolSupplier(sdf_file, removeHs=False, sanitize=False)[0]
                except:
                    continue

                try:
                    enthalpy = float(mol.GetProp('enthalpy'))
                except: print(sdf_file)

                nmr = mol.GetProp('NMR').split(',')
                nmr = {i: float(x) for i,x in enumerate(nmr)}
                
                ase_a = Rdkit_2_ASE(mol)
                old_id = row.id
                new_id = self.db.write(ase_a, str_id=str_id, enthalpy=enthalpy, opt=True, data=nmr)
                del self.db[old_id]
                
        
        

        
       
        
        
