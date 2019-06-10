import os, argparse
from DB import NMRDB
import rdkit
from rdkit import Chem
import calculators
import time, logging, glob

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--db', required=True, type=str,
                    help='database to connect')
parser.add_argument('-n', '--num_points', required=False, default=-1, type=int,
                    help='number of data point to get')
parser.add_argument('-b', '--number_of_batches', required=False, default=20, type=int,
                    help='number of batches, that is also how many seperated jobs' \
                         ' running on the queue at the same time')
parser.add_argument('-l', '--level_of_theory', required=False, default='', type=str,
                    help='level of theory to optimize or calculate SPE for data point')
parser.add_argument('-c', '--calculator', required=False, default='Summit_OPT', type=str,
                    help='calculator to use')
parser.add_argument('-r', '--restart', required=False, action='store_true',
                    help='restart calculation without initiating')

args = parser.parse_args()

logging.basicConfig(format='%(asctime)s-%(message)s', level=logging.INFO, filename='nmr.log')

try:
    DB = NMRDB(args.db)
except TypeError as e:
    raise('Cannot initiate database: %s' % e)

computers = []
if args.restart:
    logging.info('Resuming...')

    if args.calculator == 'Summit_OPT': calculator = calculators.SummitOpt()

    mol_batch_dir = glob.glob('mol_batch*')
    
    for batch_dir in mol_batch_dir:
        batch = [int(os.path.splitext(x)[0]) for x in os.listdir(batch_dir) if '.sdf' in x]
        calculator_tmp = calculator.copy()
        calculator_tmp.load_batch(batch)
        calculator_tmp.s_type = 'mol'
        calculator_tmp.path=batch_dir
        computers.append(calculator_tmp)
        
else:
    logging.info('Initiating %i data points from %s\n\n' % (args.num_points, args.db))
    
    calculator = calculators.SummitOpt()
    computers = DB.assign_computer(args.num_points, args.number_of_batches, calculator)

    logging.info('Spliting initial datapoints into %i batches' % args.number_of_batches)

for computer in computers:
    computer.connect_db(DB.db)

    if args.restart:
        jobId = computer.on_queue()
        if jobId:
            computer.jobId = jobId
        else:
            computer.restart()
    else:
        computer.launch()

while True:
    time.sleep(10)
    computers_tmp = []
    for computer in computers:
        computer.check_status()
        if computer.status == 'RUNNING' or (computer.status == 'PENDING'):
            computers_tmp.append(computer)
        if computer.status == 'TIMEOUT':
            computer.restart()
            computers_tmp.append(computer)
        if computer.status == 'Failed':
            logging.info('%s failed, please check' % computer.path)
        if computer.status == 'COMPLETED':
            logging.info('%s completed' % computer.path)
    
        computer.gather_data()
        logging.info('%s status: %s' % (computer.path, computer.status))
    
    computers = computers_tmp
    space = args.number_of_batches - len(computers)
    
    logging.info('%i space were found on the queue' % space)

    if space > 0:
        
        n_nodes = len(glob.glob('*_batch*'))
        n_mols = len([r for r in DB.db.select(filter=lambda x: x.opt is not False)])
        batch_size = int(n_mols/n_nodes)+5

        computers_new = DB.assign_computer(batch_size*space, space, calculator)
        for computer in computers_new:
            if computer.s_type == 'mol':
                computer.connect_db(DB.BDE_DB)
            elif computer.s_type == 'rad':
                computer.connect_db(DB.radical_DB)
    
            computer.launch()
            computers.append(computer)


        




 




