#!/usr/bin/python
from __future__ import print_function, absolute_import

#######################################################################
# Molecular conformer generator in progress
# genConf.py -isdf file_input.sdf -osdf file_output.sdf
# -n number_of_conformers (optional, if not specified is based
# on the nomber of rotable bonds) -rtpre rms_threshold_pre_opt(optional)
# -rtpost rms_threshold_post_opt(optional) -e energy_window (optional, Kcal/mol)
# -t number_of_threads (if not specify 1)
#######################################################################

## known issues / to-do list
## logging doesn't work properly
## need to print out rotatable bond atoms for reference
## comparison of printing vs. old full monte - more options
## comparison of printing vs. GoodVibes
## tabulation of results in addition to sdf output file
## pythonify where possible
## currently RMS but torsion_list would be better?

from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent import futures
import argparse, logging, os, sys, time, copy
logger = logging.getLogger(__name__)

# PHYSICAL CONSTANTS
GAS_CONSTANT, PLANCK_CONSTANT, BOLTZMANN_CONSTANT, SPEED_OF_LIGHT, AVOGADRO_CONSTANT, AMU_to_KG, atmos = 8.3144621, 6.62606957e-34, 1.3806488e-23, 2.99792458e10, 6.0221415e23, 1.66053886E-27, 101.325
# UNIT CONVERSION
j_to_au = 4.184 * 627.509541 * 1000.0

# version number
__version__ = "1.0.1"

# Formatted output to command line and log file
class Logger:
    # Designated initializer
    def __init__(self, filein, ext):
        # check to see if already exists
        if os.path.exists(filein+"."+ext):
            var = input("\no  Logger file %s already exists! OK to proceed? (Y/N) " % (filein+"."+ext))
            if var.lower().strip() == "n" or var.lower().strip() == "no":
                logger.error("\n   OK. Exiting gracefully ...\n"); sys.exit(1)
        # Create the log file at the input path
        self.log = open(filein+"."+ext, 'w' )

    # Write a message to the log
    def Write(self, message):
        # Print the message
        logger.info(message)
        # Write to log
        self.log.write(message + "\n")

    # Write a message only to the log and not to the terminal
    def Writeonlyfile(self, message):
        # Write to log
        self.log.write(message)

    # Write a fatal error, finalize and terminate the program
    def Fatal(self, message):
        # Print the message
        logger.error(message+"\n")
        # Write to log
        self.log.write(message + "\n")
        # Finalize the log
        self.Finalize()
        # End the program
        sys.exit(1)

    # Finalize the log file
    def Finalize(self):
        self.log.close()

dashedline = "   ------------------------------------------------------------------------------------------------------------------"
emptyline = "   |                                                                                                                     |"
normaltermination = "\n   -----------------       N   O   R   M   A   L      T   E   R   M   I   N   A   T   I   O   N      ----------------\n"
leftcol=97
rightcol=12

asciiArt = "     ___       ___                                    ___          ___          ___                   ___ \n    /  /\\     /__/\\                                  /__/\\        /  /\\        /__/\\         ___     /  /\\\n   /  /:/_    \\  \\:\\                                |  |::\\      /  /::\\       \\  \\:\\       /  /\\   /  /:/_\n  /  /:/ /\\    \\  \\:\\   ___     ___  ___     ___    |  |:|:\\    /  /:/\\:\\       \\  \\:\\     /  /:/  /  /:/ /\\\n /  /:/ /:/___  \\  \\:\\ /__/\\   /  /\\/__/\\   /  /\\ __|__|:|\\:\\  /  /:/  \\:\\  _____\\__\\:\\   /  /:/  /  /:/ /:/_\n/__/:/ /://__/\\  \\__\\:\\\\  \\:\\ /  /:/\\  \\:\\ /  /://__/::::| \\:\\/__/:/ \\__\\:\\/__/::::::::\\ /  /::\\ /__/:/ /:/ /\\\n\\  \\:\\/:/ \\  \\:\\ /  /:/ \\  \\:\\  /:/  \\  \\:\\  /:/ \\  \\:\\~~\\__\\/\\  \\:\\ /  /:/\\  \\:\\~~\\~~\\//__/:/\\:\\\\  \\:\\/:/ /:/\n \\  \\::/   \\  \\:\\  /:/   \\  \\:\\/:/    \\  \\:\\/:/   \\  \\:\\       \\  \\:\\  /:/  \\  \\:\\  ~~~ \\__\\/  \\:\\\\  \\::/ /:/\n  \\  \\:\\    \\  \\:\\/:/     \\  \\::/      \\  \\::/     \\  \\:\\       \\  \\:\\/:/    \\  \\:\\          \\  \\:\\\\  \\:\\/:/\n   \\  \\:\\    \\  \\::/       \\__\\/        \\__\\/       \\  \\:\\       \\  \\::/      \\  \\:\\          \\__\\/ \\  \\::/\n    \\__\\/     \\__\\/                                  \\__\\/        \\__\\/        \\__\\/                 \\__\\/\n  "

# algorithm to generate nc conformations
def genConf(m, nc, rms, efilter, rmspost):
    nr = int(AllChem.CalcNumRotatableBonds(m))
    #m = Chem.AddHs(m)
    Chem.AssignAtomChiralTagsFromStructure(m, replaceExistingTags=True)
    if not nc: nc = 3**nr

    print(dashedline+"\n   |    "+("FULL_MONTE search").ljust(leftcol)+("|").rjust(rightcol))
    #log.Write("   | o  "+("COMP: "+str(Params.COMP)+" degrees").ljust(leftcol)+("|").rjust(rightcol))
    #log.Write("   | o  "+("LEVL: "+str(Params.LEVL)+" force field").ljust(leftcol)+("|").rjust(rightcol))
    #log.Write("   | o  "+("DEMX: "+str(Params.DEMX)+" kcal/mol").ljust(leftcol)+("|").rjust(rightcol))
    print("   | o  "+("EWIN: "+str(efilter)+" kcal/mol").ljust(leftcol)+("|").rjust(rightcol))
    print("   | o  "+("MCNV: "+str(nr)+" ROTATABLE BONDS").ljust(leftcol)+("|").rjust(rightcol))
    #log.Write("   |    "+torstring.ljust(leftcol)+("|").rjust(rightcol))
    print("   | o  "+("STEP: "+str(nc)+" (ESTIMATED CONFORMER SPACE: "+str(nr**3)+")").ljust(leftcol)+("|").rjust(rightcol))
    print(dashedline+"\n")

    if not rms: rms = -1
    ids=AllChem.EmbedMultipleConfs(m, numConfs=nc)


    if len(ids)== 0:
        ids = m.AddConformer(m.GetConformer, assignID=True)

    diz = []
    diz2 = []
    diz3 = []
    for id in ids:
        prop = AllChem.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94s")
        ff = AllChem.MMFFGetMoleculeForceField(m, prop, confId=id)
        ff.Minimize()
        en = float(ff.CalcEnergy())
        econf = (en, id)
        diz.append(econf)

    if efilter != "Y":
        n, diz2 = energy_filter(m, diz, efilter)
    else:
        n = m
        diz2 = diz

    if rmspost != None and n.GetNumConformers() > 1:
        o, diz3 = postrmsd(n, diz2, rmspost)
    else:
        o = n
        diz3 = diz2

    return o, diz3, nr

# filter conformers based on relative energy
def energy_filter(m, diz, efilter):
    print("o  FILTERING CONFORMERS BY ENERGY CUTOFF: "+str(efilter)+" kcal/mol")
    diz.sort()
    mini = float(diz[0][0])
    sup = mini + efilter
    n = Chem.Mol(m)
    n.RemoveAllConformers()
    n.AddConformer(m.GetConformer(int(diz[0][1])))
    nid = []
    ener = []
    nid.append(int(diz[0][1]))
    ener.append(float(diz[0][0])-mini)
    del diz[0]
    for x,y in diz:
        if x <= sup:
            #print("   KEEP - Erel:", x-mini)
            n.AddConformer(m.GetConformer(int(y)))
            nid.append(int(y))
            ener.append(float(x-mini))
        else:
            #print("   REMOVE - Erel:", x-mini)
            break
    diz2 = list(zip(ener, nid))
    print("   KEEPING "+str(len(ener))+" CONFORMERS")
    return n, diz2

# filter conformers based on geometric RMS
def postrmsd(n, diz2, rmspost):
    print("o  FILTERING CONFORMERS BY RMS: "+str(rmspost))
    diz2.sort(key=lambda x: x[0])
    o = Chem.Mol(n)
    confidlist = [diz2[0][1]]
    enval = [diz2[0][0]]
    nh = Chem.RemoveHs(n)
    del diz2[0]
    for z,w in diz2:
        confid = int(w)
        p=0
        for conf2id in confidlist:
            #print(confid, conf2id)
            rmsd = AllChem.GetBestRMS(nh, nh, prbId=confid, refId=conf2id)
            if rmsd < rmspost:
                p=p+1
                break
        if p == 0:
            confidlist.append(int(confid))
            enval.append(float(z))
    diz3 = list(zip(enval, confidlist))
    print("   KEEPING "+str(len(enval))+" CONFORMERS")
    return o, diz3

# conformational search / handles parallel threads if more than one structure is defined
def csearch(log, suppl, writer, args):
    with futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        jobs = []
        nm = []
        numMol=0
        for mol in suppl:
            numMol = numMol+1
            if mol is not None:
                nm.append(mol.GetProp('_Name'))
                job = executor.submit(genConf, mol, args.nconf, args.rmspre, args.cutoff, args.rmspost)
                jobs.append(job)
            else:
                log.Write("ERROR: Impossible to generate conformers for molecule number " + str(numMol))
        widgets = ["Generating conformations; ", progressbar.Percentage(), " ", progressbar.ETA(), " ", progressbar.Bar()]

        # If there are several structures use progressbar
        if len(jobs) > 1:
            for job in futures.as_completed(jobs):
                mol,ids,nr = job.result()

        for j in range(0, len(jobs)):
            mol,ids,nr = jobs[j].result()
            name = nm[j]
            for en,id in ids:
                mol.SetProp('_Name', name)

                if args.printproperty == True:
                    mol.SetProp('ConfId', str(id))
                    mol.SetProp('ConfEnergies', str(en)+ " kcal/mol")
                    mol.SetProp('Rotatable Bonds Number', str(nr))
                writer.write(mol, confId=id)
    writer.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Molecular conformer generator')
    parser.add_argument('-isdf', required=True, help='sdf input file')
    parser.add_argument('-osdf', required=True, help='sdf output file')
    parser.add_argument('-nconf', type=int, required=False, help='number of conformers')
    parser.add_argument('-rmspre', type=float, required=False, help='rms threshold pre optimization')
    parser.add_argument('-rmspost', type=float, required=False, default=0.2, help='rms threshold post minimization')
    parser.add_argument('-cutoff', type=float, required=False, default=10.0, help='energy window')
    parser.add_argument('-printproperty', action='store_true', default=True, help='Print molecule properties (energy and rotable bond number)')
    parser.add_argument('-threads', type=int, required=False, default=1, help='number of threads')
    args = parser.parse_args()

    # Check that the input structure exists and has the correct format
    if os.path.exists(args.isdf) and os.path.splitext(args.isdf)[1] in [".sdf"]:
        inp = args.isdf; out = args.osdf
        filename = os.path.splitext(inp)[0]

    # Define input and outputs
    log = Logger(filename,"cs"); writer = Chem.SDWriter(out); suppl = Chem.SDMolSupplier(inp, removeHs=False)
    log.Write("\no  Extracting structure from "+args.isdf+" ...")

    # conformational search
    search = csearch(log, suppl, writer, args)

    end = time.strftime("%H:%M:%S", time.localtime())
    log.Write(asciiArt+end); log.Write(normaltermination); log.Finalize()
