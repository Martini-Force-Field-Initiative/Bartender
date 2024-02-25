#!/usr/bin/python

#### EXAMPLE USAGE:
# python write_bartender_inp.py --ndx BENZ_oplsaaTOcg_cgbuilder.ndx         --itp BENZ_cog.itp --out BENZ_bartender.inp  
# python write_bartender_inp.py    -n NDMBI_oplsaaTOcg_cgbuilder.ndx           -i NDMBI.itp       -o NDMBI_bartender.inp  
# python write_bartender_inp.py --ndx TOLU_oplsaaTOcg_cgbuilder_refined.ndx --itp TOLU.itp     --out TOLU_bartender.inp  
# python write_bartender_inp.py    -n PCRE_oplsaaTOcg_cgbuilder.ndx            -i PCRE.itp     --out PCRE_bartender.inp  


import sys
import argparse

# Parse the arguments
parser = argparse.ArgumentParser(description='Script to convert an AA-to-CG mapping (.ndx) and CG topology (.itp) to a Bartender input file (.inp)')
parser.add_argument('-n','--ndx'  , required=True      , type=str, help='name of the AA-to-CG mapping file (as a GROMACS .ndx file)')
parser.add_argument('-i','--itp'  , required=True      , type=str, help='name of the CG topology (as a GROMACS .itp file)')
parser.add_argument('-o','--out'  , default='out.inp'  , type=str, help='name of the outpufile')

args = parser.parse_args()

INPndx = args.ndx
INPitp = args.itp
OUTput = args.out


def check_itp_line(line):
    """
    Checks .itp lines and returns a 'status' accordingly.

    Parameters
    ----------
    line: string
        A line of the itp file received as input file.

    Returns
    --------
    status: string
        A string which tells in which section of the itp we are. 
    """

    if 'moleculetype' in line:
        return 'mtype'
    elif 'atoms' in line:
        return 'atoms'
    elif 'virtual_sitesn' in line:
        return 'virtual_sitesn'
    elif 'bonds' in line:
        return 'bonds'
    elif 'constraints' in line:
        return 'constraints'
    elif 'pairs' in line:
        return 'pairs'
    elif 'angles' in line:
        return 'angles'
    elif 'exclusions' in line:
        return 'exclusions'
    elif 'dihedrals' in line:
        return 'dihedrals'
    else:
        print(f'Something is wrong - this is the troublesome line: ', line)


def write_Bartender_inp(line,status,outfile,init_BONDS,init_ANGLES,init_DIHEDRALS,init_IMPROPERS):
    """
    Writes lines according to the Bartender input file format. 

    Parameters
    ----------
    line: string
        A line of the itp file received as input file.
    status: string
        A string which tells in which section of the itp we are. 
    outfile: string
        Name of output file. 
    init_BONDS/init_ANGLES/init_DIHEDRALS/init_IMPROPERS: bool
        Variable controlling whether the header for the BONDS/ANGLES/DIHEDRALS/IMPROPERS section should be printed (True) or it has already been (False). 

    Returns
    --------
    init_BONDS/init_ANGLES/init_DIHEDRALS/init_IMPROPERS: bool
        Variable controlling whether the header for the BONDS/ANGLES/DIHEDRALS/IMPROPERS section should be printed (True) or it has already been (False). 
    """

    aa_indices = line.split() 

    if len(aa_indices) == 0:
        pass
    elif status == 'constraints' or status == 'bonds':
        if init_BONDS:
            outfile.write(f'{aa_indices[0]},{aa_indices[1]}\n')
        else:
            outfile.write("BONDS\n")
            init_BONDS = True 
            outfile.write(f'{aa_indices[0]},{aa_indices[1]}\n')
    elif status == 'angles':
        if init_ANGLES:
            outfile.write(f'{aa_indices[0]},{aa_indices[1]},{aa_indices[2]}\n')
        else:
            outfile.write("ANGLES\n")
            init_ANGLES = True 
            outfile.write(f'{aa_indices[0]},{aa_indices[1]},{aa_indices[2]}\n')
    elif status == 'dihedrals':
        if int(aa_indices[4]) == 2:
            if init_IMPROPERS:
                outfile.write(f'{aa_indices[0]},{aa_indices[1]},{aa_indices[2]},{aa_indices[3]}\n')
            else:
                outfile.write("IMPROPERS\n")
                init_IMPROPERS = True
                outfile.write(f'{aa_indices[0]},{aa_indices[1]},{aa_indices[2]},{aa_indices[3]}\n')
        elif init_DIHEDRALS:
            outfile.write(f'{aa_indices[0]},{aa_indices[1]},{aa_indices[2]},{aa_indices[3]}\n')
        else:
            outfile.write("DIHEDRALS\n")
            init_DIHEDRALS = True 
            outfile.write(f'{aa_indices[0]},{aa_indices[1]},{aa_indices[2]},{aa_indices[3]}\n')
    else:
        pass

    return init_BONDS, init_ANGLES, init_DIHEDRALS, init_IMPROPERS


def main(INPndx, INPitp, OUTput):
    """
    Reads in the AA-to-CG mapping (.ndx) and CG topology (.itp) and produces a Bartender input file (.inp)

    Parameters
    ----------
    INPndx: string
        Name of the input NDX file (from argparse). 
    INPitp: string
        Name of the input ITP file (from argparse). 
    OUTput: string
        Name of output file (from argparse). 
    """

    print(f"- INPUT files: {INPndx}, {INPitp}")

    with open(OUTput, 'w') as outfile:
        with open(INPndx, 'r') as inpfile:
           lines = inpfile.readlines()
           bead_number = 0
           outfile.write("BEADS\n")
           for line in lines:
               if line[0] == '[':
                   bead_number = bead_number + 1
               else:
                   aa_indices = line.split() 
                   if not len(aa_indices) == 0:
                       aa_indices_for_Bartender = []
                       # First, check if any of the atom indices appears more than once (i.e., has weight < 1)
                       if len(aa_indices) != len(list(set(aa_indices))):
                           print("- INFO - index files contains atoms with weights different from 1")
                           for aa_index in aa_indices:
                               if aa_indices.count(aa_index) > 2:
                                   sys.exit('Hold on, dunno how to handle atoms with weights != from 1 or 1/2.')
                               elif aa_indices.count(aa_index) > 1:
                                   aa_indices_for_Bartender.append(aa_index)
                               else:
                                   aa_indices_for_Bartender.append(f'{str(aa_index)}/2')
                       else:
                           aa_indices_for_Bartender = aa_indices # nothing to do
                       outfile.write(f'{str(bead_number)} ')
                       aa_indices_for_Bartender = list(set(aa_indices_for_Bartender))
                       joined_indices = ",".join(aa_indices_for_Bartender)
                       outfile.write(joined_indices)
                       outfile.write('\n')

        with open(INPitp, 'r') as inpfile:
           lines = inpfile.readlines()
           status = None
           init_BONDS = False; init_ANGLES = False; init_DIHEDRALS = False; init_IMPROPERS = False
           for line in lines:
               if line[0] == '[':
                   status = check_itp_line(line)
               elif not line[0] == ';' and not line[0] == '#' and not line == '\n':
                   init_BONDS, init_ANGLES, init_DIHEDRALS, init_IMPROPERS = write_Bartender_inp(line, status, outfile,
                                                                                                 init_BONDS, init_ANGLES, init_DIHEDRALS, init_IMPROPERS)
    
    print(f"- DONE! Bartender input file written in {OUTput}")


main(INPndx, INPitp, OUTput)

