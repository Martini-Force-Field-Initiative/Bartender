#!/usr/bin/env python3

import os
import argparse


# Parse the arguments
parser = argparse.ArgumentParser(description='Produce Bartender itps that can be used in simulations.')
parser.add_argument('-l', '--list'      , required=True, type=str, help='list of all the molecules')
parser.add_argument('-r', '--ref-folder', required=True, type=str, help='folder containing the reference .itp files')
parser.add_argument('-b', '--bar-folder', required=True, type=str, help='folder containing the bartender gmx_out.itp files')

args = parser.parse_args()

LIST       = args.list        # 'list_2022ATS_90mol'
REF_FOLDER = args.ref_folder  # 'reference-COG-itps-90mol' 
BAR_FOLDER = args.bar_folder  # 'dataset-sm30-2022ATS-90mol'

VERBOSE    = False
CWD        = os.getcwd() # "~"


def check_itp_line(line):
    """
    Checks which Gromacs section of the itp file we are at, returning a 'status' accordingly.

    Parameters
    ----------
    line: string
        A line of the itp file received as input file.

    Returns
    --------
    status: string
        A string which tells in which section of the itp we are. 
    """
    if 'defaults' in line:
        return 'defaults'
    elif 'atomtypes' in line:
        return 'atomtypes'
    elif 'system' in line:
        return 'system'
    elif 'molecules' in line:
        return 'molecules'
    elif 'moleculetype' in line:
        return 'mtype'
    elif 'atoms' in line:
        return 'atoms'
    elif 'virtual_sitesn' in line:
        return 'virtual_sitesn'
    elif 'virtual_sites2' in line:
        return 'virtual_sites2'
    elif 'virtual_sites3' in line:
        return 'virtual_sites3'
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
    elif 'link' in line:
        return 'link'
    else:
        sys.exit('! ERROR ! Something is wrong. Just found the following line in the input itp file: ' + line)



def merge_bar_and_ref_itps(mol_string, out_log):
    """
    """
    if os.path.exists(os.path.join( REF_FOLDER, f'{mol_string}/{mol_string}_cog.itp')) and \
       os.path.exists(os.path.join( BAR_FOLDER, f'{mol_string}/gmx_out.itp'         )):

        with open( os.path.join( REF_FOLDER, f'{mol_string}/{mol_string}_cog.itp'), 'r') as ref_itp, \
             open( os.path.join( BAR_FOLDER, f'{mol_string}/gmx_out.itp'         ), 'r') as bar_itp, \
             open( os.path.join( BAR_FOLDER, f'{mol_string}/{mol_string}_bar.itp'), 'w') as out_itp:


                # PART 0 - first 2 lines of Bartender itp
                bar_lines = bar_itp.readlines()
                for bar_line in bar_lines[:2]:
                    out_itp.write(bar_line)

                # PART 1 - [mtype], [atoms] from reference itp
                status_ref = 'start'
                for ref_line in ref_itp:
                    if ref_line[0] == '[':
                        status_ref = check_itp_line(ref_line)
                        if status_ref in ['start','atoms','mtype']:
                            out_itp.write(ref_line)
                        continue
                    if status_ref in ['start','atoms','mtype']:
                        out_itp.write(ref_line)
                    else:
                        pass
                        if VERBOSE:
                            print(f'Skipping line {ref_line}')

                # PART 2 - [bonds], [constraints], [angles], [dihedrals] of Barteneder itp
                bar_itp.seek(0) # rewind 'bar_itp' file
                status_bar = None 
                for bar_line in bar_itp:
                    if bar_line[0] == '[':
                        status_bar = check_itp_line(bar_line)
                        if status_bar in ['bonds','constraints','angles','dihedrals']:
                            out_itp.write(bar_line)
                        continue
                    if status_bar in ['bonds','constraints','angles','dihedrals']:
                        out_itp.write(bar_line)
                    else:
                        pass
                        if not "Bartender" in bar_line and VERBOSE:
                            print(f'Skipping line {bar_line}')

                # PART 3 - [virtual_sitesn], [exclusions], etc. from reference itp
                ref_itp.seek(0) # rewind 'ref_itp' file
                status_ref = None 
                for ref_line in ref_itp:
                    if ref_line[0] == '[':
                        status_ref = check_itp_line(ref_line)
                        if status_ref in ['virtual_sitesn','virtual_sites2','virtual_sites3','exclusions']:
                            out_itp.write(ref_line)
                        continue
                    if status_ref in ['virtual_sitesn','virtual_sites2','virtual_sites3','exclusions']:
                        out_itp.write(ref_line)
                    else:
                        pass
                        if VERBOSE:
                            print(f'Skipping line {ref_line}')

        out_log.write(f'File bar.itp written for {mol_string}.\n')

    else:
        out_log.write('Either file {0} or file {1} do not exist. Skipping molecule {2}\n'.format(
                      os.path.join( REF_FOLDER, f'{mol_string}/{mol_string}_cog.itp'),
                      os.path.join( BAR_FOLDER, f'{mol_string}/gmx_out.itp'         ),
                      mol_string))


with open( os.path.join( CWD, LIST), 'r') as listfile, open( os.path.join( CWD, '3_itp_merging.log'), 'w') as out_log:
    list_of_mols = listfile.readlines()
    for mol in list_of_mols:
        mol_string = mol.split()[0]
        print(f'{mol_string}')
        merge_bar_and_ref_itps(mol_string, out_log)

