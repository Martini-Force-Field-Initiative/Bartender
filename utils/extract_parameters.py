#!/usr/bin/python

#### EXAMPLE USAGE:
# python3 extract_parameters.py -i inputs/BENZ_cog.itp           --lbl REF
# python3 extract_parameters.py -i inputs/BENZ_Bartender_out.itp --lbl BAR


import sys
import argparse

# Parse the arguments
parser = argparse.ArgumentParser(description='Script to convert an AA-to-CG mapping (.ndx) and CG topology (.itp) to a Bartender input file (.inp)')
parser.add_argument('-i','--itp'  , required=True      , type=str, help='name of the CG topology (as a GROMACS .itp file)')
parser.add_argument('-l','--lbl'  , required=True      , type=str, help='label for output files')

args = parser.parse_args()

INPitp = args.itp
INPlbl = args.lbl

OUTbonds       = f'bonds_{INPlbl}.dat'
OUTconstr = f'constraints_{INPlbl}.dat'
OUTangles      = f'angles_{INPlbl}.dat'


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


def grab_params_and_write_to_file(line,status,outbonds,outangles,outconstr):
    """
    Writes lines according to the Bartender input file format. 

    Parameters
    ----------
    line: string
        A line of the itp file received as input file.
    status: string
        A string which tells in which section of the itp we are. 

    Returns
    --------
    elements[i]: string 
        The equilibrium bond distance or angle and associated force constant. 
    """

    elements = line.split() 

    if len(elements) == 0:
        pass
    elif status == 'constraints':
        outconstr.write('{0:12.3f}\n'.format(float(elements[3])))
    elif status == 'bonds':
        outbonds.write('{0:12.3f} {1:12.3f}\n'.format( float(elements[3]), float(elements[4]) ))
    elif status == 'angles':
        outangles.write('{0:12.3f} {1:12.3f}\n'.format( float(elements[4]), float(elements[5]) ))



with open(INPitp, 'r') as inpfile, open(OUTbonds, 'a') as outbonds, open(OUTangles, 'a') as outangles, open(OUTconstr, 'a') as outconstr:
   lines = inpfile.readlines()
   status = None
   for line in lines:
       if line[0] == '[':
           status = check_itp_line(line)
       elif not line[0] == ';' and not line[0] == '#' and not line == '\n':
           grab_params_and_write_to_file(line, status, outbonds, outangles, outconstr) 

print(f"- DONE! Parameters written to: {OUTbonds}, {OUTangles}, and {OUTconstr}.")

