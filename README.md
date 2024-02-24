![Bartender](graphics/logofinal_trazo_33.png)

# Automated bonded parameters for Martini 3

If you use or wish to cite Bartender, please refer to the following manuscript:

- G.P. Pereira, R. Alessandri, M. Domínguez, R. Araya-Osorio, L. Grünewald,
  L. Borges-Araújo, S. Wu, S.J. Marrink, P.C.T. Souza, R. Mera-Adasme.
  "Bartender: Martini 3 Bonded Terms via Quantum Mechanics-based Molecular Dynamics",
  *Chemrxiv* **2024**, [10.a26434/chemrxiv-2024-62kh1](https://doi.org/10.26434/chemrxiv-2024-62kh1)

Bartender uses [xtb](https://github.com/grimme-lab/xtb). Please cite the corresponding [methods](https://xtb-docs.readthedocs.io/en/latest/xtbrelatedrefs.html#methods)

## Resources

Take a look at the [wiki](https://github.com/Martini-Force-Field-Initiative/Bartender/wiki) for
more complete documentation.

For questions, ideas, or engaging the community,
check out the [Discussions](https://github.com/Martini-Force-Field-Initiative/Bartender/discussions)

If you find a problem with Bartender or its documentation, please open a [Github Issue](https://github.com/Martini-Force-Field-Initiative/Bartender/issues)

## Install 

First, download the **tgz** file (the binary distribution)
of the latest [release](https://github.com/Martini-Force-Field-Initiative/Bartender/releases)
(Linux only) and uncompress it to a directory.

The installation instructions can be found in the file
INSTALL in the uncompressed tgz.

It should be possible to compile Bartender 
for other operating systems/architectures, but that 
has not been attempted or is supported by the Bartender team.


## Bartender use

Assuming the Bartender and xtb excecutables are in the PATH, you
can use Bartender with the command:

```
bartender [flags] Geometry.xyz BartenderInput.inp
```

The geometry file can also be in PDB or GRO format. The Bartender input format
is a simple way of specifying the bonded parameters to be obtained. You can find a
sample in the directory of the distribution, under _samples/_.

As output, Bartender will produce several files.

1. The most important one is called **gmx_out.itp** and contains the parameters for
the bonded terms, in Gromacs itp format.

1. Before actually running the parametrization, Bartender will produce a file called **Beads.pdb**
This file contains the atomistic structure given as an input, with the bead to which each atom belongs
in the residue ID field of each atom's PDB line. The b-factor field will also reflect the bead number of
each atom, interpolated from 0 to 99. Thus, if the structure is colored by b-factor, the atoms belonging to
the same bead will appear in the same color. If an atom belongs to more than one bead, the numbers assigned
will correspond to one of them. Virtual sites will be marked with Uranium atoms.

1. Plots for each fitting showing both the fitted function and the values from the MD distribution will be written
as PNG files. Each PNG's name contains information of the term type (dihedral, angle, etc) the function fitted (Hooke, Simple Periodic, etc) and the beads involved in the fit.

1. The same data in the previous point is also written as Grace XVG files (which are basically plain text and can be easily plotted with other tools), in a directory called **xvg**

### Flags

The optional flags control the way Bartender behaves. Sensible defaults have been prepared so, in
most cases, no flags are needed. Use

```
bartender -help
```

To get all the flags available and their use. We document here some of the most common flags:

*  `-charge` _int_ the total charge of the system, in a.u. (default 0)
*  `-method` _string_ The method employed in the semiempirical simulation. Valid options are gfn0, gfn1,gfn2 and gfnff (default "gfn2")
*  `-time` _int_ The total simulation time for the QM MD, in ps. If a number <0 is given, the MD will not be performed, and a previous trajectory will be used (default 10000, i.e. 10 ns).
*  `-verbose` _int_  Sets the level of verbosity
*  `-dcdSave` _filename.dcd_ Saves the trajectory produced by xtb in the more compact DCD format
*  `-owntraj` _filename_ Reads a trajectory (DCD, XTC, multiPDB or multiXYZ, identified from the file extension) instead of performing an xtb simulation.
*  `-restart` Restart from a previouls xtb simulation. The xtb.trj file must be present, and it's your responsibility to ensure that the same method is used in both runs. The already run time from the previous trajectory (rounded to ps) will be discounted from the total simulation time requested. Note that this is not a "proper" restart, in that, velocities are not taken from the previous run, only coordinates.


## Input file format

The input file has several sections. Each section is preceded by a "section header" in its own line. All section headers are written in capitalized letters. Each section comprises the lines between its section header and the next section header (excluding both). All sections need to be present, even if they contain 0 lines. The sections need to appear in the order specified here. 
Lines starting with the character '\#' will be taken as comments and ignored by Bartender. 

Follows a description of the format of each section.

### Section "Beads". Section header: BEADS
Each line of the section describes a Martini Bead in terms of the atoms that it represents. I has the form:

* BeadNumber atom1,...,atomN 

Where BeadNumber is an integer, containing the index of the bead (the same information contained in the position of the line within the section, so the first line must correspond to the BeadNumber 1, and so forth). atom1,atom2,...,atomN is a list of the 1-based indexes of the N atoms that form the bead. Alternatively, if only 1/M of the atom L wants to be assigned to a bead (needed, for instance, when an atom is split between 2 beads), the line can be written as:

* BeadNumber atom1,...,atomL/M,...,atomN 

Where M is an integer number. It is possible to do the previous with as many atoms as needed, and with different denominators for each atom.

### Section "Virtual sites". Section header: VSITES
Each line corresponds to a virtual site, a "fake" bead with a position that is determined by the position of other beads according to one of several available functions. The line has the form:

* BeadNumber  Bead1,...,BeadN, Parameter1,...,ParameterM, GromacsFunctionNumber

Where BeadNumber is the index of the virtual site. If there are L 'actual' Martini Beads, the first virtual site has an index of L+1. The numbering needs to agree with the position of the line in the section. Parameter1,..ParameterM is a list of M floating-point numbers corresponding to the parameters required by the function to obtain the position of the virtual site from the position of the Martini beads. GromacsFunctionNumber is the Gromacs number for the aforementioned function.


### Sections for bonded terms. Section headers BONDS, ANGLES, DIHEDRALS, IMPROPERS

In these sections the bonded terms to be parametrized are listed.
Each of these sections contains a line per bonded term, with the beads involved in the term,
separated by spaces. A reduced sample of these sections follows:

```
BONDS
1,2
2,3
ANGLES
1,4,5 
1,4,6
DIHEDRALS
1,2,4,5
IMPROPERS  
3,2,1,4
```

### "Goodenough" File

Though, in most cases, the default values are sufficient, The option ``` -goodenoughfile FILE``` allows the user to supply a file to indicate Bartender what are the maximum values acceptable for fitting RMSDs, and what range of values are acceptable for force constants, for each bond term type.

The file must be formatted as follows:

```
dihe k v1 v2 rmsd v3 v4
angles k v1 v2 rmsd v3
impropers k v1 v2 rmsd v3
bonds k v1 v2 rmsd v3
```

Where v1,v2,v3 and v4 are numbers (they may be floating point).

Note that the "dihe" line containes 2 numbers for the rmsd. 
These are both maximum values, for the simple-periodic and one for the
Ryckaert-Bellemans function.

## Contributing


Please open a Github issue for bugfixes/problems.
Similarly, please open a Github issue for minor code problems like typo fixes.

If you wish to make a larger contribution to the program, please tell us about it in the [Discussions](https://github.com/Martini-Force-Field-Initiative/Bartender/discussions)
or contact the lead developer at rmeraa+bartender@academicos.uta.cl



### Licenses

Bartender is open source, under LGPLv3.

The Bartender Logo is Copyright (c) Rocio Araya employing the [Dancing Script](https://github.com/impallari/DancingScript) font by Impallary Type, under the the Open Font License. The logo is under a [Creative Commons BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).


