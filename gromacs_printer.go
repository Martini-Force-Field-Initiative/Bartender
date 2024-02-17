/*
 * gromacs_printer.go, part of Bartender
 *
 *
  *
 * Copyright 2023 Raul Mera  <rmeraa{at}academicos(dot)uta(dot)cl>
 *
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *
*/

package main

import (
	"fmt"
	"math"
	"os"
)

// This seems like a good easy-to-implement compromise.
// all in kJ/mol
const bondcutoff float64 = 5000
const constraintCutoff0 float64 = 15000 //this and the following modified per Gil's request.
const constraintCutoff1 float64 = 20000
const constraintCutoff2 float64 = 50000 //this might be too high now, but can't do harmm. It could print
//unnecessary lines to the itp, but the lines will be commented out.

func PrintBonded(params map[string][]*bonded, outname string, rebcomment []bool, command ...string) {
	fout, err := os.Create(outname)
	if err != nil {
		panic("Can't create a file for the Gromacs parameters: " + err.Error())
	}

	fout.WriteString("; Topology by Bartender - www.github.com/rmera/bartender\n")
	fout.WriteString("; Please cite the Bartender reference: XXXXXXXXXXX\n;\n")
	if len(command) > 0 {
		fout.WriteString("; Bartender called as:\n")
		fout.WriteString(fmt.Sprintf("; %s\n;\n", command[0]))
	}
	//bonds
	fout.WriteString("[bonds]\n; i j  funct    length   force.c.\n")
	for _, v := range params["bonds"] {
		eq := v.params[0]
		k := v.params[1]
		rmsd := v.rmsd
		str := v.Comment()
		if k >= constraintCutoff1 {
			str += ";"
			if k >= constraintCutoff2 {
				continue
			}
		}
		//This is new. Added by Gil's request.
		if k < bondcutoff {
			k = bondcutoff
		}
		str += fmt.Sprintf("%3d %-3d 1      %5.3f     %8.2f ; rmsd: %8.2f\n", v.beads[0]+1, v.beads[1]+1, eq, k, rmsd)
		fout.WriteString(str)
	}

	//constraints
	fout.WriteString("[constraints]\n; i j  funct    length  \n")
	for _, v := range params["bonds"] {
		eq := v.params[0]
		k := v.params[1]
		rmsd := v.rmsd
		str := v.Comment()
		if k < constraintCutoff1 {
			str += ";"
			if k < constraintCutoff0 {
				continue
			}
		}
		str += fmt.Sprintf("%3d %-3d 1      %5.3f     %8.2f ; rmsd: %8.2f\n", v.beads[0]+1, v.beads[1]+1, eq, k, rmsd)
		fout.WriteString(str)
	}

	//angles
	strCos := ""
	fout.WriteString("[angles]\n; i     j       k       funct   angle   force_constant\n")
	for _, v := range params["angles"] {
		if v.functype == 1 {
			eq := v.params[0]
			k := v.params[1]
			b := v.beads
			str := fmt.Sprintf("%s%3d     %-3d     %-3d       %2d     %5.2f  %8.2f ; rmsd: %8.2f\n", v.Comment(), b[0]+1, b[1]+1, b[2]+1, v.functype, eq, k, v.rmsd)
			str = lmFixes("angles", v, str, params) //last moment fix
			fout.WriteString(str)

		}
		if v.functype == 2 {
			//Now the parameters for the Cos angle fit.
			eqcos := v.params[0]
			kcos := v.params[1]
			b := v.beads
			strCos += fmt.Sprintf("%s%3d     %-3d     %-3d       %2d     %5.2f  %8.2f ; Harmonic-Cos (Gromos96) potential rmsd: %8.2f\n", v.Comment(), b[0]+1, b[1]+1, b[2]+1, v.functype, eqcos, kcos, v.rmsd)
		}

	}
	fout.WriteString(strCos) //After processing the dihedrals, we print the alternative, Cos-Angle fit for all the angles (commented out).

	//angles
	//got to make changes here
	fout.WriteString("[angles]\n; i     j       k       funct   angle   force_constant\n")
	fout.WriteString("; ReB\n")
	for i, v := range params["reb"] {
		eq := v.params[0]
		k := v.params[1]
		b := v.beads
		if rebcomment != nil && len(rebcomment) > i && rebcomment[i] {
			v.SetComment(true)
		}
		str := fmt.Sprintf("%s%3d     %-3d     %-3d       %2d     %5.2f  %8.2f ; rmsd: %8.2f \n", v.Comment(), b[0]+1, b[1]+1, b[2]+1, v.functype, eq, k, v.rmsd)
		str = lmFixes("reb", v, str, params) //last moment fix
		fout.WriteString(str)

	}

	//dihedrals
	fout.WriteString("[dihedrals]\n; i     j       k    l       funct   phase       kd    pn\n")
	str2 := ""
	str3 := ""
	for _, v := range params["dihe"] {
		if v == nil {
			continue
		}
		if v.functype == 1 {
			eq := v.params[0]
			k := v.params[1]
			n := int(math.Round(v.params[2]))
			b := v.beads
			str := fmt.Sprintf("%s%3d     %-3d     %-3d  %-3d       %2d     %5.2f  %8.2f   %1d ; rmsd: %8.2f\n", v.Comment(), b[0]+1, b[1]+1, b[2]+1, b[3]+1, v.functype, eq, k, n, v.rmsd)
			fout.WriteString(str)
		}
		if v.functype == 3 {
			p := v.params
			b := v.beads
			str2 += fmt.Sprintf("%s%3d     %-3d     %-3d  %-3d       %2d     %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f ;; Ryckaert-Belleman's potential. rmsd: %8.2f \n", v.Comment(), b[0]+1, b[1]+1, b[2]+1, b[3]+1, v.functype, p[0], p[1], p[2], p[3], p[4], p[5], v.rmsd)
		}
		if v.functype == 11 {
			p := v.params
			b := v.beads
			str3 += fmt.Sprintf("%s%3d     %-3d     %-3d  %-3d       %2d     %5.2f %5.2f %5.2f %5.2f %5.2f  %5.2f ;; Combined bending-torsion potential. rmsd: %8.2f \n", v.Comment(), b[0]+1, b[1]+1, b[2]+1, b[3]+1, v.functype, p[0], p[1], p[2], p[3], p[4], p[5], v.rmsd)
		}

	}
	fout.WriteString(str2) //After processing the dihedrals, we print the alternative, R-B fit for all the angles
	fout.WriteString(str3) //After processing the dihedrals, we print the alternative, Combined B-T fit for all the angles (commented out).

	//improper
	fout.WriteString(";Improper \n; i     j       k    l       funct   angle       kd    \n")
	for _, v := range params["improp"] {
		eq := v.params[0]
		k := v.params[1]
		b := v.beads
		//	n := int(math.Round(param["dihe"][i][2]))
		str := fmt.Sprintf("%s%3d     %-3d     %-3d  %-3d       %2d     %5.2f  %8.2f ; rmsd: %8.2f   \n", v.Comment(), b[0]+1, b[1]+1, b[2]+1, b[3]+1, v.functype, eq, k, v.rmsd)
		fout.WriteString(str)

	}
	fout.Close()
}
