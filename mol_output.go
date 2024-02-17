/*
 * mol_output.go, part of Bartender
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
	"strings"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/traj/dcd"
	v3 "github.com/rmera/gochem/v3"
)

func MakePDB(coord *v3.Matrix, mol *chem.Molecule, beads []BeadCoord) {
	nvsites := NVsites(beads)
	binterval := 100.0 / (float64(len(beads)) - float64(nvsites))
	bfacs := make([]float64, mol.Len())
	for i, _ := range bfacs {
		bfacs[i] = -1
	}
	mol2 := chem.NewTopology(0, 1, nil)
	vsitesread := 0
	var vscoord *v3.Matrix
	if nvsites > 0 {
		vscoord = v3.Zeros(nvsites)
	}
	finalcoord := v3.Zeros(coord.NVecs() + nvsites)
	for i := 0; i < mol.Len(); i++ {
		mol2.AppendAtom(mol.Atom(i))
	}
	for i, v := range beads {
		c, in := v(coord, mol)
		if in == nil {
			vscoord.SetVecs(c, []int{vsitesread})
			at := new(chem.Atom)
			at.Symbol = "U"
			at.Name = "U"
			at.MolName = "VIR"
			at.MolID = i + 1
			mol2.AppendAtom(at)
			bfacs = append(bfacs, 0.0) //This will fail if the virtual sites are not by the end of the beads slice.
			vsitesread++
			continue
		}
		for _, w := range in {
			ifloat := float64(i)
			at := mol.Atom(w)
			at.MolID = i + 1
			if bfacs[w] == -1 {
				bfacs[w] = binterval * ifloat

			} else {
				bfacs[w] = (binterval*ifloat + bfacs[w]) / 2 //this is for atoms that are shared among beads. It will only help if the sharing beads come one after the other (say, if the atom is shared between beads 7 and 8, but not if the sharing beads are 8 and 10.
			}
		}
	}
	if vscoord != nil {
		finalcoord.Stack(coord, vscoord)
		chem.PDBFileWrite("Beads.pdb", finalcoord, mol2, bfacs)
		return
	}
	chem.PDBFileWrite("Beads.pdb", coord, mol2, bfacs)
}

// Saves the multi-xyz trajectory in the file trajname to a DCD trajectory in the file fname.
// Since this is not really  needed, it doesn't panic. Will just return an error to be printed by main.
func DCDSave(fname, trajname string) error {
	if !strings.Contains(trajname, ".trj") {
		return fmt.Errorf("Option only available for multi-xyz trajectories, such as those produced by xtb")
	}
	mol, err := chem.XYZFileRead(trajname)
	if err != nil {
		return err
	}
	dcd, err := dcd.NewWriter(fname, mol.Len())
	coord := v3.Zeros(mol.Len())
	for {
		err = mol.Next(coord)
		if err != nil {
			break
		}
		err = dcd.WNext(coord)
		if err != nil {
			return err
		}
	}
	if _, ok := err.(chem.LastFrameError); ok { //just means we read the whole thing.
		return nil
	}
	return err

}
