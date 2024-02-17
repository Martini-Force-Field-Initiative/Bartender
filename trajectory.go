/*
 * charges.go, part of Bartender
 *
 *
 * Copyright 2023 Raul Mera <rmeraa{at}academicos(dot)uta(dot)cl>
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

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/traj/dcd"
)

func OpenTraj(name string) (chem.Traj, error) {
	var traj chem.Traj
	var err error
	ext := getExtension(name)
	switch ext {
	case "pdb":
		traj, err = chem.PDBFileRead(name)
	case "xyz":
		traj, err = chem.XYZFileRead(name)
	case "dcd":
		traj, err = dcd.New(name)
	case "xtc":
		traj, err = OpenXTC(name) //now Bartender will not compile without the xdrfile libraries, which sucks.
	default:
		return nil, fmt.Errorf("Format not supported. Supported formats are multiPDB, multiXTC (x-plor/namd)-DCD and (if compiled with the xtc tag) xtc")
	}
	return traj, err
}
