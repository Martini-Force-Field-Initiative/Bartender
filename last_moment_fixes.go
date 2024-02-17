package main

/*
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
 */

import (
	"fmt"
	"slices"
)

func lmFixes(typ string, b *bonded, prevout string, params map[string][]*bonded) string {

	switch typ {
	case "angles":
		return lmFixAngles(b, prevout, params)
	case "reb":
		return lmFixAngles(b, prevout, params)
	}
	return prevout

}

func lmFixAngles(bon *bonded, prevout string, params map[string][]*bonded) string {
	b := bon.beads
	for _, v := range params["dihe"] {
		if indihedral(bon, v) && bon.functype == 10 {
			eq, modif := angleboundaries(bon.params[0])
			if modif {
				k := bon.params[1] * 2
				help := "\n;The previous angle is involved in a dihedral, and its too close to planarity,so I'm commenting it out and offering a displaced and rescaled version here\n"
				return ";" + prevout + help + fmt.Sprintf("%3d     %-3d     %-3d       %2d     %5.2f  %8.2f \n", b[0]+1, b[1]+1, b[2]+1, bon.functype, eq, k)
			}
		} else if indihedral(bon, v) {
			return ";" + prevout //we comment out any non-type 10 angle that is involved in a dihedral

		}

	}
	return prevout
}

func indihedral(a *bonded, d *bonded) bool {
	b := a.beads
	c := slices.Contains[[]int, int]
	return c(d.beads, b[0]) && c(d.beads, b[1]) && c(d.beads, b[2])

}

func angleboundaries(a float64) (float64, bool) {
	if a < 40 && a >= 0 {
		return 40.0, true
	}
	if a <= 0 && a > -10 {
		return -40.0, true
	}
	if a > 140 && a <= 180 {
		return 140.0, true
	}
	if a > 320 && a <= 360 {
		return 320, true
	}
	if a < -320 && a >= -360 {
		return -320, true
	}
	return a, false
}
