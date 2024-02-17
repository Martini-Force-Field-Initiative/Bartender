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
	"log"
)

// this is more of a procedure than a function
// it checks for some common mistakes or problems in the input file.
func initialChecks(wanted map[string][][]int) (map[string][][]int, []bool) {
	rebcomment := make([]bool, len(wanted["reb"]))
	for n, v := range wanted {
		if len(v) == 0 && n != "reb" {
			LogV(2, "WARNING: No", n, "selected. This could be a mistake.")
		}

		if n == "dihe" {
			for m, _ := range v {
				var rc []bool
				wanted, rc = checkDihe(m, wanted)
				rebcomment = append(rebcomment, rc...)
			}
		}
		if n == "angles" {
			for m, _ := range v {
				wanted = checkAngles(m, wanted)
			}
		}
		if n == "bonds" {
			for m, _ := range v {
				wanted = checkBonds(m, wanted)
			}
		}
		if n == "improp" {
			for m, _ := range v {
				wanted = checkImprop(m, wanted)
			}
		}
		if n == "reb" {
			for m, _ := range v {
				wanted = checkReB(m, wanted)
			}
		}

	}
	return wanted, rebcomment
}

//The idea is that we can add more checks/warnings to any term.

func checkReB(index int, wanted map[string][][]int) map[string][][]int {
	t := wanted["reb"][index]
	checklen("reb", index, t)
	return wanted
}

func checkImprop(index int, wanted map[string][][]int) map[string][][]int {
	t := wanted["improp"][index]
	checklen("improp", index, t)
	return wanted
}

func checkBonds(index int, wanted map[string][][]int) map[string][][]int {
	t := wanted["bonds"][index]
	checklen("bonds", index, t)
	return wanted
}

func checkAngles(index int, wanted map[string][][]int) map[string][][]int {
	t := wanted["angles"][index]
	checklen("angles", index, t)
	return wanted
}

func checkDihe(index int, wanted map[string][][]int) (map[string][][]int, []bool) {
	t := wanted["dihe"][index]
	rebcomment := make([]bool, 0)
	checklen("dihe", index, t)
	needed := [][]int{{t[0], t[1], t[2]}, {t[1], t[2], t[3]}}

	for _, v := range needed {
		if angleSearch(v, wanted["reb"]) < 0 {
			if angleSearch(v, wanted["angles"]) > 0 {
				LogV(2, "Dihedral", index, "needs angle", BeadsText(v), "which is present, but as an angle. The Restricted Bending potential (ReB) could be a safer choice here")
				LogV(2, "A ReB term it will be added, commented")
				wanted["reb"] = append(wanted["reb"], v)
				rebcomment = append(rebcomment, true)
				continue
			}
			LogV(1, "WARNING: Dihedral", index, "needs angle", BeadsText(v), "but the angle is not in the input file. Will add it and continue")
			wanted["angle"] = append(wanted["angle"], v)
			wanted["reb"] = append(wanted["reb"], v)
			rebcomment = append(rebcomment, true)

		}
	}

	return wanted, rebcomment
}

// I prefer not to upgrade to generic things until
// they are in the standard library, as this kind of
// function is likely to be implemented there when that happens.
func containsIntSlice(s [][]int, e []int) bool {
	var thesame bool = true
	for _, a := range s {
		if len(a) != len(e) {
			//this should never happen
			continue
		}
		for i, v := range a {
			if v != e[i] {
				thesame = false
				break
			}
		}
		return thesame
	}
	return false
}

// just a little struct to keep the characteristics of each term
type errhelp struct {
	n   int
	nam string
}

// Checks that a given term has the correct number of beads
// if there are too many, it just ignores the last ones and prints
// a warning.
func checklen(typ string, index int, beads []int) {
	m := map[string]*errhelp{
		"angles": {3, "angle"},
		"dihe":   {4, "dihedral"},
		"reb":    {3, "reb"},
		"improp": {4, "improper"},
		"bonds":  {2, "bond"},
	}
	if len(beads) < m[typ].n {
		log.Fatalf("%s %d has less than %d beads (%d):%v.", m[typ].nam, index, m[typ].n, len(beads), BeadsText(beads))
	}
	if len(beads) > m[typ].n {
		LogV(1, fmt.Sprintf("WARNING: %s %d has more than %d beads (%d):%v. Will ignore beads after the %dth", m[typ].nam, index, m[typ].n, len(beads), BeadsText(beads), m[typ].n))
	}
}
