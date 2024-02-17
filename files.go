/*
 * files.go, part of Bartender
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
	"bufio"
	"fmt"
	"log"
	"os"
	"regexp"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
)

// This should, at some point, replace the [][]int in the map[string][][]int currently
// returned by ParseInputGeo.
type Wanted struct {
	d       []int
	comment bool
}

func NewWanted(data []int, comment bool) *Wanted {
	w := new(Wanted)
	if len(data) > 0 {
		w.d = data
	} else {
		w.d = make([]int, 0)
	}
	w.comment = comment
	return w
}

// ParseInputGeo reads a Bartender input file and returns a map of parameters (bonds, angles, reb, dihe, and improp)
// and a 2D slice of marked atoms. The input file is expected to have a specific format. If the file cannot be opened,
// the function panics with an error message.
func ParseInputGeo(inpname string) (map[string][][]int, [][]int) {

	param := map[string][][]int{
		"bonds":  make([][]int, 0),
		"angles": make([][]int, 0),
		"reb":    make([][]int, 0),
		"dihe":   make([][]int, 0),
		"improp": nil,
	}
	var marked [][]int
	reading := ""
	finp, err := os.Open(inpname)
	if err != nil {
		panic("Failed to open Bartender input: " + err.Error())
	}
	inp := bufio.NewReader(finp)
	linenu := 0
	for {
		linenu++
		line, err := inp.ReadString('\n')
		if err != nil { //inefficient, (errs[1] can be checked once before), but clearer.
			if strings.Contains(err.Error(), "EOF") {
				break
			} else {
				p := fmt.Sprintf("Failed to read line %d in the Bartender input file: %s", linenu, err.Error())
				panic(p)
			}
		}
		if strings.HasPrefix(line, "#") {
			continue //comment line
		}
		if strings.HasPrefix(line, "BEADS") {
			reading = ""
			continue
		}
		if strings.HasPrefix(line, "VSITES") {
			reading = ""
			continue
		}

		if strings.HasPrefix(line, "BONDS") {
			reading = "bonds"
			continue
		}
		if strings.HasPrefix(line, "REB") {
			reading = "reb"
			continue
		}

		if strings.HasPrefix(line, "ANGLES") {
			reading = "angles"
			continue
		}
		if strings.HasPrefix(line, "DIHEDRALS") {
			reading = "dihe"
			continue
		}
		if strings.HasPrefix(line, "IMPROPERS") {
			reading = "improp"
			param["improp"] = make([][]int, 0, 0) //there may not always be impropers, so this one is only created here. Maybe I should do this for angles and dihedrals too.
			continue
		}
		if reading == "" {
			continue
		}
		pf := strings.ReplaceAll(line, " ", "")
		pf = strings.ReplaceAll(pf, "\n", "")
		if pf == "" {
			continue //shouldn't happen, but you know how users are.
		}
		star := strings.HasSuffix(pf, "*")
		pf = strings.TrimRight(pf, "*")
		fields := strings.Split(pf, ",")
		nums := make([]int, len(fields))
		for i, v := range fields {
			nums[i], err = strconv.Atoi(v)
			if err != nil {
				p := fmt.Sprintf("Failed to parse the %d field in the %d line of the Bartender input file: %s", i, linenu, err.Error())
				panic(p)
			}
			nums[i]-- //convert from 1-based indexes to 0-based

		}
		param[reading] = append(param[reading], nums)
		if star {
			marked = append(marked, nums)
		}

	}
	finp.Close()
	return param, marked
}

// Parses the input file, returns an slice of functions, where the nth function will yield the coordinates of the nth bead given the current state of the atomistic molecule
// and a slice of slices of float64, where the nth slice contains, for each atom, the fraction of the
// that atom that belonging to the nth bead. Thus, a bead can contain half an atom, for instance.
func ParseInputBead(inpname string) ([]BeadCoord, [][]float64) {
	beadslice := make([]BeadCoord, 0, 0)
	wslice := make([][]float64, 0, 0)
	finp, err := os.Open(inpname)
	if err != nil {
		panic("Failed to open the Bartender input file: " + err.Error())
	}
	defer finp.Close()
	inp := bufio.NewReader(finp)
	reading := false
	linenu := 0
	for {
		linenu++
		line, err := inp.ReadString('\n')
		if line == "\n" {
			continue
		}
		if err != nil { //inefficient, (errs[1] can be checked once before), but clearer.
			if strings.Contains(err.Error(), "EOF") {
				break
			} else {
				p := fmt.Sprintf("Failed to read line %d in Bartender input: %s", linenu, err.Error())
				panic(p)
			}
		}
		if strings.HasPrefix(line, "#") {
			continue //comment line
		}
		if strings.HasPrefix(line, "BEADS") {
			reading = true
			continue
		}
		if strings.HasPrefix(line, "BONDS") { //the "BEADS" part is over
			break
		}
		if strings.HasPrefix(line, "VSITES") { //the "BEADS" part is over
			break
		}
		if !reading {
			continue
		}
		//now the actual reading!
		pf := strings.ReplaceAll(strings.Fields(line)[1], " ", "")
		fields := strings.Split(pf, ",")
		indexes := make([]int, len(fields))
		wslice = append(wslice, make([]float64, len(fields)))
		m1 := len(wslice) - 1
		for i, v := range fields {
			var bead string
			var weight string = "1.0"
			wslice[m1][i] = 1.0
			bead = v
			if strings.Contains(v, "/") {
				info := strings.Split(v, "/")
				bead = info[0]
				weight = info[1]
			}
			wslice[m1][i], err = strconv.ParseFloat(weight, 64)
			if err != nil {
				p := fmt.Sprintf("Failed to parse the %d field in the %d line of the Bartender input file: %s", i, linenu, err.Error())
				panic(p)
			}
			wslice[m1][i] = 1 / wslice[m1][i]
			indexes[i], err = strconv.Atoi(bead)
			if err != nil {
				p := fmt.Sprintf("Failed to convert bead id: %s in line %d to int in the Bartender input file: %s", bead, linenu, err.Error())
				panic(p)
			}
			indexes[i]-- //to convert from 1-based indexes to 0-based indexes

		}
		beadslice = append(beadslice, Bead(indexes, wslice[m1]))

	}
	return beadslice, wslice
}

// Parses the input file, returns an slice of functions, where the nth function will yield the coordinates of the nth virtual sitethe current state of the atomistic molecule. It requires as an input the equivalent slice of functions for the regular beads.
func ParseInputVsites(inpname string, beads []BeadCoord) []BeadCoord {
	vsiteslice := make([]BeadCoord, 0, 0)
	finp, err := os.Open(inpname)
	if err != nil {
		panic("Failed to open the Bartender input file: " + err.Error())
	}
	inp := bufio.NewReader(finp)
	reading := false
	linenu := 0
	for {
		linenu++ //just to return better errors
		line, err := inp.ReadString('\n')
		if err != nil { //inefficient, (errs[1] can be checked once before), but clearer.
			if strings.Contains(err.Error(), "EOF") {
				break
			} else {
				p := fmt.Sprintf("Failed to read line %d in Bartender input: %s", linenu, err.Error())
				panic(p)
			}
		}
		if strings.HasPrefix(line, "#") {
			continue //comment line
		}
		if strings.HasPrefix(line, "VSITES") {
			reading = true
			continue
		}
		if strings.HasPrefix(line, "BONDS") { //the "VSITES" part is over
			break
		}
		if !reading {
			continue
		}
		//now the actual reading!
		//The line format for this section is:
		//beadnumber  beadcomponent1,...,beadcomponentn par1,par2,..parn gromacsfunctionnumber
		areas := strings.Fields(line)
		inf := strings.ReplaceAll(areas[1], " ", "")
		fields := strings.Split(inf, ",")
		indexes := make([]int, len(fields))
		pf := strings.ReplaceAll(areas[2], " ", "")
		parfields := strings.Split(pf, ",")
		parameters := make([]float64, len(parfields))
		for i, v := range fields {
			var bead string
			bead = v
			indexes[i], err = strconv.Atoi(bead)
			if err != nil {
				p := fmt.Sprintf("Failed to convert bead id: %s in line %d to int in the Bartender input file: %s", bead, linenu, err.Error())
				panic(p)
			}
			indexes[i]-- //to convert from 1-based indexes to 0-based indexes

		}

		//we parse the virtual site parameters now
		for i, v := range parfields {
			var par string
			par = v
			parameters[i], err = strconv.ParseFloat(par, 64)
			if err != nil {
				p := fmt.Sprintf("Failed to convert parameters: %s in line %d to int in the Bartender input file: %s", par, linenu, err.Error())
				panic(p)
			}

		}
		//finally, the function number

		funcn, err := strconv.Atoi(areas[3]) //if this panics, its the users fault :-)
		if err != nil {
			p := fmt.Sprintf("Failed to convert the Gromacs function number in line %d to int in the Bartender input file: %s", linenu, err.Error())
			panic(p)
		}

		var function BeadCoord
		switch len(indexes) {
		case 2:
			function = Vsite2(beads[indexes[0]], beads[indexes[1]], parameters[0], funcn)
		case 3:
			in := indexes
			tp := parameters
			if funcn == 4 {
				function = Vsite3(beads[in[0]], beads[in[1]], beads[in[2]], tp[0], tp[1], tp[2], funcn)

			} else {
				function = Vsite3(beads[in[0]], beads[in[1]], beads[in[2]], tp[0], tp[1], 0, funcn)
			}

		}
		vsiteslice = append(vsiteslice, function)
	}
	finp.Close()
	return vsiteslice
}

// This is a custom function that can detect and read
// the PDBs made with LigParGen, which are so wrong
// they don't actually classify as PDB at all.
// This is a pretty quick and dirty function, and I'd say quite brittle.
// It relies on LigParGen "PDBs" actually following a format where the fields
// are separated by spaces, which I am not sure is the case.
func PoorlyMadePDBFileRead(filename string) (*chem.Molecule, error) {
	fin, err := scu.NewMustReadFile(filename)
	if err != nil {
		return nil, err
	}
	symbolre := regexp.MustCompile("[a-zA-Z]+")
	first := fin.Next()
	if !strings.Contains(first, "REMARK LIGPARGEN GENERATED PDB") {
		//It it's not from LigParGen, we assume it is a reasonable PDB
		fin.Close()
		return chem.PDBFileRead(filename)
	}
	//top:=chem.NewTopology()
	//no luck, we have a LigParGen "pdb" file.
	ats := make([]*chem.Atom, 0, 10)
	tmpcoord := make([][3]float64, 0, 10)
	cont := 0
	for i := fin.Next(); i != "EOF"; i = fin.Next() {
		if !strings.HasPrefix(i, "ATOM") && !strings.HasPrefix(i, "HETATM") {
			continue
		}
		at := new(chem.Atom)
		chunks := strings.Fields(i)
		at.Symbol = symbolre.FindString(chunks[2])
		at.Name = at.Symbol
		at.MolName = chunks[3]
		at.ID = cont + 1
		at.MolID, err = strconv.Atoi(chunks[4])
		if err != nil {
			log.Printf("Couldn't obtain MolID for an atom in a LigParGen PDB: %s. Will set to 1", err.Error())
			at.MolID = 1
		}
		var coord [3]float64
		for i, c := range chunks[5:] {
			coord[i], err = strconv.ParseFloat(c, 64)
			if err != nil {
				return nil, err
			}
		}
		cont++
		tmpcoord = append(tmpcoord, coord)
		ats = append(ats, at)
	}
	coord := v3.Zeros(len(tmpcoord))
	for i, v := range tmpcoord {
		coord.Set(i, 0, v[0])
		coord.Set(i, 1, v[1])
		coord.Set(i, 2, v[2])
	}
	fin.Close()
	coords := make([]*v3.Matrix, 1)
	coords[0] = coord
	top := chem.NewTopology(0, 1, ats)
	return chem.NewMolecule(coords, top, nil)
}
