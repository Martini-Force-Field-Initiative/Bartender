/*
 * geos.go, part of Bartender
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
	"strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

// Distance *in nm*!!!!!!
func distance(a, b, aux *v3.Matrix) float64 {
	aux.Sub(a, b)
	A2nm := 1 / 10.0
	return aux.Norm(2) * A2nm
}

func TrajAn(traj chem.Traj, mol chem.Atomer, beads []BeadCoord, wanted map[string][][]int, skip int) map[string][][]float64 {
	ret := map[string][][]float64{
		"bonds":  make([][]float64, 0),
		"angles": make([][]float64, 0),
		"reb":    make([][]float64, 0),
		"dihe":   make([][]float64, 0),
		"improp": nil,
	}

	//we'll write the raw data if you give a really high verbose level.
	var a *os.File
	if verb > 3 {
		var err error
		a, err = os.Create("rawdata.out")
		if err != nil {
			LogV(3, "Error creating rawdata.out file: ", err.Error(), "Will not write raw data")
		} else {
			defer a.Close()
			bigstr := []string{}
			for _, v := range []string{"bonds", "angles", "reb", "dihe", "improp"} {
				for j := range wanted[v] {
					bigstr = append(bigstr, fmt.Sprintf("%s-%s", v, strings.TrimSpace(BeadsText(wanted[v][j]))))
				}
			}
			a.WriteString(strings.Join(bigstr, " ") + "\n")
		}
	}
	//end raw-data stuff (for this function)
	normals := make([]*v3.Matrix, len(wanted["improp"]))
	for i, _ := range normals {
		normals[i] = v3.Zeros(1)
	}
	var err error
	coord := v3.Zeros(traj.Len())
	fr := 0
	for {
		fr++
		if fr%skip != 0 {
			err = traj.Next(nil)
			if err != nil {
				break
			}
			continue
		} else {
			err = traj.Next(coord)
		}
		if err != nil {
			break
		}
		tmap := FramePar(coord, mol, beads, wanted, normals, a)
		ret = UpdateMap(tmap, ret)
	}
	if _, ok := err.(chem.LastFrameError); ok { //just means we read the whole thing.
		ret = checkImpropers(ret, wanted)
		return ret
	}
	p := fmt.Sprintf("Error reading frame %d from the QM-MD trajectory: %s", fr, err.Error())
	panic(p)
	//something wrong happened, as the error is not LastFrameError

}

// This is a bit of a hack. It relies on cutoffs to decide whether to "switch" a negative
// angle to the angle+2pi equivalent. This whole thing comes from fitting a non-periodic function
// to periodic data.
func checkImpropers(ret map[string][][]float64, wanted map[string][][]int) map[string][][]float64 {
	if ret["improp"] == nil {
		return ret
	}
	anglecutoff := 0.5 //"close enough" to either pi or -pi radians. It's about 30 degrees
	for i, v := range ret["improp"] {
		c180 := make([]int, 0)  //close to 180
		cm180 := make([]int, 0) //close to -180
		for j, w := range v {
			if math.Abs(w-math.Pi) < anglecutoff {
				c180 = append(c180, j)
			}
			if math.Abs(w+math.Pi) < anglecutoff {
				cm180 = append(cm180, j)
			}
		}
		cutoff := 0.2
		if float64(len(c180)) < cutoff*float64(len(v)) && float64(len(cm180)) < cutoff*(float64(len(v))) {
			continue //we are a-ok
		}
		LogV(2, "Angle", wanted["improp"][i], "was 'recentered' aroubnd 180 degrees")
		ret["improp"][i] = switchsmallest(v, c180, cm180)

	}
	return ret
}

func switchsmallest(v []float64, inpos, inneg []int) []float64 {
	//	criterion := func(x float64) float64 {
	//		if x > 0.0 {
	//			return x - 2*math.Pi
	//		}
	//		return x
	//	}
	//	if len(inpos) > len(inneg) {
	criterion := func(x float64) float64 {
		if x < 0.0 {
			return x + 2*math.Pi
		}
		return x
	}
	//	}
	for i, w := range v {
		v[i] = criterion(w)
	}
	return v
}

func UpdateMap(temp map[string][]float64, mmap map[string][][]float64) map[string][][]float64 {
	/*
		in temp we have, for one frame, temp["angle"]{a1,a2,a3,a4}
		in mmap we have mmap["angle"]{{a11,a12,a13...}{a21,a22,a23...}}
		so I need to append, say temp["angle"][1] to mmap[angle][1]
	*/
	for k, v := range temp {
		if v == nil || len(v) == 0 {
			continue
		}
		for l, w := range v {
			if len(mmap[k]) == 0 {
				for _, _ = range v {
					mmap[k] = append(mmap[k], make([]float64, 0, 100)) //this should happen if we are on the first frame. Note that len(mmap[k]) can never be smaller than l-2
				}
			}
			mmap[k][l] = append(mmap[k][l], w)
		}
	}
	return mmap
}

// Builds a map with al the wanted data in one frame
func FramePar(coord *v3.Matrix, mol chem.Atomer, beadsf []BeadCoord, wanted map[string][][]int, normals []*v3.Matrix, a *os.File) map[string][]float64 {
	beads := make([]*v3.Matrix, len(beadsf))
	for i, v := range beadsf {
		beads[i], _ = v(coord, mol)
	}
	aux := v3.Zeros(1)
	aux2 := v3.Zeros(1)
	ret := map[string][]float64{
		"bonds":  make([]float64, 0),
		"angles": make([]float64, 0),
		"reb":    make([]float64, 0),
		"dihe":   make([]float64, 0),
		"improp": nil,
	}
	//Distances first
	for _, v := range wanted["bonds"] {
		//	LogV(true, ret["bonds"])                                                     // distance(beads[v[0]], beads[v[1]], aux), v[0], v[1])              //////////////////////////////////////////////////////////////////////////
		ret["bonds"] = append(ret["bonds"], distance(beads[v[0]], beads[v[1]], aux)) //we don't check that v has the correct lenght. You are on your own there.
	}
	for _, v := range wanted["angles"] {
		aux.Sub(beads[v[0]], beads[v[1]])
		aux2.Sub(beads[v[2]], beads[v[1]])
		ret["angles"] = append(ret["angles"], chem.Angle(aux, aux2))
	}
	for _, v := range wanted["reb"] {
		aux.Sub(beads[v[0]], beads[v[1]])
		aux2.Sub(beads[v[2]], beads[v[1]])
		ret["reb"] = append(ret["reb"], chem.Angle(aux, aux2))
	}

	for _, v := range wanted["dihe"] {
		ret["dihe"] = append(ret["dihe"], chem.DihedralAlt(beads[v[0]], beads[v[1]], beads[v[2]], beads[v[3]]))
		//	fmt.Println("values", beads[v[0]], beads[v[1]], beads[v[2]], beads[v[3]])
		//	fmt.Println("DIHEDRAL", ret["dihe"][len(ret["dihe"])-1])
	}

	ret["improp"] = make([]float64, 0)
	for j, v := range wanted["improp"] {
		//	fmt.Println(indexes[v[0]], indexes[v[1]], indexes[v[2]], indexes[v[3]], v[0], v[1], v[2], v[3]) ///////////////////////////
		//ret["improp"] = append(ret["improp"], chem.Improper(beads[v[0]], beads[v[1]], beads[v[2]], beads[v[3]]))
		ret["improp"] = append(ret["improp"], Improper(beads[v[0]], beads[v[1]], beads[v[2]], beads[v[3]], normals, j))
	}
	if a != nil {

		bigstr := []string{}
		for _, v := range []string{"bonds", "angles", "reb", "dihe", "improp"} {
			var r2d float64
			if v != "bonds" {
				r2d = chem.Rad2Deg
			} else {
				r2d = 1.0
			}
			for j := range ret[v] {
				bigstr = append(bigstr, fmt.Sprintf("%5.2f", r2d*ret[v][j]))
			}
		}

		a.WriteString(strings.Join(bigstr, " ") + "\n")
	}

	return ret
}

// Improper calculates the improper dihedral between the points a, b,c,d
// as the angle between the plane defined by a,b,c and that defined by the plane bcd
func Improper(a, b, c, d *v3.Matrix, normals []*v3.Matrix, normalindex int) float64 {
	all := []*v3.Matrix{a, b, c, d}
	for number, point := range all {
		if point == nil {
			panic(fmt.Sprintf("Improper: Vector %d is nil", number))
		}
		pr, pc := point.Dims()
		if pr != 1 || pc != 3 {
			panic(fmt.Sprintf("Improper: Vector %d has invalid shape", number))
		}
	}
	//bma=b minus a
	amb := v3.Zeros(1)
	cmb := v3.Zeros(1)
	dmc := v3.Zeros(1)
	bmc := v3.Zeros(1)
	amb.Sub(a, b) //canged from Sub(b,a)
	cmb.Sub(c, b)
	bmc.Sub(b, c)
	dmc.Sub(d, c)
	plane1 := cross(amb, cmb)
	plane2 := cross(bmc, dmc)
	angle := chem.Angle(plane1, plane2)
	newnorm := cross(plane1, plane2)
	//The first non-zero normal we find, rules
	if normals[normalindex].Norm(2) == 0 {
		normals[normalindex] = newnorm
	} else {
		if newnorm.Dot(normals[normalindex]) < 0 {
			angle *= -1 //let's seeeeeee
		}
	}
	//	if angle < (-0.888 * math.Pi) { //around 160 degrees
	//		angle += 2 * math.Pi
	//	}
	return angle
}

// This is something that can be optimized if needed.
// cross allocates a vector on each call, it could just ask for a pre-allocated vector
func cross(a, b *v3.Matrix) *v3.Matrix {
	c := v3.Zeros(1)
	c.Cross(a, b)
	return c
}
