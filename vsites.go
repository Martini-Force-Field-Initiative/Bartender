/*
 * vsites.go, part of Bartender
 *
 *
 * Copyright 2023 Raul Mera   <rmeraa{at}academicos(dot)uta(dot)cl>
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
	"math"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	"gonum.org/v1/gonum/mat"
)

// BeadCoord is a function that, given an atomistic coordinate set
// and topology, returns the coordinate of "its" CG bead.
type BeadCoord func(*v3.Matrix, chem.Atomer) (*v3.Matrix, []int)

// NVsites eturns the number of virtual sites in a slice of BeaCoords functions
func NVsites(beads []BeadCoord) int {
	ret := 0
	for _, v := range beads {
		_, in := v(nil, nil)
		if in == nil {
			ret++
		}
	}
	return ret
}

// Bead returns a function that returns the current coordinates for the corresponding bead given the full atomistic system in
// its current state. Returns also the indexes for the atoms in that bead. If the atomistic system (both
// topology and coordinates) is nil, it returns ony the indexes.
func Bead(indexes []int, w []float64) BeadCoord {
	return func(coord *v3.Matrix, mol chem.Atomer) (*v3.Matrix, []int) {
		if coord == nil && mol == nil {
			return nil, indexes
		}
		return WCOM(coord, mol, indexes, w), indexes

	}
}

//Virtual sites.
//Procedures from
//https://manual.gromacs.org/current/reference-manual/functions/interaction-methods.html#virtual-interaction-sites

// Vsite2 Given a virtual site and functions returning coordinates for the beads that define it,
// returns a function that will take the atomistic molecule and return the updated coordinates of the virtual site.
func Vsite2(i1, i2 BeadCoord, a float64, funcn int) BeadCoord {
	return func(coord *v3.Matrix, mol chem.Atomer) (*v3.Matrix, []int) {
		if coord == nil && mol == nil {
			return nil, nil
		}
		w1, _ := i1(coord, mol)
		w2, _ := i2(coord, mol)
		ret := v3.Zeros(1)
		t := v3.Zeros(1) //these are a bit wasteful, but I wanted the function to be safe for concurrent use.
		t2 := v3.Zeros(1)
		if funcn == 1 {
			t.Scale(1-a, w1) //This might not be needed, actually. I'm should check if gonum (cont.)
			t2.Scale(a, w2)  //handles the scaling of the same receiver as in A.Scale(scalar,A) (cont.)
			ret.Add(t, t2)   //A similar thing happens in all the Vsite functions.
		} else {
			t.Sub(w2, w1)
			t2.Unit(t)
			t.Scale(a, t2)
			ret.Add(w1, t)
		}
		return ret, nil
	}
}

// Vsite3, given a virtual site and functions returning coordinates for the beads that define it,
// returns a function that will take the atomistic molecule and return the updated coordinates of the virtual site.
func Vsite3(f1, f2, f3 BeadCoord, a, b, c float64, funcn int) BeadCoord {
	return func(coord *v3.Matrix, mol chem.Atomer) (*v3.Matrix, []int) {
		if coord == nil && mol == nil {
			return nil, nil
		}
		w1, _ := f1(coord, mol)
		w2, _ := f2(coord, mol)
		w3, _ := f3(coord, mol)
		ret := v3.Zeros(1)
		if funcn == 1 { //3
			t := v3.Zeros(1)
			t2 := v3.Zeros(1)
			t3 := v3.Zeros(1)
			t.Scale(1-a-b, w1)
			t2.Scale(a, w2)
			t3.Add(t, t2)
			t2.Scale(b, w3)
			ret.Add(t3, t2)
		} else if funcn == 2 { //3fd
			t1 := v3.Zeros(1)
			t2 := v3.Zeros(1)
			ijsca := v3.Zeros(1)
			jksca := v3.Zeros(1)
			t1.Sub(w2, w1)
			t2.Sub(w3, w2)
			ijsca.Scale(1-a, t1)
			jksca.Scale(a, t2)
			t1.Add(ijsca, jksca)
			t2.Unit(t1)
			t1.Scale(b, t2)
			ret.Add(w1, t1)
		} else if funcn == 3 { //3fad. Here we don't use the "c" parameter. Just set it to 0
			d := a //just easier names
			ang := b
			ij := v3.Zeros(1)
			jk := v3.Zeros(1)
			t2 := v3.Zeros(1)
			t := v3.Zeros(1)
			ter2 := v3.Zeros(1)
			ij.Sub(w2, w1)
			jk.Sub(w3, w2)
			m := ij.Dot(jk) / ij.Dot(ij)
			t.Scale(m, ij)
			t2.Sub(jk, t)
			t.Unit(t2)
			ter2.Scale(d*math.Sin(ang), t)
			t.Unit(ij)
			t2.Scale(d*math.Cos(ang), t)
			t.Add(t2, ter2)
			ret.Add(w1, t)
		} else if funcn == 4 { //3out
			ij := v3.Zeros(1)
			ik := v3.Zeros(1)
			t2 := v3.Zeros(1)
			t := v3.Zeros(1)
			t3 := v3.Zeros(1)
			ij.Sub(w2, w1)
			ik.Sub(w3, w1)
			t.Cross(ij, ik)
			t2.Scale(c, t)
			t.Scale(b, ik)
			t3.Add(t, t2)
			t2.Scale(a, ij)
			t.Add(t2, t3)
			ret.Add(w1, t)

		}
		return ret, nil
	}

}

func Vsites22Beads(b1, b2 []int, dist float64, w1, w2 []float64) ([]int, []float64) {
	vsite := make([]int, 0, len(b1)+len(b2))
	w := make([]float64, 0, len(w1)+len(w2))
	vsite = append(vsite, b1...)
	vsite = append(vsite, b2...)
	w = append(w, w1...)
	w = append(w, w2...)

	for i, v := range w {
		if i < len(w1) {
			w[i] = v * dist
		} else {
			w[i] = v * (1 - dist)
		}
	}
	return vsite, w
}

// WCOM Obtains the geometric center for a subset of atoms from mol/coord given by indexes, where the massess are additionally
// weighted by the weights slice, if given.
// this is not a very optimized function, it allocates a lot.
func WCOM(coord *v3.Matrix, mol chem.Atomer, indexes []int, weights []float64) *v3.Matrix {
	ncoord := v3.Zeros(len(indexes))
	nmol := chem.NewTopology(0, 1)
	ncoord.SomeVecs(coord, indexes)
	nmol.SomeAtoms(mol, indexes)
	//	masses := nmol.Len()
	//	if err != nil {
	//		panic(err.Error())
	//	}
	//hopefully we always get this matrix, so we don't have to allocate a new every time.
	if weights == nil {
		weights = make([]float64, nmol.Len())
		for i, _ := range weights {
			weights[i] = 1
		}
	}
	//	massD := mat.NewDense(1, len(masses), masses)
	wD := mat.NewDense(len(weights), 1, weights)
	//	fmt.Println(len(weights), ncoord.NVecs()) ////////////////////////////////////////////////////
	//	massD.MulElem(massD, wD)
	massD := wD // Martini recommends that the centroid is used, not the COM

	gr, _ := ncoord.Dims()
	//	fmt.Println(gr) //////////////////////////////////////////////////////////////////////////////////
	tmp2 := make([]float64, gr, gr)
	for i, _ := range tmp2 {
		tmp2[i] = 1
	}
	gnOnesvector := mat.NewDense(1, gr, tmp2) //gnOnes(1, gr)
	ref := v3.Zeros(gr)
	ref.ScaleByCol(ncoord, massD)
	ref2 := v3.Zeros(1)
	ref2.Mul(gnOnesvector, ref)
	ref2.Scale(1.0/mat.Sum(massD), ref2)
	return ref2
}
