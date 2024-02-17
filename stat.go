/*
 * stat.go, part of Bartender
 *
 *
 * Copyright 2023 Raul Mera  <rmeraa{at}academicos(dot)uta(dot)cl>
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
	"sort"

	chem "github.com/rmera/gochem"
	"gonum.org/v1/gonum/stat"
)

// sidewalk walks through histo starting by the position imax, assumed
// to be the maximum**, in the direction "direction".
// After the histo starts going downhill, i.e. the values
// start decreasing, sidewalk will record  a "strike" every time
// a value is _not_ smaller than the previous (as it would always be if things were harmonic)
// after the 3rd strike (to allow for some fluctuations) it will decide that everything else
// is anharmonic, and return that point (it can actually return the first or second strike
// if they are continuous, which suggest proven that they weren't just fluctuations.
// **Be aware that sidewalk is meant to operate on _frequency_ histograms, so, the maxima there will be
// energy minima.
func sidewalk(histo []float64, imax, steps, direction int) int {
	if direction != 1 && direction != -1 {
		panic("Wrong value for the histogram climbing direction") //this is a bug, no user input should cause it
	}
	strikes := 3 //hardcoded for now, we might consider an option in the future
	strikei := make([]int, 0, 3)
	downhill := false
	downhill_pre := false
	prev := histo[imax]
	for i := 0; i <= steps; i++ {
		if i == 0 {
			continue //only to ensure that all i>0 are checked to be <=steps.
		}
		//	println("i value:", i, imin+(i*direction)) ///////////////////////////
		test := histo[imax+(i*direction)]
		if !downhill {
			if test < prev && !downhill_pre {
				downhill_pre = true
				prev = test
				continue
			}
			if test < prev {
				downhill = true
				prev = test
				continue
			}
		}
		if test >= prev {
			strikes--
			strikei = append(strikei, i)
		}
		if strikes == 0 {
			if strikei[2]-strikei[0] == 2 {
				//		println(imax+(strikei[0]*direction), direction) ////////////////////
				return imax + (strikei[0] * direction) //it means this is the actual start of the anharmonicity
			}
			if strikei[2]-strikei[1] == 1 {
				//		println(imax+(strikei[1]*direction), direction) ////////////////////
				return imax + (strikei[1] * direction) //it means this is the actual start of the anharmonicity
			}
			//	println(imax + (i * direction))
			return imax + (i * direction) //it means strikes 1 and 2 are possible outlayers, as there were
			//"non strikes" in between them. so we'll take strike 3 as the start of the anharmonic part.
		}
		prev = test
	}
	//	println("No anharmonicities! limit returned:", imax+((steps+1)*direction), direction) ////////////////////
	return imax + ((steps + 1) * direction) //this will be an index that is off-limits for the histo slice, meaning
	//no value has to be eliminated
}

// sidewalk walks through histo starting by the position imax, assumed
// to be the maximum**, in the direction "direction".
// After the histo starts going downhill, i.e. the values
// start decreasing, sidewalk will record  a "strike" every time
// a value is _not_ smaller than the previous (as it would always be if things were harmonic)
// after the 3rd strike (to allow for some fluctuations) it will decide that everything else
// is anharmonic, and return that point (it can actually return the first or second strike
// if they are continuous, which suggest proven that they weren't just fluctuations.
// **Be aware that sidewalk is meant to operate on _frequency_ histograms, so, the maxima there will be
// energy minima.
// slopetol is 1 minus the fraction of change in slope before calling it anharmonic.
// thinking of it, I might just make it a fraction of the slope, so that it is more intuitive.
func sidewalkSlope(histo []float64, imax, steps, direction int, slopetol ...float64) int {
	//The value in the next non-coment line might actually be too little tolerance,
	//so look there if you have issues. I built the thing so we can add an option for the user to set
	//the value, but the option is not currently implemented.
	slopetolerance := 0.5
	if len(slopetol) > 0 && slopetol[0] > 0 {
		slopetolerance = slopetol[0] //we allow a 30% decrease in the slope before we call it an anharmonicity
	}
	if direction != 1 && direction != -1 {
		panic("Wrong value for the histogram climbing direction") //this is a bug, no user input should cause it
	}
	strikes := 3 //hardcoded for now, we might consider an option in the future
	strikei := make([]int, 0, 3)
	downhill := false
	downhill_pre := false
	prev := histo[imax]
	prevslope := 0.0 //this is the slope of the previous point. We will just ignore the first one.
	for i := 0; i <= steps; i++ {
		if i < 2 {
			continue //we need 2 points to calculate the slope
		}
		//	println("i value:", i, imin+(i*direction)) ///////////////////////////
		test := histo[imax+(i*direction)]
		testprev := histo[imax+((i-1)*direction)]
		slope := test - testprev
		if !downhill {
			if test < prev && !downhill_pre {
				downhill_pre = true
				prev = test
				continue
			}
			if test < prev {
				downhill = true
				prev = test
				continue
			}
		}
		if test >= prev {
			strikes--
			strikei = append(strikei, i)
		} else if prevslope > 0 && slope <= slopetolerance*prevslope {
			strikes--
			strikei = append(strikei, i)
		} else {
			prevslope = slope
		}
		if strikes == 0 {
			if strikei[2]-strikei[0] == 2 {
				//		println(imax+(strikei[0]*direction), direction) ////////////////////
				return imax + (strikei[0] * direction) //it means this is the actual start of the anharmonicity
			}
			if strikei[2]-strikei[1] == 1 {
				//		println(imax+(strikei[1]*direction), direction) ////////////////////
				return imax + (strikei[1] * direction) //it means this is the actual start of the anharmonicity
			}
			//	println(imax + (i * direction))
			return imax + (i * direction) //it means strikes 1 and 2 are possible outlayers, as there were
			//"non strikes" in between them. so we'll take strike 3 as the start of the anharmonic part.
		}
		prev = test
	}
	//	println("No anharmonicities! limit returned:", imax+((steps+1)*direction), direction) ////////////////////
	return imax + ((steps + 1) * direction) //this will be an index that is off-limits for the histo slice, meaning
	//no value has to be eliminated
}

// Removes anharmonicities beyond remove_range of the data range, for each side (range, not frequency)
// by default remove_range is 0.2, so anharmonicities from 40% of the range are removed, half on each side.
func RemoveRangeAnharmonicities(histo []float64, slopetolerance float64, remove_range ...float64) []float64 {
	ran := 0.2
	if len(remove_range) > 0 {
		ran = remove_range[0]
	}
	//we first find the minimum.
	imaxleft := int(float64(len(histo)) * ran)
	imaxright := int(float64(len(histo))*1 - ran)
	//left walk
	leftlim := sidewalkSlope(histo, imaxleft, imaxleft, -1, slopetolerance)
	//right walk
	rightlim := sidewalkSlope(histo, imaxright, len(histo)-imaxright-1, 1, slopetolerance)
	//now we go ahead and set everything "off limits" to 0
	for i, _ := range histo {
		if i <= leftlim || i >= rightlim {
			histo[i] = 0
		}
	}
	return histo
}

func RemoveAnharmonicities(histo []float64) []float64 {
	//we first find the minimum.
	max := 0.0
	imax := 0
	for i, v := range histo {
		if v > max {
			max = v
			imax = i
		}
	}
	//left walk
	leftlim := sidewalk(histo, imax, imax, -1)
	//right walk
	steps_remaining := len(histo) - imax - 1
	rightlim := sidewalk(histo, imax, steps_remaining, 1)
	//now we go ahead and set everything "off limits" to 0
	for i, _ := range histo {
		if i <= leftlim || i >= rightlim {
			histo[i] = 0
		}

	}
	return histo
}

// takes a slice with values (angles, distances, dihedrals) and, from their relative abundance, obtains an energy
func IBoltzmann(inp []float64, increment, temperature float64, removeAnharmonicities bool, histcutoff float64) ([]float64, []float64) {
	sort.Float64s(inp)
	divs := make([]float64, 0, 10)
	hpoints := make([]float64, 0, 10)
	for i := inp[0] - increment; i <= inp[len(inp)-1]+increment; i += increment {
		divs = append(divs, i)
		if len(divs) >= 2 {
			hpoints = append(hpoints, (divs[len(divs)-1]+divs[len(divs)-2])/2.0)
		}
	}

	/////
	histo := make([]float64, len(divs)-1)
	histo = stat.Histogram(histo, divs, inp, nil)
	//	fmt.Println(histo) /////
	//now we invert the Maxwell-Boltzmann distribution to
	energies := make([]float64, len(histo))
	largest := 0.0
	for _, v := range histo {
		if v > largest {
			largest = v
		}
	}
	if histcutoff <= 0 {
		histcutoff = 0.0005
	}
	for i, v := range histo {
		if v <= largest*histcutoff {
			histo[i] = 0
		}
	}
	if removeAnharmonicities {
		histo = RemoveRangeAnharmonicities(histo, 0) //I changed the old behavior which ran
		//RemoveAnharmonicities here
	}
	for i, v := range histo {
		//This uses the latest goChem (devel branch), where chem.R is in kcal/molK and chem.RkJ is
		//in kJ/molK
		energies[i] = -1 * chem.RkJ * temperature * math.Log(v/largest) //math.Log is the natural log
	}
	//we remove now the points with 0 frequency
	cleanE := make([]float64, 0, len(histo))
	cleanpoints := make([]float64, 0, len(histo))
	for i, v := range histo {
		if v != 0 {
			cleanE = append(cleanE, energies[i])
			cleanpoints = append(cleanpoints, hpoints[i])
		}
	}
	LogV(2, "Points:", cleanpoints, "\nFrequencies:", histo)
	//	fmt.Println("values:", cleanpoints, "\nener:", cleanE) //////////////
	return cleanpoints, cleanE

}

type bendtor struct {
	b1  []float64
	b2  []float64
	tor []float64
}

func Newbendtor(N int) *bendtor {
	bt := new(bendtor)
	bt.tor = make([]float64, N)
	bt.b1 = make([]float64, N)
	bt.b2 = make([]float64, N)
	return bt
}

type bendtorFreq struct {
	tor []float64
	b1  []float64
	b2  []float64
	n   int
	e   float64
}

func NewbendtorFreq() *bendtorFreq {
	r := new(bendtorFreq)
	r.tor = make([]float64, 2)
	r.b1 = make([]float64, 2)
	r.b2 = make([]float64, 2)
	r.n = 0 //not reaaally needed
	return r
}

type freqs []*bendtorFreq

func (f freqs) LargestFreq() int {
	var largest int
	for _, v := range f {
		if v.n > largest {
			largest = v.n
		}
	}
	return largest
}

func (f freqs) removeZeros() freqs {
	f2 := make([]*bendtorFreq, 0, len(f))
	for _, v := range f {
		if v.n != 0 {
			f2 = append(f2, v)
		}
	}
	return freqs(f2)
}

//If any of this fails, check that I didn't copy-paste the variables wront (i.e. assigned b1 to tor, for instance).

// returns the values in f but as 3 separate, matching, slices of torsions, bend1, bend2 and energies.
func (f freqs) tbbe() ([]float64, []float64, []float64, []float64) {
	t := make([]float64, 0, len(f))
	b1 := make([]float64, 0, len(f))
	b2 := make([]float64, 0, len(f))
	e := make([]float64, 0, len(f))
	for _, v := range f {
		tor := (v.tor[0] - v.tor[1]) / 2
		bend1 := (v.b2[0] - v.b1[1]) / 2
		bend2 := (v.b2[0] - v.b2[1]) / 2
		t = append(t, tor)
		b1 = append(b1, bend1)
		b2 = append(b2, bend2)
		e = append(e, v.e)
	}
	return t, b1, b2, e

}

// CountAngles will search bt and count all the elements that have the torsion and 2 angles
// within the ranges specified in AF. The final count will be loaded to the "N" field of AF
// and the same -modified- AF will be returned.
func CountAngles(AF *bendtorFreq, bt *bendtor) *bendtorFreq {
	for i, v := range bt.tor {
		if v >= AF.tor[0] && v < AF.tor[1] {
			if bt.b1[i] >= AF.b1[0] && bt.b1[i] < AF.b1[1] {
				if bt.b2[i] >= AF.b2[0] && bt.b2[i] < AF.b2[1] {
					AF.n++
				}

			}
		}

	}
	return AF
}

//the only solution here seems to be to set up a 4D slice containing the frequencies for each combination of dihedral, angle1
//and angle2. To normalize that and then transform to a 4D slice of energies, which is then used by the fitting function.

// takes a slice with values for 2 bendings and the torsion between them. From their relative abundance, obtains an energy
// the first element in increments is the increment for the torsion, the second, for the 2 angles
func IBoltzmannBT(inpt, inpb1, inpb2, incre []float64, temperature float64) ([]float64, []float64, []float64, []float64) {
	bt := Newbendtor(len(inpb1))
	copy(bt.b1, inpb1)
	copy(bt.b2, inpb2)
	copy(bt.tor, inpt)
	sort.Float64s(inpb1)
	sort.Float64s(inpb2)
	sort.Float64s(inpt)
	allFreqs := make([]*bendtorFreq, 0, 100)
	for i := inpt[0]; i <= inpt[len(inpt)-1]+incre[0]; i += incre[0] {

		for j := inpb1[0]; j <= inpb1[len(inpb1)-1]+incre[1]; j += incre[1] {

			for k := inpb2[0]; k <= inpb2[len(inpb2)-1]+incre[2]; k += incre[2] {

				freq := NewbendtorFreq()
				freq.tor = []float64{i, i + incre[0]}
				freq.b1 = []float64{j, j + incre[1]}
				freq.b2 = []float64{k, k + incre[2]}
				freq = CountAngles(freq, bt)
				allFreqs = append(allFreqs, freq)

			}
		}
	}
	//now I need to define a type for a slice of bendtorFreq with a method that gives me the total/largest "frequencies" for the set of angles.
	F := freqs(allFreqs)
	largest := F.LargestFreq()
	F = F.removeZeros()
	for _, v := range F {
		q := float64(v.n) / float64(largest)
		v.e = -1 * chem.R * temperature * math.Log(q) //math.Log is the natural log
	}
	return F.tbbe()

}
