/*
 * stat.go, part of Bartender
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

/*To the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche*/

package main

import (
	"math"
	"math/rand"
	"sort"

	"github.com/gonum/floats"
	chem "github.com/rmera/gochem"
	"gonum.org/v1/gonum/diff/fd"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/optimize"
	"gonum.org/v1/gonum/stat"
)

// The cosine-based function for angles
func GoCosAngleFit(x, y []float64) ([]float64, float64) {
	score := func(par []float64) float64 {
		eq := par[0]
		k := par[1]
		//		b := par[2]
		//		c := par[3]
		var r2 float64 = 0.0
		for i, v := range y {
			p := 0.5 * k * math.Pow((math.Cos(x[i])-math.Cos(eq)), 2.0)
			r2 += math.Pow((v - p), 2.0)
		}
		return r2 / (2 * float64(len(x)))
	}
	ngrad := func(g, par []float64) {
		g = fd.Gradient(g, score, par, &fd.Settings{Formula: fd.Central})
	}
	nhess := func(hess *mat.SymDense, par []float64) {
		fd.Hessian(hess, score, par, nil) //&fd.Settings{Formula: fd.Central2nd})

	}
	guess := cosangleGuess(x, y)
	iterations := -1 //tells Fit to use its default
	ret, res := Fit(score, ngrad, nhess, guess, iterations)
	return ret, math.Sqrt(res * 2)

}

func cosangleGuess(x, y []float64) []*float64 {
	par := hookeGuess(x, y)
	K := (*par[1] / math.Pow(math.Sin(*par[0]), 2))
	eq := par[0]
	//	par[1] = &K
	return []*float64{eq, &K}

}

// The Harmonic function for bonds and angles
func GoHookeFit(x, y []float64) ([]float64, float64) {
	//produces a function that will return the sum of the squared residues for a y = 1/2*k*(x-eq)^2
	score := func(par []float64) float64 {
		eq := par[0]
		k := par[1]
		//		b := par[2]
		//		c := par[3]
		var r2 float64 = 0.0
		for i, v := range y {
			p := 0.5 * k * math.Pow((x[i]-eq), 2.0)
			r2 += math.Pow((v - p), 2.0)
		}
		return r2 / (2 * float64(len(x)))
	}
	ngrad := func(g, par []float64) {
		g = fd.Gradient(g, score, par, &fd.Settings{Formula: fd.Central})
	}
	nhess := func(hess *mat.SymDense, par []float64) {
		fd.Hessian(hess, score, par, nil) //&fd.Settings{Formula: fd.Central2nd})

	}
	//analytic gradient, but I am  likely to have made a mistake or 2.
	agrad := func(g, par []float64) {
		eq := par[0]
		k := par[1]
		var geq, gk float64
		for i, v := range y {
			p := 0.5 * k * math.Pow((x[i]-eq), 2.0)
			ppeq := 0.5 * k * 2 * (x[i] - eq) * x[i]
			ppk := 0.5 * math.Pow((x[i]-eq), 2.0)
			geq += 2 * (v - p) * v * ppeq
			gk += 2 * (v - p) * v * ppk

		}
		g[0] = geq / 2 * float64(len(x))
		g[1] = gk / 2 * float64(len(x))

	}
	_ = agrad ////////one can choose

	//	g3:=1
	//	g4 := 0.0
	guess := hookeGuess(x, y)
	iterations := -1 //tells fit to use its default
	ret, res := Fit(score, ngrad, nhess, guess, iterations)
	//	rmsd := math.Sqrt(res * 2)
	//	rmsd_tol := 50.0
	//if rmsd > rmsd_tol {
	//	ret, res = GoHookeFitAv(x, y)
	//	}
	return ret, math.Sqrt(res * 2)
}

func hookeGuess(x, y []float64) []*float64 {
	ret := make([]*float64, 2)
	geq := x[len(x)/2] //first approach
	gk := 1.0
	mindex := 0
	energy := y[0] + 10.0 // I only need it to be higher thanthe first element
	for i, v := range y {
		if v < energy {
			energy = v
			geq = x[i]
			mindex = i
		}
	}
	gkacc := 0.0
	for i, vy := range y {
		if i == 0 || i == mindex || i == len(y)-1 {
			continue
		}
		fprime := ((y[i-1]-vy)/(x[i-1]-x[i]) + (vy-y[i+1])/(x[i]-x[i+1])) / 2.0
		gkacc += fprime / (x[i] - geq)
	}
	gk = gkacc / float64(len(y))

	//	fmt.Println("Guess:", geq, gk) /////////////////////////////////////////////////////////
	ret[0] = &geq
	ret[1] = &gk
	LogV(2, "GUESS:", "eq", geq, "k", gk) ////////////////////////////////////////
	return ret

}

// Takes a slice of distance/angles (x) and one of energies (y). It finds the minima in y and divides
// both slices in "basins" around the minima (which can overlap). It then fits a Hooke function to each
// basin, and reports the Boltzmann-averaged equilibrium distance and force constant fo the basins found.
func GoHookeFitAv(x, y []float64, Temperature ...float64) ([]float64, float64) {
	//The following might be angle-specific and  have to be removed if we use the function for bonds
	if len(y) < 60 { //assuming steps of 1 deg. But I think it's reasonable in general.
		return GoHookeFit(x, y)
	}
	//end angle-specific thing
	T := 310.0
	if len(Temperature) > 0 {
		T = Temperature[0]
	}
	//	println("look for minima") /////////////
	//we start by finding all minima in the y slice
	pminima := make([]int, 0, 2)
	for p, _ := range y {
		m := StupidDescent(x, y, p, T)
		if m >= 0 && !isInInt(m, pminima) {
			pminima = append(pminima, m)
		}
	}
	minima := RemoveRepeatedMinima(pminima, y)

	LogV(3, "minima found", len(minima), minima) /////////////
	//if there is only one minima, better just use
	//the regular fitting function
	if len(minima) <= 1 {
		return positivize(GoHookeFit(x, y))
	}

	nminima := make([]int, 0, len(minima))
	//we now need to colect basins around each minimum
	basins := make([][2]int, 0, len(minima))

	//	println("look for basins") /////////////
	for _, v := range minima {
		f := sidewalkE(y, v, v, -1)
		l := sidewalkE(y, v, len(y)-v-1, 1)
		if f < 0 {
			f = 0
		}
		if l > len(y)-1 {
			l = len(y) - 1
		}
		minbasin := 10 //a small number, just to remove fluctuation in the anharmonic zone. // len(y) / 3
		//restricts the minimum width of the basin
		if l-f < minbasin {
			continue
		}
		ytest := y[f : l+1]

		//restricts the minimum depth of the basin
		if floats.Max(ytest)-floats.Min(ytest) <= chem.RkJ*T {
			continue
		}
		nminima = append(nminima, v)
		LogV(3, "basin depth", floats.Max(ytest)-floats.Min(ytest), len(basins)+1)
		basins = append(basins, [2]int{0, 0})
		in := len(basins) - 1
		basins[in][0] = f
		basins[in][1] = l
	}
	LogV(3, len(basins), "basins found")
	//We merge basins separated by small barriers (~RT)
	basins2 := make([][2]int, 0, len(basins))
	mins2 := make([]int, 0, len(nminima))
	for i, v := range basins {
		if i < 1 {
			basins2 = append(basins2, v)
			mins2 = append(mins2, nminima[i])
			continue
		}
		l := len(mins2) - 1 //also the last element of basins2
		merged, mergedmin, mergedbas := MergeMinima(mins2[l], nminima[i], basins2[l], v, y)
		if merged {
			mins2[l] = mergedmin
			basins2[l] = mergedbas
		} else {
			mins2 = append(mins2, nminima[i])
			basins2 = append(basins2, v)
		}

	}
	basins = basins2
	nminima = mins2
	LogV(3, len(basins), "Adjusted basins found and", len(nminima), "minima", nminima) /////////////
	for i, v := range nminima {                                                        ////
		f := basins[i][0]                ////////
		l := basins[i][1]                ////////
		LogV(4, "min", x[v], x[f], x[l]) ///////////
	} ///////////
	//if there is only one basin, then the regular GoHookefit gets to deal with the whole data set.
	if len(basins) <= 1 {
		return positivize(GoHookeFit(x, y))
	}
	//Now we fit each basin
	rets := make([][]float64, len(basins))
	res := make([]float64, len(basins))
	for j, w := range basins {
		yc := y[w[0] : w[1]+1]
		xc := x[w[0] : w[1]+1]
		LogV(3, "basin", j, xc, yc)
		//if the minimum is too high, we just move the whole basin to zero, so that the fit is not affected by it.
		if y[nminima[j]] > T*chem.RkJ {
			yc = make([]float64, len(yc))
			for i, v := range y[w[0] : w[1]+1] {
				yc[i] = v - y[nminima[j]]
			}

		}
		rets[j], res[j] = positivize(GoHookeFit(xc, yc)) // Fit(score, ngrad, nhess, guess, iterations)

		res[j] = math.Sqrt(res[j] * 2)

		LogV(3, "Fit function completed:", rets[j][0], rets[j][1], res[j]) /////////////
	}
	//Now we just Boltzmann-average everything

	weights := boltzmannw(y, nminima, T)
	eqs := make([]float64, len(basins))
	ks := make([]float64, len(basins))
	for i, v := range rets {
		eqs[i] = v[0]
		ks[i] = v[1]
	}
	LogV(3, "eqs,ks,weights:", eqs, ks, weights) /////////////
	//fmt.Println("Will return the weighted means", eqs, ks, weights) /////////////
	return []float64{stat.Mean(eqs, weights), meanKql(ks, weights) /*stat.Mean(ks, weights)*/}, stat.Mean(res, weights)

}

func meanKql(ks, ws []float64) float64 {
	n := float64(len(ks))
	if n == 1 {
		return ks[0]
	}
	//	fmt.Println("n =", n, ks, ws) ///////////////
	sum := 0.0
	for i, v := range ks {
		sum += (ws[i] / v)
	}
	x := stat.StdDev(ws, nil)
	den := n - 2*x*(n-1) //I'm trying to "divide" by a greater number the more spread  the weights are.
	//if one weight is much larger than the other one, then I should just treat it as it was close to just one
	//minima (den~1) if they are evenly spread, then I should divide by the number of minima (den~n). NOTE: I'm
	//not sure of any of this, so this is a place to look if a bug appears.

	r := math.Pow(sum, -1.0) / den
	//	println("final k", r)
	return r
}

func positivize(ret []float64, res float64) ([]float64, float64) {
	if ret[0] < 0 {
		ret[0] = math.Abs(ret[0])
	}
	return ret, res

}

// data is a NxM matrix where the _columns_ represent data sets, each with row elements
// the weights must have the same number of elements as the rows in the data matrix (N)
// it will return a slice with the weighted average for each data set
func averages(weights []float64, data ...[]float64) []float64 {
	avs := make([]float64, len(data))
	for i, v := range data {
		for j, _ := range v {
			avs[i] += data[i][j] * weights[j]
		}

	}
	return avs
}

func boltzmannw(energies []float64, minima []int, T float64) []float64 {
	Eb := make([]float64, len(minima))
	sum := 0.0
	for i, v := range minima {
		Eb[i] = math.Exp(-1 * energies[v] / T * chem.RkJ)
		sum += Eb[i]

	}
	for i, v := range Eb {
		Eb[i] = v / sum
	}
	return Eb

}

// returns the index of set that contains the angle "angle", whether it is defined in the same order (abc) or reverse (cba)
// it returns -1 if the angle is not found.
func angleSearch(angle []int, set [][]int) int {
	revangle := []int{angle[2], angle[1], angle[0]} //note that the central element remains the same as in "angle"
	//	fmt.Println(angle, revangle)                    /////////////////////
	for _, v := range [][]int{angle, revangle} {
		for j, w := range set {
			//	fmt.Println(w) ////////////
			if v[0] == w[0] && v[1] == w[1] && v[2] == w[2] {
				return j
			}
		}
	}
	return -1
}

func ManageBendingTorsion(datamap map[string][][]float64, wanted map[string][][]int, dihekey int, temperature float64, increments []float64) ([]float64, float64) {
	dbeads := wanted["dihe"][dihekey]
	angle1 := []int{dbeads[0], dbeads[1], dbeads[2]}
	angle2 := []int{dbeads[1], dbeads[2], dbeads[3]}
	akey1 := angleSearch(angle1, wanted["angles"])
	akey2 := angleSearch(angle2, wanted["angles"])
	if akey1 == -1 || akey2 == -1 {
		return nil, -100 //we'll use this negative value to signal issues in this function. In this case, you need to put the corresponding bending angles
		//if you want a torsion to be given a combined BT potential, if not, it will simply not be calculated, and it will not be considered
		//an error.
	}
	var tor, b1, b2 []float64
	tor = datamap["dihe"][dihekey]
	b1 = datamap["angles"][akey1]
	b2 = datamap["angles"][akey2]
	x1, x2, x3, y := IBoltzmannBT(tor, b1, b2, increments, temperature)
	//	fmt.Println(len(x1), len(x2), len(x3), len(y)) /////////////////
	//	for i, v := range y {                          ///////////
	//		fmt.Println(x1[i], x2[i], x3[i], v) ////////////
	//	} ////////////////////////////////////////////////////////////////////
	sin := math.Sin
	cos := math.Cos
	pow := math.Pow
	score := func(par []float64) float64 {
		k := par[0]
		//		b := par[2]
		//		c := par[3]
		var r2 float64 = 0.0
		for i, v := range y {
			acc := 0.0
			for j := 0; j <= 4; j++ {
				acc += par[j+1] * pow(cos(x1[i]), float64(j))
			}
			p := k * pow(sin(x2[i]), 3) * pow(sin(x3[i]), 3) * acc
			r2 = r2 + math.Pow((v-p), 2.0)
		}
		return r2 / (2 * float64(len(x1)))
	}
	ngrad := func(g, par []float64) {
		g = fd.Gradient(g, score, par, &fd.Settings{Formula: fd.Central})
	}
	nhess := func(hess *mat.SymDense, par []float64) {
		fd.Hessian(hess, score, par, nil)
	}
	var p1, p2, p3, p4, p5, p6 float64 = 1, 1, 1, 1, 1, 1 //yeah, not getting cute here.
	guess := []*float64{&p1, &p2, &p3, &p4, &p5, &p6}
	iterations := -1 //Fit will use its default
	ret, res := Fit(score, ngrad, nhess, guess, iterations)
	return ret, math.Sqrt(res * 2)

}

// The fit for the Restricted Bending potential (ReB).
// See: https://pubs.acs.org/doi/abs/10.1021/ct400219n
func GoReBFit(x, y []float64) ([]float64, float64) {
	//produces a function that will return the sum of the squared residues for a y = 1/2*k*(x-eq)^2
	score := func(par []float64) float64 {
		eq := par[0]
		k := par[1]
		//		b := par[2]
		//		c := par[3]
		var r2 float64 = 0.0
		for i, v := range y {
			//In the original file, in the next line we had " math.Pow((math.Cos(x[i])-eq), 2.0)"
			//I changed it, because I think that was a bug. From what I have seen, nowhere else was the angle transformed to cosine before.
			p := 0.5 * k * math.Pow((math.Cos(x[i])-math.Cos(eq)), 2.0) * (1 / math.Pow(math.Sin(x[i]), 2)) //I think "eq" should be "cos(eq)" !!!
			r2 += math.Pow((v - p), 2.0)
		}
		return r2 / (2 * float64(len(x)))
	}
	ngrad := func(g, par []float64) {
		g = fd.Gradient(g, score, par, &fd.Settings{Formula: fd.Central})
	}
	nhess := func(hess *mat.SymDense, par []float64) {
		fd.Hessian(hess, score, par, nil) //&fd.Settings{Formula: fd.Central2nd})

	}

	guess := reBGuess(x, y)
	iterations := -1 //tells Fit to use its default
	ret, res := Fit(score, ngrad, nhess, guess, iterations)
	return ret, math.Sqrt(res * 2)
}

func reBGuess(x, y []float64) []*float64 {
	cos := math.Cos
	ret := make([]*float64, 2)
	geq := 1.0
	gk := 1.0
	mindex := 0
	for i, v := range y {
		if v == 0 {
			geq = x[i]
			mindex = i
			break
		}
	}
	gkacc := 0.0
	for i, vy := range y {
		if i == 0 || i == mindex || i == len(y)-1 {
			continue
		}
		fprime := ((y[i-1]-vy)/(cos(x[i-1])-cos(x[i])) + (vy-y[i+1])/(cos(x[i])-cos(x[i+1]))) / 2.0
		gkacc += fprime / (x[i] - geq)
	}
	gk = gkacc / float64(len(y))

	//	fmt.Println("Guess:", geq, gk) /////////////////////////////////////////////////////////
	ret[0] = &geq
	ret[1] = &gk
	return ret

}

///The fit for the simple  function, such as the one used for dihedrals U = k(1+cos(nphi - phi_eq))
//where phi is the angle,  and phi_eq is the equilibrium angle, both in radians.

func GoSimplePeriodicFit(x, y []float64) ([]float64, float64) {
	//produces a function that will return the sum of the squared residues for a y = 1/2*k*(x-eq)^2
	score := func(par []float64) float64 {
		eq := par[0]
		k := par[1]
		n := par[2]
		//		b := par[2]
		//		c := par[3]
		var r2 float64 = 0.0
		for i, v := range y {
			p := k * (1 + math.Cos(n*x[i]-eq))
			r2 = r2 + math.Pow((v-p), 2.0)
		}
		return r2 / (2 * float64(len(x)))
	}
	ngrad := func(g, par []float64) {
		g = fd.Gradient(g, score, par, &fd.Settings{Formula: fd.Central})
	}
	nhess := func(hess *mat.SymDense, par []float64) {
		fd.Hessian(hess, score, par, nil)
	}
	guess := simplePeriodicGuess(x, y)
	iterations := 10000 //this is 3 orders of magnitude less than the default
	ret, res := Fit(score, ngrad, nhess, guess, iterations)
	//We try to forbid periodicity one by discarding the value, if we get it, incrementing the guess by a random number, and fitting again.
	macroiters := 30
	cont := 0
	if ret[0] > 2*math.Pi {
		//NOTE: of course the phases could be greater than 360*2, o less than 0, so we need to put something more sophisticated here.
		ret[0] = ret[0] - 2*math.Pi //hopefully this fix the problem of phases being larger than 360.
	}
	for (ret[2] > 3 || ret[2] < 1) && cont < macroiters {
		LogV(3, "Failed fit for simple periodic! will try increasing the n guess")
		if *guess[2] < 1.0 {
			*guess[2] = 1.0
		} else {
			*guess[2] += rand.Float64() //
		}
		ret, res = Fit(score, ngrad, nhess, guess, iterations)
		cont++
	}
	if cont >= macroiters {
		//LogV(1, "Simple periodic fit failed")
		ret = make([]float64, 3) //empty
		res = 225.0              //so, 15^2 Just a random large number to make it clear that things went wrong.
	}
	return ret, math.Sqrt(2 * res)
}

func simplePeriodicGuess(x, y []float64) []*float64 {
	ret := make([]*float64, 3)
	geq := 1.0
	gk := 1.0
	gn := 1.0
	//	mindex := 0
	minima := 0
	maxima := 0
	for i, v := range y {
		if v == 0 {
			geq = x[i]
		}
		if i == 0 || i == len(y)-1 {
			continue
		}
		if v > y[i-1] && v > y[i+1] {
			maxima++
		}
		if v < y[i-1] && v < y[i+1] {
			minima++
		}
	}
	//if maxima < minima {
	//	maxima = minima
	//}
	turns := float64(maxima + minima)
	xrange := x[len(x)-2] - x[1]
	gn = 1
	if turns > 1 {
		gn = math.Pi * turns / xrange
	}
	sortedy := make([]float64, len(y))
	copy(sortedy, y)
	sort.Float64s(sortedy)
	max := sortedy[len(y)-1]
	//this will be a bad guess if only a zone around a minimum is sampled!
	gk = (2 * math.Pi) * max / (2.0 * xrange) //the IBoltzmann function always sets the minimum to 0, so max *is* twice the amplitude of the cos function.

	LogV(2, "Guess:", geq, gk, gn, "xrange", xrange, "turns", turns)
	ret[0] = &geq
	ret[1] = &gk
	ret[2] = &gn
	return ret
}

func SommelierSimplePeriodicFit(x, y []float64, phase, n, guessk float64) ([]float64, float64) {
	//produces a function that will return the sum of the squared residues for a y = 1/2*k*(x-eq)^2
	score := func(par []float64) float64 {
		eq := phase
		n := n
		k := par[0]
		//		b := par[2]
		//		c := par[3]
		var r2 float64 = 0.0
		for i, v := range y {
			p := k * (1 + math.Cos(n*x[i]-eq))
			r2 = r2 + math.Pow((v-p), 2.0)
		}
		return r2 / (2 * float64(len(x)))
	}
	ngrad := func(g, par []float64) {
		g = fd.Gradient(g, score, par, &fd.Settings{Formula: fd.Central})
	}
	nhess := func(hess *mat.SymDense, par []float64) {
		fd.Hessian(hess, score, par, nil)
	}
	guess := []*float64{&guessk}
	iterations := 10000 //this is 3 orders of magnitude less than the default
	ret, res := Fit(score, ngrad, nhess, guess, iterations)
	//We try to forbid periodicity one by discarding the value, if we get it, incrementing the guess by a random number, and fitting again.
	return ret, math.Sqrt(2 * res)
}

func GoRyckBelleFit(x, y []float64) ([]float64, float64) {
	//produces a function that will return the sum of the squared residues for a y = 1/2*k*(x-eq)^2
	score := func(p []float64) float64 {
		//		b := par[2]
		//		c := par[3]

		var r2 float64 = 0.0
		cos := math.Cos
		pow := math.Pow
		for i, v := range y {
			psi := x[i] - math.Pi
			test := p[0] + p[1]*cos(psi) + p[2]*pow(cos(psi), 2) + p[3]*pow(cos(psi), 3) + p[4]*pow(cos(psi), 4) + p[5]*pow(cos(psi), 5)
			r2 = r2 + math.Pow((v-test), 2.0)
		}
		return r2 / (2 * float64(len(x)))
	}
	ngrad := func(g, par []float64) {
		g = fd.Gradient(g, score, par, &fd.Settings{Formula: fd.Central})
	}
	nhess := func(hess *mat.SymDense, par []float64) {
		fd.Hessian(hess, score, par, nil)
	}
	guess := ryckBelleGuess(x, y)
	iterations := -1 //tell the Fit function to use it's default
	ret, res := Fit(score, ngrad, nhess, guess, iterations)

	return ret, math.Sqrt(2 * res)
}

// yeah, only for consistency. I won't even try to get cute here.
func ryckBelleGuess(x, y []float64) []*float64 {
	guess := []*float64{new(float64), new(float64), new(float64), new(float64), new(float64), new(float64)}
	for i, _ := range guess {
		*guess[i] = 1.0
	}
	return guess
}

/************************
*
The  master fit function
*
*************************/

func Fit(score func([]float64) float64, grad func([]float64, []float64), hess func(*mat.SymDense, []float64), guessp []*float64, iter int) ([]float64, float64) {
	prob := optimize.Problem{Func: score, Grad: grad, Hess: hess}
	guess := make([]float64, len(guessp))
	for i, v := range guessp {
		if v == nil {
			guess[i] = 1 //as good a default guess as any, I suppose.
		} else {
			guess[i] = *v
		}
	}
	if iter <= 0 {
		iter = 10000000
	}
	//I could add flags to control these parameters also.
	conv := optimize.FunctionConverge{Iterations: iter}
	//	gradient := 0.000001
	ret, err := optimize.Minimize(prob, guess, &optimize.Settings{Converger: &conv, Concurrent: 4}, &optimize.BFGS{}) //&optimize.Newton{})

	//ret, err := optimize.Minimize(prob, guess, &optimize.Settings{Converger: &conv, Concurrent: 4}, &optimize.Newton{})
	if err == nil {
		LogV(2, ret.X, ret.F, "Iterations", ret.MajorIterations)

		return ret.X, ret.F

	}
	LogV(3, "couldn't fit, will try relaxed gradients") ///////////
	conv = optimize.FunctionConverge{Absolute: 1e-5, Iterations: iter}
	for gradient := 0.000005; gradient <= 0.01; gradient = gradient * 5.0 {

		//	LogV(4, "failed!", ret.X, ret.F, ret.MajorIterations, "will try now gradient:", gradient)

		LogV(3, "failed!", ret.X, ret.F, ret.MajorIterations, "will try now gradient:", gradient)
		//		ret, err := optimize.Minimize(prob, ret.X, &optimize.Settings{GradientThreshold: gradient, Converger: &conv, Concurrent: 4}, &optimize.Newton{})
		conv.Absolute *= 0.5 ///
		ret, err := optimize.Minimize(prob, ret.X, &optimize.Settings{GradientThreshold: gradient, Converger: &conv, Concurrent: 4}, &optimize.BFGS{})
		if err == nil {
			LogV(4, ret.X, ret.F, "Iterations", ret.MajorIterations)
			return ret.X, ret.F

		}
	}
	return guess, 999999.9 //a very large number to signal that things didn't work
	//panic(err.Error()) //we tried, but couldn't
}

//produces a function that will return the sum of the squared residues for a y = 1/2*k*(x-eq)^2
/**
		score := func(par []float64) float64 {
			eq := par[0]
			k := par[1]
			//		b := par[2]
			//		c := par[3]
			var r2 float64 = 0.0
			for i, v := range yc {
				p := 0.5 * k * math.Pow((xc[i]-eq), 2.0)
				r2 += math.Pow((v - p), 2.0)
			}
			return r2 / (2 * float64(len(xc)))
		}
		ngrad := func(g, par []float64) {
			g = fd.Gradient(g, score, par, &fd.Settings{Formula: fd.Central})
		}
		nhess := func(hess *mat.SymDense, par []float64) {
			fd.Hessian(hess, score, par, nil) //&fd.Settings{Formula: fd.Central2nd})

		}
		//analytic gradient, but I am  likely to have made a mistake or 2.
		agrad := func(g, par []float64) {
			eq := par[0]
			k := par[1]
			var geq, gk float64
			for i, v := range yc {
				p := 0.5 * k * math.Pow((xc[i]-eq), 2.0)
				ppeq := 0.5 * k * 2 * (xc[i] - eq) * xc[i]
				ppk := 0.5 * math.Pow((xc[i]-eq), 2.0)
				geq += 2 * (v - p) * v * ppeq
				gk += 2 * (v - p) * v * ppk

			}
			g[0] = geq / 2 * float64(len(xc))
			g[1] = gk / 2 * float64(len(xc))

		}
		_ = agrad ////////one can choose

		//	g3:=1
		//	g4 := 0.0
		guess := hookeGuess(xc, yc)
		iterations := -1 //tells fit to use its default
*******/
