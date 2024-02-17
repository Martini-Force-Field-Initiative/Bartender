package main

/*
 * energyminima.go part of Bartender
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
 */

import (
	"math"
	"math/rand"

	chem "github.com/rmera/gochem"
	"gonum.org/v1/gonum/stat"
)

// sidewalkE walks through histo starting by the position imin, assumed
// to be a minimum, in the direction "direction".
// After the histo starts going uphill, i.e. the values
// start increasing, sidewalkE will record  a "strike" every time
// a value is _not_ larger than the previous (as it would always be if things were harmonic)
// after the 3rd strike (to allow for some fluctuations) it will decide that everything else
// is anharmonic, and return that point (it can actually return the first or second strike
// if they are continuous, which suggest proven that they weren't just fluctuations.
func sidewalkE(histo []float64, imin, steps, direction int) int {
	if direction != 1 && direction != -1 {
		panic("Wrong value for the histogram climbing direction") //this is a bug, no user input should cause it
	}
	strikes := 3 //hardcoded for now, we might consider an option in the future
	strikei := make([]int, 0, 3)
	uphill := false
	uphill_pre := false
	prev := histo[imin]
	for i := 0; i <= steps; i++ {
		if i == 0 {
			continue //only to ensure that all i>0 are checked to be <=steps.
		}
		test := histo[imin+(i*direction)]
		//the first steps could not be uphill, as there could be a "plain" near the minimum
		//we ensure that we are actually going uphill before trying to detect the end of the basin
		if !uphill {
			if test > prev && !uphill_pre {
				uphill_pre = true
				prev = test
				continue
			}
			if test > prev {
				uphill = true
				prev = test
				continue
			}
		}
		if test <= prev {
			//		println("strike", strikes) /////////
			strikes--
			strikei = append(strikei, i)
		}
		if strikes == 0 {
			if strikei[2]-strikei[0] == 2 {
				//		println(imin+(strikei[0]*direction), direction) ////////////////////
				//		println("endsidewalk 2")
				//println("r", imin+(strikei[0]*direction)) //////
				return imin + (strikei[0] * direction) //it means this is the actual start of the anharmonicity
			}
			if strikei[2]-strikei[1] == 1 {
				//		println("endsidewalk 1") //////
				//	println("r", imin+(strikei[1]*direction)) ////////////////////
				//		println(imin+(strikei[1]*direction), direction) ////////////////////
				return imin + (strikei[1] * direction) //it means this is the actual start of the anharmonicity
			}
			//	println(imin + (i * direction))
			//	println("endsidewalk 0") //////
			//	println("r", imin+i*direction) ////////////////////
			return imin + (i * direction) //it means strikes 1 and 2 are possible outlayers, as there were
			//"non strikes" in between them. so we'll take strike 3 as the start of the anharmonic part.
		}
		prev = test
	}
	//	println("No anharmonicities! limit returned:", imin+((steps+1)*direction), direction) ////////////////////
	//println("endsidewalk noharmo") //////
	if direction < 0 {
		return 0
	} else {
		return len(histo)
	}
	//return imin + ((steps + 1) * direction) //this will be an index that is off-limits for the histo slice, meaning
	//no value has to be eliminated
}

// This is a very poor implementation, but should be good enough
func steep(x, tar []float64, point int, T float64, Delta ...int) float64 {
	var left, right float64
	RT := chem.RkJ * T
	delta := 6 //no reason at all
	if len(Delta) > 0 {
		delta = Delta[0]
	}
	if point > 0 {
		d := delta
		if d > point {
			d = point
		}
		_, left = stat.LinearRegression(x[point-d:point+1], tar[point-d:point+1], nil, false)
		left *= -1
	}
	if point < len(tar)-1 {
		d := delta
		if d+point >= len(tar) {
			d = len(tar) - point
		}
		_, right = stat.LinearRegression(x[point:point+d], tar[point:point+d], nil, false)
	}
	if left >= RT && right >= RT {
		return 0 //minimum found
	}
	//we see which one is more negative and go in that direction
	if math.Abs(left-right) <= RT {
		d := rand.Intn(1)
		if d == 0 {
			return left
		}
		return right

	}
	if left < right {
		return left //minus sign means we go to the left
	}
	return -1 * right //plus sign means we go to the right
}

func timesmostvisited(visited []int) int {
	most := 0
	//index:=-1
	for _, v := range visited {
		if v > most {
			most = v
			//		index=1
		}
	}
	return most
}

// StupidDescent returns the closest minimum to the location begin
// in the slice target. The function will return -1 if it reaches an extreme (beging or end)
// of the target slice more than once, meaning that no minima has been found.
func StupidDescent(x, target []float64, begin int, T float64) int {
	visited := make([]int, len(target))
	point := begin
	visited[point]++
	for dir := steep(x, target, point, T); dir != 0; dir = steep(x, target, point, T) {
		//println("search iteration", point, dir) //////////////////
		if dir == 0 { //redundant
			break
		}
		if dir > 0 {
			if point == len(target)-1 {
				return -1
			}
			point++
		} else {
			if point == 0 {
				return -1
			}
			point--
		}
		visited[point]++
		if timesmostvisited(visited) > 2 {
			break
		}
	}
	return point
}

func isInInt(s int, col []int) bool {
	for _, v := range col {
		if s == v {
			return true
		}
	}
	return false
}

func isCloseInInt(s int, col []int, total int) int {
	tolerance := 10
	for i, v := range col {
		if int(math.Abs(float64(v)-float64(s))) < total/tolerance {
			return i
		}
	}
	return -1
}

func RemoveRepeatedMinima(minima []int, y []float64, tol ...int) []int {
	tolerance := 10
	if len(tol) > 0 {
		tolerance = tol[0]
	}
	ret := make([]int, 0, len(minima)/2)
	exclu := make([]int, 0, len(minima)/2)

	for i, v := range minima {
		if isInInt(v, exclu) {
			continue
		}
		for _, w := range minima[i+1:] {
			//	if isInInt(w, exclu) {
			//		continue
			//	}
			if int(math.Abs(float64(v)-float64(w))) < len(y)/tolerance {
				if y[v] > y[w] {
					exclu = append(exclu, v)
				} else {
					exclu = append(exclu, w)
				}

			}

		}
	}
	for _, v := range minima {
		if !isInInt(v, exclu) {
			ret = append(ret, v)
		}
	}

	return ret
}

// returns the merged minimum and basin
func MergeMinima(min1, min2 int, basin1, basin2 [2]int, y []float64, Tolerance ...float64) (bool, int, [2]int) {
	tolerance := chem.RkJ * 300 //RT
	if len(Tolerance) > 0 {
		tolerance = Tolerance[0]
	}
	max := min1
	for i := min1; i <= min2; i++ {
		if y[i] > y[max] {
			max = i
		}
	}
	d1 := y[max] - y[min1]
	d2 := y[max] - y[min2]
	d := d1
	if d2 > d1 { //I am not sure if I should choose the smallest or the largest (now)  barrier
		d = d2
	}
	if d < tolerance {
		retmin := (min1 + min2) / 2
		retbasin := [2]int{basin1[0], basin2[1]}
		return true, retmin, retbasin
	}
	return false, 0, [2]int{-1, -1}
}

func SmallestOrLargest(E []float64, smallest bool) int {
	if E == nil {
		return -1
	}
	f := func(a, b float64) bool { return a > b }
	if smallest {
		f = func(a, b float64) bool { return a < b }
	}
	var retval float64 = E[0]
	var retindex int
	for i, v := range E {
		if f(v, retval) {
			retval = v
			retindex = i
		}
	}
	return retindex
}

func MaxWellDepth(points, E []float64) float64 {
	smallest := SmallestOrLargest(E, true)
	largest := SmallestOrLargest(E, false)
	return E[largest] - E[smallest]

}
