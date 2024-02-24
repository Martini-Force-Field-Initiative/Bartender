/*
 * fitcontrol.go, part of Bartender
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
	"math"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	"github.com/rmera/scu"
)

func f2s(f float64) string {
	return fmt.Sprintf("%5.3f", f)
}

func formatFloats(a []float64) string {
	s := make([]string, 0, len(a))
	for _, v := range a {
		s = append(s, fmt.Sprintf("%5.3f", v))
	}
	return strings.Join(s, ", ")
}
func formatInts(a []int) string {
	s := make([]string, 0, len(a))
	for _, v := range a {
		s = append(s, fmt.Sprintf("%d", v))
	}
	return strings.Join(s, ", ")
}

// fit output
type fitoutput struct {
	out    [][]string
	levels []int
}

func (f *fitoutput) addOut(level int, out ...string) {
	f.out = append(f.out, out)
	f.levels = append(f.levels, level)
}
func (f *fitoutput) addErr(level int, err error) {
	if err == nil {
		return
	}
	f.out = append(f.out, []string{err.Error()})
	f.levels = append(f.levels, level)
}
func (f *fitoutput) printOut() {
	for i, v := range f.out {
		nv := strings.Join(v, " ")
		LogV(f.levels[i], nv)
	}
}
func (f *fitoutput) String() string {
	var nv []string
	for _, v := range f.out {
		nv = append(nv, strings.Join(v, " "))
	}
	return strings.Join(nv, "\n")
}

type StatAndFitOptions struct {
	ra                bool
	increment         float64
	temperature       float64
	histcutoff        float64
	improperTolerance float64
	noplot            bool
}

type goodEnough struct {
	dihe     map[string][]float64
	bonds    map[string][]float64
	angles   map[string][]float64
	improp   map[string][]float64
	reb      map[string][]float64
	maxiters int
	iters    int
}

func (g *goodEnough) MaxIters(i int) int {
	ret := g.maxiters
	if i >= 0 {
		g.maxiters = i
	}
	return ret
}

func (g *goodEnough) Iterate() bool {
	g.iters++
	return g.iters > g.maxiters
}
func (g *goodEnough) ResetIters() int {
	ret := g.iters
	g.iters = 0
	return ret
}

func (g *goodEnough) Map(m string) map[string][]float64 {
	switch m {
	case "dihe":
		return g.dihe
	case "bonds":
		return g.bonds
	case "angles":
		return g.angles
	case "impropers":
		return g.improp
	case "improp":
		return g.improp
	case "reb":
		if len(g.reb) == 0 {
			return g.angles
		}
		return g.reb
	default:
		return nil
	}

}

func (g *goodEnough) String() string {
	return fmt.Sprintf("dihe: %v\nangles: %v\nbonds:%v\nimprop:%v\n", g.dihe, g.angles, g.bonds, g.improp)

}

func DefaultGoodEnough() *goodEnough {
	g := new(goodEnough)
	g.dihe = map[string][]float64{"rmsd": {2, 2}, "k": {50, 200}}
	g.bonds = map[string][]float64{"rmsd": {1.7}, "k": {1000, 10000}}
	g.angles = map[string][]float64{"rmsd": {1.7}, "k": {50, 200}}
	g.reb = map[string][]float64{"rmsd": {1.7}, "k": {50, 200}}
	g.improp = map[string][]float64{"rmsd": {1.7}, "k": {50, 200}}
	g.maxiters = 3
	g.iters = 0
	return g
}

func GoodEnoughFromFile(name string) *goodEnough {
	g := goodEnoughFromFileInner(name)
	if g == nil {
		LogV(1, "WARNING: couldn't open or read GoodEnough file. Will use default values")
		return DefaultGoodEnough()

	}
	return g
}

//this one just panicks and recovers if anything goes wrong.
//the wrapper above deals with the resulting fallout
//The format for the file it reads is as follows
/*
dihe k v1 v2 rmsd v3 v4
angles k v1 v2 rmsd v3
impropers k v1 v2 rmsd v3
bonds k v1 v2 rmsd v3
*/
//note that dihe has 2 rmsd values, both are the upper limit for
//simple periodic and R-B
//they all have 2 k values, a lower and upper limit (in that order)
//The order of the terms (dihe, angles, etc) is variable, and so is the k/rmsd
//order
func goodEnoughFromFileInner(name string) *goodEnough {
	g := new(goodEnough)
	f, err := scu.NewMustReadFile(name)
	if err != nil {
		panic("D:")
	}
	defer func() {
		f.Close()
		recover()
	}()
	pf := func(s string) float64 {
		f, err := strconv.ParseFloat(s, 64)
		if err != nil {
			panic(err)

		}
		return f

	}
	var i string
	for i = f.Next(); i != "EOF"; i = f.Next() {
		fi := strings.Fields(i)
		typ := fi[0]
		m := make(map[string][]float64)
		for nt := 1; ; {
			switch strings.ToLower(fi[nt]) {
			case "k":
				m["k"] = []float64{pf(fi[nt+1]), pf(fi[nt+2])}
				nt += 3
			case "rmsd":
				if typ == "dihe" {
					m["rmsd"] = []float64{pf(fi[nt+1]), pf(fi[nt+2])}
					nt += 3
				} else {
					m["rmsd"] = []float64{pf(fi[nt+1])}
					nt += 2
				}
			default:
				panic("D:") //this will be recovered
			}
			if nt > 4 {
				break //the line is over
			}
		}
		switch typ {
		case "dihe":
			g.dihe = m
		case "angles":
			g.angles = m
		case "bonds":
			g.bonds = m
		case "improp":
			g.improp = m
		case "impropers":
			g.improp = m
		case "reb":
			g.reb = m
		default:
			panic("D:") //will be recovered

		}
	}
	return g
}

func StatAndFit(F *Funcs, typ string, beads []int, n int, inp []float64, O *StatAndFitOptions, g *goodEnough) []*bonded {
	points, E := IBoltzmann(inp, O.increment, O.temperature, O.ra, O.histcutoff) //doesn't return anything for now, but prints intermediate data.
	category := CategoryName(typ)
	FO := &fitoutput{}
	params := make([]*bonded, 0, 1)
	beadstext := BeadsText(beads)
	FO.addOut(2, "Term:", category, beadstext)
	FO.addOut(2, "Points:  [", formatFloats(points), "]")
	FO.addOut(2, "Energies:  [", formatFloats(E), "]")
	switch typ {
	case "dihe":
		//some dihedrals have very shallow wells. Those are
		//better excluded. This is a very unsophisticated way
		//of doing it.
		FO.addOut(3, "Maximum Well of the dihedral potential:", f2s(MaxWellDepth(points, E)), "kJ/mol", beadstext)
		if MaxWellDepth(points, E) <= chem.RkJ*O.temperature {
			FO.addOut(2, "Dihedral %s excluded due to potential <= RT", beadstext)
			return nil
		}
		par, R2 := F.SimplePeriodicFit(points, E)
		if par[0] == 0 && par[1] == 0 && par[2] == 0 {
			FO.addOut(1, fmt.Sprintf("WARNING: Periodic fit %s failed", beadstext))
			FO.addOut(1, fmt.Sprintf("Note that the R-B fit for %s is likely still available", beadstext))
		}
		if par[0] < 0 {
			FO.addOut(2, "Dihedral corrected. Was", f2s(par[0]), "became", f2s(2*math.Pi+par[0]), "Periodicity was", f2s(par[2]), "Absolute value has been taken")
			par[0] = 2*math.Pi + par[0]
			par[2] = math.Abs(par[2])
		}

		ploterr := Plot(sperf(par), points, E, fmt.Sprintf("Simple_periodic_%s", beadstext), O.noplot)
		if verb >= 2 && ploterr != nil {
			FO.addOut(2, "Plot failed:", ploterr.Error())
		}
		par[0] = par[0] * chem.Rad2Deg
		comment := false
		if R2 > 10 {
			comment = true
		}
		b := NewBonded(n, beads, par, R2, 1, comment)
		params = append(params, b)
		FO.addOut(2, fmt.Sprintf("S. Periodic. fit for the %s between  beads %s: eq: %5.3f k: %5.3f n: %5.3f Fit RMSD: %5.3f\n", category, beadstext, par[0], par[1], par[2], R2))
		par2, R22 := F.RyckBelleFit(points, E)
		FO.addErr(3, Plot(rybef(par2), points, E, fmt.Sprintf("Ryckaert-Belleman_%s", beadstext), O.noplot))
		b = NewBonded(n, beads, par2, R22, 3, !comment) //so if simple periodic was commented, whis will not, and viceversa.
		FO.addOut(2, fmt.Sprintf("Ryckaert-Bellemans fit for the %s between  beads %s: C1: %5.3f C2: %5.3f C3: %5.3f C4 %3.5f C5 %3.5f Fit RMSD: %5.3f\n",
			category, beadstext, par2[0], par2[1], par2[2], par2[3], par2[4], R22))
		params = append(params, b)

	case "improp":
		par, R2 := F.HookeFitAv(points, E)
		FO.addErr(3, Plot(hookef(par), points, E, fmt.Sprintf("Improper_Hooke_%s", beadstext), O.noplot))

		FO.addOut(2, fmt.Sprintf("Hooke fit for the %s between  beads %s: eq: %5.3f k: %5.3f Fit RMSD: %5.3f\n", category, beadstext, par[0]*chem.Rad2Deg, par[1], R2))
		par[0] = par[0] * chem.Rad2Deg
		////we force the improper dihedrals to be 0 or 180 degrees, to avoid a discontinuity in some of the functions
		if math.Abs(par[0]-180) < O.improperTolerance {
			FO.addOut(2, "The previous equilibrium angle will be set to 180 deg! Check that it is close enough to that value\n")
			par[0] = 180
		} else if math.Abs(par[0]) < O.improperTolerance {
			FO.addOut(2, "The previous equilibrium angle will be set to 0 deg! Check that it is close enough to that value\n")
			par[0] = 0

		} else {
			FO.addOut(2, "The previous equilibrium angle is too far from 180 or 0 to set it to 180 so it will be left as-is. This could cause numerical problems in some functions. Check that it is what you want\n")
		}
		b := NewBonded(n, beads, par, R2, 2, false)

		params = append(params, b)

	case "angles":
		par, R2 := GoHookeFitAv(points, E)
		FO.addErr(3, Plot(hookef(par), points, E, fmt.Sprintf("Angle_Hooke_%s", beadstext), O.noplot))
		par = append(par, R2)

		par[0] = par[0] * chem.Rad2Deg
		b := NewBonded(n, beads, par, R2, 1, false)
		params = append(params, b)
		FO.addOut(2, fmt.Sprintf("Hooke fit for the %s between  beads %s: eq: %5.3f k: %5.3f Fit RMSD: %5.3f\n", category, beadstext, par[0], par[1], R2))

		par2, R22 := F.CosAngleFit(points, E)
		par2[0] = par2[0] * chem.Rad2Deg
		FO.addErr(3, Plot(cosanglef(par2), points, E, fmt.Sprintf("CosAngle_%s", beadstext), O.noplot))
		b = NewBonded(n, beads, par2, R22, 2, true)
		FO.addOut(3, fmt.Sprintf("Cosine Angle (Gromos96) fit for the %s between  beads %s: eq: %5.3f k: %5.3f Fit RMSD: %5.3f\n", category, beadstext, par2[0], par2[1], R22))
		params = append(params, b) //,[len(param[k])-1] = append(param[k][len(param[k])-1], par2...) //just one after the other

	case "bonds":
		par, R2 := F.HookeFitAv(points, E)
		FO.addErr(3, Plot(hookef(par), points, E, fmt.Sprintf("Bond_Hooke_%s", beadstext), O.noplot))
		b := NewBonded(n, beads, par, R2, 1, false)
		params = append(params, b)
		FO.addOut(2, fmt.Sprintf("Hooke fit for the %s between  beads %s: eq: %5.3f k: %5.3f Fit RMSD: %5.3f\n", category, beadstext, par[0], par[1], R2))
	case "reb":
		par, R2 := F.ReBFit(points, E)
		FO.addErr(3, Plot(rebf(par), points, E, fmt.Sprintf("ReB_%s", beadstext), O.noplot))
		par = append(par, R2)
		par[0] = par[0] * chem.Rad2Deg
		b := NewBonded(n, beads, par, R2, 10, false)
		params = append(params, b)
		FO.addOut(2, fmt.Sprintf("Reb fit for the %s between  beads %s: eq: %5.3f k: %5.3f Fit RMSD: %5.3f\n", category, beadstext, par[0], par[1], R2))

	}
	ok, O2 := Recomendation(typ, beadstext, params, points, E, O, g)
	if !ok {
		return StatAndFit(F, typ, beads, n, inp, O2, g)

	}
	g.ResetIters() //reset the iteration counter for the next term
	FO.printOut()  //we only print the output for the actually returning function
	return params
}

func Recomendation(typ, beads string, par []*bonded, points, E []float64, prevcall *StatAndFitOptions, g *goodEnough) (bool, *StatAndFitOptions) {
	//this is pretty much a stand-in for something more sophisticated.
	//right now it uses simple RMSD and force-constant cutoffs
	//(which are, themselves, standin for better values) but we can use
	//a classification scheme to decide when do we need to try again.
	//and also to decide what to change in the next try.
	verbo := 1
	if g.MaxIters(-1) <= 0 {
		verbo = 4
	}

	if typ == "dihe" {
		if par[0].rmsd < g.Map(typ)["rmsd"][0] || par[1].rmsd < g.Map(typ)["rmsd"][1] {
			return true, prevcall

		}
	}

	if par[0].rmsd < g.Map(typ)["rmsd"][0] {
		PrintV(verbo, "The fit for the", typ, beads, "seem OK.")
		if par[0].params[1] < g.Map(typ)["k"][0] || par[0].params[1] > g.Map(typ)["k"][1] {
			verbow := 1
			if typ == "bonds" {
				verbow = 3
			}
			PrintV(verbow, "The force constant for", typ, beads, "seems strange. You may need to scale it.")
		}
		return true, prevcall
	}
	if g.Iterate() {
		PrintV(verbo, "The max number of iterations has been reached. We'll use the previous set of parameters. ")
		PrintV(verbo, "You can re-run Bartender with the -maxiter option to keep searching for improvements.")
		return true, prevcall
	}
	PrintV(1, "The fit for the", typ, beads, "doesn't seem right")
	PrintV(1, "The current fit is", PrettyParams(typ, par))
	PrintV(2, "For points (in deg or A):", formatFloats(points))
	PrintV(2, "And energies (in kJ/mol):", formatFloats(E))
	if !prevcall.ra {
		PrintV(verbo, "I will try again with the -removeAnharmonic option")
		prevcall.ra = true
		return false, prevcall
	}
	if prevcall.histcutoff < 0.3 {
		PrintV(verbo, "I will try again with a higher histogram cutoff")
		prevcall.histcutoff += 0.05
		return false, prevcall
	}

	PrintV(verbo, "Apparently, this is the best I can do. You might want to study the plot to see if you can come with a Bartender option to fix the problems.")
	return true, prevcall

}

func PrettyParams(typ string, v []*bonded) string {
	eq := v[0].params[0]
	k := v[0].params[1]
	rmsd := v[0].rmsd
	if typ == "bonds" || typ == "angles" || typ == "improp" || typ == "reb" {
		return fmt.Sprintf("eq: %5.3f k: %5.3f RMSD: %5.3f", eq, k, rmsd)
	}

	if typ == "dihe" {
		sp := fmt.Sprintf("Simple Periodic fit\n eq: %5.3f k: %5.3f n: %5.3f RMSD: %5.3f\n", eq, k, v[0].params[2], rmsd)
		rb := fmt.Sprintf("R-B fit\n: C1: %5.3f C2: %5.3f C3: %5.3f C4 %3.5f C5 %3.5f RMSD: %5.3f\n", v[1].params[0], v[1].params[1], v[1].params[2], v[1].params[3], v[1].params[4], v[1].rmsd)
		return sp + rb
	}
	panic("Unreachable")
}

/*
func DiheHandler(F *Funcs, typ string, beadstext string, n int, points, E []float64, O *StatAndFitOptions, g *goodEnough) ([]*bonded, string) {
	FO := &fitoutput{}
	//some dihedrals have very shallow wells. Those are
	//better excluded. This is a very unsophisticated way
	//of doing it.
	FO.addOut(3, "Maximum Well of the dihedral potential:", f2s(MaxWellDepth(points, E)), "kJ/mol", beadstext)
	if MaxWellDepth(points, E) <= chem.RkJ*O.temperature {
		FO.addOut(2, "Dihedral %s excluded due to potential <= RT", beadstext)
		return nil, FO.String()
	}

	FO.addOut(2, fmt.Sprintf("S. Periodic. fit for the %s between  beads %s: eq: %5.3f k: %5.3f n: %5.3f Fit RMSD: %5.3f\n", category, beadstext, par[0], par[1], par[2], R2))
	par2, R22 := F.RyckBelleFit(points, E)
	FO.addErr(3, Plot(rybef(par2), points, E, fmt.Sprintf("Ryckaert-Belleman_%s", beadstext), O.noplot))
	b = NewBonded(n, beads, par2, R22, 3, !comment) //so if simple periodic was commented, whis will not, and viceversa.
	FO.addOut(2, fmt.Sprintf("Ryckaert-Bellemans fit for the %s between  beads %s: C1: %5.3f C2: %5.3f C3: %5.3f C4 %3.5f C5 %3.5f Fit RMSD: %5.3f\n",
		category, beadstext, par2[0], par2[1], par2[2], par2[3], par2[4], R22))
	params = append(params, b)

	par, R2 := F.SimplePeriodicFit(points, E)
	if par[0] < 0 {
		FO.addOut(2, "Dihedral corrected. Was", f2s(par[0]), "became", f2s(2*math.Pi+par[0]), "Periodicity was", f2s(par[2]), "Absolute value has been taken")
		par[0] = 2*math.Pi + par[0]
		par[2] = math.Abs(par[2])
	}

	ploterr := Plot(sperf(par), points, E, fmt.Sprintf("Simple_periodic_%s", beadstext), O.noplot)
	if verb >= 2 && ploterr != nil {
		FO.addOut(2, "Plot failed:", ploterr.Error())
	}
	par[0] = par[0] * chem.Rad2Deg
	comment := false
	if R2 > 10 {
		comment = true
	}
	b := NewBonded(n, beads, par, R2, 1, comment)
	params = append(params, b)
}
*/
