/*
 * main.go, part of Bartender
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
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	"gonum.org/v1/gonum/stat"
)

// Global variables... Sometimes, you gotta use'em
var verb int
var plotdir string    //there is no reason for this to be a global, but I got tired of passing it around
var warnings []string //this counts the total warnings top level warnings printed byt the program

// If level is larger or equal, prints the d arguments to stderr
// otherwise, does nothing.
func LogV(level int, d ...interface{}) {
	if level <= verb {
		fmt.Fprintln(os.Stderr, d...)
	}
	if level <= 1 {
		warnings = append(warnings, fmt.Sprintln(d...))
	}

}

// If level is larger or equal, prints the d arguments to stdout
// otherwise, does nothing. I could have just made one function
// with LogV, but I think this way it's more clear when reading the call
// and we can better separate "errors" from information. The plan is that you
// redirect stout and stderr to different files when calling Bartender.
func PrintV(level int, d ...interface{}) {
	if level <= verb {
		fmt.Println(d...)
	}

}

type Funcs struct {
	SimplePeriodicFit func([]float64, []float64) ([]float64, float64)
	RyckBelleFit      func([]float64, []float64) ([]float64, float64)
	CosAngleFit       func([]float64, []float64) ([]float64, float64)
	ReBFit            func([]float64, []float64) ([]float64, float64)
	HookeFitAv        func([]float64, []float64, ...float64) ([]float64, float64)
	HookeFit          func([]float64, []float64) ([]float64, float64)
}

// gets a file's extension, i.e. whatever its written after the last point/dot in the filename
func getExtension(name string) string {
	fs := strings.Split(name, ".")
	return strings.ToLower(fs[len(fs)-1])
}

func main() {
	//We need to deal with the flag casing.
	//it is absolutely absurd that they are case-sensitive by default.
	//
	skip := flag.Int("skip", 1, "Read only every N frames")
	//There will be _tons_ of flags, but they are meant not to be needed the 99% of the time.
	sommelier := flag.String("sommelier", "", "Invokes the sommelier functionality to test the new dihedral parameters in a short classical MD, with the given topology file. Requires Gromacs and Insane")

	cpus := flag.Int("cpus", -1, "the total CPUs used for the QM calculations. If a number <0 is given, all logical CPUs are used")
	restart := flag.Bool("restart", false, "Attempt a restart from a previous xtb MD")

	py := flag.Bool("pyfit", false, "Use Python for function fitting. The BTROOT environment variable, with the directory where bartender lives, must be set (requires Python3, Numpy and SciPy) [DEPRECATED]")
	refit := flag.Bool("refit", false, "Only do a re-fit for the bonded parameters from an existing trajectory. Equivalent to -time 1 -nobeadtype")
	noplot := flag.Bool("noplot", false, "Do not produce the plots that would normally be written for each parameter fitted")
	removeAnharmonic := flag.Bool("removeAnharmonic", false, "Remove anharmonicities from the frequencies histogram")
	plotDir := flag.String("plotdir", ".", "The plots and plot files will be written to the given directory (default: current directory)")
	owntraj := flag.String("owntraj", "", "Use the given trajectory for geometry analysis, instead of obtaining a GFN one. DCD, multi-PDB and multi-XYZ formats are allowed. XTC is allowed if the xdrfile library is installed")
	verbose := flag.Int("verbose", 1, "Print lots of additional information (mostly for debugging)")
	verbose2 := flag.Int("v", 1, "Print lots of additional information (mostly for debugging). Same as -verbose, which ever has the highest value will be used")

	mdtime := flag.Int("time", 10000, "the total simulation time, in ps. If a number <0 is given, the MD will not be performed, and a previous trajectory will be used (the program will crash if no such previous trajectory is present)")
	charge := flag.Int("charge", 0, "the total charge of the system, in a.u. Needed for partial charges calculation")
	dcdsave := flag.String("dcdSave", "", "If given, Bartender will save the xtb-calculated trajectory in DCD format with the filename given")
	method := flag.String("method", "gfn2", "The method employed in the semiempirical simulation. Valid options are gfn0, gfn1,gfn2 and gfnff")
	temperature := flag.Float64("temperature", 310, "The temperature for the MD simulation, in K")
	histcutoff := flag.Float64("histCutoff", 0.0005, "If a bin in the 'population histogram' for a term has a lowest population lower than histCutoff times the most populated bin, that bin's population is set to 0.")
	gefile := flag.String("goodenoughfile", "", "a file with info to decide when are fits good enough")
	bi := flag.Float64("bondIncrement", 0.001, "The bin width for the bond distance histograms, in nm")
	ai := flag.Float64("angleIcrcement", 1, "The bin width for the A-B-C angle histograms, in degrees")
	di := flag.Float64("dihedralIncrement", 10, "The bin width for the A-B-C-D dihedral histograms, in degrees")
	ii := flag.Float64("improperIncrement", 1, "The bin width for the A-B-C-D improper dihedral histograms, in degrees")
	improperTolerance := flag.Float64("improperTolerance", 10.0, "Improper dihedral equilibrium values that differ from 0 or from 180 degrees by less than this parameter will be rounded to 0 or 180 degrees") //I'm not too happy with this help text. In my defense, It's late.
	solvent := flag.String("solvent", "h2o", "The solvent for the continuum model used in QM calculations. Only some values are allowed (see code or documentation). 'vac' for in vacuo calculations")
	replicas := flag.Int("replicas", 0, "Number of replicas in a replica-exchange MD simulation, if performed. If less or equal zero, Bartender will come up with a reasonable number")
	maxtemp := flag.Float64("maxtemp", 400.0, "The maximum temperature allowed for a replica in a replica-exchange simulation, if performed")
	exfreq := flag.Int("exfreq", 2, "The frequency of attempted replica exchanges in a replica-exchange simulation, if performed")
	re := flag.Bool("re", false, "Perform a Replica-exchange MD simulation, instead of a regular MD")
	maxiters := flag.Int("maxiters", 0, "Maximum histogram-build/fit iterations for each term")
	//slopeTol := flag.Float64("slopeTol", 0.0, "If removing anharmonicities, the fraction of slope change tolerance for defining the harmonic zone")
	//The following doc string is not clear at all
	//anharmoRange := flag.Float64("anharmoRance", 0.2, "If removing anharmonicities, the fraction of the total range of points that will be considered for anharmonicity search, at each corner")
	usedcommand := strings.Join(os.Args, " ")

	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Usage:\n  %s: [flags] geomtry.pdb/.gro/.xyz bartender_input.inp \n\nFlags:\n", os.Args[0])
		flag.PrintDefaults()
	}
	flag.Parse()
	verb = *verbose
	if *verbose2 > *verbose {
		verb = *verbose2
	}
	plotdir = *plotDir
	if plotdir != "." && !*noplot {
		os.Mkdir(plotdir, 0755)
	}
	if !*noplot {
		os.Mkdir(plotdir+"/xvg", 0755)
	}
	if *refit {
		*mdtime = -1
	}
	//just in case.
	if *skip < 1 {
		*skip = 1
	}
	PrintV(1, "Bartender called as:\n", usedcommand, "\n")
	//the angle increments will be in radians
	d2r := chem.Deg2Rad
	increments := map[string]float64{
		"bonds":  *bi,
		"angles": *ai * d2r,
		"reb":    *ai * d2r,
		"dihe":   *di * d2r,
		"improp": *ii * d2r,
	}

	param := map[string][]*bonded{
		"bonds":  make([]*bonded, 0),
		"angles": make([]*bonded, 0),
		"reb":    make([]*bonded, 0),
		"dihe":   make([]*bonded, 0),
		"improp": nil,
	}
	args := flag.Args()
	if len(args) < 2 {
		log.Fatal("Bartender requires at least 2 arguments: The name of a geometry file and the name of a Bartender input file")
	}
	geoname := args[0]
	inpname := args[1]
	fmt.Printf("Use:\n  $BARTENDERPATH/bartender  [FLAGS] geometry_file input_file\n Use \"bartender -help\" to see the available flags\n\n")
	var mol *chem.Molecule
	var err error

	var g *goodEnough
	if *gefile == "" {
		g = DefaultGoodEnough()
	} else {
		g = GoodEnoughFromFile(*gefile)
	}
	LogV(3, "Good enough fitting parameters:", g.String())
	g.MaxIters(*maxiters)
	extension := getExtension(geoname)
	switch extension {
	case "gro":
		mol, err = chem.GroFileRead(geoname)
	case "pdb":
		mol, err = PoorlyMadePDBFileRead(geoname) //This reads normal PDBs and LigParGen "PDBs"
	default:
		mol, err = chem.XYZFileRead(geoname)
		mol.FillIndexes()
	}
	if err != nil {
		panic("Failed to open geometry input file: " + err.Error())
	}
	mol.SetCharge(*charge) //needed for the MD and the partial charges calculation
	wanted, marked := ParseInputGeo(inpname)
	MDEngine := MD
	//So, there are currently 2 ways of calling for a REMD simulation.
	if len(marked) != 0 || *re {
		MDEngine = REMD
	}
	beads, weights := ParseInputBead(inpname)
	vsites := ParseInputVsites(inpname, beads)
	LogV(2, wanted, "weights:", weights, "vsites (len and list):", len(vsites), vsites)
	if len(vsites) > 0 {
		beads = append(beads, vsites...)
	}
	//bonded parameters
	MakePDB(mol.Coords[0], mol, beads)
	MDS := &MDSettings{time: *mdtime, method: *method, temp: *temperature, solvent: *solvent, cpus: *cpus, replicas: *replicas, maxtemp: *maxtemp, exfreq: *exfreq, restart: *restart}

	//Here we run the calculation to get a GFN0/2/FF trajectory, or we read whatever trajectory the user wants to supply
	var mdout chem.Traj
	var datamap map[string][][]float64
	// Let's test putting the Sommelier thing here, so, when you run Sommelier, Bartender doesn't do the other things

	//This module will test the dihedral parameters
	if *sommelier != "" {
		Sommelier(geoname, *sommelier, *charge, beads, wanted, increments, *temperature)
		os.Exit(0)
	}
	var rebcomment []bool //defines if the ith ReB term is commented out, or left to its previous status
	wanted, rebcomment = initialChecks(wanted)
	//skipbonded is checked again some lines ahead. This is because I can't put the goto statement
	//in that check here, as it would jump over some variable declarations.
	if *owntraj == "" {
		mdoutname := MDEngine(mol.Coords[0], mol, MDS) //This will take a while
		_, mdout, err = chem.XYZFileAsTraj(mdoutname)
		if err != nil {
			panic("Opening of XTB trajectory failed: " + err.Error())
		}
	} else {
		mdout, err = OpenTraj(*owntraj)
		if err != nil {
			panic("Failed to open trajectory: " + err.Error())
		}
	}
	datamap = TrajAn(mdout, mol, beads, wanted, *skip)
	//select which functions will be used, Go or Python
	//HookeFit := GoHookeFit
	F := &Funcs{}
	F.HookeFitAv = GoHookeFitAv
	F.SimplePeriodicFit = GoSimplePeriodicFit
	F.RyckBelleFit = GoRyckBelleFit
	F.CosAngleFit = GoCosAngleFit
	F.ReBFit = GoReBFit
	//The following functionality will probably be removed from the final program
	//Most of the new "fancy" stuff for bimodal wells are only in the Go
	//fitting functions.
	if *py {
		//HookeFit = PyHookeFit
		F.SimplePeriodicFit = PySimplePeriodicFit
		F.RyckBelleFit = PyRyckBelleFit
		F.ReBFit = PyReBFit
		F.CosAngleFit = PyCosAngleFit
	}
	//Some checks for each term, like, check that angles have 3 beads, etc.

	PrintV(1, "All Energies in kJ/mol, distances in nm, angles in degrees\n")
	PrintV(2, "Except for the distributions printed in -verbose>=2, where angles are in radians")
	for k, v := range datamap {
		for i, w := range v {
			mean := stat.Mean(w, nil)
			PrintV(2, "\nData for:", k, BeadsText(wanted[k][i]), "Mean value:", mean, "Data points:", len(w))
			PrintV(3, "Increment used:", increments[k])

			O := &StatAndFitOptions{*removeAnharmonic, increments[k], *temperature, *histcutoff, *improperTolerance, *noplot}
			params := StatAndFit(F, k, wanted[k][i], i, w, O, g)
			param[k] = append(param[k], params...)

			if k == "dihe" {
				ia := increments["angles"]
				par3, R23 := ManageBendingTorsion(datamap, wanted, i, *temperature, []float64{increments["dihe"], ia, ia})
				//R23 should never be negative, so we'll use a negative value to signal that the fit was not obtained.
				if R23 >= 0 {
					b := NewBonded(i, wanted[k][i], par3, R23, 11, true)
					param[k] = append(param[k], b)
				} else {
					LogV(3, fmt.Sprintf("Combined bending-torsion potential for beands %s will not be obtained, for lack of bending angles in input", BeadsText(wanted[k][i])))
				}
			}
		}
	}
	PrintBonded(param, "gmx_out.itp", rebcomment, usedcommand)

	if *owntraj == "" && *dcdsave != "" {
		err = DCDSave(*dcdsave, "xtb.trj")
		if err != nil {
			LogV(0, "Couldn't save DCD trajectory: ", err.Error()) //This always gets printed, even in non-verbose mode.
		}
	}
	Done("\nYour Martini, Mr. Bond.")

}

func BeadsText(beads []int) string {
	ret := " "
	for _, v := range beads {
		ret = ret + strconv.Itoa(v+1) + "-"
	}
	return ret[:len(ret)-1]
}

func CategoryName(k string) string {
	category := strings.ReplaceAll(k, "s", "")
	//	category = strings.Title(category)
	if category == "dihe" {
		category = "dihedral"
	} else if category == "improp" {
		category = "improper dihedral"
	}
	return category

}

type bonded struct {
	ID        int
	beads     []int
	params    []float64
	rmsd      float64
	functype  int
	commented bool
}

func NewBonded(ID int, beads []int, params []float64, rmsd float64, functype int, commented bool) *bonded {
	ret := new(bonded)
	ret.beads = beads
	ret.params = params
	ret.rmsd = rmsd
	ret.functype = functype
	ret.commented = commented
	ret.ID = ID
	return ret

}

func (b *bonded) SetComment(c bool) {
	b.commented = c
}

func (b *bonded) Comment() string {
	if b.commented {
		return ";;"
	}
	return ""
}

// Settings for MD. Not all these are
// always needed.
type MDSettings struct {
	time     int
	method   string
	temp     float64
	solvent  string
	cpus     int
	replicas int
	maxtemp  float64
	exfreq   int
	restart  bool
}

func Done(message string) {
	if message == "" {
		message = "\nYour Martini, Mr. Bond."
	}
	fmt.Printf("There were %d warnings, reprinted now for your convenience:\n\n", len(warnings))
	for i, w := range warnings {
		fmt.Printf("%d: %s\n", i, w)
	}
	fmt.Println(message)
	os.Exit(0)

}
