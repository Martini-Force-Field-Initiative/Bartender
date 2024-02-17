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

//This will be a separate set of functions called "Sommelier"
//It will not be "independent" at source code level. I need several Bartender functions
//and I don't think it's worth it to put those functions in a separate library for this to
//use.
//I think it _will_ be independent regarding interface. Once control is transferred to the Sommelier function
//that function will act like "main", i.e., it will not propagate errors upwards, rather, it well fix with the ones it can
//(possibly none!) and just exit the program with an error message for the others. It will also not return control to main
//until everything here is done, and it will not return any value for main to handle.

package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"os/exec"
	"regexp"
	"strconv"
	"strings"
	str "strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
)

// yeah, yeah, it's ugly. If it makes you feel better, it's supposed to be temporary,
// while I'm still figuring out which commands do I need.
func runcq(command string, a ...interface{}) {
	w := exec.Command("sh", "-c", fmt.Sprintf(command, a...))
	w.Run()
}

func OpenGeoFile(geoname string, charge ...int) (*chem.Molecule, error) {
	var mol *chem.Molecule
	var err error
	extension := strings.ToLower(strings.Split(geoname, ".")[1])
	switch extension {
	case "gro":
		mol, err = chem.GroFileRead(geoname)
	case "pdb":
		mol, err = chem.PDBFileRead(geoname)
	default:
		mol, err = chem.XYZFileRead(geoname)
	}
	if err != nil {
		return nil, fmt.Errorf("Failed to open geometry input file: " + err.Error())
	}
	if len(charge) > 0 {
		mol.SetCharge(charge[0])
	}
	return mol, nil

}

//The following is to be removed and placed in a separate "utils" module

type topHeader struct {
	any       *regexp.Regexp
	vsitesany *regexp.Regexp
	spec      map[string]*regexp.Regexp
	set       bool
	line      string
}

func NewTopHeader() *topHeader {
	R := new(topHeader)
	R.Set()
	return R

}

func (T *topHeader) Set() {
	T.any = regexp.MustCompile(`\[\p{Zs}*.*\p{Zs}*\]`)
	T.vsitesany = regexp.MustCompile(`\[\p{Zs}*virtual_sites[01234n]?\p{Zs}*\]`)

	T.spec = map[string]*regexp.Regexp{
		"atoms":      regexp.MustCompile(`\[\p{Zs}*atoms\p{Zs}*\]`),
		"bonds":      regexp.MustCompile(`\[\p{Zs}*bonds\p{Zs}*\]`),
		"vsites1":    regexp.MustCompile(`\[\p{Zs}*virtual_sites1\p{Zs}*\]`),
		"vsites2":    regexp.MustCompile(`\[\p{Zs}*virtual_sites2\p{Zs}*\]`),
		"vsites3":    regexp.MustCompile(`\[\p{Zs}*virtual_sites3\p{Zs}*\]`),
		"vsitesn":    regexp.MustCompile(`\[\p{Zs}*virtual_sitesn\p{Zs}*\]`),
		"exclusions": regexp.MustCompile(`\[\p{Zs}*exclusions\p{Zs}*\]`),
		"molecules":  regexp.MustCompile(`\[\p{Zs}*molecules\p{Zs}*\]`),
		"dihedrals":  regexp.MustCompile(`\[\p{Zs}*dihedrals\p{Zs}*\]`),
	}
	T.set = true

}

func (T *topHeader) delcomments(line string) string {
	if str.HasPrefix(line, ";") {
		return ""
	}
	return str.Split(line, ";")[0] //remove comments
}

// Returns true if the line is a Gromacs header. It discards comments.
func (T *topHeader) Is(line string) bool {
	line = T.delcomments(line)
	return T.any.MatchString(line)
}

func (T *topHeader) IsVirtualSites(line string) bool {
	line = T.delcomments(line)
	return T.vsitesany.MatchString(line)
}

func (T *topHeader) IsDihedrals(line string) bool {
	line = T.delcomments(line)
	return T.spec["dihedrals"].MatchString(line)
}

// Returns a string indicating which Gromacs top file header
// the line is, or an empty string if the line is not a header.
func (T *topHeader) Which(line string) string {
	line = T.delcomments(line)
	if !T.any.MatchString(line) {
		return ""
	}
	for k, v := range T.spec {
		if v.MatchString(line) {
			return k
		}

	}
	return ""

}

const (
	cggro        = "sommelier.gro"
	cgtop        = "sommelier.top"
	bttop        = "gmx_out.itp"
	finaltrjname = "srunwholefit.xtc"
)

type btphaseandn struct {
	atoms []int
	phase float64
	k     float64
	n     float64
}

func (b *btphaseandn) IsPerm(t []int) bool {
	if len(t) != 4 {
		return false
	}
	for _, v := range t {
		if !scu.IsInInt(v, b.atoms) {
			return false
		}
	}
	return true
}

type btpars []*btphaseandn

func (b btpars) Find(t []int) int {
	for i, v := range b {
		if v.IsPerm(t) {
			return i
		}
	}
	return -1
}

func (b btpars) String() string {
	ret := ""
	for _, v := range b {
		ret = fmt.Sprintf("%s [", ret)
		for _, w := range v.atoms {
			ret = fmt.Sprintf("%s%d ", ret, w)
		}
		ret = fmt.Sprintf("%s ]", ret)

	}
	return ret
}

func btdihedralpars() btpars {
	top, err := scu.NewMustReadFile(bttop)
	if err != nil {
		panic(fmt.Sprint("Couldn't read Bartender parameters file", bttop))
	}
	defer top.Close()
	T := NewTopHeader()
	reading := false
	ret := make([]*btphaseandn, 0, 1)
	for fl := top.Next(); fl != "EOF"; fl = top.Next() {
		fl = strings.TrimSpace(fl)
		if strings.HasPrefix(fl, ";") {
			continue
		}
		if reading && T.Is(fl) {
			break
		}
		if !reading && T.IsDihedrals(fl) {
			reading = true
			continue
		}
		if reading {
			fi := strings.Fields(fl)
			if fi[4] != "1" {
				continue //We only care about the "function 1" lines, i.e. the simple periodic parameters
			}
			phn := new(btphaseandn)
			ats := make([]int, 4)
			for i, v := range fi[:4] {
				ats[i], err = strconv.Atoi(v)
				if err != nil {
					panic(fmt.Sprintf("Couldn't read dihedral line %s field %d, %s: %s", fl, i, v, err.Error()))
				}
				ats[i]-- //We set the atom indexes to 0-based, which BT uses internally.
			}
			phn.atoms = ats
			phn.phase, err = strconv.ParseFloat(fi[5], 64)
			if err != nil {
				panic(fmt.Sprintf("Couldn't read phase from dihedral line %s field 5, %s: %s", fl, fi[5], err.Error()))
			}
			phn.k, err = strconv.ParseFloat(fi[6], 64)
			if err != nil {
				panic(fmt.Sprintf("Couldn't read force constant from dihedral line %s field 7, %s: %s", fl, fi[7], err.Error()))
			}

			phn.n, err = strconv.ParseFloat(fi[7], 64)
			if err != nil {
				panic(fmt.Sprintf("Couldn't read multip from dihedral line %s field 7, %s: %s", fl, fi[7], err.Error()))
			}
			ret = append(ret, phn)
		}

	}
	return ret

}

// Tests the force constant of dihedral parameters for a Bartender-generated Martini3 molecule
// Using the given topn topology file (which can't contain any bonded parameters), water as a solvent
func Sommelier(atomisticgeofile, topn string, charge int, beadfuncs []BeadCoord, wanted map[string][][]int, increments map[string]float64, temperature float64) {
	fmt.Printf("\nThis is the Sommelier module of Bartender. Dihedral parameters will be tested in a (by default) short MD.\n")
	fmt.Printf("An installation of Gromacs and libxdrfile is required. Bartender needs to be compiled with xdrfile library support.\n")
	fmt.Printf("Sommelier uses Gromacs and  Insane. Please cite 10.1021/acs.jctc.5b00209 and the corresponding Gromacs Papers.\n\n")
	mol, err := OpenGeoFile(atomisticgeofile)
	if err != nil {
		panic(err.Error())
	}
	molname, err := gochem2cggro(mol.Coords[0], mol, beadfuncs, topn) //creates the gro file
	if err != nil {
		panic(err.Error())
	}
	cgmol, err := chem.GroFileRead(cggro)
	if err != nil {
		panic("Unable to read CG molecule " + err.Error())
	}
	//	fmt.Println("lens", mol.Len(), cgmol.Len()) /////////////////////
	//Now for the analysis. We need to produce this "dummy" functions for TrajAn to use.
	BTPars := btdihedralpars()
	cfuncs := make([]BeadCoord, 0, cgmol.Len())
	//	fmt.Println(cgmol.Len())
	for i := 0; i < cgmol.Len(); i++ {
		g := func(index int) BeadCoord {
			f := func(coord *v3.Matrix, mol chem.Atomer) (*v3.Matrix, []int) {
				if coord == nil && mol == nil {
					return nil, []int{i}
				}
				VecIndex := index
				view := coord.VecView(VecIndex)
				return view, []int{VecIndex}
			}
			return f
		}
		cfuncs = append(cfuncs, g(i))
	}
	_, err = makeTop(topn, molname)
	//_ = itpname ///////////////////////////////
	if err != nil {
		panic("makeTop: " + err.Error())
	}
	err = Insane(charge)
	if err != nil {
		panic("Insane: " + err.Error())
	}
	const filename string = "GMXCOMMANDS"

	DoMDFromFile(filename)
	//we might need to do something about the periodic boundary conditions
	//before continueing, but we'll see.
	mdout, err := OpenTraj(finaltrjname)
	if err != nil {
		panic("Sommelier: Failed to open trajectory: " + err.Error())
	}
	cldatamap := TrajAn(mdout, cgmol, cfuncs, wanted, 1)
	//We have to re do the fitting. At least for now, Sommelier will only check dihedral parameters.
	clparadihe := make([]*bonded, 0, 0)
	for k, v := range cldatamap {
		for i, w := range v {
			if k != "dihe" {
				continue
			}
			removeAnharmonic := false
			points, E := IBoltzmann(w, increments[k], temperature, removeAnharmonic, 0) //doesn't return anything for now, but prints intermediate data.
			btp := BTPars.Find(wanted[k][i])
			var par []float64
			var R2 float64
			if btp == -1 {
				log.Printf("Term not found in the %s file!. Wanted: %v. Got: %s", bttop, wanted[k][i], BTPars.String())
				par, R2 = GoSimplePeriodicFit(points, E)

				if par[0] < 0 {
					par[0] = 2*math.Pi + par[0]
					par[2] = math.Abs(par[2])
				}
			} else {
				bt := BTPars[btp]
				tpar, tr := SommelierSimplePeriodicFit(points, E, bt.phase, bt.n, bt.k)
				par = []float64{bt.phase, tpar[0], bt.n}
				R2 = tr //a bit silly

			}
			b := NewBonded(i, wanted[k][i], par, R2, 1, false)

			clparadihe = append(clparadihe, b)
			par2, R22 := GoRyckBelleFit(points, E)
			b = NewBonded(i, wanted[k][i], par2, R22, 3, false)

			clparadihe = append(clparadihe, b) //[len(param[k])-1] = append(param[k][len(param[k])-1], par2...) //just one after the other
			ia := increments["angles"]
			par3, R23 := ManageBendingTorsion(cldatamap, wanted, i, temperature, []float64{increments["dihe"], ia, ia})
			//R23 should never be negative, so we'll use a negative value to signal that the fit was not obtained.
			if R23 >= 0 {
				//Ill add something to the log later -_-
				b = NewBonded(i, wanted[k][i], par3, R23, 11, true)
				clparadihe = append(clparadihe, b)
			} else {
				clparadihe = append(clparadihe, nil)

			}

		}
	}
	finalpars := make(map[string][]*bonded)
	finalpars["dihe"] = clparadihe
	parname := "sommelier.itp"
	fmt.Println(finalpars)
	PrintBonded(finalpars, parname, nil)

	fmt.Printf("Fruity notes on the nose, almond aftertaste.\n")
	fmt.Printf("Sommelier run completed.\n")

}

// atsfromtop reads a gromacs top/itp file for a single molecule and uses the atom/residue names in that file
// to create and return a list of chem.Atom
func atsfromtop(topn string) ([]*chem.Atom, error) {
	H := NewTopHeader()
	top, err := scu.NewMustReadFile(topn)
	if err != nil {
		return nil, err //I might improve the whole error thing later
	}
	atoms := make([]*chem.Atom, 0, 5)
	var readingatoms bool
	for fl := top.Next(); fl != "EOF"; fl = top.Next() {
		if str.HasPrefix(fl, ";") {
			continue
		}
		l := str.Split(fl, ";")[0] //remove comments
		if readingatoms && H.Is(l) {
			break
		}
		if H.Which(l) == "atoms" {
			readingatoms = true
			continue
		}
		if !readingatoms {
			continue
		}
		f := str.Fields(l)
		if len(f) == 0 {
			break
		}
		at := new(chem.Atom)
		at.MolName = f[3]
		at.MolID = 0   //it's supposed to be only one molecule
		at.Name = f[4] //not 100% sure about this. There are disagreements in some of Riccardo's examples.
		atoms = append(atoms, at)
	}
	return atoms, nil

}

// Builds a map with al the wanted data in one frame. Returns the name of the created molecule
// taken from the top file
func gochem2cggro(coord *v3.Matrix, mol chem.Atomer, beadfuncs []BeadCoord, topn string) (string, error) {
	beads := v3.Zeros(len(beadfuncs))
	for i, v := range beadfuncs {
		vec, _ := v(coord, mol)
		beads.SetVecs(vec, []int{i})
	}

	//need to get topology info from the original top file
	//then write the pdb/gro.
	ats, err := atsfromtop(topn)
	if err != nil {
		return "", fmt.Errorf("couldn't obtain the topology for the Gro file: %s\n", err.Error())
	}
	molname := ats[0].MolName
	mol2 := chem.NewTopology(0, 1, ats)
	//the PDB is just for the user to check, so I won't check for errors.
	chem.PDBFileWrite(strings.Replace(cggro, ".gro", ".pdb", 1), beads, mol2, nil)
	return molname, chem.GroFileWrite(cggro, []*v3.Matrix{beads}, mol2)

}

func makeTop(topn, molname string) (string, error) {
	//The following makes the code dependent on this particular version of Martini3. We could
	//also read the names of the files from somewhere, but I think it's easier, at least
	//for now, to just update Bartender when new versions of Martini3 appear.
	var martinidir string = os.ExpandEnv("${BTROOT}/sommelier-files/martini3")
	const martiniff string = "martini_v3.0.0.itp"
	const martinisolv string = "martini_v3.0.0_solvents_v1.itp"
	const martiniions string = "martini_v3.0.0_ions_v1.itp"
	var topmodn string
	H := NewTopHeader()
	if str.HasSuffix(topn, ".top") {
		topmodn = str.Replace(topn, ".top", "_sommelier.top", -1)
	} else {
		topmodn = str.Replace(topn, ".itp", "_sommelier.itp", -1)
	}

	topmod, err := os.Create(topmodn)
	if err != nil {
		return "", err
	}
	defer topmod.Close()
	topold, err := scu.NewMustReadFile(topn)
	if err != nil {
		return "", err
	}
	defer topold.Close()

	topbt, err := scu.NewMustReadFile(bttop)
	if err != nil {
		return "", err
	}
	defer topbt.Close()

	//We will read the first part of the original topology field, up to before the "[ bonds ]" field and copy
	//it to the new top.
	var tl string
	for tl = topold.Next(); tl != "EOF"; tl = topold.Next() {
		if H.Which(tl) == "bonds" { //this assumes that "bonds" comes after "atoms".
			break
		}
		topmod.WriteString(tl)
	}
	//we now copy the bonding data from the gmx_out.itp file
	//written by Bartender
	for bt := topbt.Next(); bt != "EOF"; bt = topbt.Next() {
		topmod.WriteString(bt)
	}
	//Finally, we skip the bonded data from the old topology file,
	//and copy the vsites and exclusions data from there to our new file.
	reading := false
	for ; tl != "EOF"; tl = topold.Next() {
		if H.IsVirtualSites(tl) || H.Which(tl) == "exclusions" {
			topmod.WriteString(tl)
			reading = true
			continue
		}
		if !reading {
			continue
		}
		w := H.Which(tl)
		if w != "" && !(strings.Contains(w, "virtual_site") || w == "exclusions") {
			reading = false
			continue
		}
		topmod.WriteString(tl)
	}
	//We are now ready with the top/itp file for the molecule.

	//We now write the "general" top file, which includes our newly created top/itp file.
	generaltop, err := os.Create(cgtop)
	if err != nil {
		return "", err
	}
	defer generaltop.Close()
	generaltop.WriteString(fmt.Sprintf("#include \"%s/%s\"\n", martinidir, martiniff))
	generaltop.WriteString(fmt.Sprintf("#include \"./%s\"\n", topmodn))
	generaltop.WriteString(fmt.Sprintf("#include \"%s/%s\"\n", martinidir, martinisolv))
	generaltop.WriteString(fmt.Sprintf("#include \"%s/%s\"\n\n", martinidir, martiniions))

	generaltop.WriteString("[ system ]\nOne molecule\n\n")
	generaltop.WriteString(fmt.Sprintf("[ molecules ]\n%s       1\n", molname))

	return topmodn, nil

}

// Insane handles the solvation of the system using Insane
// It uses water and no salt, except for neutralization, if needed.
func Insane(charge int) error {
	var insane string = os.ExpandEnv("${BTROOT}/sommelier-files/insane/insane.py")
	var bkname string = "sommelierbak.gro"
	H := NewTopHeader()
	runcq(fmt.Sprintf("cp %s %s", cggro, bkname))
	//we might need a bigger box
	com := insane + fmt.Sprintf(" -f %s  -o %s -pbc cubic -box 7,7,7 -salt 0.0 -charge %d -sol W -p dummy.top", bkname, cggro, charge)
	runcq(com)
	dum, err := scu.NewBWFile("dummy.top")
	if err != nil {
		return err
	}
	var i string
	var watersline string
	//Insane is not very good at editing the top file, os we just gave it a dummy file
	//to write. We now need to extract the number of added water molecules
	for i, err = dum.PrevLine(); err == nil; i, err = dum.PrevLine() {
		if str.Contains(i, "W") || str.Contains(i, "NA+") || str.Contains(i, "CL-") {
			watersline += i + "\n"
		}
		if H.Which(i) == "molecules" {
			break
		}

	}
	dum.Close()
	top, err := os.OpenFile(cgtop, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		return err
	}
	top.WriteString(watersline)
	top.Close()
	runcq("rm dummy.top") //we can comment this for debugging.
	return nil

}

// DoMD runs a "Sommelier-standard" simulation protocol you can still modify the input files).
// I specified the -ntmpi 1 because otherwise you can get issues when running on GPUs.
// It's always best to use a GMXCOMMANDS file with adjustments to your machine.
func DoMD() {
	var root string = os.ExpandEnv("${BTROOT}/sommelier-files")
	runcq(fmt.Sprintf("gmx grompp -f %s/mdps/em.mdp -c sommelier.gro -p sommelier.top -o sem.tpr > semprev.out 2>&1", root))
	runcq("gmx mdrun -deffnm sem -ntmpi 1 -v > sem.out 2>&1")
	runcq(fmt.Sprintf("gmx grompp -f %s/mdps/eq.mdp -c sem.gro -p sommelier.top -o seq.tpr > seqprev.out 2>&1", root))
	runcq("gmx mdrun -deffnm seq -ntmpi 1 -v > seq.out")
	runcq(fmt.Sprintf("gmx grompp -f %s/mdps/run.mdp -c seq.gro -t seq.cpt -p sommelier.top -o srun.tpr > semdprev.out 2>&1", root))
	runcq("gmx mdrun -deffnm srun -v -ntmpi 1 > srun.out 2>&1")
	runcq("echo '0' | gmx trjconv -f srun.xtc  -o srunwhole.xtc -s sem.tpr  -pbc nojump")
	runcq("echo '2 0' | gmx trjconv -f srunwhole.xtc   -o " + finaltrjname + "  -s sem.tpr    -fit rot+tran")
	//gmx trjconv -f tempskip5whole.xtc  -o WT5skip5Fit.xtc -s npt5.tpr   -fit rot+tran -n wt.ndx

}

// DoMDFomfile attempts to open a file with the Gromacs commands for the MD protocol. If the file name
// corresponds to a non-existent file, DoMDFomFile will call DoMD to perform a "Sommelier-standard"
// MD protocol.
func DoMDFromFile(fn string) {
	var err error
	fin, err := os.Open(fn)
	if err != nil {
		log.Printf("Unable to open file %s. Will attempt to run MD with the default commands/parameters", fn)
		DoMD()
		return
	}
	bfin := bufio.NewReader(fin)
	var line string
	for line, err = bfin.ReadString('\n'); err == nil; line, err = bfin.ReadString('\n') {
		//you can not run the simulations by writing "DRY RUN" as the first line of the file.
		if strings.Contains(line, "DRY RUN") {
			return
		}
		if strings.HasPrefix(line, ";") || strings.HasPrefix(line, "#") {
			continue
		}
		line = strings.Replace(line, "\n", "", 1)
		runcq(line)
	}
	if err.Error() == "EOF" {
		return
	} else {
		log.Printf("Error while reading the file %s: %s Will attempt to run MD with the default commands/parameters", fn, err.Error())
		DoMD()
		return
	}
}
