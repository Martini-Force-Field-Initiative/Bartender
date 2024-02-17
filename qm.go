/*
 * qm.go, part of Bartender
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
	"log"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/qm"
	v3 "github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
)

// uses the ree and xtb program to run a replica-exchange simulation
func REMD(coord *v3.Matrix, mol chem.AtomMultiCharger, MD *MDSettings) string {
	if MD.method == "" {
		MD.method = "gfn0" //the default
	}
	if MD.temp == 0 {
		MD.temp = 298.0
	}
	if MD.time <= 0 {
		MD.time = 500
	}
	if MD.cpus < 0 {
		MD.cpus = runtime.NumCPU()
	}
	charge := mol.Charge()
	multi := mol.Multi()
	if MD.replicas <= 0 {
		MD.replicas = MD.cpus / 3
		if MD.replicas < 1 {
			MD.replicas = 5
		}
	}
	inpgeo := "toree.pdb"
	err := chem.PDBFileWrite(inpgeo, coord, mol, nil)
	if err != nil {
		panic("Failed to write input pdb from REMD: " + err.Error())
	}
	var reexec string = os.ExpandEnv("${BTROOT}/RE/ree")
	command := fmt.Sprintf("%s -maxt %5.3f -method %s -cpus %d -charge %d -multi %d -tref %5.3f -dielectric %3.1f -replicas %d %s %d > ree.log 2>1", reexec, MD.maxtemp, MD.method, MD.cpus, charge, multi, MD.temp, solvent2eps(MD.solvent), MD.replicas, inpgeo, MD.time)
	LogV(2, command) /////////////////////////
	remd := exec.Command("sh", "-c", command)
	err = remd.Run()
	if err != nil {
		panic("Call to the REMD utility failed: " + err.Error())
	}
	return "fulltraj.xyz"

}

func do_restart() {
	mol, traj, err := chem.XYZFileAsTraj("xtb.trj")
	if err != nil {
		log.Printf("Couldn't perform restart: %s", err.Error())
		return
	}
	prev := v3.Zeros(traj.Len())
	curr := v3.Zeros(traj.Len())
	cont := 0
	//whatever happens, the traj file will always be closed by
	//when this loop ends. traj will always close itself upon error
	//or normal termination.
	for err := traj.Next(curr); err == nil; err = traj.Next(curr) {
		prev.Copy(curr)
		cont++
	}
	if cont < 2 {
		log.Printf("Too few frames to restart: %d", cont)
	}
	fs2ps := 1.0 / 1000
	fsbetweenframes := 50.0
	psdone := int(math.Floor(float64(cont) * fsbetweenframes * fs2ps))
	if psdone < 1 {
		log.Printf("Restarts won't be performed for existing trajectories of less than 1 ps. Got: %d", psdone)
		return
	}
	//Now we need to edit the input file to subtract the already run time from the total time
	//wanted.
	inp, err := scu.NewMustReadFile("gochem.inp")
	if err != nil {
		log.Printf("Couldn't perform restart: Couldn't open gochem.inp: %s", err.Error())
		return
	}
	out, err := os.Create("gochem.inp-tmp")
	if err != nil {
		log.Printf("Couldn't perform restart: Couldn't open temp file for gochem.inp: %s", err.Error())
		return
	}
	for i := inp.Next(); i != "EOF"; i = inp.Next() {
		if strings.Contains(i, "time") {
			i = strings.Split(i, "#")[0]
			t := strings.Split(i, "=")[1] //could panic
			t = strings.Replace(t, "\n", "", -1)
			tn, err := strconv.Atoi(t)
			if err != nil {
				log.Printf("Couldn't perform restart: Failed to parse time from line: %s, Error: %s", i, err.Error())
				return
			}
			remaining := tn - psdone
			i = fmt.Sprintf(" time=%d\n", remaining)

		}
		out.WriteString(i)
	}
	//Only now we start overwriting files.
	os.Rename("xtb.trj", "xtb.trj-PREV")
	os.Rename("gochem.inp-tmp", "gochem.inp")
	chem.XYZFileWrite("gochem.xyz", prev, mol)
	return
}

// This has been tested and seems to work fine.
// note that the "notused" parameter is only there to keep the same signature as REMD
func MD(coord *v3.Matrix, mol chem.AtomMultiCharger, MD *MDSettings) string {
	var dry bool = false
	Q := new(qm.Calc)
	if MD.method == "" {
		MD.method = "gfn0" //the default
	}
	if MD.temp == 0 {
		MD.temp = 298.0
	}
	if MD.time <= 0 {
		dry = true
		MD.time = 500 //0.5 ns, rather demanding
	}
	if MD.cpus < 0 {
		MD.cpus = runtime.NumCPU()
	}
	Q.Method = MD.method
	Q.Dielectric = solvent2eps(MD.solvent)
	Q.Job = &qm.Job{MD: true}
	Q.MDTime = MD.time //simulation time (whatever unit the program uses!) it's ps for xtb
	Q.MDTemp = MD.temp
	xtb := qm.NewXTBHandle()
	xtb.SetnCPU(MD.cpus)
	err := xtb.BuildInput(coord, mol, Q)
	if err != nil {
		panic("Failed to build XTB input for the QM-MD" + err.Error()) //you know the drill
	}
	//Here we deal with restarts.
	postProcTraj := func() { return } //this function is always called after the MD ends, but it does nothing if we aren't doing a restart.
	if MD.restart {
		do_restart()
		postProcTraj = func() {
			//Bartender is unix-only anyway, and I'm lazy.
			//NOTE: Write properly,  maybe.
			cq := exec.Command("sh", "-c", "cat xtb.trj-PREV xtb.trj > xtb.trj-full")
			cq.Run()
			os.Rename("xtb.trj-full", "xtb.trj")
			os.Remove("xtb.trj-PREV")
		}
	}

	if !dry {
		err = xtb.Run(true) //we wait for the simulation to end, this will take a while!
	}
	//We will try to remove scoord files left by xtb, but if it doesn't work, it doesn't work
	//the program will just keep running.
	toremove, err := filepath.Glob("scoord*")
	if err == nil {
		for _, f := range toremove {
			_ = os.Remove(f)
		}
	}
	postProcTraj()
	//I need to check that the thing ran correctly. There is probably some file produced a the
	//end by xtb for which I can check.
	if f, err := os.Open("xtb.trj"); err == nil {
		f.Close()        //just wanted to know if it was there. Surely there is a more direct method, but I'm lazy.
		return "xtb.trj" //the name of the resulting trajectory. It's just an multi-xyz file
	} else {
		panic("Failed produce a readable XTB QM-MD trajectory " + err.Error())
	}
}

func solvent2eps(sol string) float64 {
	def := 80.0
	var solvent2eps = map[string]float64{
		"h2o":          80.0,
		"chcl3":        5.0,
		"ch2cl2":       9.0,
		"acetone":      21.0,
		"acetonitrile": 37.0,
		"methanol":     33.0,
		"toluene":      2.0,
		"thf":          7.0,
		"dmso":         47.0,
		"dmf":          38.0,
		"vac":          -1,
	}
	ret, ok := solvent2eps[sol]
	if !ok {
		log.Printf("Solvent %s not found. Will use the default dielectric, %4.2f\n", sol, def)
		return def
	}
	return ret

}
