/*
 * fit_python.go, part of Bartender
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
	"os"
	"os/exec"
	"strconv"
	"strings"
)

//we'll implement at least 2 fits. 2-degree poly, for bond/angles and a trigonometric function for  dihedrals. We'll probably implement several alternatives for dihedral angles.

var pyexec string = os.ExpandEnv("${BTROOT}/bartender_fit.py")

func PyCosAngleFit(x, y []float64) ([]float64, float64) {
	guess := cosangleGuess(x, y)
	fout, err := os.Create("task.fit")
	if err != nil {
		panic(err.Error())
	}
	fout.WriteString("cosa\n")
	fout.WriteString(fmt.Sprintln(x))
	fout.WriteString(fmt.Sprintln(y))
	fout.WriteString(fmt.Sprintf("%5.3f %5.3f\n", *guess[0], *guess[1]))
	//bartender_fit.py needs to be in the PATH!
	command := exec.Command("sh", "-c", pyexec)
	err = command.Run()
	if err != nil {
		panic(err.Error())

	}
	par, rmsd, err := FloatFileParse("cosa.out")

	if err != nil {
		panic(err.Error())
	}
	LogV(1, "Python says", par, rmsd)
	return par, rmsd
}

func PyReBFit(x, y []float64) ([]float64, float64) {
	guess := reBGuess(x, y)
	fout, err := os.Create("task.fit")
	if err != nil {
		panic(err.Error())
	}
	fout.WriteString("reb\n")
	fout.WriteString(fmt.Sprintln(x))
	fout.WriteString(fmt.Sprintln(y))
	fout.WriteString(fmt.Sprintf("%5.3f %5.3f\n", *guess[0], *guess[1]))
	//bartender_fit.py needs to be in the PATH!
	command := exec.Command("sh", "-c", pyexec)
	err = command.Run()
	if err != nil {
		panic(err.Error())

	}
	par, rmsd, err := FloatFileParse("reb.out")

	if err != nil {
		panic(err.Error())
	}
	LogV(1, "Python says", par, rmsd)
	return par, rmsd

}

func PyHookeFit(x, y []float64) ([]float64, float64) {
	guess := hookeGuess(x, y)
	fout, err := os.Create("task.fit")
	if err != nil {
		panic(err.Error())
	}
	fout.WriteString("hooke\n")
	fout.WriteString(fmt.Sprintln(x))
	fout.WriteString(fmt.Sprintln(y))
	fout.WriteString(fmt.Sprintf("%5.3f %5.3f\n", *guess[0], *guess[1]))
	//bartender_fit.py needs to be in the PATH!
	command := exec.Command("sh", "-c", pyexec)
	err = command.Run()
	if err != nil {
		panic(err.Error())

	}
	par, rmsd, err := FloatFileParse("hooke.out")

	if err != nil {
		panic(err.Error())
	}
	LogV(1, "Python says", par, rmsd)
	return par, rmsd

}

func PySimplePeriodicFit(x, y []float64) ([]float64, float64) {
	guess := simplePeriodicGuess(x, y)
	fout, err := os.Create("task.fit")
	if err != nil {
		panic(err.Error())
	}
	fout.WriteString("simple_periodic\n")
	fout.WriteString(fmt.Sprintln(x))
	fout.WriteString(fmt.Sprintln(y))
	fout.WriteString(fmt.Sprintf("%5.3f %5.3f %5.3f\n", *guess[0], *guess[1], *guess[2]))
	//bartender_fit.py needs to be in the PATH!
	command := exec.Command("sh", "-c", pyexec)
	err = command.Run()
	if err != nil {
		panic(err.Error())

	}
	par, rmsd, err := FloatFileParse("simple_periodic.out")

	if err != nil {
		panic(err.Error())
	}
	LogV(1, "Python says", par, rmsd)
	return par, rmsd

}

func PyRyckBelleFit(x, y []float64) ([]float64, float64) {
	guess := ryckBelleGuess(x, y)
	fout, err := os.Create("task.fit")
	if err != nil {
		panic(err.Error())
	}
	fout.WriteString("ryck_belle\n")
	fout.WriteString(fmt.Sprintln(x))
	fout.WriteString(fmt.Sprintln(y))
	fout.WriteString(fmt.Sprintf("%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", *guess[0], *guess[1], *guess[2], *guess[3], *guess[4], *guess[5]))
	//bartender_fit.py needs to be in the PATH!
	command := exec.Command("sh", "-c", pyexec)
	err = command.Run()
	if err != nil {
		panic(err.Error())

	}
	par, rmsd, err := FloatFileParse("ryck_belle.out")

	if err != nil {
		panic(err.Error())
	}
	LogV(1, "Python says", par, rmsd)
	return par, rmsd

}

// IndexFileParse will read a file which contains one line with integer numbers separated by spaces. It returns those numbers
// as a slice of ints, and an error or nil.
func FloatFileParse(filename string) ([]float64, float64, error) {
	parfile, err := os.Open(filename)
	if err != nil {
		return nil, -1, err
	}
	defer parfile.Close()
	indexes := bufio.NewReader(parfile)
	line, err := indexes.ReadString('\n')
	if err != nil {
		return nil, -1, err
	}
	ret1, err := floatStringParse(line)
	if err != nil {
		return nil, -1, err
	}
	line, err = indexes.ReadString('\n')
	if err != nil {
		return nil, -1, err
	}
	ret2, err := floatStringParse(line)

	return ret1, ret2[0], nil

}

func floatStringParse(str string) ([]float64, error) {
	var err error
	fields := strings.Fields(str)
	ret := make([]float64, len(fields))
	for key, val := range fields {
		ret[key], err = strconv.ParseFloat(val, 64)
		if err != nil {
			return nil, err
		}
	}
	return ret, nil
}
