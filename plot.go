package main

/*
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

import (
	"fmt"
	"image/color"
	"math"
	"os"
	"strings"

	chem "github.com/rmera/gochem"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

//This is mostly repeated code from fit_go. I might refactor things so these functions are called there to eliminate that problem.

func sperf(par []float64) func(float64) float64 {
	eq := par[0]
	k := par[1]
	n := math.Round(par[2]) //This was originally not rounded here
	//which led to discrepancies between the plots/xvgs and the numbers produced
	//from the parameters in the gmx_out.itp file.
	return func(x float64) float64 { return k * (1 + math.Cos(n*x-eq)) }
}

func hookef(par []float64) func(float64) float64 {
	eq := par[0]
	k := par[1]
	return func(x float64) float64 { return 0.5 * k * math.Pow((x-eq), 2.0) }
}

func rebf(par []float64) func(float64) float64 {
	eq := par[0]
	k := par[1]

	return func(x float64) float64 {
		//Here I also changed "eq" for "math.Cos(eq)". This seems to agree with the Gromacs
		//Manual so I think the previous behavior was a bug (see also the comment in the ReB part of the fit_go.go file)
		return 0.5 * k * math.Pow((math.Cos(x)-math.Cos(eq)), 2.0) * (1 / math.Pow(math.Sin(x), 2))
	}
}

func rybef(p []float64) func(float64) float64 {
	ret := func(x float64) float64 {
		psi := x - math.Pi
		cos := math.Cos
		pow := math.Pow
		return p[0] + p[1]*cos(psi) + p[2]*pow(cos(psi), 2) + p[3]*pow(cos(psi), 3) + p[4]*pow(cos(psi), 4) + p[5]*pow(cos(psi), 5)
	}
	return ret
}

func cosanglef(par []float64) func(float64) float64 {
	eq := par[0] * chem.Deg2Rad
	k := par[1]

	return func(x float64) float64 {
		return 0.5 * k * math.Pow((math.Cos(x)-math.Cos(eq)), 2.0)
	}
}

// plots y and f(x) vs x, as points and a line, respectively, unless given true in noplot, in which case, does nothing.
// the plot is saved to a file name.png
func Plot(f func(float64) float64, x, y []float64, name string, noplot bool) error {
	if noplot {
		return nil
	}
	name = strings.ReplaceAll(name, " ", "")
	unit_conv := 1.0
	has := strings.Contains
	if has(name, "ngle") || has(name, "Rycka") || has(name, "Cos") || has(name, "eriodic") || has(name, "mproper") {
		unit_conv = chem.Rad2Deg
	}
	fy := make([]float64, len(x))
	for i, v := range x {
		fy[i] = f(v)
	}
	pointsData := pointsPlot(x, y, unit_conv)
	funcData := pointsPlot(x, fy, unit_conv)
	p, err := plot.New()
	if err != nil {
		return err
	}

	//Here we create an xvg file with the data (unless the file doesn't open)
	//My plan is to add a Python script to plot everything with Matplotlib
	//perhaps including putting everything in one category in the same figure with subplots.
	fout, err := os.Create(plotdir + "/xvg/" + name + ".xvg")
	if err != nil {
		e := err.Error()
		LogV(3, "Failed to create file", name+".xvg", "for plotting:", e, "Will procede without it")
	} else {
		//I don't really know the xvg format, sorry.
		//feel free to change the following line.
		fout.WriteString("@    xaxis  label \"Length or angle (A or deg)\"\n@    yaxis  label \"Energy (kJ/mol)\"\n@TYPE xy\n@ view 0.15, 0.15, 0.75, 0.85\n@ legend on\n@ legend box on\n@ legend loctype view\n@ legend 0.78, 0.8\n@ legend length 2\n@ s0 legend \"Points\"\n@ s1 legend \"Fitted function\"\n")
		for i, v := range x {
			fout.WriteString(fmt.Sprintf("%f %f %f\n", v*unit_conv, y[i], fy[i]))
		}
		fout.Close()
	}
	//back to the Go plotting
	p.Title.Text = name
	p.X.Label.Text = "Angle (deg) / Dist. (A)"
	p.Y.Label.Text = "Energy (kJ/mol)"
	// Draw a grid behind the data
	p.Add(plotter.NewGrid())
	// Make a scatter plotter and set its style.
	s, err := plotter.NewScatter(pointsData)
	if err != nil {
		return err
	}
	s.GlyphStyle.Shape = draw.PyramidGlyph{}
	s.GlyphStyle.Color = color.RGBA{R: 255, A: 255}

	// Make a line plotter and set its style.
	l, err := plotter.NewLine(funcData)
	if err != nil {
		return err
	}
	l.LineStyle.Color = color.RGBA{B: 255, A: 255}
	p.Add(s, l)
	p.Legend.Add("Trajectory", s)
	p.Legend.Add("Fitted function", l)
	// Save the plot to a PNG file.
	if err := p.Save(6*vg.Inch, 6*vg.Inch, plotdir+"/"+name+".png"); err != nil {
		return fmt.Errorf("Failed to create file :%s %s %s", name+".png", "for plotting:", err.Error())
	}
	return nil
}

func pointsPlot(x, y []float64, unit float64) plotter.XYs {
	pts := make(plotter.XYs, len(x))
	for i, v := range x {
		pts[i].X = v * unit
		pts[i].Y = y[i]
	}
	return pts
}

func funcPlot(x []float64, f func(float64) float64, unit float64) plotter.XYs {
	pts := make(plotter.XYs, len(x))
	for i, v := range x {
		pts[i].X = v * unit
		pts[i].Y = f(v)
	}
	return pts
}
