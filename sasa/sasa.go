package sasa

import (
	"math"
	"os/exec"
	"strconv"
	"strings"

	"github.com/tikz/bio/pdb"
)

type SASA struct {
	Total        float64
	Side         float64
	Main         float64
	Apolar       float64
	Polar        float64
	Unknown      float64
	Residues     []*pdb.Residue
	ResiduesSASA []ResidueSASA
}

// ResidueSASA represents results for a single residue, a line in the output.
type ResidueSASA struct {
	All       float64
	RelAll    float64
	Side      float64
	RelSide   float64
	Main      float64
	RelMain   float64
	Apolar    float64
	RelApolar float64
	Polar     float64
	RelPolar  float64
}

// Run runs FreeSASA on a PDB and parses the output.
func Run(p *pdb.PDB) (sasa SASA, err error) {
	cmd := exec.Command("freesasa",
		p.PDBPath,
		"--format=rsa")

	out, err := cmd.CombinedOutput()
	if err != nil {
		return
	}

	for _, l := range strings.Split(string(out), "\n") {
		if len(l) > 0 && l[0:3] == "RES" {
			// aa := l[4:7]
			chain := string(l[8])
			pos, _ := strconv.ParseInt(strings.TrimSpace(l[9:13]), 10, 64)

			parse := func(start int, end int) float64 {
				s := strings.TrimSpace(l[start:end])
				if s == "N/A" {
					return math.NaN()
				}
				v, _ := strconv.ParseFloat(s, 64)
				return v
			}

			sasa.Residues = append(sasa.Residues, p.Chains[chain][pos])
			sasa.ResiduesSASA = append(sasa.ResiduesSASA, ResidueSASA{
				All:       parse(13, 22),
				RelAll:    parse(22, 29),
				Side:      parse(29, 35),
				RelSide:   parse(35, 41),
				Main:      parse(41, 48),
				RelMain:   parse(48, 54),
				Apolar:    parse(54, 61),
				RelApolar: parse(61, 67),
				Polar:     parse(67, 74),
				RelPolar:  parse(74, 80),
			})
		}

		if len(l) > 0 && l[0:5] == "TOTAL" {
			parse := func(s string) float64 {
				v, _ := strconv.ParseFloat(s, 64)
				return v
			}
			f := strings.Fields(l)
			sasa.Total = parse(f[0])
			sasa.Side = parse(f[1])
			sasa.Main = parse(f[2])
			sasa.Apolar = parse(f[3])
			sasa.Polar = parse(f[4])
		}
	}

	return
}
