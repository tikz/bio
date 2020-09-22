package sasa

import (
	"os/exec"
	"strconv"
	"strings"

	"github.com/tikz/bio/pdb"
)

type SASA struct {
	Total      float64
	Side       float64
	Apolar     float64
	Polar      float64
	Unknown    float64
	PerResidue map[*pdb.Residue]ResidueSASA
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

// BuriedResidues returns a list of buried residues (relative sidechain SASA < 50%)
func BuriedResidues(p *pdb.PDB) ([]*pdb.Residue, error) {
	var buried []*pdb.Residue

	cmd := exec.Command("freesasa",
		p.ID+".pdb",
		"--format=rsa")

	out, err := cmd.CombinedOutput()
	if err != nil {
		return nil, err
	}

	for _, l := range strings.Split(string(out), "\n") {
		if len(l) > 0 && l[0:3] == "RES" {
			aa := l[4:7]
			chain := string(l[8])
			pos, _ := strconv.ParseInt(strings.TrimSpace(l[9:13]), 10, 64)
			rsasa, _ := strconv.ParseFloat(strings.TrimSpace(l[35:41]), 64)
			if aa != "GLY" && rsasa < 50 {
				buried = append(buried, p.Chains[chain][pos])
			}
		}
	}

	return buried, nil
}

// RunSASA returns total, apolar and polar SASA for the given PDB path.
func RunSASA(p *pdb.PDB) (sasa SASA, err error) {
	cmd := exec.Command("freesasa",
		p.PDBPath,
		"--format=rsa")
	cmd.Dir = "bin/"

	out, err := cmd.CombinedOutput()
	if err != nil {
		return
	}

	res := make(map[*pdb.Residue]ResidueSASA)

	for _, l := range strings.Split(string(out), "\n") {
		if len(l) > 0 && l[0:3] == "RES" {
			// aa := l[4:7]
			chain := string(l[8])
			pos, _ := strconv.ParseInt(strings.TrimSpace(l[9:13]), 10, 64)
			// rsasa, _ := strconv.ParseFloat(strings.TrimSpace(l[35:41]), 64)
			res[p.Chains[chain][pos]] = ResidueSASA{}
			// if aa != "GLY" && rsasa < 50 {
			// 	buried = append(buried, p.Chains[chain][pos])
			// }
		}
	}

	return
}
