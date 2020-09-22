package sasa

import (
	"os/exec"
	"respdb/pdb"
	"strconv"
	"strings"
)

// BuriedResidues returns a list of buried residues (relative sidechain SASA < 50%)
func BuriedResidues(p *pdb.PDB) ([]*pdb.Residue, error) {
	var buried []*pdb.Residue

	cmd := exec.Command("freesasa",
		p.ID+".pdb",
		"--format=rsa")
	cmd.Dir = "bin/"

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

// SASA returns total, apolar and polar SASA for the given PDB path.
func SASA(path string) (totalSASA float64, apolarSASA float64, polarSASA float64, unknownSASA float64, err error) {
	cmd := exec.Command("freesasa",
		"--resolution=100",
		path)

	out, err := cmd.CombinedOutput()
	if err != nil {
		return 0, 0, 0, 0, err
	}

	for _, line := range strings.Split(string(out), "\n") {
		f := strings.Fields(line)
		if len(f) == 3 {
			val, _ := strconv.ParseFloat(f[2], 64)
			switch f[0] {
			case "Total":
				totalSASA = val
			case "Apolar":
				apolarSASA = val
			case "Polar":
				polarSASA = val
			case "Unknown":
				unknownSASA = val
			}
		}
	}

	return totalSASA, apolarSASA, polarSASA, unknownSASA, err
}
