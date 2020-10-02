package dssp

import (
	"os/exec"
	"strconv"
	"strings"

	"github.com/tikz/bio/pdb"
)

type DSSPResults struct {
	Residues map[*pdb.Residue]string
}

// DSSP calculates secondary structure for a given PDB.
func DSSP(p *pdb.PDB) (results DSSPResults, err error) {
	cmd := exec.Command("mkdssp", "-i", p.PDBPath)

	out, err := cmd.CombinedOutput()
	if err != nil {
		return
	}

	results.Residues = make(map[*pdb.Residue]string)

	// https://swift.cmbi.umcn.nl/gv/dssp/
	start := false
	for _, l := range strings.Split(string(out), "\n") {
		if len(l) > 17 {
			if start {
				posStr := strings.TrimSpace(l[5:10])
				if len(posStr) > 0 {
					chain := string(l[11])
					pos, _ := strconv.ParseInt(posStr, 10, 64)
					if chain, ok := p.Chains[chain]; ok {
						if res, ok := chain[pos]; ok {
							results.Residues[res] = string(l[16])
						}
					}
				}
			}
			if string(l[2]) == "#" {
				start = true
			}
		}
	}

	return
}
