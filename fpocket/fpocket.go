package fpocket

import (
	"fmt"
	"io/ioutil"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"github.com/tikz/bio/pdb"
)

// Results represents the results of a Fpocket job.
type Results struct {
	Pockets []Pocket
	Dir     string
}

// Pocket represents a single pocket from a pocketN_atm.pdb file.
type Pocket struct {
	Path      string
	DrugScore float64
	Residues  []*pdb.Residue
}

// Run runs Fpocket on a PDB and parses the results
func Run(outPath string, p *pdb.PDB) (Results, error) {
	outPath = filepath.Clean(outPath)
	outJobPath := outPath + "/" + p.ID + "_out"

	var results Results

	_, err := os.Stat(outJobPath)
	if os.IsNotExist(err) {
		// Copy PDB to output base dir since fpocket always
		// creates the results dir in the same level as the PDB.
		tmpPDB, err := p.CopyPDB(outPath)
		if err != nil {
			return results, err
		}

		// Run Fpocket
		cmd := exec.Command("fpocket", "-f", tmpPDB)

		out, err := cmd.CombinedOutput()
		if err != nil || strings.Contains(string(out), "failed") {
			return results, fmt.Errorf("fpocket: %v %s", err, string(out))
		}

		os.Remove(tmpPDB)
	}

	// Walk created folder containing pocket analysis files
	pockets, err := walkPocketDir(p, outJobPath+"/pockets")
	if err != nil {
		return results, nil
	}

	results.Pockets = pockets
	results.Dir = outJobPath
	return results, nil
}

// ResidueInPocket returns if the given residue is in a pocket, and the pocket's drug score.
func (r *Results) ResidueInPocket(res *pdb.Residue) (inPocket bool, drugScore float64) {
	for _, pocket := range r.Pockets {
		for _, pRes := range pocket.Residues {
			if pRes == res {
				inPocket, drugScore = true, pocket.DrugScore
				return
			}
		}
	}
	return
}

func walkPocketDir(p *pdb.PDB, dir string) (pockets []Pocket, err error) {
	n := 0
	err = filepath.Walk(dir, func(file string, info os.FileInfo, err error) error {
		// For each Fpocket result PDB file
		if strings.Contains(file, "_atm.pdb") {
			data, err := ioutil.ReadFile(file)
			if err != nil {
				return err
			}
			// Extract drug score
			regexScore, _ := regexp.Compile("Drug Score.*: ([0-9.]*)")
			drugScore, err := strconv.ParseFloat(regexScore.FindAllStringSubmatch(string(data), -1)[0][1], 64)
			if err != nil {
				return err
			}

			pocketPDB, err := pdb.NewPDBFromRaw(data)
			if err != nil {
				return err
			}

			var residues []*pdb.Residue
			for chain, chainPos := range pocketPDB.Chains {
				for pos := range chainPos {
					residues = append(residues, p.Chains[chain][pos])
				}
			}

			pocket := Pocket{
				Path:      dir + file,
				DrugScore: drugScore,
				Residues:  residues,
			}
			pockets = append(pockets, pocket)
			n++
		}

		return nil
	})

	return pockets, err
}
