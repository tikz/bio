package tango

import (
	"crypto/sha256"
	"encoding/hex"
	"errors"
	"io/ioutil"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/tikz/bio/pdb"
	"github.com/tikz/bio/uniprot"
)

type Tango struct {
	binPath     string
	resultsPath string
}

type ResidueResults struct {
	// http://tango.crg.es/Tango_Handbook.pdf
	// Those files will have the following columns:
	Beta             float64 // Percentage of beta-strand conformation
	BetaTurn         float64 // Percentage of beta-turn conformation
	Helix            float64 // Percentage of alpha-helical conformation
	Aggregation      float64 // Percentage of Aggregation
	HelixAggregation float64 // Percentage of Helical Aggregation
}

func generateID(seq string) string {
	hash := sha256.Sum256([]byte(seq))
	return hex.EncodeToString(hash[:])
}

func NewTango(binPath string, resultsPath string) (t *Tango, err error) {
	t = &Tango{}
	if t.binPath, err = filepath.Abs(binPath); err != nil {
		return nil, err
	}

	if t.resultsPath, err = filepath.Abs(resultsPath); err != nil {
		return nil, err
	}

	return t, nil
}

func (t *Tango) Aggregability(u *uniprot.UniProt, p *pdb.PDB) (results map[int64]*ResidueResults, err error) {
	results = make(map[int64]*ResidueResults)

	for _, chain := range p.SIFTS.UniProt[u.ID].Mappings {
		seq := u.Sequence[chain.UnpStart:chain.UnpEnd]
		residues, err := t.Run(p.ID+"-"+chain.ChainID, seq)
		if err != nil {
			return nil, err
		}

		for i, res := range residues {
			results[chain.UnpStart+int64(i)] = res
		}
	}

	return results, nil
}

func (t *Tango) Run(name string, seq string) ([]*ResidueResults, error) {

	// http://tango.crg.es/Tango_Handbook.pdf
	// The format of the sequences to be run is as follows:
	// Name Cter Nter pH Temp Ionic Sequence
	// Name = name of the sequence (less than 25 characters)
	// Cter = status of the C-terminus of the peptide (amidated Y, free N)
	// Nter = status of the N-terminus of the peptide (acetylated A, succinilated S and free N)
	// pH = pH
	// Temp = temperature in Kelvin
	// Ionic = ionic strength in M
	// sequence = sequence of the peptide in one letter code.

	outPath := t.resultsPath + "/" + name + ".txt"
	if _, err := os.Stat(outPath); os.IsNotExist(err) {
		binDir, binFile := filepath.Split(t.binPath)
		relResultsPath, err := filepath.Rel(binDir, t.resultsPath)
		if err != nil {
			return nil, err
		}

		cmd := exec.Command("./"+binFile, relResultsPath+"/"+name, "ct=N", "nt=N", "ph=7.4", "te=303", "io=0.05", "seq="+seq)
		cmd.Dir = binDir
		out, err := cmd.CombinedOutput()
		strOut := string(out)
		if err != nil || strings.Contains(strings.ToLower(strOut), "error") {
			return nil, errors.New(strOut)
		}
	}

	return t.parse(outPath)
}

func (t *Tango) parse(path string) ([]*ResidueResults, error) {
	raw, err := ioutil.ReadFile(path)
	if err != nil {
		return nil, err
	}

	var results []*ResidueResults

	for _, line := range strings.Split(string(raw), "\n")[1:] {
		fields := strings.Fields(line)
		if len(fields) >= 7 {
			toFloat := func(str string) float64 {
				f, _ := strconv.ParseFloat(str, 64)
				return f
			}
			results = append(results, &ResidueResults{
				Beta:             toFloat(fields[2]),
				BetaTurn:         toFloat(fields[3]),
				Helix:            toFloat(fields[4]),
				Aggregation:      toFloat(fields[5]),
				HelixAggregation: toFloat(fields[6]),
			})
		}
	}

	return results, nil
}
