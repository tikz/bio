package abswitch

import (
	"crypto/sha256"
	"encoding/hex"
	"errors"
	"fmt"
	"io/ioutil"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/tikz/bio/pdb"
	"github.com/tikz/bio/uniprot"
)

type AbSwitch struct {
	binPath     string
	resultsPath string
}

type ResidueResults struct {
	Gor string  `json:"gor"`
	PaH float64 `json:"pah"`
	PaE float64 `json:"pae"`
	PaC float64 `json:"pac"`
	Amb float64 `json:"amb"`
	Ins float64 `json:"ins"`
	Swi float64 `json:"swi"`
	S5s float64 `json:"s5s"`
}

func generateID(seq string) string {
	hash := sha256.Sum256([]byte(seq))
	return hex.EncodeToString(hash[:])
}

func NewAbSwitch(binPath string, resultsPath string) (ab *AbSwitch, err error) {
	ab = &AbSwitch{}
	if ab.binPath, err = filepath.Abs(binPath); err != nil {
		return nil, err
	}

	if ab.resultsPath, err = filepath.Abs(resultsPath); err != nil {
		return nil, err
	}

	return ab, nil
}

func (ab *AbSwitch) Switchability(u *uniprot.UniProt, p *pdb.PDB) (results map[int64]*ResidueResults, err error) {
	results = make(map[int64]*ResidueResults)

	for _, chain := range p.SIFTS.UniProt[u.ID].Mappings {
		seq := u.Sequence[chain.UnpStart:chain.UnpEnd]
		residues, err := ab.Run(seq)
		if err != nil {
			return nil, err
		}

		for i, res := range residues {
			results[chain.UnpStart+int64(i)] = res
		}
	}

	return results, nil
}

func (ab *AbSwitch) parse(path string) ([]*ResidueResults, error) {
	raw, err := ioutil.ReadFile(path)
	if err != nil {
		return nil, err
	}

	var results []*ResidueResults

	for _, line := range strings.Split(string(raw), "\n")[1:] {
		fields := strings.Fields(line)
		if len(fields) >= 9 {
			toFloat := func(str string) float64 {
				f, _ := strconv.ParseFloat(str, 64)
				return f
			}
			results = append(results, &ResidueResults{
				Gor: fields[1],
				PaH: toFloat(fields[2]),
				PaE: toFloat(fields[3]),
				PaC: toFloat(fields[4]),
				Amb: toFloat(fields[5]),
				Ins: toFloat(fields[6]),
				Swi: toFloat(fields[7]),
				S5s: toFloat(fields[8]),
			})
		}
	}

	return results, nil
}

func (ab *AbSwitch) Run(seq string) ([]*ResidueResults, error) {
	name := generateID(seq)
	outFile := ab.resultsPath + "/" + name + ".s5"

	if _, err := os.Stat(outFile); os.IsNotExist(err) {
		// Calculate
		dirName := filepath.Dir(ab.binPath)

		// Write fasta
		fastaFile := "abswitch_" + name + ".fasta"
		ioutil.WriteFile(dirName+"/"+fastaFile, []byte(">"+name+"\n"+seq), 0644)

		// Write cfg
		cfgFile := "abswitch_" + name + ".cfg"
		cfg := fmt.Sprintf("command=Switch5\nfasta=%s\noFile=%s", fastaFile, outFile)
		ioutil.WriteFile(dirName+"/"+cfgFile, []byte(cfg), 0644)

		defer func() {
			os.RemoveAll("bin/" + fastaFile)
			os.RemoveAll("bin/" + cfgFile)
		}()

		// Run
		cmd := exec.Command(ab.binPath, "-f", cfgFile)
		cmd.Dir = filepath.Dir(ab.binPath)
		out, err := cmd.CombinedOutput()
		strOut := string(out)
		if err != nil || !strings.Contains(strings.ToLower(strOut), "printed results") {
			return nil, errors.New(strOut)
		}
	}

	return ab.parse(outFile)
}
