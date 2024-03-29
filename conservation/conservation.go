package conservation

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"io/ioutil"
	"math"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"sync"

	"github.com/tikz/bio/pdb"
	"github.com/tikz/bio/uniprot"
)

// http://www.ebi.ac.uk/uniprot/TrEMBLstats
//                          A     C     D     E     F     G     H     I     K     L     M     N     P     Q     R     S     T     V     W     Y
var abundance = [20]float64{9.06, 1.28, 5.46, 6.23, 3.88, 7.27, 2.22, 5.54, 4.93, 9.88, 2.34, 3.79, 4.97, 3.81, 5.82, 6.79, 5.55, 6.87, 1.30, 2.88}

// Family represents a family from Pfam mapped to a sequence.
type Family struct {
	ID       string
	HMM      *HMM
	Mappings []*Mapping
}

// HMM represents a parsed model from a .hmm file
type HMM struct {
	Name         string
	Desc         string
	ConsensusAas []string
	MatchEms     [][20]float64
	MatchPs      [][20]float64
	Bitscores    []float64
	Entropies    []float64
}

// Mapping holds equivalent positions between sequence and HMM.
type Mapping struct {
	Position      int
	PositionModel int
	Bitscore      float64
	Entropy       float64
}

// Pfam represents a storage of HMM models.
type Pfam struct {
	modelsDir string
	mutex     sync.Mutex
}

// NewPfam instantiates the path for storing HMM models.
func NewPfam(modelsPath string) (*Pfam, error) {
	modelsPath = filepath.Clean(modelsPath)
	f, err := os.Stat(modelsPath)
	if os.IsNotExist(err) {
		return nil, err
	}
	if !f.IsDir() {
		return nil, fmt.Errorf("%s is not a dir", modelsPath)
	}

	return &Pfam{modelsDir: modelsPath, mutex: sync.Mutex{}}, nil
}

// Families creates a slice of Families from a UniProt.
func (pfam *Pfam) Families(unp *uniprot.UniProt) (fams []*Family, err error) {
	pfam.mutex.Lock()

	// Write temporary FASTA of sequence
	fastaPath := os.TempDir() + "/" + unp.ID + ".fasta"
	err = ioutil.WriteFile(fastaPath, []byte(">"+unp.ID+"\n"+unp.Sequence), 0644)
	if err != nil {
		return nil, fmt.Errorf("write FASTA: %v", err)
	}

	defer func() {
		os.RemoveAll(fastaPath)
		pfam.mutex.Unlock()
	}()

	for _, id := range unp.Pfam {
		var fam Family
		fam.ID = id

		hmmPath := pfam.modelsDir + "/" + id + ".hmm"

		// Get and parse HMM
		hmm, err := pfam.getHMM(id)
		if err != nil {
			return nil, fmt.Errorf("parse hmm file: %v", err)
		}
		if hmm == nil {
			continue
		}
		fam.HMM = hmm

		// Align using hmmalign
		mappings, err := align(hmmPath, fastaPath)
		if err != nil {
			return nil, fmt.Errorf("align: %v", err)
		}

		for _, m := range mappings {
			fam.Mappings = append(fam.Mappings, &Mapping{
				Position:      m[0] + 1,
				PositionModel: m[1] + 1,
				Bitscore:      hmm.Bitscores[m[1]],
				Entropy:       hmm.Entropies[m[1]],
			})
		}

		fams = append(fams, &fam)
	}

	return fams, nil
}

// ResiduesBitscore returns a map of residues to HMM bitscores.
func (pfam *Pfam) ResiduesBitscore(unpID string, fams []*Family, p *pdb.PDB) map[*pdb.Residue]float64 {
	residues := make(map[*pdb.Residue]float64)
	for _, fam := range fams {
		for _, m := range fam.Mappings {
			for _, res := range p.UniProtPositions[unpID][int64(m.Position)] {
				residues[res] = m.Bitscore
			}
		}
	}
	return residues
}

// getHMM downloads a HMM model from Pfam.
func (pfam *Pfam) getHMM(id string) (*HMM, error) {
	hmmPath := pfam.modelsDir + "/" + id + ".hmm"

	hmm, err := loadHMM(hmmPath)
	if err == nil {
		return hmm, nil
	}

	resp, err := http.Get("http://pfam.xfam.org/family/" + id + "/hmm")
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	if resp.StatusCode != 200 {
		// The file will contain HTML describing the error instead
		// of an actual HMMER3 model. This happens with dead families.
		return nil, nil
	}

	out, err := os.Create(hmmPath)
	if err != nil {
		return nil, err
	}
	defer out.Close()

	io.Copy(out, resp.Body)

	return loadHMM(hmmPath)
}

func posBitscore(matchEms [20]float64) (sum float64) {
	for i, m := range matchEms {
		sum += math.Exp(-m) * math.Log2(math.Exp(-m)/(abundance[i]/100))
	}
	return sum
}

func posEntropy(matchEms [20]float64) (sum float64) {
	for _, m := range matchEms {
		sum += math.Exp(-m) * math.Log2(math.Exp(-m))
	}
	return -sum
}

func posPs(matchEms [20]float64) (matchPs [20]float64) {
	for i, m := range matchEms {
		matchPs[i] = math.Exp(-m)
	}
	return matchPs
}

func loadHMM(path string) (*HMM, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	hmm := HMM{}
	line, modelStartLine := 0, 0
	modelStart := false

	// http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())

		if fields[0] == "NAME" {
			hmm.Name = fields[1]
		}

		if fields[0] == "DESC" {
			hmm.Desc = strings.Join(fields[1:], " ")
		}

		if fields[0] == "HMM" {
			// ""the presence of two header lines is mandatory: the parser
			// always skips the line after the HMM tag line.""
			modelStart = true
			modelStartLine = line + 2
		}

		// ""The first line in the main model section -may be- an -optional- line starting with COMPO""
		// ""the last line of the format is the “//” record separator.""
		fieldOk := fields[0] != "COMPO" && fields[0] != "//"

		// 3 lines per node: match, insert, state lines
		skip := (line-modelStartLine)%3 == 0

		if modelStart && line >= modelStartLine && fieldOk && skip {
			// ""The next K numbers for match emissions, one per symbol, in alphabetic order.""
			var matchEms [20]float64
			for i, m := range fields[1:21] {
				v, _ := strconv.ParseFloat(m, 64)
				matchEms[i] = v
			}

			hmm.ConsensusAas = append(hmm.ConsensusAas, strings.ToUpper(fields[22]))
			hmm.MatchEms = append(hmm.MatchEms, matchEms)
			hmm.MatchPs = append(hmm.MatchPs, posPs(matchEms))
			hmm.Bitscores = append(hmm.Bitscores, posBitscore(matchEms))
			hmm.Entropies = append(hmm.Entropies, posEntropy(matchEms))
		}
		line++
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return &hmm, nil
}

func align(hmmPath string, fastaPath string) (mappings [][2]int, err error) {
	cmd := exec.Command("hmmalign", "--outformat", "A2M", hmmPath, fastaPath)
	out, err := cmd.CombinedOutput()
	if err != nil {
		return nil, errors.New(string(out))
	}

	alignSeq := parseFASTA(string(out))

	modelIndex, seqIndex, alignSeqIndex := 0, 0, 0

	for alignSeqIndex < len(alignSeq) {
		alignChar := string(alignSeq[alignSeqIndex])

		switch {
		case alignChar == "-": // in model, not in seq
			modelIndex++
		case alignChar == strings.ToLower(alignChar): // in seq, not in model
			seqIndex++
		default: // match
			mappings = append(mappings, [2]int{seqIndex, modelIndex})
			modelIndex++
			seqIndex++
		}

		alignSeqIndex++
	}

	return mappings, nil
}

func parseFASTA(txt string) (sequence string) {
	for _, line := range strings.Split(string(txt), "\n") {
		if len(line) > 0 && string(line[0]) != ">" {
			sequence += line
		}
	}
	return sequence
}
