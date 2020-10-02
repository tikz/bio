package pdb

import (
	"errors"
	"fmt"
	"log"
	"math"
	"regexp"
	"strings"
)

var residueNames = [...][3]string{
	{"Alanine", "Ala", "A"},
	{"Arginine", "Arg", "R"},
	{"Asparagine", "Asn", "N"},
	{"Aspartic acid", "Asp", "D"},
	{"Cysteine", "Cys", "C"},
	{"Glutamic acid", "Glu", "E"},
	{"Glutamine", "Gln", "Q"},
	{"Glycine", "Gly", "G"},
	{"Histidine", "His", "H"},
	{"Isoleucine", "Ile", "I"},
	{"Leucine", "Leu", "L"},
	{"Lysine", "Lys", "K"},
	{"Methionine", "Met", "M"},
	{"Phenylalanine", "Phe", "F"},
	{"Proline", "Pro", "P"},
	{"Serine", "Ser", "S"},
	{"Threonine", "Thr", "T"},
	{"Tryptophan", "Trp", "W"},
	{"Tyrosine", "Tyr", "Y"},
	{"Valine", "Val", "V"},
}

// Residue represents a single residue from the PDB structure.
type Residue struct {
	Chain           string  `json:"chain"`
	StructPosition  int64   `json:"structPosition"`
	Position        int64   `json:"position"`
	UnpID           string  `json:"unp_id"`
	UnpPosition     int64   `json:"unp_position"`
	Name            string  `json:"-"`
	Name1           string  `json:"name1"`
	Name3           string  `json:"-"`
	Atoms           []*Atom `json:"-"`
	MeanBFactor     float64 `json:"mean_bfactor"`
	NormMeanBFactor float64 `json:"norm_mean_bfactor"`
}

// IsAminoacid returns true if the given letter is an aminoacid, false otherwise.
func IsAminoacid(letter string) bool {
	for _, res := range residueNames {
		if res[2] == letter {
			return true
		}
	}
	return false
}

// AminoacidNames receives a name and returns a 3-sized array of all the possible representations as a string.
func AminoacidNames(input string) (string, string, string) {
	s := strings.Title(strings.ToLower(input))
	for _, res := range residueNames {
		for _, n := range res {
			if n == s {
				return res[0], res[1], res[2]
			}
		}
	}

	return input, "Unk", "X"
}

// NewResidue constructs a new residue given a chain, position and aminoacid name.
// The name is case-insensitive and can be either a full aminoacid name, one or three letter abbreviation.
func NewResidue(chain string, pos int64, input string) *Residue {
	name, abbrv3, abbrv1 := AminoacidNames(input)

	res := &Residue{
		Chain:          chain,
		StructPosition: pos,
		Name:           name,
		Name1:          abbrv1,
		Name3:          abbrv3,
	}

	return res
}

// ExtractSeqRes parses the raw PDB for SEQRES records containing the primary sequence.
func (pdb *PDB) ExtractSeqRes(rawPDB []byte) error {
	regex, _ := regexp.Compile("SEQRES[ ]*.*?[ ]+(.*?)[ ]+([0-9]*)[ ]*([A-Z ]*)")
	matches := regex.FindAllStringSubmatch(string(rawPDB), -1)
	if len(matches) == 0 {
		return errors.New("SEQRES not found")
	}

	pdb.SeqRes = make(map[string][]*Residue)
	for _, match := range matches {
		chain := match[1]
		resSplit := strings.Split(match[3], " ")
		for i, resStr := range resSplit {
			if resStr != "" {
				res := NewResidue(chain, int64(i), resStr)
				pdb.SeqRes[chain] = append(pdb.SeqRes[chain], res)
			}
		}
	}

	return nil
}

// ExtractResidues extracts data from the ATOM and HETATM records and parses them.
func (pdb *PDB) ExtractResidues(rawPDB []byte) error {
	atoms, err := pdb.extractPDBATMRecords(rawPDB, "ATOM")
	if err != nil {
		return fmt.Errorf("extract ATOM records: %v", err)
	}

	hetatms, _ := pdb.extractPDBATMRecords(rawPDB, "HETATM")

	pdb.Atoms = atoms
	pdb.HetAtoms = hetatms

	err = pdb.ExtractPDBChains()
	if err != nil {
		return fmt.Errorf("extract PDB chains: %v", err)
	}

	return nil
}

// ExtractPDBChains parses the residue chains.
func (pdb *PDB) ExtractPDBChains() error {
	atoms := pdb.Atoms
	if len(atoms) == 0 {
		return errors.New("empty atoms list")
	}

	chains := make(map[string]map[int64]*Residue)

	var res *Residue
	for _, atom := range atoms {
		chain, chainOk := chains[atom.Chain]
		pos, posOk := chain[atom.ResidueNumber]

		if !chainOk {
			chains[atom.Chain] = make(map[int64]*Residue)
		}
		if !posOk {
			res = NewResidue(atom.Chain, atom.ResidueNumber, atom.Residue)
			res.Atoms = []*Atom{atom}
			res.calculateMeanBFactor()
			chains[atom.Chain][atom.ResidueNumber] = res
		} else {
			pos.Atoms = append(pos.Atoms, atom)
		}
	}

	pdb.Chains = chains
	for _, chain := range pdb.Chains {
		pdb.TotalLength += int64(len(chain))
	}

	pdb.calculateNormMeanBFactor()

	return nil
}

// calculateMeanBFactor calculates the mean B-factor for the residue based on all its atoms.
func (r *Residue) calculateMeanBFactor() {
	var sum float64
	for _, atom := range r.Atoms {
		sum += atom.BFactor
	}

	if len(r.Atoms) == 0 {
		log.Fatal(r.Chain, r.StructPosition)
	}
	r.MeanBFactor = sum / float64(len(r.Atoms))
}

func stddev(vals []float64, mean float64) float64 {
	var ss float64
	for _, v := range vals {
		ss += math.Pow(float64(v)-mean, 2)
	}
	return math.Pow(ss/float64(len(vals)), 0.5)
}

// calculateNormMeanBFactor calculates the z-score for residues mean B-factors.
func (pdb *PDB) calculateNormMeanBFactor() {
	var sum float64
	var n float64
	var meanBfactors []float64
	for _, residues := range pdb.Chains {
		for _, residue := range residues {
			sum += residue.MeanBFactor
			meanBfactors = append(meanBfactors, residue.MeanBFactor)
			n++
		}
	}
	mean := sum / n
	s := stddev(meanBfactors, mean)

	for _, residues := range pdb.Chains {
		for _, residue := range residues {
			residue.NormMeanBFactor = (residue.MeanBFactor - mean) / s
		}
	}
}
