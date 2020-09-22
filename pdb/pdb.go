package pdb

import (
	"fmt"
	"io/ioutil"
	"respdb/http" // TOOD: decouple this
	"time"
)

// PDB represents a single PDB entry.
type PDB struct {
	ID     string `json:"id"`     // PDB ID
	URL    string `json:"url"`    // RCSB web page URL
	PDBURL string `json:"pdbUrl"` // RCSB download URL for the PDB file
	CIFURL string `json:"cifUrl"` // RCSB download URL for the CIF file

	Title       string     `json:"title"`       // publication title
	Date        *time.Time `json:"date"`        // publication date
	Method      string     `json:"method"`      // experimental method used
	Resolution  float64    `json:"resolution"`  // method resolution
	TotalLength int64      `json:"totalLength"` // total length as sum of residues of all chains in the structure

	Atoms     []*Atom  `json:"-"`         // ATOM records in the structure
	HetAtoms  []*Atom  `json:"-"`         // HETATM records in the structure
	HetGroups []string `json:"hetGroups"` // HET groups in the structure

	// Position mapping
	SIFTS               *SIFTS           // EBI SIFTS data for residue position mapping
	SeqResOffsets       map[string]int64 `json:"seqResOffsets"`       // Chain ID to SEQRES position offsets
	ChainStartResNumber map[string]int64 `json:"chainStartResNumber"` // Chain ID to First residue number as informed in ATOM column.
	ChainEndResNumber   map[string]int64 `json:"chainEndResNumber"`   // Chain ID to Last residue number as informed in ATOM column.

	SeqRes           map[string][]*Residue           `json:"-"`      // PDB SEQRES chain ID to residue pointers
	SeqResChains     map[string]map[int64]*Residue   `json:"-"`      // PDB SEQRES chain ID and PDB ATOM position to residue in structure
	Chains           map[string]map[int64]*Residue   `json:"chains"` // PDB ATOM chain ID and position to pointer in structure
	UniProtPositions map[string]map[int64][]*Residue `json:"-"`      // UniProt ID to sequence position to residue(s) (multiple chains) in structure

	// Extra data
	// SITE records
	BindingSite map[string][]*Residue `json:"bindingSite"` // binding site identifier to residues compromising it

	// REMARK 800 site descriptions
	BindingSiteDesc map[string]string `json:"bindingSiteDesc"` // binding site identifier to description

	RawPDB []byte `json:"-"` // PDB file raw data
	RawCIF []byte `json:"-"` // CIF file raw data

	LocalPath string `json:"-"` // local path for the PDB file
}

// NewPDBFromID constructs a new instance from a UniProt accession ID and PDB ID, fetching and parsing the data.
func NewPDBFromID(pdbID string) (PDB, error) {
	pdb := PDB{ID: pdbID}

	err := pdb.Load()
	return pdb, err
}

// NewPDBFromRaw constructs a new instance from raw bytes, and only extracts ATOM records.
// This is useful for parsing PDB output files generated from external tools.
func NewPDBFromRaw(raw []byte) (*PDB, error) {
	pdb := PDB{RawPDB: raw}

	err := pdb.ExtractResidues()
	if err != nil {
		return nil, fmt.Errorf("parse: %v", err)
	}

	return &pdb, nil
}

// Load fetches and parses the necessary data.
func (pdb *PDB) Load() error {
	err := pdb.Fetch()
	if err != nil {
		return fmt.Errorf("fetch data: %v", err)
	}

	err = pdb.Parse()
	if err != nil {
		return fmt.Errorf("parse: %v", err)
	}

	return nil
}

// Parse parses the raw PDB text.
func (pdb *PDB) Parse() error {
	err := pdb.Extract()
	if err != nil {
		return fmt.Errorf("extract: %v", err)
	}

	pdb.makeMappings()

	pdb.extractSites()
	return nil
}

// Fetch downloads all external data for the entry.
func (pdb *PDB) Fetch() error {
	url := "https://www.rcsb.org/structure/" + pdb.ID
	urlCIF := "https://files.rcsb.org/download/" + pdb.ID + ".cif"
	rawCIF, err := http.Get(urlCIF)
	if err != nil {
		return fmt.Errorf("download CIF file: %v", err)
	}

	urlPDB := "https://files.rcsb.org/download/" + pdb.ID + ".pdb"
	rawPDB, err := http.Get(urlPDB)
	if err != nil {
		return fmt.Errorf("download PDB file: %v", err)
	}

	pdb.URL = url
	pdb.PDBURL = urlPDB
	pdb.CIFURL = urlCIF
	pdb.RawPDB = rawPDB
	pdb.RawCIF = rawCIF

	err = pdb.getSIFTSMappings()
	if err != nil {
		return fmt.Errorf("SIFTS: %v", err)
	}

	return nil
}

// Extract parses data from the raw PDB, raw CIF, SIFTS, and populates the entry.
func (pdb *PDB) Extract() error {
	err := pdb.ExtractSeqRes()
	if err != nil {
		return fmt.Errorf("extract SEQRES: %v", err)
	}

	err = pdb.ExtractResidues()
	if err != nil {
		return fmt.Errorf("extract PDB residues: %v", err)
	}

	err = pdb.ExtractCIFData()
	if err != nil {
		return fmt.Errorf("extract CIF data: %v", err)
	}

	return nil
}

// WriteFile writes the raw PDB contents to a file.
func (pdb *PDB) WriteFile(path string) error {
	err := ioutil.WriteFile(path, pdb.RawPDB, 0644)
	if err != nil {
		return fmt.Errorf("write PDB file: %v", err)
	}

	pdb.LocalPath = path
	return nil
}
