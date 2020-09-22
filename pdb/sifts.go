package pdb

import (
	"encoding/json"
	"errors"
	"fmt"
	"respdb/http"
	"strings"
)

// Reference: https://www.ebi.ac.uk/pdbe/api/doc/sifts.html

// SIFTS represents a valid response from the SIFTS mapping project.
type SIFTS struct {
	Pfam    map[string]*Family    `json:"Pfam"`
	UniProt map[string]*Accession `json:"UniProt"`
}

// Family represents a Pfam family.
type Family struct {
	Identifier  string
	Description string
	Name        string
	Mappings    []*Mapping
}

// Accession represents an UniProt accession.
type Accession struct {
	Identifier string     `json:"identifier"`
	Mappings   []*Mapping `json:"mappings"`
	Name       string     `json:"name"`
}

// Mapping represents position mappings between the database entry and the specific PDB.
type Mapping struct {
	PDBStart     *Position `json:"start"`
	EntityID     int64     `json:"entity_id"`
	PDBEnd       *Position `json:"end"`
	UnpStart     int64     `json:"unp_start"`
	UnpEnd       int64     `json:"unp_end"`
	ChainID      string    `json:"chain_id"`
	StructAsymID string    `json:"struct_asym_id"`
}

// Position represents the start or end position of a PDB.
type Position struct {
	ResidueNumber int64 `json:"residue_number"`
}

// getSIFTSMappings retrieves the UniProt<->PDB position mappings from the SIFTS project.
func (pdb *PDB) getSIFTSMappings() error {
	pdbID := strings.ToLower(pdb.ID)
	raw, _ := http.Get("https://www.ebi.ac.uk/pdbe/api/mappings/" + pdbID)
	pdbs := make(map[string]json.RawMessage)
	err := json.Unmarshal(raw, &pdbs)
	if err != nil {
		return fmt.Errorf("unmarshal: %v", err)
	}

	databases := make(map[string]json.RawMessage)
	err = json.Unmarshal(pdbs[pdbID], &databases)
	if err != nil {
		return fmt.Errorf("unmarshal databases keys: %v", err)
	}

	sifts := SIFTS{}
	err = json.Unmarshal(pdbs[pdbID], &sifts)
	if err != nil {
		return fmt.Errorf("unmarshal: %v", err)
	}

	// if _, ok := sifts.UniProt[pdb.UniProtID]; !ok {
	// 	return fmt.Errorf("no mappings available between %s and %s", pdb.ID, pdb.UniProtID)
	// }

	pdb.SIFTS = &sifts

	return nil
}

// GetChainMapping returns a SIFTS Mapping when passed an UniProt accession and chain name.
func (s *SIFTS) GetChainMapping(accession string, chain string) (*Mapping, error) {
	if _, ok := s.UniProt[accession]; !ok {
		return nil, errors.New("accession not in SIFTS")
	}

	for _, m := range s.UniProt[accession].Mappings {
		if m.ChainID == chain {
			return m, nil
		}
	}

	return nil, errors.New("chain not in SIFTS mappings")
}
