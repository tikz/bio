package uniprot

import (
	"encoding/json"
	"fmt"
	"strings"

	"github.com/tikz/bio/http"
)

type SIFTSBestStructures struct {
	PDBID      string  `json:"pdb_id"`
	Resolution float64 `json:"resolution"`
	Method     string  `json:"experimental_method"`
	Coverage   float64 `json:"coverage"`
}

// getSIFTSBestStructures retrieves the available structures for a given UniProt ID.
func getSIFTSBestStructures(unpID string) ([]PDB, error) {
	raw, _ := http.Get("https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/" + unpID)
	unps := make(map[string]json.RawMessage)
	err := json.Unmarshal(raw, &unps)
	if err != nil { // Empty JSON, no crystals
		return []PDB{}, nil
	}

	var structures []SIFTSBestStructures
	err = json.Unmarshal(unps[unpID], &structures)
	if err != nil {
		return nil, fmt.Errorf("unmarshal UniProt keys: %v", err)
	}

	uniquePDBIDs := make(map[string]struct{})
	var uniqueStructures []PDB

	for _, s := range structures {
		if _, ok := uniquePDBIDs[s.PDBID]; !ok && s.Method == "X-ray diffraction" {
			uniqueStructures = append(uniqueStructures, PDB{
				ID:         strings.ToUpper(s.PDBID),
				Coverage:   s.Coverage,
				Method:     s.Method,
				Resolution: s.Resolution,
			})
			uniquePDBIDs[s.PDBID] = struct{}{}
		}
	}

	return uniqueStructures, nil
}
