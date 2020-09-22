package pdb

import (
	"errors"
	"regexp"
	"strconv"
	"strings"
)

// Atom represents a single atom in the structure.
// It contains all the columns from an ATOM or HETATM record in a PDB file.
type Atom struct {
	// PDB columns for the ATOM tag
	Number        int64
	Residue       string
	Chain         string
	ResidueNumber int64
	X             float64
	Y             float64
	Z             float64
	Occupancy     float64
	BFactor       float64
	Element       string
	Charge        string
}

// extractPDBATMRecords extracts either ATOM or HETATM records.
func (pdb *PDB) extractPDBATMRecords(recordName string) ([]*Atom, error) {
	var atoms []*Atom

	r, _ := regexp.Compile("(?m)^" + recordName + ".*$")
	matches := r.FindAllString(string(pdb.RawPDB), -1)
	if len(matches) == 0 {
		return atoms, errors.New("atoms not found")
	}

	var lastRes string
	for _, match := range matches {
		var atom Atom

		// https://www.wwpdb.org/documentation/file-format-content/format23/sect9.html#ATOM
		atom.Number, _ = strconv.ParseInt(strings.TrimSpace(match[6:11]), 10, 64)
		atom.Residue = strings.TrimSpace(match[17:20])
		atom.Chain = match[21:22]
		atom.ResidueNumber, _ = strconv.ParseInt(strings.TrimSpace(match[22:26]), 10, 64)
		atom.X, _ = strconv.ParseFloat(strings.TrimSpace(match[30:38]), 64)
		atom.Y, _ = strconv.ParseFloat(strings.TrimSpace(match[38:46]), 64)
		atom.Z, _ = strconv.ParseFloat(strings.TrimSpace(match[46:54]), 64)
		atom.Occupancy, _ = strconv.ParseFloat(strings.TrimSpace(match[54:60]), 64)
		atom.BFactor, _ = strconv.ParseFloat(strings.TrimSpace(match[60:66]), 64)
		atom.Element = strings.TrimSpace(match[76:78])
		atom.Charge = strings.TrimSpace(match[78:80])

		atoms = append(atoms, &atom)

		if recordName == "HETATM" && atom.Residue != lastRes {
			lastRes = atom.Residue
			exists := false
			for _, het := range pdb.HetGroups {
				if het == lastRes {
					exists = true
				}
			}
			if !exists {
				pdb.HetGroups = append(pdb.HetGroups, lastRes)
			}
		}
	}

	return atoms, nil
}
