package pdb

import (
	"errors"
	"fmt"
	"regexp"
	"strconv"
	"strings"
	"time"
)

// ExtractCIFData parses the associated CIF file for the PDB entry.
func (pdb *PDB) ExtractCIFData(rawCIF []byte) error {
	title, err := extractCIFLine("title", "_struct.title", rawCIF)
	if err != nil {
		return err
	}

	method, err := extractCIFLine("method", "_refine.pdbx_refine_id", rawCIF)
	if err != nil {
		return err
	}

	resolutionStr, err := extractCIFLine("resolution", "_refine.ls_d_res_high", rawCIF)
	if err != nil {
		return err
	}
	resolution, err := strconv.ParseFloat(resolutionStr, 64)
	if err != nil {
		return err
	}

	date, err := extractCIFDate(rawCIF)
	if err != nil {
		return err
	}

	pdb.Title = title
	pdb.Method = method
	pdb.Resolution = resolution
	pdb.Date = date

	return nil
}

func extractCIFLine(name string, pattern string, raw []byte) (string, error) {
	r, _ := regexp.Compile("(?s)" + pattern + "[ ]*(.*?)_")
	matches := r.FindAllStringSubmatch(string(raw), -1)
	if len(matches) == 0 {
		return "", errors.New("CIF " + name + " not found")
	}

	match := matches[0][1]
	match = strings.TrimSpace(match)
	match = strings.Replace(match, "'", "", -1)

	return match, nil
}

// extractCIFDate parses the main publication date from the CIF file.
func extractCIFDate(raw []byte) (*time.Time, error) {
	r, _ := regexp.Compile("_pdbx_database_status.recvd_initial_deposition_date[ ]*([0-9]*-[0-9]*-[0-9]*)")
	matches := r.FindAllStringSubmatch(string(raw), -1)
	if len(matches) == 0 {
		return nil, errors.New("CIF date not found")
	}

	t, err := time.Parse("2006-01-02", string(matches[0][1]))
	if err != nil {
		return nil, fmt.Errorf("parse CIF date: %v", err)
	}

	return &t, nil
}
