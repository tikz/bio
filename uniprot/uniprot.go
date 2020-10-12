package uniprot

import (
	"errors"
	"fmt"
	"regexp"
	"strconv"
	"strings"

	"github.com/tikz/bio/http"
)

// UniProt contains relevant protein data for a single accession.
type UniProt struct {
	ID           string                 `json:"id"`           // accession ID
	URL          string                 `json:"url"`          // page URL for the entry
	TXTURL       string                 `json:"txtUrl"`       // TXT API URL for the entry.
	Name         string                 `json:"name"`         // protein name
	Gene         string                 `json:"gene"`         // gene code
	Organism     string                 `json:"organism"`     // organism
	Sequence     string                 `json:"sequence"`     // canonical sequence
	PDBs         []PDB                  `json:"pdbs"`         // PDBs
	Sites        []Site                 `json:"sites"`        // protein function sites
	PTMs         PTMs                   `json:"ptms"`         // post translational modifications
	Pfam         []string               `json:"pfam"`         // Pfam families accessions
	Variants     []VariantEntry         `json:"variants"`     // variants
	Publications map[string]Publication `json:"publications"` // PubMed ID to publications
	Raw          []byte                 `json:"-"`            // TXT API raw bytes.
}

// PDB represents a single available PDB structure for an UniProt.
type PDB struct {
	ID         string  `json:"id"`
	Coverage   float64 `json:"coverage"`
	Method     string  `json:"method"`
	Resolution float64 `json:"resolution"`
}

// VariantEntry represents a single variant entry extracted from the TXT.
type VariantEntry struct {
	Position  int64    `json:"position"`
	FromAa    string   `json:"fromAa"`
	ToAa      string   `json:"toAa"`
	Change    string   `json:"change"`
	Note      string   `json:"note"`
	Evidence  string   `json:"evidence"`
	PubMedIDs []string `json:"pubmedIds"`
	ID        string   `json:"id"`
	DbSNP     string   `json:"dbsnp"`
}

// Publication represents a single publication entry extracted from the TXT.
type Publication struct {
	Title   string `json:"title"`
	Authors string `json:"authors"`
	Journal string `json:"journal"`
	PubMed  string `json:"pubmed"`
	DOI     string `json:"doi"`
}

// SAS represents a single aminoacid substitution.
type SAS struct {
	Position int64  `json:"position"`
	FromAa   string `json:"fromAa"`
	ToAa     string `json:"toAa"`
}

// PTMs represents post translational modifications from the entry.
type PTMs struct {
	DisulfideBonds   []Disulfide       `json:"disulfideBonds"`
	Glycosilations   []Glycosilation   `json:"glycosilationSites"`
	ModifiedResidues []ModifiedResidue `json:"modifiedResidues"`
}

// Disulfide represents a single disulfide bond between positions.
type Disulfide struct {
	Positions [2]int64 `json:"positions"`
}

// Glycosilation represents a glycosilation site.
type Glycosilation struct {
	Position int64  `json:"position"`
	Note     string `json:"note"`
}

// ModifiedResidue represents a modified residue.
type ModifiedResidue struct {
	Position int64  `json:"position"`
	Note     string `json:"note"`
}

// Site represents either a binding, active or metal site entry.
type Site struct {
	Type     string `json:"type"`
	Position int64  `json:"position"`
	Note     string `json:"note"`
}

// NewUniProt constructs an instance from an UniProt accession ID and a list of target PDB IDs
func NewUniProt(uniprotID string) (*UniProt, error) {
	url := "https://www.uniprot.org/uniprot/" + uniprotID
	txtURL := url + ".txt"
	raw, err := http.Get(txtURL)
	if err != nil {
		return nil, fmt.Errorf("get UniProt accession %v: %v", uniprotID, err)
	}

	u := &UniProt{
		ID:     uniprotID,
		URL:    url,
		TXTURL: txtURL,
		Raw:    raw,
	}

	// Parse UniProt TXT
	err = u.extract()
	if err != nil {
		return nil, fmt.Errorf("extract PDB crystals %v: %v", uniprotID, err)
	}

	return u, nil
}

// extract parses the TXT response.
func (u *UniProt) extract() error {
	err := u.extractSequence()
	if err != nil {
		return fmt.Errorf("get seq: %v", err)
	}

	err = u.extractPDBs()
	if err != nil {
		return fmt.Errorf("extracting crystals from UniProt TXT: %v", err)
	}

	err = u.extractNames()
	if err != nil {
		return fmt.Errorf("extracting names from UniProt TXT: %v", err)
	}

	err = u.extractVariants()
	if err != nil {
		return fmt.Errorf("extracting variants from UniProt TXT: %v", err)
	}

	err = u.extractPublications()
	if err != nil {
		return fmt.Errorf("extracting publications from UniProt TXT: %v", err)
	}

	err = u.extractPTMs()
	if err != nil {
		return fmt.Errorf("extracting PTMs from UniProt TXT: %v", err)
	}

	err = u.extractFams()
	if err != nil {
		return fmt.Errorf("extracting families from UniProt TXT: %v", err)
	}

	u.extractSites()

	return nil
}

// extractPDBs parses the TXT for PDB IDs and populates UniProt.PDBs
func (u *UniProt) extractPDBs() error {
	r, _ := regexp.Compile(`(?ms)PDB; (.*?); (.*?); ([\.0-9]*).*?;.*?$`)
	rStartEnd, _ := regexp.Compile("([0-9]*)-([0-9]*)")

	matches := r.FindAllStringSubmatch(string(u.Raw), -1)
	// Parse each PDB match in TXT
	for _, m := range matches {
		coverage := 0.0
		rangeMatches := rStartEnd.FindAllStringSubmatch(m[0], -1)
		for _, rm := range rangeMatches {
			startPos, _ := strconv.Atoi(rm[1])
			endPos, _ := strconv.Atoi(rm[2])
			coverage += float64(endPos - startPos)
		}

		coverage = coverage / float64(len(u.Sequence))
		resolution, _ := strconv.ParseFloat(m[3], 64)
		if m[2] == "X-ray" {
			u.PDBs = append(u.PDBs, PDB{
				ID:         m[1],
				Coverage:   coverage,
				Method:     m[2],
				Resolution: resolution,
			})
		}
	}

	return nil
}

// extractSequence parses the canonical sequence.
func (u *UniProt) extractSequence() error {
	r, _ := regexp.Compile("(?ms)^SQ.*?$(.*?)//") // https://regex101.com/r/ZTOYaJ/1
	matches := r.FindAllStringSubmatch(string(u.Raw), -1)

	if len(matches) == 0 {
		return errors.New("canonical sequence not found")
	}

	seqGroup := matches[0][1]
	sequence := strings.ReplaceAll(seqGroup, " ", "")
	sequence = strings.ReplaceAll(sequence, "\n", "")

	u.Sequence = sequence

	return nil
}

// extractNames parses protein, gene and organism names
func (u *UniProt) extractNames() error {
	r, _ := regexp.Compile("(?m)^DE.*?Name.*?Full=(.*?)(;| {)")
	matches := r.FindAllStringSubmatch(string(u.Raw), -1)

	if len(matches) == 0 {
		return errors.New("protein name not found")
	}
	u.Name = matches[0][1]

	r, _ = regexp.Compile("(?m)^GN.*?=(.*?)[;| ]")
	matches = r.FindAllStringSubmatch(string(u.Raw), -1)

	if len(matches) != 0 {
		u.Gene = matches[0][1]
	}

	r, _ = regexp.Compile("(?m)^OS[ ]+(.*?)\\.")
	matches = r.FindAllStringSubmatch(string(u.Raw), -1)

	if len(matches) == 0 {
		return errors.New("organism name not found")
	}
	u.Organism = matches[0][1]

	return nil
}

// extractVariants parses for variant references
func (u *UniProt) extractVariants() error {
	var variants []VariantEntry

	// https://regex101.com/r/BpJ3QB/1
	r, _ := regexp.Compile("(?ms)^FT[ ]*VARIANT[ ]*([0-9]*)$(.*?)id=\"(.*?)\"")
	rPubmed, _ := regexp.Compile("PubMed:([0-9]*)")
	matches := r.FindAllStringSubmatch(string(u.Raw), -1)

	for _, variant := range matches {
		var entry VariantEntry
		pos, err := strconv.ParseInt(variant[1], 10, 64)
		if err != nil {
			return fmt.Errorf("cannot parse variant position int: %s", variant[1])
		}
		entry.Position = pos

		data := variant[2]
		s := regexp.MustCompile("\nFT \\s+")
		d := s.ReplaceAllString(data, " ")

		r, _ := regexp.Compile("(?s)/note=\"(.*?)\"")
		n := r.FindAllStringSubmatch(d, -1)
		entry.Note = n[0][1]

		r, _ = regexp.Compile("(.) -> (.) ")
		ne := r.FindAllStringSubmatch(entry.Note, -1)
		if len(ne) == 0 {
			continue
		}
		entry.FromAa = ne[0][1]
		entry.ToAa = ne[0][2]
		entry.Change = entry.FromAa + variant[1] + entry.ToAa

		r, _ = regexp.Compile("dbSNP:(rs[0-9]*)")
		nedb := r.FindAllStringSubmatch(entry.Note, -1)
		if len(nedb) > 0 {
			entry.DbSNP = nedb[0][1]
		}

		r, _ = regexp.Compile("(?s)/evidence=\"(.*?)\"")
		e := r.FindAllStringSubmatch(d, -1)
		if len(e) > 0 {
			entry.Evidence = e[0][1]
			pubmedMatches := rPubmed.FindAllStringSubmatch(entry.Evidence, -1)

			for _, match := range pubmedMatches {
				entry.PubMedIDs = append(entry.PubMedIDs, match[1])
			}
		}

		entry.ID = variant[3]
		variants = append(variants, entry)
	}

	u.Variants = variants

	return nil
}

func (u *UniProt) extractPublications() error {
	u.Publications = make(map[string]Publication)

	// https://regex101.com/r/G2FV96/1/
	rPub, _ := regexp.Compile("(?ms)^RN.*?RL.*?$")
	rTitle, _ := regexp.Compile("(?m)^RT[ ]*(.*)")
	rAuthors, _ := regexp.Compile("(?m)^RA[ ]*(.*)")
	rJournal, _ := regexp.Compile("RL[ ]*(.*)")
	rPubmed, _ := regexp.Compile("PubMed=([0-9]*)")
	rDOI, _ := regexp.Compile("DOI=(.*?);")

	clean := func(s string) string {
		s = strings.ReplaceAll(s, ";", "")
		s = strings.ReplaceAll(s, "\"", "")
		s = strings.TrimSpace(s)
		return s
	}

	pubMatches := rPub.FindAllStringSubmatch(string(u.Raw), -1)
	for _, pub := range pubMatches {
		var title string
		var authors string
		var journal string
		var pubmed string
		var doi string

		pubmedMatches := rPubmed.FindAllStringSubmatch(pub[0], -1)
		if len(pubmedMatches) > 0 {
			pubmed = pubmedMatches[0][1]
		}

		doiMatches := rDOI.FindAllStringSubmatch(pub[0], -1)
		if len(doiMatches) > 0 {
			doi = doiMatches[0][1]
		}

		titleMatches := rTitle.FindAllStringSubmatch(pub[0], -1)
		for _, titleLine := range titleMatches {
			title += titleLine[1] + " "
		}
		title = clean(title)

		authorMatches := rAuthors.FindAllStringSubmatch(pub[0], -1)
		for _, authorLine := range authorMatches {
			authors += authorLine[1] + " "
		}
		authors = clean(authors)

		journalMatches := rJournal.FindAllStringSubmatch(pub[0], -1)
		if len(journalMatches) > 0 {
			journal = journalMatches[0][1]
		}

		u.Publications[pubmed] = Publication{
			Title:   title,
			Authors: authors,
			Journal: journal,
			PubMed:  pubmed,
			DOI:     doi,
		}
	}

	return nil
}

// extractPTMs parses for post translational modifications
func (u *UniProt) extractPTMs() error {
	u.PTMs = PTMs{}

	// Glycosilation sites
	r, _ := regexp.Compile("(?ms)^FT[ ]*CARBOHYD[ ]*([0-9]*)$.*?note=\"(.*?)\"")
	matches := r.FindAllStringSubmatch(string(u.Raw), -1)

	for _, glyco := range matches {
		pos, _ := strconv.ParseInt(glyco[1], 10, 64)
		if strings.Contains(glyco[2], "N-linked") {
			u.PTMs.Glycosilations = append(u.PTMs.Glycosilations,
				Glycosilation{Position: pos, Note: glyco[2]})
		}
	}

	// Modified residues
	r, _ = regexp.Compile("(?ms)^FT[ ]*MOD_RES[ ]*([0-9]*)$.*?note=\"(.*?)\"")
	matches = r.FindAllStringSubmatch(string(u.Raw), -1)

	for _, modres := range matches {
		pos, _ := strconv.ParseInt(modres[1], 10, 64)
		u.PTMs.ModifiedResidues = append(u.PTMs.ModifiedResidues,
			ModifiedResidue{Position: pos, Note: modres[2]})
	}

	// Disulfide bonds
	r, _ = regexp.Compile("(?ms)^FT[ ]*DISULFID[ ]*([0-9]*)..([0-9]*)")
	matches = r.FindAllStringSubmatch(string(u.Raw), -1)

	for _, disulf := range matches {
		pos1, _ := strconv.ParseInt(disulf[1], 10, 64)
		pos2, _ := strconv.ParseInt(disulf[2], 10, 64)
		u.PTMs.DisulfideBonds = append(u.PTMs.DisulfideBonds,
			Disulfide{Positions: [2]int64{pos1, pos2}})
	}

	return nil
}

// extractSites parses for sites
func (u *UniProt) extractSites() {
	tags := map[string]string{
		"ACT_SITE": "active",     // https://www.uniprot.org/help/act_site
		"NP_BIND":  "nucleotide", // https://www.uniprot.org/help/np_bind
		"METAL":    "metal",      // https://www.uniprot.org/help/metal
		"BINDING":  "binding",    // https://www.uniprot.org/help/binding
		"SITE":     "site",       // https://www.uniprot.org/help/site
	}

	for tag, name := range tags {
		r, _ := regexp.Compile(`(?ms)^FT[ ]*` + tag + `[ ]*([0-9\.]*)$.*?note=\"(.*?)\"`)
		matches := r.FindAllStringSubmatch(string(u.Raw), -1)

		for _, site := range matches {
			posStr := site[1]

			if strings.Contains(posStr, "..") { // range of positions
				positions := strings.Split(posStr, "..")
				start, _ := strconv.ParseInt(positions[0], 10, 64)
				end, _ := strconv.ParseInt(positions[1], 10, 64)
				for i := start; i <= end; i++ {
					u.Sites = append(u.Sites,
						Site{
							Type:     name,
							Position: i,
							Note:     site[2],
						},
					)
				}
			} else {
				pos, _ := strconv.ParseInt(site[1], 10, 64)
				u.Sites = append(u.Sites,
					Site{
						Type:     name,
						Position: pos,
						Note:     site[2],
					},
				)
			}
		}
	}
}

// extractFams parses for Pfam families.
func (u *UniProt) extractFams() error {
	r, _ := regexp.Compile("DR[ ]*Pfam; (.*?);")
	matches := r.FindAllStringSubmatch(string(u.Raw), -1)

	for _, fam := range matches {
		u.Pfam = append(u.Pfam, fam[1])
	}

	return nil
}
