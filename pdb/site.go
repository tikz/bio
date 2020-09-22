package pdb

import (
	"regexp"
	"strconv"
	"strings"
)

func (pdb *PDB) extractSites() error {
	sites := make(map[string][]*Residue)
	r, _ := regexp.Compile("(?m)^SITE.*$")
	siteRecords := r.FindAllString(string(pdb.RawPDB), -1)

	for _, s := range siteRecords {
		// https://www.wwpdb.org/documentation/file-format-content/format23/sect7.html#SITE
		siteName := s[11:14]
		for i := 18; i < 52; i += 11 {
			residueName := s[i : i+3]
			if residueName != "   " {
				residueName = strings.TrimSpace(residueName)
				chain := strings.TrimSpace(string(s[i+4]))
				pos, _ := strconv.ParseInt(strings.TrimSpace(s[i+5:i+10]), 10, 64)
				if res, ok := pdb.Chains[chain][pos]; ok {
					if strings.ToUpper(res.Name3) == residueName {
						sites[siteName] = append(sites[siteName], res)
					}
				}
			}
		}
	}

	pdb.BindingSite = sites

	pdb.extractSitesRemarks()

	return nil
}

func (pdb *PDB) extractSitesRemarks() {
	r, _ := regexp.Compile("(?m)REMARK 800 (.*?): (.*?)$")
	remarks := r.FindAllStringSubmatch(string(pdb.RawPDB), -1)

	siteDescs := make(map[string]string)
	var identifier string
	for _, r := range remarks {
		if r[1] == "SITE_IDENTIFIER" {
			identifier = strings.TrimSpace(r[2])
		}
		if r[1] == "SITE_DESCRIPTION" {
			siteDescs[identifier] = strings.TrimSpace(r[2])
		}
	}

	pdb.BindingSiteDesc = siteDescs

}
