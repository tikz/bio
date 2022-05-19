package clinvar

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"github.com/tikz/bio/pdb"
)

const (
	summaryURL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
)

type ClinVar struct {
	summaryPath string
	SNPs        map[string][]Allele
}

// Allele represents a variant from ClinVar
type Allele struct {
	VariantID     string `json:"variantId"`
	Name          string `json:"name"`
	ClinSig       string `json:"clinSig"`
	ClinSigSimple int    `json:"clinSigSimple"`
	ProteinChange string `json:"proteinChange"`
	ReviewStatus  string `json:"reviewStatus"`
	Phenotypes    string `json:"phenotypes"`
	Chromosome    uint64 `json:"chromosome"`
	Start         uint64 `json:"start"`
	End           uint64 `json:"end"`
}

func NewClinVar(clinvarDir string) (*ClinVar, error) {
	clinvarDir = filepath.Clean(clinvarDir)
	f, err := os.Stat(clinvarDir)
	if os.IsNotExist(err) {
		return nil, err
	}
	if !f.IsDir() {
		return nil, fmt.Errorf("%s is not a dir", clinvarDir)
	}

	summaryPath := clinvarDir + "/variant_summary.txt"

	if err = downloadClinVar(summaryPath); err != nil {
		return nil, err
	}

	cv := &ClinVar{}
	cv.summaryPath = summaryPath
	cv.SNPs = make(map[string][]Allele)
	err = cv.load()

	return cv, err
}

func (cv *ClinVar) load() error {
	f, err := os.Open(cv.summaryPath)
	if err != nil {
		return err
	}

	s := bufio.NewScanner(f)
	alleles := 0
	r, _ := regexp.Compile(`\(p.([A-z]{3})([0-9]*)([A-z]{3})\)`)
	for s.Scan() {
		line := strings.Split(s.Text(), "\t")
		variantType := line[1]
		name := line[2]
		// synonymous := strings.Index(name, "=") != -1
		m := r.FindAllStringSubmatch(name, -1)
		coding := len(m) > 0
		assembly := line[16] == "GRCh38"
		// if variantType == "single nucleotide variant" && coding && !synonymous && assembly {
		if assembly {
			dbSNPID := line[9]

			var change string
			if variantType == "single nucleotide variant" && coding {
				_, _, fromAa := pdb.AminoacidNames(m[0][1])
				_, _, toAa := pdb.AminoacidNames(m[0][3])
				pos := m[0][2]
				change = fromAa + pos + toAa
			}

			clinSigSimple, _ := strconv.Atoi(line[7])
			chromosome, _ := strconv.ParseUint(line[18], 10, 64)
			start, _ := strconv.ParseUint(line[19], 10, 64)
			end, _ := strconv.ParseUint(line[20], 10, 64)
			allele := Allele{
				VariantID:     line[30],
				Name:          name,
				ClinSig:       line[6],
				ClinSigSimple: clinSigSimple,
				ProteinChange: change,
				ReviewStatus:  line[24],
				Phenotypes:    line[13],
				Chromosome:    chromosome,
				Start:         start,
				End:           end,
			}

			cv.SNPs["rs"+dbSNPID] = append(cv.SNPs["rs"+dbSNPID], allele)
			alleles++
		}
	}
	return nil
}

func (cv *ClinVar) GetVariant(dbSNPID string, proteinChange string) *Allele {
	if alleles, ok := cv.SNPs[dbSNPID]; ok {
		for _, allele := range alleles {
			if allele.ProteinChange == proteinChange {
				return &allele
			}
		}
	}
	return nil
}

// downloadClinVar downloads and decompresses the variant_summary.txt.gz file from the NCBI FTP.
func downloadClinVar(summaryPath string) error {
	_, err := os.Stat(summaryPath)
	if os.IsNotExist(err) {
		resp, err := http.Get(summaryURL)
		if err != nil {
			return err
		}
		defer resp.Body.Close()

		out, err := os.Create(summaryPath)
		if err != nil {
			return err
		}
		defer out.Close()

		gr, err := gzip.NewReader(resp.Body)
		defer gr.Close()

		_, err = io.Copy(out, gr)
		return err
	}

	return nil
}
