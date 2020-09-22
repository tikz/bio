package pdb

import (
	"io/ioutil"
	"testing"
)

func LoadTestFile(path string) ([]byte, error) {
	data, err := ioutil.ReadFile(path)
	if err != nil {
		return nil, err
	}

	return data, nil
}

func TestChains(t *testing.T) {
	raw, err := LoadTestFile("./testdata/1mso.pdb")
	if err != nil {
		t.Errorf("cannot open file: %s", err)
	}

	pdb, err := NewPDBFromRaw(raw)
	if err != nil {
		t.Error(err)
	}

	t.Logf("testing PDB chains")

	actual := pdb.TotalLength
	expected := int64(102)
	if actual != expected {
		t.Errorf("expected %d, got %d", expected, actual)
	}

	res := pdb.Chains["D"][11]
	expect := "Leucine"
	if res.Name != expect {
		t.Errorf("expected %s in D-11, got %s", expect, res.Name)
	}

	expect = "Asparagine"
	res = pdb.Chains["C"][21]
	if res.Name != expect {
		t.Errorf("expected %s in C-21, got %s", expect, res.Name)
	}

	expect = "Phenylalanine"
	res = pdb.Chains["B"][1]
	if res.Name != expect {
		t.Errorf("expected %s in B-1, got %s", expect, res.Name)
	}

}

func TestSeqRes(t *testing.T) {
	raw, err := LoadTestFile("./testdata/1mso.pdb")
	if err != nil {
		t.Errorf("cannot open file: %s", err)
	}

	pdb, err := NewPDBFromRaw(raw)
	if err != nil {
		t.Error(err)
	}

	err = pdb.ExtractSeqRes(raw)
	if err != nil {
		t.Error(err)
	}

	t.Logf("testing SEQ RES")

	res := pdb.SeqRes["A"][0]
	expected := "Glycine"
	if res.Name != expected {
		t.Errorf("expected %s in SEQRES A-1, got %s", expected, res.Name)
	}

	res = pdb.SeqRes["B"][29]
	expected = "Threonine"
	if res.Name != expected {
		t.Errorf("expected %s in SEQRES B-30, got %s", expected, res.Name)
	}

}

func TestCIF(t *testing.T) {
	raw, err := LoadTestFile("./testdata/1mso.pdb")
	if err != nil {
		t.Errorf("cannot open file: %s", err)
	}

	pdb, err := NewPDBFromRaw(raw)
	if err != nil {
		t.Error(err)
	}

	rawCIF, err := LoadTestFile("./testdata/1mso.cif")
	if err != nil {
		t.Errorf("cannot open file: %s", err)
	}

	t.Logf("Testing CIF parse")

	err = pdb.ExtractCIFData(rawCIF)
	if err != nil {
		t.Errorf("cannot extract CIF: %s", err)
	}

	res := pdb.Title
	expected := "T6 Human Insulin at 1.0 A Resolution"
	if res != expected {
		t.Errorf("expected %s, got %s", expected, res)
	}

	res = pdb.Method
	expected = "X-RAY DIFFRACTION"
	if res != expected {
		t.Errorf("expected %s, got %s", expected, res)
	}

	if pdb.Date.Day() != 19 || pdb.Date.Month() != 9 || pdb.Date.Year() != 2002 {
		t.Errorf("expected date to be 2002-09-19")
	}

	if pdb.Resolution != 1.0 {
		t.Errorf("expected %f, got %f", 1.0, pdb.Resolution)
	}
}

func TestMappings(t *testing.T) {
	raw, err := LoadTestFile("./testdata/1mso.pdb")
	if err != nil {
		t.Errorf("cannot open file: %s", err)
	}

	pdb, err := NewPDBFromRaw(raw)
	if err != nil {
		t.Error(err)
	}

	pdb.ID = "1mso"
	uniProtID := "P01308"

	err = pdb.getSIFTSMappings()
	if err != nil {
		t.Error(err)
	}

	pdb.makeMappings()

	testChain := "B"
	chainSIFT, _ := pdb.SIFTS.GetChainMapping(uniProtID, testChain)
	pdbStart := chainSIFT.PDBStart.ResidueNumber
	pdbEnd := chainSIFT.PDBEnd.ResidueNumber
	unpStart := chainSIFT.UnpStart
	unpEnd := chainSIFT.UnpEnd

	if unpStart != 25 || unpEnd != 54 || pdbStart != 1 || pdbEnd != 30 {
		t.Errorf("received unexpected mapping positions")
	}

	resInSlice := func(r *Residue, s []*Residue) bool {
		for _, rs := range s {
			if rs == r {
				return true
			}
		}
		return false
	}

	testChain = "B"
	chainSIFT, _ = pdb.SIFTS.GetChainMapping(uniProtID, testChain)
	unpStart = chainSIFT.UnpStart
	var i int64
	for i = 1; i <= 30; i++ {
		if pdb.Chains[testChain][i] != pdb.SeqResChains[testChain][i] {
			t.Errorf("chain %s: misalignment between chain and SEQRES at pos %d", testChain, i)
		}
		if !resInSlice(pdb.Chains[testChain][i], pdb.UniProtPositions[uniProtID][i+unpStart-1]) {
			t.Errorf("chain %s: misalignment between chain and UniProt seq at pos %d", testChain, i)
		}
		if !resInSlice(pdb.SeqResChains[testChain][i], pdb.UniProtPositions[uniProtID][i+unpStart-1]) {
			t.Errorf("chain %s: misalignment between SEQRES and UniProt seq at pos %d", testChain, i)
		}
	}

	testChain = "C"
	chainSIFT, _ = pdb.SIFTS.GetChainMapping(uniProtID, testChain)
	unpStart = chainSIFT.UnpStart
	for i = 1; i <= 21; i++ {
		if pdb.Chains[testChain][i] != pdb.SeqResChains[testChain][i] {
			t.Errorf("chain %s: misalignment between chain and SEQRES at pos %d", testChain, i)
		}
		if !resInSlice(pdb.Chains[testChain][i], pdb.UniProtPositions[uniProtID][i+unpStart-1]) {
			t.Errorf("chain %s: misalignment between chain and UniProt seq at pos %d", testChain, i)
		}
		if !resInSlice(pdb.SeqResChains[testChain][i], pdb.UniProtPositions[uniProtID][i+unpStart-1]) {
			t.Errorf("chain %s: misalignment between SEQRES and UniProt seq at pos %d", testChain, i)
		}
	}
}
