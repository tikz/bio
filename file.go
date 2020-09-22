package main

import (
	"encoding/gob"
	"fmt"
	"os"

	"github.com/tikz/bio/pdb"
	"github.com/tikz/bio/uniprot"
)

const (
	dataDir = "data/"
	unpDir  = dataDir + "uniprot/"
	pdbDir  = dataDir + "pdb/"
	fileExt = ".data"
)

func makeDirs() {
	os.MkdirAll(unpDir, os.ModePerm)
	os.MkdirAll(pdbDir, os.ModePerm)
}

func LoadPDB(pdbID string) (*pdb.PDB, error) {
	makeDirs()

	path := pdbDir + pdbID + fileExt
	_, err := os.Stat(path)
	if os.IsNotExist(err) {
		p, err := pdb.NewPDBFromID(pdbID)
		if err != nil {
			return nil, err
		}

		err = write(path, &p)
		if err != nil {
			return nil, fmt.Errorf("write PDB: %v", err)
		}
		return &p, nil
	}

	return readPDB(pdbID)
}

func readPDB(pdbID string) (*pdb.PDB, error) {
	path := pdbDir + pdbID + fileExt
	p := new(pdb.PDB)
	err := read(path, &p)
	if err != nil {
		return nil, fmt.Errorf("load file: %v", err)
	}

	err = p.Parse()
	return p, err
}

func LoadUniProt(unpID string) (*uniprot.UniProt, error) {
	makeDirs()

	path := unpDir + unpID + fileExt
	_, err := os.Stat(path)
	if os.IsNotExist(err) {
		u, err := uniprot.NewUniProt(unpID)
		if err != nil {
			return nil, err
		}

		err = write(path, &u)
		if err != nil {
			return nil, fmt.Errorf("write UniProt: %v", err)
		}
		return u, nil
	}

	u := new(uniprot.UniProt)
	err = read(path, &u)
	if err != nil {
		return nil, fmt.Errorf("load file: %v", err)
	}

	return u, nil
}

func write(filePath string, object interface{}) error {
	file, err := os.Create(filePath)
	if err == nil {
		encoder := gob.NewEncoder(file)
		err = encoder.Encode(object)
		if err != nil {
			return err
		}
	}
	file.Close()

	return err
}

func read(filePath string, object interface{}) error {
	file, err := os.Open(filePath)
	if err == nil {
		decoder := gob.NewDecoder(file)
		err = decoder.Decode(object)
	}
	file.Close()

	return err
}
