package foldx

import (
	"errors"
	"fmt"
	"io/ioutil"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"github.com/tikz/bio/pdb"
)

type FoldX struct {
	binFile         string
	binDir          string
	repairDir       string
	absRepairDir    string
	mutationsDir    string
	absMutationsDir string
}

// NewFoldX instantiates the required paths for FoldX.
func NewFoldX(foldxBinPath string, repairDirPath string, mutationsDirPath string) (foldx *FoldX, err error) {
	foldx = &FoldX{}
	foldx.binDir, foldx.binFile = filepath.Split(foldxBinPath)
	foldx.repairDir = filepath.Clean(repairDirPath)
	if foldx.absRepairDir, err = filepath.Abs(repairDirPath); err != nil {
		return
	}

	foldx.mutationsDir = filepath.Clean(mutationsDirPath)
	if foldx.absMutationsDir, err = filepath.Abs(mutationsDirPath); err != nil {
		return
	}

	err = func() error {
		if _, err := os.Stat(foldxBinPath); os.IsNotExist(err) {
			return err
		}
		if _, err := os.Stat(foldx.absRepairDir); os.IsNotExist(err) {
			return err
		}
		if _, err := os.Stat(foldx.absMutationsDir); os.IsNotExist(err) {
			return err
		}
		return nil
	}()
	return
}

// Repair runs the RepairPDB FoldX command on a PDB file specified by pdb.LocalPath
// and stores the resulting file in foldxDir, only if the file doesn't exist already.
// Returns the path where the repaired PDB is located.
func (foldx *FoldX) Repair(p *pdb.PDB) (outFile string, err error) {
	outFile = foldx.repairDir + "/" + p.ID + "_Repair.pdb"

	if fileNotExist(outFile) {
		tmpPDB, err := p.CopyPDB(foldx.binDir)
		if err != nil {
			return outFile, err
		}
		defer os.Remove(tmpPDB)

		cmd := exec.Command("./"+foldx.binFile,
			"--command=RepairPDB",
			"--pdb="+p.ID+".pdb",
			"--output-dir="+foldx.absRepairDir)
		cmd.Dir = foldx.binDir

		out, err := cmd.CombinedOutput()
		if err != nil {
			return outFile, err
		}

		if !strings.Contains(string(out), "run OK") || fileNotExist(outFile) {
			fmt.Println(string(out))
			return outFile, errors.New("RepairPDB failed")
		}
	}

	return outFile, nil
}

// BuildModelUniProt receives a given mutation in UniProt position and returns
// the PDB position in FoldX format, for only the first chain, i.e.: KA42I;
func (foldx *FoldX) BuildModelUniProt(repairedPath string, p *pdb.PDB, unpID string, pos int64, aa string) (float64, error) {
	residues := p.UniProtPositions[unpID][int64(pos)]
	if len(residues) == 0 {
		return 0, errors.New("no coverage")
	}
	res := residues[0]
	formattedMutant := res.Name1 + res.Chain + strconv.FormatInt(res.StructPosition, 10) + aa

	return foldx.BuildModel(repairedPath, formattedMutant)

}

func (foldx *FoldX) BuildModel(repairedPath string, formattedMutant string) (float64, error) {
	_, repairedName := filepath.Split(repairedPath)
	name := strings.Split(repairedName, "_")[0]

	destDirPath := foldx.absMutationsDir + "/" + name + "/" + formattedMutant
	diffPath := destDirPath + "/Dif_" + name + "_Repair.fxout"
	mutatedPDBPath := destDirPath + "/" + name + "_Repair_1.pdb"

	ddG, err := extractddG(diffPath)
	if fileNotExist(diffPath) || fileNotExist(mutatedPDBPath) || err != nil {
		// Create FoldX job output dir
		os.MkdirAll(destDirPath, os.ModePerm)

		// Create file containing individual list of mutations
		mutantFile := "individual_list_" + name + formattedMutant
		mutantPath := foldx.binDir + mutantFile
		writeFile(mutantPath, formattedMutant+";")

		// Copy PDB (FoldX can only open PDBs in the same dir)
		tmpPDBPath := foldx.binDir + "/" + repairedName
		_, err := exec.Command("cp", repairedPath, tmpPDBPath).Output()
		if err != nil {
			return ddG, err
		}

		// Remove files on scope exit
		defer func() {
			os.RemoveAll(tmpPDBPath)
			os.RemoveAll(mutantPath)

			// Duplicate of the original repaired PDB copied to mutation folder by FoldX
			os.RemoveAll(destDirPath + "/WT_" + name + "_Repair_1.pdb")
		}()

		cmd := exec.Command("./foldx",
			"--command=BuildModel",
			"--pdb="+name+"_Repair.pdb",
			"--mutant-file="+mutantFile,
			"--output-dir="+destDirPath)
		cmd.Dir = foldx.binDir

		out, err := cmd.CombinedOutput()
		if err != nil {
			return ddG, err
		}

		if !strings.Contains(string(out), "run OK") {
			fmt.Println(string(out))
			return ddG, errors.New("BuildModel failed")
		}
	}

	return extractddG(diffPath)
}

func extractddG(path string) (ddG float64, err error) {
	data, err := ioutil.ReadFile(path)
	if err != nil {
		return ddG, err
	}

	r, _ := regexp.Compile("pdb\t(.*?)\t")
	m := r.FindAllStringSubmatch(string(data), -1)
	if len(m) == 0 {
		return ddG, errors.New("ddG not found")
	}
	ddG, err = strconv.ParseFloat(m[0][1], 64)
	return ddG, err
}

func fileNotExist(path string) bool {
	_, err := os.Stat(path)
	return os.IsNotExist(err)
}

func writeFile(path string, contents string) {
	err := ioutil.WriteFile(path, []byte(contents), 0644)
	if err != nil {
		panic(err)
	}
}
