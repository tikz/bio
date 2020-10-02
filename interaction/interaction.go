package interaction

import (
	"github.com/tikz/bio/pdb"
)

func Chains(p *pdb.PDB, distance float64) map[*pdb.Residue][]*pdb.Residue {
	interacts := make(map[*pdb.Residue][]*pdb.Residue)
	var i1, i2 int
	for chainName1, chain1 := range p.Chains {
		for chainName2, chain2 := range p.Chains {
			if i2 > i1 && chainName1 != chainName2 {
				for _, res1 := range chain1 {
					for _, res2 := range chain2 {
						if pdb.ResiduesDistance(res1, res2) < distance {
							interacts[res1] = append(interacts[res1], res2)
							interacts[res2] = append(interacts[res2], res1)
						}
					}
				}
			}
			i2++
		}
		i1++
	}
	return interacts
}

func Hets(p *pdb.PDB, distance float64) map[*pdb.Residue][]string {
	interacts := make(map[*pdb.Residue][]string)
	for _, chain := range p.Chains {
		for _, res := range chain {
			hets := NearHets(p, res, distance)
			for het, near := range hets {
				if near {
					interacts[res] = append(interacts[res], het)
				}
			}
		}
	}
	return interacts
}

func NearHets(p *pdb.PDB, r *pdb.Residue, distance float64) map[string]bool {
	hets := make(map[string]bool)
	for _, atom := range r.Atoms {
		for _, hetAtom := range p.HetAtoms {
			hets[hetAtom.Residue] = false
			if pdb.AtomsDistance(atom, hetAtom) < distance {
				hets[hetAtom.Residue] = true
				break
			}
		}
	}
	return hets
}
