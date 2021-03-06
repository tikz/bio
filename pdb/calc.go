package pdb

import (
	"math"
)

// ResiduesDistance returns the distance between two residues as the minimum atom pair distance.
func ResiduesDistance(r1 *Residue, r2 *Residue) float64 {
	minDist := AtomsDistance(r1.Atoms[0], r2.Atoms[0])
	for _, a1 := range r1.Atoms {
		for _, a2 := range r2.Atoms {
			dist := AtomsDistance(a1, a2)
			if dist < minDist {
				minDist = dist
			}
		}
	}
	return minDist
}

// AtomsDistance returns the distance between a pair of atoms.
func AtomsDistance(a1 *Atom, a2 *Atom) float64 {
	return math.Sqrt(math.Pow(a1.X-a2.X, 2) + math.Pow(a1.Y-a2.Y, 2) + math.Pow(a1.Z-a2.Z, 2))
}

// CloseResidues returns all residues within a given distance from the specified residue.
func CloseResidues(p *PDB, r *Residue, distance float64) (residues []*Residue) {
	for _, chain := range p.Chains {
		for _, res := range chain {
			if res != r && ResiduesDistance(res, r) < distance {
				residues = append(residues, res)
			}
		}
	}
	return
}
