package pdb

import (
	"math"
	"varq/pdb"
)

// ResiduesDistance returns the distance between two residues as the minimum atom pair distance.
func ResiduesDistance(r1 *pdb.Residue, r2 *pdb.Residue) float64 {
	minDist := AtomsDistance(r1.Atoms[0], r2.Atoms[1])
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
func AtomsDistance(a1 *pdb.Atom, a2 *pdb.Atom) float64 {
	return math.Sqrt(math.Pow(a1.X-a2.X, 2) + math.Pow(a1.Y-a2.Y, 2) + math.Pow(a1.Z-a2.Z, 2))
}
