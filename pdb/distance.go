package pdb

import (
	"math"
)

// Distance returns the distance between a pair of atoms
func Distance(atom1 *Atom, atom2 *Atom) float64 {
	return math.Sqrt(math.Pow(atom1.X-atom2.X, 2) + math.Pow(atom1.Y-atom2.Y, 2) + math.Pow(atom1.Z-atom2.Z, 2))
}

// ResiduesDistance returns the distance between residues, of the closest pair of atoms.
func ResiduesDistance(res1 *Residue, res2 *Residue) float64 {
	minDist := Distance(res1.Atoms[0], res2.Atoms[0])
	for _, a1 := range res1.Atoms {
		for _, a2 := range res2.Atoms {
			dist := Distance(a1, a2)
			if dist < minDist {
				minDist = dist
			}
		}
	}

	return minDist
}
