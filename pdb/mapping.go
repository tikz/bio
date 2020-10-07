package pdb

func (pdb *PDB) makeMappings() {
	pdb.ChainStartResNumber = make(map[string]int64)
	pdb.ChainEndResNumber = make(map[string]int64)
	for c := range pdb.Chains {
		pdb.ChainStartResNumber[c] = pdb.minChainPos(c)
		pdb.ChainEndResNumber[c] = pdb.maxChainPos(c)
	}

	pdb.calculateChainsOffset()

	// SEQRES chain and position to structure residues.
	pdb.SeqResChains = make(map[string]map[int64]*Residue)

	for chain, offset := range pdb.SeqResOffsets {
		pdb.SeqResChains[chain] = make(map[int64]*Residue)
		minPos := pdb.ChainStartResNumber[chain]
		for pos, res := range pdb.Chains[chain] {
			pdb.SeqResChains[chain][pos-minPos+offset+1] = res
		}
	}

	// UniProt canonical sequence position to structure residues.
	pdb.UniProtPositions = make(map[string]map[int64][]*Residue)
	// chainMappings := pdb.SIFTS.UniProt[pdb.UniProtID].Mappings
	for unpID, unp := range pdb.SIFTS.UniProt {
		pdb.UniProtPositions[unpID] = make(map[int64][]*Residue)
		for _, m := range unp.Mappings {
			var i int64
			for i = m.UnpStart; i <= m.UnpEnd; i++ {
				seqResPos := i - m.UnpStart + m.PDBStart.ResidueNumber
				if res, ok := pdb.SeqResChains[m.ChainID][seqResPos]; ok {
					pdb.UniProtPositions[unpID][i] = append(pdb.UniProtPositions[unpID][i], res)
					res.UnpPosition = i
					res.UnpID = unpID
				}
			}
		}
	}
}

// This alignment needs to be done since the residue numbers in ATOM tags doesn't always coincide with SEQRES positions.
// TODO: see if there is a value available somewhere to skip doing this.
func (pdb *PDB) calculateChainsOffset() {
	pdb.SeqResOffsets = make(map[string]int64)
	for chain := range pdb.Chains {
		var bestOffset, bestScore int

		minPos := pdb.ChainStartResNumber[chain]
		chainLength := pdb.ChainEndResNumber[chain] - minPos
		steps := len(pdb.SeqRes[chain]) - int(chainLength)

		for offset := 0; offset < steps; offset++ {
			score := 0
			for pos, res := range pdb.Chains[chain] {
				seqResPos := pos + int64(offset) - minPos
				if res.Name1 == pdb.SeqRes[chain][seqResPos].Name1 {
					score++
				}
			}
			if score > bestScore {
				bestScore = score
				bestOffset = offset
			}
		}

		pdb.SeqResOffsets[chain] = int64(bestOffset)
		for _, res := range pdb.Chains[chain] {
			res.Position = res.StructPosition - pdb.ChainStartResNumber[chain] + int64(bestOffset) + 1
		}
	}
}

// Helpers

func (pdb *PDB) chainKeys(chain string) (k []int64) {
	for pos := range pdb.Chains[chain] {
		k = append(k, pos)
	}

	return k
}

func (pdb *PDB) minChainPos(chain string) int64 {
	ck := pdb.chainKeys(chain)
	min := ck[0]
	for _, pos := range ck {
		if pos < min {
			min = pos
		}
	}

	return min
}

func (pdb *PDB) maxChainPos(chain string) int64 {
	ck := pdb.chainKeys(chain)
	max := ck[0]
	for _, pos := range ck {
		if pos > max {
			max = pos
		}
	}

	return max
}
