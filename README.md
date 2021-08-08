# bio
Biological data formats and tool output parsers for Golang.

## PDB `tikz/bio/pdb`
https://pkg.go.dev/github.com/tikz/bio/pdb

Fetches and parses PDB files, using the Model/Chain/Residue/Atom architecture.

## UniProt `tikz/bio/uniprot`
https://pkg.go.dev/github.com/tikz/bio/uniprot

Fetches and parses UniProt TXT enties.

## Conservation (Pfam) `tikz/bio/conservation`
https://pkg.go.dev/github.com/tikz/bio/conservation

Fetches a given HMM model from Pfam, aligns it against the given sequence to get the conservation bitscore of each residue.

## FoldX `tikz/bio/foldx`
https://pkg.go.dev/github.com/tikz/bio/foldx

A wrapper for [FoldX](http://foldxsuite.crg.eu/) `RepairPDB` and `BuildModel` commands.

## SASA `tikz/bio/sasa`
https://pkg.go.dev/github.com/tikz/bio/sasa

Parses [FreeSASA](https://freesasa.github.io/) command line output.

## Fpocket `tikz/bio/fpocket`
https://pkg.go.dev/github.com/tikz/bio/fpocket

A wrapper for running [Fpocket](http://fpocket.sourceforge.net/) and parsing residues marked as pockets.


## DSSP `tikz/bio/dssp`
https://pkg.go.dev/github.com/tikz/bio/dssp

A wrapper for running [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) and parsing the assigned secondary structure to residues.

