# Threading software using double-dynamic programming

## Authors
Estelle Mariaux : estelle.mariaux@hotmail.fr

Théo Ferreira : theo.ferreira.med@gmail.com

University of Paris M2BI 2020 - 2021

## Short description

In order to model the structural representation of a protein (tertiary structure) based on its amino acid sequence (primary structure), the present project aims at reproducing the THREADER (David Jones 1998) software using the threading by double-dynamic programming.

# Installation

## Clone the repository

```
git clone https://github.com/estelleeeee/Threading_DPD.git

```
## Requirements

1. A linux distribution
2. Install the few required **python packages**:

```
conda create --name python3 python==3.8
conda activate python3
conda env update --file Threading_DPD.yml

# This command will install the following modules:
# python == 3.8
# numpy == 1.19.1
# argparse == 1.4.0

```

## Run the program

```
python src/main.py -p data/1n0a.pdb -d data/dope.par.txt -f data/1n09.fasta 
```
This line should result in a low level matrix for the first AA of your sequence, at the position 1.
```
For the fixed residue is CYS in position 1, the optimized score is -3.01
[[ 0.    0.    0.    0.    0.    0.    0.    0.    0.    0.    0.    0.  ]
 [ 0.    0.     nan   nan   nan   nan   nan   nan   nan   nan   nan   nan]
 [ 0.     nan  0.08  0.08  0.08  0.08  0.08  0.08  0.08  0.08  0.08  0.08]
 [ 0.     nan  0.08  0.08 -0.11 -0.11 -0.11 -0.11 -0.22 -0.22 -0.52 -1.31]
 [ 0.     nan  0.08  0.08 -0.11 -0.16 -0.16 -0.16 -0.27 -0.27 -0.68 -2.05]
 [ 0.     nan  0.08 -0.02 -0.19 -0.19 -0.19 -0.22 -0.45 -0.45 -0.82 -2.38]
 [ 0.     nan  0.08 -0.02 -0.19 -0.19 -0.21 -0.22 -0.45 -0.46 -0.92 -2.38]
 [ 0.     nan  0.08 -0.02 -0.19 -0.19 -0.21 -0.22 -0.45 -0.46 -0.92 -2.44]
 [ 0.     nan  0.08 -0.08 -0.37 -0.37 -0.37 -0.37 -0.62 -0.62 -0.98 -2.44]
 [ 0.     nan  0.08 -0.08 -0.38 -0.38 -0.39 -0.39 -0.64 -0.73 -1.3  -2.56]
 [ 0.     nan  0.08 -0.92 -0.92 -0.92 -0.92 -0.92 -0.92 -1.64 -3.01 -3.01]]
 ```
 And then give you the optimized score for every AA, for position 1, 4 and 7 (here are only the first 3 AA of the sequence for space management)
 ```
 For each amino acids, for the position 1, 4 or 7, we calculate the minimized score
-----
For the fixed residue is CYS in position 1, the optimized score is -3.01
For the fixed residue is CYS in position 4, the optimized score is -7.35
For the fixed residue is CYS in position 7, the optimized score is -6.56
-----
For the fixed residue is THR in position 1, the optimized score is -2.41
For the fixed residue is THR in position 4, the optimized score is -5.02
For the fixed residue is THR in position 7, the optimized score is -4.53
-----
For the fixed residue is TRP in position 1, the optimized score is -2.65
For the fixed residue is TRP in position 4, the optimized score is -4.64
For the fixed residue is TRP in position 7, the optimized score is -4.03
```
