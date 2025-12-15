# Sticky End Designer (SEDesign)

Sticky End Designer Algorithm by Hao Liu and Matthew Sample, with some comments and type hints added by Josh Evans.  
Scripts for the design of orthogonal DNA sticky end sequences, with free energy calculations based on **NUPACK 4.0**.

## What this folder contains

This folder is intended to be used as a small “toolbox” for generating sets of **mutually-orthogonal sticky ends** (minimizing undesired cross-hybridization)

## Requirements

- **Python 3** with the dependencies listed in [`SEdesign/setup.py`](SEdesign/setup.py).
- **NUPACK 4.0** installed and accessible from the environment where you run the scripts 

## Installation (recommended)

From the repo root:

```
pip install -e ./SEdesign
```

Please follow the NUPACK documentation for NUPACK installation:
<https://docs.nupack.org/start/#macoslinux-installation>

## General workflow:

There are three possible starting states compatible with our tools:

1. Starting from scratch, you have no candidate seqeucens you want to find an optimal subset of
2. You have a list of candidate sequences you want to find an optimal subset of
3. You have a list of candidate sequences and set-in-stone sequences you need to be a part of the solution (constraint seqs) that you want to find the optimal union subset

### Starting from scratch
1. Create a list of all possible sticky ends with same melting temperature using the `create_same_tm_seq_library.py` script

    a. Decide target overhang melting temperature and a tolerance (EX: 15C +- 0.5C)

    b. Choose a length for all candidate overhangs

    c. Choose a gc content 

### You have a list of candidate sequences
2. Find most orthogonal subset of candidate using `orthogonal_seq_screen_parallel.py`

    a. Determine number of output pairs (EX: 50)

    b. Determine if you want to add polyT spacer to overhangs (EX: TTTTT)

    c. Determine how much secondary strucutre is allowed for each strand (EX: 0)

This script will give you your desired number of orthogonal overhang pairs

### You have a list of candidate sequences and set-in-stone sequences
3. Sample as above, but add the contraint seqs to input to `orthogonal_seq_screen_parallel.py` and use `sequence_constraint_file` flag

## Example usage

Example usage of `create_same_tm_seq_library.py` and `orthogonal_seq_screen_parallel.py` are shown in the `./example` folder.

`create_same_tm_seq_library.py` requires no input files 

This script will:

1. create all possible sequences of a given size and gc content

2. filter out all sequences not with in the melting temperature tolerance

`orthogonal_seq_screen_parallel.py` requires a list of candidate sequences, with optionally a list of contraint sequences

This script will:

1. compute a dG matrix between all sequences

2. run a mcmc algorithm to select the best subset using ddG matrix to help determine next step



# WARNINGS
`orthogonal_seq_screen_parallel.py` has a `--dG_file` flag that was made for one purpose, to allow you to change the number of output pairs and the number of optimization steps without rerunning the computation of the dG matrix. IF the `--dG_file` flag is not changed and ANY changes are made EXCEPT for changing number of output pairs or number of optimization steps, the script will FAIL.
