import pytest
from unittest.mock import patch, mock_open, MagicMock
import os
import sys
import pickle
import pandas as pd
from pathlib import Path
from orthogonal_seq_screen_parallel import main, parse_arguments, ensure_correct_arg_input
from sequence_optimization_functions import load_ensemble, get_compseq, fill_dGs, ddG_optimization_algorithm
import numpy as np
from nupack import Model

@pytest.fixture
def mock_filesystem(tmp_path):
    # Create a temporary file for candidate sequences with 100 sequences
    example_sequences_file = "seq100.txt"
    seqs = load_ensemble(example_sequences_file)

    candidate_sequences_file = tmp_path / example_sequences_file
    candidate_sequences_file.write_text('\n'.join(seqs))
    
    # Create a dummy dG file with valid pickle data
    interim_array = np.zeros((len(seqs),len(seqs)),dtype=float)
    dGs_filled = pd.DataFrame(interim_array, index=seqs,columns=seqs)
    dG_file = tmp_path / "dGs_filled.pkl"
    dGs_filled.to_pickle(dG_file)
    
    # Path for the output sequences file
    output_sequences_file = tmp_path / "optimum_sequences.txt"

    return {
        "candidate_sequences_file": candidate_sequences_file,
        "dG_file": dG_file,
        "output_sequences_file": output_sequences_file
    }

def test_parse_arguments():
    testargs = [
        "program_name",
        "-f", "dummy_sequences.txt",
        "-n", "6",
        "-s", "100",
        "-g", "dummy_dGs.pkl",
        "-o", "dummy_output.txt"
    ]
    with patch.object(sys, 'argv', testargs):
        args, parser = parse_arguments()
        assert args.candidate_sequences_file == "dummy_sequences.txt"
        assert args.num_seq_output == 6
        assert args.optimization_steps == 100
        assert args.dG_file == "dummy_dGs.pkl"
        assert args.output_sequences_file == "dummy_output.txt"

def test_ensure_correct_arg_input(mock_filesystem):
    args = MagicMock()
    args.candidate_sequences_file = str(mock_filesystem["candidate_sequences_file"])
    args.num_seq_output = 6
    args.optimization_steps = 100
    args.dG_file = str(mock_filesystem["dG_file"])
    args.output_sequences_file = str(mock_filesystem["output_sequences_file"])
    args.force_overwrite = False
    parser = MagicMock()

    candidate_sequences_file, num_seq_output, optimization_steps, dG_file, output_sequences_file = ensure_correct_arg_input(args, parser)
    
    assert candidate_sequences_file == mock_filesystem["candidate_sequences_file"]
    assert num_seq_output == 6
    assert optimization_steps == 100
    assert dG_file == mock_filesystem["dG_file"]
    assert output_sequences_file == mock_filesystem["output_sequences_file"]


def test_load_ensemble(mock_filesystem):
    seqs = load_ensemble(mock_filesystem["candidate_sequences_file"])
    assert len(seqs) == 200


def test_get_compseq(mock_filesystem):
    seq = 'ATCCGT'
    compseq = get_compseq(seq)
    assert compseq == 'ACGGAT'
    

def test_fill_dGs_with_dG_file(mock_filesystem):
    seqs = load_ensemble(mock_filesystem["candidate_sequences_file"])
    dGs = fill_dGs(seqs, None, mock_filesystem["dG_file"])
    print(seqs)
    print(dGs.columns.tolist())
    assert dGs.shape == (200, 200)
    # How to asset all sequence is seqs is in dGs.columns and dGs.index but order dosent matter
    assert set(dGs.columns.tolist()) == set(seqs)
    assert set(dGs.index.tolist()) == set(seqs)
    assert dGs.values.min() == 0.0


def test_fill_dGs_no_dG_file(mock_filesystem):
    model = Model(material='dna', celsius=10)
    seqs = load_ensemble(mock_filesystem["candidate_sequences_file"])
    tmp_path = Path.cwd()
    dGs = fill_dGs(seqs, model, tmp_path / 'new_file.pkl')
    assert os.path.exists(tmp_path / 'new_file.pkl')
    assert dGs.shape == (200, 200)
    assert set(dGs.columns.tolist()) == set(seqs)
    assert set(dGs.index.tolist()) == set(seqs)    


def test_ddG_optimization_algorithm(mock_filesystem):
    model = Model(material='dna', celsius=10)
    seqs = load_ensemble(mock_filesystem["candidate_sequences_file"])
    tmp_path = Path.cwd()
    dGs = fill_dGs(seqs, model, tmp_path / 'new_file.pkl')
    
    best_subset, best_distance, subset_df = ddG_optimization_algorithm(dGs, 6, 100)

def test_main(mock_filesystem):
    testargs = [
        "program_name",
        "-f", str(mock_filesystem["candidate_sequences_file"]),
        "-n", "2",
        "-s", "10",
        "-g", str(mock_filesystem["dG_file"]),
        "-o", str(mock_filesystem["output_sequences_file"]),
        "--force_overwrite"
    ]
    with patch.object(sys, 'argv', testargs):
        result = main()
        
        assert result == 0
        # Check if the output file was created
        assert mock_filesystem["output_sequences_file"].exists()
        # Check if the output file contains the expected number of sequences
        with open(mock_filesystem["output_sequences_file"], 'r') as f:
            lines = f.readlines()
            assert len(lines) == 2  # should match num_seq_output

if __name__ == "__main__":
    pytest.main()
