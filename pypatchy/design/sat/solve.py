# solves for a specific nC, nS pair
import argparse
import logging
import os
from pathlib import Path

from pypatchy.design.sat.design_rule import Polysat
from pypatchy.util import get_input_dir

logging.getLogger()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SAT solver')
    parser.add_argument('solve_spec', type=str, help='Solve specification file name')
    parser.add_argument("-c", "--num_colors", type=int, help="Number of colors used in the solution", required=True)
    parser.add_argument("-s", "--num_species", type=int, help="Number of species used in the solution", required=True)
    parser.add_argument("-b", "--bad_solutions", type=str, nargs="+", help="")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose mode")

    args = parser.parse_args()

    cpu_count_os = os.cpu_count()

    logger = logging.getLogger()

    solve_spec_path: Path = get_input_dir() / "topologies" / args.solve_space

    solve_spec = Polysat
