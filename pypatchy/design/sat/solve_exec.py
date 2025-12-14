#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 11:16:24 2019
@author: joakim
"""
import argparse
import datetime
import sys
import uuid
from multiprocessing import Pool

from pypatchy.structure import Binding

from .solution import PolysatSolution
from ..solve_utils import *
from .design_rule import Polysat
from ..solve_params import *


def smartEnumerate(xMax, yMax):
    return sorted(
        ((x, y) for x in range(1, xMax + 1) for y in range(1, yMax + 1)),
        key=lambda i: sum(i)
    )

results = {}
finalResult = None
ongoing = 0

def solve(solve_params: SolveParams) -> Union[None, PolysatSolution]:
    """
    Args:
        solve_params: parameters for this solve attempt

    Returns:
        a SATSolution object, or None if no solution is found
    """
    tstart = datetime.now()
    soln = None

    logger = solve_params.get_logger()

    logger.info(f"{solve_params.nC} colors and {solve_params.nS} cube types")
    # try:
    if solve_params.crystal:
        raise Exception("This isn't public yet")
    else:
        mysat = Polysat(solve_params)
    mysat.init()
    logger.info(
        f"Constructed sat problem with {mysat.problem.num_variables()} variables and {mysat.problem.num_clauses()} clauses.")
    good_soln = mysat.find_solution()

    if good_soln:
        logger.info(f"WORKING RULE: {good_soln.decRuleNew()}")
    else:
        logger.info("Could not find a working solution!")

    logger.info(f"Solved for {solve_params.nS} types and {solve_params.nC} colors in {datetime.now() - tstart}")
    return good_soln


def solve_worker(params: tuple[int, int, int, dict, str], attach_parent=True):
    """
    worker for multiprocessing structure solve

    """
    (idx, nS, nC, data, topname) = params

    solve_params = SolveParams(
        topname,
        nColors=nC,
        nSpecies=nS,
        **data)

    while True:  # gonna lose my mind if Python doesn't eventually add do-while loops
        logger_fp = get_log_dir() / "SAT" / f"{solve_params.get_logger_name()}_{uuid.uuid1()}.log"
        if not logger_fp.is_file():
            break
            
    logger = setup_logger(solve_params.get_logger_name(), logger_fp)
    if attach_parent:
        # add superlogger handler
        for main_logger_handler in logging.getLogger(f"{topname}_main").handlers:
            logger.addHandler(main_logger_handler)

    logging.getLogger(f"{topname}_main").info(f"Attaching process to solve "
                                              f"nS={nS}, nC={nC} "
                                              f"logged to file `{str(logger_fp)}`")

    logger.info(f"------------------COMPUTING: nS={nS}, nC={nC}--------------")

    solution = solve(solve_params)
    if isinstance(solution, PolysatSolution):
        solution.printToLog()
        if solution.has_coord_map():
            solution.exportScene(topname)
            # return when we find solution


def solve_multi(solve_params: dict,
                solve_name: str,
                checks: list[tuple[int, int]],
                num_cpus: int):
    """
    batches solve of structures

    """
    worker_params = [(idx, nS, nC, solve_params, solve_name)
                     for (idx, (nS, nC)) in enumerate(checks)]
    logging.getLogger(f"{solve_name}_main").info(
        f"Starting solve pool for {solve_name} with {num_cpus} processes, over {len(checks)} specs.")
    if num_cpus > 1:
        with Pool(processes=num_cpus) as pool:
            pool.map(solve_worker, worker_params)
    else:
        for p in worker_params:
            solve_worker(p)


def solve_from_spec(topname: str, nS: int, nC: int, verbose: bool):
    solve_spec_path: Path = get_input_dir() / "topologies" / (topname + ".json")
    print(f"Loading solver specifications + topology from {str(solve_spec_path)}")
    main_logger = setup_logger(f"{solve_spec_path.stem}_main")
    if verbose:
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        main_logger.addHandler(handler)
    try:
        with solve_spec_path.open('r') as f:
            data = json.load(f)
            main_logger.info("Topology file: ")
            main_logger.info(str(solve_spec_path))
        solve_worker((0, nS, nC, data, topname))
    except FileNotFoundError:
        main_logger.error(f"No file {str(solve_spec_path)}")


def get_max_colors(solveSpec: dict) -> int:
    """
    Given a solve spec, returns the number of colors in the fully-addressable
     rule for the provided topology
    Args:
        solveSpec: a solve spec dict

    Returns: int, number of colors required for fully-addressable rule

    """
    if "topologies" in solveSpec:
        return sum([
            get_max_colors(top)
            for top in solveSpec["topologies"]
        ])
    else:
        nBind = len(solveSpec['bindings'])
        if 'extraConnections' in solveSpec:
            return nBind + len(solveSpec['extraConnections'])
        else:
            return nBind



def get_max_species(solveSpec: Union[dict, list[Binding]]) -> int:
    """
    Given a solve spec, returns the number of cube types (species) in the fully
    addressable rule for the provided topology
    Args:
        solveSpec: a solve spec dict

    Returns: int, number of cube types reqd for fully addressable rule
    """
    if isinstance(solveSpec, list) or isinstance(solveSpec, set):
        return max(itertools.chain.from_iterable([[
            x for (i, x) in enumerate(b) if i % 2 == 0]
            for b in solveSpec])) + 1
    else:
        if "topologies" in solveSpec:
            return sum([
                get_max_species(top)
                for top in solveSpec["topologies"]
            ])
        else:
            return get_max_species(solveSpec["bindings"])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SAT solver')
    parser.add_argument('solve_spec', type=str, help='Solve specification file name')
    parser.add_argument('--nS', type=int, help='Minimum number of species')
    parser.add_argument('--nC', type=int, help='Minimum number of colors')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')  # Add this line

    args = parser.parse_args()

    if args.solve_spec.find(".") == -1:
        args.solve_spec = args.solve_spec + ".json"

    try:
        with open(get_input_dir() / "topologies" / args.solve_spec, 'r') as f:
            filestr = f.read()

            data = json.loads(filestr)
            solve_params = SolveParams(
                args.solve_spec[args.solve_spec.rfind(os.sep) + 1:args.solve_spec.find(".")],
                topology=data['bindings'],
                nColors=args.nC,
                nSpecies=args.nS,
                **data)
            setup_logger(solve_params.get_logger_name())
            if args.verbose:
                handler = logging.StreamHandler(sys.stdout)
                handler.setLevel(logging.DEBUG)
                formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
                handler.setFormatter(formatter)
                solve_params.get_logger().addHandler(handler)

            solve_params.get_logger().info("Topology file: ")
            solve_params.get_logger().info(filestr)

            solve_name = args.solve_spec[:args.solve_spec.find(".")]

            solution = solve(solve_params)
            if isinstance(solution, PolysatSolution):
                solution.printToLog(logger=solve_params.get_logger())
                solution.exportScene(solve_name)
            else:
                solve_params.get_logger().info("No solution found!")

            solve_params.get_logger().info("Done!")

    except FileNotFoundError:
        logging.error(f"No file `{get_input_dir() / 'topologies' / args.solve_spec}`")
