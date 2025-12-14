from __future__ import annotations

import argparse
import copy
from datetime import datetime # this is a crime
import itertools
import logging
import multiprocessing
import sys
import threading
import time
import uuid
from multiprocessing import Pool
from multiprocessing.managers import DictProxy
from pathlib import Path

import json
from typing import Union, Callable, Any, Optional

from .design_rule import Polysat
from .solution import SATSolution
from ..solve_params import load_solve_params, SolveParams
from ...util import get_input_dir, get_log_dir, get_output_dir
from ..solve_utils import setup_logger
from .solve_exec import get_max_colors, get_max_species, solve_multi


def multisolve(solve_spec: str,
               minS: int = 1,
               minC: int = 1,
               maxS: int = 0,
               maxC: int = 0,
               diff_limit: int = -1,
               num_cpus: int = 1,
               verbose: bool = False):
    if solve_spec.find(".") == -1:
        solve_spec = solve_spec + ".json"

    solve_spec_path: Path = get_input_dir() / "topologies" / solve_spec
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
            filestr = f.read()
            main_logger.info("Topology file: ")
            main_logger.info(filestr)
            data = json.loads(filestr)
            if maxC == 0:
                maxC = get_max_colors(data)
            if maxS == 0:
                maxS = get_max_species(data)
        assert minS <= maxS, f"Minimum num species {minS} is greater than maximum number of species {maxS}!"
        assert minC <= maxC, f"Minimum num colors {minC} is greater than maximum number of colors {maxC}!"
        checks = list(itertools.product(range(minS, maxS + 1), range(minC, maxC + 1)))
        if diff_limit > -1:
            checks = [(s, c) for s, c in checks if abs(s - c) <= diff_limit]

        checks.sort(key=lambda x: x[0] + x[1])

        solve_multi(data,
                    solve_spec_path.stem,
                    checks,
                    num_cpus)
        main_logger.info("Done!")

    except FileNotFoundError as e:
        if e.filename == solve_spec_path:
            main_logger.error(f"No file {str(solve_spec_path)}")
        else:
            raise e

class MultiSolver:
    """
    class to manage parallel solves with multiple solvers running with different nCs and nSs
    """

    # mapping of nC, nS pairs to Polysat (if solve in progress) or a list of solutions (if solve complete)
    data: dict[tuple[int, int]: Union[Polysat, SATSolution, False, None]]
    __refresh_queue_flag: bool
    __check_queue: list[tuple[int, int]]
    __minC: int
    __maxC: int
    __minS: int
    __maxS: int

    # spec to solve for, with dummy values for nC and nS (will fill in during execution)
    spec: SolveParams

    # function to generate nS, nC pairs
    check_gen_func: Callable[[MultiSolver], list[tuple[int, int]]]

    minS = property(fget=lambda self: self.__minS)
    num_cpus: int

    @minS.setter
    def minS(self, s: int):
        self.__refresh_queue_flag = s != self.__minS
        self.__minS = s

    maxS = property(fget=lambda self: self.__maxS)

    @maxS.setter
    def maxS(self, s: int):
        self.__refresh_queue_flag = s != self.__maxS
        self.__maxS = s

    minC = property(fget=lambda self: self.__minC)

    @minC.setter
    def minC(self, c: int):
        self.__refresh_queue_flag = c != self.__minC
        self.__minC = c

    maxC = property(fget=lambda self: self.__maxC)

    @maxC.setter
    def maxC(self, c: int):
        self.__refresh_queue_flag = c != self.__maxC
        self.__maxC = c

    # todo: class member function to generate + queue nC+nS pairs

    def __init__(self, solve_spec: str, **kwargs):
        self.logger = setup_logger(f"{solve_spec}_main")

        self.__refresh_queue_flag = True
        self.__check_queue = []

        if "num_cpus" in kwargs:
            self.num_cpus = kwargs["num_cpus"]
        else:
            self.num_cpus = 1

        self.spec = load_solve_params(solve_spec,
                                      -1,
                                      -1)

        # allow but don't require a custom check function (e.g. knight solve)
        if "gen_func" not in kwargs:
            if "diff_limit" in kwargs:
                diff_limit = kwargs["diff_limit"]
            else:
                diff_limit = -1
            def generate_checks(s: MultiSolver):
                checks = list(itertools.product(range(s.__minS, s.__maxS + 1), range(s.__minC, s.__maxC + 1)))
                checks = [(c, s) for s, c in checks
                          # nightmare filter statement
                          if (abs(s - c) <= diff_limit if diff_limit > -1 else True) and (s, c) not in self.data
                          ]
                return sorted(checks, key=lambda x: x[0] + x[1])
            self.check_gen_func = generate_checks
        else:
            # todo check that this is callable
            self.check_gen_func = kwargs["gen_func"]

        self.data = dict()
        # TODO: find and load exiting data

        if "minS" in kwargs:
            self.__minS = kwargs["minS"]
        else:
            self.__minS = len(self.spec.topologies) * max(1, len(self.spec.nanoparticles))
        if "minC" in kwargs:
            self.__minC = kwargs["minC"]
        else:
            self.__minC = 1
        if "maxS" in kwargs:
            self.__maxS = kwargs["maxS"]
        else:
            self.__maxS = get_max_species(self.spec.topology.bindings_list)
        if "maxC" in kwargs:
            self.__maxC = kwargs["maxC"]
        else:
            # can use shorthand here
            self.__maxC = len(self.spec.topology.bindings_list)

        if "verbose" in kwargs and kwargs["verbose"]:
            handler = logging.StreamHandler(sys.stdout)
            handler.setLevel(logging.DEBUG)
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)

    def get_check_queue(self) -> list[tuple[int, int]]:
        """
        Returns: a queue (as ordered list) containing nC, nS pairs for which solve hasn't started yet
        """
        # use a multiprocessing.Lock() to prevent changing this while running solves?
        # queue will need to be refreshed if nC, nS have been changed since last call to get_check_queue
        if self.__refresh_queue_flag:
            self.__check_queue = self.check_gen_func(self)
            self.__refresh_queue_flag = False # unset refresh queue flag
        return self.__check_queue

    def n_checks(self) -> int:
        return len(self.get_check_queue())

    def run(self, callback: Callable[[int, int, Optional[Union[SATSolution, bool]]], Any] = None):
        """
        Runs the SAT Solver
        :param callback: a callback function that will execute once each result is found for each nC, nS pair
        """
        queue = self.get_check_queue()
        worker_params = [copy.deepcopy(self.spec) for _ in range(len(queue))]

        for i, (nC, nS) in enumerate(queue):
            worker_params[i].nC = nC
            worker_params[i].nS = nS

        logging.getLogger(f"{self.spec.name}_main").info(
            f"Starting solve pool for {self.spec.name} with {self.num_cpus} processes, over {len(queue)} specs.")

        if self.num_cpus > 1:
            manager = multiprocessing.Manager()  # multiprocessing manager
            self.status = manager.dict()
            for nC, nS in queue:
                self.status[(nC, nS)] = "queued"
            worker_args = [(p, self.status) for p in worker_params]
            stop_flag = threading.Event()
            monitor_thread = threading.Thread(
                target=monitor_status,
                args=(self.spec.name, self.status),
                kwargs={"interval": 30, "stop_flag": stop_flag},
                daemon=True
            )
            monitor_thread.start()
            with Pool(processes=self.num_cpus) as pool:
                # generator (will
                for (nC, nS), result in pool.imap_unordered(solve_single_problem, worker_args):
                    if callback:
                        callback(nC, nS, result)
                    self.data[(nC, nS)] = result
            stop_flag.set()
            monitor_thread.join()

        else:
            for p in worker_params:
                (nC, nS), result = solve_single_problem((p, None))
                self.data[(nC, nS)] = result
                if callback is not None:
                    callback(nC, nS, result)

def monitor_status(name: str, status_dict: DictProxy, interval=60, stop_flag=None):
    time.sleep(1) # wait for solve processes to start
    stop_process = False
    while not stop_process:
        counts = {"queued": 0, "unsat": 0, "running": 0, "timeout": 0, "solved": 0}
        for v in status_dict.values():
            if v in counts:
                counts[v] += 1
            else:
                counts["solved"] += 1
                 # assume that all values for status_dict.values() are either a solution identifier or a known code
        logging.getLogger(f"{name}_main").info(
            f"Multisolve Progress: {counts['solved']} solved, {counts['unsat']} unsat,"
            f" {counts['timeout']} timeod out, {counts['running']} running, total {counts['solved']}"
        )

        write_status_table(status_dict, get_output_dir() / "SAT" / f"{name}_{datetime.now().strftime('%Y-%m-%d')}.txt")

        # do this at the end so it prints one final time
        stop_process = stop_flag is not None or stop_flag.is_set()
        time.sleep(interval)

def write_status_table(status_dict: DictProxy, out_path: Path):
    table = dict(status_dict)

    nCs = sorted({key[0] for key in table})
    nSs = sorted({key[1] for key in table})

    # Build header
    header_row = ["nC \\ nS"] + [f"{nS:>9}" for nS in nSs]
    rows = [" | ".join(header_row)]
    rows.append("-" * len(rows[0]))

    # Build each row
    for nC in nCs:
        row = [f"{nC:>9}"]
        for nS in nSs:
            status = table.get((nC, nS), "")
            row.append(f"{status:^9}")
        rows.append(" | ".join(row))

    with open(out_path, "w") as f:
        f.write("\n".join(rows) + "\n")


def solve_single_problem(args: tuple[SolveParams, Union[DictProxy, None]]) -> tuple[tuple[int, int], Union[SATSolution, bool, None]]:
    p, status_dict = args
    key = (p.nC, p.nS)
    if status_dict is not None:
        status_dict[key] = "running"
    while True:
        logger_fp = get_log_dir() / "SAT" / f"{p.get_logger_name()}_{uuid.uuid1()}.log"
        if not logger_fp.is_file():
            break
    logger = setup_logger(p.get_logger_name(), logger_fp)

    # add parent logger
    for handler in logging.getLogger(f"{p.name}_main").handlers:
        logger.addHandler(handler)

    logger.info(f"------------------COMPUTING: nS={p.nS}, nC={p.nC}--------------")
    tstart = datetime.now()
    logger.info("Setting up solver...")
    if p.crystal:
        raise Exception("This is not public yet.")
    else:
        mysat = Polysat(p)
    mysat.init()
    logger.info(f"Constructed SAT problem with {mysat.problem.num_variables()} variables and {mysat.problem.num_clauses()} clauses.")
    result = mysat.find_solution()

    if status_dict is not None:
        status_dict[key] = "done"

    if result:
        assert isinstance(result, SATSolution)
        logger.info(f"WORKING RULE: {result.decRuleNew()}. Saved to {str(mysat.solution_info_dir(result))}")
        if status_dict is not None:
            status_dict[key] = mysat.solution_info_dir(result).stem
    elif result is None:
        if status_dict is not None:
            status_dict[key] = "timeout"
        logger.info("Could not find a working solution!")
    else:
        if status_dict is not None:
            status_dict[key] = "unsat"
        logger.info(f"SAT problem {mysat.problem.describe_brief()} is not satisfiable!")

    logger.info(f"Finished in {datetime.now() - tstart}")
    return (p.nC, p.nS), result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SAT solver')
    parser.add_argument('solve_spec', type=str, help='Solve specification file name')
    parser.add_argument('--minS', type=int, default=1, help='Minimum number of species')
    parser.add_argument('--minC', type=int, default=1, help='Minimum number of colors')
    parser.add_argument('--maxS', type=int, default=0, help='Maximum number of species')
    parser.add_argument('--maxC', type=int, default=0, help='Maximim number of colors')

    parser.add_argument('--num_cpus', type=int, default=1, help='Number of CPUs available')

    parser.add_argument('--diff_limit', type=int, default=-1,
                        help="Limit on difference between nC and nS. The solver will skip specs where"
                             " |nC-nS| > diff_limit")

    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')

    args = parser.parse_args()

    multisolve(args.solve_spec,
               args.minS,
               args.minC,
               args.maxS,
               args.maxC,
               args.diff_limit,
               args.num_cpus,
               args.verbose)