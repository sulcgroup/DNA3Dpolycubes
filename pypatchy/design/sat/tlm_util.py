import logging
import re
import tempfile
from math import sqrt
from pathlib import Path
from typing import Union, Any, Callable

import libtlm
import numpy as np
import scipy
from scipy.stats import lognorm, truncnorm

from .solution import PolysatSolution
from ..solve_utils import toPolycube, CrystalloydRating, gen_map_crystalloyd_onto, setup_logger
from ...polycubeutil import tlm_data
from ...polycubeutil.polycube_structure import PolycubeStructure

from pypatchy.polycubeutil.tlm_record_set import TLMRecords
from ...polycubeutil.polycubesRule import RULE_ORDER, PolycubesRule
from ...polycubeutil.tlm_data import TLMHistoryRecord
from ...structure import TypedStructure
from ...util import get_log_dir


def records_polycubes(records: list[TLMHistoryRecord], rule: PolycubesRule) -> list[list[PolycubeStructure]]:
    return [[toPolycube(rule, pc)
             for pc in record.getPolycubes()] for record in records]


def run_tlm(i: int, tlm_params: libtlm.TLMParameters, n_sims: int, save_to_dir: Path, record_interval: Union[int, None]=None) -> list[TLMRecords]:
    tlm_params.write_to_file(str(get_log_dir() / "SAT" / f"last_settings_{i}.json"))
    results: list = [None for _ in range(n_sims)]
    logger = setup_logger(f"run_tlm_{i}")
    logging.getLogger()
    if not record_interval:
        record_interval = int(max(1e6, tlm_params.num_steps // int(1e6)))
    with tempfile.TemporaryDirectory() as td:
        # run sims and save data to temporary directory
        libtlm.runSimsAndSave(tlm_params,
                              n_sims,
                              record_interval,  # avoid spamming stout
                              td)
        logger.info(f"Execution complete!")
        for fp in Path(td).glob("sim*.json"):
            match = re.search(r"sim(\d+)\.json", fp.name)
            assert match
            isim: int = int(match.group(1))  # Extracted integer
            records = TLMRecords(fp)
            logger.info(f"Cached data from simulation {isim} in file {str(save_to_dir / f'sim_{isim}.pickle')}")
            records.export(save_to_dir / f"sim_{isim}.pickle")
            last_record = records.records[max(records.records.keys())]  # off by one: TODO FIX
            results[isim] = last_record
            # do i need to explicitly dealloc memory?
    assert not any([result is None for result in results])
    return results


def spinoff_logger(parent_logger: logging.Logger, name: str) -> logging.Logger:
    """
    logging helper function written by chatGPT
    todo: align with logging semantics in solve_exec.solve_worker
    """
    child_logger = logging.getLogger(name)
    child_logger.setLevel(parent_logger.level)

    # if not child_logger.handlers:
    #     for h in parent_logger.handlers:
    #         child_logger.addHandler(h)

    return child_logger


def get_max_nL(topologies: list[TypedStructure]) -> int:
    return max([len(top) for top in topologies])


def score_sim_group(analyze_results: list[Union[list[CrystalloydRating], None]],
                    score_func: Callable,
                    n_sim_cubes: int,
                    logger: logging.Logger,
                    min_score: float = 0.5) -> float:
    """
    scoring function for scikit-optimize bayesian optimization
    bayesian optimization tries to MINIMIZE the score, so for this function (unlike Polysat.analyze_crystalloyd_scores)
    a higher score is WORSE

    """
    n_melts = len([True for r in analyze_results if not r])  # probably a nicer way to do this
    # todo: make threshold settable
    if n_melts > len(analyze_results) // 2:
        return n_melts * 1e5
    # consider each polycube sepeerrately
    top0_scores = []
    top1_scores = []
    # flatten list of analyzed crystalloyds

    for sim, ratings in enumerate(analyze_results):
        if ratings:  # if scoring is not None
            # see function C(S, N) in doc "Multifarious crystals et al.", section "Scoring a Simulation"
            top0_scores.append(sum([
                score_func(rating) * (rating.pc.num_particles() ** 2) if score_func(rating) >= min_score else 0
                for rating in ratings if rating.topology_idx == 0
            ]) / (n_sim_cubes ** 2))
            top1_scores.append(sum([
                score_func(rating) * (rating.pc.num_particles() ** 2) if score_func(rating) >= min_score else 0
                for rating in ratings if rating.topology_idx == 1
            ]) / (n_sim_cubes ** 2))
            logger.info(f"Simulation {sim}: C_1(S) = {top0_scores[-1]}, C_2(S) = {top1_scores[-1]}")
            assert top0_scores[-1] + top1_scores[-1] <= 1.

    # now compute simulation group score (see "Scoring a Simulation Group" in "Multifarious Crystals et al.")
    score = (2 * sqrt(sum(top0_scores) * sum(top1_scores))) / len(analyze_results)
    logger.info(f"Score")
    return 1 / score if score > 0 else 1e8  # return a large value


def analyze_crystalloyds_from_record(isim: int,
                                     record: tlm_data.TLMHistoryRecord,
                                     sat_solution: PolysatSolution,
                                     logger: logging.Logger,  # todo: get rid of this so we can parallelize
                                     hyperparams: dict,
                                     topologies: list[TypedStructure],
                                     nS: int) -> Union[bool, list[CrystalloydRating]]:
    """
    compares all crystalloyd candidtes in the provided record with the provided topologies
    """
    compare_crystalloyds(
        isim=isim,
        record=record,
        rule=sat_solution.rule,
        typing_map=sat_solution.nanoparticle_map,
        hyperparams=hyperparams,
        logger=logger,
        topologies=topologies,
        nS=nS
    )

def compare_crystalloyds(isim: int,
                         record: tlm_data.TLMHistoryRecord,
                         rule: PolycubesRule,
                         typing_map: dict,
                         logger: logging.Logger,  # todo: get rid of this so we can parallelize
                         hyperparams: dict,
                         topologies: list[TypedStructure],
                         nS: int) -> Union[bool, list[CrystalloydRating]]:
    """
    tests polycubes from a simul
    """

    nL = get_max_nL(topologies)
    logger.info(f"Testing polycubes from simulation #{isim}....")
    # convert polycube data to objects
    assert record is not None and record.energy() is not None
    polycube_data: list[PolycubeStructure] = [toPolycube(rule, pc)
                                              for pc in record.getPolycubes()]
    logger.info(f"Sim #{isim} has {len(polycube_data)} polycubes")
    # get the biggest polycube
    n_unit_cell_cutoff = hyperparams["n_unit_cell_cutoff"]
    # should already be sorted, but might as well
    polycube_data = sorted(polycube_data, key=len, reverse=True)

    ratings: list[CrystalloydRating] = list()

    # for each polycube in the record
    for ipc, pc in enumerate(polycube_data):
        logger.info(f"Testing polycube {ipc} of size {len(pc)} from simulation {isim} at timestep {record.stepCount()}")
        # test if polycube satisfies crystal bindings
        # we are not using my SAT-test :((((
        # if not self.polycube_satisfies(solution_vars, pc):
        # gen blank map
        npmap = {
            i: 0 for i in range(nS)
        }
        # update w/ rule nps from soln
        npmap.update(typing_map)
        rating = map_crystalloyd_to_topology(pc, npmap, topologies)
        if rating.topology_idx == -1:
            logger.info(f"Crystal num. {ipc} for simulation no. {isim} "
                        f"does not match expected crystal topology!")
        else:
            # if polycube is too small
            if len(pc) < n_unit_cell_cutoff * nL:
                logger.info(f"Polycube size {len(pc)} is less than {n_unit_cell_cutoff} unit cells "
                            f"(unit cell size = {nL}), not big enough to test for"
                            f" crystalization...")
                # if this polycube is the largest in the simulation, yell bc we didn't find a polycube big enough
                if ipc == 0:
                    logger.warning(f"Simulation did not produce any polycube large enough to test for "
                                   f"crystalloydization! Largest polycube size was {len(pc)}, "
                                   f"cutoff was {n_unit_cell_cutoff} x {nL} ="
                                   f" {n_unit_cell_cutoff * nL}. Conclusion: rule does "
                                   f"not crystallize properly.")
                    return [CrystalloydRating(pc)]  # return a crystalloyd with no data, only polycube
                # if we did find good polycubes, note how many n continue
                else:
                    logger.info(f"Simulation {isim} had {ipc + 1}, analysis of polycubes complete!")
                    break
                # if the polycube has formed a crystal other than the desired one, return false
                # use SAT-solver approach
            else:  # polycube is large enough
                # we don't actually want to reject rules here, it may be valid if we fuck with the
                ratings.append(rating)
    return ratings


def map_crystalloyd_to_topology(pc: PolycubeStructure, np_map: dict[int, int],
                                topologies: list[TypedStructure]) -> CrystalloydRating:
    mapping = None
    for (top_idx, typed_structure) in enumerate(topologies):

        # TODO this is going to explode pretty badly when I introduce multifarious crystals
        # mapping can be an aggregate, what we are weeding out here are the off target structures
        mapping = gen_map_crystalloyd_onto(pc,
                                           typed_structure,
                                           np_map)
        if mapping:
            mapping.topology_idx = top_idx
            return mapping
    # if mapping has not been found
    # we don't actually want to reject rules here, it may be valid if we fuck with the
    # interaction strengths
    return CrystalloydRating(pc)

def construct_trunc_norm(n1: float, n2:float):
    if n1 > n2:
        n1, n2 = n2,n1
    # Mean and standard deviation
    mean = (n1 + n2) / 2
    std = (n2 - n1) / 6  # Keeps the shape nice and centered

    # Calculate truncation limits in standard normal units
    a, b = (n1 - mean) / std, (n2 - mean) / std

    # Create truncated normal distribution
    return truncnorm(a, b, loc=mean, scale=std)

def fit_lognorm_from_bounds(lower: float, upper: float, frac_within: float = 0.9):
    """
    via chatGPT. SHOULD construct a lognorm distribution where frac_within of values will be between lower and upper
    """
    assert upper > lower > 0.
    # Convert percentage to quantiles
    q1 = (1 - frac_within) / 2
    q2 = 1 - q1

    # Use log-normal quantiles: we solve for mu, sigma in log-space
    log_lower = np.log(lower)
    assert not np.isnan(log_lower)
    log_upper = np.log(upper)
    assert not np.isnan(log_upper)

    # Estimate mu and sigma assuming symmetrical quantiles
    mu = (log_lower + log_upper) / 2
    sigma = (log_upper - log_lower) / (2 * scipy.stats.norm.ppf(q2))

    # Create the distribution
    dist = lognorm(s=sigma, scale=np.exp(mu))
    return dist


def score_crystalloyd_out_node_frac(result: CrystalloydRating) -> float:
    """
    Compute the fraction of the polycube structure cubes that have fewer than 6 connections
    """
    n_border_nodes = 0
    for node in result.pc.graph.nodes:
        # TODO: better handling for sparse crystal unit cells
        is_border_node = len(result.pc.graph.edges(node)) < (len(result.unit_cell.graph.edges(result.mapping[node]))
                                                             if result.unit_cell is not None else len(RULE_ORDER))
        if is_border_node:
            n_border_nodes += 1

    n_particles = result.pc.num_particles()
    assert n_particles >= n_border_nodes
    # Raw interior count
    n_interior = n_particles - n_border_nodes

    # Raw interior fraction
    raw_frac = n_interior / n_particles

    # Normalize by expected interior fraction of a compact cube-like shape
    # Surface area ~ N^(2/3), volume = N => border ~ N^(2/3)
    expected_border = n_particles ** (2 / 3)
    expected_interior = n_particles - expected_border
    expected_frac = expected_interior / n_particles

    # Normalize: 1 = as compact as ideal cube, 0 = all surface
    normalized_score = raw_frac / expected_frac if expected_frac > 0 else 0.0
    normalized_score = min(max(normalized_score, 0.0), 1.0)  # Clamp to [0, 1]

    return normalized_score


def score_crystalloyd_surf_cubes_frac(result: CrystalloydRating) -> float:
    """
    similar to the above out_node_frac method
    Compute the fraction of the polycube structure cubes that have fewer than 6 adjacent cubes
    """
    n_border_cubes = 0
    for cube in result.pc.particles():
        for d in RULE_ORDER:
            if not result.pc.has_cube_at(cube.position() + d):
                n_border_cubes += 1
                break

    n_particles = result.pc.num_particles()
    assert n_particles >= n_border_cubes
    # Raw interior count
    n_interior = n_particles - n_border_cubes

    # Raw interior fraction
    raw_frac = n_interior / n_particles

    # Normalize by expected interior fraction of a compact cube-like shape
    # Surface area ~ N^(2/3), volume = N => border ~ N^(2/3)
    expected_border = n_particles ** (2 / 3)
    expected_interior = n_particles - expected_border
    expected_frac = expected_interior / n_particles

    # Normalize: 1 = as compact as ideal cube, 0 = all surface
    normalized_score = raw_frac / expected_frac if expected_frac > 0 else 0.0
    normalized_score = min(max(normalized_score, 0.0), 1.0)  # Clamp to [0, 1]

    return normalized_score


def score_crystalloyd_stupid(result: CrystalloydRating) -> float:
    """
    uses a "stupid" algorithm to score polycube
    just measures crystalloyd compactness
    """
    mins, maxs = result.pc.minmaxs()
    bb_vol = np.product(maxs - mins)  # bounding box volume
    return result.pc.num_particles() / bb_vol  # fraction of bounding box which is occupied by crystalloyd
