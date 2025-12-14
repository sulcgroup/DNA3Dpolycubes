import copy
import itertools
import json
import logging
import math
import os.path
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union, Any

from libtlm import TLMParameters

from ..structure import Structure, GenericTypedStructure, read_topology
from ..util import get_input_dir


@dataclass
class CrystalloydizationTestParams:
    """
    dataclass that holds wide-scale info about how to test a crystal, for the sat solver
    """
    n_replicas: int = field(default=24)
    n_unit_cells: int = field(default=400)  # can apply even to finite size assemblies
    # simulation step count - I strongly advise against leaving this as default!
    n_steps: int = field(default=int(1e6))
    # historical record interval
    record_interval: int = field(default=1e3)
    density: float = field(default=0.1)

    # min, max, interval temperatures
    Tmin: float = field(default=0)
    Tmax: float = field(default=1)

    torsion: bool = field(default=True)

    def get_total_type_counts(self, cube_type_counts: list[int]):
        return [n * self.n_unit_cells for n in cube_type_counts]

    def __post_init__(self):
        self.n_steps = int(self.n_steps)


class MultifariousTypeRequirement:
    NO_MULTIFARIOUS = -1  # topology is not multifarious
    NO_REQUIREMENTS = 0
    ANY_SHARED = 1
    ALL_SHARED = 2


class SolveParams:
    """
    A class that holds the parameters for a SAT solver execution
    It include number of species, number of colors, etc.
    """
    # structure name
    name: str
    # structure topology

    # todo: flag non-continuous connections seperately
    topologies: Union[GenericTypedStructure, list[GenericTypedStructure]]
    # number of species to solve for
    nS: int
    # number of colors to solve for
    nC: int
    # forbid self interact?
    forbid_self_interact: bool
    # positions of nanoparticles
    nanoparticles: Union[dict[int, int], list[int]]
    # number of dimensions (currently locked to 3)
    nDim: int
    # if the structure we're trying to solve for is a crystal
    crystal: bool
    # maximum alternative attempts to solve problem
    maxAltTries: int
    # number of solutions to find (should always be 1?)
    nSolutions: int
    # number of seconds to run the SAT solver for before giving up
    solver_timeout: int
    # structure_ids: list[set[int]]

    tlm_params: Union[CrystalloydizationTestParams, None]

    # how should multifarious behavior be treated?
    multifarious_behavior: int

    # other stuff
    pc_analysis_hyperparams: Union[None, dict[str, Any]]

    # directories containing solutions to pre-forbid
    # note that these are NOT ALL solutions for THIS sat problem! they must be validated
    # by Polysat later
    # they aren't all even guaranteed to be real directories as of end of __init__!
    pre_forbid_solutions: list[Path]

    # for crystals only. if true, the SAT solver will use a two-step method to design the rule,
    # desiging the "bindings" and "extraConnections" seperately
    segregate_bindings: bool
    # todo: god any name except this one
    # if true, do not allow the same to appear in both internal and external bindings
    segregate_colors: bool

    # external topology files
    forbidden_topologies: list[GenericTypedStructure]

    def __init__(self, name, **kwargs):
        self.name = name
        # default behavior
        if "topologies" not in kwargs:
            self.topologies = [read_topology(kwargs)]
            self.nanoparticles = dict()
            self.structure_ids = [
                0 for _ in range(len(self.topologies))
            ]
            self.bindings = kwargs["bindings"]
            if "extraConnections" in kwargs:
                self.extraConnections = kwargs["extraConnections"]
            else:
                self.extraConnections = []
            self.nanoparticles = kwargs["nanoparticles"]

        else:
            assert len(kwargs["topologies"]) == 2
            # a list of two
            pt_idxs = set()
            self.bindings = []
            self.extraConnections = []
            self.nanoparticles = dict()
            self.structure_ids = []
            self.topologies = []
            for topology in kwargs["topologies"]:
                self.topology = read_topology(topology, pt_idxs)
                pt_idxs.update(self.topologies[-1].vertices())

        # fix nanoparticles if num np types == 1
        if len(self.nanoparticles) > 0 and len(set(self.nanoparticles.values())) == 1:
            self.nanoparticles = list(self.nanoparticles.keys())

        assert "topologies" not in kwargs or len(kwargs["topologies"]) == 1 or self.topology.is_multifarious()

        self.maxAltTries = 10 if "maxAltTries" not in kwargs else kwargs["maxAltTries"]
        if not isinstance(self.maxAltTries, int):  # typechecking bc this is easier than editing typescript code
            self.maxAltTries = int(self.maxAltTries)
        self.nS = kwargs["nSpecies"] if "nSpecies" in kwargs else None
        self.nC = kwargs["nColors"] if "nColors" in kwargs else None
        # Use 1 solution at a time(glucose) unless otherwise specified
        # TODO: option "dynamic" to choose nSolutions depending on topology size
        self.nSolutions = kwargs["nSolutions"] if "nSolutions" in kwargs else 1
        self.nDim = kwargs["nDim"] if "nDim" in kwargs else 3
        self.torsion = kwargs["torsion"] if "torsion" in kwargs else True
        self.crystal = bool(kwargs["crystal"]) if "crystal" in kwargs else len(self.extraConnections) > 0
        self.solver_timeout = kwargs["solve_timeout"] if "solve_timeout" in kwargs else 21600
        self.forbid_self_interact = bool(kwargs["forbid_self_interact"]) if "forbid_self_interact" in kwargs else False
        # default timeout: 6 hrs, which is probably both too long and not long enough

        # if we're testing a crystal, initialize parameters to pass to tlm
        if self.crystal:
            if "segregate_bindings" in kwargs:
                self.segregate_bindings = kwargs["segregate_bindings"]
            else:
                self.segregate_bindings = False # default to false
            if self.segregate_bindings and "segregate_colors" in kwargs:  # god any name but this
                self.segregate_colors = kwargs["segregate_colors"]
            else:
                self.segregate_colors = False
            if "tlm_params" in kwargs and kwargs["tlm_params"]:
                self.tlm_params = CrystalloydizationTestParams(torsion=self.torsion, **kwargs["tlm_params"])
            else:  # json file can have "none" or "false", both are false-y
                self.tlm_params = None
            if self.crystal and "crystalloyd_test_hyperparams" in kwargs and kwargs["crystalloyd_test_hyperparams"]:
                if "crystalloyd_test_hyperparameters" in kwargs:  # backwards compatibility
                    kwargs["crystalloyd_test_hyperparams"] = kwargs["crystalloyd_test_hyperparameters"]
                self.pc_analysis_hyperparams = {
                    # this is bs, don't do this
                    "n_unit_cell_cutoff": int(math.sqrt(self.tlm_params.n_unit_cells)),
                    "num_search_iters": 6,
                    "score_func": "out_node_frac",
                    "min_good_score": 0.6,
                    "bayes_chunk_size": 8,
                    "method": "isopotent"
                }
                if "crystalloyd_test_hyperparams" in kwargs:
                    if kwargs["crystalloyd_test_hyperparams"]:
                        self.pc_analysis_hyperparams.update(kwargs["crystalloyd_test_hyperparameters"])
                    else:
                        self.pc_analysis_hyperparams = None
        else:
            self.tlm_params = None
            # use hpyerparams to store finite size params
            self.pc_analysis_hyperparams = {
                "pcs_n_batches": 20,
                "pcs_batch_size": 4
            }
            if "polycubes_batches" in kwargs:
                self.pc_analysis_hyperparams["pcs_n_batches"] = kwargs["polycubes_batches"]
            if "batch_size" in kwargs:
                self.tlm_params["pcs_batch_size"] = kwargs["batch_size"]

        if "mf_reqs" in kwargs:
            if kwargs["mf_reqs"] == "all_shared":
                self.multifarious_behavior = MultifariousTypeRequirement.ALL_SHARED
            elif kwargs["mf_reqs"] == "any_shared":
                self.multifarious_behavior = MultifariousTypeRequirement.ANY_SHARED
            elif kwargs["mf_reqs"] == "none":
                self.multifarious_behavior = MultifariousTypeRequirement.NO_REQUIREMENTS
            else:
                raise ValueError(f"Invalid multifarious type requirement {kwargs['mf_reqs']}")
        else:
            if isinstance(self.topology, list):
                # default to ANY_SHARED since a structure with no multifarious sharing requirements isn't
                # really multifarious
                self.multifarious_behavior = MultifariousTypeRequirement.ANY_SHARED
            else:
                self.multifarious_behavior = MultifariousTypeRequirement.NO_MULTIFARIOUS

        # can pass a list of file paths to solutions to forbid from the outset
        if "pre_forbid" in kwargs:
            self.pre_forbid_solutions = [Path(soln_path) for soln_path in kwargs["pre_forbid"]]
        else:
            self.pre_forbid_solutions = []

        if "forbidden_topologies" in kwargs:
            self.forbidden_topologies = []
            for top in kwargs["forbidden_topologies"]:
                self.forbidden_topologies.append(read_topology(top))


    # topology: GenericTypedStructure = property(lambda self: GenericTypedStructure(bindings=self.bindings + self.extraConnections,
    #                                                                                 types={i: (1 if i in self.nanoparticles else 0) for i in self.get_locations()}))

    @property
    def topology(self):
        if len(self.topologies) > 1:
            return sum(self.topologies)
        else:
            return self.topologies[0]

    def get_locations(self) -> set[int]:
        # can ignore crystal bindings
        return set(itertools.chain.from_iterable([(x[0], x[2]) for x in self.bindings]))

    def num_locations(self) -> int:
        return len(self.get_locations())

    def has_nanoparticles(self):
        return len(self.nanoparticles) > 0

    def brief_descr(self):
        return f"{self.nC},{self.nS}"

    def get_logger(self):
        logger_name = self.get_logger_name()
        assert logger_name in logging.Logger.manager.loggerDict
        return logging.getLogger(logger_name)

    def get_logger_name(self):
        return f"{self.name}_{self.nS}S_{self.nC}C"

    def use_test_lattice(self) -> bool:
        """
        return true if the input params require solutions to be validated on a thernal lattice model,
        false otherwise
        """
        return self.tlm_params is not None and self.pc_analysis_hyperparams is not None


def construct(name: str, min_nC: int, max_nC: int, min_nS: int, max_nS: int, max_diff: Union[int, None],
              **kwargs) -> list[SolveParams]:
    """
    what
    """
    paramsets = [SolveParams(
        name,
        nSpecies=nS,
        nColors=nC,
        **kwargs)
        for nC, nS in itertools.product(range(min_nC, max_nC), range(min_nS, max_nS))
        if max_diff is None or abs(nS - nC) <= max_diff]
    paramsets.sort(key=lambda p: p.nS + p.nC)
    return paramsets

def load_solve_params(spec_file_path: Union[Path, str], nColors: int, nSpecies: int) -> SolveParams:
    if isinstance(spec_file_path, str):
        spec_file_path = Path(spec_file_path)
    if not spec_file_path.suffix == ".json":
        spec_file_path = spec_file_path.with_suffix(".json")
    spec_file_path = spec_file_path.expanduser()
    if not spec_file_path.is_file():
        spec_file_path = get_input_dir() / "topologies" / spec_file_path
    assert spec_file_path.is_file(), f"No such file as {str(spec_file_path)}"
    with spec_file_path.open('r') as f:
        spec = json.load(f)
    return SolveParams(name=spec_file_path.name[:spec_file_path.name.find(".")],
                       nSpecies=nSpecies, nColors=nColors, **spec)
