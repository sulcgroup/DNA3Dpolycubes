"""
Polycube SAT specification adapted from Lukas chrystal one.

"Atoms" are "movable and rotatable", have 6 slots
"Positions" are fixed in the crystal, have 6 slots and bind according to spec
The problem has 2 parts:
A. find color bindings and colorings of position slots where each patch neightbors according to crystal model have
    colors that bind
B. find colorings of atoms s.t. all crystal positions are identical to (some) species rotation. The species definitions
    must not allow for bad 2-binds

indexes:
- colors:   1...c...#c (variable number)
- atoms:    1...s...#s (variable number)
- slots:    0...p...5=#p-1 (bindings places on atoms - 0,1,2 on one side, 3,4,5 on the other)
- position: 1...l...16=#l (number of positions in the crystal)
- rotation: 1...r...6=#r possible rotations of an species
- condition: 1...d...#d conditions to avoid bad crystal
- qualification: 0..#q (0 for good crystal, one more for each bad one)

(boolean) variables:
- B(c1, c2): color c1 binds with c2 (n=#c*#c)
- F(l, p, c): patch p at position l has color c (n=#l*#p*#c)
- C(s, p, c): patch p on species s has color c (n=#s*#p*#c)
- P(l, s, r): position l is occupied by species s with rotation r (n=#l*#s*#r)

encoding functions:
- rotation(p, r) = patch that p rotates to under rotation r

See: https://arxiv.org/pdf/2207.06954v1.pdf
"""
import datetime
import functools
import json
import pickle
import random
import uuid
from enum import Enum
from threading import Timer

from libtlm import TLMParameters
from pysat.formula import CNF
from pysat.solvers import Glucose4

from ..solve_params import *
from ..solve_utils import match_crystalloyd_to_unit_cell
from .tlm_util import *
from ...polycubeutil.polycubesRule import PolycubesRule, PolycubesPatch, get_orientation
from ...polycubeutil.polycubesRule import RULE_ORDER
from .polycube_sat_problem import PolycubeSATProblem
from .sat_problem import SATClause, interrupt
from ...polycubeutil.tlm_data import TLMHistoryRecord
from ...structure import TypedStructure, GenericTypedStructure
from ...util import tlm_data_dir

RELSAT_EXE = 'relsat'


# exit/error/response codes

# solut
class SolverResponse(Enum):
    SOLN_UNSAT = 0
    SOLN_UNDEFINED = 1
    SOLN_TIMEOUT = 2
    SOLN_ERR = 3


class Polysat:
    """
    class that facilitates a solve attempt for a specific structure and specific nC/nS
    """
    # ----------------------------------- basic stuff ---------------------------------------------
    # bindings
    # bindings: dict[tuple[int, int]: tuple[int, int]]  # ngl I have no idea why they're like this
    # internal_bindings: list[tuple[int, int, int, int]]

    # map of locations of nanoparticles (a list of locations if we only have one np type, a dict otherwise
    # todo:  use this + target_structure and change target_structure to TypedStructure
    target_structure: GenericTypedStructure

    # logger for SAT solver
    logger: logging.Logger = property(fget=lambda self: self.input_params.get_logger() )
    # timeout for SAT solver
    solver_timeout: Union[int, None]

    # whether patches are torsional
    torsionalPatches: bool

    # allo_clauses: Union[None, list[SATClause]]
    input_params: SolveParams

    problem: PolycubeSATProblem
    bindings = property(lambda self: self.target_structure.bindings_list)

    data_dir: Path
    __soln_info_dirs: dict[str, Path]
    # idxs should correspond to problem parts called forbidden{idx} in self.problem
    # however this isn't crazy important to get right, it's mostly aesthetics if we're being real
    forbidden_solutions: list[PolysatSolution]

    def __init__(self,
                 params: SolveParams):
        """
        constructor
        Parameters:
            params: a SolveParams object
        """
        self.data_dir = tlm_data_dir()
        self.solver_timeout = params.solver_timeout
        self.target_structure = GenericTypedStructure(bindings=[*params.bindings, *params.extraConnections],
                                                      types={i: 1 if i in params.nanoparticles else 0 for i in
                                                             params.topology.vertices()})

        # save solve parameter set
        self.input_params = params

        self.__soln_info_dirs = dict()

        if params.has_nanoparticles() and isinstance(params.nanoparticles, dict):
            raise Exception("Not yet tested!!! this will 100% crash something")
            # if the nanoparticle system has multiple nanoparticle types, we need to
            # add our dummy nanoparticle type and calcuate the number of nanoparticle types
            nNPT = max(params.nanoparticles.values())
            self.np_locations = {
                l: params.nanoparticles[l] if l in params.nanoparticles else self.nNPT
                for l in self.target_structure.vertices()
            }
            nNPT += 1  # add additional nanoparticle type for "no nanoparticle"
        elif isinstance(params.nanoparticles, list) and len(params.nanoparticles):
            self.np_locations = params.nanoparticles
            nNPT = 1
        else:
            self.np_locations = []
            nNPT = 0

        # nL = number of locations
        # self.internal_bindings = copy.deepcopy(params.topology.bindings_list)
        self.problem = PolycubeSATProblem(
            params.nS,
            (params.nC + 1) * 2,
            params.nDim,
            params.topology.num_vertices(),
            nNPT,
            params.torsion
        )
        self.forbidden_solutions = []

    def init(self, strict_counts=True):
        # generate basic constraints
        self.generate_constraints()

        # Solution must use all particles
        if strict_counts:
            self.problem.add_constraints_all_particles()

        self.problem.add_constraints_all_patches_except([0], [1])

        # A color cannot bind to itself
        # this is probably redundant w/ fix_color_interaction
        self.problem.add_constraints_no_self_complementarity()

        # Make sure color 0 binds to 1 and nothing else
        # self.fix_color_interaction(0, 1)

        # Fix interaction matrix, to avoid redundant solution
        # in other words, c2 should bind with c3, c4 with c5, etc.
        # if we don't do this then when we forbid a solution it will just find an alternative
        # arrangement of color vars

        # if we have torsional patches,
        if self.problem.nD == 3 and self.input_params.torsion:
            # i think this just saves time
            self.problem.add_constraints_fixed_blank_orientation()

        if self.input_params.multifarious_behavior != MultifariousTypeRequirement.NO_MULTIFARIOUS:
            if self.input_params.multifarious_behavior == MultifariousTypeRequirement.ANY_SHARED:
                self.problem.gen_req_mf_any_shared(self.input_params.topology)
            elif self.input_params.multifarious_behavior == MultifariousTypeRequirement.ALL_SHARED:
                self.problem.gen_req_mf_all_shared(self.input_params.topology)

        self.logger.info(f"Constructed {self.problem.describe_brief()}")

        # if any pre-forbidden solutions have been provided, try to forbid them
        if len(self.input_params.pre_forbid_solutions) > 0:
            self.logger.info(f"Forbidden solution paths have been provided!")
            for forbid_path in self.input_params.pre_forbid_solutions:
                self.forbid_solution_from_candidate(forbid_path)

    def forbid_solution_from_candidate(self, forbid_path: Path):
        """

        """
        # raw path is allowed but unusual
        if not forbid_path.is_dir():
            forbid_path = forbid_path.expanduser()
            if not forbid_path.is_dir():
                if not str(forbid_path).startswith("/"):
                    if (tlm_data_dir() / forbid_path).is_dir():
                        forbid_path = tlm_data_dir() / forbid_path
                    elif (__path__ / forbid_path).is_dir():  # completely unnessecary ngl
                        forbid_path = __path__ / forbid_path
                    else:
                        self.logger.warning(f"Could not find forbid-solution path {forbid_path}")
                        return
        if (forbid_path / "problem_info.json").is_file():
            with (forbid_path / "problem_info.json").open("r") as f:
                forbid_problem_info = json.load(f)
                if not self.matches_problem(forbid_problem_info):
                    self.logger.info(f"Forbidden solution at path {str(forbid_path)} does not match problem")
                    return
                else:
                    self.logger.info(f"Forbidding solution at path {str(forbid_path / 'solution.pickle')}...")
                    with (forbid_path / "solution.pickle").open("rb") as solution_in:
                        solns = pickle.load(solution_in)
                        for forbidden_solution in solns:
                            self.forbidden_solutions.append(forbidden_solution)
                            self.problem.forbidSolution(forbidden_solution)
        else:
            self.logger.warning(f"Directory {str(forbid_path)} does not contain a problem_info.json file!")

    def matches_problem(self, forbid_problem_info: dict) -> bool:
        if (self.problem.nL != forbid_problem_info["nL"]
                or self.problem.nC != forbid_problem_info["nC"]
                or self.problem.nS != forbid_problem_info["nS"]
                or self.problem.torsionalPatches != forbid_problem_info["torsion"]
                or self.problem.nNPT != forbid_problem_info["nNPT"]):
            return False
        else:
            return all([tuple(binding) in [tuple(b) for b in self.bindings] for binding in forbid_problem_info["bindings"]])

    def solution_info_dir(self, soln: PolysatSolution) -> Path:
        key = str(soln.rule)
        if key not in self.__soln_info_dirs:
            while True:  # every `while True` is a case for python do-while loops
                uid = uuid.uuid4().hex[:8]  # you will run out of memory on your machine before you hit the limit
                path = self.data_dir / f"soln_{uid}"
                try:
                    path.mkdir()
                    break
                except FileExistsError:
                    continue  # Try again with a new UUID
            self.__soln_info_dirs[key] = path
        return self.__soln_info_dirs[key]

    def is_bound_patch(self, p: int, l: int) -> int:
        """
        :param p: a patch index [0-6)
        :param l: a location index [0-self.nL)

        :return: true if patch p at location l is involved in binding, false otherwise

        """
        return (p, l) in self.bindings or (p, l) in self.bindings.values()

    def all_bindings(self) -> list:
        return list(self.bindings.keys()) + list(self.bindings.values())

    def generate_constraints(self):
        # do topology last
        # clauses iv & ix
        if self.input_params.crystal and self.input_params.segregate_bindings:
            generate_clauses(self.problem,
                             GenericTypedStructure(bindings=self.input_params.bindings,
                                                   types={
                                                       k: 1 if k in self.np_locations else 0 for k in
                                                       self.target_structure.vertices()
                                                   }))
            if self.input_params.segregate_colors:
                # if we're using seperate color sets for internal / external connection
                raise Exception()  # todo
        else:
            # old/default behavior
            generate_clauses(self.problem, self.input_params.topology)

    def run_sat_solve(self, nSolutions: int, solver_timeout: Union[int, None]) -> Union[SolverResponse,
                                                                                        list[PolysatSolution]]:
        # TODO: REFACTPOR
        solutions = []
        forbid_start = len(self.forbidden_solutions)  # number of solutions we have already run at start
        attempts_limit = self.input_params.maxAltTries * 1000  # TODO: SETTABLE
        while len(solutions) < nSolutions:
            result = self.run_glucose()

            # if not unsat or timeout
            if isinstance(result, PolysatSolution):
                self.save_solution(result) # save all solutions, good or bad
                if self.validate_solution(result):
                    solutions.append(result)
                    # if we're not going to test on lattice grid, save solution
                    # update: assume always use some version of lattice
                    # if not self.input_params.use_test_lattice():
                    #     self.save_solution(result)
                self.forbid_solution(result)
                if len(self.forbidden_solutions) - forbid_start - len(solutions) > attempts_limit:
                    self.logger.info(f"After {len(self.forbidden_solutions) - len(solutions) - forbid_start} we still "
                                     f"have found {len(solutions)} which does not also satisfy other topologies. Giving up...")
                    return SolverResponse.SOLN_TIMEOUT  # it's not technically, but go with it
            elif len(solutions) > 0:
                break
            else:
                return result

        return solutions

    def forbid_solution(self, solution: PolysatSolution):
        self.logger.info(f"Forbidding solution {solution.rule}...")
        ppname = self.problem.forbidSolution(solution)
        self.forbidden_solutions.append(solution)
        self.logger.info(f"Solution at path {self.solution_info_dir(solution)} forbidden by SAT problem part {ppname}. "
                         f"We now have forbidden {len(self.forbidden_solutions)} solutions")

    def validate_solution(self, soln: PolysatSolution) -> bool:
        #
        for forbidden in self.forbidden_solutions:
            # if solution rule is equivalent to a forbidden rule
            if soln.rule.to_canonical() == str(forbidden.rule):
                self.logger.info(f"Rule {soln.rule} is synonymous with an existing rule!")
                return False  # todo: should we check if nanoparticles match as well?
        return True

    def run_glucose(self, assumptions: list[int] = [], problem: Union[None, PolycubeSATProblem] = None) -> Union[
        SolverResponse, PolysatSolution]:
        """
        Uses Glucose to solve the SAT problem
        Returns either a PolysatSolution object or a SolverResponse code specifying what didn't work
        """
        if problem is None:
            problem = self.problem  # allow for trying problems other than the "standard" one
        formula = CNF(from_string=problem.output_cnf())
        tstart = datetime.datetime.now()
        self.logger.info(f"Starting solve with Glucose4, timeout {self.solver_timeout}")
        with Glucose4(bootstrap_with=formula.clauses) as m:
            # if the solver has a timeout specified
            if self.solver_timeout:
                timer = Timer(self.solver_timeout, interrupt, [m])
                timer.start()
                solved = m.solve_limited(assumptions=assumptions,
                                         expect_interrupt=True)
                timer.cancel()
            else:
                solved = m.solve(assumptions=assumptions)
            if solved:
                self.logger.info("Solved!")
                # pysat returns solution as a list of variables
                # which are "positive" if they're true and "negative" if they're
                # false.
                model = m.get_model()
                # we can pass the model directly to the PolysatSolution constructor because it will
                # just check for positive variables to be present
                return PolysatSolution(problem, frozenset(model), self.target_structure)
            else:
                # if the solve solver timed out
                if self.solver_timeout and (datetime.datetime.now() - tstart).seconds > self.solver_timeout:
                    return SolverResponse.SOLN_TIMEOUT
                # if the solver failed but didn't time out, conclude that the problem
                # is not satisfiable
                else:
                    return SolverResponse.SOLN_UNSAT

    def readSolutionAsRule(self, sol: str) -> PolycubesRule:
        colorCounter = 1
        # need to map B variable indeces to colors
        colorMap: dict[int, int] = {}

        # initialize blank rule
        rule = PolycubesRule(nS=self.nS)
        # color pairing variables (B) that are true
        B_vars = re.findall(r'B\((\d+),(\d+)\)', sol)
        # iterate values for B
        for c1, c2 in B_vars:  # color c1 binds with c2
            # print("Color {} binds with {}".format(c1, c2))
            assert (c1 not in colorMap or c2 not in colorMap)
            # map colors
            if int(c1) < 2 or int(c2) < 2:
                colorMap[c1] = 0
                colorMap[c2] = 0
            else:
                colorMap[c1] = colorCounter
                colorMap[c2] = -colorCounter
                colorCounter += 1

        # patch color (C) variable matches
        C_matches = re.findall(r'C\((\d+),(\d+),(\d+)\)', sol)
        for s, p, c in C_matches:  # Patch p on species s has color c
            patch_direction = RULE_ORDER[p]
            if not rule.particle(p).has_patch(patch_direction):
                rule.add_particle_patch(s, PolycubesPatch(uid=None,  # idx will be assigned in add method
                                                          color=colorMap[c],
                                                          direction=patch_direction,
                                                          orientation=get_orientation(p, 0)
                                                          )
                                        )
        oMatches = re.findall(r'O\((\d+),(\d+),(\d+)\)', sol)
        if len(oMatches) > 0:
            for s, p, o in oMatches:  # Patch on species l has orientation o
                # print("Patch {} on species {} has orientation {}".format(p, s, o))
                rule.particle(s).get_patch_by_diridx(p).set_align_rot(o)

        return rule

    def readSolutionRuleFromPath(self, path: Union[str, Path]) -> PolycubesRule:
        with open(path) as f:
            sol = f.read()

        return self.readSolutionAsRule(sol)

    def find_solution(self) -> Union[PolysatSolution, None, False]:
        """
        big method! go go go!
        """
        nTries = 0
        good_soln = None
        while nTries < self.input_params.maxAltTries and not good_soln:
            # run SAT
            result: list[PolysatSolution] = self.run_sat_solve(self.input_params.nSolutions,
                                                               self.input_params.solver_timeout)

            # check response
            if result == SolverResponse.SOLN_TIMEOUT:
                self.logger.info("SAT solver timed out!")
                break  # assume more iterations will not help (???)
            elif result == SolverResponse.SOLN_UNSAT:
                self.logger.info('Sorry, no solution')
                good_soln = False # set explicit false to indicate that problem is UNSAT
                break  # assume more iterations will not help
            elif len(result) > 0:  # found solutions!
                # if self.input_params.segregate_bindings:
                #     raise Exception("I am still working on this!")
                #     for i, soln in result:
                #         unit_cell_good = self.test_multifarious_finite_size(soln)
                #         if unit_cell_good:
                #             # create working copy of problem
                #             # todo: parallelize here
                #             next_itr_problem = copy.deepcopy(self.problem)
                #             next_itr_problem.hardcode_rule(soln)
                #
                #         else:
                #             self.problem.forbidSolution(soln)
                # else:
                self.logger.info(f"Found {len(result)} possible solutions:")
                # for now let's use Joakim's rule formulation
                self.logger.info("\n".join("\t" + soln.decRuleNew() for soln in result))
                # search for good solution
                # if we have specified crystalloyd test info
                good_soln = self.test_solutions(result)
                # if returned none, all solutions failed polycube tests, so forbid them and try again
                if good_soln is None:
                    for soln in result:
                        self.logger.info(f"Solution {soln.decRuleNew()} was found"
                                         f" invalid by Polycubes, forbidding.")
                        self.problem.forbidSolution(soln)
                    nTries += 1

            else:
                self.logger.error("Undefined problem of some sort!!!")
                break

        return good_soln

    def save_solution(self, soln: PolysatSolution):
        """
        Saves a polyccubes SAT solution to a file
        """
        soln_pickle_fn = self.solution_info_dir(soln) / f"solution.pickle"
        assert not soln_pickle_fn.is_file(), \
            f"File {str(soln_pickle_fn)} already exists - what are you DOING?"
        # write solution
        # todo: better approach than pickle
        with soln_pickle_fn.open("wb") as f:
            pickle.dump(soln, f)
        self.logger.info(f"Cached solution {soln.decRuleNew()} to file {str(soln_pickle_fn)}")
        # write solution info txt file
        # TODO: the way that it saves problems is misleading
        #  vis a vis `forbidden{n}` values, fix in future
        with (self.solution_info_dir(soln) / "soln_info.txt").open("w") as f:
            f.write(soln.describe() + "\n")
            f.write(f"Topology connections: {len(self.bindings)}\n\n")
            f.write(f"Solves problem:\n{self.problem.describe()}")
        # write problem info
        self.problem.save_problem_info(self.solution_info_dir(soln) / "problem_info.json",
                                       self.bindings)

    def test_solutions(self, sat_solutions: list[PolysatSolution]) -> Union[PolysatSolution, None]:
        """
        ok this is complecated because what makes a solution "good" depends on what kind of structure we're
        trying to design here
        """
        assert not self.input_params.crystal, "Use subclass CrystalPolysat to solve crystals!"
        # for each solution in our list
        for soln in sat_solutions:
            # if we're trying to design a multifarious system
            if self.input_params.multifarious_behavior != MultifariousTypeRequirement.NO_MULTIFARIOUS:
                # test rule in polycubes
                if self.test_multifarious_finite_size(soln):
                    return soln
            # if our design target isn't multifrious or a crystal but does have nanoparticles:
            elif self.input_params.has_nanoparticles():
                # if the solve target has nanoparticles, you can still use polycubes
                # but warn the user that this may create false positives
                if self.test_type_specific_finite_size(soln):
                    return soln
            else:
                # loop each solution in the results
                if self.test_finite_size(soln):
                    return soln
        return None



    def test_multifarious_finite_size(self, sat_solution: PolysatSolution) -> bool:
        """
        tricky
        """
        # i mean i guess odds are better than not that we are testing a crystalloyd, in a sense...
        mf_types = self.input_params.topology.iter_components()
        # there are sort of two ways i can approch this: either give the lattice model all the cubes for both structures
        # and then classify whatever it spits out, or give it the particles for each structure
        # doin the first one
        type_counts = [0 for _ in mf_types]
        n_iters = self.input_params.pc_analysis_hyperparams["num_search_iters"]
        polycube_records: list[tlm_data.TLMHistoryRecord] = libtlm.testFiniteSize(str(sat_solution.rule),
                                                                                  n_iters)

        for record_idx, record in enumerate(polycube_records):
            assert record.numPolycubes() == 1, "Record has more or less than one polycube!"
            top_idx = -1
            for i, topology in enumerate(mf_types):
                # use the same method we use to identify a unit cell match, even if this isn't technically a crystalloyd
                pc = toPolycube(sat_solution.rule, record.getPolycubes()[0])
                try:
                    rot_idx, coord_matrix, uids, ptypes = next(match_crystalloyd_to_unit_cell(pc,
                                                                                              topology,
                                                                                              sat_solution.nanoparticle_map))
                    # todo: info message w/ return values?
                    top_idx = i
                    break
                except StopIteration:
                    pass
            if top_idx == -1:
                self.logger.info(f"Finite-size unit cell test of rule `{sat_solution.rule} found a topology that ")
                # todo: should it keep going to see if it finds more misassemblies or is one too many?
                # since we are about to forbid solution and move on, th
                self.solution_info_dir(sat_solution) / f"misassembly.json"
                return False
        self.logger.info(f"Finite-size testing of {sat_solution.rule} found {type_counts} of our topologies over "
                         f"{self.input_params.pc_analysis_hyperparams['num_search_iters']}.")
        return True

    def test_type_specific_finite_size(self, sat_solution: PolysatSolution) -> bool:
        """

        """
        raise

    def test_finite_size(self, sat_solution: PolysatSolution) -> bool:
        """
        this is old code but it should still work
        """
        self.logger.info(f"Testing rule {sat_solution.decRuleNew()}")
        tlm_params = libtlm.TLMParameters(
            self.input_params.torsion,
            False,
            0,
            0,
            str(sat_solution.rule),
            sat_solution.type_counts(),  # unused
            self.problem.nL + 1, # we need to offer one additional step after nL to see if it keeps going
            1
        )
        # hate this naming convention
        cube_type_types = sat_solution.nanoparticle_map
        for ct in sat_solution.rule:
            if ct.type_id() not in cube_type_types:
                cube_type_types[ct.type_id()] = 0
        # todo: make this less hardcoded
        # these are more properly hp
        batches = self.input_params.pc_analysis_hyperparams["pcs_n_batches"]
        batch_size = self.input_params.pc_analysis_hyperparams["pcs_batch_size"]
        for i in range(batches):
            records = run_tlm(i, tlm_params, batch_size, self.solution_info_dir(sat_solution), 1)
            for ri, record in enumerate(records):
                # we still want to cache polycube even if it's going to be doa
                pc = toPolycube(sat_solution.rule, record.getPolycubes()[0])
                pc.save_polycube(self.solution_info_dir(sat_solution) / f"polycube_{i*batch_size+ri}.json")
                # low-hanging fruit: if the number of cubes doesn't line up with any topology, return false
                if record.numCubeInstances() not in [top.num_vertices() for top in self.input_params.topologies]:
                    return False
                else:
                    # convert polycube now to make faster
                    pc_as_typed = pc.to_typed_structure(cube_type_types)
                    found_match = False
                    # iter topologies (for multifarious targets)
                    for top in self.input_params.topologies:
                        # if a mapping exists from top <-> pc
                        # (again, screen w/ num verts since the homomorphism code will eat a lot of time if the sizes are different)
                        if top.num_vertices() == pc_as_typed.num_vertices() and top.homomorphism(pc_as_typed):
                            found_match = True
                            break
                    if not found_match:
                        return False
        return True


        #
        # if libtlm.isBoundedAndDeterministic(sat_solution.decRuleNew(),
        #                                     self.input_params.torsion,
        #                                     100):
        #     self.logger.info(f"{sat_solution.decRuleNew()} works!! We can stop now!")
        #     return True
        # else:
        #     self.logger.info(f"{sat_solution.decRuleNew()} found to be unbounded and/or nondeterministic.")
        #     return False

    @functools.cached_property
    def targets_as_typed(self) -> list[TypedStructure]:
        return list(self.target_as_typed_structure().iter_components())

    def target_as_typed_structure(self) -> TypedStructure:
        """
        returns the target structure as a typed structure
        """
        if len(self.np_locations) == 0:
            # if no nanoparticle info is present
            # edge case
            return GenericTypedStructure(graph=self.target_structure.graph,
                                         types={k: 0 for k in self.target_structure.vertices()})
        elif isinstance(self.np_locations, list):
            # if we have one nanoparticle type
            return GenericTypedStructure(graph=self.target_structure.graph,
                                         types={
                                             k: 1 if k in self.np_locations else 0 for k in
                                             self.target_structure.vertices()
                                         })
        else:
            assert isinstance(self.np_locations, dict)
            # todo test
            np_type_map = {k: -1 for k in self.target_structure.vertices()}
            # todo: make sure "type -1" doesn't cause problems
            np_type_map.update(self.np_locations)
            return GenericTypedStructure(graph=self.target_structure.graph,
                                         types={k: 0 for k in self.target_structure.vertices()})


def generate_clauses(problem: PolycubeSATProblem, topology: TypedStructure, skip_clauses: list[str] = []):
    """
            Adds all constraints from the original paper (clauses i - x)
            Returns: a list of CNF clauses expressing SAT clauses i - x in the original paper, plus
            a few more that are implicit in the paper)

            """
    # todo: incorporate skip_clauses in a way that doesn't suck
    problem.generate_BCO()
    # clause i
    if "B" not in skip_clauses:
        problem.gen_legal_color_bindings()
    # clause ii
    if "C" not in skip_clauses:
        problem.gen_legal_species_coloring()
    # ADD CRYSTAL and COLORS:
    # clause ???
    if "F" not in skip_clauses:
        problem.gen_legal_patch_coloring()
    # clause ???
    if "A" not in skip_clauses:
        problem.gen_legal_position_patch_orientation()
    # clause ix?
    if "D" not in skip_clauses:
        problem.gen_hard_code_orientations()
    # clause iii
    if "P" not in skip_clauses:
        problem.gen_legal_species_placement()
    # clause iix, or viii if you know how to use roman numerals correctly
    if "O" not in skip_clauses:
        problem.gen_legal_species_orientation()
    # clause v
    if "rotC" not in skip_clauses:
        problem.gen_legal_species_position_coloring()
    # clause x
    if "rotO" not in skip_clauses:
        problem.gen_legal_species_patch_orientation()
    # 2D-specific mystery clauses
    if "lock" not in skip_clauses:
        problem.gen_lock_patch_orientation()
    if "BFi" not in skip_clauses and "DAi" not in skip_clauses:
        problem.gen_forms_desired_structure(topology.bindings_list)
    if "E" not in skip_clauses:
        problem.fix_empties(topology)

    # if problem has nanoparticles
    if problem.nNPT > 0:  # use this instead of num_particle_types because negative design stuff, trust me it's useful
        # if the nanoparticle data is provided as a list, that means single np type
        if topology.num_particle_types() <= 2:
            problem.gen_nanoparticle_singleparticle([k for k, v in topology.get_particle_types().items() if v != 0])
        else:
            problem.gen_nanoparticle_multiparticle({k: v for k, v in topology.get_particle_types().items() if v != 0})
