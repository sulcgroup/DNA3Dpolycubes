from __future__ import annotations

import json
import logging
import os
import re

import pandas as pd

from .sat_problem import SATProblem
from ...design.sat.polycube_sat_problem import PolycubeSATProblem
from ...polycubeutil.polycube_structure import PolycubeStructure
from ...polycubeutil.cube import PolycubesStructureCube
from ...structure import *
from ...util import get_output_dir

from ..solve_utils import compute_coordinates, to_xyz, rot_idx_to_quaternion, quaternion_inverse


class SATSolution:
    # list of variables that are true
    raw: frozenset[int]
    # each key is a var type, value is a tuple where first item is variable parameters,
    # second is variable number (in raw)
    vars_by_type: dict[str, set[ tuple[int, tuple[int, ...]]]]

    # variables mapped for var number to var info
    var_types: dict[int, tuple[str, tuple[int, ...]]]

    def __init__(self, problem: SATProblem, sat_results: frozenset[int]):
        self.raw = sat_results
        self.vars_by_type = dict()
        self.var_types = dict()
        for vname, vnum in problem.list_vars():
            if vnum in sat_results:
                m = re.match(r"([A-Za-z]+)\(([\d\s,]*)\)", vname)
                if m:
                    label: str = m.group(1)  # 'FooBar'
                    numbers: tuple[int, ...] = tuple([int(x) for x in re.findall(r"\d+", m.group(2))])
                    if label not in self.vars_by_type:
                        self.vars_by_type[label] = set()
                    self.vars_by_type[label].add( (vnum, numbers) )
                    self.var_types[vnum] = (label, numbers)
                else:
                    raise Exception(f"Invalid sat variable name {vname}")

    def get_vars(self, var_name: str) -> list[int]:
        return list([
            i for i, _ in self.vars_by_type[var_name]
        ])

    def dump_vars(self, fp: Path):
        """
        saves all variables that are true to a text file
        :param fp: path to file to save
        """
        with fp.open("w") as f:
            for v in self.raw:
                if v > 0:
                    f.write(f"{v} ")

    def dump_var_names(self, fp: Path):
        """
        saves solution to json form
        """
        vars_list = [
            [self.var_types[v][0], list(self.var_types[v][1])] for v in self.raw if v > 0
            ]
        with fp.open("w") as f:
            json.dump(vars_list, f)


class PolysatSolution(SATSolution):

    num_species: int
    num_colors: int

    def get_C_vars(self) -> list[int]:
        return self.get_vars("C")

    def get_O_vars(self) -> list[int]:
        return self.get_vars("O")

    # var C(s,p,c) = patch p on species s has color c
    C_vars: list[str] = property(get_C_vars)
    # var O(s,p,o) = patch p on species s has orientation o
    O_vars: list[str] = property(get_O_vars)

    rule: PolycubesRule

    # spacial map is a list of tuples where the first value is a particle species, the second is a rotation
    # todo: replace with typedstructure
    spacial_map: list[tuple[int, int]]
    top: Structure
    # map of cube types -> nanoparticle types
    nanoparticle_map: dict

    coordinate_map: Union[dict[int, np.array], None]

    def __init__(self, sat_problem: PolycubeSATProblem, sat_results: frozenset[int], target_structure: Structure):
        """
        Constructor.
        Args:
            solver: an instance of Polysat, used to contextualize SAT solver data
            sat_results: a frozenset object containing variables (as numbers) in the solution that are True
        """
        # save solver parameters
        super().__init__(sat_problem, sat_results)
        self.num_species = sat_problem.nS
        self.num_colors = int(sat_problem.nC / 2 - 1)

        self.top = target_structure

        # loop through variables
        colorMap = {}
        self.rule = PolycubesRule(nS=sat_problem.nS)
        colorCounter = 1
        # species identity, rotation tuple for each position in nL
        self.spatial_map = [(-1, -1) for _ in range(sat_problem.nL)]
        self.nanoparticle_map = {}
        # check that we have P(l,s,r) for some s,r for all L
        assert all([
            any([
                sat_problem.P(l, s, r) in sat_results
                for s, r in itertools.product(range(sat_problem.nS), range(sat_problem.nR))])
            for l in
            range(sat_problem.nL)]), "Missing some spacial solution informatiuon ( P(l,s,r) has no s,r for some l)"

        # iter sat problem variables
        # TODO: read these from SATSolution
        for vname, vnum in sat_problem.list_vars():
            # if variable is True in solution
            if vnum in sat_results:
                # check if this variable is a color match spec
                m = re.match(r"B\((\d*),\s*(\d*)\)", vname)
                if m:
                    c1, c2 = m.groups()
                    assert (c1 not in colorMap or c2 not in colorMap)
                    if int(c1) < 2 or int(c2) < 2:
                        colorMap[int(c1)] = 0
                        colorMap[int(c2)] = 0
                    else:
                        colorMap[int(c1)] = colorCounter
                        colorMap[int(c2)] = -colorCounter
                        colorCounter += 1
                    continue

                # check if this variable is a position specifier
                m = re.match(r"P\((\d*),\s*(\d*),\s*(\d*)\)", vname)
                if m:
                    location, species, rotation = m.groups()
                    species = int(species)
                    self.spatial_map[int(location)] = (species, int(rotation))
                    continue

                # check for nanoparticles
                if sat_problem.nNPT > 0:
                    # if there's 1 nanoparticle type
                    if sat_problem.nNPT == 1:
                        m = re.match(r"N\((\d*)\)", vname)
                        if m:
                            self.nanoparticle_map[int(m.groups()[0])] = 1
                            continue
                    # if there's more than 1 nanoparticle type
                    else:
                        m = re.match(r"N\((\d*),\s*(\d*)\)", vname)
                        if m:
                            species, nptype = m.groups()
                            self.nanoparticle_map[int(species)] = int(nptype)
                            continue

        if sat_problem.nNPT !=0 and len(self.nanoparticle_map) == 0:
            raise MalformedSolutionException()

        if not all(species != -1 and rotation != -1 for species, rotation in self.spatial_map):
            raise MalformedSolutionException()

        # apply patch color data
        for _, (s, p, c) in self.vars_by_type["C"]:
            if colorMap[c]:
                particle_type_idx = int(s)
                patch_direction = RULE_ORDER[int(p)]
                if not self.rule.particle(particle_type_idx).has_patch(patch_direction):
                    self.rule.add_particle_patch(particle_type_idx, PolycubesPatch(
                        uid=None,  # idx will be assigned in add method
                        color=colorMap[c],
                        direction=patch_direction,
                        orientation=get_orientation(int(p), 0)
                    ))

        # if applicable, apply patch orientation data
        if "O" in self.vars_by_type:
            for _, (s, p, o) in self.vars_by_type["O"]:  # Patch on species l has orientation o
                # print("Patch {} on species {} has orientation {}".format(p, s, o))
                species_idx = int(s)
                direction_idx = int(p)
                if self.rule.particle(species_idx).has_patch(direction_idx):
                    rotation = int(o)
                    self.rule.particle(species_idx).get_patch_by_diridx(direction_idx).set_align_rot(rotation)

        # construct map of location indexes to coordinates in 3-space....
        if not target_structure.is_multifarious() and not target_structure.is_crystal():
            self.coord_map = compute_coordinates(target_structure.bindings_list)
        else:
            self.coord_map = None

    def decRuleNew(self):
        return str(self.rule)

    def hasNanoparticles(self):
        return len(self.nanoparticle_map) > 0

    def numNanoparticleTypes(self):
        return np.unique(self.nanoparticle_map.values()).size

    def printToLog(self, logger: logging.Logger = logging.root):
        logger.info("-----------------")
        logger.info(f"\tNum Species: {self.num_species}")
        logger.info(f"\tNum Colors: {self.num_colors}")
        logger.info(f"\tRule: {self.decRuleNew()}")
        logger.info(f"\tSpatial Map: {self.spatial_map}")
        if self.hasNanoparticles():
            logger.info(f"\tNum Nanoparticle Types: {self.numNanoparticleTypes()}")
            for nptype in range(self.numNanoparticleTypes()):
                logger.info(
                    f"\t\tType {nptype + 1}: {[i for i, ct in enumerate(self.rule.particles()) if self.nanoparticle_map[i] == nptype + 1]}")

        logger.info("-----------------")

    def type_counts(self) -> list[int]:
        tyoe_counts = [0 for _ in range(self.num_species)]
        for species, _ in self.spatial_map:
            tyoe_counts[species] += 1
        return tyoe_counts

    def to_polycube(self) -> PolycubeStructure:
        cubes = []
        for i in range(len(self.spatial_map)):
            ctid, ctrot = self.spatial_map[i]
            cubes.append(
                PolycubesStructureCube(
                    uid=ctid,
                    cube_position=self.coord_map[i],
                    cube_rotation=ctrot,
                    cube_type=self.rule.particle(ctid),
                    state=[True],  # ignore for now TODO come back
                )
            )

        return PolycubeStructure(
            rule=self.rule,
            cubes=cubes
        )

    def exportScene(self, modelname="scene"):
        if modelname.startswith("/"):
            p = Path(modelname)
        elif not isinstance(modelname, Path):
            if modelname.startswith("~"):
                p = Path(os.path.expanduser(modelname))
            else:
                p = Path(modelname)

        elif modelname.endswith(".json"):
            p = get_output_dir() / "SAT" / modelname
        else:
            p = get_output_dir() / "SAT" / f"{modelname}_{self.num_species}S_{self.num_colors}C.json"
        try:
            self.to_polycube().save_polycube(p)
        except FileNotFoundError as e:
            print(e.strerror)

    def cubeToJSON(self, i: int) -> dict:
        """
        Exports spacial data on a cube
        Args:
            i: the location index (in self.spacial_map) of the cube to export

        Returns:
            a dict of data for the spacial location of a cube
        """
        ctid, ctrot = self.spatial_map[i]

        return {
            "position": to_xyz(self.coord_map[i]),
            "rotation": {
                k: float(v) for k, v in zip(
                    ["w", "x", "y", "z"],
                    quaternion_inverse(rot_idx_to_quaternion(ctrot)))
            },
            "state": [True, *([False] * 12)],  # ignore for now TODO come back to
            "type": ctid
        }

    def has_coord_map(self):
        return self.coord_map is not None

    def describe(self) -> str:
        """
        outputs human-readable description
        """
        description = f"Num Species: {self.num_species}" \
                      f"\nNum Colors: {self.num_colors}" \
                      f"\nRule: {self.decRuleNew()}" \
                      f"\n Num locations: {len(self.spatial_map)}"
        if self.hasNanoparticles():
            # todo: may need to alter code for systems w/ multiple nanoparticle types
            description += f"\nNum Nanoparticle Types: {self.numNanoparticleTypes()}\n"
            for nptype in range(self.numNanoparticleTypes()):
                description += f"\n\tType {nptype+1}: {[i for i, ct in enumerate(self.rule.particles()) if i in self.nanoparticle_map and self.nanoparticle_map[i] == nptype+1]}\n"
            if self.top.is_multifarious():
                # species locations
                spatial_map: list[int] = [i for i, _ in self.spatial_map]
                structures = list(self.top.iter_components())
                if len(structures) != 2:
                    raise ValueError("Multifarious not implemented for num structures != 2")
                c1, c2 = structures
                type_counts = [[0 for _ in range(self.num_species)] for c in [c1, c2]]
                for location, species in enumerate(spatial_map):
                    if location in c1.vertices():
                        type_counts[0][species] += 1
                    elif location in c2.vertices():
                        type_counts[1][species] += 1
                    else:
                        raise Exception("Somehow found location that isn't in either substructure of multifarious topology")
                description += str(f"Type Counts {type_counts}\n")
            locs_df = pd.DataFrame(self.spatial_map, columns=["Species", "Rotation"])
            locs_df.index.set_names("Location", inplace=True)
            description += "Location map:\n"
            description += locs_df.to_string()
        return description

    def mf_type_counts(self, tops: Iterable[Structure]) -> list[list[int]]:

        spatial_map: list[int] = [i for i, _ in self.spatial_map]
        type_counts = [[0 for _ in range(self.num_species)] for top in tops]
        for location, species in enumerate(spatial_map):
            for i, top in enumerate(tops):
                if location in top.vertices():
                    type_counts[i][species] += 1
                elif location in top.vertices():
                    type_counts[i][species] += 1
                else:
                    raise Exception()
        # todo: more checks
        assert all(sum(type_counts[i]) == top.num_vertices() for i, top in enumerate(tops))
        return type_counts


class MalformedSolutionException(BaseException):
    """
    the main purpose of this is to let the constructor throw catchable exceptions
    """
    def __init__(self, msg: str = ""):
        self.msg = msg

    def __str__(self):
        return self.msg
