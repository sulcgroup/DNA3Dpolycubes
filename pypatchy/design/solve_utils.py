import dataclasses
import itertools
import logging
from datetime import datetime
from typing import Union, Generator

import libtlm
import networkx as nx
import numpy as np

from scipy.spatial.transform import Rotation as R

from ..polycubeutil.cube import PolycubesStructureCube
from ..polycubeutil.polycube_structure import PolycubeStructure
from ..structure import TypedStructure, Structure
from ..util import get_log_dir, enumerateRotations, getSignedAngle
from ..polycubeutil.polycube_util import rotation_mapping_to_matrix
from ..polycubeutil.polycubesRule import RULE_ORDER, PolycubesRule
# DON'T IMPORT ANYTHING FROM design_rule !!!

logging.basicConfig(filename=get_log_dir() / "SAT" / datetime.now().strftime("log_%Y-%m-%d-%H:%M.txt"))
logging.root.setLevel(logging.INFO)

def to_diridx(x: int) -> libtlm.DirIdx:
    return [libtlm.LEFT, libtlm.RIGHT,
            libtlm.BOTTOM, libtlm.TOP,
            libtlm.BACK, libtlm.FRONT][x]


def setup_logger(logger_name, file_path=None):
    if file_path is None:
        file_path = get_log_dir() / "SAT" / f"{logger_name}.log"
    logger = logging.getLogger(str(logger_name))
    logger.setLevel(logging.INFO)

    # create a file handler
    handler = logging.FileHandler(file_path)
    handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(handler)

    return logger


## Utility functions
def patchRotToVec(i: int, rot):
    """ Get vector indicating patch rotation, given rotation state and face index

    Args:
        i (int): Index of patch [0...5]
        rot (int): Rotation state [0,1,2,3] = North , East, South, West

    Returns:
        vector: Patch rotation vector
    """
    v = getFACE_ROTATIONS()[i]
    axis = RULE_ORDER[i]
    angle = rot * np.pi / 2
    r = R.from_rotvec(angle * axis)
    return r.apply(v).round()


def getFACE_ROTATIONS():
    """

    """
    return [
        np.array([0, -1, 0]), np.array([0, 1, 0]),
        np.array([0, 0, -1]), np.array([0, 0, 1]),
        np.array([-1, 0, 0]), np.array([1, 0, 0])
    ]


def getFlatFaceRot():
    """
    currently unused? for 2D?
    """
    return [1, 1, 2, 0, 0, 0]


def getIndexOf(elem, array):
    """
    like 90% sure this can be replaced with some kinda call to `find`
    """
    for i, e in enumerate(array):
        if (e == elem).all():
            return i
    return -1


def coordsFromFile(path):
    with open(path) as f:
        return [[int(c) for c in line.strip('()\n').split(',')] for line in f]


def patchVecToRot(i: int, v: np.array) -> int:
    """ Get rotation state, given face index and patch rotation vector

    :param i: Index of patch [0...5]
    :param v: Patch rotation vector

    :return: Rotation state [0,1,2,3] = North , East, South, West
    """
    angle = getSignedAngle(
        getFACE_ROTATIONS()[i],
        v,
        RULE_ORDER[i]
    )
    return int((angle * (2 / np.pi) + 4) % 4)


def forbidden_symmetries():
    # begin with an arbitrary reflective symmetry
    asym = {
        0: 1,
        1: 0,
        2: 2,
        3: 3,
        4: 4,
        5: 5
    }
    # apply rotations to our seed symmetry
    forbidden = [
        {
            idx_key: rot_dict[asym[idx_key]] for idx_key in rot_dict
        }
        for sym_key, rot_dict in enumerateRotations().items()
    ]
    # perform tests (TODO: remove once I'm more confident)
    assert len(forbidden) == 24
    assert all([
        not any([
            all(
                f[key] == sym[key] for key in range(6)
            ) for sym in enumerateRotations().values()]
        ) for f in forbidden]
    )
    return forbidden


def rotation_matrix_to_quaternion(rotation):
    """
    chatGPT wrote this code
    Convert a 3x3 rotation matrix to a quaternion.
    Is there a neater way to do this? unclear!!!
    """

    assert rotation.shape == (3, 3)
    trace = np.trace(rotation)

    if trace > 0:
        s = 0.5 / np.sqrt(trace + 1.0)
        w = 0.25 / s
        x = (rotation[2, 1] - rotation[1, 2]) * s
        y = (rotation[0, 2] - rotation[2, 0]) * s
        z = (rotation[1, 0] - rotation[0, 1]) * s
    else:
        if rotation[0, 0] > rotation[1, 1] and rotation[0, 0] > rotation[2, 2]:
            s = 2.0 * np.sqrt(1.0 + rotation[0, 0] - rotation[1, 1] - rotation[2, 2])
            w = (rotation[2, 1] - rotation[1, 2]) / s
            x = 0.25 * s
            y = (rotation[0, 1] + rotation[1, 0]) / s
            z = (rotation[0, 2] + rotation[2, 0]) / s
        elif rotation[1, 1] > rotation[2, 2]:
            s = 2.0 * np.sqrt(1.0 + rotation[1, 1] - rotation[0, 0] - rotation[2, 2])
            w = (rotation[0, 2] - rotation[2, 0]) / s
            x = (rotation[0, 1] + rotation[1, 0]) / s
            y = 0.25 * s
            z = (rotation[1, 2] + rotation[2, 1]) / s
        else:
            s = 2.0 * np.sqrt(1.0 + rotation[2, 2] - rotation[0, 0] - rotation[1, 1])
            w = (rotation[1, 0] - rotation[0, 1]) / s
            x = (rotation[0, 2] + rotation[2, 0]) / s
            y = (rotation[1, 2] + rotation[2, 1]) / s
            z = 0.25 * s

    q = np.array([w, x, y, z])
    return q


def quaternion_inverse(q):
    # q is a numpy array of four elements [w, x, y, z]

    norm_squared = np.dot(q, q)  # Compute the square of the norm

    if norm_squared == 0:
        raise ValueError("The input quaternion is zero, which cannot be inverted.")

    # Compute the conjugate of the quaternion
    conjugate = np.array([q[0], -q[1], -q[2], -q[3]])

    # Compute the inverse by dividing the conjugate by the norm squared
    inverse = conjugate / norm_squared

    return inverse


def rot_idx_to_quaternion(idx: int) -> np.ndarray:
    return rot_map_to_quat(enumerateRotations()[idx])


def rot_map_to_quat(rot_map: dict[int, int]) -> np.ndarray:
    return rotation_matrix_to_quaternion(rotation_mapping_to_matrix(rot_map))


def countParticlesAndBindings(topology):
    """
    Returns the number of bindings and particles in a topology
    Doesn't take into account external bindings (for crystals) or nanoparticles
    Args:
        topology:

    Returns: a tuple of the number of particles and the number of bindings
    """
    pidsa = [x[0] for x in topology]
    pidsb = [x[2] for x in topology]
    particles = pidsa + pidsb
    return max(particles) + 1, len(topology)


def parseHexRule(hexRule):
    ruleset = []
    faces = []
    for i in range(0, len(hexRule), 2):
        if i % 12 == 0 and i != 0:
            ruleset.append(faces)
            faces = []
        face_hex = hexRule[i:i + 2]
        face_int = int(face_hex, 16)
        face_bin = bin(face_int)[2:].zfill(8)
        face = {}
        sign = int(face_bin[0], 2)
        face['color'] = int(face_bin[1:6], 2) * (-1 if sign else 1)
        face['orientation'] = int(face_bin[6:8], 2)
        face['conditional'] = ""  # hex conditionals can't be expressed in hex (TODO: yet?)
        faces.append(face)
    ruleset.append(faces)
    return ruleset


def ruleToHex(ruleset):
    hexRule = ''
    for rule in ruleset:
        for face in rule:
            sign = bin(face['color'] < 0)[2:]
            color = bin(abs(face['color']))[2:].zfill(5)
            orientation = bin(abs(face['orientation']))[2:].zfill(2)
            binStr = sign + color + orientation
            hexStr = hex(int(binStr, 2))[2:].zfill(2)
            hexRule += hexStr
    return hexRule


def to_xyz(vector):
    return {k: int(v) for k, v in zip(["x", "y", "z"], vector)}


def compute_coordinates(topology: frozenset) -> dict[int, np.array]:
    """
    Args:
        topology:

    Returns:

    """
    top_queue = list(topology)

    # Initialize the first coordinate at (0,0,0) and create dict to hold coordinates
    coord_dict = {top_queue[0][0]: np.zeros((3,))}

    # top_queue = top_queue[1:]  # pop zeroth element

    giveupcount = 0

    while len(top_queue) > 0 and giveupcount < 1e7:
        queue_len = len(top_queue)
        for j in range(queue_len):
            i = queue_len - j - 1
            loc1, dir1, loc2, dir2 = top_queue[i]
            assert -1 < dir1 < len(RULE_ORDER), f"Invalid direction index {dir1}"
            assert -1 < dir2 < len(RULE_ORDER), f"Invalid direction index {dir2}"
            if loc1 in coord_dict and loc2 in coord_dict:  # unclear how this happened
                top_queue = top_queue[:i] + top_queue[i + 1:]
                continue
            # assert loc1 not in coord_dict or loc2 not in coord_dict, \
            #     f"Both locations {loc1} and {loc2} have already been handled!"
            assert (RULE_ORDER[dir1] + RULE_ORDER[dir2] == np.zeros((3,))).all(), \
                f"Particles at {loc1} and {loc2} are not bound by opposite directional patches!!! " \
                f"Directions are {RULE_ORDER[dir1]} and {RULE_ORDER[dir2]}"
            if loc1 in coord_dict:
                # Compute the new coordinate by adding the direction vector to the current coordinate
                coord_dict[loc2] = coord_dict[loc1] + RULE_ORDER[dir1]
            elif loc2 in coord_dict:
                # Compute the new coordinate by adding the direction vector to the current coordinate
                coord_dict[loc1] = coord_dict[loc2] + RULE_ORDER[dir2]
            else:
                continue
            top_queue = top_queue[:i] + top_queue[i + 1:]

    if len(top_queue) > 0:
        raise Exception("Unable to construct a coordinate map! Perhaps topology is not connected?")

    assert len(coord_dict) == len(set(itertools.chain.from_iterable([(b[0], b[2]) for b in topology])))

    # Ensure minimum distance of 1 between structures by sorting and offsetting coordinates
    sorted_locs = sorted(coord_dict.keys())
    for i in range(1, len(sorted_locs)):
        diff = coord_dict[sorted_locs[i]] - coord_dict[sorted_locs[i - 1]]
        if np.abs(diff).max() < 1:
            coord_dict[sorted_locs[i]] += np.sign(diff) * (1 - np.abs(diff))

    return coord_dict


def build_graphs(topology):
    """
    warning: the following function was written by ChatGPT
    Args:
        topology:

    Returns:

    """
    # Initialize an empty graph and list to hold all graphs
    G = nx.Graph()
    graph_list = []

    for binding in topology:
        loc1, _, loc2, _ = binding
        # Add an edge to the graph for this binding
        G.add_edge(loc1, loc2)

    # Generate subgraphs for each connected component
    for component in nx.connected_components(G):
        subgraph = G.subgraph(component).copy()
        graph_list.append(subgraph)

    return graph_list


def toPolycube(rule: PolycubesRule, pc: libtlm.TLMPolycube) -> PolycubeStructure:
    return PolycubeStructure(rule=rule, cubes=[
        {
            "uid": tlm_cube_info.getUID(),
            "cube_position": tlm_cube_info.getPosition(),
            "cube_rotation": tlm_cube_info.getRotation(),
            "state": tlm_cube_info.getState(),
            "cube_type": rule.particle(tlm_cube_info.getTypeID())
            # "personalName": f"cube_{tlm_cube_info.getUID()}"
        } for tlm_cube_info in pc.getCubes()
    ])


def match_unit_cell_at_pt(
        crystalloyd: PolycubeStructure,
        tiling_structure_dimensions: tuple[int, int, int],
        coord_matrix: np.ndarray,
        np_species_map: dict[int: int],
        ptypes: np.ndarray,
        uids: np.ndarray,
        rmap: dict[int: int],
        rot_idx: int,
        tran: np.ndarray) -> Union[tuple[int, np.ndarray, np.ndarray, np.ndarray], bool]:
    """
    attempts to find a rotation that maps the unit cell (with coordinates defined in coord_matrix)
    to the crystalloyd, with the provided translation
    """
    is_type_map_good = True

    # need to make this return a generator as well bc there may be multiple valid rotmaps
    # at a given translation and we need to check them all
    # compute rotation matrix from mapping
    rotmat = rotation_mapping_to_matrix(rmap)
    # apply transformation to coordinate matrix
    cm = ((coord_matrix - (np.array(tiling_structure_dimensions) // 2)) @ rotmat) + np.array(tran)[np.newaxis, :]
    # if the polycube has a cube at every position which map to the tiling structure
    pos_map_good = all([crystalloyd.has_cube_at(pos) for pos in cm])
    if pos_map_good:
        # now we need to test that it's possible to create a type-mapping from tiling_structure onto self
        # for type map:
        # start with assumption that map is good

        # check each particle in the tile
        for (pos, tile_vert_type_id) in zip(cm, ptypes):
            cube = crystalloyd.cubeAtPosition(pos)
            # if the types don't map correctly
            if np_species_map[cube.get_cube_type().type_id()] != tile_vert_type_id:
                is_type_map_good = False
                break
        # if all particle types mapped properly, we found a good one!
        if is_type_map_good:
            # this return is such a mess i'm almost tempted to make it a real class
            return rot_idx, cm, uids, ptypes
    return False


def match_unit_cell_from_igraph_mapping(
    crystalloyd: PolycubeStructure,
    unit_cell: TypedStructure,
    mapping: list[int],
    np_species_map: dict[int, int]) -> Union[tuple[int, np.ndarray, np.ndarray, np.ndarray], bool]:
    """
    Given a vertex mapping from crystalloyd (subgraph) to unit cell (supergraph),
    verify spatial alignment and type mapping, and return the closest discrete rotation index.
    Code by chatGPT

    Returns:
        (rot_idx, unit_coords, unit_uids, unit_types) if valid, else False
    """

    # --- Convert igraph indices back to node labels ---
    g_crys = crystalloyd.to_igraph()
    g_unit = unit_cell.to_igraph()
    crys_labels = g_crys.vs['uid']
    unit_labels = g_unit.vs['uid']

    crys_nodes = [crys_labels[i] for i in range(len(mapping))]

    # --- Get positions ---

    crys_coords = np.array([crystalloyd.get_particle(uid).position() for uid in crys_nodes])
    unit_coords, unit_nodes, _ = unit_cell.matrix_etc()

    # --- Kabsch alignment ---
    A = crys_coords - crys_coords.mean(axis=0)
    B = unit_coords - unit_coords.mean(axis=0)
    H = A.T @ B
    U, _, Vt = np.linalg.svd(H)
    R_kabsch = Vt.T @ U.T
    if np.linalg.det(R_kabsch) < 0:
        Vt[2, :] *= -1
        R_kabsch = Vt.T @ U.T

    # --- Match to known discrete rotations ---
    best_rot_idx = None
    best_error = np.inf

    for rot_idx, rmap in enumerateRotations().items():
        R_discrete = rotation_mapping_to_matrix(rmap)
        error = np.linalg.norm(R_kabsch - R_discrete)
        if error < best_error:
            best_rot_idx = rot_idx
            best_error = error

    if best_error > 1e-2:
        return False  # No close-enough rotation found

    # --- Check alignment after applying the best discrete rotation ---
    A_discrete = crys_coords - crys_coords.mean(axis=0)
    aligned = (A_discrete @ rotation_mapping_to_matrix(enumerateRotations()[best_rot_idx])) + unit_coords.mean(axis=0)

    if not np.allclose(aligned, unit_coords, atol=1e-3):
        return False

    # --- Check type mapping ---
    for crys_uid, unit_uid in zip(crys_nodes, unit_nodes):
        crys_type = crystalloyd.get_particle(crys_uid).get_type().type_id()
        unit_type = unit_cell.particle_type(unit_uid)
        if np_species_map[crys_type] != unit_type:
            return False

    uids = np.array(unit_nodes)
    types = np.array([unit_cell.cubeAtUID(uid).get_cube_type().type_id() for uid in uids])
    return best_rot_idx, unit_coords, uids, types


@dataclasses.dataclass
class CrystalloydRating:
    pc: PolycubeStructure # when this is instantiated we are not confident that this is a crystalloyd
    # next two properties will be none if the polycube has not been successfully mapped to a crystalloyd
    # unit cell
    unit_cell: Union[TypedStructure, None] = None
    mapping: Union[dict[int, int], None] = None
    topology_idx: int = None  # index of the topology, as a multifarious structure

    score: float = -1  # property for classing crystalloyd score. -1 = not scored


def match_crystalloyd_to_unit_cell(crystalloyd: PolycubeStructure,
                                   unit_cell: TypedStructure,
                                   np_species_map: dict[int: int],
                                   allow_small: bool = True) -> Generator[tuple[int, np.ndarray, np.ndarray, np.ndarray], None, None]:
    """
    matches the provided crystalloyd to the unit cell
    Parameters:
        crystalloyd a crystalloyd polycube strucutre
        unit_cell unit cell to match to the crystalloyd
        np_species_map a mapping of cube type IDs in crystalloyd to particle types in unit_cell
    Returns: a tuple (rotation idx, translation vector)
    """
    # if the crystalloyd candidate we're testing is smaller than the unit cell,

    if len(unit_cell) > len(crystalloyd):
        # we can override to force ignore small crystalloyd candidates. otherwise:
        if allow_small:
            # find isomorphic subgraph of unit_cell.graph and crystalloyd.graph (both are nx.MultiDiGraph objects)

            g_unit = unit_cell.to_igraph()
            g_crys = crystalloyd.to_igraph()

            matcher = g_unit.get_subisomorphisms_vf2(
                other=g_crys,
                color1=g_unit.vs['type'],
                color2=[np_species_map[t] for t in g_crys.vs['type']]
            )

            for sub_iso in matcher:
                match = match_unit_cell_from_igraph_mapping(crystalloyd, unit_cell, sub_iso, np_species_map)
                if match:
                    yield match

        else:
            return
    else:
        # step 1: find point in self to start looking for core
        crystalloyd_centroid: np.ndarray = crystalloyd.centroid()
        # step 2: find box size of unit cell. this is harder than expected
        tiling_structure_dimensions = unit_cell.dimension_sizes()

        # step 3: construct coordinate matrix + particle unique ids + particle types
        coord_matrix, uids, ptypes = unit_cell.matrix_etc()

        # we loop through rotations, which are defined as mappings [6]->[6] (dir idxs mapped to other dir idxs)
        # Bounding box check for valid centerpoints
        crystalloyd_min = crystalloyd.matrix().min(axis=0)
        crystalloyd_max = crystalloyd.matrix().max(axis=0)

        # helper function,

        translation_options = list(itertools.product(*[list(range(int(i), int(j))) for i, j
                                                       in zip(crystalloyd_min, crystalloyd_max)]))
        # sort translation options by distnace from crystalloyd centroid to point
        translation_options.sort(key=lambda pt: sum((crystalloyd_centroid - pt) ** 2))
        # try them until you find a good one. list will often be v long but we don't generally
        # need to go that long
        for tran in translation_options:
            # search rotations
            for rot_idx, rmap in enumerateRotations().items():
                # gonna lose it
                match = match_unit_cell_at_pt(crystalloyd, tiling_structure_dimensions, coord_matrix, np_species_map,
                                              ptypes, uids, rmap, rot_idx, tran)
                if match:
                    yield match


def gen_map_crystalloyd_onto(crystalloyd: PolycubeStructure,
                             unit_cell: TypedStructure,
                             species_np_map: dict[int, int]) -> Union[CrystalloydRating, bool]:
    """
    generates a mapping that maps each cube in self onto the unit cell of the structure
    returns false if no such mapping is possibles
    """

    # step 1: find a single copy of the unit cell that can map onto the crystalloyd
    for rot_idx, coord_matrix, uids, ptypes in match_crystalloyd_to_unit_cell(crystalloyd, unit_cell, species_np_map):
        structure_mapping = try_crystalloyd_map_rotran(crystalloyd,
                                                       unit_cell,
                                                       species_np_map,
                                                       rot_idx,
                                                       coord_matrix,
                                                       uids,
                                                       )
        if structure_mapping:
            return CrystalloydRating(crystalloyd, unit_cell, structure_mapping)
    return False


def try_crystalloyd_map_rotran(crystalloyd: PolycubeStructure,
                               unit_cell: TypedStructure,
                               species_np_map: dict[int:int],
                               rot_idx: int,
                               coord_matrix: np.ndarray,
                               uids: np.ndarray) -> Union[bool, dict[int, int]]:
    assert len(uids) > 0
    rot_map = enumerateRotations()[rot_idx]
    # we cant use the "StructuralHomomorphism" class here bc our homomorphism isn't injective
    # may be worth object-orientifying anyway tho
    # first populate mapping w/ unit cell
    structure_mapping: dict[int, int] = {
        crystalloyd.cubeAtPosition(pos).get_uid(): uid for pos, uid in zip(coord_matrix, uids)
    }
    # generate queue of cubes that have at least one adjacent position with a missing cube
    edge_cubes_queue: list[PolycubesStructureCube] = [
        crystalloyd.cubeAtPosition(pos) for pos in coord_matrix if
        any([crystalloyd.has_cube_at(pos + d) and crystalloyd.cubeAtPosition(pos + d).get_uid() not in structure_mapping
             for d in RULE_ORDER])
    ]
    # iterative mapping algorithm
    # completely by accident i have made this a breadth first "search"
    while len(structure_mapping) < crystalloyd.num_particles():
        # this algoritm is going to suck, theoretical max O(inf)
        # otoh premature optimization is the root of all evil
        # start by picking a random cube that we have already mapped
        cube = edge_cubes_queue[0]
        edge_cubes_queue = edge_cubes_queue[1:]
        # for (u, v, k) in unit_cell.graph.out_edges(structure_mapping[cube.get_uid()]):
        for (u, v, k) in unit_cell.graph.edges:
            # we are only looping edges originating at cube, so if u is not cube skip
            if u != structure_mapping[cube.get_uid()]:
                continue  # there has to be a better way t
            # get direction index of (u,v,k), post rot mapping
            dirIdx = rot_map[unit_cell.graph.edges[u, v, k]["dirIdx"]]
            # if our cube has a neighbor  in diection,
            # todo: edge between?? we can maybe use polycube structure graph here??
            if crystalloyd.has_cube_at(cube.position() + RULE_ORDER[dirIdx]):
                neighbor: PolycubesStructureCube = crystalloyd.cubeAtPosition(cube.position() + RULE_ORDER[dirIdx])
                # if the neighboring cube isn't in our structure mapping yet
                if neighbor.get_uid() not in structure_mapping:
                    # if the types don't match up, that's enough to ruin everything
                    if species_np_map[neighbor.get_type()] != unit_cell.particle_type(v):
                        # break out of all loops up to top
                        return False

                    # add neighbor to structure mapping
                    structure_mapping[neighbor.get_uid()] = v
                    # if neighbor has any cubes not in our queue yet
                    if any([
                        crystalloyd.has_cube_at(neighbor.position() + d)
                        and crystalloyd.cubeAtPosition(neighbor.position() + d) not in edge_cubes_queue
                        for d in RULE_ORDER
                    ]):
                        edge_cubes_queue.append(neighbor)
                    else:
                        pass  # only here so i can be sure this conditional works
    return structure_mapping


class TypedTopology(TypedStructure):
    special_types: dict[int, int]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if "nanoparticles" in kwargs and len(kwargs["nanoparticles"]) > 0:
            self.special_types = {int(key): val for key, val in kwargs["nanoparticles"].items()}
        else:
            self.special_types = dict()

    def particle_type(self, particle_id: int) -> int:
        return self.special_types[particle_id] if particle_id in self.special_types else -1

    def get_particle_types(self) -> dict[int, int]:
        return self.special_types

    def get_topology(self) -> libtlm.TLMTopology:
        """
        converts to C++ library object
        """
        return libtlm.TLMTopology({
            (
                (i, to_diridx(di)),
                (j, to_diridx(dj))
            ) for i, di, j, dj in self.bindings_list},
            self.get_particle_types())

