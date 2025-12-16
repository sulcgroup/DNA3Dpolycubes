This package contains utility functions for working with polycubes in Python. It provides tools for creating, manipulating, and analyzing polycubes, which are three-dimensional shapes formed by joining cubes face-to-face.

# cube.py
Defines `PolycubesStructureCube`, the runtime representation of an individual cube in a polycube structure. Each cube stores its lattice position, orientation (as a `Rotation` object), and its cube-type definition from a `PolycubeRuleCubeType`. The class provides:

- flexible constructors that accept raw values, dicts (JSON), quaternions, rotation matrices, or rotation indices  
- helpers for converting between global and local patch directions (`typedir`), checking for patches, and retrieving rotated `PolycubesPatch` objects  
- state tracking for rule-driven activation/interaction logic  
- export to a JSON-friendly dict format for serialization  

This class is the building block used by higher-level polycube structure and conversion tools. 

# polycubesRule.py
Defines the full data model for **polycube interaction rules**, including cube types, patches, dynamic state variables, and logical activation effects. It provides:

- low-level primitives such as direction indices, rotation utilities, and orientation helpers  
- the `PolycubesPatch` class for directional, optionally torsional patches with activation/state metadata  
- the `PolycubeRuleCubeType` class, representing a cube species with patches, state variables, and rule-based effects that determine which patches are active  
- the `PolycubesRule` class, a container for multiple cube types with flexible constructors (string-based rule syntax, JSON dictionaries, or empty templates)

The module supports patch-condition evaluation, state-transition graph construction, canonicalization (for comparing rules up to symmetries and color relabeling), SAT-encoding export, and various editing utilities. It acts as the foundational rule language for designing polycube-based self-assembly systems.

# polycube_structure.py
Implements `PolycubeStructure`, a full spatial and topological representation of a polycube assembly built from `PolycubesRule` cube types. It loads structures from JSON or cube objects, positions and orients cubes on an integer lattice, and constructs a directed binding graph based on patch complementarity, alignment, and cube state variables.

The class provides:

- fast neighbor detection via a spatial hash, enabling automatic bond inference  
- conversion utilities (to JSON, to typed structures, to substructures) and geometric transforms  
- cluster, cycle, and type-count queries  
- graph visualization helpers and compatibility with the broader Patchy scene/structure interface

Together, this module forms the concrete structural layer enabling rule-based design, validation, and analysis of polycube self-assembly systems.

# polycube_util.py
Provides standalone utility functions for constructing, comparing, and transforming polycube rules and structures without relying on the polycube C++ backend. It includes:

- `make_colors_consecutive()` for renumbering patch colors in a rule to a compact, consecutive range  
- `get_fully_addressable_rule()` for generating a polycube rule in which every cube position and every adjacency receives a unique cube type and unique patch color (useful for debugging or enumerating all possible bindings)  
- geometric and structural helpers such as `coord_equal()` for rotationally/permutationally checking equivalence of coordinate sets, and `rotation_mapping_to_matrix()` for converting discrete direction mappings into rotation matrices  
- a few convenience imports and small utilities used across polycube tooling

