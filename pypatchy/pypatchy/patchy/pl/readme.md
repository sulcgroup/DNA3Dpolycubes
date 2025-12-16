Code for specifc subsets of patchy particles. What exactly distinguishes a "PL" patchy particle is not well defined and none of us know what "pl" means. However the category is useful for the purpose of our code structure, especially in distinguishing between patch particles and particles on the integer lattice.

# patchy_exports.py
Provides functionality for exporting patchy-particle simulation scenes to **glTF/GLB** files for visualization in 3D tools such as Blender or web-based viewers. It includes geometry generators for spheres, cubes, cones, and skewed “torsional” cones used to visually represent patches. The main function, `to_gltf()`, builds a full glTF scene from a `PLPSimulation`, creating buffers, accessors, materials (colored per particle type or patch), meshes, and nodes corresponding to each particle and its orientation. The exporter supports optional centering, selectable particle shapes, and patch visualization modes. Overall, this module serves as the bridge between PyPatchy simulation data and portable 3D visualization formats.

# patchyio.py
Implements the I/O backbone for reading and writing patchy-particle simulations across several different file formats (Flavio [Romano], Lorenzo [Rovigatti], Josh [Evans], Subhajit [Roy], raspberry [the fruit], and MGL-like formats). It defines an abstract `PLBaseWriter` class and multiple concrete writer subclasses, each responsible for:

- loading particle types, patches, and interaction matrices  
- writing topology (`.top`), configuration (`.dat`), particles/patches specification files, and input-file directives  
- reading full simulation snapshots into `PLPSimulation` objects  
- managing absolute/relative paths, conf indexing, and format quirks

The module also includes an MGL→patchy converter, utilities for exporting interaction matrices, and a registry (`get_writer`, `register_writer`) for selecting writers based on server or user settings. Overall, it provides the serialization, deserialization, and format-interoperability layer that the rest of PyPatchy relies on for simulation setup and replay.

# pljsonrw.py
Implements JSON encoding/decoding utilities for PL-style patchy-particle data structures. It provides small helper functions to serialize and reconstruct patches, particle types, particle sets, multidentate source maps, and conversion settings. The central class, `PLJSONEncoder`, extends Python’s `json.JSONEncoder` to properly handle numpy types, PLPatch/PLPatchyParticle objects, full `PLParticleSet` instances, and multidentate-conversion metadata. This module serves as the low-level JSON I/O layer used by higher-level ensemble serialization tools. 

# plparticle.py
Defines the `PLPatchyParticle` class, which represents both **particle types** and **particle instances** within the PL patchy-particle framework. Each particle contains a list of patches, a radius, a position, and an orientation described by orthonormal vectors `a1`, `a3` (with `a2` derived by cross product). The class provides:

- helpers for managing patches, querying patch positions/orientations, and instantiating particles  
- rotation and normalization utilities for aligning particles or applying random orientations  
- export helpers for MGL/XYZ formats and compatibility with downstream simulation writers  
- geometric logic such as rotated patch placement and matching particle types for structural comparison

It acts as the core geometric and structural representation of particles used throughout the simulation and conversion pipeline. 

# plparticleset.py
Defines the `PLParticleSet` class, the core container for PL-format particle types and their patches. It manages:

- interaction matrices for patch–patch binding  
- multidentate conversion via `MultidentateConvertSettings` and `PLMultidentateSourceMap`, which map unidentate patches to multiple “teeth” and generate new patch geometries, colors, and energy scalings  
- particle lookup, normalization checks, and flexible equality tests  
- patch grouping utilities (important for multidentate behavior and geometry alignment)

Supporting classes (`PLSourceMap`, `PLMultidentateSourceMap`) store mappings between original and converted particle sets, while `decode_particles_dict()` reconstructs particle sets from JSON-like structures.

Overall, this module implements the structural and interaction-level logic required for both standard and multidentate patchy-particle models. 

# plpatch.py
Defines the `PLPatch` class, the fundamental representation of a patch on a PL-format patchy particle. Each patch stores a type ID, color, interaction strength, relative position, and orientation vectors (`a1`, `a2`, with `a3` derived). The class provides utilities for rotation, binding checks, serialization to/from text formats, and structural equality. 

The module also contains helper functions for constructing interaction matrices and for assigning colors to patches based on an underlying interaction graph, including a one-to-one bipartite coloring algorithm and a fallback assignment strategy. Collectively, these tools define how patches are represented, validated, and converted within the PL patchy-particle framework. 

# plpatchylib.py
Contains conversion utilities between various particle representations—MGL scenes, PL particle sets, and polycube structures. It provides:

- `polycube_rule_to_PL()` and `PL_to_rule()` for converting between PolycubesRule objects and PL-format particle sets  
- `polycube_to_pl()` for generating full `PLPSimulation` scenes from polycube assemblies, including multidentate expansion, orientation mapping, and box-size computation  
- `mgl_particles_to_pl()` and `mgl_to_pl()` for translating MGL particles and scenes into PL-format particles, handling color remapping, patch orientation recovery, and type alignment  
- the `MGLPLSourceMap` class for maintaining mapping information between MGL and PL particles  
- `load_pl_particles()` for reconstructing PL particle sets from serialized dictionaries or writer-based formats

Overall, this module forms the interoperability layer that allows patchy-particle workflows to bridge between MGL, PL, and polycube-based representations. 

# plpotential.py
Implements some **pairwise interaction potentials** used for PL-format patchy-particle simulations. It defines:

- a shared abstract interface (`PLPotential`, `PLPatchyPotential`) for computing particle–particle and patch–patch energies  
- **excluded-volume potentials**  
  - `PLLRExclVolPotential`: Lorenzo-style LJ-based repulsion with optional spherical attraction  
  - `PLFRExclVolPotential`: Flavio-style smoothed Lennard-Jones repulsion  
- **patchy interaction potentials**  
  - `PLLRPatchyPotential`: Lorenzo’s detailed patchy potential using an interaction matrix and exponential distance modulation  
  - `PLFRPatchyPotential`: Flavio’s non-torsional patch binding model  
  - `PLFRTorsionalPatchyPotential`: the torsional extension incorporating angular alignment terms and narrow-type presets  

The module also provides periodic-boundary distance helpers and narrow-type parameter presets. Together, these classes supply the physical interaction models governing how particles bind, repel, and align during simulations. 
Note that these classes are used to evaluate energies in simulation data produced elsewhere. They do not perform dynamics or Monte Carlo moves themselves, and should not be used to do so without rigorous additional testing.

# plscene.py
Defines `PLPSimulation`, the central in-memory representation of a patchy-particle configuration in PL format. It extends the general `Scene` class with spatial indexing via `CellLists`, physical interaction evaluation via attached `PLPotential` objects, and numerous geometric and structural utilities. Key features include:

- box management, inboxing, translation, and rotation of entire configurations  
- full and neighbor-list–accelerated energy evaluation, overlap detection, and random particle/cluster placement  
- graph construction (`compute_scene_graph`, `compute_scene_multigraph`) to identify bound particles, patch–patch bonds, and structural clusters  
- methods for extracting subscenes, merging scenes, and converting to multidentate representations  
- configuration export/import helpers for MGL and general PL workflows  

In practice, `PLPSimulation` is the workhorse object for manipulating and analyzing static or quasi-static patchy-particle states.