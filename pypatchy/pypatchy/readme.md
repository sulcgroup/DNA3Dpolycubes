PyPatchy Top Level Source Module
================================

# [analysis](analysis)
Module for analyzing patchy particle simulation data.

# [patchy](patchy)
Primary submodule. Used for setting up and managing patchy particle simulations.

# [polycubeutil](polycubeutil)
Module for working with structures built out of multiple instances of cube types, on a 3D integer lattice.

# cell_lists.py
Implements a lightweight spatial partitioning system used to accelerate neighbor searches in patchy-particle simulations. It defines a `Cell` structure representing cubic subregions of the simulation box, and a `CellLists` class that divides the simulation volume into a regular 3D grid, assigns particles to cells based on their positions, and efficiently yields nearby cells or particles for interaction checks. The module supports periodic wrapping of neighbor queries and handles dynamic resizing when the simulation box or cell size changes.

# configure.py
 Provides a utility function for initializing a default PyPatchy configuration environment in the user’s home directory. It creates a `~/.pypatchy/settings.cfg` file with basic ANALYSIS and SETUP entries, populates the local `spec_files` directory by copying bundled specification files (excluding test data), and ensures that required subdirectories such as `input/targets` and `output/logs` exist. The script is primarily intended as a one-time setup helper for new installations and can be run directly as a standalone script.

# interaction_matrix.py
Defines a lightweight data structure for representing pairwise interaction strengths between color types used in patchy-particle and crystal-design simulations. It provides an `InteractionMatrix` class that stores symmetric interactions, supports querying and updating strengths, converts itself into array or graph form, and analyzes interaction structure (e.g., one-to-one connectivity). A small factory function builds a default matrix where each color only interacts with its negative counterpart. 

# origami_deabstractor.py
Defines an abstract base class that converts a high-level patchy-particle scene into a fully realized DNA origami representation suitable for oxDNA simulations. It manages mappings between particle types and their DNA counterparts, assigns and generates complementary sticky-end sequences, positions and orients DNA structures according to a simulation scene, and exports monomer files for downstream use. Subclasses implement a `convert()` method that drives the full translation from patchy particles to DNA structures.

# patchy_base_particle.py
Defines the abstract foundation for all patchy-particle types and particle instances used throughout the simulation framework. It includes base classes for particle *types* (which define patch geometry and metadata), patch types (which define patch positions, colors, and rotational behavior), particle *instances* (which carry UID, type, position, rotation, and patch access), and a `BaseParticleSet` container for managing collections of particle and patch types. Together these classes provide the minimal shared interface for geometric queries, rotation/alignment, patch access, and type registration that higher-level patchy and DNA-origami modules build upon.

# scene.py
Defines the abstract `Scene` class, which represents a collection of patchy particles along with the rules needed to query their configuration, types, and binding interactions. It manages particle storage, provides utilities for iterating over bound particles or bound patch pairs, exposes geometric helpers such as centroids and bounding boxes, and delegates all system-specific logic—like particle types, binding rules, and configuration export—to subclasses. This class serves as the backbone for any simulation or structural environment that contains patchy particles.

# server_config.py
Defines the configuration layer that connects PyPatchy to the computational environment where oxDNA simulations are executed. It provides tools to load server-specific JSON config files, enumerate available configurations, and represent them as a `PatchyServerConfig` dataclass containing paths, SLURM settings, input-parameter presets, and execution flags. The module standardizes how simulation jobs are prepared—especially on SLURM clusters—by generating sbatch headers, handling CUDA MPS batching, and exposing convenience accessors for job parameters. 

# slurmlog.py
Implements a convenience wrapper around lists of `SlurmLogEntry` objects, providing indexed, date-based, and metadata-based access to logged SLURM jobs. It maintains entries sorted by submission time, supports efficient binary-search lookup for date ranges, and offers filtering utilities (by job ID, type, subject, or arbitrary metadata keys). The class also handles insertion of new log entries in chronological order and can export the log as a list of dictionaries for serialization or analysis.

# slurm_log_entry.py
Defines the data model for an individual SLURM job record used in PyPatchy’s logging system. A `SlurmLogEntry` stores submission time, job ID, job type, associated simulation metadata, script and log file paths, freeform notes, and optional user-defined metadata. It supports converting entries to dictionaries, loading the contents of the corresponding log file, and producing a formatted string representation. An abstract `LogEntryObject` base class defines the required interface for simulation objects that appear within log entries.

# structure.py
Implements the core graph-based representation of polycube structures and their transformations. A `Structure` stores a directed multigraph encoding bindings between cube-like units, and provides extensive tooling for analyzing, transforming, slicing, tiling, rotating, and comparing these structures. It can compute coordinates, extract substructures, detect crystallinity, evaluate homomorphisms via igraph, and generate structural mappings.

Several subclasses extend this model:  
• `FiniteLatticeStructure` adds lattice-indexed coordinate utilities.  
• `StructuralMapping` and `StructuralHomomorphism` represent vertex/direction mappings between structures.  
• `TypedStructure` and `GenericTypedStructure` incorporate particle-type annotations and enable type-constrained homomorphisms.  

The module also provides helpers for reading topologies from JSON, computing empty patch sites, and supporting periodic tiling operations. Overall, it forms the computational backbone for reasoning about 3-D polycube assemblies, symmetry, and structural compatibility.

# util.py
Provides a collection of stateless helper functions used throughout the PyPatchy codebase, covering vector math, rotations, quaternion operations, color generation, and small filesystem/configuration utilities. It also includes functions for working with simulation directories, reading spec files, and enumerating rotational symmetries. Two small custom exceptions support missing-file handling and validating simulation directory structure. Overall, the module serves as a central toolbox for common mathematical and I/O tasks across the project.

# vis_util.py
Provides a small visualization helper for assigning consistent colors to particle types. Using a golden-angle hue progression, it maps each particle-type index to an RGB value via HSV conversion, ensuring visually distinct yet harmonized colors for plots or renderings of polycube or patchy-particle structures. The module is intentionally minimal but designed to match the color conventions used elsewhere in the codebase.


