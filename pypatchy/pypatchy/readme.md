PyPatchy Top Level Source Module
================================

# [analysis](analysis)
Module for analyzing patchy particle simulation data.

# [patchy](patchy)
Primary submodule. Used for setting up and managing patchy particle simulations.

# [polycubeutil](polycubeutil)
Module for working with structures built out of multiple instances of cube types, on a 3D integer lattice.

(Note: The following descriptions were produced with the aid of Claude AI)

# cell_lists.py

This file implements **spatial partitioning using cell lists** for efficient neighbor searching in molecular simulations.

## Classes:

### 1. **`Cell`** (dataclass)
Represents a single cubic cell in the spatial grid

**Attributes:**
- `idxs`: NumPy array of cell indices (x, y, z) in the grid
- `startcoords`: NumPy array of cell's minimum corner coordinates
- `endcoords`: NumPy array of cell's maximum corner coordinates  
- `particles`: List of particles contained in this cell

**Methods:**
- `__post_init__()`: Validates that endcoords >= startcoords

### 2. **`CellLists`**
Manages a 3D grid of cells for spatial partitioning of a simulation box

**Attributes:**
- `cells_shape`: Tuple (nx, ny, nz) number of cells in each dimension
- `cells`: Dict mapping (x,y,z) indices to Cell objects
- `particle_cells`: Dict mapping particle UIDs to their containing Cell
- `cell_size`: Size of cubic cells
- `cells_matrix_dimensions`: NumPy array [nx, ny, nz]
- `_box_size`: NumPy array of simulation box dimensions

**Methods:**
- `__init__()`: Initializes with None/0 values
- `get_cell()`: Returns Cell for a given position (ndarray), particle UID (int), particle object, or cell index tuple
- `interaction_cells()`: Generator yielding the 27 neighboring cells (3×3×3 including self) with periodic boundaries
- `interaction_particles()`: Generator yielding all particles in neighboring cells that could interact with given particle
- `box_size()` / `set_box_size()`: Get/set simulation box size (triggers cell recomputation if size changes)
- `compute_cell_size()`: Computes cell size from either explicit value, number of cells, or number of particles
- `get_cell_size()`: Returns current cell size
- `apportion_cells()`: Divides simulation box into grid of cells based on cell_size
- `apportion_cell_particles()`: Assigns particles to cells based on their positions

## Notes:

- Uses periodic boundary conditions (modulo operation in `interaction_cells()`)
- Optimizes particle-particle interaction searches from O(n²) to O(n) by only checking nearby cells
- Cell size can be computed automatically based on particle count using heuristic: n_cells ≈ (n_particles/2)^(1/3)/2

# configure.py

This file handles **initial setup and configuration** for the pypatchy package by creating config files and directory structure in the user's home directory.

## Functions:

### **`create_default_config()`**
Sets up the pypatchy environment in `~/.pypatchy/`

**What it does:**

1. **Creates config file** at `~/.pypatchy/settings.cfg` with default values:
   - `ANALYSIS` section: `simulation_data_dir = "/scratch/jrevan21/patchysimulations"`
   - `SETUP` section: `server_config = "agave_classic"`
   - Only creates if file doesn't already exist (won't overwrite)

2. **Copies spec files**: 
   - Recursively copies subdirectories from package's `spec_files/` to `~/.pypatchy/spec_files/`
   - Excludes `test_files` subdirectory
   - Uses `importlib.resources` to access package resources
   - Handles `FileNotFoundError` if resource directory doesn't exist

3. **Creates directory structure**:
   - `~/.pypatchy/input/targets/`
   - `~/.pypatchy/output/logs/`

**Script execution:**
- Can be run directly: `python configure.py`
- Executes `create_default_config()` when run as main

# interaction_matrix.py

This file defines an **interaction matrix data structure** for managing pairwise patch color interactions in patchy particle simulations.

## Module-level functions:

### **`make_int_mat_for_colors(maxColor: int)`**
Factory function that creates an interaction matrix where each color c (1 to maxColor) interacts only with its negative -c, all with strength 1.0

## Classes:

### **`InteractionMatrix`**
Wraps a dictionary mapping color pairs to interaction strengths, with helper methods for querying and manipulating interactions

**Attributes:**
- `__interactions`: Private dict mapping `(color1, color2)` tuples to float strengths
  - Always stores with smaller color first: `(c1, c2)` where `c2 > c1`

**Methods:**

- `__init__()`: Accepts dict or list of `((c1, c2), strength)` tuples, normalizes to store with smaller index first

- `intmap()`: Returns bidirectional dictionary with both `(c1, c2)` and `(c2, c1)` entries

- `__getitem__()`: Gets interaction strength for a color pair, returns 0.0 if not defined

- `__setitem__()`: Sets interaction strength for a color pair (automatically normalizes order)

- `__contains__()`: 
  - For tuple `(c1, c2)`: checks if interaction exists
  - For int `c`: checks if color interacts with anything

- `num_interactions()`: Returns number of unique interactions defined

- `is_one_to_one()`: Returns True if every connected component in interaction graph has ≤2 nodes

- `num_colors()`: Returns total number of distinct colors (positive + negative)

- `to_array()`: Converts to 2D NumPy array where indices are colors and values are strengths

- `self_interacts()`: Checks if a color can interact with itself

- `__iter__()`: Iterates over `((c1, c2), strength)` items

- `get_interacting_colors()`: Generator yielding all colors that interact with given color c

- `graph()`: Returns NetworkX graph where nodes are colors and edges represent interactions with strength attributes

## Notes:
- Symmetric storage ensures `(c1, c2)` and `(c2, c1)` refer to same interaction

# structure.py

This file defines **graph-based representations of molecular structures** for patchy particle assemblies, particularly for polycube/lattice-based systems.

## Key type definitions:

- `Binding`: Type alias for `tuple[int, int, int, int]` representing (node1, direction1, node2, direction2)

## Utility functions:

- `flip()`: Reverses a binding tuple to swap node perspective
- `get_nodes_overlap()`: Finds nodes common to all provided cycles

## Main classes:
# structure.py - Structure class

This class defines a **graph-based representation of molecular structures** for polycube/lattice-based patchy particle systems.

## Class:

### **`Structure`**
Represents molecular structures using NetworkX MultiDiGraph where nodes are particles and edges are directional bonds

**Attributes:**
- `graph`: NetworkX MultiDiGraph storing the structure topology

**Key properties:**
- `bindings_list`: Cached property returning set of unique bindings as `(u, du, v, dv)` tuples
  - Each binding appears once despite bidirectional edge storage in graph
  - Computed from graph edges, validates that edge count = 2 × binding count

**Initialization:**
- `__init__()`: Accepts either:
  - `bindings`: List of `(n1, d1, n2, d2)` tuples (creates bidirectional edges)
  - `graph`: Existing NetworkX MultiDiGraph
  - Empty initialization if neither provided

**Core methods:**
- `to_json()`: Serializes bindings_list to nested list format
- `vertices()`: Returns list of node IDs as integers
- `filter_edges()`: Returns edges satisfying a predicate function
- `rotate()`: Applies rotation mapping to structure, returns StructuralHomomorphism
- `num_dimensions()`: Returns 3 (hardcoded, has TODO for dynamic)
- `dimension_sizes()`: Computes max length along each dimension for periodic structures
- `get_empties()`: Returns list of `(vertex_id, direction)` tuples for unbound patches
- `draw()`: Visualizes structure using matplotlib with NetworkX layout
- `substructures()`: Generator yielding all connected substructures (2 to n-1 nodes)
- `substructure()`: Extracts connected substructure containing specified nodes, remaps IDs starting from 0
- `homomorphism()` / `homomorphisms()`: Finds structural homomorphisms to another structure using igraph's VF2 algorithm with rotation enumeration
- `edge_exists()` / `bi_edge_exists()`: Check if directional edge exists from vertex
- `positions_connected()`: Tests if edge exists between two vertices
- `is_connected()` / `is_multifarious()`: Check graph connectivity
- `iter_components()`: Generator yielding connected components as Structure objects
- `num_vertices()`: Returns node count
- `splice()`: Combines two structures by connecting at specified bindings, remapping vertex IDs to avoid conflicts
- `slice()`: Returns structure containing only edges orthogonal to specified axis
- `remap_vert_ids()`: Creates copy with vertex IDs remapped to avoid collisions, returns StructuralHomomorphism
- `matrix()`: Returns N×3 NumPy array of particle coordinates
- `neighbor()`: Returns neighbor vertex in given direction (or None)
- `is_crystal()`: Checks if structure has cycles along single axes (indicates crystallinity)
- `to_igraph()`: Converts to igraph format with validation checks
- `is_bindings_euclidian()`: Checks if all bindings connect opposite directions

**Container methods:**
- `__contains__()`: Check if vertex ID or binding exists
- `__len__()`: Returns number of vertices
- `__str__()`: Returns formatted string with vertex and binding counts

## Notes:

- Contains several commented-out lines that can be deleted (lines 13, 24, 70, 181)
- Docstring references external documentation URL
- Many methods note ChatGPT assistance in their implementation

### 2. **`TypedStructure`** (ABC, extends `Structure`)
Abstract class for structures where particles have types

**Abstract methods:**
- `particle_type()`: Returns type ID for a particle
- `set_particle_types()`: Updates particle type assignments

### 3. **`GenericTypedStructure`** (extends `TypedStructure`)
Concrete implementation with dictionary-based particle typing

**Key attributes:**
- `particle_types`: Dict mapping vertex IDs to particle type IDs

**Key methods:**
- `particle_type()`: Returns type for a given particle ID
- `set_particle_types()`: Updates particle type mappings
- `remap_vert_ids()`: Reassigns vertex IDs to avoid collisions, returns homomorphism
- `get_particle_types()`: Returns the particle types dictionary
- `num_particle_types()`: Counts distinct particle types
- `substructures()` / `substructure()`: Return typed substructures
- `iter_components()`: Generator yielding connected components as typed structures
- `splice()`: Merges structures while preserving particle types
- `rotate()`: Abstract method (pass statement only)
- `slice()`: Extracts edges orthogonal to specified axis
- `to_igraph()`: Converts to igraph format with type vertex attributes

### 4. **`StructuralHomomorphism`**
Represents a mapping between two structures, preserving bindings

**Attributes:**
- `source`: Source Structure
- `target`: Target Structure
- `mapping`: Dict mapping source vertex IDs to target vertex IDs
- 

## Module-level functions:

- `read_topology()`: Loads structure from JSON file with bindings and particle types
- `calcEmptyFromTop()`: Computes empty patches from binding list

## Design notes:

- Structures use bidirectional edge representation (each bond stored as two directed edges)
- Direction indices (0-5) represent patches on cubic lattice faces
- Heavy use of NetworkX for graph algorithms and igraph for advanced analysis

# origami_deabstractor.py

This file defines an **abstract base class for converting coarse-grained patchy particle simulations into detailed DNA origami structures**.

## Class:

### **`OrigamiDeabstractor`** (ABC)
Handles conversion from abstract patchy particle representations to nucleotide-level DNA structures

**Class attributes:**
- `scene`: Scene object containing patchy particles
- `check_strict_rc`: Boolean for strict reverse complement checking (default True)
- `particle_type_map`: Dict mapping patchy particle type IDs to DNA particle templates
- `sticky_book`: Set tracking which sticky ends have been added (or None if tracking disabled)
- `dna_particles`: Dict mapping particle UIDs to positioned DNA particle instances
- `sticky_length`: Length of sticky end sequences (can be None for dynamic calculation)
- `color_sequences`: Dict mapping colors (int or str) to DNA sequences
- `bondcount`: Counter for bonds
- `scale_factor`: Converts simulation distance units to DNA model units (can be None)

**Methods:**

- `__init__()`: Initialize with optional sticky_length, spacer_length, and expected_num_edges

- `set_track_conf_stickies()` / `is_track_conf_stickies()`: Enable/disable and check sticky end tracking

- `get_dna_origami()`: Get DNA particle template for a patchy particle type (accepts type ID, name, or particle object)

- `get_full_conf()`: Merges all DNA particles into single DNAStructure, optionally clearing clusters

- `get_sticky_length()`: Returns sticky length (hardcoded value or average of assigned sequences)

- `color_sequence()`: Gets/generates DNA sequence for a color
  - If color exists: returns stored sequence
  - If matching color exists: returns reverse complement
  - Otherwise: generates random sequence

- `get_color_match()`: Returns complementary color
  - For strings: "dark" prefix toggle (e.g., "darkred" ↔ "red")
  - For ints: negation (e.g., 5 ↔ -5)
  - Respects `color_match_overrides` dict

- `ready_to_position()`: Checks if all particle types have DNA mappings

- `position_particles()`: Clones DNA templates, aligns with patchy particles, applies scale factor and positions

- `get_particles()`: Returns list of positioned DNA particles (calls `position_particles()` if needed)

- `get_dna_particle()`: Gets positioned DNA particle for a specific patchy particle instance

- `dump_monomer()`: Exports single DNA particle type to .top and .dat files with sticky ends added

- `dump_monomers()`: Exports all particle types to separate files, optionally filtering to only types present in scene

- `convert()`: Abstract method for subclasses to implement conversion logic

- `assign_particles()`: Assigns DNA particle template(s) to patchy particle type(s), calls `link_patchy_particle()`

# patchy_base_particle.py

This file defines **abstract base classes for patchy particle systems**, establishing the core hierarchy for particle types, patch types, particle sets, and particle instances.

## Classes:

### 1. **`PatchyBaseParticleType`** (ABC)
Abstract base class for particle type definitions (templates)

**Attributes:**
- `_type_id`: Unique identifier for particle type
- `_patches`: List of patches on this particle type

**Methods:**
- `__init__()`: Constructs particle type with uid and patches list
- `type_id()` / `set_type_id()`: Get/set type ID
- `name()`: Abstract method returning particle name
- `num_patches()`: Returns number of patches
- `patches()`: Returns list of patches
- `patch()`: Get patch by int index or by position (np.ndarray)
- `radius()`: Abstract method for particle radius along a normal vector

### 2. **`BasePatchType`** (ABC)
Abstract base class for patch definitions

**Attributes:**
- `_uid`: Unique patch identifier
- `_color`: Patch color/type (flexible type for forward compatibility)
- `_key_points`: List of important 3D points (position, orientation vectors, etc.)

**Methods:**
- `__init__()`: Initialize with uid and color
- `get_id()` / `set_id()`: Get/set patch ID
- `add_key_point()`: Append a key point, returns its index
- `get_key_point()` / `num_key_points()`: Access key points
- `position()`: Abstract method returning patch position
- `color()` / `set_color()`: Get/set color
- `colornum()`: Abstract method returning numeric color
- `can_bind()`: Abstract method checking binding compatibility
- `rotate()`: Applies rotation matrix to all key points in-place
- `has_torsion()`: Abstract method for torsional constraints

### 3. **`BaseParticleSet`**
Container managing collections of particle types and patch types (not abstract)

**Attributes:**
- `_particle_types`: List of particle types
- `_patch_types`: List of patch types

**Methods:**
- `__init__()`: Optionally initialize with particle list
- `particles()` / `particle()`: Access all particles or by index/type_id
- `num_particle_types()`: Count of particle types
- `patches()` / `patch()`: Access patch types
- `num_patches()`: Count of patch types
- `add_particle()` / `add_particles()`: Add particle types (auto-assigns IDs)
- `add_patch()` / `add_patches()`: Add patch types
- `__len__` / `__iter__` / `__getitem__`: Container interface

### 4. **`PatchyBaseParticle`** (ABC)
Abstract base class for particle instances in simulation

**Attributes:**
- `_uid`: Unique particle instance ID
- `_type_id`: Reference to particle type
- `_position`: 3D position as NumPy array

**Methods:**
- `__init__()`: Initialize with uid, type_id, and position
- `get_type()` / `set_type()`: Get/set type ID
- `get_uid()` / `set_uid()`: Get/set unique ID
- `rotation()`: Abstract method for particle orientation
- `position()` / `set_position()`: Get/set position
- `patches()` / `patch()` / `num_patches()`: Abstract methods for patch access
- `translate()`: Move particle by translation vector
- `rotate()`: Abstract method for rotation
- `rotation_from_to()`: Computes rotation matrix to align this particle with another using SVD superposition, trying all permutations to find best match

## Notes:
- Design uses type-instance pattern: types are templates, instances are actual particles in simulation
- Heavy use of abstract methods allows different simulation backends to implement specific behaviors

# scene.py

This file defines an **abstract base class for simulation scenes** that manage collections of particle instances and their interactions.

## Class:

### **`Scene`** (ABC)
Abstract container for particle instances in a simulation

**Attributes:**
- `_particles`: List of all particle instances in the scene

**Concrete methods:**

- `__init__()`: Initializes empty particle list

- `particles()`: Returns list of all particles
- `num_particles()`: Returns particle count
- `particle_type_counts()`: Returns dict mapping particle type IDs to their counts in scene

- `add_particle()` / `add_particles()`: Add single particle or iterable of particles to scene

- `get_particle()`: Get particle by index, with assertion that particle UID matches index

- `iter_bound_particles()`: Generator yielding all pairs of bound particles (uses `particles_bound()`)

- `iter_binding_patches()`: Generator yielding bound patch pairs between two specific particles (uses `patches_bound()`)

- `centroid()`: Computes average position of all particles

- `minmaxs()`: Returns array of (min, max) coordinate pairs for each axis (bounding box)

- `particle_coords()`: Returns stacked NumPy array of all particle positions

**Abstract methods (must be implemented by subclasses):**

- `num_particle_types()`: Returns number of distinct particle types
- `particle_types()`: Returns BaseParticleSet containing type definitions
- `set_particle_types()`: Sets the particle type definitions
- `get_conf()`: Returns Configuration in oxDNA_analysis_tools format
- `particles_bound()`: Checks if two particles are bound
- `patches_bound()`: Checks if two specific patches on two particles are bound

## Design notes:

- Provides interface for different simulation backends to implement
- Manages particle instances while delegating type management to subclasses
- Bonding/interaction queries are abstract to allow backend-specific implementations

# server_config.py

This file defines **server/cluster configuration management** for running oxDNA simulations on different computational environments.

## Module-level functions:

### **`list_server_configs()`**
Generator that iterates through all JSON files in the server configs directory and yields loaded `PatchyServerConfig` objects

### **`load_server_settings(settings_name: str)`**
Loads a server configuration from JSON file by name and returns `PatchyServerConfig` instance

### **`get_server_config()`**
Returns the currently active server configuration from user settings

## Class:

### **`PatchyServerConfig`** (dataclass)
Encapsulates all settings for running simulations on a specific server or cluster

**Attributes:**
- `name`: Configuration identifier
- `oxdna_path`: Path to oxDNA executable/installation
- `patchy_format`: Format specification for patchy particle interactions
- `slurm_bash_flags`: Dict of SLURM scheduler flags (e.g., memory, time, partition)
- `slurm_includes`: List of bash commands to include in scripts (e.g., "module load cuda")
- `input_file_params`: Default oxDNA input file parameters (ParamSet)
- `absolute_paths`: Whether to use absolute paths in job scripts (default False)
- `is_slurm`: Whether this is a SLURM-managed cluster (default False)
- `cuda_mps`: CUDA Multi-Process Service for GPU sharing (bool or SimulationManager, default False)
- `no_oxpy_check`: Skip oxpy validation before running (default False)

**Methods:**

- `__post_init__()`: Converts input_file_params to ParamSet if needed, initializes SimulationManager if cuda_mps enabled

- `write_sbatch_params()`: Writes SLURM batch script header with:
  - Shebang (`#!/bin/bash`)
  - SBATCH flags (long flags with `--`, short flags with `-`)
  - Job name and output file path
  - Include lines (module loads, etc.)

- `is_batched()`: Returns whether simulations use batching (checks cuda_mps)

- `is_server_slurm()`: Returns whether server uses SLURM scheduler

- `get_slurm_bash_flags()`: Returns SLURM flags dict (asserts if not SLURM server)

- `get_slurm_n_tasks()`: Extracts number of tasks from SLURM flags (checks "n" or "ntasks", defaults to 1)

- `to_dict()`: Serializes configuration to dictionary (converts bools explicitly)

# slurmlog.py

This file defines a **container and query interface for SLURM job log entries** with efficient searching and filtering capabilities.

## Class:

### **`SlurmLog`**
Wrapper for a sorted list of `SlurmLogEntry` objects with accessor and filter methods

**Attributes:**
- `log_list`: List of `SlurmLogEntry` objects, sorted chronologically by submission date
- `id_map`: Dict mapping job IDs to log entries for O(1) lookup

**Initialization:**
- `__init__()`: Takes variable number of `SlurmLogEntry` objects, automatically sorts by submission date and builds id_map

**Indexing methods:**

- `__getitem__()`: Supports multiple access patterns:
  - `int`: Returns entry at that index
  - `datetime.date`: Returns all entries from that date
  - `slice`: Supports slicing by int or date ranges

**Binary search methods:**

- `idx_begin()`: Finds index of first entry matching given date using modified binary search, returns None if out of bounds

- `idx_end()`: Finds index of last entry matching given date using modified binary search, returns None if out of bounds

- `find_idx()`: Standard binary search to find insertion point for a datetime (note: written by ChatGPT per comment)

**Lookup methods:**

- `by_id()`: O(1) lookup by SLURM job ID using id_map

**Filter methods (all return new `SlurmLog` objects):**

- `by_type()`: Filter by job type string or list of type strings

- `by_subject()`: Filter by LogEntryObject (simulation) that job relates to

- `by_other()`: Generic filter by arbitrary key-value pairs in additional_metadata

**List operations:**

- `append()`: Add new entry maintaining chronological sort (uses binary search for insertion if needed)

- `to_list()`: Convert all entries to list of dicts

- `empty()`: Check if log is empty

- `__len__()`: Returns number of entries

- `__iter__()`: Iterate over entries

- `__repr__()`: String representation with double newlines between entries

# slurm_log_entry.py

This file defines **individual SLURM job log entries** and the interface for objects that can be logged.

## Classes:

### 1. **`LogEntryObject`** (ABC)
Abstract base class for objects that can be associated with log entries (e.g., simulations)

**Abstract method:**
- `to_dict()`: Must return dict with str, int, or float values for serialization

### 2. **`SlurmLogEntry`**
Represents a single SLURM job submission with metadata

**Attributes:**
- `job_submit_date`: DateTime when job was submitted
- `job_id`: SLURM job ID (integer)
- `simulation`: LogEntryObject reference (simulation/analysis this job relates to)
- `script_path`: Path to bash/SLURM script file
- `log_path`: Path to SLURM output log file
- `notes`: Optional user notes/comments
- `job_type`: String categorizing the job (e.g., "oxdna", "analysis")
- `additional_metadata`: Dict for arbitrary extra metadata

**Methods:**

- `__init__()`: Constructor with flexible date handling
  - Accepts `datetime` object or string in "YYYY-MM-DD" format for `start_date`
  - Converts string paths to Path objects
  - Defaults: `notes=""`, `additional_metadata={}`, `start_date=datetime.now()`

- `get_log_txt()`: Reads and returns contents of log file at `log_path`, returns error message if file not found

- `to_dict()`: Serializes entry to dict
  - Converts paths to strings
  - Formats date as "YYYY-MM-DD"
  - Calls `simulation.to_dict()` for nested serialization

- `__str__()`: Returns formatted multi-line string representation with all fields indented and formatted

## Design notes:

- Flexible design allows arbitrary metadata storage while maintaining structure
- Bridges SLURM job management with simulation tracking
- Date handling accommodates both programmatic and user input


# structure.py

This file defines **graph-based representations of molecular structures** for patchy particle assemblies, particularly for polycube/lattice-based systems.

## Key type definitions:

- `Binding`: Type alias for `tuple[int, int, int, int]` representing (node1, direction1, node2, direction2)

## Utility functions:

- `flip()`: Reverses a binding tuple to swap node perspective
- `get_nodes_overlap()`: Finds nodes common to all provided cycles

## Main classes:

### 1. **`Structure`**
Base class for molecular structure graphs using NetworkX MultiDiGraph

**Key attributes:**
- `graph`: NetworkX MultiDiGraph where nodes are particles and edges are bonds with direction indices

**Key methods:**
- `bindings_list`: Cached property returning set of unique bindings (each edge appears once despite bidirectional storage)
- `vertices()`: Returns list of node IDs
- `filter_edges()`: Returns edges satisfying a predicate function
- `rotate()`: Applies rotation mapping to structure, returns `StructuralHomomorphism`
- `num_dimensions()`: Returns 3 (hardcoded for now)
- `dimension_sizes()`: Computes max length along each dimension for periodic structures
- `get_empties()`: Returns list of unbound patches (node, direction pairs)
- `draw()`: Visualizes structure using matplotlib
- `to_json()`: Serializes bindings list to nested list format

(File truncated at line 208 - many additional methods exist but weren't fully visible)

### 2. **`TypedStructure`** (ABC, extends `Structure`)
Abstract class for structures where particles have types

**Abstract methods:**
- `particle_type()`: Returns type ID for a particle
- `set_particle_types()`: Updates particle type assignments

### 3. **`GenericTypedStructure`** (extends `TypedStructure`)
Concrete implementation with dictionary-based particle typing

**Key attributes:**
- `particle_types`: Dict mapping vertex IDs to particle type IDs

**Key methods:**
- `particle_type()`: Returns type for a given particle ID
- `set_particle_types()`: Updates particle type mappings
- `remap_vert_ids()`: Reassigns vertex IDs to avoid collisions, returns homomorphism
- `get_particle_types()`: Returns the particle types dictionary
- `num_particle_types()`: Counts distinct particle types
- `substructures()` / `substructure()`: Return typed substructures
- `iter_components()`: Generator yielding connected components as typed structures
- `splice()`: Merges structures while preserving particle types
- `rotate()`: Abstract method (pass statement only)
- `slice()`: Extracts edges orthogonal to specified axis
- `to_igraph()`: Converts to igraph format with type vertex attributes

## Module-level functions:

- `read_topology()`: Loads structure from JSON file with bindings and particle types
- `calcEmptyFromTop()`: Computes empty patches from binding list

## Design notes:

- Structures use bidirectional edge representation (each bond stored as two directed edges)
- Direction indices (0-5) represent patches on cubic lattice faces
- Heavy use of NetworkX for graph algorithms and igraph for advanced analysis