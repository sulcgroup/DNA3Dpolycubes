# ipy_oxDNA
This library provides classes and functions for simulating DNA nanostructures using oxDNA. It acts as a wrapper for oxpy and includes an API for editing oxDNA structures pythonically. 

## oxdna_simulation.py
Provides a large, modular interface for building, running, and managing **oxDNA molecular dynamics simulations** in Python. At the top level is the `Simulation` class, which wraps all functionality needed to construct input files, launch simulations (locally, via subprocess, or on SLURM), manage forces and observables, and interact with oxDNA output. The module includes:

- builder components (`BuildSimulation`, `BuildSimulationFromStructure`) for assembling topology, configuration, and input files  
- runtime interfaces (`OxpyRun`, `SlurmRun`) for executing simulations via `oxpy` or SLURM scripts  
- helper subsystems for forces (`SimForces`), sequence-dependent parameterization (`SequenceDependant`), proteins, and runtime observables  
- a high-level `SimulationManager` that schedules many simulations across CPUs/GPUs with memory-aware allocation  

Overall, this module provides the operational infrastructure enabling automated oxDNA workflows inside PyPatchy, from file preparation to full simulation orchestration and analysis hooks. 

## defaults/defaults.py
Contains the logic for loading and evaluating **default oxDNA input templates** used throughout the simulation system. The core class, `DefaultInput`, loads a preset JSON file from the `ipy_oxdna` package, supports dynamic expressions of the form `f(x,y) = ...`, and evaluates them using user-provided keyword arguments. After evaluation, it exposes a clean, fully resolved dictionary suitable for writing to `input`/`input.json`.

The module also defines:

- `default_input_exist()` — checks whether a named default template exists  
- `get_default_input()` — returns a `DefaultInput` instance for a given preset  
- Predefined parameter dictionaries (`SEQ_DEP_PARAMS`, `NA_PARAMETERS`, `RNA_PARAMETERS`) used for sequence-dependent simulations  
- Exceptions (`IncompleteInputError`, `MissingParamError`) that ensure all dynamic expressions are fully evaluated before input generation  

In short, this file provides the machinery that turns high-level input presets into ready-to-use oxDNA simulation parameters. 

## defaults/inputs/
Directory containing default input file JSONs. The files are sets of common oxDNA operations like relaxation and production simulations, with both CPU and GPU variants. These files are loaded and evaluated by the `DefaultInput` class in `defaults.py`.

## dna_structure/dna_structure.py
Defines a comprehensive toolkit for representing, manipulating, and exporting **DNA structures** compatible with oxDNA and OxView. The core class, `DNAStructure`, holds multiple `DNAStructureStrand` objects, each storing base identities, positions, orientations, and globally unique base IDs. The module supports:

- sequence editing, strand splicing (`nick`, prepend/append), strand merging, and UID remapping  
- geometric transformations, helix construction (`construct_strands`, `gen_helix_coords`, `helixify`), and validation  
- IO utilities to load DNA structures from `.top`/`.dat` or `.oxview` files, and export structures back into oxDNA or glTF formats  
- cluster and coloration assignment for visual or analytical grouping  
- base-level indexing, sequence extraction, bonding checks, and mass calculations  

Together, these tools provide a high-level, editable in-memory representation of DNA suitable for building, analyzing, and visualizing complex origami and duplex assemblies. 

## utils/force.py
Implements the full set of **external forces** available to oxDNA simulations within this wrapper. It provides:

- an enumeration (`ForceType`) describing each supported force, its required parameters, and default values  
- the `Force` class, a dictionary-like container that validates parameters, supports spatial transformations, and serializes cleanly to JSON-compatible formats  
- convenience generators such as `zip_strands()` for pairing nucleotides with mutual traps  
- factory functions (`morse`, `skew_force`, `mutual_trap`, `string`, `harmonic_trap`, `rotating_harmonic_trap`, `repulsion_plane`, `repulsion_sphere`) to produce force dictionaries in legacy style  
- loaders for force specifications (`load_forces_from_txt`, `load_forces_from_json`)  

Together, this module defines the interface for specifying per-particle or multi-particle external fields used in controlled oxDNA simulations (e.g., traps, pulling forces, COM restraints, repulsive boundaries). 

## utils/observable.py
Defines a system for specifying **oxDNA observables**—quantities to be recorded during simulation runs. It replaces an older `Observable` class with a cleaner interface built around two components:

- **`Observable`**: a container specifying the output file name, print interval, and one or more observable columns. It exports to the JSON format expected by oxDNA.
- **`ObservableColumn`**: a single measurement type (e.g., `distance`, `hb_list`, `potential_energy`) with its corresponding parameters.

The module also provides convenience functions (`distance`, `hb_list`, `particle_position`, `potential_energy`, `force_energy`, `kinetic_energy`, `pair_energy`) that construct correctly formatted observable dictionaries in one call. These observables integrate with the simulation builder to generate `observables.json` files automatically.

Overall, the module offers a uniform, extensible mechanism for defining simulation-time measurements in oxDNA workflows. 

## utils/util.py
Contains general-purpose helper functions used across the codebase. It includes:

- `rotation_matrix()` — computes a 3×3 rotation matrix for a given axis/angle  
- `process_path()` — normalizes relative or user-prefixed paths into `Path` objects  
- unit-conversion dictionaries and functions (`si_units`, `ox_units`) for translating between SI and oxDNA/oxRNA units  
- `generate_distinct_colors()` — produces visually distinct, colorblind-friendly colors for plotting  

These utilities support filesystem handling, numerical operations, visualization, and physical-unit conversions needed throughout simulation and analysis workflows. 

