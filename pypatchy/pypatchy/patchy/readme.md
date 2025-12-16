# pl
Sub-package handling patchy particle simulations which utilize the `PLPatchyParticle` class. The term "PL" doesn't really come from anything. More information is available in the package readme.

# analysis_lib.py
Provides a large collection of concrete `AnalysisPipelineStep` subclasses used throughout the PyPatchy analysis system. These steps load raw simulation data (trajectories, energies, cluster files, bond graphs), convert them into structured graph or dataframe forms, and perform higher-level analyses such as cluster classification, polycube-structure matching, yield computation, and statistical aggregation across simulations. The module integrates tightly with oxDNA/oxpy readers, PyPatchy’s particle/patch representations, and igraph-based subgraph isomorphism tools. Collectively, it forms the core library of reusable operations that users compose into end-to-end analysis pipelines. 

# base_param_val.py
Defines the data structures used to represent simulation and analysis parameters throughout PyPatchy. The core class, `ParameterValue`, stores a single named parameter and a primitive value (int/float/bool/str), supporting equality, comparison, hashing, and readable string forms. More complex grouped parameters are handled by `ParamValueGroup`, which aggregates multiple `ParameterValue` objects under a shared name while still behaving like a single parameter for comparison and pipeline operations. Custom exceptions (`InvalidParameterError`, `InvalidParameterTypeError`) help ensure parameters are well-formed. These classes form the foundation for parameter handling in simulation specifications and analysis pipelines. 

# dna_particle.py
Defines `DNAParticle`, a wrapper around an oxDNA `DNAStructure` that enables representing DNA origami components as patchy-particle analogs. The class tracks sets of strands corresponding to patch sites, maps those strands to patchy-particle patches, and provides geometric utilities for computing patch centroids, normals, 3'/5' ends, and alignment transforms. It supports scaling, rotation, and translation so that DNA structures can be positioned to match patchy-particle instances during origami-to-patchy conversion workflows. The module also includes `PatchyOriRelation`, a small dataclass for storing the mapping and alignment quality between a patchy particle and a DNAParticle. 

# ensemble_parameter.py
Implements the machinery for representing, constructing, and managing variable parameters within a simulation ensemble. It provides a flexible `parameter_value` factory that converts raw JSON-like dictionaries into appropriately typed `ParameterValue` or `ParamValueGroup` objects, including special handling for particle-set parameters, multidentate-conversion settings, and staged-assembly specifications. The `EnsembleParameter` class stores all allowed values for a given parameter (e.g., temperatures, particle sets, staging modes) and supports lookup, iteration, and directory-name generation. `EnsembleParameters` then manages collections of these parameters, ensuring unique keys and enabling easy updates or deletions. This module is central to configuring and enumerating simulation ensembles.

# mgl.py
Provides support for reading, representing, and manipulating **MGL-format** patchy-particle scenes. It defines `MGLParticle` and `MGLPatch` as lightweight implementations of the `PatchyBaseParticle` and `BasePatchType` interfaces, where particle “types” correspond to color labels rather than geometric templates. `MGLScene` loads an entire MGL file into a `Scene`, offering utilities for recentring, type extraction, patch–patch binding checks, and computing particle rotations based on SVD alignment. The module also includes `MGLParticleSet` for managing unique particle-type definitions and a `load_mgl()` function that parses MGL text files into full simulation-ready scenes. Overall, it acts as a bridge between MGL geometry files and the patchy-particle framework. 

# param_val_types.py
Defines the specialized `ParameterValue` subclasses used to represent complex simulation parameters in PyPatchy ensembles. These include:

- **`ParticleSetParam`**, which wraps a full `PLParticleSet` (mapping particle types, colors, and interactions) so that entire particle sets can be treated as single parameter values.  
- **`MDTConvertParams`**, a structured wrapper for multidentate–patch conversion settings, validating and exposing conversion parameters such as number of teeth, torsion, and radius.  
- **`StagedAssemblyParam`** and **`StageInfoParam`**, which encode multi-stage assembly protocols, including stage timing, particle-addition methods, and per-stage input-file parameters.

Together, these classes allow JSON or dictionary-based ensemble specifications to be converted into strongly typed, internally consistent objects suitable for simulation setup and pipeline logic.

# particle_adders.py
Defines the mechanisms used to introduce particles into staged-assembly simulations. It provides the abstract base class `StageParticleAdder`, which specifies an interface for retrieving particle counts to be added during a simulation stage. Concrete implementations include:

- **`RandParticleAdder`** — adds user-specified numbers of particle types, typically for randomized or unstructured assembly stages.  
- **`FromPolycubeAdder`** — generates particle additions based on one or more polycube structures, supporting replication and distance-scaling options for multi-copy insertion.  
- **`FromConfAdder`** (placeholder) — intended to derive particle additions from existing configuration files.

These adders integrate with staged-assembly parameters, allowing simulation protocols to introduce particles dynamically and reproducibly across stages.

# patchy_origami_convert.py
Contains the high-level machinery that converts a patchy-particle scene (MGL or PL format) into a complete DNA-origami representation suitable for oxDNA simulations. Built atop the `OrigamiDeabstractor` base class, `PatchyOrigamiConverter` handles:

- **Geometric alignment** between patchy particles and their DNA-origami counterparts (via SVD-superimposition, Hungarian assignment for multidentate patches, and patch-normal matching).  
- **Sticky-end generation and bonding**, automatically extending strands at 3' ends to form interparticle linkages and optionally adding unbound sticky ends.  
- **Scaling and spacing logic**, determining how patchy-particle distances map into oxDNA units using spacer lengths, sticky-end lengths, and inferred interparticle distances.  
- **Export utilities**, including writing `.top/.dat` pairs, generating combined or per-particle `.oxview` visualizations, and exporting sticky-end sequences to Excel.

Additional helpers compute patch-pairing orderings, render scene multigraphs, and iterate over generated staples. This module is one of the central components enabling conversion from abstract patchy models to physically realizable DNA nanostructures. 

# patchy_scripts.py
Collects a variety of high-level helper routines for working with patchy-particle simulation files, converting between file formats, and constructing or manipulating interaction specifications. It includes:

- Utilities for **converting unidentate (UDT) particle definitions to multidentate (MDT)** forms, reading and writing flavored patchy-particle files, and switching between Flavio and Lorenzo file formats.  
- Tools such as `int_mat_to_keywise` for transforming interaction matrices into color-keyed dictionaries.  
- Convenience functions for adding standard patchy interaction potentials or selecting the best configuration from a trajectory based on cluster contiguity.

These routines are primarily scripting conveniences layered on top of the more formal classes elsewhere in the library and are often used for preprocessing or format translation tasks. 

# patchy_sim_observable.py
Defines a lightweight, now-deprecated wrapper (`PatchySimObservable`) for specifying oxDNA observables within patchy-particle simulations. It stores configuration such as file names, print intervals, nonlinear sampling parameters, and observable column specifications, and provides helpers to write these observables into oxDNA input files or input dictionaries. The module also includes `observable_from_file()`, which loads observable definitions from JSON spec files and constructs modern `Observable` objects from `ipy_oxdna`. Overall, the file serves as a compatibility and I/O utility layer for observable definitions used in analysis steps.

# sim_ensemble_json.py
Defines the JSON serialization/deserialization logic for complete PatchySimulationEnsemble configurations. At its core is `SimEnsembleJsonEncoder`, a subclass of `PLJSONEncoder` that translates diverse PyPatchy objects—including parameter sets, ensemble parameters, staged-assembly specifications, polycube structures, server configs, observables, and full analysis pipelines—into JSON-safe dictionaries. It also performs the inverse operation: reconstructing complex objects and pipeline steps from JSON via a type-dispatch system.

The module is the glue enabling ensembles, analysis workflows, and particle/type specifications to be saved, transferred, and restored intact, making it central to reproducible simulations and analysis automation within the PyPatchy ecosystem. 

# simulation_ensemble.py
Implements the central class **`PatchySimulationEnsemble`**, which manages the full lifecycle of patchy-particle simulations: setup, execution, staging, data organization, observable handling, metadata tracking, logging, and interfacing with the analysis pipeline. It orchestrates parameter sets, stages of assembly, input/output file generation, job submission (including SLURM awareness), and provides high-level utilities for retrieving configurations, computing interaction potentials, determining simulation progress, and running analysis steps. Supporting functions help locate and build ensembles from JSON specs or metadata files, enabling fully reproducible multi-simulation workflows.
We use the term "ensemble" here to refer to a collection of related simulations that vary over specified parameters, such as temperature or particle set. This class is the primary interface for users to define, run, and analyze patchy-particle simulation campaigns within PyPatchy. 
**The use of "ensemble" should not be confused with statistical ensembles in physics, and we regret the confusion.**

 # simulation_specification.py
defines the core parameter container class **`ParamSet`**, which represents the complete set of parameters defining a single patchy-particle simulation. A `ParamSet` stores both simple and grouped `ParameterValue` objects, exposes dictionary-style access to parameter values, and provides methods for merging, updating, hashing, iteration, and generating filesystem folder names. It also implements meaningful equality and order-independent comparison for use in ensemble management. The module aliases `PatchySimulation` to `ParamSet` for backward compatibility, provides a `NoSuchParamError` for missing-parameter handling, and includes `get_param_set()` for loading parameter sets from JSON spec files. Overall, it is the fundamental specification layer used by higher-level ensemble and simulation orchestration code.

# stage.py
Implements the **Stage** class, the core unit of PyPatchy’s staged-assembly simulation framework. Each `Stage` represents a contiguous segment of an oxDNA simulation with its own input parameters, duration, particle additions, box size, and file structure. A stage knows its predecessor and successor, computes start/end times, manages particle insertions (random, polycube-based, or configuration-based), and applies interaction potentials and box-size adjustments. It also generates input files, writes topology/conf data, and integrates with the ensemble’s parameter system.

The module further defines several staging-related exceptions (`IncompleteStageError`, `NoStageTrajError`, etc.) used to signal issues in multi-stage simulations. Overall, `stage.py` provides the operational backbone for constructing, validating, and executing multi-step patchy-particle assembly protocols. 

# vis_lib.py
Provides a suite of visualization utilities for inspecting simulation results and analysis outputs from PyPatchy. It focuses on high-level plotting workflows built on seaborn, matplotlib, and networkx, enabling users to:

- Visualize **time-series analysis outputs** (e.g., yield curves, cluster-size metrics) across ensemble parameters using faceted line plots.  
- Compare results across **multiple ensembles** with shared parameter spaces.  
- Display cluster graphs at particular timepoints, including particle-type coloration.  
- Render derived metrics such as total graph size over time, energy curves, and histogram summaries of final cluster sizes.

The module also defines lightweight figure-wrapper classes (`PolycubesFigure`, `YieldCurveFigure`) to standardize plot creation. Overall, it supplies the plotting layer for exploratory analysis of patchy-particle assembly behavior.

Please note that many of these functions are deprecated and may be very buggy. We invite any users who encouter issues to [reach out directly](mailto:jrevan21@asu.edu).