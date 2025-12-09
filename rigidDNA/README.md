# RigidDNA: Rigid Body Simulation for DNA

This project implements a rigid body dynamics simulation specifically designed for DNA strands. It groups particles into "clusters" (acting as rigid bodies) and simulates their interaction through springs, damping, and repulsion. This is particularly useful for relaxing initial DNA configurations or simulating large-scale motions where base-pair level detail can be rigidified.

## Table of Contents
- [Compilation](#compilation)
- [Usage](#usage)
- [Input Parameters](#input-parameters)
- [Algorithm Detailed Explanation](#algorithm-detailed-explanation)
- [Technical Appendix: Mathematical Formulation](#technical-appendix-mathematical-formulation)

## Compilation

The simulation core is written in C++. To compile:

1. Navigate to the source directory:
   ```bash
   cd rigid_body_sim
   ```
2. Run `make`:
   ```bash
   make
   ```
   This will produce the `RigidBodySim` executable.

## Usage

After compilation, you can run the simulation using an input parameter file.

**Basic usage:**
```bash
./rigid_body_sim/RigidBodySim <path_to_input_file>
```

**Overriding files via command line:**
```bash
./rigid_body_sim/RigidBodySim <input_file> <topology_file> <conf_file>
```

### Example
An example is provided in the `example/` directory.
```bash
cd example
../rigid_body_sim/RigidBodySim input
```
This will read `input`, load `output.top` and `output.dat`, run the simulation, and write the result to `last_conf.dat`.

## Input Parameters

The input file uses a `key=value` format (with optional comments starting with `#`).

| Parameter | Description |
|-----------|-------------|
| `steps` | Number of simulation steps to run. |
| `dt` | Time step size. |
| `k` | Spring constant for connections between clusters. |
| `b` | Damping coefficient (linear and rotational). |
| `r0` | Equilibrium length for the connecting springs. |
| `repulsion` | Force constant for soft repulsion between clusters. |
| `repulsion_offset` | Offset distance for repulsion check. |
| `topology` | Path to the topology file (defines particles and clusters). |
| `conf_file` | Path to the initial configuration file (positions/velocities). |
| `last_conf` | Path where the final configuration will be written. |

## Algorithm Detailed Explanation

The simulation treats groups of particles as **Rigid Bodies** (Clusters). Instead of integrating equations of motion for every particle, it integrates the motion of the center of mass (COM) and the orientation (Quaternion) of each cluster.

### 1. Initialization
- **Clustering**: Particles are grouped based on `cluster_id` provided in the topology file.
- **Rigid Body Properties**: For each cluster, the code calculates:
  - **Center of Mass (COM)**: $\mathbf{R}_{cm} = \frac{1}{M} \sum m_i \mathbf{r}_i$
  - **Inertia Tensor ($I_{body}$)**: Computed relative to the COM in the initial frame.
  - **Radius**: The maximum distance from COM to any particle in the cluster (used for collision detection).
  - **Orientation ($q$)**: Initialized to identity $[1, 0, 0, 0]$.

### 2. Force Computation
Forces are accumulated on the COM of each cluster, and torques are accumulated based on the impact point relative to the COM.

- **Spring Connections**:
  - Acts between bonded neighbors (3' connection `n3`) if they belong to *different* clusters.
  - **Formula**: $\mathbf{F}_{spring} = -k (|\mathbf{r}_{ij}| - r_0) \hat{\mathbf{r}}_{ij}$
  - Generates both force on the COM and torque: $\boldsymbol{\tau} = (\mathbf{r}_{contact} - \mathbf{R}_{cm}) \times \mathbf{F}$.

- **Damping**:
  - Drag force proportional to the relative velocity between connected points.
  - **Formula**: $\mathbf{F}_{damp} = -b (\mathbf{v}_i - \mathbf{v}_j)$
  - Helps dissipate energy, useful for relaxation.
  - **Rotational Damping**: A torque is applied to dampen angular velocity: $\boldsymbol{\tau}_{damp} = -(b \cdot R^2) \boldsymbol{\omega}$.

- **Repulsion**:
  - Soft repulsive force if clusters overlap (distance $d < R_1 + R_2 + \text{offset}$).
  - **Formula**: $\mathbf{F}_{rep} = k_{rep} (1 - \frac{d}{R_{limit}}) \hat{\mathbf{r}}_{12}$.
  - Prevents rigid bodies from passing through each other.

### 3. Integration Scheme
The simulation uses a semi-implicit Euler or symplectic-like integration step for rigid body dynamics:

1.  **Update Momentum**:
    - $\mathbf{P}(t+\Delta t) = \mathbf{P}(t) + \mathbf{F} \Delta t$
    - $\mathbf{L}(t+\Delta t) = \mathbf{L}(t) + \boldsymbol{\tau} \Delta t$

2.  **Update Velocities**:
    - **Linear**: $\mathbf{v} = \mathbf{P} / M$
    - **Angular**:
        - Convert Angular Momentum to body frame: $\mathbf{L}_{body} = q^{-1} \mathbf{L}_{world} q$
        - Apply inverse inertia: $\boldsymbol{\omega}_{body} = I_{body}^{-1} \mathbf{L}_{body}$
        - Convert back to world/angular velocity: $\boldsymbol{\omega} = q \boldsymbol{\omega}_{body} q^{-1}$ (or handled via quaternion algebra).

3.  **Update State**:
    - **Position**: $\mathbf{R}_{cm} += \mathbf{v} \Delta t$
    - **Orientation**: Integration of quaternion using $\boldsymbol{\omega}$.
        - $\dot{q} = \frac{1}{2} \boldsymbol{\omega} q$
        - $q(t+\Delta t) = q(t) + \dot{q} \Delta t$ (Normalized).

### 4. Particle Update (Reconstruction)
Finally, the positions and velocities of individual constituent particles are reconstructed from their cluster's state:
- $\mathbf{r}_i(t) = \mathbf{R}_{cm}(t) + q(t) \mathbf{r}_{i, rel} q^{-1}(t)$
- $\mathbf{v}_i(t) = \mathbf{v}_{cm}(t) + \boldsymbol{\omega}(t) \times \mathbf{r}_i(t)$

This allows the output file (`last_conf`) to remain compatible with standard viewing or analysis tools that expect particle data.

## Technical Appendix: Mathematical Formulation

### Overview
This appendix defines the precise mathematical formulation used in the codebase.
The simulation employs a rigid body dynamics approach to coarse-grain the DNA system. Groups of particles are defined as rigid "clusters". The motion of these clusters is governed by the symplectic integration of their center of mass (COM) translational degrees of freedom and their rotational degrees of freedom (represented by quaternions).

### Rigid Body Kinematics
For a cluster $C$ consisting of particles $\{i\}$, the gross properties are defined as:

$$
    M_C = \sum_{i \in C} m_i
$$
$$
    \mathbf{R}_{cm} = \frac{1}{M_C} \sum_{i \in C} m_i \mathbf{r}_i
$$

The inertia tensor $\mathbf{I}_{body}$ is computed in the body-fixed frame (initially aligned with the world frame) relative to the COM:
$$
    \mathbf{I}_{body} = \sum_{i \in C} m_i \left( |\tilde{\mathbf{r}}_i|^2 \mathbf{1} - \tilde{\mathbf{r}}_i \otimes \tilde{\mathbf{r}}_i \right)
$$
where $\tilde{\mathbf{r}}_i = \mathbf{r}_i(t=0) - \mathbf{R}_{cm}(t=0)$.

The state of the rigid body at time $t$ is described by the tuple $(\mathbf{R}_{cm}(t), \mathbf{P}(t), q(t), \mathbf{L}(t))$, where $\mathbf{P}$ is the linear momentum, $\mathbf{L}$ is the angular momentum, and $q$ is the unit quaternion representing orientation.

### Force Fields
The simulation incorporates three primary interaction terms: harmonic bond potentials, viscous damping, and soft-core repulsion.

#### Connection Forces (Springs)
Interactions between specific particles $i$ and $j$ (typically representing varying strands or specific binding sites) spanning different clusters are modeled using a harmonic potential.

The force $\mathbf{F}_{spring}$ exerted on particle $i$ by particle $j$ is:
$$
    \mathbf{F}_{spring}^{(i)} = -k_{spring} (|\mathbf{r}_{ij}| - r_0) \hat{\mathbf{r}}_{ij}
$$
where $\mathbf{r}_{ij} = \mathbf{r}_i - \mathbf{r}_j$, $k_{spring}$ is the stiffness constant, and $r_0$ is the equilibrium bond length.

The corresponding force on the cluster COM is $\mathbf{F}_{spring}^{(i)}$, and the torque is:
$$
    \boldsymbol{\tau}_{spring} = (\mathbf{r}_i - \mathbf{R}_{cm}) \times \mathbf{F}_{spring}^{(i)}
$$

#### Viscous Damping
To facilitate relaxation and energy dissipation, a damping force proportional to the relative velocity of connected particles is applied.

$$
    \mathbf{F}_{damp}^{(i)} = -b_{damp} (\mathbf{v}_i - \mathbf{v}_j)
$$
Here, the velocity of a constituent particle $\mathbf{v}_i$ is derived from the rigid body motion:
$$
    \mathbf{v}_i = \mathbf{v}_{cm} + \boldsymbol{\omega} \times (\mathbf{r}_i - \mathbf{R}_{cm})
$$
This force also contributes to both the net force and net torque on the cluster.

#### Soft-Core Repulsion
To prevent unphysical overlap between rigid bodies, a soft repulsive force is applied between the centers of mass of any two clusters $A$ and $B$ that come within a threshold distance $R_{limit}$.

$R_{limit} = R_A + R_B + \delta_{offset}$, where $R_A$ is the bounding radius of cluster $A$.

If $d_{AB} = |\mathbf{R}_A - \mathbf{R}_B| < R_{limit}$:
$$
    \mathbf{F}_{rep}^{(A)} = k_{rep} \left( 1 - \frac{d_{AB}}{R_{limit}} \right) \hat{\mathbf{n}}_{AB}
$$
where $\hat{\mathbf{n}}_{AB}$ is the unit vector pointing from $B$ to $A$.

#### Rotational Damping
A global rotational damping torque is applied to stabilize angular velocities:
$$
    \boldsymbol{\tau}_{damp}^{rot} = - (b_{damp} \cdot R_{cluster}^2) \boldsymbol{\omega}
$$

### Time Integration
The system evolves via a symplectic-like integration scheme with time step $\Delta t$.

#### Momentum Update
$$
\begin{aligned}
    \mathbf{P}(t + \Delta t) &= \mathbf{P}(t) + \mathbf{F}_{total} \Delta t \\
    \mathbf{L}(t + \Delta t) &= \mathbf{L}(t) + \boldsymbol{\tau}_{total} \Delta t
\end{aligned}
$$

#### Velocity and Angular Velocity Update
Linear velocity is simply $\mathbf{v}(t + \Delta t) = \mathbf{P}(t + \Delta t) / M$.

Angular velocity $\boldsymbol{\omega}$ is computed by transforming $\mathbf{L}$ to the body frame to utilise the constant inertia tensor:
$$
    \mathbf{L}_{body} = q^{-1} \mathbf{L}_{world} q
$$
$$
    \boldsymbol{\omega}_{body} = \mathbf{I}_{body}^{-1} \mathbf{L}_{body}
$$
$$
    \boldsymbol{\omega}_{world} = q \boldsymbol{\omega}_{body} q^{-1}
$$

#### Position and Orientation Update
Positions update linearly:
$$
    \mathbf{R}_{cm}(t + \Delta t) = \mathbf{R}_{cm}(t) + \mathbf{v}(t + \Delta t) \Delta t
$$

Orientation updates by integrating the quaternion derivative $\dot{q} = \frac{1}{2}\boldsymbol{\omega}q$:
$$
    q(t + \Delta t) = \text{Norm}\left( q(t) + \frac{1}{2} \boldsymbol{\omega}_{world} q(t) \Delta t \right)
$$
