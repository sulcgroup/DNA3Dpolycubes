#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// --- Vector3 Class ---
struct Vector3 {
  double x, y, z;

  Vector3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

  Vector3 operator+(const Vector3 &other) const {
    return Vector3(x + other.x, y + other.y, z + other.z);
  }
  Vector3 operator-(const Vector3 &other) const {
    return Vector3(x - other.x, y - other.y, z - other.z);
  }
  Vector3 operator*(double s) const { return Vector3(x * s, y * s, z * s); }
  Vector3 operator/(double s) const { return Vector3(x / s, y / s, z / s); }
  Vector3 &operator+=(const Vector3 &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }
  Vector3 &operator-=(const Vector3 &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }

  double dot(const Vector3 &other) const {
    return x * other.x + y * other.y + z * other.z;
  }
  Vector3 cross(const Vector3 &other) const {
    return Vector3(y * other.z - z * other.y, z * other.x - x * other.z,
                   x * other.y - y * other.x);
  }
  double norm() const { return std::sqrt(x * x + y * y + z * z); }
  double norm2() const { return x * x + y * y + z * z; }
  Vector3 normalized() const {
    double n = norm();
    return (n > 1e-10) ? *this / n : Vector3(0, 0, 0);
  }
};

Vector3 operator*(double s, const Vector3 &v) { return v * s; }

// --- Matrix3 Class (for Inertia Tensor) ---
struct Matrix3 {
  double m[3][3];

  Matrix3() {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        m[i][j] = 0;
  }

  static Matrix3 identity() {
    Matrix3 res;
    res.m[0][0] = res.m[1][1] = res.m[2][2] = 1;
    return res;
  }

  Vector3 operator*(const Vector3 &v) const {
    return Vector3(m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
                   m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
                   m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z);
  }

  Matrix3 operator+(const Matrix3 &other) const {
    Matrix3 res;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        res.m[i][j] = m[i][j] + other.m[i][j];
    return res;
  }

  // Inverse using determinant (assuming symmetric/invertible)
  Matrix3 inverse() const {
    double det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
                 m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                 m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    Matrix3 res;
    if (std::abs(det) < 1e-10)
      return res; // Return zero matrix if singular

    double invDet = 1.0 / det;
    res.m[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invDet;
    res.m[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invDet;
    res.m[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invDet;
    res.m[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invDet;
    res.m[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invDet;
    res.m[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invDet;
    res.m[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invDet;
    res.m[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invDet;
    res.m[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invDet;
    return res;
  }
};

// --- Quaternion Class ---
struct Quaternion {
  double w, x, y, z;

  Quaternion(double w = 1, double x = 0, double y = 0, double z = 0)
      : w(w), x(x), y(y), z(z) {}

  Quaternion normalized() const {
    double n = std::sqrt(w * w + x * x + y * y + z * z);
    return (n > 1e-10) ? Quaternion(w / n, x / n, y / n, z / n)
                       : Quaternion(1, 0, 0, 0);
  }

  Vector3 rotate(const Vector3 &v) const {
    // v' = q * v * q_inv
    // Optimized implementation
    double t2 = w * x;
    double t3 = w * y;
    double t4 = w * z;
    double t5 = -x * x;
    double t6 = x * y;
    double t7 = x * z;
    double t8 = -y * y;
    double t9 = y * z;
    double t10 = -z * z;
    return Vector3(
        2 * ((t8 + t10) * v.x + (t6 - t4) * v.y + (t3 + t7) * v.z) + v.x,
        2 * ((t4 + t6) * v.x + (t5 + t10) * v.y + (t9 - t2) * v.z) + v.y,
        2 * ((t7 - t3) * v.x + (t2 + t9) * v.y + (t5 + t8) * v.z) + v.z);
  }

  // Integration: q_new = q + 0.5 * w_vec * q * dt
  Quaternion integrate(const Vector3 &omega, double dt) const {
    // dq/dt = 0.5 * omega * q
    // omega as quaternion (0, wx, wy, wz)
    double nw = -0.5 * (omega.x * x + omega.y * y + omega.z * z);
    double nx = 0.5 * (omega.x * w + omega.y * z - omega.z * y);
    double ny = 0.5 * (omega.y * w + omega.z * x - omega.x * z);
    double nz = 0.5 * (omega.z * w + omega.x * y - omega.y * x);

    return Quaternion(w + nw * dt, x + nx * dt, y + ny * dt, z + nz * dt)
        .normalized();
  }
};

// --- Particle Struct ---
struct Particle {
  int id;
  int strand_id;
  std::string base_name;
  int n3, n5; // Neighbors
  int cluster_id;

  Vector3 pos;
  Vector3 a1, a3;
  Vector3 L, v; // Angular momentum, Linear velocity

  // Relative to Cluster
  Vector3 rel_pos;
  Vector3 rel_a1, rel_a3; // In body frame
};

// --- Cluster Struct ---
struct Cluster {
  int id;
  std::vector<int> particle_ids;
  double mass;
  double radius;

  Vector3 com;
  Quaternion orientation; // Orientation of body frame relative to world
  Vector3 P;              // Linear Momentum
  Vector3 L;              // Angular Momentum

  Matrix3 I_body; // Inertia Tensor in body frame
  Matrix3 I_inv_body;

  Vector3 force;
  Vector3 torque;

  // Computed state
  Vector3 v;
  Vector3 omega;
};

// --- Simulation Class ---
class Simulation {
public:
  std::vector<Particle> particles;
  std::vector<Cluster> clusters;

  // Parameters
  int steps;
  double k_spring;
  double b_damp;
  double repulsion_k;
  double repulsion_offset;
  double r0_spring;  // Equilibrium length for spring
  double dt = 0.005; // Time step

  std::string last_conf_file;
  std::string topology_file;
  std::string conf_file;
  Vector3 box;

  void readInput(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      std::cerr << "Error opening input file: " << filename << std::endl;
      exit(1);
    }
    std::string line;
    while (std::getline(file, line)) {
      if (line.empty() || line[0] == '#')
        continue;
      std::stringstream ss(line);
      std::string key, val;
      if (std::getline(ss, key, '=') && std::getline(ss, val)) {
        if (key == "steps")
          steps = std::stoi(val);
        else if (key == "k")
          k_spring = std::stod(val);
        else if (key == "b")
          b_damp = std::stod(val);
        else if (key == "repulsion")
          repulsion_k = std::stod(val);
        else if (key == "repulsion_offset")
          repulsion_offset = std::stod(val);
        else if (key == "r0")
          r0_spring = std::stod(val);
        else if (key == "last_conf")
          last_conf_file = val;
        else if (key == "dt")
          dt = std::stod(val);
        else if (key == "topology")
          topology_file = val;
        else if (key == "conf_file")
          conf_file = val;
      }
    }
  }

  void readTopology(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      std::cerr << "Error opening topology file: " << filename << std::endl;
      exit(1);
    }
    int n_particles, n_strands, n_clusters;
    file >> n_particles >> n_strands >> n_clusters;

    particles.resize(n_particles);
    clusters.resize(n_clusters);
    for (int i = 0; i < n_clusters; ++i)
      clusters[i].id = i;

    for (int i = 0; i < n_particles; ++i) {
      particles[i].id = i;
      file >> particles[i].strand_id >> particles[i].base_name >>
          particles[i].n3 >> particles[i].n5 >> particles[i].cluster_id;
      if (particles[i].cluster_id >= 0 &&
          particles[i].cluster_id < n_clusters) {
        clusters[particles[i].cluster_id].particle_ids.push_back(i);
      }
    }
  }

  void readConf(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      std::cerr << "Error opening conf file: " << filename << std::endl;
      exit(1);
    }
    std::string line;
    // Skip metadata lines or parse them
    std::getline(file, line); // t = 0
    std::getline(file, line); // b = ...
    if (line.substr(0, 4) == "b = ") {
      std::stringstream ss(line.substr(4));
      ss >> box.x >> box.y >> box.z;
    }
    std::getline(file, line); // E = ...

    for (int i = 0; i < particles.size(); ++i) {
      if (!std::getline(file, line))
        break;
      std::stringstream ss(line);
      Particle &p = particles[i];
      ss >> p.pos.x >> p.pos.y >> p.pos.z >> p.a1.x >> p.a1.y >> p.a1.z >>
          p.a3.x >> p.a3.y >> p.a3.z >> p.L.x >> p.L.y >> p.L.z >> p.v.x >>
          p.v.y >> p.v.z;
    }
  }

  void initRigidBodies() {
    for (auto &cluster : clusters) {
      if (cluster.particle_ids.empty())
        continue;

      // 1. Calculate COM
      Vector3 com(0, 0, 0);
      double total_mass = 0;
      for (int pid : cluster.particle_ids) {
        com += particles[pid].pos; // Assuming mass=1
        total_mass += 1.0;
      }
      cluster.com = com / total_mass;
      cluster.mass = total_mass;

      // 2. Calculate Radius and Inertia
      cluster.radius = 0;
      Matrix3 I;
      for (int pid : cluster.particle_ids) {
        Vector3 r = particles[pid].pos - cluster.com;
        double dist = r.norm();
        if (dist > cluster.radius)
          cluster.radius = dist;

        // Inertia: m * (r^2 * Identity - r * r^T)
        double r2 = r.dot(r);
        I.m[0][0] += r2 - r.x * r.x;
        I.m[0][1] += -r.x * r.y;
        I.m[0][2] += -r.x * r.z;
        I.m[1][0] += -r.y * r.x;
        I.m[1][1] += r2 - r.y * r.y;
        I.m[1][2] += -r.y * r.z;
        I.m[2][0] += -r.z * r.x;
        I.m[2][1] += -r.z * r.y;
        I.m[2][2] += r2 - r.z * r.z;
      }
      // Add particle intrinsic inertia? Assuming point particles for now for
      // simplicity, or rather, the cluster inertia dominates.

      // 3. Set Orientation (Identity initially, as we define body frame =
      // current frame)
      cluster.orientation = Quaternion(1, 0, 0, 0);
      cluster.I_body = I; // Since orientation is identity, I_body = I_current
      cluster.I_inv_body = I.inverse();

      // 4. Store relative positions
      for (int pid : cluster.particle_ids) {
        particles[pid].rel_pos = particles[pid].pos - cluster.com;
        particles[pid].rel_a1 = particles[pid].a1; // Since R=I
        particles[pid].rel_a3 = particles[pid].a3;
      }

      // 5. Initialize Momentum (Zero for relaxation, or from particles?)
      // User said "relax DNA strand quickly". Zero velocity is best for
      // minimization/relaxation start.
      cluster.P = Vector3(0, 0, 0);
      cluster.L = Vector3(0, 0, 0);
      cluster.v = Vector3(0, 0, 0);
      cluster.omega = Vector3(0, 0, 0);
    }
  }

  void computeForces() {
    // Reset forces
    for (auto &cluster : clusters) {
      cluster.force = Vector3(0, 0, 0);
      cluster.torque = Vector3(0, 0, 0);
    }

    // 1. Connection Forces (Damped Harmonic Oscillator)
    for (const auto &p_i : particles) {
      // Check 3' neighbor
      if (p_i.n3 != -1) {
        int j_idx = p_i.n3;
        if (j_idx < 0 || j_idx >= particles.size())
          continue; // Should not happen
        const Particle &p_j = particles[j_idx];

        // Only if in different clusters
        if (p_i.cluster_id != p_j.cluster_id) {
          Vector3 r_ij = p_i.pos - p_j.pos;
          double dist = r_ij.norm();
          Vector3 dir = r_ij.normalized();

          // Spring Force: F = -k(|r| - r0) * dir
          // Force on i
          Vector3 f_spring = -k_spring * (dist - r0_spring) * dir;

          // Damping Force: F = -b * (v_i - v_j)
          Vector3 v_i = clusters[p_i.cluster_id].v +
                        clusters[p_i.cluster_id].omega.cross(
                            p_i.pos - clusters[p_i.cluster_id].com);
          Vector3 v_j = clusters[p_j.cluster_id].v +
                        clusters[p_j.cluster_id].omega.cross(
                            p_j.pos - clusters[p_j.cluster_id].com);
          Vector3 f_damp = -b_damp * (v_i - v_j);

          Vector3 total_f = f_spring + f_damp;

          // Apply to Cluster i
          clusters[p_i.cluster_id].force += total_f;
          clusters[p_i.cluster_id].torque +=
              (p_i.pos - clusters[p_i.cluster_id].com).cross(total_f);

          // Apply to Cluster j (Newton's 3rd Law)
          clusters[p_j.cluster_id].force -= total_f;
          clusters[p_j.cluster_id].torque +=
              (p_j.pos - clusters[p_j.cluster_id].com).cross(-1.0 * total_f);
        }
      }
    }

    // 2. Repulsion Forces
    for (int i = 0; i < clusters.size(); ++i) {
      for (int j = i + 1; j < clusters.size(); ++j) {
        Cluster &c1 = clusters[i];
        Cluster &c2 = clusters[j];

        Vector3 r_12 = c1.com - c2.com;
        double d = r_12.norm();
        double limit = c1.radius + c2.radius + repulsion_offset;

        if (d < limit && d > 1e-10) {
          double f_mag = repulsion_k * (1.0 - d / limit);

          Vector3 f_vec = f_mag * r_12.normalized();

          c1.force += f_vec;
          c2.force -= f_vec;
          // No torque for COM-COM force
        }
      }
    }

    // 3. Global Rotational Damping
    // To help relaxation, we damp the angular velocity of each cluster.
    // Torque_damp = -b * omega.
    // We scale b by radius^2 to have consistent units/magnitude with linear
    // damping if desired, or just use b_damp as a simple coefficient. Given b
    // is 0.2 and radius ~ 1-5, b * r^2 is reasonable.
    for (auto &cluster : clusters) {
      if (cluster.particle_ids.empty())
        continue;
      // Using b_damp directly or scaled. Let's use b_damp * radius^2 for better
      // scaling.
      double damping_factor =
          b_damp *
          (cluster.radius > 1.0 ? cluster.radius * cluster.radius : 1.0);
      cluster.torque -= damping_factor * cluster.omega;
    }
  }

  void integrate() {
    for (auto &cluster : clusters) {
      if (cluster.particle_ids.empty())
        continue;

      // Update Momentum
      cluster.P += cluster.force * dt;
      cluster.L += cluster.torque * dt;

      // Update Velocity
      cluster.v = cluster.P / cluster.mass;

      // Update Omega
      // I_inv_world = R * I_inv_body * R^T
      // But actually we can just convert L to body frame, apply I_inv_body,
      // then convert back? L_body = R^T * L_world omega_body = I_inv_body *
      // L_body omega_world = R * omega_body

      Vector3 L_body = cluster.orientation.rotate(Vector3(
          cluster.L.x, -cluster.L.y, -cluster.L.z)); // Conjugate rotate? No.
      // Rotate vector v by q: v' = q v q*.
      // Inverse rotate: v = q* v' q.
      // Let's use a helper for inverse rotate or just conjugate.
      Quaternion q_conj(cluster.orientation.w, -cluster.orientation.x,
                        -cluster.orientation.y, -cluster.orientation.z);
      L_body = q_conj.rotate(cluster.L);

      Vector3 omega_body = cluster.I_inv_body * L_body;
      cluster.omega = cluster.orientation.rotate(omega_body);

      // Update Position
      cluster.com += cluster.v * dt;

      // Update Orientation
      cluster.orientation = cluster.orientation.integrate(cluster.omega, dt);
    }
  }

  void updateParticles() {
    for (auto &cluster : clusters) {
      for (int pid : cluster.particle_ids) {
        Particle &p = particles[pid];

        // Position
        p.pos = cluster.com + cluster.orientation.rotate(p.rel_pos);

        // Orientation vectors
        p.a1 = cluster.orientation.rotate(p.rel_a1);
        p.a3 = cluster.orientation.rotate(p.rel_a3);

        // Velocity (for output)
        p.v = cluster.v + cluster.omega.cross(p.pos - cluster.com);

        // Angular momentum (approximate as particle rotating with cluster)
        // This is just for output format, maybe not physically rigorous for
        // point particles but usually oxDNA particles have inertia. We'll leave
        // L as 0 or some value? Let's just output 0 for L to be safe, or keep
        // initial? User said "start MD from there". MD needs velocities. I'll
        // set L = 0 for particles as they are part of rigid body. Or better,
        // L_particle = I_particle * omega. Assuming I_particle is identity or
        // similar.
        p.L = Vector3(0, 0, 0);
      }
    }
  }

  void writeOutput(const std::string &filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
      std::cerr << "Error opening output file: " << filename << std::endl;
      return;
    }

    file << "t = " << steps * dt << "\n"; // Approx time
    file << "b = " << box.x << " " << box.y << " " << box.z << "\n";
    file << "E = 0 0 0\n"; // Placeholder

    file << std::fixed << std::setprecision(10);
    for (const auto &p : particles) {
      file << p.pos.x << " " << p.pos.y << " " << p.pos.z << " " << p.a1.x
           << " " << p.a1.y << " " << p.a1.z << " " << p.a3.x << " " << p.a3.y
           << " " << p.a3.z << " " << p.L.x << " " << p.L.y << " " << p.L.z
           << " " << p.v.x << " " << p.v.y << " " << p.v.z << "\n";
    }
  }

  void run() {
    initRigidBodies();

    for (int step = 0; step < steps; ++step) {
      computeForces();
      integrate();
      updateParticles();
    }

    writeOutput(last_conf_file);
  }
};

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input> [topology] [conf]"
              << std::endl;
    return 1;
  }

  Simulation sim;
  sim.readInput(argv[1]);

  if (argc >= 3)
    sim.topology_file = argv[2];
  if (argc >= 4)
    sim.conf_file = argv[3];

  if (sim.topology_file.empty() || sim.conf_file.empty()) {
    std::cerr << "Error: Topology and configuration files must be specified in "
                 "input file or command line."
              << std::endl;
    return 1;
  }

  sim.readTopology(sim.topology_file);
  sim.readConf(sim.conf_file);

  sim.run();

  std::cout << "Simulation completed. Output written to " << sim.last_conf_file
            << std::endl;

  return 0;
}
