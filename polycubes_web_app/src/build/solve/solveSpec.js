import { RULE_ORDER } from "../utils";
/**
 * Data class for SAT solve specifications, should make easier to translate to/from json
 */
export class SolveSpec {
    // const params
    // num species (can be single value or range)
    nS;
    // num colors
    nC;
    // num locations
    nL;
    // num dimensions (usually 3)
    nD;
    // num patches per particle
    nP;
    // are patches torsional
    torsion;
    // num possible orientations for a tortional patch
    nO;
    // crystal topology
    bindings;
    nanoparticles;
    extraConnections;
    maxAltTries;
    forbid_self_interact;
    crystal;
    constructor(topology, nDim = 3, tortionalPatches = true) {
        this.bindings = compute_topology(topology);
        this.nC = undefined;
        this.nS = undefined;
        // Read number of particles from the topology
        this.nL = countParticles(topology); //: Numper of particle positions in the crystal lattice
        this.nD = nDim; //: Number of dimensions
        this.nP = RULE_ORDER.length; //: Number of patches on a single particle
        this.torsion = tortionalPatches; //tortionalPatches and nD > 2 // only tortion for
        if (this.torsion) {
            this.nO = 4; //: Number of possible orientations for (const a patch, N,S,W,E
        }
        // extra connections
        this.extraConnections = [];
        // whether this is a crystal (usual equivalent to this.extraConnections.length > 0)
        this.crystal = undefined;
        // do we allow patches that interact with other patches on the same particle?
        this.forbid_self_interact = false;
        // nanoparticles?
        this.nanoparticles = new Map();
        // this.nNPTypes = 0
        this.maxAltTries = 10;
    }
    /**
     * warning: code translated from python by chatGPT
     */
    getEmpties() {
        const ids = new Set();
        for (const [[i,], [j,]] of this.bindings.entries()) {
            ids.add(i);
            ids.add(j);
        }
        const patches = new Set();
        for (const [[i, dPi], [j, dPj]] of this.bindings.entries()) {
            patches.add([i, dPi]);
            patches.add([j, dPj]);
        }
        const empty = [];
        ids.forEach(i => {
            for (let dPi = 0; dPi < 6; dPi++) {
                if (!patches.has([i, dPi])) {
                    empty.push([i, dPi]);
                }
            }
        });
        return empty;
    }
    setExtraConnections(bindings) {
        this.extraConnections = bindings;
        this.crystal = bindings.length > 0;
    }
    set_nanoparticles(nps) {
        this.nanoparticles = nps;
    }
    set_nanoparticle(loc_idx, nptype) {
        this.nanoparticles.set(loc_idx, nptype);
    }
    /**
     * most important result of this function is "one" or "more than one"
     */
    num_nanoparticle_types() {
        return (new Set(this.nanoparticles.values())).size;
    }
    has_nanoparticles() {
        return this.nanoparticles.size > 0;
    }
    assign_nS(s) {
        this.nS = s;
    }
    assign_nC(c) {
        this.nC = c;
    }
    toJSON() {
        const { extraConnections, crystal, nanoparticles, nS, nC, nL, nD, nP, forbid_self_interact, ...export_vars } = this;
        if (this.forbid_self_interact) {
            export_vars["forbid_self_interact"] = this.forbid_self_interact;
        }
        if (this.crystal != void 0) {
            export_vars["crystal"] = this.crystal;
        }
        if (this.extraConnections.length > 0) {
            export_vars["extraConnections"] = this.extraConnections;
        }
        if (this.nanoparticles) {
            export_vars["nanoparticles"] = this.nanoparticles;
        }
        if (this.nC != void 0) {
            export_vars["nC"] = this.nC;
        }
        if (this.nS != void 0) {
            export_vars["nS"] = this.nS;
        }
        return export_vars;
    }
    /**
     * manually written deepcopy method because this language is a goddamn joke
     */
    cpy() {
        let copy = new SolveSpec([], this.nD, this.torsion);
        copy.bindings = new Map(this.bindings);
        copy.nS = this.nS;
        copy.nC = this.nC;
        copy.nP = this.nP;
        copy.nO = this.nO;
        copy.nL = this.nL;
        copy.extraConnections = this.extraConnections.map((b) => [...b]);
        copy.crystal = this.crystal;
        copy.forbid_self_interact = this.forbid_self_interact;
        copy.nanoparticles = new Map(this.nanoparticles);
        return copy;
    }
}
function compute_topology(bindings) {
    /**
     Accepts an array of integer tuples bindings, of format [particle_id1,patch1,particle_id_2,patch2], where particle_id1 uses patch1 to bind to particle_id2 on patch2
     Each interacting pair is only to be listed once
     */
    let top = new Map();
    for (const [p1, s1, p2, s2] of bindings) {
        top.set([p1, s1], [p2, s2]);
    }
    return top;
}
function countParticles(topology) {
    let particles;
    particles = topology.map(x => x[0]).concat(topology.map(x => x[2]));
    return Math.max(...particles) + 1;
}
//# sourceMappingURL=solveSpec.js.map