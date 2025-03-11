import {RULE_ORDER} from "../utils";

export type Binding = [number, number, number, number]

/**
 * Data class for SAT solve specifications, should make easier to translate to/from json
 */
export class SolveSpec {
    // const params
    // num species (can be single value or range)
    protected nS: number | [number, number];
    // num colors
    protected nC: number | [number, number];

    // num locations
    protected nL: number;
    // num dimensions (usually 3)
    protected nD: number;
    // num patches per particle
    protected nP: number;
    // are patches torsional
    public torsion: boolean;
    // num possible orientations for a tortional patch
    protected nO: number;
    // crystal topology
    bindings: Map<[number, number], [number, number]>;

    protected nanoparticles: Map<number, number>;
    extraConnections: Binding[];

    maxAltTries: number;

    public forbid_self_interact: boolean;
    crystal: boolean | undefined;

    constructor(topology: Binding[],
                nDim = 3, tortionalPatches = true) {
        this.bindings = compute_topology(topology);
        this.nC = undefined;
        this.nS = undefined;

        // Read number of particles from the topology
        this.nL = countParticles(topology);  //: Numper of particle positions in the crystal lattice
        this.nD = nDim;  //: Number of dimensions
        this.nP = RULE_ORDER.length;   //: Number of patches on a single particle
        this.torsion = tortionalPatches; //tortionalPatches and nD > 2 // only tortion for
        if (this.torsion) {
            this.nO = 4   //: Number of possible orientations for (const a patch, N,S,W,E
        }
        // extra connections
        this.extraConnections = [];
        // whether this is a crystal (usual equivalent to this.extraConnections.length > 0)
        this.crystal = undefined;

        // do we allow patches that interact with other patches on the same particle?
        this.forbid_self_interact = false;

        // nanoparticles?
        this.nanoparticles = new Map<number, number>();
        // this.nNPTypes = 0

        this.maxAltTries = 10;
    }

    /**
     * warning: code translated from python by chatGPT
     */
    public getEmpties(): [number, number][] {
        const ids = new Set<number>();
        for (const [[i,], [j, ]] of this.bindings.entries()) {
            ids.add(i);
            ids.add(j);
        }

        const patches = new Set<[number, number]>();
        for (const [[i, dPi], [j, dPj]] of this.bindings.entries()) {
            patches.add([i, dPi]);
            patches.add([j, dPj]);
        }

        const empty: [number, number][] = [];
        ids.forEach(i => {
            for (let dPi = 0; dPi < 6; dPi++) {
                if (!patches.has([i, dPi])) {
                    empty.push([i, dPi]);
                }
            }
        });

        return empty;
    }


    public setExtraConnections(bindings: Binding[]) {
        this.extraConnections = bindings;
        this.crystal = bindings.length > 0;
    }

    public set_nanoparticles(nps: Map<number, number>) {
        this.nanoparticles = nps;
    }

    public set_nanoparticle(loc_idx: number, nptype: number) {
        this.nanoparticles.set(loc_idx, nptype);
    }

    /**
     * most important result of this function is "one" or "more than one"
     */
    public num_nanoparticle_types(): number {
        return (new Set(this.nanoparticles.values())).size;
    }

    public has_nanoparticles(): boolean {
        return this.nanoparticles.size > 0;
    }

    public assign_nS(s: number | [number, number]) {
        this.nS = s;
    }

    public assign_nC(c: number | [number, number]) {
        this.nC = c;
    }

    public toJSON() {
        const {
            extraConnections,
            crystal,
            nanoparticles,
            nS,
            nC,
            nL,
            nD,
            nP,
            forbid_self_interact, ...export_vars
        } = this;
        if (this.forbid_self_interact) {
            export_vars["forbid_self_interact"] = this.forbid_self_interact;
        }
        if (this.crystal != void 0) {
            export_vars["crystal"] = this.crystal
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
    cpy(): SolveSpec {
        let copy = new SolveSpec([], this.nD, this.torsion);
        copy.bindings = new Map(this.bindings);
        copy.nS = this.nS;
        copy.nC = this.nC;
        copy.nP = this.nP;
        copy.nO = this.nO
        copy.nL = this.nL;
        copy.extraConnections = this.extraConnections.map((b) => [...b]);
        copy.crystal = this.crystal;
        copy.forbid_self_interact = this.forbid_self_interact;
        copy.nanoparticles = new Map(this.nanoparticles);
        return copy;
    }
}

function compute_topology(bindings: Binding[]): Map<[number, number], [number, number]> {
    /**
     Accepts an array of integer tuples bindings, of format [particle_id1,patch1,particle_id_2,patch2], where particle_id1 uses patch1 to bind to particle_id2 on patch2
     Each interacting pair is only to be listed once
     */
    let top = new Map<[number, number], [number, number]>();
    for (const [p1, s1, p2, s2] of bindings) {
        top.set([p1, s1], [p2, s2]);
    }
    return top;
}

function countParticles(topology) {
    let particles;
    particles = topology.map(x => x[0]).concat(topology.map(x => x[2]));
    return Math.max(...particles) + 1
}