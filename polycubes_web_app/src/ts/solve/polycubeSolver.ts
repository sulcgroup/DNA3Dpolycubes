/*
Polycube SAT specification adapted from Lukas chrystal one.

"Atoms" are "movable and rotatable", have 6 slots
"Positions" are fixed in the crystal, have 6 slots and bind according to spec
The problem has 2 parts:
A. find color bindings and colorings of position slots where each patch neightbors according to crystal model have
    colors that bind
B. find colorings of atoms s.t. all crystal positions are identical to (some) species rotation. The species definitions
    must not allow for (const bad 2-binds

indexes:
- colors:   1...c...//c (variable number)
- atoms:    1...s...//s (variable number)
- slots:    0...p...5=//p-1 (bindings places on atoms - 0,1,2 on one side, 3,4,5 on the other)
- position: 1...l...16=//l (number of positions in the crystal)
- rotation: 1...r...6=//r possible rotations of an species
- condition: 1...d...//d conditions to avoid bad crystal
- qualification: 0..//q (0 for (const good crystal, one more for (const each bad one)

(boolean) variables:
- B(c1, c2) { color c1 binds with c2 (n=#c*#c)
- F(l, p, c) { patch p at position l has color c (n=#l*#p*#c)
- C(s, p, c) { patch p on species s has color c (n=//s*#p*#c)
- P(l, s, r) { position l is occupied by species s with rotation r (n=#l*#s*#r)

encoding functions:
- rotation(p, r) = patch that p rotates to under rotation r
*/
import * as THREE from "three"

import {combinations, FACE_ROTATIONS, getSignedAngle, range, RULE_ORDER} from "../utils"
import {DynamicPolycubeCubeType, PolycubeCubeType, ruleToDec} from "../rule";
import * as Logic from "logic-solver";
import {Binding, SolveSpec} from "./solveSpec";
import {Vector3} from "three";
import {getSysBehavior, SystemBehavior, SystemBehaviorType} from "../analysis";
import {PolycubePatch} from "../patch";

declare global {
    interface Map<K, V> {
        setdefault(key: K, default_value: V | null): V;
    }
}


//Polycube SAT Solver
export class polysat extends SolveSpec {


    // map of possible rotations (make global const?)
    private rotations: Map<any, any>;
    // num possible rotations
    private nR: number;
    //  map of sat vars
    private variables: Map<string, number>;
    private BCO_varlen: number;

    solver: Logic.Solver;

    // private basic_sat_clauses: any;

    constructor(topology: Binding[] = [],
                empty: [number, number][] = [],
                nCubeTypes: number = -1,
                nColors: number = -1, nDim=3, tortionalPatches=true) {
        super(topology, nDim, tortionalPatches);

        // Different color coding, color n binds not to -n but
        // to another m, also ignore 0 and 1.
        this.nC = nColors;

        this.nS = nCubeTypes;  //: Number of distinct particle types for (const the solver

        //problem specification:
        this.rotations = enumerateRotations(this.nD);
        this.nR = this.rotations.size;

        this.variables = new Map<string, number>(); // map of variable names (string) to variable indexes (int)
        // this.basic_sat_clauses = null;       //  string of s basic sat clause
        // this.additional_sat_clauses = null;  //  some additional conditions
        this.BCO_varlen = null;              //  the number of clauses that determine B and C

        this.solver = new Logic.Solver();
    }

    public init(){
        this.generate_constraints();

        // Solution must use all particles
        this.add_constraints_all_particles();

        // Solution must use all patches, except color 0 which should not bind
        this.add_constraints_all_patches_except(0);

        // A color cannot bind to itself
        this.add_constraints_no_self_complementarity();

        // Make sure color 0 binds to 1 and nothing else
        this.fix_color_interaction(0, 1);

        // Fix interaction matrix, to avoid redundant solution
        for (const c of range(2, this.num_colors() - 1, 2)) {
            this.fix_color_interaction(c, c + 1);
        }

        if (this.torsion) {
            this.add_constraints_fixed_blank_orientation();
        }

        for (const [particle, patch] of this.getEmpties()) {
            this.fix_slot_colors(particle, patch, 1);
            //console.log("Particle {} patch {} should be empty".format(particle, patch))
        }

        if (this.forbid_self_interact) {
            this.add_constraints_no_self_interact_patches();
        }
        if (this.has_nanoparticles()) {
            if (this.num_nanoparticle_types() > 1){
                this.add_constraints_nanoparticle_multitype();
            }
            else {
                this.add_constraints_nanoparticle_singletype();
            }
        }
    }

    public variable(name: string, ...args: any[]) : number{
        let str_key = `${name}(${args.join(",")})`;
        return this.solver.getVarNum(str_key, false);
        // return this.variables.setdefault(str_key, this.variables.size + 1)
    }

    // TODO: more thorougly integrate logic-solver
    public add_clause(...clause: number[][]){
        this.solver.require(clause);
    }

    public num_species() : number {
        return <number>this.nS;
    }

    public num_colors() : number {
        return (<number>this.nC + 1) * 2;
    }

    B (c1: number, c2: number) : number {
        // color c1 binds with c2
        if (c2 < c1) {
            [c1, c2] = [c2, c1];
        }
        //console.assert(0 <= c1 <= c2 && c1 <= c2 < this.nC, `c1 <= c2 out of bounds 0 <= ${c1 <= c2} < ${this.nC}`);
        //print >> sys.stderr, 'B({c1},{c2})'.format(c1=c1, c2=c2)
        return this.variable("B", c1, c2);
    }

    D (p1: number, o1: number, p2: number, o2: number) : number {
        // patch p1, orientation c1 binds with patch p2, orientation c2 //
        if (p2 < p1) {
            [o1, o2] = [o2, o1];
            [p1, p2] = [p2, p1];
        }
        //console.assert(0 <= p1 <= p2 && p1 <= p2 < this.nP, `p1 <= p2 out of bounds 0 <= ${p1 <= p2} < ${this.nP}`);
        console.assert(0 <= o1 && o1 < this.nO, `o1 out of bounds 0 <= ${o1} < ${this.nO}`);
        console.assert(0 <= o2 && o2 < this.nO, `o2 out of bounds 0 <= ${o2} < ${this.nO}`);
        return this.variable("D", p1, o1, p2, o2);
    }

    F (l: number, p: number, c: number) : number {
        // patch p at position l has color c //
        console.assert(0 <= l && l < this.nL, `l out of bounds 0 <= ${l} < ${this.nL}`);
        console.assert(0 <= p && p < this.nP, `p out of bounds 0 <= ${p} < ${this.nP}`);
        console.assert(0 <= c && c < this.num_colors(), `c out of bounds 0 <= ${c} < ${this.num_colors()}`);
        return this.variable("F", l, p, c)
    }

    A (l: number, p: number, o: number) : number {
        // patch p at position l has orientation o //
        console.assert(0 <= l && l < this.nL, `l out of bounds 0 <= ${l} < ${this.nL}`);
        console.assert(0 <= p && p < this.nP, `p out of bounds 0 <= ${p} < ${this.nP}`);
        console.assert(0 <= o && o < this.nO, `o out of bounds 0 <= ${o} < ${this.nO}`);
        return this.variable("A", l, p, o);
    }

    C (s: number, p: number, c: number) : number {
        // patch p on species s has color c //
        console.assert(0 <= s && s < this.num_species(), `s out of bounds 0 <= ${s} < ${this.num_species()}`);
        console.assert(0 <= p && p < this.nP, `p out of bounds 0 <= ${p} < ${this.nP}`);
        console.assert(0 <= c && c < this.num_colors(), `c out of bounds 0 <= ${c} < ${this.num_colors()}`);
        return this.variable("C", s, p, c);
    }

    O (s: number, p: number, o: number) : number{
        // patch p on species s has orientation o //
        console.assert(0 <= s && s < this.num_species(), `s out of bounds 0 <= ${s} < ${this.num_species()}`);
        console.assert(0 <= p && p < this.nP, `p out of bounds 0 <= ${p} < ${this.nP}`);
        console.assert(0 <= o && o < this.nO, `o out of bounds 0 <= ${o} < ${this.nO}`);
        return this.variable("O", s, p, o);
    }

    public P (l: number, s: number, r: number) : number {
        // position l is occupied by species s with rotation r //
        console.assert(0 <= l && l < this.nL, `l out of bounds 0 <= ${l} < ${this.nL}`);
        console.assert(0 <= s && s < this.num_species(), `s out of bounds 0 <= ${s} < ${this.num_species()}`);
        console.assert(0 <= r && r < this.nR, `r out of bounds 0 <= ${r} < ${this.nR}`);
        return this.variable("P", l, s, r);
    }

    /**
     * Single nanoparticle version of speies-has-nanoparticle clause
     * @param s
     * @constructor
     */
    public N_single(s: number) : number {
        return this.variable("N", s);
    }

    public N_multi(s: number, i: number) : number {
        return this.variable("N", s, i);
    }

    rotation (p: number, r: number) : number {
        // patch that p rotates to under rotation r //
        console.assert(0 <= p && p < this.nP, `p out of bounds 0 <= ${p} < ${this.nP}`);
        console.assert(0 <= r && r < this.nR, `r out of bounds 0 <= ${r} < ${this.nR}`);
        // console.assert(all(len(set(rotations[r].keys())) == nP for (const r in rotations))
        // console.assert(all(len(set(rotations[r].values())) == nP for (const r in rotations))
        console.assert(this.rotations.size == this.nR)
        console.assert(this.rotations.has(r))
        console.assert(this.rotations.get(r).has(p))
        return this.rotations.get(r).get(p)
    }

    orientation (p: number, r: number, o: number) : number {
        // new orientation for (const patch p with initial orientation o after getting rotated by r //
        console.assert(0 <= p && p < this.nP, `p out of bounds 0 <= ${p} < ${this.nP}`);
        console.assert(0 <= r && r < this.nR, `r out of bounds 0 <= ${r} < ${this.nR}`);
        console.assert(0 <= o && o < this.nO, `o out of bounds 0 <= ${o} < ${this.nO}`);
        console.assert(this.rotations.size == this.nR)
        console.assert(this.rotations.has(r))
        console.assert(this.rotations.get(r).has(p))

        // Calculate patch that p rotates to:
        let p_rot = this.rotations.get(r).get(p)
        // Calculate vector corresponding to o
        let v = patchRotToVec(p, o)
        // Which p has the same vector as o?
        let p_temp = RULE_ORDER.findIndex(e=>e.equals(v))
        //assert(p != p_temp)
        // How would that p get rotated?
        let p_temp_rot = this.rotations.get(r).get(p_temp)
        //assert(p_rot != p_temp_rot)
        // And what vector does that correspond to?
        let v_rot = RULE_ORDER[p_temp_rot]
        // And what orientation value does that give us?
        return patchVecToRot(p_rot, v_rot)
    }

    check_settings() {
        console.assert(this.bindings.size == (this.nL * this.nP) / 2.0)
        console.assert(new Set(this.bindings.values()).size == this.bindings.size);
        //console.assert((new Set(this.bindings) | new Set(this.bindings.values())).size == this.nL * this.nP)
        //console.assert(Math.min([s for (const s, _ in this.bindings] + [s for (const _, s in this.bindings] + [s for (const s, _ in this.bindings.values()] + [s for (const _, s in this.bindings.values()]) == 0)
        //console.assert(Math.max([s for (const s, _ in this.bindings] + [s for (const s, _ in this.bindings.values()]) == this.nL - 1)
        //console.assert(Math.max([s for (const _, s in this.bindings] + [s for (const _, s in this.bindings.values()]) == this.nP - 1)
        for (const p of range(this.nP)) {
            for (const r of range(this.nR)) {
                this.rotation(p, r);
            }
        }
    }

    exactly_one(vs: number[]) : number[][] {
        // returns a list of constraints implementing "exacly one of vs is true" //
        console.assert(vs.every(v=>v>0));
        if (vs.length == 1){
            return [vs];
        }
        vs.sort((a,b)=>{return a-b}); // Sort numerically
        let constraints = [vs];
        for (const v1 of vs) {
            for (const v2 of vs) {
                if (v2 >= v1) {
                    break;
                }
                constraints.push([-v1, -v2]);
            }
        }
        console.assert((new Set(constraints)).size == (vs.length * (vs.length-1)) / 2 + 1)
        return constraints;
    }

    generate_constraints() {
        // make sure B, C and O vars are first:
        for (const c1 of range(this.num_colors())) {
            for (const c2 of range(this.num_colors())) {
                this.B(c1, c2)
            }
        }
        for (const s of range(this.num_species())) {
            for (const p of range(this.nP)) {
                for (const c of range(this.num_colors())) {
                    this.C(s, p, c);
                }
            }
        }
        if (this.torsion) {
            for (const s of range(this.num_species())) {
                for (const p of range(this.nP)) {
                    for (const o of range(this.nO)) {
                        this.O(s, p, o);
                    }
                }
            }
        }

        //console.log('c settings: nS=%d nC=%d nP=%d ' % (nS, nC, nP) )
        //console.log('c Last B and C var number: %s' % variables.length)
        // this.basic_sat_clauses = [];
        //this.add_clause('c settings: nS=%d nC=%d nP=%d ' % (this.nS, this.nC, this.nP) )
        //this.add_clause('c Last B and C var number: %s' % this.variables.size)
        this.BCO_varlen = this.variables.size;

        // BASIC THINGS:
        // - Legal color bindings:
        // "Each color has exactly one color that it binds to"
        // 	forall c1 exactly one c2 s.t. B(c1, c2)
        for (const c1 of range(this.num_colors())) {
            this.add_clause(...this.exactly_one([...range(this.num_colors())].map(c2=>this.B(c1, c2))))
            //print >> sys.stderr, [B(c1, c2) for (const c2 of range(nC) if c2 != c1]
        }

        // - Legal species patch coloring (unnecesay, implied by "Legal species coloring in positions" and "Legal position
        //   patch coloring") {
        // "Each patch on every species has exactly one color"
        //   forall s, forall p, exactly one c p.t. C(s, p, c)
        for (const s of range(this.num_species())) {
            for (const p of range(this.nP)) {
                this.add_clause(...this.exactly_one([...range(this.num_colors())].map(c=>this.C(s, p, c))));
            }
        }

        // - Legal species patch orientation
        // "Each patch on every species has exactly one orientation"
        //   forall s, forall p, exactly one o p.t. O(s, p, o)
        if (this.torsion) {
            for (const s of range(this.num_species())) {
                for (const p of range(this.nP)) {
                    this.add_clause(...this.exactly_one([...range(this.nO)].map(o=>this.O(s, p, o))));
                }
            }
        }


        // ADD CRYSTAL and COLORS:
        // - Legal position patch coloring:
        // "Every position patch has exactly one color"
        // 	for (const all l, p exactly one c st. F(l, p, c)
        for (const l of range(this.nL)) {
            for (const p of range(this.nP)) {
                this.add_clause(...this.exactly_one([...range(this.num_colors())].map(c=>this.F(l, p, c))))
            }
        }

        // - Legal position patch orientation:
        // "Every position patch has exactly one orientation"
        // 	for (const all l, p exactly one o st. A(l, p, o)
        if (this.torsion) {
            for (const l of range(this.nL)) {
                for (const p of range(this.nP)) {
                    this.add_clause(...this.exactly_one([...range(this.nO)].map(o=>this.A(l, p, o))))
                }
            }
        }

        // - Forms desired crystal:
        // "Specified binds have compatible colors"
        // 	forall (l1, p1) binding with (l2, p2) from crystal spec:
        // 		forall c1, c2: F(l1, p1, c1) and F(l2, p2, c2) => B(c1, c2)
        for (const [[l1, p1], [l2, p2]] of this.bindings) {
            for (const c1 of range(this.num_colors())) {
                for (const c2 of range(this.num_colors())) {
                    this.add_clause([-this.F(l1, p1, c1), -this.F(l2, p2, c2), this.B(c1, c2)])
                }
            }
        }

        // - Forms desired crystal:
        // "Specified binds have compatible orientations"
        // 	forall (l1, p1) binding with (l2, p2) from crystal spec:
        // 		forall o1, o2: A(l1, p1, o1) and A(l2, p2, o2) => D(c1, c2)
        if (this.torsion) {
            for (const [[l1, p1], [l2, p2]] of this.bindings) {
                for (const o1 of range(this.nO)) {
                    for (const o2 of range(this.nO)) {
                        this.add_clause([-this.A(l1, p1, o1), -this.A(l2, p2, o2), this.D(p1, o1, p2, o2)])
                    }
                }
            }
        }

        // Hard-code patch orientations to bind only if they point in the same direction
        if (this.torsion) {
            for (const p1 of range(this.nP)) {
                for (const p2 of range(this.nP)) {
                    if (p2 >= p1) {
                        break
                    }
                    for (const o1 of range(this.nO)) {
                        for (const o2 of range(this.nO)) {
                            const v1 = patchRotToVec(p1, o1)
                            const v2 = patchRotToVec(p2, o2)
                            // Do they point in the same global direction?
                            // And do the patches face each other?
                            if (v1.equals(v2) && p2%2==0 && p2+1 == p1) {
                                this.add_clause([this.D(p1, o1, p2, o2)])
                            }
                            else {
                                this.add_clause([-this.D(p1, o1, p2, o2)])
                            }
                        }
                    }
                }
            }
        }

        // - Legal species placement in positions:
        // "Every position has exactly one species placed there with exactly one rotation"
        //   forall l: exactly one s and r p.t. P(l, s, r)
        for (const l of range(this.nL)) {
            let ps = [];
            for (const r of range(this.nR)) {
                for (const s of range(this.num_species())) {
                    ps.push(this.P(l, s, r))
                }
            }
            this.add_clause(...this.exactly_one(ps))
        }

        // - Legal species coloring in positions:
        // "Given a place, species and its rotation, the patch colors on the position and (rotated) species must be the same"
        //   for (const all l, s, r:
        //       P(l, s, r) => (forall p, c: F(l, p, c) <=> C(s, rotation(p, r), c))
        for (const l of range(this.nL)) {
            for (const s of range(this.num_species())) {
                for (const r of range(this.nR)) {
                    // forall part
                    for (const p of range(this.nP)) {
                        for (const c of range(this.num_colors())) {
                            const p_rot = this.rotation(p, r) // Patch after rotation
                            // Species 's' rotated by 'r' gets color 'c' moved from patch 'p' to 'p_rot':
                            this.add_clause([
                                -this.P(l, s, r), // EITHER no species 's' at position 'l' with rot 'r'
                                -this.F(l, p, c), // OR no patch 'p' at position 'l' with color 'c'
                                this.C(s, p_rot, c) // OR patch 'p_rot' on species 's' DOES have the color 'c'
                            ]);
                            this.add_clause([
                                -this.P(l, s, r), // EITHER no species 's' at position 'l' with rot 'r'
                                this.F(l, p, c), // OR there is a patch 'p' at position 'l' with color 'c'
                                -this.C(s, p_rot, c) // OR there is no patch 'p_rot' on species 's' with the color 'c'
                            ]);
                        }
                    }
                }
            }
        }


        // - Legal species patch orientation in positions:
        // "Given a place, species and its rotation, the patch orientations on the position and (rotated) species must be correct"
        //   for (const all l, s, r:
        //       P(l, s, r) => (forall p, c: F(l, p, c) <=> C(s, rotation(p, r), c))
        if (this.torsion) {
            for (const l of range(this.nL)) {
                for (const s of range(this.num_species())) {
                    for (const r of range(this.nR)) {
                        for (const p of range(this.nP)) {
                            for (const o of range(this.nO)) {
                                const p_rot = this.rotation(p, r) // Patch after rotation
                                const o_rot = this.orientation(p, r, o) // Patch orientation after rotation
                                // Species 's' rotated by 'r' gets orientation 'o' of patch 'p' changed to 'o_rot' at the new path 'p_rot':
                                //console.log("Species {} rotated by {}: patch {}-->{}, orientation {}-->{}".format(s, r, p, p_rot, o, o_rot))
                                this.add_clause([ 
                                    -this.P(l, s, r), // EITHER no species 's' at position 'l' with rot 'r'
                                    -this.A(l, p, o), // OR no patch 'p' at position 'l' with orientation 'o'
                                    this.O(s, p_rot, o_rot) // OR patch 'p_rot' on species 's' has the orientation 'o_rot'
                                ])
                                this.add_clause([
                                    -this.P(l, s, r), // EITHER no species 's' at position 'l' with rot 'r'
                                    this.A(l, p, o), // OR there is a patch 'p' at position 'l' with orientation 'o'
                                    -this.O(s, p_rot, o_rot) // OR there is no patch 'p_rot' on species 's' with the orientation 'o_rot'
                                ])
                            }
                        }
                    }
                }
            }
        }
        
        if (this.nD == 2) {
            // Lock patch orientation if 2D
            for (const s of range(this.num_species())) {
                for (const p of range(this.nP)) {
                    if (p>3) {
                        // Patch p on species s is empty
                        this.add_clause([this.C(s, p, 1)])
                    }
                    const o = getFlatFaceRot()[p]
                    // Patch p has orientation 'o'
                    this.add_clause([this.O(s, p, o)])
                }
            }
        }
    }


    output_cnf(constraints: any[]) {
        // Outputs a CNF formula //
        let num_vars = Math.max(...this.variables.values());
        let num_constraints = constraints.length;
        let outstr = `p cnf ${num_vars} ${num_constraints}\n`;
        for (const c of constraints) {
            outstr += c.join(' ') + ' 0\n';
        }
        return outstr;
    }


    // load_solution_from_lines(lines) {
    //     // loads solution from sat solution output in s string//
    //     let satline;
    //     if (lines.length > 1) {
    //         console.assert(lines[0].strip() == 'SAT')
    //         satline = lines[1].strip().split()
    //     }
    //     else {
    //         satline = lines[0].strip().split()
    //     }
    //
    //     //line = myinput.readline().strip()
    //     console.assert(line == 'SAT')
    //     let sols = satline.map(v=>Math.floor(v))
    //     console.assert(sols[-1] == 0)
    //     sols.pop(); //Remove last element
    //     console.assert(sols.length <= this.variables.size)
    //
    //     let vnames = []
    //     for (const [vname, vnum] in this.variables) {
    //         if (sols.includes(vnum)) {
    //             vnames.push(vname);
    //         }
    //     }
    //     return vnames;
    // }

    add_constraints_from_vnames(vnames: string) {
        let constraints = []
        for (const vname of vnames) {
            if (!(vname in this.variables)) {
                console.error("Trying to add variables that have not been defined, probably incompatible problem formulation?")
            }
            constraints.push(this.variables[vname])
        }
        this.add_clause(constraints);
    }

    // convert_solution(solution: Logic.Solution) : PolycubeCube[] {
    //     // loads solution from minisat sol in myinput handle, writes variable names to output //
    //     let lines = solution.split(' ');
    //     console.assert(lines[0] == 'SAT')
    //     let sols: any[] = lines.slice(1);
    //     //console.assert(sols[-1] == 0)
    //     sols.pop();
    //     console.assert(sols.length <= this.variables.size)
    //
    //     let output = '';
    //
    //     for (const vname of [...this.variables.keys()].sort()) {
    //         const vnum = this.variables.get(vname);
    //         if (vnum > sols.length) {
    //             break;
    //         }
    //         if (sols[vnum-1] > 0) {
    //             output += vname+'\n';
    //         }
    //     }
    //     return output;
    // }

    add_constraints_all_particles() {
        for (const s of range(this.num_species())) {
            let ps = [];
            for (const r of range(this.nR)){
                for (const l of range(this.nL)) {
                    ps.push(this.P(l,s,r));
                }
            }
            this.add_clause(ps);
        }
    }

    add_constraints_all_patches_except(forbidden: number) {
        for (const c of range(this.num_colors())) {
            if (c != forbidden) {
                let cs = [];
                for (const s of range(this.num_species())){
                    for (const p of range(this.nP)) {
                        cs.push(this.C(s, p, c));
                    }
                }
                this.add_clause(cs)
            }
            // Do not use forbidden color
            for (const p of range(this.nP)) {
                for (const s of range(this.num_species())) {
                    this.add_clause(
                            [-this.C(s, p, forbidden)]
                    )
                }
            }
        }
    }

    add_constraints_fixed_blank_orientation() {
        for (const p of range(this.nP)) {
            for (const s of range(this.num_species())) {
                this.add_clause([
                    -this.C(s, p, 1), // Either patch p on species s isn't empty
                     this.O(s, p, 0)  // Or the patch is oriented up
                ]);
            }
        }
    }

    add_constraints_no_self_complementarity(above_color=0) {
        for (const c of range(above_color,this.num_colors())) {
            this.add_clause([-this.B(c,c)]);
        }
    }

    public add_constraints_nanoparticle_singletype() {
        // for each location
        for (const l of range(this.nL)){
            // if location l has a nanoparticle
            let hasParticle = l in this.nanoparticles;
            // eval now to save time
            for (const r of range(this.nR)){
                for (const s of range(this.num_species())){
                    if (hasParticle){
                        this.add_clause([this.N_single(s), -this.P(l, s, r)])
                    } else {
                        this.add_clause([-this.N_single(s), -this.P(l, s, r)])
                    }
                }
            }
        }
    }

    public add_constraints_nanoparticle_multitype() {
        // generate disjointitiy clause - each species s has exactly one nanoparticle type
        for (const s of range(this.num_species())) {
            this.add_clause(...this.exactly_one([...range(this.num_nanoparticle_types())].map((t)=>this.N_multi(s, t))));
        }

        // occupancy clause - each position must have its correct nanoparticle type
        for (const l of range(this.nL)) {
            // if location l has a nanoparticle
            let nptype = this.nanoparticles[l] // eval now to save time
            for (const r of range(this.nR)) {
                for (const s of range(this.num_species())){
                    this.add_clause([this.N_multi(s, nptype), -this.P(l, s, r)]);
                }
            }
        }
    }

    add_constraints_no_self_interact_patches() {
        for (const s of range(this.num_species())){
            // for each pair of colors
            for (const [c1, c2] of combinations([...range(this.num_colors())], 2)){
                if (c1 == c2){
                    continue; // colors should never be self-interacting so we can safely ignore
                }
                // either the colors do not interact
                let var_no_inter = -this.B(c1, c2);
                for (const p of range(this.nP)){
                    // this clause was written by chatGPT so I'm only like 70% sure I trust it
                    // The clause says: either var_no_interact is True (colors do not interact)
                    // or at least one patch is not c1 or at least one patch is not c2 for the species
                    this.add_clause([var_no_inter, -this.C(s, p, c1), -this.C(s, p, c2)]);
                }
            }
        }
    }

    fix_position_species(l: number, s: number) {
        // fix position l to be occupied by species s
        this.add_clause([...range(this.nR)].map(r=>this.P(l,s,r)))
    }

    fix_slot_colors(loc: number, sid: number, cid: number) {
        this.add_clause([this.F(loc, sid, cid)]);
    }

    fix_color_interaction(c1: number, c2: number) {
        this.add_clause([this.B(c1,c2)])
    }

    run_solve() : string | {status: SystemBehaviorType} | {status: SystemBehaviorType, rule: string} {
        let tries = 0;
        let status : SystemBehaviorType;
        let decRule: string;
        while (tries < this.maxAltTries) {
            if (minSol && minSol[0]+minSol[1] < this.num_species() + (this.num_colors() / 2)) {
                return 'skipped';
            }
            // try to solve the thing
            let result: Logic.Solution | false = this.solver.solve();
            // if the thing is solved
            if (result) {
                let rule = this.readSolution(result);
                rule.sort((a,b)=>{return patchCount(b)-patchCount(a)});
                decRule = ruleToDec(rule); // translate rule to dec
                status = getSysBehavior(rule);
                console.assert(false, "Gotta write bounded and deterministic code!!!!!")
                if (status.short == "UNB" || status.short == "NDB") {
                    // it did not work! let's not try this again
                    console.log(`https://akodiat.github.io/polycubes/?decRule=${decRule}"`);
                    this.solver.forbid(result);
                }
                else {
                    // it worked! we can stop now.
                    minSol = [this.num_species(), this.num_colors()];
                    return {'status': status, 'rule': decRule};
                }
            } else {
                return {'status': status};
            }
        }
        return {'status': status, 'rule': decRule};
    }

    public readSolution(sol: Logic.Solution) : PolycubeCubeType[] {
        let colorCounter = 1;
        let colorMap = new Map();
        let solution_true_vars: string[] = sol.getTrueVars();
        let ruleMap : Map<number, Map<number, {color, orientation}>> = new Map();
        let bMatches = solution_true_vars.filter(/B\((\d+),(\d+)\)/g.test);
        for (const m of bMatches) {  // color c1 binds with c2
            let b1 = Number(m[1]);
            let b2 = Number(m[2]);
            // console.log("Color {} binds with {}".format(c1, c2))
            console.assert(!colorMap.has(b1) || !colorMap.has(b2));
            if (b1 < 2 || b2 < 2) {
                colorMap.set(b1, 0);
                colorMap.set(b2, 0);
            } else {
                colorMap.set(b1, colorCounter);
                colorMap.set(b2, -colorCounter);
                colorCounter += 1;
            }
        }
        let cMatches = solution_true_vars.filter(/C\((\d+),(\d+),(\d+)\)/g.test);
        for (const m of cMatches) {  // Patch p on species s has color c
            let s = Number(m[1]);
            let p = Number(m[2]);
            let c = Number(m[3]);
            //console.log("Patch {} on species {} has color {}".format(p, s, c))
            if (!ruleMap.has(s)) {
                ruleMap.set(s, new Map());
            }
            if (!ruleMap.get(s).has(p)) {
                ruleMap.get(s).set(p, {color: null, orientation: null});
            }
            ruleMap.get(s).get(p).color = colorMap.get(c);
        }
        let hasOrientation = false;
        let oMatches = solution_true_vars.filter(/O\((\d+),(\d+),(\d+)\)/g.test);
        for (const m of oMatches) {  // Patch on species l has orientation o
            let s = Number(m[1]);
            let p = Number(m[2]);
            let o = Number(m[3]);
            //console.log("Patch {} on species {} has orientation {}".format(p, s, o))
            hasOrientation = true;
            ruleMap.get(s).get(p).orientation = o;
        }
        if (!hasOrientation) {
            console.log("Found no orientation values")
            for (const patches of ruleMap.values()) {
                for (const [i, p] of patches) {
                    p.orientation = getFlatFaceRot()[Number(i)];
                }
            }
        }
        return [...ruleMap.values()].map((ct, i)=> {
            return new DynamicPolycubeCubeType(i, [...ct.values()].map(patch => {
                return new PolycubePatch(patch.color, undefined, patch.orientation);
            }))
        });
    }
}

////// Helper functions:
function getFlatFaceRot() : number[] {
    return [1,1,2,0,0,0]
}

function patchRotToVec(i: number, rot: number) : Vector3{
    /** Get vector indicating patch rotation, given rotation state and face index

    Args:
        i (Math.floor) { Index of patch [0...5]
        rot (Math.floor) { Rotation state [0,1,2,3] = North , East, South, West

    Returns:
        vector: Patch rotation vector
    */

    let v = FACE_ROTATIONS[i].clone();
    let axis = RULE_ORDER[i].clone();
    let angle = rot * Math.PI/2;
    let q = new THREE.Quaternion().setFromAxisAngle(axis, angle);
    return v.clone().applyQuaternion(q).round();
}

function patchVecToRot(i: number, v: Vector3) : number{
    /** Get rotation state, given face index and patch rotation vector

    Args:
        i (Math.floor) { Index of patch [0...5]
        v (vector) { Patch rotation vector

    Returns:
        Math.floor: Rotation state [0,1,2,3] = North , East, South, West
    */
    let angle = getSignedAngle(
        FACE_ROTATIONS[i],
        v,
        RULE_ORDER[i]
    )
    return Math.floor((angle * (2/Math.PI)+4) % 4)
}


function enumerateRotations(dim=3) : Map<number, Map<number, number>> {
    if (dim==2) {
        return new Map([
            [0, new Map([[0, 0], [1, 1], [2, 2], [3, 3], [4, 4], [5, 5]])],
            [1, new Map([[0, 1], [1, 0], [2, 3], [3, 2], [4, 4], [5, 5]])],
            [2, new Map([[0, 2], [1, 3], [2, 1], [3, 0], [4, 4], [5, 5]])],
            [3, new Map([[0, 3], [1, 2], [2, 0], [3, 1], [4, 4], [5, 5]])]
        ]);
    } else {
        return new Map([
            [0,  new Map([[0, 0], [1, 1], [2, 2], [3, 3], [4, 4], [5, 5]])],
            [1,  new Map([[0, 0], [1, 1], [2, 3], [3, 2], [4, 5], [5, 4]])],
            [2,  new Map([[0, 0], [1, 1], [2, 4], [3, 5], [4, 3], [5, 2]])],
            [3,  new Map([[0, 0], [1, 1], [2, 5], [3, 4], [4, 2], [5, 3]])],
            [4,  new Map([[0, 1], [1, 0], [2, 2], [3, 3], [4, 5], [5, 4]])],
            [5,  new Map([[0, 1], [1, 0], [2, 3], [3, 2], [4, 4], [5, 5]])],
            [6,  new Map([[0, 1], [1, 0], [2, 4], [3, 5], [4, 2], [5, 3]])],
            [7,  new Map([[0, 1], [1, 0], [2, 5], [3, 4], [4, 3], [5, 2]])],
            [8,  new Map([[0, 2], [1, 3], [2, 0], [3, 1], [4, 5], [5, 4]])],
            [9,  new Map([[0, 2], [1, 3], [2, 1], [3, 0], [4, 4], [5, 5]])],
            [10, new Map([[0, 2], [1, 3], [2, 4], [3, 5], [4, 0], [5, 1]])],
            [11, new Map([[0, 2], [1, 3], [2, 5], [3, 4], [4, 1], [5, 0]])],
            [12, new Map([[0, 3], [1, 2], [2, 0], [3, 1], [4, 4], [5, 5]])],
            [13, new Map([[0, 3], [1, 2], [2, 1], [3, 0], [4, 5], [5, 4]])],
            [14, new Map([[0, 3], [1, 2], [2, 4], [3, 5], [4, 1], [5, 0]])],
            [15, new Map([[0, 3], [1, 2], [2, 5], [3, 4], [4, 0], [5, 1]])],
            [16, new Map([[0, 4], [1, 5], [2, 0], [3, 1], [4, 2], [5, 3]])],
            [17, new Map([[0, 4], [1, 5], [2, 1], [3, 0], [4, 3], [5, 2]])],
            [18, new Map([[0, 4], [1, 5], [2, 2], [3, 3], [4, 1], [5, 0]])],
            [19, new Map([[0, 4], [1, 5], [2, 3], [3, 2], [4, 0], [5, 1]])],
            [20, new Map([[0, 5], [1, 4], [2, 0], [3, 1], [4, 3], [5, 2]])],
            [21, new Map([[0, 5], [1, 4], [2, 1], [3, 0], [4, 2], [5, 3]])],
            [22, new Map([[0, 5], [1, 4], [2, 2], [3, 3], [4, 0], [5, 1]])],
            [23, new Map([[0, 5], [1, 4], [2, 3], [3, 2], [4, 1], [5, 0]])]
        ]);

    }
}

export function make_solver(spec: SolveSpec): polysat {
    let sat = Object.assign(new polysat(), spec);
    sat.init();
    return sat;
}

function patchCount(cube) :  number{
    return cube.filter(face=>face.color!=0).length;
}

let minSol;
