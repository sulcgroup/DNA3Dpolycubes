// three.js imports
import * as THREE from "three";
import { Matrix4, Quaternion } from "three";
// import seedrandom library
import seedrandom from "seedrandom";
import { ruleToDec } from "./rule";
import { allRotations, getSignedAngle, randOrdering, randstr, range, RULE_ORDER, strToVec, vecToStr } from "./utils";
import { AssemblyComponentList, AssemblyMethod, cubeSetColorMap, PiecewiseStage } from "./solve/assembly";
import { Polycube } from "./polycube";
export class PolycubeSystem extends Polycube {
    // important ones!
    // map of string reps of coords -> cube objs
    /// ... or that's what i thought at least; i'm actually not convinced at this point that this does ANYTHING
    // current maximum coordinate value
    // list of keys for possuble moves
    // moveKeys: any[];
    // center of mass of the entire polycube
    // max allowed number of cubes
    nMaxCubes;
    // max allowed coord
    maxCoord;
    maxPieceSize;
    // private pieces: PolycubeSystem[];
    // map showing which colors appear in which pieces
    // private pieceColorMap: {};
    // ASSEMBLY PARAMETERS!!!!!
    // "stochatic" or "seeded"
    // staging order is a list of sets
    // each set is a stage
    // for unstaged assemblies, staging.size == 1
    staging;
    // for piecewise-assembly systems, we need an additional array for piecewise stages
    piecewise_staging;
    // index in `staging` of current stage of assembly
    current_stage;
    assembly_method;
    // step counter
    steps;
    // graphics materials for object colors
    objGroup;
    confMap;
    mismatches_allowed;
    // system-specific rng
    rng_seed;
    rng;
    /**
     @param rules: the rules to use to construct this system
     @param scene: root object for THREE graphics
     @param nMaxCubes: the maximum number of cubes to use to construct this system
     @param maxCoord: the maximum distance along any axis from the center that the system will tolarate
     @param assemblyMode
     'ordered', and 'seeded'
     @param buildConfmap: ?????
     @param torsion: true to make torsion active, false otherwise
     * @param envbrowser
     * @param maxPieceSize
     * @param rand_seed
     */
    constructor(nMaxCubes = 1000, maxCoord = 100, buildConfmap = true, torsion = true, maxPieceSize = 1, rand_seed = "abcd") {
        super();
        /**
        @type THREE.Vector3
        @desc something something graphics? @todo fill in real descriptio
         */
        /**
        @type int
        @desc the maximum number of cubes tat the system will model before it stops adding more cubes
         */
        this.nMaxCubes = nMaxCubes;
        // maximum coordinate possible in this build space
        this.maxCoord = maxCoord;
        // minimum and maximum coordinates of this structure - will be updated
        // as structure assembles
        this.current_stage = 0;
        this.torsion = torsion;
        this.staging = [new AssemblyComponentList(0)];
        this.maxPieceSize = maxPieceSize;
        this.mismatches_allowed = true;
        if (this.maxPieceSize > 1) {
            this.piecewise_staging = [new PiecewiseStage(0)];
        }
        this.steps = 0;
        if (buildConfmap) {
            this.confMap = new Map();
        }
        this.reset_rng(rand_seed);
    }
    /**
     * adds a cube type
     * @param cubeType
     */
    addCubeType(cubeType, update_piecewise = true) {
        super.addCubeType(cubeType);
        if (this.staging.length == 1) {
            this.staging[0].add(cubeType);
        }
        if (this.isPiecewise() && update_piecewise) {
            this.refresh_piecewise_staging();
        }
    }
    /**
     * toggles whether a cube type is allowed to be used
     * for simplicity, we'll treat this for now as adding/removing from
     * @param i
     */
    toggleCubeType(i) {
        let ct = this.cube_types[i];
        // Check if `i` exists in the first set
        if (this.staging[0].in(ct)) {
            // If it exists, remove it
            this.staging[0].removeItem(ct);
        }
        else {
            // If it doesn't exist, add it
            this.staging[0].add(ct);
        }
    }
    /**
     * tests if the provided cube type is allowed at the provided stage of assembly
     */
    cube_type_allowed(cube_type_idx, stage_idx = undefined) {
        if (stage_idx === undefined) {
            stage_idx = this.current_stage;
        }
        return cube_type_idx in this.staging[stage_idx];
    }
    /**
     * returns the number of cube types that are accessable in any stage of assembly
     */
    countAllowedCubeTypes() {
        // Create a set to hold the unique values
        let unionSet = new Set();
        // Iterate over each set in staging
        for (let set of this.staging) {
            // Add each number in the set to the union set
            set.forEach(num => unionSet.add(num));
        }
        // The size of the union set is the number of unique values
        console.assert(unionSet.size <= this.numCubeTypes());
        return unionSet.size;
    }
    allowedCubeTypes(stage = undefined) {
        if (typeof stage == "undefined") {
            stage = this.current_stage;
        }
        return this.staging[stage].copy();
    }
    numStages() {
        return this.staging.length;
    }
    cubeTypeIndex(r) {
        for (let i = 0; i < this.cube_types.length; i++) {
            if (this.cube_types[i].typeName == r.typeName) {
                return i;
            }
        }
        return -1;
    }
    removeCubeType(ct) {
        // i hate this language so much
        if (typeof ct == "number") {
            this.removeCubeType(this.getCubeType(ct));
        }
        else {
            // treat ct as type id
            // loop cube types
            for (let i = 0; i < this.cube_types.length; i++) {
                // we cannot assume that if ct is a PolycubeCube, ct.typeid is the index in cube_types
                if (this.cube_types[i].equals(ct)) {
                    // joke language
                    let new_cube_types = [...this.cube_types.slice(0, i),
                        ...this.cube_types.slice(i + 1, this.cube_types.length)];
                    this.set_rule(new_cube_types);
                    break;
                }
            }
            // gotta remove cube types from staging
            for (let i_stage = 0; i_stage < this.staging.length; i_stage++) {
                for (let i_ct = 0; i_ct < this.staging[i_stage].size(); i_ct++) {
                    let staging_ct = this.staging[i_stage].get(i_ct);
                    if (staging_ct.equals(ct)) {
                        this.staging[i_stage].removeItem(staging_ct);
                        break;
                    }
                }
            }
            // if there are any allosteric cubes
            if (this.cube_types.length > 0) {
                if (this.cube_types.map(ct => ct.hasAllostery()).reduce((a, b) => a || b)) {
                    // refresh piecewise staging
                    this.refresh_piecewise_staging();
                }
            }
        }
    }
    getAssemblyMethod() {
        return this.assembly_method;
    }
    setAssemblyMethod(a) {
        if (typeof a === 'string') {
            switch (a) {
                case "Ordered":
                case "Stepwise":
                    this.assembly_method = AssemblyMethod.StepwiseStaging;
                    break;
                case "Stochastic":
                    this.assembly_method = AssemblyMethod.Stochastic;
                    break;
                // ... handle other string cases if necessary ...
            }
        }
        else {
            this.assembly_method = a;
        }
    }
    setMisMatchesAllowed(bNewVal) {
        this.mismatches_allowed = bNewVal;
    }
    clone() {
        let copy = super.clone();
        copy.nMaxCubes = this.nMaxCubes;
        copy.maxCoord = this.maxCoord;
        if (this.confMap !== undefined) {
            copy.confMap = new Map();
        }
        copy.torsion = this.torsion;
        copy.maxPieceSize = this.maxPieceSize;
        copy.steps = this.steps;
        if (typeof this.piecewise_staging == "undefined") {
            copy.piecewise_staging = undefined;
        }
        else {
            copy.staging = this.staging.map(a => a.copy());
        }
        copy.mismatches_allowed = this.mismatches_allowed;
        return copy;
    }
    isPolycubeSystem() {
        return true;
    }
    isPiecewise() {
        return this.maxPieceSize > 1 && this.cube_types.some(r => r.hasAllostery());
    }
    setMaxPieceSize(n) {
        this.maxPieceSize = n;
    }
    /**
     * Resets staging, so there's only one assembly stage which contains all particles
     */
    clear_staging() {
        this.staging = [new AssemblyComponentList(0, [...this.cube_types])];
    }
    // adds initial "seed" particle
    seed() {
        let i;
        if (this.cube_types.length == 0 || this.countAllowedCubeTypes() == 0) {
            console.error("No cube types!!!");
            return;
        }
        if (this.countAllowedCubeTypes() > 0) {
            if (this.isPiecewise()) {
                // choose random seed piece
                let seed_piece = this.piecewise_staging[0].all_pieces.r(this.rng);
                seed_piece.forEachCube((c, k) => {
                    this.addParticle(c.instantiate(c.position, this.genCubeName()));
                });
                this.steps += seed_piece.numCubes();
                console.assert((this.numMoves() == 0) == (this.numMovePositions() == 0));
            }
            else {
                let seed_cube_type = this.staging[0].r(this.rng);
                let seed_cube = this.instantiateCube(seed_cube_type, this.genCubeName());
                let rotidx = Math.floor(allRotations().length * this.rng()); // rotate seed randomly
                let mat4 = new Matrix4();
                mat4.setFromMatrix3(allRotations()[rotidx]);
                let q = new Quaternion();
                q.setFromRotationMatrix(mat4);
                seed_cube.rotate(q);
                this.addParticle(seed_cube);
                this.steps++;
                console.assert((this.numMoves() == 0) == (this.numMovePositions() == 0));
            }
        }
        else {
            //this will likely occur mainly if not only if there's no rules
            console.log("No rule!");
        }
    }
    reset_rng(newseed) {
        if (newseed != void 0) {
            this.rng_seed = newseed;
        }
        this.rng = seedrandom(this.rng_seed);
    }
    reset(reset_random = true) {
        if (reset_random) {
            this.reset_rng();
        }
        this.steps = 0;
        this.current_stage = 0;
        this.cubes = new Map();
        this.moves = new Map();
        if (this.confMap) {
            this.confMap = new Map;
        }
        this.matches = 0;
        this.mismatches = 0;
        this.connections = [];
    }
    regenerate() {
        this.reset();
        this.seed();
        //this.processMoves();
    }
    set_rule(rules = [], reset_random = true) {
        // if no change has been made to rule
        if (ruleToDec(rules) == ruleToDec(this.cube_types)) {
            return;
        }
        super.set_rule(rules);
        // reset staged-assembly pieces
        if (this.isPiecewise()) {
            this.refresh_piecewise_staging();
        }
    }
    /**
     * refreshes staging for piecewise assembly
     */
    refresh_piecewise_staging() {
        // construct piecewise staging list
        this.piecewise_staging = this.staging.map(s => new PiecewiseStage(s.index));
        // iter stages
        for (let iStage = 0; iStage < this.numStages(); iStage++) {
            // construct mapping of color -> pieces for this stage
            // let stageColorMap = new Map<number, AssemblyComponentList<PolycubeSystem>>();
            let pieces = [];
            // loop allowed cube types
            for (let ct of this.allowedCubeTypes(iStage)) {
                // add single-cube pieces
                let p = new PolycubeSystem(null, 1, false, true, 1, randstr(8, this.rng));
                p.set_rule([ct]);
                p.seed();
                pieces.push(p);
            }
            // i've replaced the previous algorithm for piecewise assembly on account of what it's indefensibly awful
            // better algorithm: recursive breadth-first tree search
            console.assert(this.maxPieceSize > 1);
            pieces = this.iterPiecewiseAssembly(pieces, this.allowedCubeTypes(iStage).options(), this.maxPieceSize);
            this.addPieces(iStage, ...pieces);
        }
        this.reset_rng();
    }
    numSteps() {
        return this.steps;
    }
    addPiece(p, stage = 0) {
        this.piecewise_staging[stage].add(p);
    }
    addPieces(stage = 0, ...pieces) {
        this.piecewise_staging[stage].addAll(pieces);
    }
    /**
     * recursive piece-generation function
     * @param pieces
     * @param cubeOptions
     * @param maxDepth
     */
    iterPiecewiseAssembly(pieces, cubeOptions, maxDepth) {
        console.assert(pieces.length > 0);
        console.assert(cubeOptions.length > 0);
        let cubeColors = cubeSetColorMap(cubeOptions);
        // list for next lvl of iteration
        let new_pieces = [];
        // for each polycube we provide on this level
        pieces.map(piece => {
            // for each move
            piece.allMoves().map(m => {
                let matching_colors = this.matchingColors(m.color);
                // if no options
                if (!matching_colors.some(c => cubeColors.has(c))) {
                    return [];
                }
                // iter colors compatible with move
                matching_colors.forEach(c => {
                    cubeColors.get(c).forEach(ct => {
                        // iter patches on cube
                        let cubes = [];
                        // construct all the ways we can put a cube of type ct on the polycube
                        // at move m
                        ct.patches.forEach((patch, i) => {
                            // skip no color patches
                            if (!patch.color) {
                                return;
                            }
                            // if this is one of the patches compatible with move
                            if (this.compatibleColors(m.color, patch.color)) {
                                // Rotate rule so that the matching face has
                                // the same direction:
                                let cube = this.instantiateCube(ct, this.genCubeName(), m.pos).rotateInstanceFromTo(RULE_ORDER[patch.direction], m.direction);
                                // get direction index
                                let i = RULE_ORDER.indexOf(m.direction);
                                // if torsion is off, match patch align directions
                                if (this.torsion) {
                                    let theta = getSignedAngle(m.alignDir, cube.patch_by_local(i).alignDir, RULE_ORDER[i]);
                                    // apply rotatiion
                                    cube.rotateInstanceAroundAxis(m.direction, theta);
                                    console.assert(this.compatibleColors(m.color, cube.patch(m.direction).color), "Incompatible colors");
                                    console.assert(!this.torsion || m.alignDir.distanceTo(cube.patch_by_local(i).alignDir) < 1e-5, `Angle mismatch between ${vecToStr(m.alignDir)} and ${vecToStr(cube.patch_by_local(i).alignDir)}`);
                                    // final checks for stuff
                                    // console.assert(cube.patches[RULE_ORDER.indexOf(move.direction)].color == -move.color)
                                    if (this.canPlaceCube(cube)) {
                                        console.assert(this.compatibleColors(m.color, cube.patch_by_local(i).color));
                                        // Return the rotated rule b
                                        cubes.push(cube);
                                    }
                                }
                                else {
                                    // iter possible directions
                                    [...Array(4).keys()].forEach(rotIdx => {
                                        let theta = Math.floor(rotIdx) * Math.PI / 2;
                                        // apply rotatiion
                                        cube.rotateInstanceAroundAxis(m.direction, theta);
                                        console.assert(!this.torsion || m.alignDir.distanceTo(cube.patch_by_local(i).alignDir) < 1e-5, `Angle mismatch between ${m.alignDir} and ${cube.patch_by_local(i).alignDir}`);
                                        cube.position = m.pos;
                                        // final checks for stuff
                                        // console.assert(cube.patches[RULE_ORDER.indexOf(move.direction)].color == -move.color)
                                        if (this.canPlaceCube(cube)) {
                                            console.assert(this.compatibleColors(m.color, cube.patch_by_local(i).color));
                                            // Return the rotated rule b
                                            cubes.push(cube);
                                        }
                                    });
                                }
                            }
                        });
                        // iter cubes always will be 1-4 copies of the same cube
                        cubes.forEach(cube => {
                            // clone piece
                            let p2 = piece.clone();
                            p2.addParticle(cube);
                            // if new piece doesn't have an equivalent in mew_pieces
                            if (!new_pieces.some(p2.isEquivelent)) {
                                new_pieces.push(p2);
                            }
                        });
                    });
                });
            });
        });
        let lvl = Math.max(...pieces.map(p => p.numCubes()));
        if (lvl >= maxDepth) {
            return new_pieces;
        }
        else
            return [...new_pieces, ...pieces];
    }
    setMaxCubes(n) {
        this.nMaxCubes = n;
    }
    /**
     * todo: complexity!
     * @param c
     */
    matchingColors(c) {
        return [-c];
    }
    countColors() {
        return this.listColors().length;
    }
    listColors(absvals = true) {
        let allRules = [].concat.apply([], this.cube_types.map(x => {
            return x.patch_colors();
        }));
        return [...new Set(allRules)];
    }
    /**
    Finds the first name of format R[number] that is not in use, and returns it
     */
    nameNewCubeType() {
        let i = 1;
        while (true) { //BAD PRACTICE :D
            if (!this.cube_types.some(r => r.typeName == `CT${i}`)) {
                return `CT${i}`;
            }
            i++;
        }
    }
    /**
    @param move the move spot to check
    @param cube the rule to check
    Checks if a rule fits in a given move
     */
    fitCube(move, cube) {
        console.assert(!this.hasCube(move.pos), "Already a cube here!!!");
        // Traverse rule and move faces in random order
        let rb = randOrdering(RULE_ORDER.length, this.rng);
        // Check each face in cube type
        // not great stochastically, gonna live with it
        for (let j of rb) {
            // skip non-patches
            if (!cube.patch_by_local(j).color) {
                continue;
            }
            // skip inactive patches
            if (!cube.patch_is_active(j)) {
                continue;
            }
            console.assert(j == cube.patch_by_local(j).direction, "Invalid patch somehow");
            // If we find an equal color
            if (this.compatibleColors(move.color, cube.patch_by_local(j).color)) {
                // "fit" the cube by rotating it so the patch with teh matching color is oriented in the
                // same direction as the move
                let fit_cube = this.fitCubeFace(move, cube, j);
                if (fit_cube) { // check for non-fits (usually from mismatch problems)
                    return fit_cube;
                }
            }
        }
        // Return false if we didn't find any matching faces
        return false;
    }
    fitCubeFace(move, cube, face) {
        // Rotate cube so that the matching face has
        // the same direction:
        cube = cube.clone();
        console.assert(this.compatibleColors(move.color, cube.patch(face).color), `Move color ${move.color} does not match patch color ${cube.patch(move.direction).color}`);
        cube.rotateInstanceFromTo(RULE_ORDER[face], move.direction);
        console.assert(move.direction.distanceTo(RULE_ORDER[cube.patch(move.direction).direction]) < 1e-5);
        console.assert(this.compatibleColors(move.color, cube.patch(move.direction).color), `Move color ${move.color} does not match patch color ${cube.patch(move.direction).color}`);
        let theta;
        // if torsion is off, match patch align directions
        let align_rot_axis = move.direction;
        if (this.torsion) {
            var patch_curr_align = cube.patch(move.direction).alignDir;
            theta = getSignedAngle(patch_curr_align, move.alignDir, align_rot_axis);
        }
        else {
            // otherwise randomize rotation
            theta = Math.floor(this.rng() * 4) * Math.PI / 2;
        }
        // apply rotation
        console.assert(RULE_ORDER[cube.patch(move.direction).direction] == move.direction);
        cube = cube.rotateInstanceAroundAxis(align_rot_axis, theta);
        console.assert(this.compatibleColors(move.color, cube.patch(move.direction).color), `Move color ${move.color} does not match patch color ${cube.patch(move.direction).color}`);
        console.assert(RULE_ORDER[cube.patch(move.direction).direction] == move.direction);
        let patch_new_align = cube.patch(move.direction).alignDir;
        console.assert(!this.torsion || move.alignDir.distanceTo(patch_new_align) < 1e-5, `Angle mismatch between ${vecToStr(move.alignDir)} and ${vecToStr(patch_new_align)}`);
        cube.position = move.pos;
        // final checks for stuff
        // console.assert(cube.patches[RULE_ORDER.indexOf(move.direction)].color == -move.color)
        if (this.canPlaceCube(cube)) {
            console.assert(this.compatibleColors(move.color, cube.patch(move.direction).color));
            // Return the rotated rule b
            return cube;
        }
        else {
            return false;
        }
    }
    /**
     * Determines if a cube fits in a location
     * @param c a cube, already rotated and translated
     */
    canPlaceCube(c) {
        let [connectionCount, mismatchCount] = this.countConnections(c);
        return connectionCount > 0 && (this.mismatches_allowed || mismatchCount == 0);
    }
    /**
    @param move the move spot to check
    @param piece the piece to check
    Checks if a piece fits in a given move
    @todo: optimize this by pre-mapping possible pieces for any given color
     */
    pieceFits(move, piece) {
        // Traverse move faces and piece faces in random order
        let piece_moves = piece.allMoves();
        let rc = randOrdering(piece.moves.size, this.rng);
        // for each move on the piece
        for (let ric = 0; ric < piece.moves.size; ric++) {
            let j = rc[ric];
            let piece_move = piece.moveFromIdx(j);
            let valid_piece = this.testPieceMove(move, piece, piece_move);
            if (valid_piece) {
                return valid_piece;
            }
        }
        // Return false if we didn't find any matching faces
        return false;
    }
    /**
    @param move
    @param piece
    @param i the index of the face of the move (in this) being tested
    @param j the index of the face of the piece move being tested
    @param k the index of the piece move being tested
     */
    testPieceMove(move, piece, piece_move) {
        // If we find an equal color
        if (!this.compatibleColors(move.color, piece_move.color)) {
            return false;
        }
        // first, find the location of piece_move in the global space
        let global_translation = move.pos.clone();
        //            let piece_translated = piece.translate(t_initial);
        //reposition the piece so that the cube that's the parent of piece_move[j]
        // is at 0,0,0
        let piece_translation = piece_move.pos.clone().add(piece_move.direction).negate();
        let piece_t = piece.translate(piece_translation);
        // Rotate the piece so that the matching face has
        // the same direction:
        let q = new THREE.Quaternion(); // create one and reuse it
        q.setFromUnitVectors(piece_move.direction, move.direction.clone().negate());
        let t_rot = new THREE.Matrix4().makeRotationFromQuaternion(q);
        let piece_rotated = piece_t.transform(t_rot);
        let rotated_move_key = piece_move.pos.clone().add(piece_translation).applyQuaternion(q).round(); //however move keys will of course be different
        let piece_move_rotated = piece_rotated.getMoveCell(rotated_move_key);
        //now, reposition piece_rotated so that it lines up in the global space
        // moves should index the same way after applying a transformation - key word "should"
        // reindex j. use applyQuaternion(q) to avoid translating face vector
        let j_prime = RULE_ORDER.indexOf(piece_move.direction.clone().applyQuaternion(q).round());
        //            console.assert(rotated_move_key.clone().applyMatrix4(t).equals(move.pos));
        console.assert(piece_move_rotated != void 0);
        console.assert(this.compatibleColors(move.color, piece_move_rotated.speciesFit[j_prime].color));
        let theta;
        if (this.torsion) {
            theta = -getSignedAngle(move.alignDir, piece_move_rotated.speciesFit[j_prime].alignDir, move.direction);
            //                console.assert(move[i].alignDir.distanceTo(rule_rotated.alignments[i]) < 1e-5);
        }
        else {
            // ...and a random rotation:
            theta = Math.floor(this.rng() * 4) * Math.PI / 2;
        }
        if (theta) {
            q = new THREE.Quaternion().setFromAxisAngle(move.direction, theta);
            piece_rotated = piece_rotated.rotate(q);
            //            console.assert(this.compatibleColors(move[i].color, rule_rotated.colors[i]));
            // Return the rotated rule b
        }
        // check for overlaps
        piece_rotated = piece_rotated.translate(global_translation);
        if ([...piece_rotated.cubes.keys()].some((cube_position) => this.hasCube(cube_position))) {
            return false;
        }
        return piece_rotated;
    }
    /**
     @param move
     @param cubeType the index of the rule to test at this location
     */
    tryProcessSingleCubeMove(move, cubeType) {
        let cube = this.instantiateCube(cubeType, this.genCubeName(), move.pos.clone());
        // let cube: PolycubeCubeType | false = cubeType.instantiate(move.pos.clone(), this.genCubeName());
        // console.assert(this.hasMove(movekey))
        cube = this.fitCube(move, cube);
        if (cube) { // if a rule is identified that fits the move
            // loop faces of the identified cube type
            // for (let i=0; i < RULE_ORDER.length; i++) {
            // 	// identify corresponding face in move object
            //     let neigb: null | {color, alignDir} = this.getMoveCell(move.pos).speciesFit[i]
            //     if (neigb != null && cube.patch_by_local(i).color) {
            // 		// Josh Note: I think this next line of code is a bug. If we're checking for matches,
            // 		// it should be this.compatibleColors(neighb.color, rule[i].color)
            //         // if (neigb.color == rule[i].color && neigb.alignDir.equals(rule[i].alignDir)) {
            // 		let color_match = this.compatibleColors(neigb.color, cube.patch_by_local(i).color)
            // 		let compatible_align = neigb.alignDir.distanceTo(cube.patch_by_local(i).alignDir) < 1e-5;
            // 		if (color_match && compatible_align){
            //             this.matches++;
            //             let parent_position: Vector3 = cube.position.clone();
            // 			parent_position.add(RULE_ORDER[i]);
            // 			console.assert(parent_position instanceof Vector3)
            //             this.addConnection(cube, this.getCube(parent_position) , neigb.color);
            //
            //         } else {
            //             this.mismatches++;
            //         }
            //     }
            // }
            this.addParticle(cube);
            return cube;
        }
        return false;
    }
    tryProcessPieceMove(move, piece) {
        let fitted_piece = this.pieceFits(move, piece);
        if (fitted_piece) { // if a rule is identified that fits the move
            [...fitted_piece.cubes.values()].forEach((c) => {
                let cube = this.instantiateCube(c.getType(), this.genCubeName(), c.position);
                cube.rotate(c.rotation);
                this.addParticle(cube);
                if (this.hasMove(c.position)) {
                    for (let i = 0; i < RULE_ORDER.length; i++) {
                        let neigb = this.getMoveCell(c.position).speciesFit[i];
                        if (neigb != null) {
                            if (this.compatibleColors(neigb.color, cube.patch_from_cube_type(i).color) && neigb.alignDir.equals(cube.patch_from_cube_type(i).alignDir)) {
                                this.matches++;
                                let parent_position = cube.position.clone().add(RULE_ORDER[i]);
                                this.addConnection(cube, this.getCube(parent_position), neigb.color);
                            }
                            else {
                                this.mismatches++;
                            }
                        }
                    }
                    // Remove processed move
                    this.removeMove(c.position);
                }
            });
        }
        return fitted_piece;
    }
    getValidPieces(move, stage) {
        if (stage == void 0) {
            stage = this.current_stage;
        }
        console.assert(move != void 0);
        let validPieces = new AssemblyComponentList(stage);
        if (move.color in this.piecewise_staging[stage].colormap) {
            this.piecewise_staging[stage].colormap[move.color].forEach(p => {
                if (!validPieces.in(p)) {
                    validPieces.add(p);
                }
            });
        }
        return validPieces;
    }
    run() {
        let hasNext;
        do {
            hasNext = this.step();
        } while (hasNext && this.numCubes() <= this.nMaxCubes);
    }
    step() {
        if (this.numMovePositions() == 0) {
            return false;
        }
        console.assert(this.numMoves() > 0);
        let result = false; //gotta deal with scoping problems
        // Pick a random move
        let n_moves = this.numMoves();
        let move = this.moveFromIdx(Math.floor(this.rng() * n_moves));
        console.assert(move != null);
        if (!this.isPiecewise()) {
            // list allowed cube types
            let cube_types = this.allowedCubeTypes();
            // shuffle list
            cube_types = cube_types.shuffled(this.rng);
            // Check if we have a rule that fits this move
            result = false;
            while (cube_types.size() > 0) {
                let ct = cube_types.pop();
                let cube = this.tryProcessSingleCubeMove(move, ct);
                if (cube) {
                    result = true;
                    this.steps++;
                    break;
                }
            }
            // if we didnd't successfully place a cube, remove the face of this move
            if (!result) {
                this.removeMove(move);
            }
            else {
                // if we did place a cube, the entire cell is now occupied, so remove the entire move
                this.removeMove(move.pos);
            }
        }
        else {
            let validPieces = this.getValidPieces(move);
            validPieces = validPieces.shuffled(this.rng);
            // let pieceIdxs = randOrdering(validPieces.length);
            // Check if we have a rule that fits this move
            result = false;
            while (validPieces.size() > 0) {
                let piece = validPieces.pop();
                let placed_piece = this.tryProcessPieceMove(move, piece);
                if (placed_piece) {
                    result = true;
                    this.steps += placed_piece.numCubes();
                    this.removeMove(move);
                    break;
                }
            }
            // move keys will be automatically removed in piecewise assembly
            if (this.hasMove(move)) {
                // Remove processed move
                this.removeMove(move);
            }
        }
        console.assert((this.numMoves() == 0) == (this.numMovePositions() == 0));
    }
    isProcessingComplete() {
        return this.moves.size == 0 || this.cubes.size > this.nMaxCubes;
    }
    // randomCubeType(stage?: number | undefined) : PolycubeCube {
    // 	if (typeof stage == "undefined"){
    // 		stage = this.current_stage;
    // 	}
    // 	return this.staging[stage].r(this.rng);
    // }
    //
    // getCubeObjects() : Mesh[] {
    // 	return Object.entries(this.cubes).map((key, value: any) => {
    // 		console.assert(! (value.cube_vis instanceof undefined))
    // 		return value.cube_vis;
    // 	});
    // }
    resetMoves() {
        this.moves.clear();
        this.cubes.forEach((_, move_key) => this.updateMoves(this.cubes.get(move_key)));
    }
    /**
     * lists possible orders in which we can add cube types
     */
    getAllPossibleMoveOrders(stage) {
        if (stage == void 0) {
            stage = 0;
        }
        let cts = this.staging[stage];
        let num_cts = cts.size();
        let num_assembly_orders = num_cts ** this.maxPieceSize;
        //probably a better algorithm for this but i am too lazy to find it
        let assembly_order = [...range(num_assembly_orders)].map(() => [...range(this.maxPieceSize)].map(() => new AssemblyComponentList(stage)));
        // iter location indexes
        for (let iPos = 0; iPos < this.maxPieceSize; iPos++) {
            let ct_idx = 0; // can NOT assume this indexes this.cube_types
            // iter assembly orders
            for (let i = 0; i < Math.pow(this.countAllowedCubeTypes(), this.maxPieceSize); i++) {
                console.assert(assembly_order[i][iPos].size() == 0);
                assembly_order[i][iPos].add(cts.get(ct_idx));
                if (i % Math.pow(num_cts, iPos) == 0) {
                    ct_idx++;
                }
                if (ct_idx == this.countAllowedCubeTypes()) {
                    ct_idx = 0;
                }
            }
        }
        return assembly_order;
    }
    transform(t) {
        let copy = this.clone();
        console.assert(copy.cubes.size == this.cubes.size);
        copy.reset();
        [...this.cubes.values()].forEach((cube, _) => {
            let new_position = cube.position.clone().applyMatrix4(t).round();
            let q = new THREE.Quaternion().setFromRotationMatrix(t);
            cube.rotate(q);
            copy.cubes.set(vecToStr(new_position), cube);
            copy.cubes.get(vecToStr(new_position)).position = new_position.clone();
        });
        copy.moves = new Map();
        [...this.moves.keys()].forEach((v, i) => {
            let key = strToVec(v);
            key = key.applyMatrix4(t).round();
            copy.moves.set(vecToStr(key), this.getMoveCell(v).transform(t));
        });
        console.assert(copy.numCubes() == this.numCubes());
        return copy;
    }
    rotate(rot) {
        if (rot instanceof THREE.Quaternion) {
            return this.transform(new THREE.Matrix4().makeRotationFromQuaternion(rot));
        }
        else if (rot instanceof THREE.Matrix3) {
            return this.transform(new THREE.Matrix4().setFromMatrix3(rot));
        }
        else {
            return this.transform(rot);
        }
    }
    translate(transl) {
        return this.transform(new THREE.Matrix4().setPosition(transl));
    }
    getCurrentTopology() {
        let bindings = [];
        let empty = [];
        let donePairs = []; // Keep track so that only one bond per pair is saved
        let coords = [...this.confMap.keys()];
        // For each position
        coords.forEach((key, i) => {
            let current = this.cubes[key];
            current.connections.forEach((conn, dPi) => {
                if (!conn) {
                    return;
                }
                let other;
                if (conn.cube_1 == current) {
                    other = conn.cube_2;
                }
                else {
                    other = conn.cube_1;
                }
                let j = coords.findIndex(c => c == other.position);
                if (!donePairs.includes([i, j].sort().toString())) {
                    // add the connection of i and j by dP, dP * -1
                    bindings.push([
                        // Particle {} patch {} 
                        i, dPi,
                        // with Particle {} patch {}
                        j, dPi + (dPi % 2 == 0 ? 1 : -1) // equivelant of RULE_ORDER.idxOf(dP * -1)
                    ]);
                    //console.log(`Particle ${i} patch ${dPi} with particle ${j} patch ${dPi + (dPi % 2 == 0 ? 1 : -1)}`);
                    donePairs.push([i, j].sort().toString());
                }
                else {
                    // if there is no connection between coords[i] and coords[i] + dP...
                    // add (i, dPi) to list of non-connections
                    //... right?
                    empty.push([i, dPi]);
                }
            });
        });
        return [bindings, empty];
    }
    genCubeName() {
        let cube_name;
        do {
            cube_name = randstr(4, this.rng);
        } while (this.cube_name_map.has(cube_name));
        return cube_name;
    }
    /**
     * Exports a polycube as a json file
     */
    export() {
        let json = {};
        json['cube_types'] = this.cube_types;
        json['cubes'] = [...this.cubes.entries()].map(([key, c]) => {
            return {
                "position": {
                    "x": c.position.x,
                    "y": c.position.y,
                    "z": c.position.z
                },
                "rotation": {
                    "w": c.rotation.w,
                    "x": c.rotation.x,
                    "y": c.rotation.y,
                    "z": c.rotation.z
                },
                "type": c.type_id(),
                "personal_name": c.getPersonalName(),
                "state": c.getState()
            };
        });
        json['tstep'] = this.steps;
        return JSON.stringify(json, null, 4);
    }
    updateMove(cube, face) {
        let movePos = cube.position.clone().add(RULE_ORDER[face.direction]);
        if (Math.abs(movePos.x) > this.maxCoord ||
            Math.abs(movePos.y) > this.maxCoord ||
            Math.abs(movePos.z) > this.maxCoord) {
            // Neigbour outside of bounding box, stopping here
            return;
        }
        super.updateMove(cube, face);
    }
}
//# sourceMappingURL=polycubeSystem.js.map