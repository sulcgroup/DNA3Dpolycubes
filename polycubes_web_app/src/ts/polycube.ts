import {Move, PolycubeCellMove} from "./move";
import {PolycubeCubeType, ruleToDec} from "./rule";
import * as THREE from "three";
import {Quaternion, Vector3} from "three";
import {allRotations, RULE_ORDER, strToVec, vecToStr} from "./utils";
import {DefaultCubeInst, PatchWrapper, PolycubeCubeInst, PolycubeSystemConnection} from "./cube_instance";
import {PolycubePatch} from "./patch";


// base-class for polycube. exclude piecewise and staging info
export class Polycube {
    // list of all possible moves
    protected moves: Map<string, PolycubeCellMove>;
    // list of keys for moves, which are depleted as the program checks possible moves, and
    // reset at a new step
    // protected moveKeys: Vector3[];
    // cubes
    protected cubes: Map<string, DefaultCubeInst>;
    protected cube_name_map: Map<string, DefaultCubeInst>; // maps cube personal names to cube instances


    protected currMinCoord: number;
    protected currMaxCoord: number;

    cube_types: PolycubeCubeType[];
    protected cubeTypeCount: Map<number, number>
    protected connections: PolycubeSystemConnection[];
    // torsion or no?
    protected torsion: boolean;
    // number of correct color matches
    protected matches: number;
    // number of incorrect color matches
    protected mismatches: number;

    constructor() {
        /**
         @type dict
         @desc a dict where keys are vector strings and values are move objects
         */
        this.moves = new Map<string, PolycubeCellMove>();
        this.cubes = new Map<string, PolycubeCubeInst<{patch: PolycubePatch, srcdiridx: number}>>();
        this.cube_types = [];
        this.cubeTypeCount = new Map<number, number>();


        this.cube_name_map = new Map()


        this.currMinCoord = 0;
        this.currMaxCoord = 0;

        // list of connections between two blocks
        this.connections = [];

        // semi-deprectated count of successful connections
        this.matches = 0;

        // semi-deprecated count of unsuccessful connections
        this.mismatches = 0;
    }

    clone() {
        let copy = new Polycube();
        copy.set_rule(this.cube_types);

        [...this.cubes.values()].forEach((v) => {
            copy.cubes.set(vecToStr(v.position), v.clone());
        });
        [...this.moves.keys()].forEach((v, i) => {
            copy.moves.set(v, this.getMoveCell(v).clone());
        });
        return copy;
    }

    set_rule(rules : PolycubeCubeType[] = []) {
        this.cubes = new Map();
        this.cube_types = [];
        rules.forEach(r=>{this.addCubeType(r)});
        // todo: check for existing cubes?
        // console.assert(this.numCubes() == 0);
        this.cubeTypeCount = new Map();
        this.cube_types.forEach(ct=>{
            this.cubeTypeCount.set(ct.type_id, 0);
        });
    }

    public hasTorsion() : boolean{
        return this.torsion;
    }

    public numMismatches() : number {
        return this.mismatches;
    }


    /**
     * todo: better
     * @param c1
     * @param c2
     */
    compatibleColors(c1: number, c2: number) :boolean{
        return c1 == -c2;
    }

    public addCubeType(cubeType: PolycubeCubeType) {
        this.cube_types.push(cubeType);
    }

    /**
     * applies the function f to each move in this polycube
     * @param f function to apply to moves
     */
    public forEachMove(f: (PolycubeMove, string) => void) {
        this.moves.forEach(f);
    }

    /**
     * @return the number of move positions
     */
    public numMovePositions() : number {
        return this.moves.size;
    }

    /**
     * accessor for cube types
     * @param i cube type index
     * @return cube type with index i
     */
    public getCubeType(i: number) : PolycubeCubeType{
        return this.cube_types[i];
    }

    /**
     * why does this pass by value???
     * @return a list of the cube types in self (pass by value)
     */
    public listCubeTypes() : PolycubeCubeType[]{
        return this.cube_types.map((ct: PolycubeCubeType) => ct.clone(ct.type_id));
    }

    /**
     *
     */
    getCubeCoords() : Vector3[] {
        return [...this.cubes.keys()].map((s: string) => strToVec(s));
    }

    /**
     *
     */
    public cubeCoordsStrs() : string[] {
        return [...this.cubes.keys()];
    }

    /**
     *
     */
    public allMoves(): Move[] {
        // Convert the values of the 'moves' Map into an array of Move arrays.
        const arrayOfMoveArrays: Move[][] = [...this.moves.values()].map(m => m.getMoves());
        // Flatten the array of Move arrays into a single array of Moves.
        return arrayOfMoveArrays.reduce((acc, moveArray) => acc.concat(moveArray), []);
    }

    /**
     * removes a move
     * @param key Move or THREE.Vector3 object to remove
     */
    public removeMove(key: Vector3 | Move){
        if (key instanceof Move){
            // remove move from move cell
            let movegroup = this.moves.get(vecToStr(key.pos));
            movegroup.clearFace(key.direction);
            // if move cell has no more moves now
            if (movegroup.isEmpty()){
                // delete it
                return this.moves.delete(vecToStr(key.pos))
            }
        }
        else {
            // remove move position
            this.moves.delete(vecToStr(key));
        }
    }

    /**
     *
     * @param key
     */
    public moveFromIdx(key: number) : Move {
        // mostly for random indexing, since these indexes won't stay still for long
        let moves = this.allMoves();
        console.assert(key < moves.length)
        return moves[key]
    }

    /**
     *
     */
    public numMoves() : number {
        return this.allMoves().length
    }

    /**
     *
     * @param position
     */
    hasMove(position: string | Vector3 | Move) : boolean {
        if (position instanceof Vector3 && position.isVector3){
            return this.moves.has(vecToStr(position))
        }
        else if (position instanceof  Move){
            return this.hasMove(position.pos);
        }
        else {
            console.assert(position instanceof String)
            return this.moves.has(<string>position)
        }
    }

    /**
     *
     * @param f
     */
    forEachCube(f: (PolycubeCubeInst, string) => void){
        this.cubes.forEach(f);
    }

    /**
     *
     * @param position
     */
    getMoveCell(position: string | Vector3 | number) : PolycubeCellMove {
        if (position instanceof Object && position.isVector3){
            return this.getMoveCell(vecToStr(position));
        }
        else if (typeof position == "string"){
            return this.moves.get(position);
        }
        else {
            console.assert(typeof position == "number")
            throw Error();
        }
    }

    /**
     *
     * @param position
     */
    hasCube(position: string | Vector3) : boolean {
        if (position instanceof Object && position.isVector3){
            return this.cubes.has(vecToStr(position));
        }
        else {
            console.assert(typeof position == "string")
            return this.cubes.has(<string>position);
        }
    }

    /**
     *
     * @param position
     */
    getCube(position: string | Vector3 | number) : PolycubeCubeInst<{patch: PolycubePatch, srcdiridx: number}> {
        let c;
        if (position instanceof Object && position.isVector3){
            c = this.cubes.get(vecToStr(position));
        }
        else if (typeof position == "string"){
            c = this.cubes.get(position);
        }
        else {
            console.assert(position instanceof Number)
            // @ts-ignore
            c = this.cubes.get(this.cubes.keys()[position]);
        }
        console.assert(c != void 0, `Invalid position type ${typeof position}`);
        return c;
    }

    numCubes() : number{
        return this.cubes.size;
    }

    numCubeTypes() :number{
        return this.cube_types.length;
    }

    public getCubes() : PolycubeCubeInst<{patch: PolycubePatch, srcdiridx, localdir: number}>[]{
        return [...this.cubes.values()];
    }

    /**
     * @return the system rule (list of cube types) in decimal serial form
     */
    getRuleStr() : string {
        return ruleToDec(this.cube_types);
    }

    /**
     @param c1
     @param c2
     @param color the color of the connection (should be able to infer but guess not)
     */
    addConnection(c1: PolycubeCubeInst<PatchWrapper>,
                  c2: PolycubeCubeInst<PatchWrapper>, color: number) {
        // create connection object
        let conn = new PolycubeSystemConnection(c1, c2, color);
        this.connections.push(conn);

        // locate patch 1 idx
        let patch1idx = RULE_ORDER.indexOf(c2.position.clone().sub(c1.position));
        // assign connection object
        c1.connections[patch1idx] = conn;
        // update state
        c1.set_patch_bound(patch1idx);

        // locate patch 2 idx
        let patch2idx = RULE_ORDER.indexOf(c1.position.clone().sub(c2.position));
        // assign connection object
        c2.connections[patch2idx] = conn;
        // update state
        c2.set_patch_bound(patch2idx);

        for (const cube of Array(c1, c2)) {
            if (cube.updateState()){
                this.updateMoves(cube);
            }
        }
    }

    public countConnections(c: PolycubeCubeInst<PatchWrapper>) : [number, number] {
        let connectionCount = 0;
        let mismatchCount = 0;
        // loop directions, relative to self
        RULE_ORDER.forEach((dir, i)=>{
            if (!c.hasPatch(i)){
                return;  // ignore zeroes
            }
            if (!c.patch_is_active(i)){
                return;
            }
            // compute adjacent position
            let adj = c.position.clone().add(dir);
            // if cube in adjacent position
            if (this.hasCube(adj)){
                let c2 = this.getCube(adj);
                let c2patch = c2.patch(dir.clone().negate());
                if (c2patch.color == 0){
                    return; // ignore zeroes
                }
                if (!c2.patch_is_active(RULE_ORDER.indexOf(dir.clone().negate()))){
                    return;
                }
                if (this.compatibleColors(c2patch.color, c.patch(i).color)){
                    if (!this.torsion || c2patch.alignDir.equals(c.patch(i).alignDir)){
                        connectionCount++;
                    }
                    else {
                        mismatchCount++;
                    }
                }
                else {
                    mismatchCount++;
                }
            }
        });
        return [connectionCount, mismatchCount];
    }

    /**
     * creates a new cube instance from the given cube type
     * override for subclasses!
      */
    instantiateCube(cube_type: PolycubeCubeType, name: string, position: Vector3 = new Vector3(0,0,0)) : DefaultCubeInst {
        return new PolycubeCubeInst(
            cube_type,
            position,
            Quaternion.prototype.identity(),
            name
        );
        // 	console.assert(position instanceof Vector3);
        // 	// copy object
        // 	let copy = this.clone(instance_name, false);
        // 	copy.rotation = new THREE.Quaternion().identity();
        //
        // 	// copy connections
        // 	// TODO: do connections actually do anything?
        // 	copy.connections = RULE_ORDER.map(x => false);
        //
        // 	// copy position
        // 	copy.position = position;
        // 	//https://stackoverflow.com/a/8084248
        //
        // 	// create new personal name
        // 	copy.personalName = instance_name;
        // 	return copy;
    }

    /**
     Need both rule and ruleIdx to determine color as the rule might be rotated
     * @param cube_instance
     */
    addParticle(cube_instance: PolycubeCubeInst<PatchWrapper>) {
        // Go through all non-zero parts of the rule and add potential moves
        // Josh Note: the system moves through the rule in a non-random order. is this an issue?

        let position = cube_instance.position;
        let [connections, mismatches] = this.countConnections(cube_instance);
        // console.assert(connections > 0);
        // todo: move the following assert elsewhere
        // console.assert(this.mismatches_allowed || mismatches == 0)

        this.cubes.set(vecToStr(position), cube_instance);

        this.currMinCoord = Math.min(this.currMinCoord, ...position.toArray());
        this.currMaxCoord = Math.max(this.currMaxCoord, ...position.toArray());

        // update moves for newly added particle
        this.updateMoves(cube_instance);

        this.cubeTypeCount.set(cube_instance.type_id(), this.cubeTypeCount.get(cube_instance.type_id()) + 1);
        this.mismatches += mismatches;
        this.matches += connections;
    }

    /**
     updates the moves map after a cube is added or has its state change
     */
    updateMoves(cube: PolycubeCubeInst<PatchWrapper>) {
        cube.iter_patches((patch)=>this.updateMove(cube, patch.patch));
        console.assert((this.numMoves() == 0) == (this.numMovePositions() == 0))
    }
    updateMove(cube: PolycubeCubeInst<PatchWrapper>, patch: PolycubePatch){
        // check if patch has colorx
        if (!patch.has_color()) {
            return; // nothing to do here
        }
        // find move position
        let movePos: Vector3 = cube.position.clone().add(RULE_ORDER[patch.direction])
        // check if move is outside area limit

        // if there's a cube where the move would go
        if (this.hasCube(movePos)){
            // There is already a cube at pos,
            // no need to add this neigbour to moves
            return;
        }

        // direction of face on move should be opposite of patch direction
        let move_direction: Vector3 = RULE_ORDER[patch.direction].clone().negate();
        console.assert(move_direction.length() == 1, "Distance between cube and connected cube should always be 1 unit.")
        // compute direction index of move direction
        let dirIdx = RULE_ORDER.indexOf(move_direction);
        console.assert(dirIdx >= 0)
        //if the face is not active
        if (!cube.patch_is_active(patch.direction))
        {
            // if move key is in moves, remove it
            if (this.hasMove(movePos) && this.getMoveCell(movePos).speciesFit[dirIdx]) {
                // face in moves, should be removed now that it's inactive
                this.getMoveCell(movePos).speciesFit[dirIdx] = null;
            }
            return;
        }

        // if move key isn't currently in moves
        if (!this.getMoveCell(movePos)) {
            //if there's not already a move in this position in the moves list,
            // add a new move to the moves list for this
            this.moves.set(vecToStr(movePos), new PolycubeCellMove(movePos));
        }
        //Make sure we haven't written anything here before:
        let move_cell = this.getMoveCell(movePos);
        // if there's a move in this position+direction w/ same color and align, it was just added earlier & is fine
        if (move_cell.speciesFit[dirIdx] != null){
            console.assert((move_cell.speciesFit[dirIdx].color == patch.color) &&
                move_cell.speciesFit[dirIdx].alignDir.distanceTo(patch.alignDir) < 1e-5,
                `There's somehow already a move in position ${vecToStr(movePos)} facing direction ${dirIdx}`)
        }
        this.getMoveCell(movePos).setFace(
            move_direction,
            patch.color,
            patch.alignDir);
    }

    getReasonableOverlapTransforms(){
        let transforms = [];
        for (let x = this.currMinCoord; x <= this.currMaxCoord; x++) {
            for (let y = this.currMinCoord; y <= this.currMaxCoord; y++){
                for (let z = this.currMinCoord; z <= this.currMaxCoord; z++){
                    transforms.push(new THREE.Vector3(x,y,z));
                }
            }
        }
        return transforms;
    }

    /**
     * computes if two polycubes are equivalent (made up of same cubes in same positions and rotations)
     * @param other
     */
    isEquivelent(other: Polycube) : boolean{
        if (this.cubes.size != other.cubes.size){
            return false;
        }
        if (this.cubeTypeCount.size == other.cubeTypeCount.size && this.cubeTypeCount.toString() != other.cubeTypeCount.toString()){
            return false;
        }
        // todo: more filtering
        let all_rotations = allRotations();
        let all_translations = this.getReasonableOverlapTransforms();
        let order = all_rotations.length * all_translations.length; // gimme some idea of how fucked this is
        // iterate through all possible rotation-combination pairs
        return all_rotations.some(rot => {
            return all_translations.some(tran => {
                // if any rot, tran pair applied to other produces a polycube which is equoivalent to
                // this, the polycubes are equivalent
                let t = new THREE.Matrix4();
                t.setFromMatrix3(rot);
                t.setPosition(tran);
                // check match cube-by-cube
                this.forEachCube((cube: PolycubeCubeInst<PatchWrapper>, _:string) => {
                    // apply transform
                    let position = cube.position.clone().applyMatrix4(t);
                    // if no cube is at the position transformed, return false
                    if (!other.hasCube(position))
                        return false;
                    let other_cube = other.getCube(position);
                    // check that type IDs match
                    if (other_cube.type_id != cube.type_id) // todo: skip / alternative to typechecking?
                        return false
                    // aaahhhhhh
                    // check that rotations match
                    if (cube.rotation.clone().multiply(new Quaternion().setFromRotationMatrix(t)).equals(other_cube.rotation))
                        return false
                    // todo: check that states match? two polycubes with non-matching q
                });
                return true;
            })
        })
    }
}