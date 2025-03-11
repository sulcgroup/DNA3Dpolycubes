import * as THREE from "three";
import {BoxGeometry, BufferGeometry, Group, Material, Mesh, MeshLambertMaterial, Quaternion, Vector3} from "three";
import * as $ from "jquery";
import {cubeTypeColor, DynamicPolycubeCubeType, FACE_NAMES, patchColor, PolycubeCubeType} from "./rule";
import {FACE_ROTATIONS, RULE_ORDER, selectColor, vecToStr} from "./utils";
import faceShape from "../ui/face.svg";
import {PolycubePatch} from "./patch";

export type PatchWrapper = {patch: PolycubePatch, srcdiridx: number};

export class PolycubeCubeInst<T extends PatchWrapper> {
    protected cube_type: PolycubeCubeType
    position: Vector3
    rotation: Quaternion
    protected personalName: string

    connections: (boolean | PolycubeSystemConnection)[];

    protected state: boolean[];

    /**
     * very simple constructor which assigns stuff
     * @param ct
     * @param pos
     * @param rot
     * @param name
     */
    constructor(ct: PolycubeCubeType, pos: Vector3, rot: Quaternion, name: string) {
        console.assert(ct != undefined);
        this.cube_type = ct;
        this.position = pos;
        this.rotation = rot;
        this.personalName = name;
        this.state = Array(this.cube_type.state_size()).map(_=>false);
        this.state[0] = true;
        this.connections = Array(RULE_ORDER.length).map(_=>false);
    }

    public getPersonalName() : string{
        return this.personalName;
    }

    public getState() : boolean[]{
        return this.state;
    }

    public getType() : PolycubeCubeType {
        return this.cube_type;
    }

    /**
     * sets a patch to bound
     * @param idx LOCAL INDEX OF PATCH
     */
    public set_patch_bound(idx: number) {
        this.state[this.patch_by_local(idx).state_var] = true;
    }

    public type_id() :number {
        return this.cube_type.type_id;
    }

    /**
     * patch indexing in this class is.... nightmarinsh
     * returns the patch at the given direction, indexed using the orientation
     * of the type cube
     * @param idx
     */
    public patch_from_cube_type(idx: number | Vector3): PolycubePatch {
        if (idx instanceof Vector3){
            return this.patch_from_cube_type(RULE_ORDER.indexOf(idx))
        }
        else {
            return this.cube_type.patches[idx];
        }
    }

    /**
     * given a direction (either a vector or an index) on the rotated cube,
     * returns the patch information for that direction
     * @param direction
     */
    public patch(direction: number | Vector3) : PolycubePatch {
        if (direction instanceof Vector3){
            let direction_on_type: Vector3 = direction.clone().applyQuaternion(this.rotation.clone().invert() );
            direction_on_type.round();
            console.assert(direction_on_type.length() == 1);
            let uncorrected_ct_patch = this.patch_from_cube_type(direction_on_type);
            return this.correct_patch(uncorrected_ct_patch).patch;
        }
        else {
            console.assert(typeof direction == "number");
            return this.patch(RULE_ORDER[direction]);
        }
    }
    /**
     * returns the
     * @param idx
     */
    public patch_by_local(idx: number | Vector3) : PolycubePatch {
        if (idx instanceof Vector3){
            return this.patch_by_local(RULE_ORDER.indexOf(idx));
        }
        else {
            return this.correct_patch(this.patch_from_cube_type(idx)).patch;
        }
    }
    /**
     * returns mapping of patch GLOBAL (aka CUBE TYPE)
     * direction indexes to LOCAL (aka ROTATED CUBE DIRIDX) patch objects
     * so we do NOT expect key == val.direction in all or even most circumstances!!!
     */
    public get_patch_map() : Map<number, T> {
        let mapping = new Map<number, T>();
         this.cube_type.patches.forEach((p, i) => {
                if (p.color) {
                    mapping.set(i, this.correct_patch(p));
                    // assertion?
                }
            });
         return mapping;
    }

    /**
     * adds visualization info if applicable and corrects rotation for patches from cube type
     * @param p patch object, pulled from cube type
     * @protected
     */
    protected correct_patch(p: PolycubePatch): T {
        return {patch: p.rotate(this.rotation), srcdiridx: p.direction} as T;
    }

    public iter_patches(f: (T) => void) {
        [...this.get_patch_map().values()].forEach(f);
    }

    public clone() : PolycubeCubeInst<T> {
        return new PolycubeCubeInst<T>(this.getType(), this.position.clone(), this.rotation.clone(), this.getPersonalName());
    }

    rotate(q: Quaternion) {
        console.assert(q !== void 0)
        console.assert(this.rotation !== void 0)
        // cache rotation, important for later
        // TODO: use this exclusively
        this.rotation = q.multiply(this.rotation)
    }

    /**
     * @param idx global patch idx
     */
    public patch_is_active(idx: number) : boolean {
        if (!this.cube_type.patch_allosterically_controlled(idx)){
            return true;
        }
        else if (this.get_patch_map().get(idx).patch.activation_var > 0){
            return this.state[this.get_patch_map().get(idx).patch.activation_var];
        }
        else {
            return !this.state[this.get_patch_map().get(idx).patch.activation_var];
        }
    }

    /**
     * given a direction on the rotated cube, returns true if there is a nonzero patch on that location
     * @param idx
     */
    public hasPatch(idx: number) : boolean{
        let direction_on_type: Vector3 = RULE_ORDER[idx].clone().applyQuaternion(this.rotation.clone().invert() );
        direction_on_type.round();
        return this.cube_type.patches[RULE_ORDER.indexOf(direction_on_type)].color != 0;
    }

    //
    /**
     * rotates instance around axis, returns self
     * https://stackoverflow.com/a/25199671
     * @param vFrom
     * @param vTo
     */
    rotateInstanceFromTo(vFrom: Vector3, vTo: Vector3) {
        let quaternion = new THREE.Quaternion(); // create one and reuse it
        quaternion.setFromUnitVectors(vFrom, vTo);
        this.rotate(quaternion);
        return this;
    }

    rotateInstanceAroundAxis(axis: Vector3, angle: number) {
        let quaternion = new THREE.Quaternion(); // create one and reuse it
        quaternion.setFromAxisAngle(axis, angle);
        this.rotate(quaternion);
        return this;
    }

    public updateState() :boolean {
        let changes:boolean[];
        [changes, this.state] = this.cube_type.reeval_conditionals(this.state)
        return changes.reduce((a,b) => {return a || b;});
    }
}

type PatchVisual = {
    // the large square of the patch
    big: THREE.Mesh;
    // the small square part of the patch that indicates orientation
    small: THREE.Mesh;
    material: MeshLambertMaterial;
    env_vis_rep: JQuery;
}


export class VisualPolycubeCubeInst extends PolycubeCubeInst<{patch: PolycubePatch, srcdiridx: number, vis: PatchVisual}> {
    envvis: JQuery;
    cube_vis: Mesh;

    // array of objects indexed identically to patches in this.cube_type
    patch_instances: Map<number,PatchVisual>

    constructor(ct: PolycubeCubeType,
                pos: Vector3,
                rot: Quaternion,
                name: string) {
        super(ct, pos, rot, name);
        this.patch_instances = new Map<number, PatchVisual>();
        this.getType().patches.filter(p => p.color).forEach(p => {
            this.patch_instances.set(p.direction, {
                big: undefined,
                small: undefined,
                material :undefined,
                env_vis_rep: undefined
            })
        });
    }

    protected correct_patch(p: PolycubePatch): {patch: PolycubePatch, srcdiridx: number, vis: PatchVisual} {
        console.assert(!this.is_drawn() || p.color == 0 || (this.patch_instances.has(p.direction) && this.patch_instances.get(p.direction).big !== undefined))
        return {
            patch: super.correct_patch(p).patch,
            srcdiridx: p.direction,
            vis: this.patch_instances.get(p.direction)
        };
    }

    rotate(q: Quaternion) {
        console.assert(!this.is_drawn(), "Please do not rotate a cube which has been drawn!!!");
        super.rotate(q);
    }

    is_drawn() : boolean {
        return this.cube_vis !== undefined;
    }

    draw(connector_geometry: BoxGeometry, cube_mat: MeshLambertMaterial, cube_geo: BufferGeometry) : Group {

        let cube_vis_group = new THREE.Group();

        // add binding sites on cube faces
        cube_vis_group.position.copy(this.position);
        cube_vis_group.name = "Cube";

        // add env browser elements
        this.envvis = $("<table class='env-vis-box'>");
        this.envvis.attr('name', this.personalName);
        let header_label_colspan;
        let header_val_colspan;
        if (this.cube_type.hasAllostery()) {
            header_label_colspan = 2;
            header_val_colspan = 3;
        }
        else {
            header_label_colspan = 1;
            header_val_colspan = 1;
        }
        // cube type name
        let row = $("<tr>")
        this.envvis.append(row);
        row.append($(`<td colspan="${header_label_colspan}"><b>Cube Type</b>:</td>`));
        row.append($(`<td colspan="${header_val_colspan}">${this.cube_type.typeName}</td>`));
        // personal name
        row = $("<tr>")
        this.envvis.append(row);
        row.append($(`<td colspan="${header_label_colspan}"><b>Cube Name</b>:</td>`));
        row.append($(`<td colspan="${header_val_colspan}">${this.personalName}</td>`))
        // cube position
        row = $("<tr>");
        this.envvis.append(row);
        row.append($(`<td colspan="${header_label_colspan}"><b>Position</b>:</td>`));
        row.append($(`<td colspan="${header_val_colspan}"> ${vecToStr(this.position)}</td>`))
        // patch info column headers
        row = $("<tr>");
        this.envvis.append(row);
        row.append("<td style=\"border-bottom: 2px solid black;\">Dir.</td>");
        row.append("<td style=\"border-bottom: 2px solid black;\">Ori.</td>");
        row.append("<td style=\"border-bottom: 2px solid black;\">Color</td>");
        if (this.getType().hasAllostery()){
            row.append("<td style=\"border-bottom: 2px solid black;\">State Var.</td>");
            row.append("<td style=\"border-bottom: 2px solid black;\">Act. Var.</td>");
        }

        // init patch graphics
        this.iter_patches((p: {patch: PolycubePatch, srcdiridx: number, vis: PatchVisual}) => {
            // skip non-patches
            if (p.patch.has_color()) {
                let patch_color = patchColor(p.patch.color);
                // init patch visual

                // ---------------- patch info in env browser ----------------- //
                // patch dir
                let face_label = $(`<tr class='face-label' name='${FACE_NAMES[p.patch.direction]}'>`);
                face_label.append($(`<td style="font-weight: bold;">${FACE_NAMES[p.patch.direction]}:</td>`))
                // patch align
                let face_rot = FACE_ROTATIONS[p.patch.direction].angleTo(p.patch.alignDir) * (2 / Math.PI);
                face_label.append($(`<td style="width: 35px;"><img src=${faceShape} class="face-img rot${face_rot}" style="cursor: auto;"/></td>`));
                face_label.append($(`<td>${p.patch.color}</td>`));
                if (this.getType().hasAllostery()) {
                    // state variable
                    face_label.append($(`<td>${p.patch.state_var}</td>`));
                    // activation variable
                    face_label.append($(`<td>${p.patch.activation_var}</td>`));
                }

                this.envvis.append(face_label);

                // construct color material
                let mat = new THREE.MeshLambertMaterial({
                    color: patch_color
                });

                // ------------------- draw scene stuff ----------------- //
                p.vis = {
                    material: mat,
                    big: new THREE.Mesh(connector_geometry, mat),
                    small: new THREE.Mesh(connector_geometry, mat),
                    env_vis_rep: face_label

                }
                p.vis.material.transparent = true;
                p.vis.material.opacity = this.state[p.patch.activation_var] ? 1.0 : 0.5;
                if (p.patch.color >= 0) {
                    p.vis.material.emissive = p.vis.material.color.clone().addScalar(-0.5);
                } else {
                    p.vis.material.color.addScalar(-0.1);
                }
                cube_vis_group.add(p.vis.small);
                cube_vis_group.add(p.vis.big);
                p.vis.big.position.add(
                    RULE_ORDER[p.patch.direction].clone().multiplyScalar(0.4)
                );
                // dimension of patch
                let dim = RULE_ORDER[p.patch.direction].clone();
                dim.setX(Math.abs(dim.x)).setY(Math.abs(dim.y)).setZ(Math.abs(dim.z));
                // flatten visual mesh on appropriate axis
                p.vis.big.scale.copy(
                    new THREE.Vector3(1, 1, 1).sub(dim.multiplyScalar(0.65))
                );
                p.vis.small.scale.copy(p.vis.big.scale);
                p.vis.small.position.copy(p.vis.big.position);
                p.vis.small.position.add(p.patch.alignDir.clone().multiplyScalar(0.2));
                p.vis.small.scale.multiplyScalar(0.5);
            }
            console.assert(this.patch_instances.has(p.srcdiridx));
            this.patch_instances.set(p.srcdiridx, p.vis);

        });

        this.envvis.append($(`<tr>
			<td style="border-bottom: 2px solid black;">Index</td>
			<td style="border-bottom: 2px solid black;">Value</td>
			<td colspan="2" style="border-bottom: 2px solid black;">Conditional</td>
		</tr>`));
        this.envvis.css('border-color', cubeTypeColor(this.type_id()))

        if (this.cube_type instanceof DynamicPolycubeCubeType){
            this.state.forEach((statevar, i) =>{
                let staterow = $(`<tr class='state-label'>`);
                staterow.append($(`<td>${i}</td>`));
                staterow.append($(`<td name="state${i}">${statevar}</td>`));

                if (i != 0) {
                    staterow.append($(`<td colspan="2">Dynamic Effect</td>`));
                }
                else {
                    staterow.append($(`<td>`));
                }
                this.envvis.append(staterow);
            });
        }

        // draw main cube last for complecated reasons
        this.cube_vis = new THREE.Mesh(cube_geo, cube_mat);
        // cube.bindings_vis = [];
        cube_vis_group.add(this.cube_vis);

        return cube_vis_group;
    }

    /**
     *
     */
    refresh_patch_active_visual(p: {patch: PolycubePatch, vis: PatchVisual}) {
        if (!this.patch_is_active(p.patch.direction)) {
            let bg = `repeating-linear-gradient(45deg, ${patchColor(p.patch.color)}, ${selectColor(p.patch.color)} 1%, white 1%, white 2%)`;
            p.vis.env_vis_rep.css('background-image', bg);
            p.vis.env_vis_rep.css('background-color', '');
            (p.vis.big.material as Material).opacity = 0.5;
        } else {
            p.vis.env_vis_rep.css('background-image', '');
            p.vis.env_vis_rep.css('background-color', patchColor(p.patch.color));
            (p.vis.big.material as Material).opacity = 1;
        }
    }

    public clone(): VisualPolycubeCubeInst {
        console.assert(!this.is_drawn());
        // there is no better way to do this than copy-pasting the line from PolycubeCubeInst::clone
        // because Javascript was designed by an illiterate homophobic macaque
        return new VisualPolycubeCubeInst(
            this.getType(),
            this.position.clone(),
            this.rotation.clone(),
            this.getPersonalName());
    }

    set_state(state: boolean[]) {
        this.state = state;
    }
}

export class PolycubeSystemConnection {
    cube_1: PolycubeCubeInst<any>;
    cube_2: PolycubeCubeInst<any>;
    color: number;

    /**
     @param cube_1 the first cube for this connection
     @param cube_2 the second cube for this connection
     @param color the color of the connection
     */
    constructor(cube_1, cube_2, color) {
        this.cube_1 = cube_1;
        this.cube_2 = cube_2
        this.color = Math.abs(color)
    }
}

export type DefaultCubeInst = PolycubeCubeInst<any>;
