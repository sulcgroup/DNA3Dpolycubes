import {FACE_ROTATIONS, getSignedAngle, resize_array, RULE_ORDER, selectColor, strToVec} from "./utils"
import * as THREE from "three"
import {BoxGeometry, BufferGeometry, Group, MeshLambertMaterial, Quaternion, Vector3} from "three"
import {PolycubePatch} from "./patch";

export const FACE_NAMES = ["left", "right", "bottom", "top", "back", "front"];

export function cubeTypeColor(i: number){
	return selectColor(i);
}

export function patchColor(i: number) : string{
	return selectColor( Math.abs(i) - 1);
}

export class PolycubeCubeType {
	readonly typeName: string;
	public readonly type_id: number;
	patches: PolycubePatch[];
	// cube instance vars

	protected _state_size: number;

	constructor(
			typeid: number,
			patches: PolycubePatch[],
			type_name: string | undefined,
			stateSize = 0) {
		console.assert(typeof typeid == "number")
		console.assert(patches.length == 6)
		this.type_id = typeid;
		if (typeof type_name == "undefined"){
			type_name = `CT${typeid}`;
		}
		this.typeName = type_name;
		this.patches = patches;
		// add 1 to state size to incorporate tautology variable
		this._state_size = stateSize
	}

	/**
	 * copies this cube type
	 * @param newName
	 * @param newtype
	 */
	clone(newtype: number){
		return new PolycubeCubeType(
            newtype,
            this.patches.map(patch => {
                return new PolycubePatch(
                    patch.color,
                    patch.direction,
                    patch.alignDir,
                    patch.state_var,
                    patch.activation_var
                );

            }),
            `CT${newtype}`,
            this.state_size());
	}

	/**
	 *
	 * @returns {number} number of state variables (not counting tautology variable)
	 */
	state_size() : number{
		return this._state_size;
	}

	// // https://stackoverflow.com/a/25199671
	// /**
	//  * rotates cube type. utility of this. unclear.
	//  * @param q
	//  */
	rotate(q: Quaternion) : PolycubeCubeType{
		console.assert(q !== void 0)
		// console.assert(this.rotation !== void 0)
		let connections = [];
		let new_patches = [];
		// loop patches
		this.patches.forEach((patch, i) => {
			console.assert(i == patch.direction);
			let oldDir = patch.direction;
			let newDir = RULE_ORDER[patch.direction].clone().applyQuaternion(q).round();
			let iNewFace = RULE_ORDER.indexOf(newDir);
			// if ('connections' in this && this.connections !== void 0) {
			// 	connections[iNewFace] = this.connections[i];
			// }
			new_patches[iNewFace] = patch.rotate(q);
		});
		let rotated_cube = this.clone(this.type_id);
		rotated_cube.patches = new_patches;
		return rotated_cube;
	}

	hasAllostery() : boolean {
		return false;
	}

	patch_color(i: number) : number{
		return this.patches[i].color
	}

	face(i: number) : PolycubePatch {
		return this.patches[i];
	}

	toDec() : string {
		return this.patches.map((patch,i)=>{
			if (patch.color === 0) {
				return '';
			} else {
				let theta = getSignedAngle(FACE_ROTATIONS[i],
					patch.alignDir,
					RULE_ORDER[i]);
				let orientation = Math.round(theta * (2 / Math.PI) + 4) % 4;
				if (orientation == 0 && patch.state_var == 0 && patch.activation_var == 0){
					return `${patch.color}`
				}
				else if (patch.state_var == 0 && patch.activation_var == 0){
					return `${patch.color}:${orientation}`
				}
				else {
					return `${patch.color}:${orientation}:${patch.state_var}:${patch.activation_var}`;
				}
			}
		}).join('|');
	}

	num_colored_faces() : number{
		return this.patches.map(patch=>patch.color).reduce((count, color) => count + Math.min(Math.abs(color), 1));
	}

	colored_faces_idxs() : number[] {
		return this.patches.filter(p=>p.color).map(p=>p.direction);
	}

	patch_colors() : Set<number> {
		return new Set(this.patches.filter(patch=> patch.color != 0).map(patch=>patch.color));
	}

	public getPatchCoords() {
		return this.patches.filter(patch=>{
			return patch.color !== 0;
		}).map(patch=>RULE_ORDER[patch.direction]);
	}

	public patch_allosterically_controlled(n: number) : boolean{
		throw Error();
	}

	public reeval_conditionals(state: boolean[]) : [boolean[], boolean[]] {
		return [[], []]
	}

	/**
	 *
	 * @param other
	 */
	public equals(other: PolycubeCubeType) : boolean {
		return (this.type_id == other.type_id)
			&& (this.typeName == other.typeName)
			&& (this.state_size() == other.state_size())
			&& this.patches.map((p, i) => {
				let other_patch = other.patches[i];
				return (p.color == other_patch.color)
					&& (RULE_ORDER[p.direction].distanceTo(RULE_ORDER[other_patch.direction]) < 1e-5) // prob a freebie
					&& (p.alignDir.distanceTo(other_patch.alignDir) < 1e-5)
					&& (p.state_var == other_patch.state_var)
					&& (p.activation_var == other_patch.activation_var);
			}).reduce((a,b)=> a && b);
	}
}
// export class StaticPolycubeCubeType extends PolycubeCubeType {
// 	conditionals: string[];
// 	/**
// 	 * The Static Polycube Cube is a subclass of PolycubeCube that has state size 12
// 	 * State 0 is the tautology varaible (always 1)
// 	 * States 1 - 6 are binding on patches (corresponding with RULE_ORDER in an unrotated cube)
// 	 * States 7 - 12 are activation of patches (corresponding with RULE_ORDER in an unrotated cube)
// 	 * @param type_id
// 	 * @param name
// 	 * @param patches
// 	 * @param conditionals
// 	 */
// 	constructor(type_id, patches?, type_name?: string, conditionals=Array(RULE_ORDER.length).fill("")) {
// 		super(type_id, patches, type_name, 2 * RULE_ORDER.length); // constructor will add tautology variable
//
// 		// set non-allosteric activation variables to true
// 		for (let i = 0; i < RULE_ORDER.length; i++){
// 			if (!conditionals[i].trim()){
// 				this.state[i + RULE_ORDER.length + 1] = true;
// 			}
// 		}
// 		// conditionals which determine values for state vars 7 - 12
// 		// does NOT need to be reindexed with rotations b/c patches keep their
// 		// state and activation var idxs when they're rotated
// 		this.conditionals = conditionals;
// 		this.reeval_conditionals();
// 	}
//
// 	hasAllostery() {
// 		return this.patches.some((p,i) => {
// 			return this.patch_allosterically_controlled(i);
// 		} );
// 	}
//
// 	patch_allosterically_controlled(i){
// 		if (this.patches[i].activation_var  == 0){
// 			return false;
// 		}
// 		let conditional = this.conditionals[this.patches[i].activation_var - RULE_ORDER.length - 1];
// 		return (conditional.trim().length > 0) && (conditional != "(true)");
// 	}
//
// 	patch_is_active(idx){
// 		return !this.patch_allosterically_controlled(idx) || this.state[this.patches[idx].activation_var];
// 	}
//
// 	set_patch_bound(idx){
// 		this.state[this.patches[idx].state_var] = true;
// 	}
//
// 	/**
// 	 re-evaluates the activation states of the faces, based on which faces are bound. returns
// 	 true if any activation state has changed, false otherwise
// 	 */
// 	reeval_conditionals() : boolean[] {
// 		// init state change list
// 		let state_changed = Array(RULE_ORDER.length).fill(false);
//
// 		// copy state variables 1-6 to temporary variable b
// 		// b should start at 0 for backwards compatibility reasons
// 		let b = this.state.slice(1, 1 + RULE_ORDER.length);
//
// 		for (let i_conditional = 0; i_conditional < RULE_ORDER.length; i_conditional++){
// 			if (this.conditionals[i_conditional].trim()) {
// 				let cond_var_idx = 1 + RULE_ORDER.length + i_conditional;
// 				// save current value for comparison
// 				let prev_var_val = this.state[cond_var_idx];
// 				this.state[cond_var_idx] = eval(this.conditionals[i_conditional]);
// 				// check if state has changed
// 				state_changed[cond_var_idx] = prev_var_val != this.state[cond_var_idx];
// 			}
// 		}
//
// 		// // iterate conditionals for activations
//
// 		// if any states changed and this cube has a visualization in three.js, update visualizer
// 		if (state_changed.some(Boolean) && 'personalName' in this) {
// 			// loop patches
// 			this.patches.forEach((patch, i_patch) => {
// 				// if the patch has a color and is allosterically controlled
// 				if (patch.color != 0 && this.patch_allosterically_controlled(i_patch) && state_changed[patch.activation_var]) {
// 					// grab var for binding visualization
// 					let activation = this.state[patch.activation_var];
//
// 					// loop binding visualization components
// 					(patch.visual.big.material as Material).opacity = (activation ? 1 : 0.5)
// 				}
// 			});
// 		}
// 		return state_changed;
// 	}
//
// 	clone(newName? : string | undefined, newtype: number | false = false) {
// 		let copy = super.clone(newName, newtype) as StaticPolycubeCubeType;
// 		copy.conditionals = [...this.conditionals];
// 		return copy;
// 	}
//
// 	instantiate(position: Vector3, instance_name: string) {
// 		let instance = toStatic(super.instantiate(position, instance_name));
// 		instance.rotation = new THREE.Quaternion().identity();
//
// 		instance.conditionals = [...this.conditionals];
// 		instance.reeval_conditionals();
//
// 		return instance;
// 	}
//
// 	draw(connector_geometry: BoxGeometry, cube_mat: MeshLambertMaterial, cube_geo: BufferGeometry) : Group  {
// 		let g = super.draw(connector_geometry, cube_mat, cube_geo);
// 		this.state.forEach((statevar, i) =>{
// 			let staterow = $(`<tr class='state-label'>`);
// 			staterow.append($(`<td>${i}</td>`));
// 			staterow.append($(`<td name="state${i}">${statevar}</td>`));
// 			if (i - (RULE_ORDER.length + 1) >= 0)
// 			{
// 				let conditional = this.conditionals[i - (RULE_ORDER.length + 1)];
// 				staterow.append($(`<td colspan="2">${conditional}</td>`));
// 			}
// 			else {
// 				staterow.append($("<td>"));
// 			}
// 			this.envvis.append(staterow);
// 		});
// 		return g;
// 	}
//
// 	toDec() {
// 		return this.patches.map((patch,i)=>{
// 			let conditional = this.conditionals[patch.activation_var - (RULE_ORDER.length + 1)];
// 			if (patch.color === 0) {
// 				return '';
// 			} else {
// 				let orientation = Math.round(getSignedAngle(FACE_ROTATIONS[i], patch.alignDir, RULE_ORDER[i])*(2/Math.PI)+4)%4;
// 				return `${patch.color}:${orientation}:${conditional != "(true)" ? conditional : ""}`;
// 			}
// 		}).join('#')
// 	}
//
// 	rotate(q: Quaternion): PolycubeCubeType {
// 		return toStatic(super.rotate(q));
// 	}
// }

// todo: change this class name to the term we used in the paper once i remember what it was / can be bothered
type Effect = {
	sources: number[];
	target: number;
	// todo: activation energy
}

export class DynamicPolycubeCubeType extends PolycubeCubeType {
	effects: Effect[];
	constructor(type_id, patches?, name?, state_size=0, effects=[]) {
		super(type_id, patches, name, state_size);
		this.effects = effects;
	}

	hasAllostery() :boolean{
		return this.patches.some((p,i) => {
			return this.patch_allosterically_controlled(i);
		} );
	}

	patch_allosterically_controlled(i: number):boolean{
		return this.patches[i].activation_var != 0;
	}

	/**
	 * reevaluates state passed from a cube type instance
	 * @param state
	 * @return a tuple where the first element is CHANGES (each element is true if the reevaluation has changed the state var) and the second is the NEW STATE
	 */
	public reeval_conditionals(state: boolean[]) : [boolean[], boolean[]] {
		let change;
		let newState = [...state]
		do {
			change = false;
			this.effects.forEach(effect=>{
				let can_fire = effect.sources.every(s=>newState[s]);
				if (can_fire && !newState[effect.target]){
					change = true;
					newState[effect.target] = true;
				}
			});
		} while (change);
		let changes = state.map((x,i)=>{
			return x !== newState[i]
		});
		state = newState
		return [changes, state];
	}

	clone(newtype: number) {
		let copy = toDynamic(super.clone(newtype));
		copy.effects = this.effects.map((e: Effect) => {
			return JSON.parse(JSON.stringify(e));
		});
		return copy;
	}

	// instantiate(position: Vector3, instance_name: string) {
	// 	let instance : DynamicPolycubeCubeType = toDynamic(super.instantiate(position, instance_name));
	// 	instance.rotation = new THREE.Quaternion().identity();
	//
	// 	instance.effects = [...this.effects];
	// 	instance.reeval_conditionals();
	// 	return instance;
	// }


	update_state_size() {
		let new_state_size = this.state_size();
		this.patches.forEach(p=>{
			new_state_size = Math.max(new_state_size, p.state_var);
			new_state_size = Math.max(new_state_size, Math.abs(p.activation_var));
		});
		this._state_size = Math.max(new_state_size, ...this.effects.map(e=>{return e.target;}));
		// resize_array(this.state, new_state_size);
	}

	toDec() : string {
		let decstr = super.toDec();
		if (this.effects.length > 0){
			decstr += "@" + this.effects.map(e=> {
				return `[${e.sources.join(",")}]>${e.target}`;
			}).join(";")
		}
		return decstr;
	}

	// rotate(q: Quaternion): PolycubeCubeType {
	// 	return toDynamic(super.rotate(q));
	// }
}

function toDynamic(c: PolycubeCubeType) : DynamicPolycubeCubeType {
	let newCube = new DynamicPolycubeCubeType(
		c.type_id,
		c.patches,
		c.typeName,
		c.state_size()
	);
	Object.assign(newCube, c);
	return newCube;
}

// function toStatic(c: PolycubeCubeType) : StaticPolycubeCubeType {
// 	let newCube = new StaticPolycubeCubeType(
// 		c.type_id,
// 		c.patches,
// 		c.personalName
// 	);
// 	Object.assign(newCube, c);
// 	return newCube;
// }

/**

 */
export function ruleToDec(cubeTypes: PolycubeCubeType[]) : string{
	if (cubeTypes.length > 0) {
		return cubeTypes.map(ct => ct.toDec()).join('_');
	}
	else {
		return "";
	}
}

/**
 */
export function parseDecRule(rulesStr: string) : PolycubeCubeType[]{
	// split by underscore
    return rulesStr.split('_').map((s, j)=>{
		let patch_conditionals = RULE_ORDER.map(x=>"");
		var state_size = 1;
		let delimeter;
		if (s.match(/^([^|]*\|){5}[^|]*@?.*?$/)){
			delimeter = "|";
		}
		else {
			delimeter = "#";
		}
		var face_props: [string, string, string, string] | [string, string, string];
		let patches = s.split("@")[0].split(delimeter).map((face, i) => {
			let state_var: string = "0";
			let activation_var: string = "0";
			let i_color: number = 0;
			let r = FACE_ROTATIONS[i].clone();
			let i_orientation: number  = 0;
			let logic = "";
			if (!(face.trim() == "")) {
				// @ts-ignore
				face_props = face.split(":");
				let color = face_props[0];
				if (face_props.length > 1) {
					let patchOriIdx: string = face_props[1];
					// color, orientation, state var, activation var
					if (face_props.length == 4){
						[color, patchOriIdx, state_var, activation_var] = face_props;
						var i_state_var = parseInt(state_var);
						var i_activation_var = parseInt(activation_var);
						state_size = Math.max(state_size, Math.abs(i_state_var), Math.abs(i_activation_var));
					}
					if (face_props.length == 3) {
						[color, patchOriIdx, logic] = face_props;
						// var b for querying in inline js statements in static formulation
						// (DEPRECATED!!!)
						let b = Array(RULE_ORDER.length).fill(false);
						if (logic.trim().length > 0) {
							try {
								eval(logic);
							} catch (e) {
								throw `Malformed rule - malformed logical statement ${logic[i]}`;
							}
						}
					}
					i_orientation = parseInt(patchOriIdx);
					//[colors[i], orientation] = face.slice(0, logic_start).split(':').map(v=>parseInt(v));
					if (isNaN(i_orientation)) {
						throw "Malformed rule";
					}
					// if this string is in the dynamic format
					if (face_props.length == 3) {
						patch_conditionals[i] = logic;
					}
				}
				i_color = parseInt(color);
				if (isNaN(i_color)){
					throw "Malformed rule";
				}
				r.applyAxisAngle(RULE_ORDER[i], i_orientation * Math.PI / 2).round();
			}

			return new PolycubePatch(i_color,
				i, r,
				i_state_var,
				i_activation_var);
		});
		console.assert(patches.length == RULE_ORDER.length)
		// if (!patch_conditionals.some(x=>x.length > 0)){
		let effects;
		if (s.search("@") > -1) {
			effects = s.split("@")[1].split(";").map(e_str => {
				let [sources, target] = e_str.split(">");
				const regex = /-?\d+/g;
				const matches = sources.match(regex);

				// Convert string matches to numbers and return the resulting array
				return {
					target: Number(target),
					sources: matches.map(Number)
				};
			});
		}
		else {
			effects = [];
			state_size = 1;
		}
		 let ct = new DynamicPolycubeCubeType(j, patches, `CT${j}`, state_size, effects)
		ct.update_state_size()
		return ct;
		// }
		// else {
			// return new StaticPolycubeCubeType(j, patches, `CT${j}`, patch_conditionals);
		// }
    });
}

/**
https://stackoverflow.com/a/45054052 
*/
export function parseHexRule(ruleStr) {
    let ruleSize = RULE_ORDER.length;
    let rule = [];
    for (let i=0; i < ruleStr.length; i+=2*ruleSize) {
        let cube_type = [];
        //console.log("Rule ",(i/(2*ruleSize))+1);
        for (let j = 0; j<ruleSize; j++) {
            let face = ruleStr.substring(i+(2*j), i+(2*j) + 2);
            let binStr = (parseInt(face, 16).toString(2)).padStart(8, '0');
            let sign = parseInt(binStr[0], 2);
            let color = parseInt(binStr.substring(1,6),2);
            let orientation = parseInt(binStr.substring(6,8),2);

            let r = FACE_ROTATIONS[j].clone();
            r.applyAxisAngle(RULE_ORDER[j], orientation*Math.PI/2);
            r.round();
            cube_type.push( {'color': color * (sign ? -1:1), 'alignDir': r} );
        }
        rule.push(cube_type);
    }
    //added this line of code to switch to new method of storing rules
  	rule = rule.map((x,i)=>{
		  return new DynamicPolycubeCubeType(i, x.map((patch, i) => {
			  return new PolycubePatch(patch.color, i, patch.alignDir, 0, 0);
		  }), `CT${i}`);
	});
    return rule;
}
// Accepts either ~ or | as face separators
// and either . or : as rotation separators
// export function parseDecRuleOld(ruleStr) {
// 	return ruleStr.split('_').map((s,j)=>{
// 		return new StaticPolycubeCube(
// 			`CT${j}`,
// 			s.split(/[|.]/).map((face,i)=>{
// 				let color = 0;
// 				let orientation = 0;
// 				if (face !== '') {
// 					let faceVal = face.split(/[r:]/).map(v=>parseInt(v));
// 					color = faceVal[0]
// 					if (faceVal.length > 1){
// 						orientation = faceVal[1];
// 					}
// 				}
// 				let r = FACE_ROTATIONS[i].clone();
// 				r.applyAxisAngle(RULE_ORDER[i], orientation*Math.PI/2).round();
// 				return new PolycubePatch(color, undefined, r, i + 1, i + 1 + RULE_ORDER.length);
// 			})
// 		);
// 	});
// }
/**
@deprecated */
export function ruleToHex(rule) {
    const ruleSize = 6;
    let ruleStr = "";
    for (let i=0; i< rule.length; i++) {
        for (let j = 0; j<ruleSize; j++) {
            const face = rule[i].face(j);
            const sign = face.color < 0 ? "1" : "0";
            const color = Math.abs(face.color).toString(2).padStart(5,'0');
            let ori_num = (getSignedAngle(FACE_ROTATIONS[j], face.alignDir, RULE_ORDER[j])*(2/Math.PI)+4)%4
            let ori_str = ori_num.toString(2).padStart(2,'0');
            const binStr = sign + color + ori_str;
            const hexStr = parseInt(binStr,2).toString(16).padStart(2,'0');
            ruleStr += hexStr;
        }
    }
    return ruleStr;
}

interface PatchJson {
	"color": number;
	"dir": {
		"x": number;
		"y": number;
		"z": number;
	};
	"alignDir"?: {
		"x": number;
		"y": number;
		"z": number;
	};
	"orientation"?: string;
	"activation_var"?: number;
	state_var?: number;
}

export interface CubeTypeJson {
	typeName: string;
	type_id: number;
	state?: Array<any>;
	patches: Array<PatchJson>;
	effects?: Array<any>;
	state_size?: number;
}

export function cubeTypeFromJSON(jsonobj: CubeTypeJson, name) {
	// init blank patches in case some are not in JSON file
	let patches = RULE_ORDER.map((d, i)=> {
		return new PolycubePatch(0, d.clone(), FACE_ROTATIONS[i].clone());
	});

	jsonobj.patches.forEach((patchobj, i) => {
		let idx;
		// if the patch direction is specified,
		if ("dir" in patchobj) {
			let patchDir = new THREE.Vector3(
				patchobj.dir.x,
				patchobj.dir.y,
				patchobj.dir.z);
			idx = RULE_ORDER.indexOf(patchDir);
		}
		else {
			idx = i;
		}
		console.assert(patches[idx].color == 0);
		if ("alignDir" in patchobj) {
			patches[idx] = new PolycubePatch(patchobj.color, RULE_ORDER[idx].clone(), new THREE.Vector3(
				patchobj.alignDir.x,
				patchobj.alignDir.y,
				patchobj.alignDir.z));
		}
		else {
			patches[idx] = new PolycubePatch(patchobj.color, undefined, strToVec(patchobj.orientation));
		}
		patches[idx].state_var = patchobj.state_var;
		patches[idx].activation_var = patchobj.activation_var;
	});
	// TODO: import effects
	// let effects = jsonobj.effects if "effects" in jsonobj else [];

	if (("effects" in jsonobj) && ("state_size" in jsonobj)){
		return new DynamicPolycubeCubeType(jsonobj.type_id, patches, name, jsonobj.state_size, jsonobj.effects)
	}
	else {
		return new DynamicPolycubeCubeType(jsonobj.type_id, patches, name, 1, [])
	}
	// else {
	// 	return new StaticPolycubeCubeType(jsonobj.type_id, patches, name);
	// }
}

