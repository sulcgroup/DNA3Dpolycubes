import { FACE_ROTATIONS, getSignedAngle, RULE_ORDER, selectColor, strToVec } from "./utils";
import * as THREE from "three";
import { PolycubePatch } from "./patch";
export const FACE_NAMES = ["left", "right", "bottom", "top", "back", "front"];
export function cubeTypeColor(i) {
    return selectColor(i);
}
export function patchColor(i) {
    return selectColor(Math.abs(i) - 1);
}
export class PolycubeCubeType {
    typeName;
    type_id;
    patches;
    // cube instance vars
    _state_size;
    constructor(typeid, patches, type_name, stateSize = 0) {
        console.assert(typeof typeid == "number");
        console.assert(patches.length == 6);
        this.type_id = typeid;
        if (typeof type_name == "undefined") {
            type_name = `CT${typeid}`;
        }
        this.typeName = type_name;
        this.patches = patches;
        // add 1 to state size to incorporate tautology variable
        this._state_size = stateSize;
    }
    /**
     * copies this cube type
     * @param newName
     * @param newtype
     */
    clone(newtype) {
        return new PolycubeCubeType(newtype, this.patches.map(patch => {
            return new PolycubePatch(patch.color, patch.direction, patch.alignDir, patch.state_var, patch.activation_var);
        }), `CT${newtype}`, this.state_size());
    }
    /**
     *
     * @returns {number} number of state variables (not counting tautology variable)
     */
    state_size() {
        return this._state_size;
    }
    // // https://stackoverflow.com/a/25199671
    // /**
    //  * rotates cube type. utility of this. unclear.
    //  * @param q
    //  */
    rotate(q) {
        console.assert(q !== void 0);
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
    hasAllostery() {
        return false;
    }
    patch_color(i) {
        return this.patches[i].color;
    }
    face(i) {
        return this.patches[i];
    }
    toDec() {
        return this.patches.map((patch, i) => {
            if (patch.color === 0) {
                return '';
            }
            else {
                let theta = getSignedAngle(FACE_ROTATIONS[i], patch.alignDir, RULE_ORDER[i]);
                let orientation = Math.round(theta * (2 / Math.PI) + 4) % 4;
                if (orientation == 0 && patch.state_var == 0 && patch.activation_var == 0) {
                    return `${patch.color}`;
                }
                else if (patch.state_var == 0 && patch.activation_var == 0) {
                    return `${patch.color}:${orientation}`;
                }
                else {
                    return `${patch.color}:${orientation}:${patch.state_var}:${patch.activation_var}`;
                }
            }
        }).join('|');
    }
    num_colored_faces() {
        return this.patches.map(patch => patch.color).reduce((count, color) => count + Math.min(Math.abs(color), 1));
    }
    colored_faces_idxs() {
        return this.patches.filter(p => p.color).map(p => p.direction);
    }
    patch_colors() {
        return new Set(this.patches.filter(patch => patch.color != 0).map(patch => patch.color));
    }
    getPatchCoords() {
        return this.patches.filter(patch => {
            return patch.color !== 0;
        }).map(patch => RULE_ORDER[patch.direction]);
    }
    patch_allosterically_controlled(n) {
        throw Error();
    }
    reeval_conditionals(state) {
        return [[], []];
    }
    /**
     *
     * @param other
     */
    equals(other) {
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
            }).reduce((a, b) => a && b);
    }
}
export class DynamicPolycubeCubeType extends PolycubeCubeType {
    effects;
    constructor(type_id, patches, name, state_size = 0, effects = []) {
        super(type_id, patches, name, state_size);
        this.effects = effects;
    }
    hasAllostery() {
        return this.patches.some((p, i) => {
            return this.patch_allosterically_controlled(i);
        });
    }
    patch_allosterically_controlled(i) {
        return this.patches[i].activation_var != 0;
    }
    /**
     * reevaluates state passed from a cube type instance
     * @param state
     * @return a tuple where the first element is CHANGES (each element is true if the reevaluation has changed the state var) and the second is the NEW STATE
     */
    reeval_conditionals(state) {
        let change;
        let newState = [...state];
        do {
            change = false;
            this.effects.forEach(effect => {
                let can_fire = effect.sources.every(s => newState[s]);
                if (can_fire && !newState[effect.target]) {
                    change = true;
                    newState[effect.target] = true;
                }
            });
        } while (change);
        let changes = state.map((x, i) => {
            return x !== newState[i];
        });
        state = newState;
        return [changes, state];
    }
    clone(newtype) {
        let copy = toDynamic(super.clone(newtype));
        copy.effects = this.effects.map((e) => {
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
        this.patches.forEach(p => {
            new_state_size = Math.max(new_state_size, p.state_var);
            new_state_size = Math.max(new_state_size, Math.abs(p.activation_var));
        });
        this._state_size = Math.max(new_state_size, ...this.effects.map(e => { return e.target; }));
        // resize_array(this.state, new_state_size);
    }
    toDec() {
        let decstr = super.toDec();
        if (this.effects.length > 0) {
            decstr += "@" + this.effects.map(e => {
                return `[${e.sources.join(",")}]>${e.target}`;
            }).join(";");
        }
        return decstr;
    }
}
function toDynamic(c) {
    let newCube = new DynamicPolycubeCubeType(c.type_id, c.patches, c.typeName, c.state_size());
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
export function ruleToDec(cubeTypes) {
    if (cubeTypes.length > 0) {
        return cubeTypes.map(ct => ct.toDec()).join('_');
    }
    else {
        return "";
    }
}
/**
 */
export function parseDecRule(rulesStr) {
    // split by underscore
    return rulesStr.split('_').map((s, j) => {
        let patch_conditionals = RULE_ORDER.map(x => "");
        var state_size = 1;
        let delimeter;
        if (s.match(/^([^|]*\|){5}[^|]*@?.*?$/)) {
            delimeter = "|";
        }
        else {
            delimeter = "#";
        }
        var face_props;
        let patches = s.split("@")[0].split(delimeter).map((face, i) => {
            let state_var = "0";
            let activation_var = "0";
            let i_color = 0;
            let r = FACE_ROTATIONS[i].clone();
            let i_orientation = 0;
            let logic = "";
            if (!(face.trim() == "")) {
                // @ts-ignore
                face_props = face.split(":");
                let color = face_props[0];
                if (face_props.length > 1) {
                    let patchOriIdx = face_props[1];
                    // color, orientation, state var, activation var
                    if (face_props.length == 4) {
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
                            }
                            catch (e) {
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
                if (isNaN(i_color)) {
                    throw "Malformed rule";
                }
                r.applyAxisAngle(RULE_ORDER[i], i_orientation * Math.PI / 2).round();
            }
            return new PolycubePatch(i_color, i, r, i_state_var, i_activation_var);
        });
        console.assert(patches.length == RULE_ORDER.length);
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
        let ct = new DynamicPolycubeCubeType(j, patches, `CT${j}`, state_size, effects);
        ct.update_state_size();
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
    for (let i = 0; i < ruleStr.length; i += 2 * ruleSize) {
        let cube_type = [];
        //console.log("Rule ",(i/(2*ruleSize))+1);
        for (let j = 0; j < ruleSize; j++) {
            let face = ruleStr.substring(i + (2 * j), i + (2 * j) + 2);
            let binStr = (parseInt(face, 16).toString(2)).padStart(8, '0');
            let sign = parseInt(binStr[0], 2);
            let color = parseInt(binStr.substring(1, 6), 2);
            let orientation = parseInt(binStr.substring(6, 8), 2);
            let r = FACE_ROTATIONS[j].clone();
            r.applyAxisAngle(RULE_ORDER[j], orientation * Math.PI / 2);
            r.round();
            cube_type.push({ 'color': color * (sign ? -1 : 1), 'alignDir': r });
        }
        rule.push(cube_type);
    }
    //added this line of code to switch to new method of storing rules
    rule = rule.map((x, i) => {
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
    for (let i = 0; i < rule.length; i++) {
        for (let j = 0; j < ruleSize; j++) {
            const face = rule[i].face(j);
            const sign = face.color < 0 ? "1" : "0";
            const color = Math.abs(face.color).toString(2).padStart(5, '0');
            let ori_num = (getSignedAngle(FACE_ROTATIONS[j], face.alignDir, RULE_ORDER[j]) * (2 / Math.PI) + 4) % 4;
            let ori_str = ori_num.toString(2).padStart(2, '0');
            const binStr = sign + color + ori_str;
            const hexStr = parseInt(binStr, 2).toString(16).padStart(2, '0');
            ruleStr += hexStr;
        }
    }
    return ruleStr;
}
export function cubeTypeFromJSON(jsonobj, name) {
    // init blank patches in case some are not in JSON file
    let patches = RULE_ORDER.map((d, i) => {
        return new PolycubePatch(0, d.clone(), FACE_ROTATIONS[i].clone());
    });
    jsonobj.patches.forEach((patchobj, i) => {
        let idx;
        // if the patch direction is specified,
        if ("dir" in patchobj) {
            let patchDir = new THREE.Vector3(patchobj.dir.x, patchobj.dir.y, patchobj.dir.z);
            idx = RULE_ORDER.indexOf(patchDir);
        }
        else {
            idx = i;
        }
        console.assert(patches[idx].color == 0);
        if ("alignDir" in patchobj) {
            patches[idx] = new PolycubePatch(patchobj.color, RULE_ORDER[idx].clone(), new THREE.Vector3(patchobj.alignDir.x, patchobj.alignDir.y, patchobj.alignDir.z));
        }
        else {
            patches[idx] = new PolycubePatch(patchobj.color, undefined, strToVec(patchobj.orientation));
        }
        patches[idx].state_var = patchobj.state_var;
        patches[idx].activation_var = patchobj.activation_var;
    });
    // TODO: import effects
    // let effects = jsonobj.effects if "effects" in jsonobj else [];
    if (("effects" in jsonobj) && ("state_size" in jsonobj)) {
        return new DynamicPolycubeCubeType(jsonobj.type_id, patches, name, jsonobj.state_size, jsonobj.effects);
    }
    else {
        return new DynamicPolycubeCubeType(jsonobj.type_id, patches, name, 1, []);
    }
    // else {
    // 	return new StaticPolycubeCubeType(jsonobj.type_id, patches, name);
    // }
}
//# sourceMappingURL=rule.js.map