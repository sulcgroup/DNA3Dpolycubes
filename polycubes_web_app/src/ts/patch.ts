import * as THREE from "three";
import {Quaternion, Vector3} from "three";
import {FACE_ROTATIONS, RULE_ORDER} from "./utils";

export class PolycubePatch {
    color: number; // patch color - determines interaction behavior
    alignDir: Vector3; // alignment of patch - restricts interaction behavior
    direction: number; // orientation of patch. currently specifies index in RULE_ORDER but we would like to do fwd-compatibility w/ klossar
    state_var: number; // binding sets this cube internal state var to true
    activation_var: number; // if this cube inernal state var is true, patch can bind

    constructor(color: number, direction: number | THREE.Vector3, alignment: THREE.Vector3 | number | undefined = undefined, state_var: number = 0, activationvar: number = 0) {
        this.color = color;
        if (direction instanceof THREE.Vector3) {
            direction = RULE_ORDER.indexOf(direction);
        }
        if (typeof alignment == "number") {
            let default_alignment = FACE_ROTATIONS[direction].clone();
            default_alignment.applyAxisAngle(RULE_ORDER[direction],
                alignment * Math.PI / 2).round();
            alignment = default_alignment;
        }
        this.alignDir = alignment;
        console.assert(typeof direction == "number")
        this.direction = direction;
        let dir_vector: Vector3 = RULE_ORDER[direction];
        if (this.has_torsion()) {
            console.assert(this.alignDir.dot(dir_vector) < 1e-6)
        }
        this.state_var = state_var;
        this.activation_var = activationvar;
    }

    rotate(q: Quaternion) : PolycubePatch {
        let newFaceAlign = this.alignDir.clone().applyQuaternion(q).round();
        let newFaceDir = RULE_ORDER[this.direction].clone().applyQuaternion(q).round();
        return new PolycubePatch(this.color, newFaceDir, newFaceAlign, this.state_var, this.activation_var);
    }

    has_color() {
        return this.color != 0;
    }

    has_torsion() {
        return typeof this.alignDir != "undefined";
    }

    abscolor() {
        return Math.abs(this.color);
    }
}