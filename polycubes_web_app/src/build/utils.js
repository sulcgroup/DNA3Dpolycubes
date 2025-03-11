import "seedrandom";
import * as THREE from "three";
import { Vector3 } from "three";
import { saveAs } from "file-saver";
/**
 * from stack overflow https://stackoverflow.com/a/32055229
 * @param arr
 * @param newSize
 * @param defaultValue
 */
export function resize_array(arr, newSize, defaultValue = undefined) {
    while (newSize > arr.length)
        arr.push(defaultValue);
    arr.length = newSize;
}
// Modified from https://stackoverflow.com/a/8273091
export function* range(start, stop, step) {
    if (typeof stop === 'undefined') {
        stop = start;
        start = 0;
    }
    if (typeof step === 'undefined') {
        step = 1;
    }
    let iterationCount = 0;
    if (!((step > 0 && start >= stop) || (step < 0 && start <= stop))) {
        for (let i = start; (step > 0 ? i < stop : i > stop); i += step) {
            iterationCount++;
            yield i;
        }
    }
    return iterationCount;
}
// https://stackoverflow.com/questions/4856717/javascript-equivalent-of-pythons-zip-function
// thank you stach overflow for desperately trying to turn JS into a serious language
export const zip = rows => rows[0].map((_, c) => rows.map(row => row[c]));
export function selectColor(number, saturation = 50, value = 65) {
    const hue = number * 137.508; // use golden angle approximation
    return `hsl(${hue % 360},${saturation}%,${value}%)`;
}
/**
 * Emulating https://docs.python.org/3.8/library/stdtypes.html#dict.setdefault
 * If key is in the dictionary, return its value.
 * If not, insert key with a value of default and return default. default defaults to None.
 * @param {*} key
 * @param {*} default_value
 */
Map.prototype.setdefault = function (key, default_value) {
    if (!this.has(key)) {
        this.set(key, default_value);
    }
    return this.get(key);
};
/**
 * custom array equals function
 */
Object.defineProperty(Array.prototype, 'aeq', {
    value: function (a, eq_fn) {
        // loop elements of self
        return a instanceof Array && this.every((v, i) => {
            // if element in this isn't another array
            if (!(v instanceof Array)) {
                return eq_fn == void 0 ? v == a[i] : eq_fn(v, a[i]);
            }
            else {
                return v.aeq(a[i], eq_fn);
            }
        });
    }
});
// permutation code from https://stackoverflow.com/a/30551462
Object.defineProperty(Array.prototype, 'nPr', {
    value: function (r) {
        if (!r)
            return [];
        let a = this;
        return this.reduce(function (memo, cur, i) {
            var others = a.slice(0, i).concat(a.slice(i + 1)), perms = others.nPr(r - 1), newElms = !perms.length ? [[cur]] :
                perms.map(function (perm) { return [cur].concat(perm); });
            return memo.concat(newElms);
        }, []);
    }
});
String.prototype.escapeHTML = function () {
    return this.replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;').replaceAll('"', '&quot;').replaceAll("'", '&#039;');
};
export function getNc(rule) {
    return (new Set([].concat.apply([], rule.map(r => r.map(f => { return Math.abs(f.color); }))))).size - 1;
}
export function getNt(rule) {
    return rule.length;
}
export function vectorAbs(v) {
    return new THREE.Vector3(Math.abs(v.x), Math.abs(v.y), Math.abs(v.z));
}
export function saveString(text, filename) {
    let element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
    element.setAttribute('download', filename);
    element.style.display = 'none';
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
}
export function vecToStr(v) {
    return `(${v.x},${v.y},${v.z})`;
}
export function strToVec(s) {
    let v_components = s.slice(1, -1).split(',').map(x => parseInt(x));
    return new THREE.Vector3(v_components[0], v_components[1], v_components[2]);
}
// function parseKlossString(ruleStr) {
//     return JSON.parse(ruleStr).map(patches => patches.map(function(p) {
//         let color = p[0];
//         let pos = new THREE.Vector3(p[1], p[2], p[3]);
//         let q = new THREE.Quaternion(p[4], p[5], p[6], p[7]);
//         return new Patch(color, pos, q);
//     }));
// }
// function patchySpecToKloss(particlesStr, patchesStr) {
//     // Remove whitespace
//     particlesStr = particlesStr.replaceAll(' ', '');
//     patchesStr = patchesStr.replaceAll(' ', '');
//
//     let getScalar = (name, s) => {
//         let m = s.match(new RegExp(`${name}=(-?\\d+)`));
//         if (m) {
//             return parseFloat(m[1]);
//         }
//         return false
//     }
//     let getArray = (name, s) => {
//         let m = s.match(new RegExp(`${name}=([\\,\\d\\.\\-\\+]+)`));
//         if (m) {
//             return m[1].split(',').map(v=>parseFloat(v));
//         }
//         return false
//     }
//     let particles = [];
//     let currentParticle;
//     for (const line of particlesStr.split('\n')) {
//         const particleID = line.match(/particle_(\d+)/)
//         if (particleID) {
//             if (currentParticle) {
//                 particles.push(currentParticle);
//             }
//             currentParticle = {'id': parseInt(particleID[1])}
//         }
//         const type = getScalar('type', line);
//         if (type !== false) {
//             currentParticle['type'] = type
//         }
//         const patches = getArray('patches', line);
//         if (patches !== false) {
//             currentParticle['patches'] = patches;
//         }
//     }
//     particles.push(currentParticle);
//
//     let patches = new Map();
//
//     let currentId;
//     for (const line of patchesStr.split('\n')) {
//
//         console.log(line)
//         const patchID = line.match(/patch_(\d+)/)
//         if (patchID) {
//             currentId = parseInt(patchID[1]);
//             patches.set(currentId, {});
//         }
//         const color = getScalar('color', line);
//         if (color !== false) {
//             if (Math.abs(color) <= 20) {
//                 console.warn("Self-interactive patches!");
//             }
//             //patches.get(currentId)['color'] = color - Math.sign(color) * 20;
//             patches.get(currentId)['color'] = color;
//         }
//         for (const k of ['position', 'a1', 'a2']) {
//             const a = getArray(k, line);
//             if (a) {
//                 const v = new THREE.Vector3().fromArray(a);
//                 patches.get(currentId)[k] = v;
//             }
//         }
//     }
//
//     for (const particle of particles) {
//         particle['patches'] = particle['patches'].map(id=>patches.get(id));
//     }
//
//     const species = particles.map(particle=>{
//         return particle['patches'].map(patch=>{
//             let p = new Patch(
//                 patch.color,
//                 patch.position,
//                 new THREE.Quaternion()
//             );
//             p.q.copy(rotateVectorsSimultaneously(
//                 p.alignDir, p.dir,
//                 patch.a2, patch.a1
//             ));
//             return p;
//         })
//     });
//
//     return species
// }
// function polycubeRuleToKloss(cubeRule) {
//     let klossRule = [];
//     cubeRule.forEach(species=>{
//         let klossSpecies = [];
//         species.forEach((patch,i)=>{
//             if (patch.color !== 0) {
//                 let p = new Patch(
//                     patch.color,
//                     RULE_ORDER[i].clone().multiplyScalar(0.5),
//                     new THREE.Quaternion()
//                 );
//
//                 p.q.copy(rotateVectorsSimultaneously(
//                     p.alignDir, p.dir,
//                     patch.alignDir, RULE_ORDER[i]
//                 ));
//
//                 console.assert(
//                     RULE_ORDER[i].distanceTo(p.dir) < 1e-4,
//                     `Misdirected polycube rule: ${RULE_ORDER[i].toArray()} !== ${p.dir.toArray()}`
//                 );
//                 console.assert(
//                     patch.alignDir.distanceTo(p.alignDir) < 1e-4,
//                     `Misaligned polycube rule: ${patch.alignDir.toArray()} !== ${p.alignDir.toArray()}`
//                 );
//                 klossSpecies.push(p);
//             }
//         });
//         klossRule.push(klossSpecies);
//     });
//     return klossRule;
// }
export function getSignedAngle(v1, v2, axis) {
    console.assert(v1 instanceof Vector3);
    console.assert(v2 instanceof Vector3);
    console.assert(axis instanceof Vector3);
    let s = v1.clone().cross(v2);
    let c = v1.clone().dot(v2);
    let a = Math.atan2(s.length(), c);
    if (!s.equals(axis)) {
        a *= -1;
    }
    return a;
}
//https://stackoverflow.com/a/55248720
//https://robokitchen.tumblr.com/post/67060392720/finding-a-quaternion-from-two-pairs-of-vectors
export function rotateVectorsSimultaneously(u0, v0, u2, v2) {
    const q2 = new THREE.Quaternion().setFromUnitVectors(u0, u2);
    const v1 = v2.clone().applyQuaternion(q2.clone().conjugate());
    const v0_proj = v0.projectOnPlane(u0);
    const v1_proj = v1.projectOnPlane(u0);
    let angleInPlane = v0_proj.angleTo(v1_proj);
    if (v1_proj.dot(new THREE.Vector3().crossVectors(u0, v0)) < 0) {
        angleInPlane *= -1;
    }
    const q1 = new THREE.Quaternion().setFromAxisAngle(u0, angleInPlane);
    return new THREE.Quaternion().multiplyQuaternions(q2, q1);
}
function getRotationFromSpecies(speciesA, speciesB) {
    const decRuleB = speciesB.toDec();
    for (const r of allRotations()) {
        const q = new THREE.Quaternion().setFromRotationMatrix(new THREE.Matrix4().setFromMatrix3(r));
        if (speciesA.rotate(q).toDec() === decRuleB) {
            return q;
        }
    }
    throw "Matching rotation not found";
}
export function randstr(length = 8, prng) {
    const characters = 'abcdefghijklmnopqrstuvwxyz';
    const charactersLength = characters.length;
    let result = '';
    for (let i = 0; i < length; i++) {
        if (prng) {
            result += characters.charAt(Math.floor(prng() * charactersLength));
        }
        else {
            result += characters.charAt(Math.floor(Math.random() * charactersLength));
        }
    }
    return result;
}
// From stackoverflow/a/12646864
export function shuffleArray(a, rng) {
    for (let i = a.length - 1; i > 0; i--) {
        // let j = Math.floor(Math.random() * (i + 1));
        let j = Math.floor(Math.random() * (i + 1));
        let temp = a[i];
        a[i] = a[j];
        a[j] = temp;
    }
}
export function randOrdering(length, rng) {
    let l = new Array(length);
    for (let i = 0; i < length; i++) {
        l[i] = i;
    }
    shuffleArray(l, rng);
    return l;
}
// this may contain redundancies
export function angle_between_all(vectors, v1) {
    if (v1 == void 0) {
        v1 = vectors[0];
        vectors = vectors.slice(1);
    }
    if (vectors.length == 1) {
        let q = new THREE.Quaternion(); // create one and reuse it
        q.setFromUnitVectors(v1, vectors[0]);
        return [q];
    }
    else {
        return vectors.map(v2 => {
            let q = new THREE.Quaternion(); // create one and reuse it
            q.setFromUnitVectors(v1, v2);
            return q;
        }).concat(angle_between_all(vectors.slice(1), vectors[0]));
    }
}
export function patchCount(cube) {
    return cube.filter(face => face.color != 0).length;
}
function getCenterOfMass(coords) {
    let tot = new THREE.Vector3();
    coords.forEach(c => tot.add(c));
    tot.divideScalar(coords.length);
    return tot;
}
export function rotCoords(coords, r) {
    return coords.map(c => c.clone().applyMatrix3(r));
}
export function coordEqual(coordsA, coordsB) {
    // Copy coords so as not to modify originals
    let cA = coordsA.map(c => c.clone());
    let cB = coordsB.map(c => c.clone());
    if (cA.length != cB.length) {
        console.log(`Coords not equal (different lengths)`);
        return false;
    }
    let comA = getCenterOfMass(cA).round();
    let comB = getCenterOfMass(cB).round();
    cA.forEach(c => c.sub(comA));
    cB.forEach(c => c.sub(comB));
    for (const r of allRotations()) {
        if (compCols(cA, rotCoords(cB, r))) {
            return true;
        }
    }
    console.log(`Coords not equal (Could not find fitting rotation)`);
    return false;
}
// Compare matrices, ignoring column order
export function compCols(coordsA, coordsB) {
    console.assert(coordsA.length == coordsB.length);
    let s1 = new Set(coordsA.map(c => c.toArray().toString()));
    let s2 = new Set(coordsB.map(c => c.toArray().toString()));
    // Check if sets are equal
    if (s1.size !== s2.size)
        return false;
    for (const c of s1) {
        if (!s2.has(c)) {
            return false;
        }
    }
    return true;
}
// function saveString(text, filename) {
//     let element = document.createElement('a');
//     element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
//     element.setAttribute('download', filename);
//     element.style.display = 'none';
//     document.body.appendChild(element);
//     element.click();
//     document.body.removeChild(element);
// }
export function allRotations() {
    return [
        new THREE.Matrix3().set(1, 0, 0, 0, 1, 0, 0, 0, 1),
        new THREE.Matrix3().set(0, -1, 0, 1, 0, 0, 0, 0, 1),
        new THREE.Matrix3().set(-1, 0, 0, 0, -1, 0, 0, 0, 1),
        new THREE.Matrix3().set(0, 1, 0, -1, 0, 0, 0, 0, 1),
        new THREE.Matrix3().set(0, 0, 1, 0, 1, 0, -1, 0, 0),
        new THREE.Matrix3().set(-1, 0, 0, 0, 1, 0, 0, 0, -1),
        new THREE.Matrix3().set(0, 0, -1, 0, 1, 0, 1, 0, 0),
        new THREE.Matrix3().set(1, 0, 0, 0, 0, -1, 0, 1, 0),
        new THREE.Matrix3().set(1, 0, 0, 0, -1, 0, 0, 0, -1),
        new THREE.Matrix3().set(1, 0, 0, 0, 0, 1, 0, -1, 0),
        new THREE.Matrix3().set(0, -1, 0, 0, 0, 1, -1, 0, 0),
        new THREE.Matrix3().set(0, -1, 0, -1, 0, 0, 0, 0, -1),
        new THREE.Matrix3().set(0, -1, 0, 0, 0, -1, 1, 0, 0),
        new THREE.Matrix3().set(0, 0, 1, 1, 0, 0, 0, 1, 0),
        new THREE.Matrix3().set(0, 1, 0, 1, 0, 0, 0, 0, -1),
        new THREE.Matrix3().set(0, 0, -1, 1, 0, 0, 0, -1, 0),
        new THREE.Matrix3().set(0, 0, -1, 0, -1, 0, -1, 0, 0),
        new THREE.Matrix3().set(0, 0, 1, 0, -1, 0, 1, 0, 0),
        new THREE.Matrix3().set(-1, 0, 0, 0, 0, 1, 0, 1, 0),
        new THREE.Matrix3().set(-1, 0, 0, 0, 0, -1, 0, -1, 0),
        new THREE.Matrix3().set(0, 1, 0, 0, 0, -1, -1, 0, 0),
        new THREE.Matrix3().set(0, 1, 0, 0, 0, 1, 1, 0, 0),
        new THREE.Matrix3().set(0, 0, -1, -1, 0, 0, 0, 1, 0),
        new THREE.Matrix3().set(0, 0, 1, -1, 0, 0, 0, -1, 0)
    ];
}
export function signed_axis_angle(Va, Vb, norm) {
    // despise this library.
    var angle = Va.clone().angleTo(Vb);
    var cross = Va.clone().cross(Vb);
    if (norm.clone().dot(cross) < 0) { // Or > 0
        angle = -angle;
    }
    return angle;
}
export function opposite(i) {
    if (i < 0 || i >= RULE_ORDER.length) {
        throw new Error('Index out of bounds');
    }
    return RULE_ORDER.indexOf(RULE_ORDER[i].clone().negate());
}
export function saveArrayBuffer(buffer, filename) {
    saveAs(new Blob([buffer], { type: 'application/octet-stream' }), filename);
}
export const RULE_ORDER = [
    new THREE.Vector3(-1, 0, 0),
    new THREE.Vector3(1, 0, 0),
    new THREE.Vector3(0, -1, 0),
    new THREE.Vector3(0, 1, 0),
    new THREE.Vector3(0, 0, -1),
    new THREE.Vector3(0, 0, 1),
];
RULE_ORDER.indexOf = function (v, fromIndex) {
    return this.findIndex(k => v.equals(k), fromIndex);
};
export const FACE_ROTATIONS = [
    new THREE.Vector3(0, -1, 0),
    new THREE.Vector3(0, 1, 0),
    new THREE.Vector3(0, 0, -1),
    new THREE.Vector3(0, 0, 1),
    new THREE.Vector3(-1, 0, 0),
    new THREE.Vector3(1, 0, 0),
];
// From: https://html-online.com/articles/get-url-parameters-javascript/
export function getUrlVars() {
    const vars = {};
    const searchParams = new URLSearchParams(window.location.search);
    searchParams.forEach((value, key) => {
        vars[key] = value;
    });
    return vars;
}
export function getUrlParam(param, defaultVal) {
    let vars = getUrlVars();
    return param in vars ? vars[param] : defaultVal;
}
function last(arr) { return arr[arr.length - 1]; }
function* numericCombinations(n, r, loc = []) {
    const idx = loc.length;
    if (idx === r) {
        yield loc;
        return;
    }
    for (let next of range(idx ? last(loc) + 1 : 0, n - r + idx)) {
        yield* numericCombinations(n, r, loc.concat(next));
    }
}
export function* combinations(arr, r) {
    for (let idxs of numericCombinations(arr.length, r)) {
        yield idxs.map(i => arr[i]);
    }
}
//# sourceMappingURL=utils.js.map