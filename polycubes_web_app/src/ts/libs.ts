// file for useful high-level functions
// avoid importing from this file in any but even higher-level functions!!!!

import {saveAs} from "file-saver";

import {getUrlVars, saveArrayBuffer, saveString} from "./utils";
import {GLTFExporter} from "three/examples/jsm/exporters/GLTFExporter";

// import {addSpecies, log, update} from "./edit";

// function getCubeTypeCount(rule, assemblyMode = 'seeded', maxPieceSize = 3) {
//     let test_system = new PolycubeSystem(rule,
//         100,
//         100,
//         false,
//         true,
//          1,
//          getRandSeed());
//     test_system.seed();
//     let processed: boolean | string = false;
//     while (test_system.moveKeys.length > 0 && test_system.cubeMap.size < test_system.nMaxCubes) {
//         processed = test_system.step();
//         if (processed == 'oub') {
//             console.warn("Getting cube type count for unbounded rule");
//             break;
//         }
//     }
//     return test_system.cubeTypeCount;
// }

//should REALLY REALLY have nTries > 15
//may also be wise to have step limit greater than 100
//really this whole function needs to be rewritten to actually prove whether arbitrary
//systems are bounded and deterministic
//it's possible that I, as a non-mathmetician, have just casually suggested solving a
//millenium-prize level mathematical question in this javascript applet
// export function isBoundedAndDeterministic(rule, nTries = 15, assemblyMode = 'seeded') {
//     let oldCoords;
//     while (nTries--) {
//         let system = new PolycubeSystem(rule, undefined, 100, 100, assemblyMode);
//         system.seed();
//         let processed: boolean | string = false;
//         while (!processed) {
//             processed = system.step();
//             if (processed == 'oub') {
//                 return '∞';
//             }
//         }
//         let strCoords = [...system.cubeMap.values()];
//         if (oldCoords && !coordEqual(oldCoords, strCoords)) {
//             return '?';
//         }
//         oldCoords = strCoords;
//     }
//     return true;
// }

// export function simplify2(rules, onUpdate) {
//     rules = simplify(rules);
//
//     let correctCoords = getCoords(rules);
//
//     const rotations = allRotations().map(
//         r => [r, new THREE.Quaternion().setFromRotationMatrix(
//             new THREE.Matrix4().setFromMatrix3(r))]
//     );
//
//     let alreadyTried = new Set();
//     let triedStr = (a, b) => `${[a, b].sort()}`;
//
//     let updatedRuleset = true;
//     while (updatedRuleset) {
//         updatedRuleset = false;
//         for (const [iA, cA] of rules.entries()) {
//             const coordsA = getPatchCoords(cA);
//             for (const [iB, cB] of rules.entries()) {
//                 if (iA < iB) {
//                     const coordsB = getPatchCoords(cB);
//                     // If they have the same number of patches
//                     if (coordsA.length === coordsB.length && !alreadyTried.has(triedStr(iA, iB))) {
//                         for (const [r, q] of rotations) {
//                             if (compCols(coordsA, rotCoords(coordsB, r))) {
//                                 console.log(`Cube type ${iA} is similar to ${iB}`);
//                                 let colorMap = new Map();
//                                 //let rotMap = new Map();
//                                 //const q = new THREE.Quaternion().setFromRotationMatrix(new THREE.Matrix4().setFromMatrix3(r));
//                                 const rotatedB = cB.rotate(q);
//                                 for (let i = 0; i < 6; i++) {
//                                     if (rotatedB.colors[i] !== 0) {
//                                         console.assert(cA.colors[i] !== 0);
//                                         colorMap.set(rotatedB.colors[i], cA.colors[i]);
//                                         //rotMap.set(rotatedB[i].color, getSignedAngle(cB[i].alignDir, cA[i].alignDir, RULE_ORDER[i]));
//                                     }
//                                 }
//                                 let newRules = [];
//                                 // Clone rule and create new ruleset
//                                 rules.map(r => r.clone(false)).forEach((c, i) => {
//                                     if (i !== iB) {
//                                         c.toOld().forEach((p, i) => {
//                                             // Replace all of the old color with the new
//                                             const color = p.color;
//                                             if (colorMap.has(color)) {
//                                                 c.colors[i] = colorMap.get(color);
//                                                 c.alignments[i] = cA.alignments[i];
//                                             } else if (colorMap.has(-color)) {
//                                                 c.colors[i] = -colorMap.get(-color);
//                                                 p.alignDir.applyQuaternion(q.invert());
//                                                 c.alignments[i] = p.alignDir;
//                                             }
//                                         });
//                                     } else {
//                                         c.toOld().forEach((p, i) => {
//                                             c.colors[i] = 0
//                                             c.alignments[i] = FACE_ROTATIONS[i];
//                                         })
//                                     }
//                                     newRules.push(c);
//                                 });
//
//                                 const simplifiedNewRule = simplify(newRules); //remove empty slots
//
//                                 // Accept new ruleset if we still get the same chape
//                                 if (
//                                     (isBoundedAndDeterministic(simplifiedNewRule, undefined, 'seeded') === true) &&
//                                     coordEqual(getCoords(simplifiedNewRule, 'seeded'), correctCoords)
//                                 ) {
//                                     updatedRuleset = true;
//                                     rules = newRules;
//                                     console.log(`Changing  ${iB} to ${iA} did work: ${ruleToDec(simplifiedNewRule)}`);
//                                     alreadyTried.add(triedStr(iA, iB));
//                                     if (onUpdate) {
//                                         onUpdate(ruleToDec(simplifiedNewRule));
//                                     }
//                                     break;
//                                 } else {
//                                     console.log(`Changing  ${iB} to ${iA} did not work: ${ruleToDec(simplifiedNewRule)}\nBounded & Deterministic = ${isBoundedAndDeterministic(simplifiedNewRule)}\nEqual coords = ${coordEqual(getCoords(simplifiedNewRule), correctCoords)}`);
//                                     alreadyTried.add(triedStr(iA, iB));
//                                     //break;
//                                 }
//                             }
//                         }
//                     }
//                 }
//                 if (updatedRuleset) {
//                     break;
//                 }
//             }
//             if (updatedRuleset) {
//                 break;
//             }
//         }
//     }
//     return simplify(rules);
// }

// function reduce(r1start = 0, r2start = 0) {
//     //TODO: verify these sort functions
//     let rules_ordered_by_free_space = [...window.system.cube_types].sort(
//         (a, b) => a.num_colored_faces() - b.num_colored_faces()
//     )
//     // loop rules ny first checking the rule with the most free space against the one with the most colors, and contining as such
//     for (let i = r1start; i < window.system.cube_types.length; i++) {
//         for (let j = r2start; j < window.system.cube_types.length; j++) {
//             let r1 = rules_ordered_by_free_space[i]
//             let r2 = rules_ordered_by_free_space[window.system.cube_types.length - j - 1]
//             if (r1 == r2 || r1.num_colored_faces() > RULE_ORDER.length - r2.num_colored_faces()) {
//                 continue;
//             }
//             let r1_faces = r1.colored_faces_idxs()
//             let r2_faces = [...RULE_ORDER.keys].filter(i => r2.colors[i] == 0)
//             for (let i_fc = Math.min(r1.num_colored_faces(), RULE_ORDER.length - r2.num_colored_faces()); i_fc >= 2; i_fc--) {
//                 let r1_faceperm = r1_faces.nPr(i_fc)
//                 let r2_faceperm = r2_faces.nPr(i_fc)
//                 // |r1_facep| should equal |r2_facep|
//                 for (let i_perm = 0; i_perm < r1_faceperm.length; i_perm++) {
//                     let perm1 = r1_faceperm[i_perm]
//                     let perm2 = r2_faceperm[i_perm]
//                     // if the angle between the empty faces on r2 is the same as
//                     // the angle between the colored faces on r1...
//                     if (angle_between_all(perm1.map(fidx => RULE_ORDER[fidx])).equals(angle_between_all(perm2.map(fidx => RULE_ORDER[fidx])), (a, b) => a.equals(b))) {
//                         let r2_incumbants = r2.colors.keys().filter(x => !(x in r2))
//                         // replace the empty faces on r2 with the colored faces on r1
//                         for (let fidx = 0; fidx < i_fc; fidx++) {
//                             let fidx1 = perm1[fidx];
//                             let fidx2 = perm2[fidx];
//                             r2.colors[fidx2] = r1.colors[fidx1]
//                             r2.alignments[fidx2] = r1.alignments[fidx1].clone()
//                             // create a new conditional so that the patches currently on r1 will not
//                             // be active if any of the patches currently on r2 are active
//                             let new_conditional = "(" + r2_incumbants.map(idx => `!p[${idx}]`).join(" & ") + ")";
//                             // if the patch being copied
//                             if (r1.conditionals[fidx1] != "(true)") {
//                                 new_conditional = `(${new_conditional} & ${r1.conditionals[fidx1]})`
//                             }
//                             (r2 as StaticPolycubeCube).conditionals[fidx2] = new_conditional;
//                         }
//
//                         // update conditionals on faces previously on r2 so that they don't
//                         // cross-talk with r1
//
//                         r2_incumbants.forEach(fidx => {
//                             if (!(fidx in perm2)) {
//                                 // set this patch to inactivate if any binding occurs on any of the new patches
//                                 let new_conditional = "(" + perm2.map(idx=>[`!p[${idx}]`]).join(" & ") + ")";
//                                 // if the incumbent patch already has a conditional
//                                 if ((r2 as StaticPolycubeCube).conditionals[fidx] != "(true)") {
//                                     new_conditional = `(${new_conditional} & ${(r2 as StaticPolycubeCube).conditionals[fidx]})`
//                                 }
//                                 (r2 as StaticPolycubeCube).conditionals[fidx] = new_conditional;
//                             }
//                         });
//                         // remove r1 from system
//                         window.system.cube_types = window.system.cube_types.splice(i, i + 1)
//                         // continue reducing within new function call to avoid
//                         // messing with indexing
//                         return reduce(i, j);
//                     }
//                 }
//             }
//         }
//     }
// }

/**
 * Computes a random seed
 * if if the random seed is specified in the url, use that. otherwise, find a random seed from the document
 */
export function getRandSeed(fromDoc=false, seedlen=32) {
    if ("randseed" in getUrlVars()) {
        return getUrlVars()["randseed"];
    } else{
        let randseed;
        if (fromDoc) {
            let rstart = Math.round(Math.random() * document.documentElement.innerHTML.length - seedlen);
            let randseed = document.documentElement.innerHTML.substring(rstart, rstart + seedlen);
        }
        else {
            const characters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
            randseed = [...Array(seedlen)].map(()=>characters[Math.floor(Math.random() * characters.length)]).join("");
        }
        console.log(`Seeded RNG with string "${randseed.escapeHTML()}"`);
        return randseed;
    }

}
// function assessVariability(num_replications, fname="stats"){
//     if (num_replications == void 0){
//         num_replications = Math.pow(2, window.system.cube_types.length); //made-up number
//     }
//     let shape_types = [];
//     let shape_names = [];
//     let shape_counts = [];
//     for (let i = 0; i < num_replications; i++)
//     {
//         let test_model = window.system.clone();
//         test_model.assemblyMode = 'stochastic'
//         test_model.resetRules(test_model.cube_types, false);
//         seedrandom(i);
//         test_model.seed();
//         while (test_model.moveKeys.length > 0 && test_model.cubeMap.size < test_model.nMaxCubes) {
//             test_model.step();
//         }
//         test_model.locmap = test_model.getLocationMap();
//         let idx = shape_types.findIndex(s=>s.isEquivelent(test_model));
//         if (idx > -1){
//             shape_counts[idx]++;
//         }
//         else {
//             shape_types.push(test_model);
//             shape_names.push(`Shape${shape_types.length}`);
//             shape_counts.push(1);
//         }
//     }
//     let json = {};
//     json['rules'] = window.system.cube_types;
//     json['shape_types'] = [];
//     for (let idx = 0; idx < shape_types.length; idx++){
//         json['shape_types'].push({
//             'name': shape_names[idx],
//             'shape': Object.fromEntries(shape_types[idx].locmap),
//             'occurrances': shape_counts[idx]
//         });
//     }
//     /**
//      stolen from
//      https://ourcodeworld.com/articles/read/189/how-to-create-a-file-and-generate-a-download-with-javascript-in-the-browser-without-a-server
//      because incredibly, there is no built-in javascript function for this*/
//     json = JSON.stringify(json, null, 4);
//     var element = document.createElement('a');
//     element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(json));
//     element.setAttribute('download', `${fname}.json`);
//
//     element.style.display = 'none';
//     document.body.appendChild(element);
//
//     element.click();
//
//     document.body.removeChild(element);
// }

function exportGLTF(name = 'scene') {
    // Instantiate an exporter
    let exporter = new GLTFExporter();
    let options = { 'forceIndices': true };
    // Parse the input and generate the glTF output
    exporter.parse(window.system.objGroup, function (result) {
        if (result instanceof ArrayBuffer) {
            saveArrayBuffer(result, 'scene.glb');
        }
        else {
            let output = JSON.stringify(result, null, 2);
            console.log(output);
            saveString(output, `${name}.gltf`);
        }
    }, () => { }, options);
}

export function downloadStructure(){
    // let rule = JSON.stringify(window.system.cube_types, null, 4);
    let struct = window.system.export() as BlobPart;
    saveAs(new Blob(
        [struct],
        {type: "text/plain;charset=utf-8"}), "ruleset.json");
}
