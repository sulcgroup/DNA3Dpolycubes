// import * as THREE from "three";
// import {VRButton} from "three/examples/jsm/webxr/VRButton";
// import {parseHexRule} from "./rule";
// let rotation = 0;
//
// let rig = new THREE.PerspectiveCamera();
// rig.add(window.camera);
// window.scene.add(rig);
// rig.lookAt(window.system.centerOfMass);
//
// document.body.appendChild(VRButton.createButton(window.renderer));
// window.renderer.vr.enabled = true;
// window.renderer.setAnimationLoop(function(){
//     rotation += 0.001;
//     rig.position.x = Math.sin(rotation) * 5;
//     rig.position.z = Math.cos(rotation) * 5;
//     rig.lookAt(new THREE.Vector3(0,0,0));
//     window.renderer.render(window.scene, window.camera);
// });
//
// const controller = (window.renderer.vr as any).getController(0);
//
// const selectListener = (event) => {
//     if(firstRule) { // Ignore first click
//         firstRule = false;
//         return;
//     }
//     let maxRuleSize = 20;
//     let ruleSize = Math.round(Math.random()*maxRuleSize)+1;
//     let hexRule = "";
//     while(ruleSize--) {
//         hexRule += (Math.abs(Math.random()*0xFFFFFFFFFFFF<<0)).toString(16);
//     }
//
//     // Parse rule
//     let rules = parseHexRule(hexRule);
//     window.system.set_rule(rules);
// };
//
// controller.addEventListener('select', selectListener);
// let firstRule = true;
//
//
//
