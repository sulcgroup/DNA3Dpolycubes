// import * as THREE from "three";
// import {onWindowResize, render} from "./view";
// import "./init";
// import {CustomTransformControls} from "./libs/controls/TransformControls";
//
// window.transform = new CustomTransformControls(window.camera, window.canvas);
// window.transform.setSpace('local');
// window.transform.addEventListener('change', ()=>{
//     render();
// });
//
// window.transform.addEventListener('dragging-changed', event=>{
//     window.orbit.enabled = !event.value;
//     if (!event.value && window.transform.object) {
//         window.transform.object.patch.update(undefined,
//             window.transform.object.position,
//             window.transform.object.quaternion);
//         window.system.set_rule([]);
//     }
// });
//
// window.scene.add(window.transform);
//
// window.raycaster = new THREE.Raycaster();
// window.mouse = new THREE.Vector2()
// function onMouseMove(event) {
// 	window.mouse.x = ( event.clientX / window.innerWidth ) * 2 - 1;
// 	window.mouse.y = - ( event.clientY / window.innerHeight ) * 2 + 1;
//     window.requestAnimationFrame(rayrender);
// }
//
// const hoverscale = 1.2;
//
// function rayrender() {
// 	// update the picking ray with the camera and mouse position
// 	window.raycaster.setFromCamera( window.mouse, window.camera );
//
//     if (window.orbit.enabled) {
//         // calculate objects intersecting the picking ray
//         const intersects = window.raycaster.intersectObjects(window.system.getCubeObjects(), true);
//
//         if (intersects.length > 0) {
//             let patchGroup = intersects[0].object.parent;
//
//             if (window.transform.object) {
//                 window.transform.object.children.forEach(c=>{
//                     c.scale.divideScalar(hoverscale);
//                 });
//             }
//
//             window.transform.attach(patchGroup);
//
//             patchGroup.children.forEach(c=>{
//                 c.scale.multiplyScalar(hoverscale);
//             });
//         }
//     }
//
// 	window.renderer.render( window.scene, window.camera );
// }
//
// window.canvas.addEventListener('mousedown', ()=>{
//     window.canvas.focus();
//     console.log("Canvas has focus")
//     if (window.orbit.enabled && window.transform.object) {
//         window.transform.object.children.forEach(c=>{
//             c.scale.divideScalar(hoverscale);
//         });
//         window.transform.detach();
//     }
// }, false)
//
// window.canvas.addEventListener('mousemove', onMouseMove, false);
//
// window.canvas.addEventListener('keydown', event=>{
//
//     console.log(event.key)
//
//     switch (event.key) {
//
//         case 'q':
//             window.transform.setSpace(window.transform.space === 'local' ? 'world' : 'local');
//             break;
//
//         case 'Shift':
//             window.transform.setTranslationSnap(100);
//             window.transform.setRotationSnap(THREE.MathUtils.degToRad(15));
//             break;
//
//         case 't':
//             window.transform.setMode('translate');
//             break;
//
//         case 'r':
//             window.transform.setMode('rotate');
//             break;
//
//         // deprecated bc it's not typescript compatible and frankly i don't care
//         // case 'c':
//         //     const position = camera.position.clone();
//         //
//         //     camera = camera.isPerspectiveCamera ? cameraOrtho : cameraPersp;
//         //     camera.position.copy(position);
//         //
//         //     orbit.object = camera;
//         //     transform.camera = camera;
//         //
//         //     camera.lookAt(orbit.target.x, orbit.target.y, orbit.target.z);
//         //     onWindowResize();
//         //     break;
//
//         case 'v':
//             const randomFoV = Math.random() + 0.1;
//             const randomZoom = Math.random() + 0.1;
//
//             window.cameraPersp.fov = randomFoV * 160;
//             window.cameraOrtho.bottom = -randomFoV * 5;
//             window.cameraOrtho.top = randomFoV * 5;
//
//             window.cameraPersp.zoom = randomZoom * 5;
//             window.cameraOrtho.zoom = randomZoom * 5;
//             onWindowResize();
//             break;
//
//         case '+':
//             window.transform.setSize(window.transform.size + 0.1);
//             break;
//
//         case '-':
//             window.transform.setSize(Math.max(window.transform.size - 0.1, 0.1));
//             break;
//
//         case 'x':
//             window.transform.showX = !window.transform.showX;
//             break;
//
//         case 'y':
//             window.transform.showY = !window.transform.showY;
//             break;
//
//         case 'z':
//             window.transform.showZ = !window.transform.showZ;
//             break;
//
//         case ' ':
//             window.transform.enabled = !window.transform.enabled;
//             break;
//
//         case 'Escape':
//             window.transform.reset();
//             break;
//     }
//
//     const stepAngle = Math.PI/12;
//     switch (event.code) {
//         case 'Numpad0':
//             window.orbit.reset();
//             break;
//         case 'Numpad1':
//             if (event.ctrlKey || event.metaKey) {
//                 window.orbit.setToAxis(new THREE.Vector3(-1, 0, 0));
//                 break;
//             }
//             else {
//                 window.orbit.setToAxis(new THREE.Vector3(1, 0, 0));
//                 break;
//             }
//         case 'Numpad2':
//             window.orbit.stepAroundAxis(new THREE.Vector3(-1, 0, 0), stepAngle);
//             break;
//         case 'Numpad3':
//             if (event.ctrlKey || event.metaKey) {
//                 window.orbit.setToAxis(new THREE.Vector3(0, -1, 0));
//                 break;
//             }
//             else {
//                 window.orbit.setToAxis(new THREE.Vector3(0, 1, 0));
//                 break;
//             }
//         case 'Numpad4':
//             window.orbit.stepAroundAxis(new THREE.Vector3(0, 1, 0), stepAngle);
//             break;
//         // case 'Numpad5':
//         //     api.switchCamera();
//         //     break;
//         case 'Numpad6':
//             window.orbit.stepAroundAxis(new THREE.Vector3(0, -1, 0), stepAngle);
//             break;
//         case 'Numpad7':
//             if (event.ctrlKey || event.metaKey) {
//                 window.orbit.setToAxis(new THREE.Vector3(0, 0, -1));
//                 break;
//             }
//             else {
//                 window.orbit.setToAxis(new THREE.Vector3(0, 0, 1));
//                 break;
//             }
//         case 'Numpad8':
//             window.orbit.stepAroundAxis(new THREE.Vector3(1, 0, 0), stepAngle);
//             break;
//         case 'Numpad9':
//             if (event.ctrlKey || event.metaKey) {
//                 window.orbit.setToAxis(new THREE.Vector3(0, 0, 1));
//                 break;
//             }
//             else {
//                 window.orbit.setToAxis(new THREE.Vector3(0, 0, -1));
//                 break;
//             }
//     }
// });
//
// window.canvas.addEventListener('keyup', event=>{
//     if (event.shiftKey) {
//         window.transform.setTranslationSnap(null);
//         window.transform.setRotationSnap(null);
//         window.transform.setScaleSnap(null);
//     }
// });