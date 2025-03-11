import "./init";
import {parseHexRule} from "./rule";

let rotation = 0;
rotateCamera();

regenerate()

window.system.setMaxCubes(200);

// Regenerate when out of bounds (nMaxCubes reached)
window.addEventListener('oub', function(e) {
    regenerate();
}, false);

// Regenerate when there are no more cubes to add
window.addEventListener('movesProcessed', async function(e) {
    await sleep(2000);
    regenerate();
}, false);

function rotateCamera() {
    rotation += 0.005;
    window.camera.position.y = Math.sin(rotation) * 20;
    window.camera.position.x = Math.cos(rotation) * 20;
    window.camera.lookAt(window.system.centerOfMass);
    requestAnimationFrame(rotateCamera.bind(this));
}

function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

function regenerate() {
    let maxRuleSize = 8;
    let ruleSize = Math.round(Math.random()*maxRuleSize)+1;
    let hexRule = "";
    while(ruleSize--) {
        hexRule += (Math.abs(Math.random()*0xFFFFFFFFFFFF<<0)).toString(16);
    }
    let argstr = "?hexRule="+hexRule;
    window.history.pushState(null, null, argstr);
    
    // Parse rule
    let rules = parseHexRule(hexRule);
    window.system.set_rule(rules);
}
