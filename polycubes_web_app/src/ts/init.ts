import * as $ from "jquery";
import {FACE_ROTATIONS, getUrlParam, getUrlVars, RULE_ORDER} from "./utils";
import {DynamicPolycubeCubeType, parseDecRule, parseHexRule} from "./rule";
import  {VisualPolycubeSystem} from "./view";
import {getRandSeed} from "./libs";
import {PolycubePatch} from "./patch";

declare global {
    export interface Window {
        system: VisualPolycubeSystem;
        canvas: HTMLCanvasElement;
    }
}
export function createPolycubeSystem() : VisualPolycubeSystem {
    let cube_types = {};

    // Parse rule
    let vars = getUrlVars();
	if ("rule" in vars) {
        cube_types = parseDecRule(vars["rule"]);
    } else if ("hexRule" in vars) {
        cube_types = parseHexRule(vars["hexRule"]);
    } else if ("decRule" in vars) {
        cube_types = parseDecRule(decodeURIComponent(vars["decRule"]));
    } else {
        let defaultRule = [[1,1,1,1,1,1],[-1,0,0,0,0,0]];
        cube_types = defaultRule.map((cube_type, i) => {
            return new DynamicPolycubeCubeType(i,
                cube_type.map((c, j) => {
                    return new PolycubePatch(c, j, FACE_ROTATIONS[j], 0, 0);
            }),
                );
        });
    }

    let nMaxCubes = parseInt(getUrlParam("nMaxCubes", 100));
    let maxCoord = parseInt(getUrlParam("maxCoord", 100));

    let system = new VisualPolycubeSystem(
        window.scene,
        nMaxCubes,
        maxCoord,
        true,
        true,
        $("#env-browser"), 3,
                getRandSeed());
    system.set_rule(cube_types);
    // check for staged assembly info
    if (("groups" in vars) || ("staging" in vars) || ("stages" in vars)){
        // TODO
    }
    else {
        system.clear_staging();
        $("#stagedAssemblyToggle").prop("checked", false);
    }

    window.orbit.target = system.centerOfMass;

    return system
}


// function createKlossSystem() {
//     let vars = getUrlVars();
//     let rule;
//     if ("rule" in vars) {
//         rule = parseKlossString(vars["rule"]);
//     } else if ("hexRule" in vars) {
//         rule = polycubeRuleToKloss(parseHexRule(vars["hexRule"]));
//     } else if ("decRule" in vars) {
//         rule = polycubeRuleToKloss(parseDecRule(vars["decRule"]));
//     } else {
//         let defaultRule = "[[[1,-1.25,0,0,-0.5,1,0,0],[1,-1,-1,0,-1,0,0,0],[1,0.5,0,0,1,1,-0,0],[1,0,-0.5,0,-1,0,1,0],[1,0,0.5,0,0,-1,0,1],[1,0,0,-0.5,-0,1,-1,0],[1,0,0,0.5,1,0,0,1]],[[-1,-0.5,0,0,-1,1,0,0]]]";
//         rule = parseKlossString(defaultRule);
//     }
//
//     system = new KlossSystem(rule, scene, 100);
//     orbit.target = system.centerOfMass;
//
//     system.addParticle(new THREE.Vector3(), new THREE.Quaternion(), system.rule[0], 0);
//     system.processMoves();
//     render();
// }
