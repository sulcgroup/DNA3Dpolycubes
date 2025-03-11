import {allRotations, compCols, coordEqual, FACE_ROTATIONS, randstr, rotCoords} from "./utils";
import {Matrix3, Matrix4, Quaternion, Vector3} from "three";
import {parseDecRule, PolycubeCubeType, ruleToDec} from "./rule";
import {PolycubeSystem} from "./polycubeSystem";

function simplify(rule: PolycubeCubeType[]) : PolycubeCubeType[] {
    let colors = new Set([].concat.apply([], rule.map(r=>r.patches.map(f=>{return f.color}))));
    let newRule: PolycubeCubeType[] = [];
    rule.forEach((cubeType, iCube)=>{
        let allZero = true;
        cubeType.patches.forEach((face, iFace)=>{
            let c = face.color
            if (!colors.has(c*-1)) {
                face.color = 0;
            }
            if (face.color == 0) {
                face.alignDir = FACE_ROTATIONS[iFace];
            }
            else {
                allZero = false;
            }
        })

        if (!allZero || iCube == 0) {
            newRule.push(cubeType);
        }
    });

    let colorset = Array.from(new Set([].concat.apply([],
        rule.map(r=>r.patches.map(f=>{
            return Math.abs(f.color)
        }))))).filter(x => x != 0)
    newRule.forEach(cubeType=>{
        cubeType.patches.forEach(face=>{
            let c = face.color;
            if (c != 0) {
                face.color = colorset.indexOf(Math.abs(c)) + 1;
                if (c < 0) {
                    face.color *= -1;
                }
            }
        })
    })
    return newRule;
}

export function getCoords(rules: PolycubeCubeType[],
                          nMaxCubes=100,
                          torsion=true,
                          piece_size = 1 // pls do not
                         ): Vector3[] | SystemBehaviorType {
    let sys = new PolycubeSystem(nMaxCubes,
        undefined,
        true,
        torsion,
        piece_size,
        randstr(12))
    sys.set_rule(rules);
    sys.seed();
    // run system
    while (!sys.isProcessingComplete()){
        sys.step();
    }
    // if system has moves remaining, it's unbounded
    if (sys.numMovePositions()){
        return SystemBehavior.UNBOUNDED;
        // todo: support for crystals
    }
    else {
        return sys.getCubeCoords();
    }
}

export function simplify2(rule: PolycubeCubeType[], onUpdate) {
    // first-iteration simplification function
    rule = simplify(rule);

    let correctCoords: Vector3[] = getCoords(rule) as Vector3[];
    console.assert("length" in correctCoords, "Cannot simplify non-bounded system"); // joke of a language

    // construct quaternions from rotations, for some reason
    const rotations: [Matrix3, Quaternion][] = allRotations().map(
        r=>[r,
            new Quaternion().setFromRotationMatrix(new Matrix4().setFromMatrix3(r))
        ]
    );

    let alreadyTried = new Set();
    let triedStr = (a,b)=>`${[a,b].sort()}`;

    let updatedRuleset = true;
    while (updatedRuleset) {
        updatedRuleset = false;
        for (const [iA, cA] of rule.entries()) {
            const coordsA = cA.getPatchCoords();
            for (const [iB, cB] of rule.entries()) {
                if (iA < iB) {
                    const coordsB = cB.getPatchCoords();
                    // If they have the same number of patches
                    if (coordsA.length === coordsB.length && !alreadyTried.has(triedStr(iA,iB))) {
                        for (const [r, q] of rotations) {
                            if (compCols(coordsA, rotCoords(coordsB, r))) {
                                console.log(`Cube type ${iA} is similar to ${iB}`);
                                let colorMap = new Map();
                                //let rotMap = new Map();
                                const rotatedB = cB.rotate(q);
                                for (let i=0; i<6; i++) {
                                    if (rotatedB[i].color !== 0) {
                                        console.assert(cA[i].color !== 0);
                                        colorMap.set(rotatedB[i].color, cA[i].color);
                                        //rotMap.set(rotatedB[i].color, getSignedAngle(cB[i].alignDir, cA[i].alignDir, ruleOrder[i]));
                                    }
                                }
                                let newRule = [];
                                // Clone rule and create new ruleset
                                parseDecRule(ruleToDec(rule)).forEach((c,i)=>{
                                    if (i !== iB) {
                                        c.patches.forEach((p,i)=>{
                                            // Replace all of the old color with the new
                                            const color = p.color;
                                            if(colorMap.has(color)) {
                                                p.color = colorMap.get(color);
                                                p.alignDir = cA[i].alignDir;
                                            }
                                            else if(colorMap.has(-color)) {
                                                p.color = -colorMap.get(-color);
                                                p.alignDir.applyQuaternion(q.invert());
                                            }
                                        });
                                    } else {
                                        c.patches.forEach((p,i)=>{
                                            p.color = 0;
                                            p.alignDir = FACE_ROTATIONS[i];
                                        })
                                    }
                                    newRule.push(c);
                                });

                                const simplifiedNewRule = simplify(newRule); //remove empty slots

                                // Accept new ruleset if we still get the same chape
                                let is_bounded_deterministic = getSysBehavior(simplifiedNewRule).short == SystemBehavior.DETERMINISTIC.short;
                                if (is_bounded_deterministic &&
                                    coordEqual(getCoords(simplifiedNewRule) as Vector3[], correctCoords)
                                ) {
                                    updatedRuleset = true;
                                    rule = newRule;
                                    console.log(`Changing  ${iB} to ${iA} did work: ${ruleToDec(simplifiedNewRule)}`);
                                    alreadyTried.add(triedStr(iA,iB));
                                    if (onUpdate) {
                                        onUpdate(ruleToDec(simplifiedNewRule));
                                    }
                                    break;
                                } else {
                                    console.log(`Changing  ${iB} to ${iA} did not work: ${ruleToDec(simplifiedNewRule)}\nBounded & Deterministic = ${is_bounded_deterministic}\nEqual coords = ${coordEqual(getCoords(simplifiedNewRule) as Vector3[], correctCoords)}`);
                                    alreadyTried.add(triedStr(iA,iB));
                                    //break;
                                }
                            }
                        }
                    }
                }
                if (updatedRuleset) {
                    break;
                }
            }
            if (updatedRuleset) {
                break;
            }
        }
    }
    return simplify(rule);
}

export const SystemBehavior = {
    DETERMINISTIC: {
        symbol: "✓",
        short: "DET",
        description: "The system consistently forms the same finite-size shape."
    },
    FINITE_SHAPES: {
        symbol: ":",
        short: "FIN",
        description: "Similar to NON_DETERMINISTIC_BOUNDED but the number of distinct shapes formed is less than sqrt(num tests)."
    },
    UNBOUNDED: {
        symbol: "∞",
        short: "UNB",
        description: "The system grows forever with no pattern."
    },
    NON_DETERNINISTIC_BOUNDED: {
        symbol: "?",
        short: "NDB",
        description: "The system forms finite-size shapes but not always the same one."
    },
    CRYSTAL: {
        symbol: "#",
        short: "CRY",
        description: "The system forms a repeating pattern."
    },
    SOME_NONDET: {
        symbol: "×",
        short: "NDT",
        description: "The system is nondeterministic, no further information is available"
    }
} as const;

// Type for the system behavior
export type SystemBehaviorType = typeof SystemBehavior[keyof typeof SystemBehavior];


export function getSysBehavior(rule: PolycubeCubeType[],
                                   nTries=15,
                                   nMaxCubes = 1000,
                                   torsion:boolean = true,
                                   piece_size=1, // please don't
                                   stopIfNondet=false) : SystemBehaviorType{
    let coordGroups: Vector3[][] = [];
    let count = nTries;
    while (count--) {
        let sysCoords = getCoords(rule, nMaxCubes, torsion, piece_size);
        if ("symbol" in sysCoords) { // check if returned unbounded, this language is a fuckin joke
            console.assert(sysCoords.short == "UNB") // todo: crystal considerations
            return SystemBehavior.UNBOUNDED;
        }
        let foundGroup = false;
        let i =0;
        while (!foundGroup && i < coordGroups.length){
            let toCheck = coordGroups[i]
            // if there's a better search algorithm here I can't think of it
            if (coordEqual(toCheck, sysCoords)){
                foundGroup = true;
            }
            i++;
        }
        if (!foundGroup){
            coordGroups.push(sysCoords);
            if (stopIfNondet){
                return SystemBehavior.SOME_NONDET;
            }
            else if (coordGroups.length >= Math.sqrt(nTries)){ // can't confirm finite shapes, but can refute
                return SystemBehavior.NON_DETERNINISTIC_BOUNDED;
            }
        }


    }
    if (coordGroups.length == 1){
        return SystemBehavior.DETERMINISTIC;
    }
    else if (coordGroups.length < Math.sqrt(nTries)){
        return SystemBehavior.FINITE_SHAPES;
    }
    else {
        return SystemBehavior.NON_DETERNINISTIC_BOUNDED;
    }
}
