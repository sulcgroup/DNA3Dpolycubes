import "libs/minisat.js";
import { find_solution } from "./polycubeSolver";
onmessage = function (e) {
    let [topology, empty, nCubeTypes, nColors, nDim, torsionalPatches] = e.data;
    postMessage(find_solution(topology, e, nCubeTypes, nColors, nDim, torsionalPatches));
};
//# sourceMappingURL=solveWorker.js.map