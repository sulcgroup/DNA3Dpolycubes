import {polysat} from "./polycubeSolver";

interface SolveWorkerMsgEvent extends MessageEvent {
    data: polysat;
}

onmessage = function(e: SolveWorkerMsgEvent) {
    postMessage((e.data as polysat).run_solve());
}