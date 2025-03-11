import {parseDecRule, ruleToDec} from "./rule";
import {simplify2} from "./analysis";

onmessage = function(e) {
    let inputRule = parseDecRule(e.data);
    let result = ruleToDec(simplify2(inputRule, rule=>postMessage(rule)));
    console.log("Got result: "+result);
    postMessage(result);
}