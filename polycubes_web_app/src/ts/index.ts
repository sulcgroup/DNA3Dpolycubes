import * as THREE from "three";
import * as $ from "jquery";
import {createPolycubeSystem} from "./init";
import {addSpecies, update, W_LIMIT, stepBack, pause, stepFwd, resume, updateRuleVis, import_system} from "./edit";
import {toggleFog, render} from "./view";
import {downloadStructure, getRandSeed} from "./libs";
import ClickEvent = JQuery.ClickEvent;
import {cubeTypeFromJSON, CubeTypeJson} from "./rule";
import {drawCutterWindow} from "./tangible_export";

// https://stackoverflow.com/questions/10730362/get-cookie-by-name
function getCookie(name) {
    const value = `; ${document.cookie}`;
    const parts = value.split(`; ${name}=`);
    if (parts.length === 2) return parts.pop().split(';').shift();
}

// add listeners to DOM objects
$("#addSpeciesBtn").on("click", function() { addSpecies();});
$("#reseedRandomBtn").on("click", function (){
    window.system.reset_rng(getRandSeed());
    // window.system.reset_rng("69");
    window.system.regenerate();
    resume(1, $("#play"));
});
$("#toggle-grid").on("click", function (){
    window.system.toggleGrid();
});
$("#stagedAssemblyToggle").on("click", function (){
});

$("#reset").on("click", ()=>{update(); resume(500);});

$("#back").on("click", stepBack);
$("#pause").on("click", pause);
$("#stepfwd").on("click", stepFwd);
$("#play").on("click", function (){
    resume(500, $("#play"));
});
$("#ff").on("click", function (){
    resume(100, $("#ff"));
});
$("#fff").on("click", function (){
    resume(1, $("#fff"));
});

$("#stepto").on("change", function (){
    window.system.setMaxCubes(<number>$("#stepto").val())
});

$("#download-lcut").on("click", ()=>{
    if ($("#lcut-window").length == 0) {
        drawCutterWindow(window.system.getCubes());
    }
    else {
        $("#lcut-window").remove();
    }
})

$("#toggle-fog").on("click", ()=>{toggleFog();});

$("#torsion").on("click", function (){
    window.system.setTorsion($("#torsion").prop("checked"));
    update();
});

$("#toggle-rule").on("click", function (){
    $('#display-rule-panel').toggle();
    updateRuleVis();
});

$("#mismatches").on("click", function (e) {
    window.system.setMisMatchesAllowed($("#mismatches").prop("checked"));
    $("#mismatch_count_field").toggle();
    update();
});

$("#envToggle").on("click", function (e: ClickEvent){
    let source = $(e.target);
    source.html(source.html() == "Show Environment" ? "Hide Environment" : "Show Environment");
    $("#env-browser").toggle();
});

$("#toggle-lineview").on("click", function(){
    window.system.toggleLineView();
});

$("#pieceSize").on("change", function (){
    window.system.setMaxPieceSize(<number>$("#pieceSize").val());
    update();
});

$("#uploadsys").on("change", function (e){
    var reader = new FileReader();
    reader.onload = function(e){
        if (typeof this.result === "string") {
            import_system(JSON.parse(this.result));
        }
        else {
            console.log("Bad system file I guess")
        }
    }
    reader.readAsText((e.target as HTMLInputElement).files[0])
});

$("#uploadstructurebtn").on("click", function (){
    $('#uploadsys').trigger('click');
});

$("#uploadrules").on("change", function (e){
    var reader = new FileReader();
    reader.onload = function(e){
        if (typeof this.result === "string") {
            let rule_data = JSON.parse(this.result);
            if (Array.isArray(rule_data)) {
                rule_data = rule_data.map(
                    (r, i) => {
                        if (!("type_id" in <CubeTypeJson>r)) {
                            (<CubeTypeJson>r).type_id = i + window.system.cube_types.length;
                        }
                        return cubeTypeFromJSON(<CubeTypeJson>r, r["typeName"])
                    });
            }
            else {
                rule_data = Object.entries(rule_data).map(
                    ([n, r], i) => {
                        if (!("type_id" in <CubeTypeJson>r)) {
                            (<CubeTypeJson>r).type_id = i + window.system.cube_types.length;
                        }
                        return cubeTypeFromJSON(<CubeTypeJson>r, r["typeName"])
                    });
            }
            let rule = Object.entries(rule_data).map(([n, r]) => cubeTypeFromJSON(<CubeTypeJson>r, r["typeName"]));
            // add all species
            rule.forEach(addSpecies);
            update();
        }
        else {
            console.log("Bad rule file I guess")
        }

    }
    reader.readAsText((e.target as HTMLInputElement).files[0])
});


$("#uploadrulebtn").on("click", function () {
    $("#uploadrules").trigger('click');
    // $("#uploadrulebtn").val();
});

$("#downloadstructure").on("click", downloadStructure);


function saveCanvasImage(){
    window.canvas.toBlob(function(blob){
        var a = document.createElement('a');
        var url = URL.createObjectURL(blob);
        a.href = url;
        a.download = `${window.system.getRuleStr()}.png`;
        a.click();
    }, 'image/png', 1.0);
}

document.addEventListener("keydown", event => {
    if (event.key == 's' && event.ctrlKey) {
        event.preventDefault();
        window.system.export()
    } else if (event.key == 'p') {
        saveCanvasImage();
    }
});
// window.addEventListener('movesProcessed', function(e) {
//     if (!window.transform || !window.transform.object) {
//         fitCamera();
//     }
// }, false);

$(function(){
    // read cookie and set oxdna-dir
    // not currently working, not a priority
    // $('[name="oxdna-dir"]').val(getCookie("oxDNA_home"));
    // $('[name="asu-id"]').val(getCookie("asu_id"));
    $("#rule-display-canvas").attr("width", `${W_LIMIT}px`);
    $("#env-browser").hide();
    $("#patchy-export-panel").hide();
    $("#display-rule-panel").hide();
    $("#graphui").hide();

    window.system = createPolycubeSystem();

    for (let i = 0; i < window.system.numCubeTypes(); i++){
        addSpecies(window.system.cube_types[i]);
    }

    update();
    resume(1);
});
