import * as $ from "jquery";
import * as SVG from "@svgdotjs/svg.js";
import * as THREE from "three";
import editIcon from "../ui/edit.svg";
import deleteIcon from "../ui/delete.svg";
import faceIcon from "../ui/face.svg";
import duplicateIcon from "../ui/duplicate.svg";
import {
	cubeTypeColor,
	cubeTypeFromJSON,
	CubeTypeJson,
	DynamicPolycubeCubeType,
	FACE_NAMES,
	parseDecRule,
	patchColor,
	PolycubeCubeType,
	ruleToDec,
} from "./rule"

import {FACE_ROTATIONS, RULE_ORDER, saveString, selectColor, signed_axis_angle} from "./utils";

import {render, VisualPolycubeSystem} from "./view";
import {PolycubePatch} from "./patch";
import {VisualPolycubeCubeInst} from "./cube_instance";

fetch(faceIcon)
	.then(response => response.text())
	.then(svgContent => {
		// Use `svgContent` to create or modify your SVG graphics
	})
	.catch(error => console.error('Error loading SVG:', error));

// import {render, scene} from "./view";

interface PatchOrientationElement extends HTMLImageElement {
	value: number;
}


export function addSpecies(new_species?: PolycubeCubeType | undefined) {
	if (new_species == undefined) {
		// new_species = [];
		let ctName = window.system.nameNewCubeType();
		new_species =  new DynamicPolycubeCubeType(
			window.system.numCubeTypes(),
			FACE_NAMES.map((x, i) => new PolycubePatch(0, i, FACE_ROTATIONS[i])),
			ctName
		);
	}
	if (window.system.cube_types.indexOf(new_species) == -1)
	{
		window.system.addCubeType(new_species);
	}
	let ruleset = document.getElementById("ruleset");
	let ruleField = document.createElement("fieldset");
	let e_autoupdate = document.getElementById("autoUpdate") as HTMLInputElement;
	$(ruleField).addClass("rule-display");
	$(ruleField).append($(`<legend name="${new_species.typeName}">${new_species.typeName}</legend>`));
	ruleField.style.borderColor = cubeTypeColor(new_species.type_id)
	for (let i = 0; i < RULE_ORDER.length; i++) {
		let face = document.createElement("span");
		//face.faceIdx = i;
		let color_field = document.createElement("input") as HTMLInputElement;
		color_field.type = "number";

		color_field.value = String(new_species.patch_color(i));
		let color = parseInt(color_field.value)
		if (color != 0) {
			face.style.backgroundColor = patchColor(parseInt(color_field.value))
		}
		color_field.addEventListener("change", updateCubeColor.bind(
			event, color_field, new_species, i)
		);
		if (window.system.isPolycubeSystem()) {
			face.title = FACE_NAMES[i];
		}
		if (window.system.isPolycubeSystem()) { //under most operational circumstances (as of this writing)
			let rotation = document.createElement("img") as PatchOrientationElement;
			rotation.value = FACE_ROTATIONS[i].angleTo(new_species.patches[i].alignDir) * (2 / Math.PI);

			rotation.src = faceIcon;
			rotation.className = `face-img rot${rotation.value}`;
			rotation.addEventListener("click", updatePatchRot.bind(
				event, rotation, new_species, i)
			);
			face.appendChild(rotation);
			face.appendChild(color_field);
		}
		ruleField.appendChild(face);
	}

	ruleField.appendChild(document.createElement("br"));
	let control_button_size = 30;
	
	//add button to delete rule
	let remove = document.createElement("img");
	remove.src = deleteIcon;
	remove.classList.add("control_button");
	remove.height = control_button_size;
	remove.alt = "Remove";
	remove.addEventListener("click", removeCubeType.bind(
		event, new_species, ruleField)
	);
	ruleField.appendChild(remove);
	
	//add button to duplicate rule
	let duplicate = document.createElement("img");
	duplicate.src = duplicateIcon;
	duplicate.height = control_button_size;
	duplicate.alt = "Duplicate";
	duplicate.classList.add("control_button");
	duplicate.addEventListener("click", duplicateCubeType.bind(
		event, new_species, ruleField)
	);
	ruleField.appendChild(duplicate);

	//add button to edit rule
	let edit = document.createElement("img");
	edit.src = editIcon; // Use the imported SVG file
	edit.height = control_button_size;
	edit.alt = "Edit";
	edit.classList.add("control_button");
	edit.addEventListener("click", editCubeType.bind(
		event, new_species, ruleField)
	);
	ruleField.appendChild(edit);

	let togglePresent = $(`<svg xmlns="http://www.w3.org/2000/svg"
     width="${control_button_size}" height="${control_button_size}"
     viewBox="0 0 256 256">
  <circle id="big" cx="128" cy="128" r="110" stroke="black" stroke-width="12" fill="gray"/>
</svg>`)
	togglePresent.addClass("control_button");
	togglePresent.attr("alt", "Toggle Present");
	togglePresent.click(function() {
		let circ = $(this).children("circle");
		let active = circ.attr("fill") == "gray";
		window.system.toggleCubeType(window.system.cubeTypeIndex(new_species));
		if (active) {
			circ.attr("fill", "white");
		} 
		else {
			circ.attr("fill", "gray");
		}
		update();
	});
	$(ruleField).append(togglePresent);
	ruleset.appendChild(ruleField);
}


function updateCubeColor(e, new_species: PolycubeCubeType, faceIdx: number) {
	// let ctidx = window.system.cube_types.indexOf(new_species);
	let c;
	let value = parseInt(e.value);
	if (value != 0) {
		c = patchColor(value);
	}
	else {
		c = "White";
	}
	e.parentElement.style.backgroundColor = c;
	new_species.patches[faceIdx].color = value;
	if ((document.getElementById("autoUpdate") as HTMLInputElement).checked) {
		update();
	}
}

function updatePatchRot(e: any, new_species: PolycubeCubeType, faceIdx: number) {
	$(e).removeClass(`rot${e.value}`);
	// let ruleIdx = window.system.cube_types.indexOf(new_species);
	e.value = (parseInt(e.value) + 1) % 4;
	$(e).addClass("rot" + e.value);
	let r = FACE_ROTATIONS[faceIdx].clone();
	r.applyAxisAngle( RULE_ORDER[faceIdx], e.value * Math.PI / 2);
	new_species.patches[faceIdx].alignDir = r.round();
	if ((document.getElementById("autoUpdate") as HTMLInputElement).checked) {
		update();
	}
}

function removeCubeType(rule: PolycubeCubeType, ruleField) {
	let ncts = window.system.numCubeTypes();
	ruleField.parentNode.removeChild(ruleField);
	window.system.removeCubeType(rule);
	console.assert(ncts - 1 == window.system.numCubeTypes(), "Appear to have removed more than one cube types?")
	$('#ruleset').children().each(function(idx){
		$(this).css('borderColor', cubeTypeColor(idx));
	});
	if ((document.getElementById("autoUpdate") as HTMLInputElement).checked) {
		update();
		resume();
	}
}

function duplicateCubeType(rule: PolycubeCubeType) {
	// let ruleIdx = window.system.cube_types.indexOf(rule);
	let copy = rule.clone(window.system.numCubeTypes());
	window.system.addCubeType(copy);
	addSpecies(copy);
	let cubeMaterial = new THREE.MeshLambertMaterial({
		color: cubeTypeColor(window.system.numCubeTypes()-1)
	});
	window.system.particleMaterials.set(copy.type_id, cubeMaterial);
	if ((document.getElementById("autoUpdate") as HTMLInputElement).checked) {
		update();
	}
}

/**
creates a new subwindow to edit a rule in detail 
 */
function editCubeType(cube_type, ruleField) {
	if ($("#editbox").length > 0) {
		return;
	}
	let editbox = document.createElement("fieldset");
	let staged_name_change = cube_type.typeName;
	// make array for staged rotation changes
	let staged_rotations = [];
	// make array for staged color changes
	let staged_colors = [];
	// make array for staged allostery changes
	let staged_allostery_alterations;
	if ("effects" in cube_type) {
		staged_allostery_alterations = {
			"patches": [],
			"effects": []
		};
	}
	else {
		staged_allostery_alterations = [];
	}
	editbox.id = "editbox";
	let legend = document.createElement("legend");
	legend.innerHTML = `Editing ${cube_type.typeName}`;
	editbox.appendChild(legend);
	let name_field = document.createElement("input")
	name_field.value = cube_type.typeName;
	// event for cube type name change
	name_field.addEventListener("change", e=>{
		staged_name_change = (e.target as HTMLInputElement).value;
		legend.innerHTML = `Editing ${staged_name_change}`;
	});
	editbox.appendChild(name_field);
	let rule_faces = document.createElement("table");
	$(rule_faces).append(
		$("<tr><td>Face</td><td>Rot.</td><td>Color</td><td>S.V</td><td>A.V</td></tr>")
	)
	for (let i = 0; i < RULE_ORDER.length; i++) {
		let face = document.createElement("tr");
		let face_label = document.createElement("td");
		face_label.innerHTML = `${i}:${FACE_NAMES[i]}:`;
		face.appendChild(face_label)
		//under most operational circumstances (as of this writing)
		if (window.system.isPolycubeSystem()) {
			let rotation = document.createElement("img") as PatchOrientationElement;
			rotation.value = FACE_ROTATIONS[i].angleTo(cube_type.patches[i].alignDir) * (2 / Math.PI);

			rotation.src = faceIcon;
			rotation.className = `face-img rot${rotation.value}`;
			rotation.addEventListener("click", e=>{
					let rot = (parseInt((e.target as HTMLInputElement).value) + 1) % 4;
					$(e.target).removeClass(`rot${$(e.target).val()}`);
					$(e.target).addClass(`rot${rot}`);
					$(e.target).val(rot)
					staged_rotations[i] = rot;
				}
			);
			face.appendChild(rotation);
		} //TODO: if kloss ever gets updated, add an alternative field here
		
		let color_input = document.createElement("input");
		let color_cell = document.createElement("td");
		color_input.style.display = "table-cell";
		color_input.type = "number";
		color_input.value = cube_type.patches[i].color;
		let color = parseInt(color_input.value)
		if (color != 0) {
			face.style.backgroundColor = patchColor(parseInt(color_input.value))
		}
		// function on color change
		color_input.addEventListener("change", e=>{
			let value = parseInt((e.target as HTMLInputElement).value);
			let c;
			if (value != 0) {
				if (!(Math.abs(value) in window.system.colorMaterials)) {
					window.system.colorMaterials[value] = new THREE.MeshLambertMaterial({
						color: patchColor(value)
					});
				}

				c = patchColor(Math.abs(value) - 1);
			}
			else {
				c = "White";
			}
			(e.target as HTMLInputElement).parentElement.style.backgroundColor = c;
			staged_colors[i] = value;
		});
		if (window.system.isPolycubeSystem()) {
			face.title = FACE_NAMES[i];
		}
		color_cell.appendChild(color_input);
		face.appendChild(color_cell);

		// if cube type is static formulation
		if ("conditionals" in cube_type) {
			let logic = document.createElement("input");
			let logic_cell = document.createElement("td");
			logic.style.display = "table-cell";
			logic.type = "text";
			logic.value = cube_type.conditionals[i];
			logic.addEventListener("change", e => {
				staged_allostery_alterations[i] = (e.target as HTMLInputElement).value;
			});
			logic_cell.appendChild(logic);
			face.append(logic_cell)
		}
		// dynamic formualtion
		else {
			let state_var = $("<td>").append($("<input>").val(
				cube_type.patches[i].state_var
			).attr({
				"type": "number",
				"step": "1"
			}).on("input", function (){
				if (!(i in staged_allostery_alterations.patches)){
					staged_allostery_alterations.patches[i] = {};
				}
				staged_allostery_alterations.patches[i].state_var = parseInt((this as HTMLInputElement).value);
			}));
			let activation_var = $("<td>").append($("<input>").val(
				cube_type.patches[i].activation_var
			).attr({
				"type": "number",
				"step": "1"
			}).on("input", function (){
				if (!(i in staged_allostery_alterations.patches)){
					staged_allostery_alterations.patches[i] = {};
				}
				staged_allostery_alterations.patches[i].activation_var = parseInt((this as HTMLInputElement).value);
			}));
			face.appendChild(state_var.get(0))
			face.appendChild(activation_var.get(0))
		}
		rule_faces.appendChild(face);
	}
	editbox.appendChild(rule_faces);

	function add_effect_ui(e, i) {
		return $("<tr>").append(
			$("<td>").append($("<input>")
				.val(e.sources.join(","))
				.on("input", function () {
					// validate input
					let pattern = /^-?[0-9]+(,-?[0-9]+)*$/;
					if (pattern.test((this as HTMLInputElement).value)) {
						if (!(i in staged_allostery_alterations.effects)){
							staged_allostery_alterations.effects[i] = {}
						}
						let sources = (this as HTMLInputElement).value.split(",").map(x => parseInt(x));
						staged_allostery_alterations.effects[i].sources = sources;
					}
				})
			),
			$("<td>").append($("<input>")
				.val(e.target)
				.attr({
					"type": "number",
					"step": "1",
					"min": 0
				})
				.on("input", function () {
					if (!(i in staged_allostery_alterations.effects)){
						staged_allostery_alterations.effects[i] = {}
					}
					staged_allostery_alterations.effects[i].target = (this as HTMLInputElement).value;
				})
			)
		);
	}

	// if the cube is dynamic formulation, add a box for effects
	if ("effects" in cube_type) {
		let n_ct_fx = cube_type.effects.length;
		$(editbox).append($("<div>")
			.css({
				"overflow": "hidden",
				"overflow-y": "scroll",
				"height": "150px"
			})
			.append($("<table id='ct-fx-list'>").append(
					$("<tr><td>Source</td><td>Target</td></tr>"),
					...cube_type.effects.map(add_effect_ui)
				),
				$("<button>+</button>").on("click", function (){
					$("#ct-fx-list").append(add_effect_ui({
						"sources": [],
						"target": 0
					}, n_ct_fx));
					n_ct_fx++;
				})
			)
		);
	}

	let apply_button = document.createElement("button")
	apply_button.innerHTML = "Apply";
	editbox.appendChild(apply_button);

	// function to exec when the "apply" button is clicked
	apply_button.addEventListener("click", e=>{
		$(`legend[name=${cube_type.name}]`).text(staged_name_change);
		cube_type.typeName = staged_name_change;
		let flag_made_change = false;
		// loop rotations
		staged_rotations.forEach((rot_key, rot)=>{
			// grab old rotation
			let rval_old = $(ruleField).children(`[title=${FACE_NAMES[rot_key]}]`).children("img").val();
			// if new rotation is different than old rotation
			if (rot != rval_old){
				// apply rotation
				$(ruleField).children(`[title=${FACE_NAMES[rot_key]}]`).children("img").removeClass(`rot${rval_old}`);
				let r = FACE_ROTATIONS[rot_key].clone();
				r.applyAxisAngle(RULE_ORDER[rot_key], staged_rotations[rot_key] * Math.PI / 2);
				cube_type.patches[rot_key].alignDir = r.round();
				$(ruleField).children(`[title=${FACE_NAMES[rot_key]}]`).children("img").val(staged_rotations[rot_key]);
				$(ruleField).children(`[title=${FACE_NAMES[rot_key]}]`).children("img").addClass(`rot${staged_rotations[rot_key]}`)
				flag_made_change = true; // set flag
			}
		});
		// loop colors
		staged_colors.forEach((color_int, color_key)=>{
			// if new color is different from old color
			if (cube_type.patches[color_key].color != color_int) {
				// set color
				cube_type.patches[color_key].color = staged_colors[color_key];
				$(ruleField).children(`[title=${FACE_NAMES[color_key]}]`).children("input").val(staged_colors[color_key]);
				let color = staged_colors[color_key] != 0 ? patchColor(Math.abs(staged_colors[color_key]) - 1) : "white";
				$(ruleField).children(`[title=${FACE_NAMES[color_key]}]`).css("background-color", color);
				flag_made_change = true; // set flag
			}
		});
		// update static fomrulatiuon cube
		if ("conditionals" in cube_type){
			for (const conditional_key in staged_allostery_alterations){
				cube_type.conditionals[conditional_key] = staged_allostery_alterations[conditional_key].replace("/\s/g", '');
			}
		}
		// update dynamic formulation cube
		else {
			let flag_changed_allo = false
			staged_allostery_alterations.patches.forEach((allo_updates, idx)=>{
				if ("state_var" in allo_updates){
					// set state var
					cube_type.patches[idx].state_var = allo_updates.state_var;
					flag_changed_allo = true; // set flag
				}
				if ("activation_var" in allo_updates){
					// set activation var
					cube_type.patches[idx].activation_var = allo_updates.activation_var;
					flag_changed_allo = true; // set flag
				}
			});
			let offset=0; // correct effects indexes when we remove them
			staged_allostery_alterations.effects.forEach((e, i)=>{
				if (i in cube_type.effects) {
					if ("target" in e && cube_type.effects[i-offset].target != e.target) {
						if (e.target > 0) {
							cube_type.effects[i-offset].target = e.target;
						}
						// setting an effect to target state 0 deletes it
						else {
							cube_type.effects.splice(i-offset, -1);
							offset++;
						}
						flag_changed_allo = true; // set flag
					}

					if ("sources" in e && cube_type.effects[i-offset].sources != e.sources) {
						cube_type.effects[i-offset].sources = e.sources;
						flag_made_change = true;
					}
				}
				else {
					cube_type.effects.push(e);
				}
			});
			if (flag_changed_allo){
				flag_made_change = true;
				cube_type.update_state_size();
			}
		}
		// skip update unless change was made
		if (flag_made_change) {
			window.system.refresh_piecewise_staging();
			update();
		}
		document.body.removeChild(editbox);
	});
	let cancel_button = document.createElement("button");
	cancel_button.innerHTML = "Cancel";
	editbox.appendChild(cancel_button);
	cancel_button.addEventListener("click", e=>{
		document.body.removeChild(editbox);
	});
	document.body.appendChild(editbox);
	$(editbox).css("top", Math.min($(ruleField).offset().top, window.innerHeight - $(editbox).height() - 25))
}

export function clearRules() {
	document.getElementById("ruleset").innerText = "";
	window.system.cube_types.forEach(addSpecies);
	if ((document.getElementById("autoUpdate") as HTMLInputElement).checked) {
		window.system.regenerate();
		render();
	}
}

export function rgbToHex(rgb) {
	return "#" + ((1 << 24) + (rgb.r * 255 << 16) + (rgb.g * 255 << 8) + rgb.b * 255).toString(16).slice(1);
}

function toggleRuleSet() {
	let ruleDiv = document.getElementById("ruleset");
	let hideToggle = document.getElementById("hideToggle");
	if (ruleDiv.style.height == "0px") {
		ruleDiv.style.height = "100%"
		hideToggle.innerHTML = "Hide";
	} else {
		ruleDiv.style.height = "0px"
		hideToggle.innerHTML = "Show";
	}
}


export function exportRule() : string{
	return ruleToDec(window.system.listCubeTypes());
}

// map of 2d grid coordinates for faces in the cutout
// index is face enumeration, numbers are x and y translation of face
export const FACE_SVG_COORDS = [
	[2, 2], //left
	[0, 2], //right
	[1, 0], //bottom
	[1, 2], //top
	[1, 3], //back
	[1, 1] //front
]

export const FACE_ROT_OFFSETS = [
	1, //left
	1, //right
	0, //bottom
	0, //top
	1, //back
	-1 //front
]

export let W_LIMIT = 800;

let big_cum_x;
let big_cum_y;

let row_h = 325;

export function updateRuleVis() {
	let rule_display = $("#rule-display-canvas");
	rule_display.empty();
	var svg: SVG.Container = SVG.SVG("#rule-display-canvas") as SVG.Container;
	big_cum_x = 0;
	big_cum_y = 0;
	let w_max = W_LIMIT;
	window.system.cube_types.forEach((ct)=>{
		drawCubeType(svg, ct.type_id);
		if (big_cum_x > W_LIMIT){
			w_max = Math.max(w_max, big_cum_x)
			big_cum_x = 0;
			big_cum_y += row_h;
		}
	});
	rule_display.attr("width", `${w_max}`);
	rule_display.attr("height", `${big_cum_y + row_h}`);
}


const FACE_SQUARE_W = 64;
const FACE_FIG_PAD = 12;

function drawCubeType(svg: SVG.Container, ct_idx: number){
	let cubeType: PolycubeCubeType = window.system.cube_types[ct_idx];
	let f;
	let x = big_cum_x;
	let y = big_cum_y;
	let group = svg.group();
	let arrow_end = svg.marker(10, 10, function (add) {
		add.path("M 0 0 L 10 5 L 0 10 z").fill("black");
	});
	// move cube type group
	group.transform({translateX: x, translateY: y});
	for (f = 0; f < RULE_ORDER.length; f++){
		[x, y] = FACE_SVG_COORDS[f];
		const cube_type_color = selectColor(ct_idx, 57, 48)
		group.rect(FACE_SQUARE_W, FACE_SQUARE_W).move(x * FACE_SQUARE_W, y * FACE_SQUARE_W).fill(cube_type_color);
		if (cubeType.patches[f].color != 0){
			const patch_color = patchColor(cubeType.patches[f].color);
			const patch_angle = -signed_axis_angle(FACE_ROTATIONS[f], cubeType.patches[f].alignDir, RULE_ORDER[f])
			const face_rot_idx = patch_angle * (2 / Math.PI);
			const rotation = (face_rot_idx + FACE_ROT_OFFSETS[f]) * 90
			group.path("m 30,85 0,-45 10,0 0,-15 22,0 0,15 10,0 0, 45 z")
				.fill(patch_color)
				.stroke({'color': cubeType.patches[f].color > 0 ? "black" : "white", 'width': 3.5})
				.transform({
					scale: .5
				})
				.transform({
					translateX: x * FACE_SQUARE_W - 20,
					translateY: y * FACE_SQUARE_W - 20,
				}, true)
				.transform({rotate: rotation}, true);
			
		}
	}
	// cumulative x coordinate
	let x_cum = FACE_SQUARE_W * 3 + 2 * FACE_FIG_PAD
	function activates(x1: number, y1: number,
					   x2: number, y2: number, deactivates=false, arc_r=1.25){
		let radius = Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2)) * arc_r;
		let color = deactivates ? "red" : "green";
		group.path(`M ${x1} ${y1} A ${radius} ${radius} 0 0 1 ${x2} ${y2}`)
			.stroke(color)
			.fill("none")
			.marker("end", arrow_end);
	}

	function coordspace(n: number) : number{
		return (n + 0.5) * FACE_SQUARE_W;
	}

	// only to be used in dynamic model cubes
	function state_box_intr_coords(s: number, side = "left") {
		console.assert(s > -cubeType.state_size() - 1);
		console.assert(s <  + cubeType.state_size() + 1);
		let x = state_box_x + (state_box_w * (side == "left" ? 0 : 1));
		x += 6 * (side == "left" ? -1 : 1);
		let y = (s + 0.5) * state_box_w + row_center;
		// y += states_interact_counter[s + cubeType.state_size()] * 10 * (states_interact_counter % 2 == 1 ? 1 : -1);
		// states_interact_counter[s + cubeType.state_size()]++;
		return [x, y]
	}

	function sets_state(x1: number, y1: number, state_idx: number, x2=undefined){
		let x3, y3;
		let approach_dir = x1 < state_box_x ? "left" : "right";
		[x3, y3] = state_box_intr_coords(state_idx, approach_dir);
		if (!x2) {
			x2 = (x1 + x3) / 2;
		}
		let y2 = (y1 + y3) / 2;
		group.path(`M ${x1} ${y1} 
						S ${x2} ${y1} ${x2} ${y2}
						S ${x2} ${y3} ${x3} ${y3}`)
			.stroke("green")
			.fill("none")
			.marker("end", arrow_end);
	}

	let x1, x2, y1, y2;

	// if static formulation
	if ("effects" in cubeType) {
		if ((cubeType as DynamicPolycubeCubeType).effects.length > 0) {
			var state_box_w = row_h / (2 * cubeType.state_size() + 3);
			var row_center = row_h / 2 - 1.5*state_box_w;
			x_cum += + 2 * FACE_FIG_PAD;
			var state_box_x = x_cum;
			x_cum += state_box_w;
			// loop state variables
			for (let s = -cubeType.state_size(); s < cubeType.state_size() + 1; s++) {
				// draw state var
				group.rect(state_box_w, state_box_w)
					.move(state_box_x, row_center + s * state_box_w)
					.stroke("black")
					.fill("white");
				group.text(`${s}`)
					.move(state_box_x + state_box_w / 2, row_center + (s + 0.5) * state_box_w)
					.font({size: 10, anchor: 'middle'});
			}
			// let states_interact_counter = Array(cubeType.state_size() * 2 + 1).fill(0);

			for (let f = 0; f < RULE_ORDER.length; f++) {
				if (cubeType.patches[f].color != 0) {
					if (cubeType.patches[f].state_var != 0) { // don't clutter the diagram with patches that have no state var
						[x1, y1] = FACE_SVG_COORDS[f];
						sets_state(coordspace(x1), coordspace(y1), cubeType.patches[f].state_var);
					}
					if (cubeType.patches[f].activation_var != 0) { // don't clutter the diagram with patches that have no state var
						[x1, y1] = state_box_intr_coords(cubeType.patches[f].activation_var);
						[x2, y2] = FACE_SVG_COORDS[f];
						activates(x1, y1, coordspace(x2), coordspace(y2), cubeType.patches[f].activation_var < 0);
					}
				}
			}
			let dest_x, dest_y;
			x_cum += 2 * FACE_FIG_PAD;
			(cubeType as DynamicPolycubeCubeType).effects.forEach((e, j) => {
				let origins_coords_list = e.sources.map(s=>{return state_box_intr_coords(s, "right");});
				[dest_x, dest_y] = state_box_intr_coords(e.target, "right");

				let e_draw_size, e_draw_y;
				if (e.sources.length == 1){
					e_draw_size = state_box_w;
					e_draw_y = origins_coords_list[0][1] - state_box_w/2;
				}
				else {
					e_draw_size = state_box_w * 0.75 * origins_coords_list.length;
					e_draw_y = origins_coords_list.reduce((a,p)=>{
						return a + (p[1] / origins_coords_list.length);
					}, 0);
				}
				group.path(`M ${x_cum} ${e_draw_y} 
							h ${e_draw_size * 0.25}
							a ${e_draw_size / 2} ${e_draw_size / 2} 0 0 1 0 ${e_draw_size}
							h ${-e_draw_size * 0.25}
							z`)
					.stroke("black")
					.fill("white");
				origins_coords_list.forEach((xy,k)=>{
					[x1, y1] = xy;
					let x2 = x_cum - 5;
					let y2 = e_draw_y + (k+1) * e_draw_size / (origins_coords_list.length + 1);
					activates(x1, y1, x2, y2, e.sources[k] < 0, 4);
				});
				x_cum += e_draw_size * 0.75;
				// activates(x_cum, e_draw_y + e_draw_size/2, dest_x, dest_y, e.target < 0,0.75);</
				y1 = e_draw_y + e_draw_size / 2;
				let x_mid = x_cum + Math.sqrt(Math.abs(e_draw_y + e_draw_size / 2 - dest_y));
				sets_state(x_cum, y1, e.target, x_mid);
				x_cum += 22;
				if (e.sources.length == 1){
					x_cum -= 10;
				}
			});
		}
	}

	group.text(cubeType.typeName)
		.move(96, 275)
		.font({
			size: 20,
			anchor: 'middle'
		})
	big_cum_x += x_cum;
}


function updateSystemGraph(){
	// wip
}


export function import_system(jssys: {"cube_types": any, "cubes": any[]}) {
	// import species from json tag "cube_types"
	let rule = Object.entries(jssys.cube_types).map(
		([n, r], i) => {
			if (!("type_id" in <CubeTypeJson>r)){
				(<CubeTypeJson>r).type_id = i;
			}
			return cubeTypeFromJSON(<CubeTypeJson>r, r["typeName"])
		});
	window.system = new VisualPolycubeSystem(
		window.scene,
		jssys.cubes.length,
		100, //whatever
		false, //
		true,
		window.system.envbrowser); //new window.system
	// add all species
	rule.forEach(addSpecies);
	jssys.cubes.forEach(c => {
			let ctidx = c.type;
			let cube_position = new THREE.Vector3(c.position.x, c.position.y, c.position.z);
			let cube_rotation = new THREE.Quaternion(c.rotation.x, c.rotation.y, c.rotation.z, c.rotation.w);
			let cinst: VisualPolycubeCubeInst = window.system.instantiateCube(rule[ctidx], c["name"], cube_position);
			cinst.rotate(cube_rotation);
			cinst.set_state(c.state);
			window.system.addParticle(cinst);
		}
	);

	render();
}

function uploadSystem(source: any) {
	var reader = new FileReader();
	reader.onload = function(e){
		if (typeof this.result === "string") {
			import_system(JSON.parse(this.result));
		}
		else {
			console.log("Bad system file I guess")
		}
	}
	reader.readAsText(source.files[0])
}

function download_topology() {
    let [topology, _] = window.system.getCurrentTopology();
    let out = {
        nDim: 3,
        bindings: topology,
        torsion: window.system.hasTorsion()
    };
    saveString(
        JSON.stringify(out),
        'topology.json'
    );
}


let SELECTED_BTN_BG = "grey";
let UNSELECTED_BTN_BG = "white";

/** 
 */
 export function update(){
	window.system.regenerate();
	if ($("#display-rule-panel").is(":visible")){
		updateRuleVis();
	}
	render();
}

/**
the simulation can't actually "step back", so this function simply runs the simulation one step before current
*/
export function stepBack() {
	$("#stepto").val(window.system.numSteps()-1);
    window.system.regenerate();
    resume(0);
}

/**
stops automatically adding blocks
 */
export function pause() {
	$(".stepctrlbtn").css("background-color", UNSELECTED_BTN_BG);
	window.clearTimeout(stepTimeoutID)
	$("#pause").css("background-color", SELECTED_BTN_BG);
}

/**
adds one block
 */
export function stepFwd() {
	doStep();
}

/**

 */
let stepTimeoutID; // Timeout for steps

export function resume(step_length = 0, src: undefined | JQuery = undefined) {
	$(".stepctrlbtn").css("background-color", UNSELECTED_BTN_BG);
	window.clearTimeout(stepTimeoutID); // Clear any existing timeout to prevent duplicates

	// Function to execute doStep and schedule the next invocation
	function scheduleNextStep() {
		doStep();
		// Schedule the next execution of doStep after the current one completes
		if (!window.system.isProcessingComplete()){
			stepTimeoutID = window.setTimeout(scheduleNextStep, step_length);
		}
	}

	// Start the execution loop
	stepTimeoutID = window.setTimeout(scheduleNextStep, step_length);
	console.log(`Setting timer ${stepTimeoutID} with interval ${step_length}`);

	if (src) {
		$(src).css("background-color", SELECTED_BTN_BG);
	}
}



function doStep() {
	window.system.step();
	$("#stepcount").html(`Step: ${window.system.numSteps()}`);
	let step_to = parseInt(<string>$("#stepto").val())
	if (window.system.numMovePositions() == 0 || ( (step_to <= window.system.numSteps()) && (step_to != 0))) {
		pause();
	}
	render();
}

// (window as any).addSpecies = addSpecies;

function addRule(rule_str, fmt = 0) {
	let rule;
	rule = parseDecRule(rule_str);
	rule.forEach(function (r) {
		addSpecies(r);
	});
	update();
}