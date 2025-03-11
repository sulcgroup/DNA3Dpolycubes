// // DEPRECATED!!! use pypatchy ffs!!!
// import * as $ from "jquery";
// import PolycubeSystem from "./polycubeSystem"
//
// let export_groups = [];
// let export_system = undefined;
//
// let input_file_params;
// input_file_params = null;
//
// const inputFileName = 'input';
// const topFileName = 'init.top';
// const confFileName = 'init.conf';
// const particlesFileName = 'particles.txt';
// const patchesFileName = 'patches.txt';
//
// class ExportGroup {
// 	private exportGroupName: any;
// 	private particle_type_levels: any;
// 	private num_assemblies: number;
// 	private num_duplicates: number;
// 	private particle_density: number;
// 	private narrow_types: any[];
// 	private temperatures: any;
// 	constructor(name){
// 		this.exportGroupName = name;
// 		this.num_assemblies = 10;
// 		this.particle_type_levels = getCubeTypeCount(export_system.cube_types);
// 		this.temperatures = {};
// 		$('[name="temperatures"]').val().split(",").forEach(t=>{this.temperatures[parseFloat(t)] = true;});
// 		this.num_duplicates = 1;
// 		this.particle_density = 0.1;
// 		this.narrow_types = []
// 	}
//
// 	visualize(){
// 		let workaround = this; //livid
// 		let vis = $("<div>");
// 		vis.addClass("patchy-export-group");
// 		vis.attr("name", this.exportGroupName);
// //		let numcols = this.rule_levels.length+1;
// //		vis.css("grid-template-columns", [...Array(this.rule_levels.length + 1).keys()].map(i=>`${Math.round(100 / numcols)}%`).join(" "));
//
// 		//name
// 		vis.append($('<label for="export-group-name">Name</label>'));
// 		let name_field = $('<input name="export-group-name"></input>');
// 		name_field.val(this.exportGroupName);
// 		name_field.change(e=>{
// 			this.exportGroupName = $(e.target).val();
// 		});
// 		vis.append(name_field);
//
// 		//num assemblies
// 		vis.append($("<label for='num-assemblies' text-align:right;'>Num Assemblies:</label>"));
// 		let assemblies_field = $("<input name='num-assemblies' type='number'>");
// 		assemblies_field.css("text-align", "left");
// 		assemblies_field.css("border", "1px solid black");
// 		assemblies_field.val(this.num_assemblies);
// 		vis.append(assemblies_field);
// 		assemblies_field.change(e=>{
// 			this.num_assemblies = parseInt($(e.target).val());
// 			this.particle_type_levels.map((r, i) =>{
// 				let rule_level_labels = vis.children(".rule-container").children("label")
// 				rule_level_labels.text(function(i,t)
// 				{
// 					let num = workaround.particle_type_levels[i] * assemblies_field.val();
// 					return	`${export_system.cube_types[i].name}: ${num}`;
// 				});
// 			});
// 		});
//
// 		//num duplicates
// 		vis.append($("<label for='num-duplicates' style='text-align:right;'>Num Duplicates:</label>"));
// 		let duplicates_field = $("<input name='num-duplicates' type='number'>");
// 		duplicates_field.css("text-align", "left");
// 		duplicates_field.css("border", "1px solid black");
// 		duplicates_field.val(this.num_duplicates);
// 		duplicates_field.change(e=>{
// 			this.num_duplicates = parseInt(<string>$(e.target).val());
// 		});
// 		vis.append(duplicates_field);
//
// 		vis.append($("<div style='width:100%;'></div>"))
//
// 		//cube type levels
// 		this.particle_type_levels.forEach((rl, i)=>{
// 			let rule_container = $("<div class='rule-container'>");
// 			let rule_level_input = $('<input type="range" min="0" class="rule-level">');
// 			rule_level_input.attr("name",`${export_system.cube_types[i].typeName}-count`);
// 			rule_level_input.attr("type", "range");
// 			rule_level_input.attr("min", 0);
// 			rule_level_input.attr("max", $('[name="max-rule-particles"]').val());
// 			rule_level_input.val(this.particle_type_levels[i])
// 			let rule_level_label = $(`<label class="group-temperature-label" for="${export_system.cube_types[i].typeName}-count">`);
// 			rule_level_label.text(`${export_system.cube_types[i].typeName}: ${this.particle_type_levels[i] * assemblies_field.val()}`);
// 			rule_container.append(rule_level_input);
// 			rule_container.append($("<br>"));
// 			rule_container.append(rule_level_label);
// 			rule_level_input.change(e=>{
// 				this.particle_type_levels[i] = parseInt(<string>$(e.target).val());
// 				rule_level_label.text(`${export_system.cube_types[i].typeName}: ${this.particle_type_levels[i] * assemblies_field.val()}`);
// 			});
// 			rule_level_input.css("-webkit-appearance", "slider-vertical");
// 			vis.append(rule_container);
// 		});
//
// 		//temperatures
// 		let temp_check_buttons = $('<div class="temperature-grouping" name="group-temps">');
// 		temp_check_buttons.append($(`<label for="group-temps" style="margin-top: auto;">Temperatures</label>`));
// 		vis.append(temp_check_buttons);
// 		setTimeout(function(){this.updateTemperatures()}, 500);
// 		//remove export group
// 		let rm_btn = $('<img src="ui/delete.svg" class="control_button">');
// 		rm_btn.css("width", 30);
// 		rm_btn.css("height", 30);
// 		rm_btn.click(function (e){
// 			let idx = export_groups.findIndex(x=>x==workaround);
// 			export_groups.splice(idx, 1);
// 			$(this).parent().remove();
// 		});
// 		vis.append(rm_btn);
// 		return vis;
// 	}
//
// 	updateTemperatures(){
// 		let new_temperatures = $('[name="temperatures"]').val().split(",").map(parseFloat).filter(Boolean);
// 		new_temperatures.forEach(t=>{
// 			if ($(`[name="${this.exportGroupName}"] > .temperature-grouping > [name="T${t}-span"]`).length == 0){
// 				this.addTemperatureControl(t);
// 			}
// 		});
// 		$(`[name="${this.exportGroupName}"] > .temperature-grouping`).children('input').filter((i, e)=>{ //order for jquery is idx, elmnt
// 			return !($(e).attr("name") in new_temperatures);
// 		}).remove();
// 	}
//
// 	addTemperatureControl(t, on=true){
// 		this.temperatures[t] = on;
// 		let tinput = $("<input type='checkbox'>");
// 		tinput.attr("name", t);
// 		if (on){
// 			tinput.attr('checked','checked');
// 		}
// 		let w = this;//ugh
// 		tinput.change(e=>{
// 			w.temperatures[$(this).attr("name")] = $(this).attr('checked');
// 		});
// 		let tspan = $(`<br><span name="T${t}-span">`);
// 		tspan.append(tinput);
// 		tspan.append(t);
// 		$(`[name="${this.exportGroupName}"] > .temperature-grouping`).append(tspan);
// 	}
//
// 	getDynamicModeFiles(zip){
// 		let [patchesStr, particlesStr, patchCounter, stateTransitionsStr] = genDynamicModePatchesAndParticleFiles(export_system.cube_types);
// 	   	let [topStr, confStr] = this.getTopAndConfStrs();
// 	   	let fileNames = [topFileName,
// 	   					particlesFileName,
// 	   					patchesFileName,
// 	   					...Object.keys(stateTransitionsStr).map(k=>`particle_${k}_statetransition.txt`)
// 	   					];
// 	   	let fileStrs = [topStr,
// 	   					particlesStr,
// 	   					patchesStr,
// 	   					...Object.keys(stateTransitionsStr).map(k=>stateTransitionsStr[k])
// 						];
// 		this.genSimDirs(zip, fileNames, fileStrs);
// 	}
//
// 	getStaticModeFiles(zip){
//
// 	    let [patchesStr, particlesStr, patchCounter] = genStaticModePatchesAndParticleFiles(export_system.cube_types);
// //	    const interactionMatrixFileName = 'interactions.txt'
// //
// //	    let interactionsStr = "";
// //	    for (let c1 = -export_system.countColors(); c1 <= export_system.countColors(); c1++){
// //			let color_1 = (c1 + 20) * Math.sign(c1);
// //			for (let c2 = -export_system.countColors(); c2 <= export_system.countColors(); c2++){
// //				color_2 = (c2 + 20) * Math.sign(c2);
// //				interactionsStr += `patchy_eps[${color_1}][${color_2}]=${c1==0||c2==0 ? 0 : 1}\n`
// //			}
// //		}
//
// 	   	let [topStr, confStr] = this.getTopAndConfStrs();
//
// 	    let fileNames = [topFileName,
// 	   					particlesFileName,
// 	   					patchesFileName
// 	   					];
// 	   	let fileStrs = [topStr,
// 	   					particlesStr,
// 	   					patchesStr];
// 	   	this.genSimDirs(zip, fileNames, fileStrs);
//
// 	}
//
// 	getTopAndConfStrs(){
// 		let topStr, confStr;
//
// 	    if (this.num_assemblies > 1) {
// 	        const count = this.particle_type_levels.map(n=>n*this.num_assemblies);
//
// 	        let total = 0;
// 	        let top = [];
// 	        count.forEach((c, typeID)=>{
// 	            total+=c;
// 	            for(let i=0; i<c; i++) {
// 	                top.push(typeID);
// 	            }
// 	        });
// 	        topStr = `${total} ${export_system.cube_types.length}\n${top.join(' ')}`;
// 	    } else {
// 	        [topStr, confStr] = generateTopAndConfig(export_system.cube_types, assemblyMode);
// 	    }
// 	    return [topStr, confStr];
// 	}
//
// 	genSimDirs(zip, fileNames, fileStrs){
// 	    let user = $('[name="asu-id"]').val();
// 	    let oxDNA_dir = $('[name="oxdna-dir"]').val();
// 		for (let dup = 0; dup < this.num_duplicates; dup++) {
// 	        let dupfolder = zip.folder(`${this.exportGroupName}_duplicate_${dup}`);
// 	        for (const narrow_type of $('[name="narrow-types"]').val().split(",")) {
// 	            let ntfolder = dupfolder.folder(`nt${narrow_type}`);
// 	            let temperatures = $('[name="temperatures"]').val().split(",").map(parseFloat).filter(t=>this.temperatures[t]);
// 	            for (const temperature of temperatures) {
// 	                let folder = ntfolder.folder(`T_${temperature}`);
//
// 					let inputStr = generateInputFile({
// 						'particle_types_N': export_system.cube_types.length,
// 						'patch_types_N': export_system.cube_types.map(ct=> {
// 							return ct.patches.filter(p=>p.color != 0).length
// 						}).reduce((a,b) => a + b),
// 						'narrow_type': narrow_type,
// 						'T': temperature,
// 						'topology': topFileName,
// 						'conf_file': confFileName,
// 						'patchy_file': 'patches.txt',
// 						'particle_file': 'particles.txt'
// 					});
//
// 	                let simulateStr = `addqueue -c "${this.exportGroupName} T=${temperature} nt${narrow_type} - 1 week" ${oxDNA_dir}/build/bin/oxDNA ${inputFileName}`;
// 					let slurm_args = {
// 						'p': 'sulccpu1',
// 						'q': 'sulcgpu1',
// 						'n': '1',
// 						't': '8-00:00',
// 						'job-name': `"${$('[name="export-model-name"]').val()}_${this.exportGroupName}_T${temperature}_nt${narrow_type}"`,
// 						'o': `${this.exportGroupName}_T${temperature}_nt${narrow_type}.out`,
// 						'e': `${this.exportGroupName}_T${temperature}_nt${narrow_type}.out`
// 						//'mail-type': 'BEGIN,END,FAIL',
// 						//'mail-user': `${user}@asu.edu`
// 					};
// 					let submit_slurmStr = "#!/bin/bash\n";
// 					submit_slurmStr += Object.entries(slurm_args).map((pair) => {
// 						let [k, v] = pair
// 						return "#SBATCH " + (k.length == 1 ? `-${k} ${v}` : `--${k}=${v}`);
// 					}).join('\n');
// 					submit_slurmStr += `\n
// module add cuda/10.2.89
// module add gcc/8.4.0
// echo $SLURM_JOB_ID >> jobid.txt
//
// if [ ! -f init.conf ]; then
// 	${oxDNA_dir}/build/bin/confGenerator input ${$('[name="particle-density"]').val()}
// fi
// ${oxDNA_dir}/build/bin/oxDNA ${inputFileName}
// `;
//
// 					fileNames.forEach((fn, i)=>{
// 						folder.file(fn, fileStrs[i]);
// 					});
// 	                folder.file(inputFileName, inputStr);
// 	                folder.file('simulate.sh', simulateStr);
// 	                folder.file('submit_slurm.sh', submit_slurmStr)
// 	            }
// 	        }
// 	    }
// 	}
// }
//
// function downloadPatchySimFiles() {
// 	let zip = new JSZip();
// 	export_groups.forEach(getDentalRadius() > 0 ? g=>g.getDynamicModeFiles(zip) : g=>g.getStaticModeFiles(zip));
// 	zip.file(
//         	'simulateAll.sh', `
// for var in */nt*/T_*
// do
//     cd $var
//     sbatch submit_slurm.sh
//     cd ../../..
// done`
// 	    );
// 	let setup = setupToJSON();
// 	zip.file("patchy_export_setup.json",  JSON.stringify(setup, null, 4));
// 	zip.generateAsync({type:"blob"})
//     	.then(function(content) {
//         	saveAs(content, `${$('[name="export-model-name"]').val()}.zip`); //FileSaver.js
//     });
// }
//
// function setupToJSON() {
// 	return {
// 		'export_name': $('[name="export-model-name"]').val(),
// 		'oxDNA_directory': $('[name="oxdna-dir"]').val(),
// 		'asu_id': $('[name="asu-id"]').val(),
// 		'narrow_types': $('[name="narrow-types"]').val().split(",").map(i=>parseInt(i)),
// 		'temperatures': $('[name="temperatures"]').val().split(",").map(parseFloat),
// 		'max_cubes': parseInt($('[name="max-rule-particles"]').val()),
// 		'particle_density': parseFloat($('[name="particle-density"]').val()),
// 		'dental_radius': getDentalRadius(),
// 		'cube_types': export_system.cube_types,
// 		// 'cube_types': export_system.cube_types.map(r=>{
// 		// 	return {
// 		// 		'name': r.name,
// 		// 		'patches': r.patches.map(patch=>{
// 		// 			let p = {
// 		// 				"color": patch.color,
// 		// 				"alignment": patch.alignDir,
// 		// 			};
// 		// 			if ("state_var" in patch){
// 		// 				p["state_var"]
// 		// 			}
// 		//
// 		// 			return p;
// 		// 		}),
// 		// 		'conditionals': r.conditionals,
// 		// 	}
// 		// }),
// 		'input_params': input_file_params,
// 		'export_groups': export_groups
// 	};
// }
//
// function getDentalRadius() {
// 	return parseFloat($('[name="dental-radius"]').val());
// }
//
// function downloadPatchySimFilesTemplate(){
// 	let setup_json = JSON.stringify(setupToJSON(), null, 4);
// 	saveAs(new Blob([setup_json], {type: "text/plain;charset=utf-8"}), "patchy_export_setup.json");
// }
//
// function uploadPatchySimFilesTemplate(source){
// 	var reader = new FileReader();
// 	reader.onload = function(e){
// 		let json = JSON.parse(this.result);
// 		let cubetypes = json['cube_types'].map(function(r){
// 			r = Object.assign(new StaticPolycubeCube(), r);
// 			r.patches = r.patches.map(patch=>{
// 				let p = Object.assign(new PolycubePatch(), patch);
// 				p.alignDir = Object.assign(new THREE.Vector3(), patch.alignDir);
// 				return p;
// 			})
// 			return r;
// 		});
// 		export_system = new PolycubeSystem(cubetypes, false);
// 		$('[name="export-model-name"]').val(json["export_name"]);
// 		$('[name="oxdna-dir"]').val(json["oxDNA_directory"]);
// 		$('[name="asu-id"]').val(json["asu_id"]);
// 		$('[name="narrow-types"]').val(json["narrow_types"].join(","));
// 		$('[name="temperatures"]').val(json["temperatures"].join(","));
// 		$('[name="max-rule=particles"]').val(json["max_cubes "]);
// 		$('[name="particle-density"]').val(json["particle_density"]);
// 		if ("dental_radius" in json){
// 			$('[name="dental-radius"]').val(json["dental_radius"]);
// 		}
// 		export_groups = json['export_groups'].map(eg => Object.assign(new ExportGroup(), eg));
// 		$('[name="export-rules"]').val(ruleToDec(export_system.cube_types));
// 		export_groups.forEach(eg=>{
// 			$('#groups').append(eg.visualize());
// 		});
// 		updateInputFileInput(json["input_params"]);
// 	};
// 	reader.readAsText(source.files[0])
// }
//
// // this function is from w3schools, unclear why it isn't built in
// function getCookie(cname) {
// 	  let name = cname + "=";
// 	  let decodedCookie = decodeURIComponent(document.cookie);
// 	  let ca = decodedCookie.split(';');
// 	  for(let i = 0; i <ca.length; i++) {
// 	    let c = ca[i];
// 	    while (c.charAt(0) == ' ') {
// 	      c = c.substring(1);
// 	    }
// 	    if (c.indexOf(name) == 0) {
// 	      return c.substring(name.length, c.length);
// 	    }
// 	  }
// 	  return "";
// 	}
//
// function updateCookies(){
//
// 	document.cookie = `oxDNA_home=${$('[name="oxdna-dir"]').val()};`;
// 	document.cookie = `asu_id=${$('[name="asu-id"]').val()};`;
// }
//
// function updateTemperatures(){
// 	export_groups.forEach(g=>{
// 		g.updateTemperatures();
// 	})
// }
//
// function updateNarrowTypes(source){
// 	export_groups.forEach(eg=>{
// 		eg.narrow_types = $(source).val().split(",").filter(Boolean).map(parseInt);
// 	})
// }
//
// function updateParticleDensity(source) {
// 	export_groups.forEach(eg=>{
// 		eg.particle_density = $(this).val();
// 	})
// }
//
// function addExportGroup(){
// 	if (export_system == void 0){
// 		log("No system used for export.");
// 		return;
// 	}
// 	let eg = new ExportGroup(`EG${export_groups.length}`);
// 	export_groups.push(eg);
// 	$("#groups").append(eg.visualize());
// }
//
// function updateExportSystem(){
// 	if (export_system == void 0){
// 		export_system = new PolycubeSystem(parseDecRule($('[name="export-rules"]').val()), false);
// 	}
// 	else {
// 		let old_num_rules = export_system.cube_types.length;
// 		export_system.rules = parseDecRule($('[name="export-rules"]').val());
// 		if (export_system.cube_types.length != old_num_rules) {
// 			$("#groups").empty();
// 			export_groups.forEach(g=>{
// 				$("#groups").append(g.visualize());
// 			});
// 		}
// 	}
// }
//
// function updateMaxRules(){
// 	$('input.rule-level').attr("max", $('[name="max-rule-particles"]').val());
// }
//
// function updateInputFileInput(input_params){
// 	input_file_params = input_params;
// 	Object.keys(input_params["OBSERVABLES"]).forEach((key)=>{
// 		input_file_params["OBSERVABLES"][key] = Object.assign(new Observable(), input_params["OBSERVABLES"][key]);
// 	});
// 	let input_params_div = $("#input-file-view");
// 	input_params_div.empty();
// 	input_params_div.append($("<legend>Input File</legend>"));
// 	Object.keys(input_params).forEach(key=>{
// 		input_params_div.append($(`<h1>${key}</h1>`));
// 		let paramslist = $(`<ul id="${key.toLowerCase().replaceAll(/[\s*\/]/g, '-')}-params-list" class="input-file-params-list"></ul>`);
// 		Object.keys(input_params[key]).forEach(paramName => {
// 			let item = $(`<li class="input-file-param-list-item"><label for="${paramName}">${paramName}</label></li>"`);
// 			if (input_params[key][paramName] != null){
// 				let ipt = $(`<input type="text" id="param-${paramName}"/>`);
// 				ipt.change(function (){
// 					input_file_params[key][paramName] = this.value;
// 				});
// 				ipt.val(input_params[key][paramName]);
// 				item.append(ipt);
// 			}
// 			else {
// 				item.append($("<span>Dependant</span>"));
// 			}
// 			paramslist.append(item);
// 		});
// 		input_params_div.append(paramslist);
// 	});
// }
//
// function generateInputFile(dependancies={}){
// 	// TODO: better parameter name than "dependencies
// 	// loop categories
// 	let str = "";
//
// 	Object.entries(input_file_params).forEach(([key,value]) => {
// 		str += "################################\n";
// 		str += "####" + " ".repeat((24 - key.length)/2) + key + " ".repeat((24 - key.length)/2) + "####\n"
// 		str += "################################\n";
// 		for (const [k2, v2] of Object.entries(value)) {
// 			let val = v2;
// 			if (val == null){
// 				console.assert(dependancies[k2] != undefined);
// 				val = dependancies[k2];
// 			}
// 			if (v2 instanceof Observable){
// 				str += `${k2} = ${val.to_string()}\n`
// 			}
// 			else {
// 				str += `${k2} = ${val}\n`
// 			}
// 		}
// 	});
// 	return str;
// }
//
// class Observable{
// 	constructor(name, print_every, cols) {
// 		this.obsName = name;
// 		this.print_every = print_every;
// 		this.cols = cols;
// 	}
//
// 	to_string(){
// 		// TODO: make more generic
// 		return `{
// 	name = ${this.obsName}
// 	print_every = ${this.print_every}
// 	col_1 = {
// 		type = ${cols[0].type}
// 		show_types = ${cols[0].show_types}
// 	}
// }`
// 	}
// }
//
//
// function generateTopAndConfig(rule, assemblyMode='seeded') {
//     let sys = new PolycubeSystem(rule, undefined, 100, 100, assemblyMode, true);
//     sys.seed();
//     let processed;
//     while (!sys.isProcessingComplete()) {
//         processed = sys.processMoves();
//         if (processed == 'oub') {
//             console.warn("Getting config for unbounded rule");
//             break;
//         }
//     }
//     let top = [];
//     let conf = [];
//     let nParticles = 0;
//     let max = new THREE.Vector3();
//     let min = new THREE.Vector3();
//     let mean = new THREE.Vector3();
//     for (const [key, c] of sys.confMap) {
//         const p = sys.cubeMap.get(key);
//         for (let i=0; i<3; i++) {
//             max.setComponent(i, Math.max(max.getComponent(i), p.getComponent(i)));
//             min.setComponent(i, Math.min(min.getComponent(i), p.getComponent(i)));
//         }
//         mean.add(p);
//         nParticles++;
//     }
//     let box = max.clone().sub(min).multiplyScalar(3);
//     mean.divideScalar(nParticles);
//
//     for (const [key, c] of sys.confMap) {
//         const p = sys.cubeMap.get(key).sub(mean).add(box.clone().divideScalar(2));
//         const a1 = new THREE.Vector3(1, 0, 0).applyQuaternion(c.q);
//         const a3 = new THREE.Vector3(0, 0, 1).applyQuaternion(c.q);
//         const vF = (v) => v.toArray().map(
//             n=>Math.round(n*100)/100 // Round to 2 decimal places
//         ).join(' ');
//         conf.push(`${vF(p)} ${vF(a1)} ${vF(a3)} 0 0 0 0 0 0`);
//         top.push(`${c.ruleIdx}`);
//     }
//     let confStr = `t = 0\nb = ${box.toArray().join(' ')}\nE = 0 0 0\n` + conf.join('\n');
//     let topStr = `${nParticles} ${rule.length}\n` + top.join(' ');
//     return [topStr, confStr];
// }
//
// function getPatchString(patch){
// 	let patchDirection = RULE_ORDER[patch.idx].clone();
// 	let dentalRadius = getDentalRadius();
// 	if (dentalRadius == 0){
//         return [
//             `patch_${patch.id} = {`,
//             `  id = patch_${patch.id}`,
//             `  color = ${patch.color}`,
//             `  strength = ${patch.strength}`,
//             `  position = ${patchDirection.clone().divideScalar(2).toArray()}`,
//             `  a1 = ${patchDirection.toArray()}`,
//             `  a2 = ${patch.align.toArray()}`,
//             `  allostery_conditional = ${patch.allosteric_control}`,
//             '}', ''
//         ].join('\n');
// 	 }
// 	 else {
// 		// clone patch position vector
// 		let patchPosition = patchDirection.clone().divideScalar(2);
// 		// get patch alignment
// 		let patchAlign = patch.align.clone();
// 		// get rotation index (have to rotate in opposite directions for positive and negative colors)
// 		let rotIdx = (patch.color > 0 ? patch.toothIdx : 3 - patch.toothIdx);
// 		// apply tooth-appropriate rotation to patch alignment
// 		patchAlign.applyAxisAngle(patchDirection, rotIdx * Math.PI / 2);
// 		// apply alignment/rotation to position
// 		patchPosition = patchPosition.add(patchAlign.multiplyScalar(dentalRadius));
// 		// apply alignment/rotation to a1
// 		a1 = patchDirection.clone().add(patchAlign.multiplyScalar(dentalRadius));
// 		let lines = [
//             `patch_${patch.id} = {`,
//             `  id = patch_${patch.id}`,
//             `  color = ${patch.color}`,
//             `  strength = ${patch.strength}`,
//             `  position = ${patchPosition.toArray().map(x=>x.toFixed(3))}`,
//             `  a1 = ${a1.toArray().map(x=>x.toFixed(3))}`];
// 		if ('activation_var' in patch){
// 			lines.push(`  activation_var = ${patch.activation_var}`);
// 		}
// 		if ('state_var' in patch){
// 			lines.push(`  state_var = ${patch.state_var}`);
// 		}
//         return [...lines,'}', ''].join('\n');
// 	}
// }
//
// function getParticleString(typeID, patches, stateSize=0) {
// 	lines = [
//         `particle_${typeID} = {`,
//         `  type = ${typeID}`,
//         `  patches = ${patches}`];
//    if (stateSize > 1){
// 		lines.push(`state_size = ${stateSize}`);
// 	}
//    return [...lines,
//         '}',''
//     ].join('\n');
// }
//
// function getPatchySimFiles(rule, user, nAssemblies=1, name='sim',
//     oxDNA_dir = '/home/jrevan21/oxDNA_torsion',
//     temperatures = ['0.01'],
//     confDensity = 0.2,
//     narrow_types=['0'],
//     duplicates = 1,
// ) {
//     let zip = new JSZip();
//
//     let particlesStr = "";
//     let patchesStr = "";
//
//     let patchCounter = 0;
//     if (getDentalRadius() == 0){
// 	    [patchesStr, particlesStr, patchCounter] = genStaticModePatchesAndParticleFiles(rule);
// 	}
// 	else {
// 		[patchesStr, particlesStr, patchCounter] = genDynamicModePatchesAndParticleFiles(rule);
// 	}
//
//     const particlesFileName = 'particles.txt'
//     const patchesFileName = 'patches.txt'
//
//     const cubeTypeCount = getCubeTypeCount(rule);
//
//     let topStr, confStr;
//     const topFileName = 'init.top';
//     const confFileName = 'init.conf';
//     if (nAssemblies > 1) {
//         const count = cubeTypeCount.map(n=>n*nAssemblies);
//
//         let total = 0;
//         let top = [];
//         count.forEach((c, typeID)=>{
//             total+=c;
//             for(let i=0; i<c; i++) {
//                 top.push(typeID);
//             }
//         });
//         topStr = `${total} ${rule.length}\n${top.join(' ')}`;
//     } else {
//         [topStr, confStr] = generateTopAndConfig(rule, assemblyMode);
//     }
//
//     const inputFileName = 'input';
//
//     for (let dup = 0; dup < duplicates; dup++) {
//         let dupfolder = zip.folder(`duplicate_${dup}`);
//         for (const narrow_type of narrow_types) {
//             let ntfolder = dupfolder.folder(`nt${narrow_type}`);
//             for (const temperature of temperatures) {
//                 let folder = ntfolder.folder(`T_${temperature}`);
//
// 				let inputStr = generateInputFile({
// 					'particle_types_N': rule.length,
// 					'patch_types_N': patchCounter,
// 					'narrow_type': narrow_type,
// 					'T': temperature,
// 					'topology': topFileName,
// 					'conf_file': confFileName
// 				});
//
//                 // let inputStr = generateInputFile(rule.length, patchCounter, oxDNA_dir, patchesFileName,
//                 //     particlesFileName, topFileName, temperature, narrow_type
//                 // );
//
//                 let simulateStr = `addqueue -c "${name} T=${temperature} nt${narrow_type} - 1 week" ${oxDNA_dir}/build/bin/oxDNA ${inputFileName}`;
// 				let slurm_args = {
// 					'p': 'sulccpu1',
// 					'q': 'sulcgpu1',
// 					'n': '1',
// 					't': '8-00:00',
// 					'job-name': `"${name}_T${temperature}_nt${narrow_type}"`,
// 					'o': `${name}_T${temperature}_nt${narrow_type}.out`,
// 					'e': `${name}_T${temperature}_nt${narrow_type}.out`
// 					//'mail-type': 'BEGIN,END,FAIL',
// 					//'mail-user': `${user}@asu.edu`
// 				};
// 				let submit_slurmStr = "#!/bin/bash\n";
// 				submit_slurmStr += Object.entries(slurm_args).map((pair) => {
// 					let [k, v] = pair
// 					return "#SBATCH " + (k.length == 1 ? `-${k} ${v}` : `--${k}=${v}`);
// 				}).join('\n');
// 				submit_slurmStr += `\n
// module add cuda/10.2.89
// module add gcc/8.4.0
// if [ ! -f init.conf ]; then
// 	${oxDNA_dir}/build/bin/confGenerator input ${confDensity}
// fi
// ${oxDNA_dir}/build/bin/oxDNA ${inputFileName}`;
//
//                 folder.file(particlesFileName, particlesStr);
//                 folder.file(patchesFileName, patchesStr);
//                 folder.file(topFileName, topStr);
//                 folder.file(inputFileName, inputStr);
//                 folder.file('simulate.sh', simulateStr);
//                 folder.file('submit_slurm.sh', submit_slurmStr)
//
//                 if (confStr) {
//                     folder.file(confFileName, confStr);
//                 }
//             }
//         }
//     }
//
//     zip.file(
//         'simulateAll.sh', `
// for var in */nt*/T_*
// do
//     cd $var
//     sbatch submit_slurm.sh
//     cd ../../..
// done`
//     );
//
//     let complClusters = cubeTypeCount.flatMap((n,idx) => {
//         let cluster = [];
//         for (let i=0; i<n; i++) {
//             cluster.push(idx)
//         }
//         return cluster;
//     });
//
//     zip.file('data.py',
// `completeCluster = [${complClusters}]
// rule = '${ruleToDec(rule)}'`
//     )
//
//     zip.generateAsync({type:"blob"})
//     .then(function(content) {
//         saveAs(content, `${name}.zip`); //FileSaver.js
//     });
// }
//
// function genStaticModePatchesAndParticleFiles(rule){
// 	let patchesStr = "";
// 	let particlesStr = "";
// 	let patchCounter = 0;
// 	rule.forEach((cubeType, typeID)=>{
//         let patches = [];
//         let internal_index = 0;
//         cubeType.patches.forEach((patch, i)=>{
//             if(patch.color != 0) {
//                 // Needs to be > 20 to not be self complementary
//                 let color = patch.color + 20 * Math.sign(patch.color);
// //                patchesStr += getPatchStr(patchCounter, color, i, patch.alignDir);
//                 patches.push({
// 					'id': patchCounter,
// 					'idx': i,
// 					'internal_index': internal_index,
// 					'color': color,
// 					'align': patch.alignDir,
// 					'strength': 1
// 				});
// 				internal_index++;
//                 patchCounter++;
//             }
//         });
//         patches.forEach(p=>{
// 			p.allosteric_control = cubeType.conditionals[p.idx]
// 			p.allosteric_control = p.allosteric_control.slice(1,-1)
// 			p.allosteric_control = p.allosteric_control.replaceAll(/b\[\d+\]/g, s=>{
// 					let ctrl_idx = parseInt(s.slice(s.indexOf('[')+1, s.indexOf(']')));
// 					return patches.find(pt=>pt.idx == ctrl_idx).internal_index
// 				}
// 			);
// 			if (!p.allosteric_control.trim()){
// 				p.allosteric_control = "(true)";
// 			}
// 		});
//         patchesStr += patches.map(getPatchString).join("\n");
//         particlesStr += getParticleString(typeID, patches.map(p=>p.id));
//     });
//     return [patchesStr, particlesStr, patchCounter];
// }
//
// const STATE_TRANSITION_SUBDIV = 4;
// function genDynamicModePatchesAndParticleFiles(rule) {
// 	let patchesStr = "";
// 	let particlesStr = "";
// 	let stateTransitionMaps = {};
// 	let patchCounter = 0;
//
// 	const toothCount = 4;
// 	rule.forEach((cubeType, typeID)=>{
// 		let activationVarCounter = 0;
// 		let activationConditionals = []; // indexed list of conditionals controlling activation vars
// 		let stateVarCounter = 1;
// 		let stateVarMap = {}; // maps patches to state vars
//         let patches = [];
//         let internal_index = 0;
//         cubeType.patches.forEach((patch, i)=>{
//             if(patch.color != 0) {
// 				let allosteric_control = cubeType.conditionals[i];
//                 // Needs to be > 20 to not be self complementary
// //                patchesStr += getPatchStr(patchCounter, color, i, patch.alignDir);
// 				[...Array(toothCount).keys()].forEach(x=>{
// 					let color = patch.color;
// 					color = color > 0 ? toothCount * color + x : toothCount * color - x;
// 					color += 20 * Math.sign(color);
// 					let patchDict = {
// 						'id': patchCounter,
// 						'idx': i,
// 						'internal_index': toothCount * internal_index + x,
// 						'color': color,
// 						'align': patch.alignDir,
// 						'toothIdx': x,
// 						'strength': 1 / toothCount
// 					};
// 					if (allosteric_control != "(true)"){
// 						patchDict['activation_var'] = activationVarCounter;
// 					}
// 					patches.push(patchDict);
// 					internal_index++;
// 	                patchCounter++;
// 				});
// 				if (allosteric_control != "(true)"){
// 					// index activation vars
// 					activationConditionals.push(allosteric_control);
// 					activationVarCounter++;
//
// 					// index state vars
// 					[... allosteric_control.matchAll(/b\[\d+\]/g)].forEach(p=>{
// 						let sz = p.toString();
// 						let ctrl_idx = parseInt(sz.slice(sz.indexOf('[')+1, sz.indexOf(']')));
// 						if (!(ctrl_idx in stateVarMap)){
// 							stateVarMap[ctrl_idx] = stateVarCounter++;
// 						}
// 					});
// 				}
//             }
//         });
//
//         // assign state vars
//         patches.forEach(p=>{
// 			if (p.idx in stateVarMap){
// 				p['state_var'] = stateVarMap[p.idx];
// 			}
// 		});
//
//         // reindex activation conditionals
//        	let state_size = stateVarCounter + activationVarCounter;
//         if (state_size > 1){
// 	        let numStates = Math.pow(2, state_size);
//
// 	        let stateTransitionList = [];
// 	        // loop states, discarding states where var 0 is false b/c var 0 is the tautology variable
// 	        [...Array(numStates).keys()].forEach(statenum=>{
// 				// ignore states where state[0] == 0
// 			/*	if (statenum % 2 == 0){
// 					return;
// 				}*/
// 				let varActivations = [];
// 				// create binding array that corresponds to state
// 				let b = [...Array(RULE_ORDER.length).keys()].map(i => {
// 					// if patch i corresponds to a state variable
// 					if (i in stateVarMap){
// 						return statenum / Math.pow(2, stateVarMap[i]) >= 1;
// 					}
// 				});
// 				// default target state is current state (no change)
// 				let targetState = statenum;
// 				// loop activation conditionals
// 				activationConditionals.forEach((c, i) => {
// 					let newActivation = eval(c);
// 					// if activation is true
// 					if (newActivation) {
// 						// set bit at i + stateVarCounter to 1
// 						targetState |= 1 << stateVarCounter + i;
// 					}
// 					else {
// 						// set bit at i + stateVarCounter to 0
// 						targetState &=  ~(1 << stateVarCounter + i);
// 					}
// 				});
// 				var stateLine = Array(STATE_TRANSITION_SUBDIV).fill(targetState).join(" ");
// 				// add a comment specifying the actual state
//
// 				var state = [...Array(state_size).keys()].map(i => (statenum >> i) % 2);
// 				stateLine += " # state=[" + state.join(" ") + "]";
// 				stateTransitionList.push(stateLine);
// 			});
// 			stateTransitionMaps[typeID] = stateTransitionList.join("\n");
//
// 		}
//
//         // replace tiny values with 0 so as to not screw up oxview
//         patchesStr += patches.map(getPatchString).join("\n");
//         particlesStr += getParticleString(typeID, patches.map(p=>p.id), state_size);
//     });
//     return [patchesStr, particlesStr, patchCounter, stateTransitionMaps];
// }
//
//
// const DEFAULT_INPUT_FILE_PARAMS = {
// 	"PROGRAM_PARAMETERS": {
// 		"backend": "CPU",
// 		"backend_precision": "double",
// 		"ensemble": "NVT",
// 		"delta_translation": 0.1,
// 		"delta_rotation": 0.1,
// 		"narrow_type": null
// 	},
// 	"SIM_PARAMETERS": {
// 		"newtonian_steps": 103,
// 		"diff_coeff": 0.1,
// 		"thermostat": "john",
// 		"sim_type": "MD",
// 		"dt": 0.001,
// 		"verlet_skin": 0.05,
// 		"no_stdout_energy": 0,
// 		"restart_step_counter": 1,
// 		"energy_file": "energy.dat",
// 		"print_conf_interval": 1e7,
// 		"print_energy_every": 1e5,
// 		"time_scale": "linear",
// 		"PATCHY_alpha": 0.12,
// 		"no_multipatch": 1,
// 		"steps": 5e9,
// 		"check_energy_every": 10000,
// 		"check_energy_threshold": 1.e-4,
// 		"T": null,
// 		"refresh_vel": 1
// 	},
// 	"PATCHY_SETUP": {
// 		"interaction_type": "PatchyShapeInteraction",
// 		"plugin_search_path": "/home/jrevan21/alloassembly/contrib/romano",
// 		"shape": "sphere",
//
// 		"particle_types_N": null,
// 		"patch_types_N":  null,
//
// 		"patchy_file": null,
// 		"particle_file": null,
// 		"same_type_bonding": 1,
// 		"use_torsion": 1,
// 		"interaction_tensor": 0,
//
// 		"PATCHY_radius": 0.5,
//
// 		"lastconf_file": "last_conf.dat"
// 	},
// 	"INPUT_OUTPUT": {
// 		"topology": null,
// 		"conf_file": null,
// 		"trajectory_file": "trajectory.dat"
// 	},
// 	"OBSERVABLES": {
// 		"data_output_1": new Observable(
// 			name="clusters.txt",
// 			print_every= 1e9,
// 			cols=[{
// 				type: "PLClusterTopology",
// 				show_types: 1
// 			}]
// 		)
// 	}
// }
//
// updateInputFileInput(DEFAULT_INPUT_FILE_PARAMS);