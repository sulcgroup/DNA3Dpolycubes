import * as SVG from "@svgdotjs/svg.js";
import * as $ from "jquery";
import { Queue } from "queue-typescript";
import { FACE_ROTATIONS, RULE_ORDER, signed_axis_angle } from "./utils";
import { FACE_ROT_OFFSETS } from "./edit";
const LETTERS = "ABCDFGHJKMPZ"; // skip "I" and "L" because they look too much like "J",
// "F" because it looks too much like "E"
let LCUT_TAB_LABELS = [
    // edges are in order top, right, bottom, left
    // indices will be used to index `LETTERS`
    [5, 11, 9, 7], //left
    [6, 3, 2, 10], //right
    [8, 11, 1, 10], //bottom
    [0, 7, 4, 3], //top
    [4, 9, 8, 2], //back
    [1, 5, 0, 6] //front
];
const PX_PER_CM = 50;
const LABEL_TEXT_FONT = { family: 'Helvetica', size: 24, style: "bold", anchor: "middle" };
const TAB_LABEL_FONT = { family: 'Helveticad', style: "bold", size: 18, anchor: 'middle' };
// default laser cut settings
const CUT_DEFAULTS = {
    // cardboard we currently use is 12 in x 12 in
    board_w_cm: 30.48,
    board_h_cm: 30.48,
    board_thickness_cm: 0.635,
    face_w_cm: 8,
    patch_w_cm: 4.7,
    patch_h_cm: 3.15,
    n_tabs: 2
};
// settable
let settings = CUT_DEFAULTS;
function canvasDimensions() {
    return [settings.board_w_cm * PX_PER_CM, settings.board_h_cm * PX_PER_CM];
}
function faceCounterPerBoard() {
    return [Math.floor(settings.board_w_cm / (settings.face_w_cm - settings.board_thickness_cm)),
        Math.floor(settings.board_h_cm / (settings.face_w_cm - settings.board_thickness_cm))];
}
function cubeFaceInnerW() {
    /**
     * @return the width of the inner portion of a cube face (sans tabs) in pixels
     */
    return (settings.face_w_cm - 2 * settings.board_thickness_cm) * PX_PER_CM;
}
/**
 * opens a dialogue box to display lasercutter templates
 * @param cubes cubes to draw
 */
export function drawCutterWindow(cubes) {
    // construct laser cutter window
    $("body").append($("<div id='lcut-window'>"));
    let lcutWindow = $("#lcut-window");
    // add preview window
    lcutWindow.append($("<div id='lcut-preview'>"));
    $("#lcut-preview").append($("<div id='lcut-preview-inner' style='padding: 0 0 0; margin: 0 0 0 0; height: 100%'>"));
    // generate actual designs
    let boardCount = 0;
    for (let it of getCubesLCutSvg(cubes)) {
        $("#lcut-preview-inner").append(it.node);
        boardCount++;
    }
    // add a few labels to describe it
    // num boards
    lcutWindow.append($(`<p><b>Num. Boards</b>:${boardCount}</p>`));
    // counts of each color (reqd to print patches
    let cubeColorCounts = new Map();
    for (let cube of cubes) {
        cube.iter_patches((patch) => {
            if (!cubeColorCounts.has(patch.patch.color)) {
                cubeColorCounts.set(patch.patch.color, 0);
            }
            cubeColorCounts.set(patch.patch.color, cubeColorCounts.get(patch.patch.color) + 1);
        });
    }
    let color_info = $("<p>");
    lcutWindow.append(color_info);
    for (let [color, count] of cubeColorCounts.entries()) {
        color_info.append($(`<p style="display: inline-block;"><b>Color ${color}</b>: ${count}</p>`));
    }
    // upload cut settings button
    lcutWindow.append($("<button id='lcut-upload-button' class='lcut-ui-btn'>Upload Cut Settings</button>"));
    lcutWindow.append($(`<input type="file"
                   accept="json"
                   style="display:none;"
                   id="upload-lcut-settings">`)); // ghost input
    // trigger our ghost input when button is pressed
    $("#lcut-upload-button").on('click', (e) => {
        $("#upload-lcut-settings").trigger('click');
    });
    // add handler for input
    $("#upload-lcut-settings").on('change', (e) => {
        var reader = new FileReader();
        reader.onload = function (e) {
            if (typeof this.result === "string") {
                settings = JSON.parse(this.result);
            }
            else {
                console.log("Malformed laser cutter settings json");
            }
        };
        reader.readAsText(e.target.files[0]);
        // close lasercutter window
        $("#lcut-close-button").trigger('click');
        // and reopen it
        $("#download-lcut").trigger('click');
    });
    // download cut files button
    lcutWindow.append($("<button id='lcut-download-button' class='lcut-ui-btn'>Download Cut Files</button>"));
    $("#lcut-download-button").on('click', (e) => {
        // loop svg elements
        let i = 1; // only to name files
        for (let canvas of $("#lcut-preview-inner").children()) {
            const blob = new Blob([canvas.outerHTML], { type: "image/svg+xml;charset=utf-8" }); // Create a blob from the SVG data
            const url = URL.createObjectURL(blob); // Create a URL for the blob
            // Create a download link for the SVG
            const downloadLink = document.createElement('a');
            downloadLink.href = url;
            downloadLink.download = `lasercutter_cutouts_pt${i}.svg`; // Name the file
            document.body.appendChild(downloadLink);
            downloadLink.click(); // Programmatically click the link to trigger the download
            document.body.removeChild(downloadLink); // Clean up
            i++;
        }
        // close laser cutter window
        $("#lcut-close-button").trigger('click');
    });
    // button to close this window
    lcutWindow.append($("<button id='lcut-close-button' class='lcut-ui-btn'>Close</button>"));
    $("#lcut-close-button").on('click', (e) => {
        // remove lasercutter window element
        $("#lcut-window").remove();
    });
    // done generating controls
}
/**
 * Takes a list of cubes and yields SVG elements until it's displayed all cubes
 * each svg element is one material board
 * @param cubes
 * @param numTabs
 */
let STROKE_CUT = {
    // width: .001,
    width: .001,
    color: "red"
};
let STROKE_DRAW = {
    // width: .001,
    width: .001,
    color: "blue"
};
function* getCubesLCutSvg(cubes) {
    let [nfaces_x, nfaces_y] = faceCounterPerBoard();
    console.assert(nfaces_x > 0);
    console.assert(nfaces_y > 0);
    // compute width and height
    let svgObj;
    // one-dimensional canvas indexer, flattened from 2D indexer
    // index in canvas, should be 0 <= canvasIdx < ncubes_x * ncubes_y
    let canvasIdx = 0;
    let [canvasWidth, canvasHeight] = canvasDimensions();
    let svgDimensions = {
        height: '95%', // Fill the available height
        viewBox: `0 0 ${canvasWidth} ${canvasHeight}`, // Set the viewBox
        preserveAspectRatio: 'xMidYMid meet', // Maintain aspect ratio, centered
        display: 'inline-block'
    };
    svgObj = SVG.SVG();
    svgObj.attr(svgDimensions);
    let cubes_queue = new Queue(...cubes);
    // continue doing this until we run out of cubes
    while (cubes_queue.length > 0) {
        let cube = cubes_queue.dequeue();
        // can't use forEach here bc generator bullshit
        for (let patchIdx = 0; patchIdx < RULE_ORDER.length; patchIdx++) {
            // compute x and y index coordinates in canvas
            let patch_idx_x = canvasIdx % nfaces_x;
            let patch_idx_y = Math.floor(canvasIdx / nfaces_y);
            // Draw the rectangle with finger joints + labels
            let face_group = drawCubeLCutFace(cube, patchIdx, svgObj);
            face_group.translate(
            // i can't figure out why the x translation should be like this
            (patch_idx_x * (settings.face_w_cm - settings.board_thickness_cm)) * PX_PER_CM, (patch_idx_y * (settings.face_w_cm - settings.board_thickness_cm) + settings.board_thickness_cm) * PX_PER_CM);
            canvasIdx++; // increment positional 1d indexer on canvas
            // if we just passed off the canvas,
            if (canvasIdx == nfaces_y * nfaces_x) {
                yield svgObj;
                svgObj = SVG.SVG();
                svgObj.attr(svgDimensions);
                canvasIdx = 0; // reset canvas index
            }
        }
    }
    if (canvasIdx != nfaces_x * nfaces_y) {
        // yield partially filled cut template at end, but don't double count if last template is filled
        yield svgObj;
    }
}
// Function to generate the path data for a single square
/**
 * Function to generate the path data for a single square
 * Use the inner (minus tabs) width for the square and extrude tabs outward from that
 * @param cube cube to draw face of
 * @param patch_idx index of patch in RULE_ORDER
 * @param canvas svg canvas or upper-level container to use as parent in svg doc heirarchy
 */
function drawCubeLCutFace(cube, patch_idx, canvas) {
    let pathData = "M 0,0 ";
    let tab_labels = LCUT_TAB_LABELS[patch_idx].map(i => LETTERS[i]);
    const tabHeight = settings.board_thickness_cm * PX_PER_CM;
    const tabWidth = cubeFaceInnerW() / settings.n_tabs;
    const slotAdjustment = tabWidth * 0.0025;
    let face_group = canvas.group();
    // top edge
    for (let i = 0; i < settings.n_tabs; i++) {
        let x = tabWidth / 2 + slotAdjustment;
        if (i == 0) {
            x -= slotAdjustment / 2;
        }
        pathData += `l ${x},0 `; // move right
        pathData += `l 0,${-tabHeight} `; // extrude tab
        pathData += `l ${tabWidth / 2 - slotAdjustment},0 `; // move right
        pathData += `l 0,${tabHeight} `; // intrude tab
        // add tab label text
        face_group.text(tab_labels[0])
            .center(tabWidth * (i + 0.75), 0)
            .fill('none')
            .stroke(STROKE_DRAW)
            .font(TAB_LABEL_FONT);
    }
    // right edge
    for (let i = 0; i < settings.n_tabs; i++) {
        let y = tabWidth / 2 + slotAdjustment;
        if (i == 0) {
            y -= slotAdjustment / 2;
        }
        pathData += `l 0,${y} `; // move down
        pathData += `l ${tabHeight},0 `; // extrude tab
        pathData += `l 0,${tabWidth / 2 - slotAdjustment} `; // move down
        pathData += `l ${-(tabHeight)},0 `; // intrude tab
        // add tab label text
        face_group.text(tab_labels[1])
            .center(cubeFaceInnerW(), tabWidth * (i + 0.75))
            .fill('none')
            .stroke(STROKE_DRAW)
            .font(TAB_LABEL_FONT);
    }
    // bottom edge
    for (let i = 0; i < settings.n_tabs; i++) {
        let x = tabWidth / 2 + slotAdjustment;
        if (i == 0) {
            x -= slotAdjustment / 2;
        }
        pathData += `l -${x},0 `; // move left
        pathData += `l 0,${tabHeight} `; // extrude tab
        pathData += `l -${tabWidth / 2 - slotAdjustment},0 `; // move left
        pathData += `l 0,${-(tabHeight)} `; // intrude tab
        // add tab label text
        face_group.text(tab_labels[2])
            .center(cubeFaceInnerW() - tabWidth * (i + 0.75), cubeFaceInnerW())
            .fill('none')
            .stroke(STROKE_DRAW)
            .font(TAB_LABEL_FONT);
    }
    // left edge
    for (let i = 0; i < settings.n_tabs; i++) {
        let y = tabWidth / 2 + slotAdjustment;
        if (i == 0) {
            y -= slotAdjustment / 2;
        }
        pathData += `l 0,-${y} `; // move up
        pathData += `l -${tabHeight},0 `; // extrude tab
        pathData += `l 0,-${tabWidth / 2 - slotAdjustment} `; // move up
        pathData += `l ${tabHeight},0 `; // intrude tab
        // add tab label text
        face_group.text(tab_labels[3])
            .center(0, cubeFaceInnerW() - tabWidth * (i + 0.75))
            .fill('none')
            .stroke(STROKE_DRAW)
            .font(TAB_LABEL_FONT);
    }
    // end tabs path string
    pathData += "Z";
    // add face border path
    face_group.path(pathData)
        .fill('none')
        .stroke(STROKE_CUT);
    // need to have rotation even for non color faces so we know where to write cube text
    let rotation = 0;
    // if patch has color
    let patch = cube.patch(patch_idx);
    if (patch.color != 0) {
        // rotation is about to become very tricky
        const patch_angle = -signed_axis_angle(FACE_ROTATIONS[patch_idx], patch.alignDir, RULE_ORDER[patch_idx]);
        const face_rot_idx = patch_angle * (2 / Math.PI);
        rotation = (1 + face_rot_idx + FACE_ROT_OFFSETS[patch_idx]) * 90;
        // need to orient text labels later so they don't overlap cutout
        // construct cutout for interaction site surface
        face_group.ellipse(settings.patch_w_cm * PX_PER_CM, settings.patch_h_cm * PX_PER_CM)
            .move(cubeFaceInnerW() / 2 - settings.patch_w_cm / 2 * PX_PER_CM, cubeFaceInnerW() / 2 - settings.patch_h_cm / 2 * PX_PER_CM)
            .rotate(rotation)
            .stroke(STROKE_CUT)
            .fill("none");
        // draw text slightly above the ellipse `Color: ${patch_color}`
        face_group.text(`Color: ${patch.color}`)
            .center(settings.face_w_cm * PX_PER_CM / 2, settings.board_thickness_cm * PX_PER_CM + LABEL_TEXT_FONT.size * 1.25)
            .rotate(rotation + 180, // use overall face w here, not inner
        (settings.face_w_cm - settings.board_thickness_cm) * PX_PER_CM / 2, (settings.face_w_cm - settings.board_thickness_cm) * PX_PER_CM / 2)
            .stroke(STROKE_DRAW)
            .fill('none')
            .font(LABEL_TEXT_FONT);
    }
    // draw cube name
    face_group.text(cube.getPersonalName())
        .center(settings.face_w_cm * PX_PER_CM / 2, settings.board_thickness_cm * PX_PER_CM + LABEL_TEXT_FONT.size * 1.25)
        .rotate(rotation, (settings.face_w_cm - settings.board_thickness_cm) * PX_PER_CM / 2, (settings.face_w_cm - settings.board_thickness_cm) * PX_PER_CM / 2)
        .stroke(STROKE_DRAW)
        .fill("none")
        .font(LABEL_TEXT_FONT);
    return face_group;
}
//# sourceMappingURL=tangible_export.js.map