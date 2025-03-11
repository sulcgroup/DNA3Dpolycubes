// IMPORTANT: THESE ARE VERY DIFFERENT FROM JOAKIM'S ASSEMBLY METHODS
import { shuffleArray, strToVec } from "../utils";
export var AssemblyMethod;
(function (AssemblyMethod) {
    // Stochastic method includes one-pot and non-stepwise staging
    // it assumes that each stage continues until no moves remain
    AssemblyMethod[AssemblyMethod["Stochastic"] = 0] = "Stochastic";
    // stepwise-staging method is used for piece construction, may have future applications
    // it assumes that stage numbers are also step numbers
    AssemblyMethod[AssemblyMethod["StepwiseStaging"] = 1] = "StepwiseStaging";
})(AssemblyMethod || (AssemblyMethod = {}));
export class AssemblyComponentList {
    [Symbol.iterator]() {
        let counter = 0;
        return {
            next: () => {
                if (counter < this._items.length) {
                    return {
                        done: false,
                        value: this._items[counter++]
                    };
                }
                else {
                    return {
                        done: true,
                        value: undefined
                    };
                }
            }
        };
    }
    index;
    _items;
    constructor(idx, items = []) {
        this.index = idx;
        this._items = items;
    }
    r(rng) {
        return this._items[Math.floor(rng() * this.size())];
    }
    // Method to add an item to the list
    add(item) {
        this._items.push(item);
    }
    in(item) {
        return this._items.some(a => item == a);
    }
    // Method to remove an item from the list
    removeItem(item) {
        this._items = this._items.filter(listItem => listItem !== item);
    }
    // Method to get the size of the list
    size() {
        return this._items.length;
    }
    forEach(f) {
        this._items.forEach(f);
    }
    copy() {
        return new AssemblyComponentList(this.index, [...this._items]);
    }
    pop() {
        return this._items.pop();
    }
    shuffled(rng) {
        let items_shuffled = [...this._items];
        shuffleArray(items_shuffled, rng);
        return new AssemblyComponentList(this.index, items_shuffled);
    }
    get(idx) {
        console.assert(idx >= 0);
        console.assert(idx < this.size());
        return this._items[idx];
    }
    options() {
        // todo: pass by value!
        return [...this._items];
    }
}
export class PiecewiseStage extends AssemblyComponentList {
    colormap;
    all_pieces;
    constructor(idx) {
        super(idx);
        this.colormap = new Map();
        this.all_pieces = new AssemblyComponentList(idx);
    }
    /**
     * Adds a polycube to the assembly component list
     * @param p
     */
    add(p) {
        console.assert(this !== void 0);
        // for each move in the polycube
        p.forEachMove((_, m) => {
            // find move cell
            let move = p.getMoveCell(strToVec(m));
            // for each face in the move
            move.forEachFace(f => {
                // if the face is not null (is a move face)
                if (f) {
                    // if the face color isn't in the color map
                    if (!(-f.color in this.colormap)) {
                        // make a new list to contain polycubes with this face color move
                        this.colormap[-f.color] = new AssemblyComponentList(this.index);
                    }
                    // add polycube to list to contain polycubes with this face color move
                    this.colormap[-f.color].add(p);
                }
            });
        });
        // add piece to overall list
        this.all_pieces.add(p);
    }
    contains(p) {
        return this.all_pieces.in(p);
    }
    hasColor(c) {
        return c in this.colormap;
    }
    addAll(pieces) {
        // todo: optimzie?
        pieces.forEach(p => this.add(p));
    }
}
/**
 * I've definately written this function before
 * @param cube_types
 */
export function cubeSetColorMap(cube_types) {
    let cube_map = new Map();
    cube_types.forEach(ct => {
        ct.patches.forEach(p => {
            if (!p.color) {
                return;
            }
            if (!cube_map.has(p.color)) {
                cube_map.set(p.color, []);
            }
            cube_map.get(p.color).push(ct);
        });
    });
    return cube_map;
}
//# sourceMappingURL=assembly.js.map