let CPM = require("Artistoo/build/artistoo-cjs")


/*  ----------------------------------
    CONFIGURATION SETTINGS
    ----------------------------------
*/
let config = {

    // Grid settings
    ndim: 2,
    field_size: [200, 200],
    chemokine_res: 5,

    // CPM parameters and configuration
    conf: {
        // Basic CPM parameters
        torus: [true, true],                        // Should the grid have linked borders?
        seed: 1,                            // Seed for random number generation.
        T: 10,                              // CPM temperature
        D: 0.1,                             // Diffusion parameter
        SECR: 5,                            // Chemokine secrection rate
        DECAY: 0.99,                        // Chemokine decay
        LAMBDA_CH: [0, 0, 500],            // Importance of chemokine for each cell kind 

        // Constraint parameters. 
        // Mostly these have the format of an array in which each element specifies the
        // parameter value for one of the cellkinds on the grid.
        // First value is always cellkind 0 (the background) and is often not used.

        // Adhesion parameters:
        J: [[0, 100, 10], [100, 10, -1], [10, -1, 0]],

        // VolumeConstraint parameters
        LAMBDA_V: [0, 1000, 5],                     // VolumeConstraint importance per cellkind
        V: [0, 10, 500],                            // Target volume of each cellkind

        // PerimeterConstraint parameters
        LAMBDA_P: [0, 1, 2],                        // PerimeterConstraint importance per cellkind
        P: [0, 5, 260],                             // Target perimeter of each cellkind

        // ActivityConstraint parameters
        LAMBDA_ACT: [0, 0, 300],                // ActivityConstraint importance per cellkind
        MAX_ACT: [0, 0, 30],                    // Activity memory duration per cellkind
        ACT_MEAN: "geometric"                   // Is neighborhood activity computed as a
        // "geometric" or "arithmetic" mean?

    },

    // Simulation setup and configuration
    simsettings: {

        // Cells on the grid
        NRCELLS: [10, 1],                       // Number of cells to seed for all
        // non-background cellkinds.
        // Runtime etc
        BURNIN: 500,
        RUNTIME: 1000,
        RUNTIME_BROWSER: "Inf",

        // Visualization
        CANVASCOLOR: "eaecef",
        CELLCOLOR: ["00ff00", "000000"],
        ACTCOLOR: [false, true],                    // Should pixel activity values be displayed?
        SHOWBORDERS: [true, false],             // Should cellborders be displayed?
        zoom: 2,                            // zoom in on canvas with this factor.

        // Output images
        SAVEIMG: false,                 // Should a png image of the grid be saved
        // during the simulation?
        IMGFRAMERATE: 1,                    // If so, do this every <IMGFRAMERATE> MCS.
        SAVEPATH: "output/img/ActModel",    // ... And save the image in this folder.
        EXPNAME: "ActModel",                    // Used for the filename of output images.

        // Output stats etc
        STATSOUT: { browser: false, node: true }, // Should stats be computed?
        LOGRATE: 10                         // Output stats every <LOGRATE> MCS.

    }
}
/*  ---------------------------------- */


let sim, meter, gm


function initialize() {
    let custommethods = {
        postMCSListener: postMCSListener,
        drawCanvas: drawCanvas
    }

    sim = new CPM.Simulation(config, custommethods)
    sim.g = new CPM.Grid2D([sim.C.extents[0] / config.chemokine_res, sim.C.extents[1] / config.chemokine_res], config.torus, "Float32")
    sim.gi = new CPM.CoarseGrid(sim.g, config.chemokine_res)

    sim.C.add(new CPM.ChemotaxisConstraint({
        CH_FIELD: sim.gi,
        LAMBDA_CH: config.conf.LAMBDA_CH
    }
    ))

    gm = new CPM.GridManipulator(sim.C)
}

function postMCSListener() {
    killCels()

    let centroids = this.C.getStat(CPM.CentroidsWithTorusCorrection)
    for (let cid in centroids) {
        if (this.C.cellKind(cid) === 1) {
            let c = [Math.floor(centroids[cid][0] / config.chemokine_res), Math.floor(centroids[cid][1] / config.chemokine_res)]
            this.g.setpix(c, this.C.conf["SECR"])
        }
    }

    for (let i = 1; i <= config.chemokine_res; i++) {
        this.g.diffusion(this.C.conf["D"])
    }
    this.g.multiplyBy(this.C.conf["DECAY"])
}

function killCels() {
    let conncomp = sim.C.getStat(CPM.CellNeighborList)
    for (let cid of sim.C.cellIDs()) {
        if (sim.C.cellKind(cid) === 1) {
            for (let ocid in conncomp[cid]) {
                if (sim.C.cellKind(ocid) === 2) {
                    gm.killCell(cid)
                }
            }
        }
    }
}

function drawCanvas() {
    // Add the canvas if required
    if (!this.helpClasses["canvas"]) { this.addCanvas() }
    this.Cim.drawField(this.gi)
    this.Cim.drawCells(1, "00ff00")
    this.Cim.drawActivityValues(2)
    this.Cim.drawCellBorders(2, "000000")

}

initialize()

sim.run()
