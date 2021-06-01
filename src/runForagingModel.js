let CPM = require("Artistoo/build/artistoo-cjs")


/*  ----------------------------------
    CONFIGURATION SETTINGS
    ----------------------------------
*/
let config = {

    // Grid settings
    ndim: 2,
    field_size: [400, 400],
    CHEMOKINE_RES : 5,

    // CPM parameters and configuration
    conf: {
        // Basic CPM parameters
        torus: [false, false],                        // Should the grid have linked borders?
        seed: 1,                            // Seed for random number generation.
        T: 10,                              // CPM temperature
        D: 0.1,                             // Diffusion parameter
        SECR: 3,                            // Chemokine secrection rate
        DECAY: 0.999,                        // Chemokine decay
        LAMBDA_CH: [0, 0, 500],            // Importance of chemokine for each cell kind 

        // Constraint parameters. 
        // Mostly these have the format of an array in which each element specifies the
        // parameter value for one of the cellkinds on the grid.
        // First value is always cellkind 0 (the background) and is often not used.

        LAMBDA_CONNECTIVITY : [0, 0, 100], // (Soft(Local))ConnectivyConstraint lambda parameter for the strength of the penalty of breaking connectiviy

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
        RUNTIME: 10000,
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
        SAVEPATH: "output/img/ForagingModel",    // ... And save the image in this folder.
        EXPNAME: "ForagingModel",                    // Used for the filename of output images.

        // Output stats etc
        STATSOUT: { browser: false, node: false }, // Should stats be computed?
        LOGRATE: 10,                         // Output stats every <LOGRATE> MCS.
        DEBUG : false,
        FINAL_OUTPUT : true
    }
}
/*  ---------------------------------- */


let sim, gm
let livelihood, maxLivelihood, foodIncrement, livelihoodDecay

class GatheredFood {
    constructor(centroid, currentMCS) {
        this.centroid = centroid
        this.respawnTime = currentMCS + clip(Math.floor(Math.random()*100), 20, 100)
    }

    getRespawnTime() {
        return this.respawnTime
    }

    getCentroid() {
        return this.centroid
    }
}

function initialize() {
    let custommethods = {
        postMCSListener : postMCSListener
    }

    // Foraging parameters
    livelihood = maxLivelihood = 1000
    foodIncrement = 200
    livelihoodDecay = -0.5

    eatenFood = []
    // Respawn food at random location or original?
    respawnFoodAtRandom = true

    mainCellKind = 2
    foodCellKind = 1

    sim = new CPM.Simulation(config, custommethods)
    sim.g = new CPM.Grid2D([sim.C.extents[0] / config.CHEMOKINE_RES, sim.C.extents[1] / config.CHEMOKINE_RES], config.torus, "Float32")
    sim.gi = new CPM.CoarseGrid(sim.g, config.CHEMOKINE_RES)

    sim.C.add(new CPM.ConnectivityConstraint({
        LAMBDA_CONNECTIVITY : config.conf.LAMBDA_CONNECTIVITY
    }))

    sim.C.add(new CPM.ChemotaxisConstraint({
        CH_FIELD: sim.gi,
        LAMBDA_CH: config.conf.LAMBDA_CH
    }
    ))

    gm = new CPM.GridManipulator(sim.C)
}

function postMCSListener() {
    // Decay life every Monte Carlo step
    mutateLivelihood(livelihoodDecay)
    
    if (killCels()) {
        // Cell died, return from simulation
        config.simsettings.RUNTIME = -1
    }
    
    // findFoodToRespawn()

    chemotaxisMCS(this)
}

function chemotaxisMCS(context) {
    let centroids = context.C.getStat(CPM.CentroidsWithTorusCorrection)
    for (let cid in centroids) {
        if (context.C.cellKind(cid) === foodCellKind) {
            let c = [Math.floor(centroids[cid][0] / config.CHEMOKINE_RES), Math.floor(centroids[cid][1] / config.CHEMOKINE_RES)]
            context.g.setpix(c, context.C.conf["SECR"])
        }
    }

    for (let i = 1; i <= config.CHEMOKINE_RES; i++) {
        context.g.diffusion(context.C.conf["D"])
    }
    context.g.multiplyBy(context.C.conf["DECAY"])
}

function findFoodToRespawn() {
    for (let food of eatenFood) {
        if (food.getRespawnTime() == sim.time) {
            if (config.simsettings.DEBUG) {
                console.log("Respawning food cell!")
            }
            // TODO: seedCell(At) crashes on the length of a centroid..
            if (respawnFoodAtRandom) {
                gm.seedCell(foodCellKind)
            } else {
                // cellkind, position [a, b]
                console.log(food.getCentroid().length   )
                gm.seedCellAt(foodCellKind, food.getCentroid())
            }
        }
    }
}

function clip(x, lb, ub) {
    return Math.max(lb, Math.min(ub, x))
}

function mutateLivelihood(amount) {
    livelihood = clip(livelihood + amount, 0, maxLivelihood)
}

// Returns true if the main cell died
function killCels() {
    for (let cid of sim.C.cellIDs()) {
        if (sim.C.cellKind(cid) === foodCellKind) {
            let conncomp = sim.C.getStat(CPM.CellNeighborList)
            for (let ocid in conncomp[cid]) {
                if (sim.C.cellKind(ocid) === mainCellKind) {
                    if(config.simsettings.DEBUG) { console.log(`Food found, increasing livelihood by ${foodIncrement}`) }
                    mutateLivelihood(foodIncrement)
                    
                    // Remember location of removed food before killing the food cell
                    // Get centroid of the food cell eaten
                    let centroid = sim.C.getStat(CPM.CentroidsWithTorusCorrection)[cid]
                    // Memorise
                    eatenFood.push(new GatheredFood(centroid, sim.time))
                    // Kill
                    gm.killCell(cid)
                }
            }
        } else if (sim.C.cellKind(cid) === mainCellKind && livelihood <= 0) {
            // TODO: check neighbors of main cell for food instead of vice versa
            if (config.simsettings.DEBUG) { console.log("Killed the cell due to starvation!") }
            gm.killCell(cid)
            return true
        }
    }
    return false
}

initialize()

sim.run()

if (config.simsettings.FINAL_OUTPUT) {
    console.log(sim.time)
}