let CPM = require("Artistoo/build/artistoo-cjs")
const fs = require('fs')

/*  ----------------------------------
    CONFIGURATION SETTINGS
    ----------------------------------
*/
let config = {

    // Grid settings
    ndim: 2,
    field_size: [400, 400],
    CHEMOKINE_RES: 5,

    // Boolean indicating whether gathered food respawns
    RESPAWN_FOOD : true,

    // CPM parameters and configuration
    conf: {
        // Basic CPM parameters
<<<<<<< HEAD
        torus: [false, false],                        // Should the grid have linked borders?
        //seed: 1,                            // Seed for random number generation.
=======
        torus: [false, false],              // Should the grid have linked borders?
        seed: 1,                            // Seed for random number generation.
>>>>>>> eed185822e20162785e6173487092fcadc28ecbc
        T: 10,                              // CPM temperature
        D: 0.1,                             // Diffusion parameter
        SECR: 3,                            // Chemokine secrection rate
        DECAY: 0.999,                       // Chemokine decay
        LAMBDA_CH: [0, 0, 500],             // Importance of chemokine for each cell kind 

        // Constraint parameters. 
        // Mostly these have the format of an array in which each element specifies the
        // parameter value for one of the cellkinds on the grid.
        // First value is always cellkind 0 (the background) and is often not used.

        CONNECTIVITY: [false, false, true], // ConnectivyConstraint boolean parameter

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
        IMGFRAMERATE: 10,                    // If so, do this every <IMGFRAMERATE> MCS.
        SAVEPATH: "output/img/ForagingModel",    // ... And save the image in this folder.
        EXPNAME: "ForagingModel",                    // Used for the filename of output images.

        // Output stats etc
        STATSOUT: { browser: false, node: true }, // Should stats be computed?
        LOGRATE: 10,                         // Output stats every <LOGRATE> MCS.
        DEBUG: false,
        FINAL_OUTPUT: false
    }
}
/*  ---------------------------------- */
try {
	const fileConfig = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'))
	config["conf"]["MAX_ACT"] = Object.values(fileConfig["conf"]["MAX_ACT"])

	config["conf"]["LAMBDA_ACT"] = Object.values(fileConfig["conf"]["LAMBDA_ACT"])
	config["conf"]["P"] = Object.values(fileConfig["conf"]["P"])
	config["conf"]["LAMBDA_P"] = Object.values(fileConfig["conf"]["LAMBDA_P"])
	config["conf"]["V"] = Object.values(fileConfig["conf"]["V"])
	config["conf"]["LAMBDA_V"] = Object.values(fileConfig["conf"]["LAMBDA_V"])
    config["conf"]["LAMBDA_CH"] = Object.values(fileConfig["conf"]["LAMBDA_CH"])
    
  } catch (err) {
	console.error(err)
}

try{
    config["simsettings"]["SAVEIMG"] = process.argv[3]
}
catch (err){
    // No value, leave default
}
let sim, gm
let livelihood, maxLivelihood, foodIncrement, livelihoodDecay, startX, startY, endX, endY, distanceToFood

class GatheredFood {
    constructor(centroid, currentMCS) {
        this.centroid = centroid.map(function(e) {
            return Math.round(e);
        })
        this.respawnTime = currentMCS + clip(Math.floor(Math.random() * 500), 300, 500)
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
        postMCSListener: postMCSListener,
<<<<<<< HEAD
        logStats: logStats
=======
        drawCanvas: drawCanvas
>>>>>>> eed185822e20162785e6173487092fcadc28ecbc
    }

    // Foraging parameters
    livelihood = maxLivelihood = 1000
    foodIncrement = 200
    livelihoodDecay = -0.5

    eatenFood = []
    // Respawn food at random location or original?
    respawnFoodAtRandom = false

    mainCellKind = 2
    foodCellKind = 1

    sim = new CPM.Simulation(config, custommethods)
    sim.g = new CPM.Grid2D([sim.C.extents[0] / config.CHEMOKINE_RES, sim.C.extents[1] / config.CHEMOKINE_RES], config.torus, "Float32")
    sim.gi = new CPM.CoarseGrid(sim.g, config.CHEMOKINE_RES)

    sim.C.add(new CPM.ConnectivityConstraint({
        CONNECTED: config.conf.CONNECTIVITY
    }))

    sim.C.add(new CPM.ChemotaxisConstraint({
        CH_FIELD: sim.gi,
        LAMBDA_CH: config.conf.LAMBDA_CH
    }
    ))

    gm = new CPM.GridManipulator(sim.C)
    for (let cid of sim.C.cellIDs()) {
        if (sim.C.cellKind(cid) === mainCellKind) {
            let centroids = sim.C.getStat(CPM.CentroidsWithTorusCorrection)
            startX = centroids[cid][0]
            startY = centroids[cid][1]
        }
    }
}

function postMCSListener() {
    // Decay life every Monte Carlo step
    mutateLivelihood(livelihoodDecay)

    if (killCels()) {
        
        // Cell died, return from simulation


        config.simsettings.RUNTIME = -1
    }

    chemotaxisMCS()

    if (config.RESPAWN_FOOD) {
        findFoodToRespawn()
    }
}
function drawCanvas() {
    // Add the canvas if required
    if (!this.helpClasses["canvas"]) { this.addCanvas() }
    this.Cim.clear( "FFFFFF" )
    this.Cim.drawField(this.gi)

    let cellcolor=( this.conf["CELLCOLOR"] || [] ), actcolor=this.conf["ACTCOLOR"], 
    nrcells=this.conf["NRCELLS"], cellkind, cellborders = this.conf["SHOWBORDERS"]
    for( cellkind = 0; cellkind < nrcells.length; cellkind ++ ){

        // draw the cells of each kind in the right color
        if( cellcolor[ cellkind ] !== -1 ){
            this.Cim.drawCells( cellkind+1, cellcolor[cellkind] )
        }
        
        // Draw borders if required
        if(  this.conf.hasOwnProperty("SHOWBORDERS") && cellborders[ cellkind  ] ){
            let bordercol = "000000"
            if( this.conf.hasOwnProperty("BORDERCOL") ){
                bordercol = this.conf["BORDERCOL"][cellkind] || "000000"
            }
            this.Cim.drawCellBorders( cellkind+1, bordercol )
        }
        
        // if there is an activity constraint, draw activity values depending on color.
        if( this.C.conf["LAMBDA_ACT"] !== undefined && this.C.conf["LAMBDA_ACT"][ cellkind + 1 ] > 0 ){ //this.constraints.hasOwnProperty( "ActivityConstraint" ) ){
            let colorAct
            if( typeof actcolor !== "undefined" ){
                colorAct = actcolor[ cellkind ] || false
            } else {
                colorAct = false
            }
            if( ( colorAct ) ){
                this.Cim.drawActivityValues( cellkind + 1 )//, this.constraints["ActivityConstraint"] )
            }			
        }
        

    }
}

<<<<<<< HEAD
function logStats() {
		
    // compute centroids for all cells
    let allcentroids; 
    let torus = false;
    for( let d = 0; d < this.C.grid.ndim; d++ ){
        if( this.C.grid.torus[d] ){
            torus = true;
        }
    }
    if( torus ){
        allcentroids = this.C.getStat( CPM.CentroidsWithTorusCorrection );
    } else {
        allcentroids = this.C.getStat( CPM.Centroids );
    }
    for( let cid of this.C.cellIDs() ){
        if (this.C.cellKind(cid) !== mainCellKind) {
            continue;
        }
    
        let thecentroid = allcentroids[cid];
        
        // eslint-disable-next-line no-console
        console.log( this.time + ", " + thecentroid.join(", ") );
        
    }

}

function chemotaxisMCS(context) {
=======
function chemotaxisMCS() {
>>>>>>> eed185822e20162785e6173487092fcadc28ecbc
    // TODO: this crashes after the food cells are reseeded in findFoodToRespawn() (probably something with mixed grids?)
    let centroids = sim.C.getStat(CPM.Centroids)
    for (let cid in centroids) {
        if (sim.C.cellKind(cid) === foodCellKind) {
            let c = [Math.floor(centroids[cid][0] / config.CHEMOKINE_RES), Math.floor(centroids[cid][1] / config.CHEMOKINE_RES)]
            sim.g.setpix(c, sim.C.conf["SECR"])
        }
    }

    for (let i = 1; i <= config.CHEMOKINE_RES; i++) {
        sim.g.diffusion(sim.C.conf["D"])
    }
    sim.g.multiplyBy(sim.C.conf["DECAY"])
}

function findFoodToRespawn() {
    for (var i = 0; i < eatenFood.length; i++) {
        if (eatenFood[i].getRespawnTime() === sim.time) {
            if (config.simsettings.DEBUG) {
                console.log("Respawning food cell!")
                console.log(eatenFood)
            }
            if (respawnFoodAtRandom) {
                gm.seedCell(foodCellKind)
            } else {
                // cellkind, position [a, b]
                gm.seedCellAt(foodCellKind, eatenFood[i].getCentroid())
            }
            // Remove food from array
            eatenFood.splice(i, 1)
            i--
            if (config.simsettings.DEBUG) {
                console.log("Food respawned!")
                console.log(eatenFood)
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
                    if (config.simsettings.DEBUG) { console.log(`Food found, increasing livelihood by ${foodIncrement}`) }
                    mutateLivelihood(foodIncrement)

                    // Remember location of removed food before killing the food cell
                    // Get centroid of the food cell eaten
                    let centroid = sim.C.getStat(CPM.Centroids)[cid]
                    // Memorise
                    eatenFood.push(new GatheredFood(centroid, sim.time))
                    // Kill
                    gm.killCell(cid)
                }
            }
        } else if (sim.C.cellKind(cid) === mainCellKind && livelihood <= 0) {
            // TODO: check neighbors of main cell for food instead of vice versa
            if (config.simsettings.DEBUG) { console.log("Killed the cell due to starvation!") }
            if(typeof endX == "undefined" || typeof endY == "undefined"){
                // Undefined final position because cell is alive, calculate now
                let centroids = sim.C.getStat(CPM.Centroids)
                endX = centroids[cid][0]
                endY = centroids[cid][1]
            }
            distanceToFood = calc_distance_to_nearest_food()
            gm.killCell(cid)
            return true
        }
    }
    return false
}

initialize()

console.log("t, x, y");
sim.run()

function euclidean_distance(x1, x2, y1, y2){
    return Math.sqrt(Math.pow(x1-x2, 2) + Math.pow(y1-y2, 2))
}

function calc_distance_to_nearest_food(){
    let _distance_to_food = 100000
    let _cellX, _cellY
    let centroids = sim.C.getStat(CPM.Centroids)
    for (let cid of sim.C.cellIDs()) {
        if (sim.C.cellKind(cid) === mainCellKind) {
            _cellX = centroids[cid][0]
            _cellY = centroids[cid][1]
        }
    }
    for (let cid of sim.C.cellIDs()) {
        if (sim.C.cellKind(cid) === foodCellKind) {
            let distance = euclidean_distance(_cellX, centroids[cid][0], _cellY, centroids[cid][1])
            _distance_to_food = Math.min(_distance_to_food, distance)
        }
    }
    return _distance_to_food
}

if (config.simsettings.FINAL_OUTPUT) {
    if(typeof distanceToFood == "undefined"){
        distanceToFood = calc_distance_to_nearest_food()
    }
    let distance_traveled = Math.sqrt(Math.pow(startX-endX, 2) + Math.pow(startY-endY, 2))
    console.log(sim.time+","+livelihood+","+(-distanceToFood))
}
