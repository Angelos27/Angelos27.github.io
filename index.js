var U_FIELD = 0;
var V_FIELD = 1;
var S_FIELD = 2;



// these are used only to call setObstacle to adjuct the angle via the slider

var globalX = 0.6;
var globalY = 0.5;

//var cnt = 0;

function cX(x) {
    return x * cScale;
}

function cY(y) {
    return canvas.height - y * cScale;
}

// ----------------- start of simulator ------------------------------

class Fluid {
    constructor(density, numX, numY, h) {
        // Initializes the fluid properties and simulation grid.

        this.density = density; 
        // The density of the fluid, which affects how it behaves under  gravity).
        this.numX = numX + 2; 
        // Number of grid cells in the x-direction, including 2 extra for boundary cells.
        this.numY = numY + 2; 
        // Number of grid cells in the y-direction, including 2 extra for boundary cells.
        this.numCells = this.numX * this.numY; 
        // Total number of cells in the simulation grid (including boundaries).
        this.cellSize = h; 
        // Size of each grid cell (spatial resolution of the grid).

        // Fluid properties stored as 1D arrays for efficient computation.

        this.horizontalVelocity = new Float32Array(this.numCells); 
        // Horizontal velocity component (u) at each grid cell.
        this.verticalVelocity = new Float32Array(this.numCells); 
        // Vertical velocity component (v) at each grid cell.
        this.newU = new Float32Array(this.numCells); 
        // Temporary array to store updated horizontal velocity during advection.
        this.newV = new Float32Array(this.numCells); 
        // Temporary array to store updated vertical velocity during advection.
        this.p = new Float32Array(this.numCells); 
        // Pressure values at each grid cell (used to correct velocities).
        this.cellState = new Float32Array(this.numCells); 
        // Solid mask: 0 for solid cells (walls), 1 for fluid cells.
        this.m = new Float32Array(this.numCells); 
        // Scalar field values (e.g., smoke density).
        this.newM = new Float32Array(this.numCells); 
        // Temporary array for storing updated scalar field values.
        this.m.fill(1.0); 
        // Initialize the scalar field (e.g., smoke) with a uniform value of 1.0.

        this.totalDivergence = 0;// Total Divergence
        this.divFreq = 0; //Total divergence additions
    }

    integrate(dt, gravity) {
        // Applies external forces (like gravity) to the velocity field.

        var n = this.numY; 
        // `n` is the stride (number of rows) in the grid for easy indexing.

        for (var i = 1; i < this.numX; i++) {
            // Iterate over interior cells in the x-direction (excluding boundaries).

            for (var j = 1; j < this.numY - 1; j++) {
                // Iterate over interior cells in the y-direction (excluding boundaries).

                if (this.cellState[i * n + j] != 0.0 && this.cellState[i * n + j - 1] != 0.0) {
                    // Apply gravity only if the current cell and the one below it are fluid cells.
                    this.verticalVelocity[i * n + j] += gravity * dt; 
                    // Add gravity to the vertical velocity, scaled by the time step.
                }
            }
        }
    }

    solveIncompressibility(numIters, dt) {
        // Ensures that the fluid velocity field is incompressible (divergence-free).
        var n = this.numY; 
        // `n` is the stride for accessing grid cells in a 1D array.
        var pressureCoef = this.density * this.cellSize / dt; 
        // Pressure coefficient based on fluid density, grid spacing, and time step.

        for (var iter = 0; iter < numIters; iter++) {
            // Perform the specified number of iterations to solve for pressure and correct velocities.
            for (var i = 1; i < this.numX - 1; i++) {
                for (var j = 1; j < this.numY - 1; j++) {
                    if (this.cellState[i * n + j] == 0.0) continue; 
                    // Skip solid cells, as they do not contribute to the fluid velocity.
                    
                    var sx0 = this.cellState[(i - 1) * n + j]; 
                    var sx1 = this.cellState[(i + 1) * n + j];
                    var sy0 = this.cellState[i * n + j - 1]; 
                    var sy1 = this.cellState[i * n + j + 1]; 
                    // Gather the solid mask values of neighboring cells (left, right, bottom, top).
                    
                    var s = sx0 + sx1 + sy0 + sy1; 
                    // Sum of neighboring solid mask values; determines fluid connectivity.
                    
                    if (s == 0.0) continue; 
                    // Skip if surrounded entirely by solids.
                    
                    var divergence = this.horizontalVelocity[(i + 1) * n + j] - this.horizontalVelocity[i * n + j] + 
                              this.verticalVelocity[i * n + j + 1] - this.verticalVelocity[i * n + j]; 
                    // Compute the velocity divergence for the current cell.

                    this.totalDivergence += divergence; // Add the divergence to total
                    this.divFreq += 1; // Frequecy count

                    var pressureAdjust = -divergence / s;  // Calculate pressure adjustment to reduce divergence.

                    pressureAdjust *= scene.overRelaxation; 
                    // Apply over-relaxation to accelerate convergence (a numerical optimization).
                    this.p[i * n + j] += pressureCoef * pressureAdjust; 
                    // Update pressure value for the current cell.
                    this.horizontalVelocity[i * n + j] -= sx0 * pressureAdjust; 
                    this.horizontalVelocity[(i + 1) * n + j] += sx1 * pressureAdjust; 
                    this.verticalVelocity[i * n + j] -= sy0 * pressureAdjust; 
                    this.verticalVelocity[i * n + j + 1] += sy1 * pressureAdjust; 
                    // Correct the velocity components based on the calculated pressure.
                }
            }
        }
    }

    extrapolate() {
        // Copies fluid velocity values to boundary cells to maintain consistent behavior.

        var n = this.numY;

        for (var i = 0; i < this.numX; i++) {
            this.horizontalVelocity[i * n + 0] = this.horizontalVelocity[i * n + 1]; 
            // Copy horizontal velocity from the second row to the first boundary row.
            this.horizontalVelocity[i * n + this.numY - 1] = this.horizontalVelocity[i * n + this.numY - 2]; 
            // Copy horizontal velocity from the second-to-last row to the last boundary row.
        }

        for (var j = 0; j < this.numY; j++) {
            this.verticalVelocity[0 * n + j] = this.verticalVelocity[1 * n + j]; 
            // Copy vertical velocity from the second column to the first boundary column.
            this.verticalVelocity[(this.numX - 1) * n + j] = this.verticalVelocity[(this.numX - 2) * n + j]; 
            // Copy vertical velocity from the second-to-last column to the last boundary column.
        }
    }

    fieldSample(x, y, field) {
        // Samples a value from a given field (u, v, or scalar) at position (x, y) using bilinear interpolation.

        var n = this.numY; 
        var h = this.cellSize;
        var h1 = 1.0 / h; // The reciprocal of the cell size (used to normalize positions to grid space).
        var h2 = 0.5 * h; // Half the cell size (used for staggered grids or offsets).
         
        x = Math.max(Math.min(x, this.numX * h), h); 
        y = Math.max(Math.min(y, this.numY * h), h); 
        // Clamp x and y to stay within grid boundaries.
        // This prevents accessing invalid grid cells, which could cause errors or undefined behavior.

        var dx = 0.0, dy = 0.0, f;

        switch (field) {
            case U_FIELD: f = this.horizontalVelocity; dy = h2; break; 
            case V_FIELD: f = this.verticalVelocity; dx = h2; break; 
            case S_FIELD: f = this.m; dx = h2; dy = h2; break; 
            // Choosing the correct field and offsets based on the type of field requested.
        }
        // Compute the indices of the four surrounding grid cells in the x-direction.
        var x0 = Math.min(Math.floor((x - dx) * h1), this.numX - 1);
        // `tx` is the fractional distance between `x0` and the sampling position.
        var tx = ((x - dx) - x0 * h) * h1;
        // `x1` is the index of the grid cell to the right of the sampling position.
        var x1 = Math.min(x0 + 1, this.numX - 1);
        // `x1` is the index of the grid cell to the right of the sampling position.

        // Compute the indices of the four surrounding grid cells in the y-direction.
        var y0 = Math.min(Math.floor((y - dy) * h1), this.numY - 1); 
        // `y0` is the index of the grid cell below the sampling position.
        var ty = ((y - dy) - y0 * h) * h1;
        // `ty` is the fractional distance between `y0` and the sampling position.
        var y1 = Math.min(y0 + 1, this.numY - 1);
        // `y1` is the index of the grid cell above the sampling position.

        // Compute the interpolation weights.
        var sx = 1.0 - tx, sy = 1.0 - ty;

        // Perform bilinear interpolation.
        // Interpolates between the values of the four nearest grid cells:
        var val = sx * sy * f[x0 * n + y0] +  // Bottom-left corner contribution.
              tx * sy * f[x1 * n + y0] +  // Bottom-right corner contribution.
              tx * ty * f[x1 * n + y1] +  // Top-right corner contribution.
              sx * ty * f[x0 * n + y1];   // Top-left corner contribution.
        // Perform bilinear interpolation using the values from the four nearest grid cells.

        return val; 
        // Return the interpolated value.
    }

    avgU(i, j) {
        // Computes the average horizontal velocity (u) around a cell (i, j).

        var n = this.numY;
        return (this.horizontalVelocity[i * n + j - 1] + this.horizontalVelocity[i * n + j] +
                this.horizontalVelocity[(i + 1) * n + j - 1] + this.horizontalVelocity[(i + 1) * n + j]) * 0.25;
    }

    avgV(i, j) {
        // Computes the average vertical velocity (v) around a cell (i, j).

        var n = this.numY;
        return (this.verticalVelocity[(i - 1) * n + j] + this.verticalVelocity[i * n + j] +
                this.verticalVelocity[(i - 1) * n + j + 1] + this.verticalVelocity[i * n + j + 1]) * 0.25;
    }

    advectVel(dt) {
            // Updates the velocity field by advecting it using the existing velocity field.
            // Copy the current velocity fields into temporary storage
            this.newU.set(this.horizontalVelocity);
            this.newV.set(this.verticalVelocity);
            var n = this.numY, h = this.cellSize, h2 = 0.5 * h;
        
            // Loop through the grid, excluding boundary cells
            for (var i = 1; i < this.numX; i++) {
                for (var j = 1; j < this.numY; j++) {
        
                    // Check if the current cell and its left neighbor are fluid (not boundary)
                    if (this.cellState[i * n + j] != 0.0 && this.cellState[(i - 1) * n + j] != 0.0 && j < this.numY - 1) {
                        // Compute the physical position of the cell center
                        var x = i * h, y = j * h + h2;
        
                        // Retrieve velocity components at the current cell
                        var u = this.horizontalVelocity[i * n + j], v = this.avgV(i, j);
        
                        // Backtrace the position by moving against the velocity field
                        x -= dt * u;
                        y -= dt * v;
        
                        // Sample the velocity field at the backtraced position
                        this.newU[i * n + j] = this.fieldSample(x, y, U_FIELD);
                    }
                    // Check if the current cell and its below neighbor are fluid (not boundary)
                    if (this.cellState[i * n + j] != 0.0 && this.cellState[i * n + j - 1] != 0.0 && i < this.numX - 1) {
                        // Compute the physical position of the cell center
                        var x = i * h + h2, y = j * h;
        
                        // Retrieve velocity components at the current cell
                        var u = this.avgU(i, j), v = this.verticalVelocity[i * n + j];
                        // Backtrace the position by moving against the velocity field
                        x -= dt * u;
                        y -= dt * v;
                        // Sample the velocity field at the backtraced position
                        this.newV[i * n + j] = this.fieldSample(x, y, V_FIELD);
                    }
                }
            }
            // Update the velocity fields with the newly computed values
            this.horizontalVelocity.set(this.newU);
            this.verticalVelocity.set(this.newV);
        }    
    
    advectSmoke(dt) {
        // Updates the scalar field (e.g., smoke density) by advecting it along the velocity field.
        // This simulates the movement of smoke particles by following the flow.
    
        // Copy the current smoke field into a temporary buffer (newM) for updated values.
        this.newM.set(this.m);
    
        // Retrieve grid parameters: n (number of cells in y-direction), cell size h, and half cell size h2.
        var n = this.numY, h = this.cellSize, h2 = 0.5 * h;
    
        // Loop over the grid cells, excluding boundary cells to avoid sampling issues.
        for (var i = 1; i < this.numX - 1; i++) {
            for (var j = 1; j < this.numY - 1; j++) {
                // Process only cells that are marked as fluid (non-zero cell state).
                if (this.cellState[i * n + j] != 0.0) {
                    // Compute the average horizontal velocity at the cell face:
                    // Average between the current cell and the right neighbor.
                    var u = (this.horizontalVelocity[i * n + j] + this.horizontalVelocity[(i + 1) * n + j]) * 0.5;
    
                    // Compute the average vertical velocity at the cell face:
                    // Average between the current cell and the top neighbor.
                    var v = (this.verticalVelocity[i * n + j] + this.verticalVelocity[i * n + j + 1]) * 0.5;
    
                    // Determine the cell's physical center position:
                    // i * h + h2 gives the x-coordinate; j * h + h2 gives the y-coordinate.
                    // Then, backtrace the position along the velocity field (semi-Lagrangian method)
                    // by subtracting the velocity components scaled by the time step dt.
                    var x = i * h + h2 - dt * u;
                    var y = j * h + h2 - dt * v;
    
                    // Sample the smoke field (S_FIELD) at the backtraced position (x, y)
                    // using bilinear interpolation for smoothness.
                    this.newM[i * n + j] = this.fieldSample(x, y, S_FIELD);
                }
            }
        }
        // Update the main smoke field with the advected values stored in newM.
        this.m.set(this.newM);
    }
    

    simulate(dt, gravity, numIters) {
        // Executes the entire simulation process for one time step.

        this.integrate(dt, gravity); // Step 1: Apply external forces.
        this.p.fill(0.0); // Reset pressure field to zero.
        this.solveIncompressibility(numIters, dt); // Step 2: Solve for incompressible velocity field.
        this.extrapolate(); // Step 3: Update boundary values.
        this.advectVel(dt); // Step 4: Advect velocity field.
        this.advectSmoke(dt); // Step 5: Advect scalar field (e.g., smoke).
    }
}


var scene = 
{
    gravity : -9.81,
    dt : 1.0 / 120.0,
    numIters : 500,
    frameNr : 0,
    overRelaxation : 1.9,
    obstacleX : 0.0,
    obstacleY : 0.0,
    obstacleRadius: 0.15,
    paused: false,
    sceneNr: 0,
    showObstacle: false,
    showStreamlines: false,	
    showPressure: false,
    showSmoke: false,
    smokeMode: false,
    fluid: null,
    circle: true,
    square: false,
    rainDrop: false,
    aeroFoil: false,
    showAngleSlide: false,
    draw: false,
    count:0,
    totalPressure:0,
    lowestMinP : 100000000,
    highestMaxP : -100000000


};

function setupScene(sceneNr = 0) 
{
    // Set the current scene number (determines the type of simulation)
    scene.sceneNr = sceneNr;
    
    // Define the radius of the obstacle in the simulation
    scene.obstacleRadius = 0.20;
    
    // Over-relaxation factor for the numerical solver (used to speed up convergence)
    scene.overRelaxation = 1.9;

    // Set the simulation timestep (how much time advances per frame)
    scene.dt = 1.0 / 90; // 90 FPS simulation speed

    // Set the number of solver iterations (more iterations improve accuracy)
    scene.numIters = 80;

    // Retrieve the simulation resolution from the UI input field
    var res = document.getElementById("resolution").value;

    // Define the height of the simulation domain (fixed to 1.0 units)
    var domainHeight = 1.0;

    // Compute the domain width based on aspect ratio (simWidth/simHeight)
    var domainWidth = domainHeight / simHeight * simWidth;

    // Calculate the size of each computational cell based on the resolution
    var h = domainHeight / res; 

    // Determine the number of cells in the X and Y directions
    var numX = Math.floor(domainWidth / h);
    var numY = Math.floor(domainHeight / h);

    // Define the density of the fluid (in kg/m³, assuming water)
    var density = 1000.0;

    // Initialize the fluid simulation with computed grid dimensions
    f = scene.fluid = new Fluid(density, numX, numY, h);

    // Store the number of Y-direction cells for quick access
    var n = f.numY;

    if (sceneNr == 0) {   		// tank

        for (var i = 0; i < f.numX; i++) {
            for (var j = 0; j < f.numY; j++) {
                var s = 1.0;	// fluid
                if (i == 0 || i == f.numX-1 || j == 0)
                    s = 0.0;	// solid
                f.cellState[i*n + j] = s
            }
        }
        scene.gravity = -9.81;
        scene.showPressure = true;
        scene.showSmoke = false;
        scene.showStreamlines = false;
        setObstacle(0.6,0.5, true)
    }
    // If the selected scene is "Vortex Shedding" (Scene 1 or 3)
    else if (sceneNr == 1 || sceneNr == 3) { 
        
        // Retrieve the velocity slider element from the UI
        var slider = document.getElementById("speedRange");
        
        // Attach an event listener to update velocity when the slider is moved
        slider.oninput = function(){
            var output = document.getElementById("velOutout");
            output.innerHTML = slider.value;
            setupScene(1); // Reinitialize scene on slider change
        }
        
        // Set the inlet velocity based on slider input (scaled by a factor of 10)
        var inVel = slider.value / 10; 
        
        // Loop through each grid cell and define the fluid region
        for (var i = 0; i < f.numX; i++) {
            for (var j = 0; j < f.numY; j++) {
                var s = 1.0; // Default: Fluid cell
                
                // Set solid boundary conditions at the domain edges:
                if (i == 0 || j == 0|| j == f.numY - 1)
                    s = 0.0; // Mark as solid

                f.cellState[i * n + j] = s; // Assign cell type to the grid

                // Set the velocity at the leftmost column (fluid inflow)
                if (i == 1) {
                    f.horizontalVelocity[i * n + j] = inVel; // Assign horizontal velocity
                }
            }
        }

        if(scene.smokeMode){
            // Define the pipe region within the fluid domain
            var pipeH = 0.1 * f.numY; // Pipe height as a fraction of domain height
            var minJ = Math.floor(0.5 * f.numY - 0.5 * pipeH); // Lower bound of pipe
            var maxJ = Math.floor(0.5 * f.numY + 0.5 * pipeH); // Upper bound of pipe

            // Set pipe constraints in the fluid model (disables flow inside pipe)
            for (var j = minJ; j < maxJ; j++){
                f.m[j] = 0.0;
            }
        }else{
            var numPipes = 10; // Number of repeated pipes
            var pipeSpacing = Math.floor(f.numY / (numPipes + 1)); // Space between pipes
            var pipeH = 0.05 * f.numY; // Height of each pipe
        
            for (var p = 1; p <= numPipes; p++) {
                var centerJ = p * pipeSpacing; // Calculate pipe position
                var minJ = Math.max(0, Math.floor(centerJ - 0.5 * pipeH));
                var maxJ = Math.min(f.numY - 1, Math.floor(centerJ + 0.5 * pipeH));
        
               for (var j = minJ; j < maxJ; j++) {
                    f.m[j] = 0.0; // Mark pipe region
                }
            }
        }
        
        // Place an obstacle at the center of the flow
        setObstacle(0.6, 0.5, true);

        // Disable gravity for this simulation since we are modeling air/water flow
        scene.gravity = 0.0;

        // Configure visualization settings
        scene.showPressure = false; // Do not display pressure contours
        scene.showSmoke = true; // Enable smoke visualization for flow tracking
        scene.showStreamlines = false; // Disable streamlines visualization
        scene.showVelocities = false; // Do not show velocity vectors

        // If Scene 3 is selected, modify simulation properties
        if (sceneNr == 3) {
            scene.dt = 1.0 / 120.0; // Increase time resolution for more accuracy
            scene.numIters = 100; // Maintain high iteration count
            scene.showPressure = true; // Enable pressure visualization
            scene.showSmoke = false; // Disable smoke rendering
        }
    }
    // If Scene 2 is selected: "Paint" mode (passive scalar diffusion)
    else if (sceneNr == 2) { 
        setObstacle(0.6, 0.5, true);
        scene.gravity = 0.0; // No gravity in this mode
        scene.overRelaxation = 1.0; // Use normal relaxation instead of over-relaxation
        scene.showPressure = false; // No pressure visualization
        scene.showSmoke = true; // Enable smoke visualization (color mixing)
        scene.showStreamlines = false; // Disable streamlines
        scene.showVelocities = false; // Disable velocity field visualization
        scene.obstacleRadius = 0.1; // Set smaller obstacle for paint mode
    }

    // Sync UI checkbox states with the scene settings
    document.getElementById("streamButton").checked = scene.showStreamlines;
    document.getElementById("pressureButton").checked = scene.showPressure;
    document.getElementById("smokeButton").checked = scene.showSmoke;
    
}

// draw -------------------------------------------------------

function setColor(r, g, b) {
    // Set the fill and stroke color for the canvas context
    // Converts normalized RGB values (0 to 1) to 0-255 range
    c.fillStyle = `rgb(${Math.floor(255 * r)}, ${Math.floor(255 * g)}, ${Math.floor(255 * b)})`;
    c.strokeStyle = `rgb(${Math.floor(255 * r)}, ${Math.floor(255 * g)}, ${Math.floor(255 * b)})`;
}
function getSciColor(val, minVal, maxVal) {
    // Ensure the value is within the given range to prevent out-of-bounds errors
    val = Math.min(Math.max(val, minVal), maxVal - 0.0001);
    
    // Compute the difference between max and min values for normalization
    var d = maxVal - minVal;
    
    // Normalize val to a range of 0 to 1
    val = d == 0.0 ? 0.5 : (val - minVal) / d;
    
    // Define the segment width for the color scale
    var m = 0.25;
    
    // Determine which segment of the color scale the value falls into
    var num = Math.floor(val / m);
    
    // Compute the fractional position within the current segment
    var s = (val - num * m) / m;
    
    // Variables for RGB components
    var r, g, b;

    // Assign RGB values based on the segment number to create a gradient
    switch (num) {
        case 0: r = 0.0; g = s; b = 1.0; break;  // Transition from blue to cyan
        case 1: r = 0.0; g = 1.0; b = 1.0 - s; break;  // Transition from cyan to green
        case 2: r = s; g = 1.0; b = 0.0; break;  // Transition from green to yellow
        case 3: r = 1.0; g = 1.0 - s; b = 0.0; break;  // Transition from yellow to red
    }
    
    // Return the computed color as an RGBA array (alpha is set to 255 for full opacity)
    return [255 * r, 255 * g, 255 * b, 255];
}

function drawCanvas() {
    // Clear the entire canvas before rendering the new frame
    c.clearRect(0, 0, canvas.width, canvas.height);
    
    var f = scene.fluid; // Retrieve the fluid object from the scene
    var n = f.numY; // Number of grid cells in the Y direction
    var cellScale = 1.1;
    var h = f.cellSize; // Grid cell size

    // Initialize min and max pressure values for visualization
    var minP = f.p[0];
    var maxP = f.p[0];

    // Find the min and max pressure values for color mapping
    for (var i = 0; i < f.numCells; i++) {
        minP = Math.min(minP, f.p[i]);
        maxP = Math.max(maxP, f.p[i]);
    }

    //minP = Math.min(...f.p);
    //maxP = Math.max(...f.p);


    // Get image data from canvas (used for direct pixel manipulation)
    var id = c.getImageData(0, 0, canvas.width, canvas.height);
    var color = [255, 255, 255, 255]; // Default white color

    // Loop through each cell in the fluid grid
    for (var i = 0; i < f.numX; i++) {
        for (var j = 0; j < f.numY; j++) {
            if (scene.showPressure) {
                if (f.cellState[i * n + j] == 0.0) {
                    color = [255, 255, 255, 255]; // Solid obstacles are black
                } else {
                    var p = f.p[i * n + j]; // Get pressure value
                    var s = f.m[i * n + j]; // Smoke density (if enabled)
                    
                    color = getSciColor(p,0.75*minP, maxP); // Convert pressure to color scale
                    if (scene.showSmoke) {
                        // Darken the color based on smoke density
                        color[0] = Math.max(0.0, color[0] - 255 * s);
                        color[1] = Math.max(0.0, color[1] - 255 * s);
                        color[2] = Math.max(0.0, color[2] - 255 * s);
                    }
                }
            } else if (scene.showSmoke) {
                var s = f.m[i * n + j]; // Get smoke density
                color = [255 * s, 255 * s, 255 * s, 255]; // Grayscale color
                if (f.cellState[i * n + j] == 0.0) {
                    color = [0, 0, 255, 255]; // Solid obstacles displayed as blue
                } else if (scene.sceneNr == 2) {
                    color = getSciColor(s, 0.0, 1.0); // Apply scientific color scale
                }
            } else if (f.cellState[i * n + j] == 0.0) {
                color = [0, 0, 255, 255]; // Solid obstacles displayed as blue

            }else{
                color = [255, 255, 255, 255];
            }
            
            

            var x = Math.floor(cX((i-1.4)* h)); // Convert grid position to canvas coordinates
            var y = Math.floor(cY((j + 1.1) * h));
            var cx = Math.floor(cScale * cellScale * h) + 1;
            var cy = Math.floor(cScale * cellScale * h) + 1;

            // Apply color to the pixel buffer
            for (var yi = y; yi < y + cy; yi++) {
                var p = 4 * (yi * canvas.width + x);
                for (var xi = 0; xi < cx; xi++) {
                    id.data[p++] = color[0]; // Red component
                    id.data[p++] = color[1]; // Green component
                    id.data[p++] = color[2]; // Blue component
                    id.data[p++] = 255;      // Alpha (fully opaque)
                }
            }
        }
    }

    // Render updated image to canvas
    c.putImageData(id, 0, 0);

    // Draw streamlines if enabled
    if (scene.showStreamlines) {
        c.strokeStyle = "#0000FF"; // Blue color for streamlines
        for (var i = 1; i < f.numX; i += 13) {
            for (var j = 1; j < f.numY - 1; j += 5) {
                var x = (i + 0.5) * f.cellSize;
                var y = (j + 0.5) * f.cellSize;
                c.beginPath();
                c.moveTo(cX(x), cY(y));
                for (var n = 0; n < 15; n++) {
                    var u = f.fieldSample(x, y, U_FIELD); // Sample velocity field in X direction
                    var v = f.fieldSample(x, y, V_FIELD); // Sample velocity field in Y direction
                    l = Math.sqrt(u * u + v * v); // Compute velocity magnitude
                    x += u * 0.02 / (1 + l); // Adjust position based on velocity
                    y += v * 0.02 / (1 + l);
                    if (x > f.numX * f.cellSize) break;
                    c.lineTo(cX(x), cY(y));
                }
                c.stroke();
                c.lineWidth = 1;
            }
        }
    }
    //scene.numIters*f.numCells*scene.frameNr
    //Display the average diveegence of the whole grid 
    var d =`Avegrage Divergence: ${(f.totalDivergence/(f.divFreq)).toFixed(9)}`

    // Display min/max pressure on screen if 
    if (scene.showPressure) {
        var s = `Pressure (min)/(mid)/(max): ${minP.toFixed(0)}, ${((minP+ maxP)/2).toFixed(0)}, ${maxP.toFixed(0)} N/m`;
        if(scene.paused){
            s = `Pressure at point (${scene.obstacleX.toFixed(4)},${scene.obstacleY.toFixed(4)}): ${(scene.totalPressure/scene.count).toFixed(0)} N/m`;
        }
        
        d =`Avegrage Divergence: ${(f.totalDivergence/(f.divFreq)).toFixed(9)}` 

        if(scene.showSmoke){
        c.fillStyle = "#FFFFFF";
        }else{
        c.fillStyle = "#000000";
        }
        c.font = "20px Arial bold";
        c.fillText(s, 10, 30); // Draw text at top-left corner of canvas

    }
    c.font = "20px Arial bold";
    c.fillText(d, 1100, 30)
}

function setObstacle(x, y, reset) {

    // Initialize velocity components for the obstacle
    var vx = 0.0; // Horizontal velocity
    var vy = 0.0; // Vertical velocit
    // If the obstacle is not being reset, calculate its velocity based on movement
    if (!reset) {
        vx = (x - scene.obstacleX) / scene.dt; // Compute horizontal velocity
        vy = (y - scene.obstacleY) / scene.dt; // Compute vertical velocity
    }
    
    // Update obstacle position in the scene
    scene.obstacleX = x;
    scene.obstacleY = y;
    // Reference the fluid simulation object from the scene
    var f = scene.fluid;
    var n = f.numY; // Number of cells in the Y direction
    
    // Retrieve the obstacle radius from the user interface (slider value)
    scene.obstacleRadius = document.getElementById("shapeSize").value / 100;
    
    // Retrieve the slider element controlling the obstacle size
    var slider = document.getElementById("shapeSize");
    
    // Update the displayed value when the slider is adjusted
    slider.oninput = function() {
        var output = document.getElementById("sliderShapeSize");
        output.innerHTML = slider.value;
    };

    if (scene.square) {  
        var squareWidth = scene.obstacleRadius; // Set the square's width based on user input
    
        // Iterate through the entire fluid grid, excluding the boundary cells
        for (var i = 1; i < f.numX - 2; i++) {
            for (var j = 1; j < f.numY - 2; j++) {
                
                f.cellState[i * n + j] = 1.0; // Mark all cells as fluid by default
    
                // Calculate distance of the current cell from the obstacle's center
                dx = (i + 0.5) * f.cellSize - x;
                dy = (j + 0.5) * f.cellSize - y;
    
                // Check if the current cell lies within the square obstacle's boundaries
                if (dx > -squareWidth / 2 && dx < squareWidth / 2 &&
                    dy > -squareWidth / 2 && dy < squareWidth / 2) {
    
                    f.cellState[i * n + j] = 0.0; // Mark the cell as part of the solid obstacle
    
                    // Assign material properties to the obstacle if it's in scene 2, else set default properties
                    if (scene.sceneNr == 2) 
                        f.m[i * n + j] = 0.5 + 0.5 * Math.sin(0.1 * scene.frameNr);
                    else 
                        f.m[i * n + j] = 1.0;
    
                    // Apply velocity to the obstacle cells to match obstacle movement
                    f.horizontalVelocity[i * n + j] = vx;
                    f.horizontalVelocity[(i + 1) * n + j] = vx;
                    f.verticalVelocity[i * n + j] = vy;
                    f.verticalVelocity[i * n + j + 1] = vy;
                }
            }
        }   
    }
    
    if (scene.rainDrop) {  
        var r = scene.obstacleRadius; // Set the radius of the raindrop based on user input
    
        // Iterate through the entire fluid grid, excluding the boundary cells
        for (var i = 1; i < f.numX - 2; i++) {
            for (var j = 1; j < f.numY - 2; j++) {
    
                f.cellState[i * n + j] = 1.0; // Mark all cells as fluid by default
    
                // Calculate the distance of the current cell from the obstacle’s center
                dx = (i + 0.5) * f.cellSize - x;
                dy = (j + 0.5) * f.cellSize - y;    
                // Check if the current cell falls within the raindrop-shaped obstacle
                // The equation here represents an approximation of a falling water droplet shape
                if (dx * dx + (dy / -((dx - r) / (2 * r))) * (dy / -((dx - r) / (2 * r))) < r * r) { 
                    f.cellState[i * n + j] = 0.0; // Mark the cell as part of the solid obstacle
    
                    // Assign material properties if in scene 2, allowing dynamic variations
                    if (scene.sceneNr == 2) 
                        f.m[i * n + j] = 0.5 + 0.5 * Math.sin(0.1 * scene.frameNr);
                    else 
                        f.m[i * n + j] = 1.0;
    
                    // Apply velocity to obstacle cells to match obstacle movement
                    f.horizontalVelocity[i * n + j] = vx;
                    f.horizontalVelocity[(i + 1) * n + j] = vx;
                    f.verticalVelocity[i * n + j] = vy;
                    f.verticalVelocity[i * n + j + 1] = vy;
                }
            }
        }
    }
    if (scene.aeroFoil) {  
        // Get the user-defined angle of attack (AoA) for the airfoil and convert it from degrees to radians
        var angle = document.getElementById("angle").value * -(Math.PI / 180);
        var slider = document.getElementById("angle");
    
        // Update the displayed angle value dynamically when the slider is adjusted
        slider.oninput = function(){
            var output = document.getElementById("sliderAngle");
            output.innerHTML = slider.value;
        }
    
        // Get the user-defined curvature parameter (K) for the airfoil camber line
        var k = document.getElementById("K").value / 100;
        var aeroCurveSlider = document.getElementById("K");
    
        // Update the displayed curvature value dynamically when the slider is adjusted
        aeroCurveSlider.oninput = function(){
            var output = document.getElementById("k");
            output.innerHTML = aeroCurveSlider.value;
        }
    
        // Define the radius scale of the airfoil shape
        var r = scene.obstacleRadius * 4;
    
        // Loop through all grid cells in the fluid domain, excluding boundary cells
        for (var i = 1; i < f.numX - 2; i++) {
            for (var j = 1; j < f.numY - 2; j++) {
    
                f.cellState[i * n + j] = 1.0; // Default state is fluid (1.0)
    
                // Apply a linear transformation (rotation) to align the airfoil at the given angle
                // The transformation rotates the grid coordinates (i, j) to match the airfoil orientation
                rotatedI = i * Math.cos(angle) + j * Math.sin(angle);
                rotatedJ = i * Math.sin(angle) - j * Math.cos(angle);
    
                rotatedX = x * Math.cos(angle) + y * Math.sin(angle);
                rotatedY = x * Math.sin(angle) - y * Math.cos(angle);
    
                // Compute distances from the rotated airfoil position
                dx = (rotatedI + 0.5) * f.cellSize - rotatedX;
                dy = (rotatedJ + 0.5) * f.cellSize - rotatedY;
                
                //The equation below defines the shape of the cambered airfoil profile.
                //It is derived from thin airfoil theory, approximating the camber line.
                if (dy < -k * dx * ((1 - (1 / r) * dx) * (1 - (1 / r) * dx)) + 
                          0.15 * Math.sqrt(r * dx) * Math.sqrt(2 - (1 / r) * dx) * (1 - (1 / r) * dx) && 
                    dy > -k * dx * ((1 - (1 / r) * dx) * (1 - (1 / r) * dx)) - 
                          0.15 * Math.sqrt(r * dx) * Math.sqrt(2 - (1 / r) * dx) * (1 - (1 / r) * dx)) {
    
                    f.cellState[i * n + j] = 0.0; // Mark the cell as part of the airfoil (solid obstacle)
    
                    // Assign material properties based on the scene configuration
                    if (scene.sceneNr == 2) 
                        f.m[i * n + j] = 0.5 + 0.5 * Math.sin(0.1 * scene.frameNr); // Dynamic material property
                    else 
                        f.m[i * n + j] = 1.0; // Default material property
    
                    // Assign velocity components to the obstacle cells to simulate movement
                    f.horizontalVelocity[i * n + j] = vx;
                    f.horizontalVelocity[(i + 1) * n + j] = vx;
                    f.verticalVelocity[i * n + j] = vy;
                    f.verticalVelocity[i * n + j + 1] = vy;
                }
            }
        }
    }
    if (scene.circle) {  // Check if the circular obstacle mode is enabled

        // Create an array to track previous solid obstacle states
        var previousSolid = [f.numCells]; 
        previousSolid.fill(0);  // Initialize all cells as not being solid before

    
        // Define the radius of the circular obstacle, scaled down by half
        var r = scene.obstacleRadius / 2;
    
        // Loop through all fluid grid cells, excluding boundary cells
        for (var i = 1; i < f.numX - 2; i++) {
            for (var j = 1; j < f.numY - 2; j++) {
    
                // If drawing is disabled, reset the cell state to fluid (1.0)
                if (!scene.paused){
                    
                    if (!scene.draw) {
                        f.cellState[i * n + j] = 1.0;
                    } 
                    // If drawing is enabled but this cell was not previously solid, set it as fluid
                    else if (previousSolid[i * n + j] == 0) {
                        f.cellState[i * n + j] = 1.0;
                    }
                }
                // Compute distances from the center of the circular obstacle
                dx = (i + 0.5) * f.cellSize - x;
                dy = (j + 0.5) * f.cellSize - y;
    
                // Check if the current cell is inside the circular region
                // Equation: (dx^2 + dy^2 < r^2) represents a filled circle
                if (dx * dx + dy * dy < r * r) {
                    
                    if(!scene.paused){
                    f.cellState[i * n + j] = 0.0; // Mark the cell as a solid obstacle
                    }else{
                        scene.count =+ 1;
                        scene.totalPressure =+ f.p[i * n + j];
                    }

                    if (scene.draw){ 
                        f.cellState[i * n + j] = 0.0; 
                        previousSolid[i * n + j] = 1; // Store as previously solid
                    }
                    // Assign material properties based on the scene configuration
                    if (scene.sceneNr == 2) 
                        f.m[i * n + j] = 0.5 + 0.5 * Math.sin(0.1 * scene.frameNr); // Dynamic material property
                    else 
                        f.m[i * n + j] = 1.0; // Default material property
                    // Assign velocity components to the obstacle cells to simulate movement
                    f.horizontalVelocity[i * n + j] = vx;
                    f.horizontalVelocity[(i + 1) * n + j] = vx;
                    f.verticalVelocity[i * n + j] = vy;
                    f.verticalVelocity[i * n + j + 1] = vy;
                }
            }   
        }
    }
    // Enable obstacle visibility in the scene
    scene.showObstacle = true;
}

function resetDiv(){
    var f = scene.fluid;
    f.totalDivergence = 0; // Resets total divegence to 0
    f.divFreq = 0; // Resets the count to 0
}


// interaction -------------------------------------------------------

var mouseDown = false;

function startDrag(x, y) {
    let bounds = canvas.getBoundingClientRect();

    let mx = x - bounds.left - canvas.clientLeft + 10;
    let my = y - bounds.top - canvas.clientTop;
    mouseDown = true;

    x = mx / cScale;
    y = (canvas.height - my) / cScale;

    setObstacle(x,y, true);
}

function drag(x, y) {
    if (mouseDown) {
        let bounds = canvas.getBoundingClientRect();
        let mx = x - bounds.left - canvas.clientLeft + 10;
        let my = y - bounds.top - canvas.clientTop;
        x = mx / cScale;
        y = (canvas.height - my) / cScale;
        setObstacle(x,y, false);
        globalX = x;
        globalY = y;
    }
}

function endDrag() {
    mouseDown = false;
}

function checkRes(res){
    if (res > 150){
        alert("Warning!\nThe recommended resolution is less than 150.\nPress Wind Tunnel if you wish to continue.");
    }else if (res > 0){
        setupScene(1)
    }else{
       alert("Invald input. Please try again.")
    }
}

function controlBox(){
    alert("Controls:\nP - Pause the simulation\nM - Run one iteration\nA - Toggle Pressure\nS - Toggle Smoke\nD - Toggle Draw Mode\nWIND TUNNEL - If anything goes wrong\nMouse - For everything else");
}


canvas.addEventListener('mousedown', event => {
    startDrag(event.x, event.y);
});

canvas.addEventListener('mouseup', event => {
    endDrag();
});

canvas.addEventListener('mousemove', event => {
    drag(event.x, event.y);
});

canvas.addEventListener('touchstart', event => {
    startDrag(event.touches[0].clientX, event.touches[0].clientY)
});

canvas.addEventListener('touchend', event => {
    endDrag()
});

canvas.addEventListener('touchmove', event => {
    event.preventDefault();
    event.stopImmediatePropagation();
    drag(event.touches[0].clientX, event.touches[0].clientY)
}, { passive: false});

document.addEventListener('keydown', event => {
    switch(event.key) {
        case 'p': scene.paused = !scene.paused; 
        document.getElementById("sliderShapeSize").value = "1";
        document.getElementById('shapeSize').value = "1";break;
        case 'm': scene.paused = false; run(); scene.paused = true; break;
        case 's': document.getElementById('smokeButton').click(); break;
        case 'd': document.getElementById('drawButton').click(); break
        case 'a': document.getElementById('pressureButton').click(); break;
    }
});



// main -------------------------------------------------------

// Function to run the fluid simulation
function run() {
    // Check if the scene is not paused before running the simulation
    if (!scene.paused)  
        // Perform fluid simulation with time step, gravity, and iteration count
        scene.fluid.simulate(scene.dt, scene.gravity, scene.numIters);
    
    // Increment the frame counter
    scene.frameNr++;
}

// Function to update the scene continuously
function update() {
    run();  // Run the simulation step
    drawCanvas();      // Render the updated scene
    requestAnimationFrame(update); // Schedule the next frame update
}

// Initialize the scene with a given parameter (e.g., level of detail or mode)
setupScene(1);

// Start the update loop
update();
