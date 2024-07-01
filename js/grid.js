class Grid3D {
    constructor(nx, ny, nz, w, h, d) {
        this.nx = nx; // Number of cells in x direction
        this.ny = ny; // Number of cells in y direction
        this.nz = nz; // Number of cells in z direction
        this.w = w;   // Domain width
        this.h = h;   // Domain height
        this.d = d;   // Domain depth

        this.cells = new Array(nx * ny * nz).fill().map(() => new Cell());

        this.init();
    }

    init() {
        for (let i = 0; i < this.nx; i++) {
            for (let j = 0; j < this.ny; j++) {
                for (let k = 0; k < this.nz; k++) {
                    let idx = i + j * this.nx + k * this.nx * this.ny;
                    let c = this.cells[idx];
                    this.computeNeighbors(i, j, k, c);
                }
            }
        }
    }

    //for the moment, compute just the face-adjacent neighbours
    computeNeighbors(i, j, k, cell) {
        let idx = i + j * this.nx + k * this.nx * this.ny;

        // Add neighbors in x direction
        if (i > 0) {
            cell.halfNeighbors.push(this.cells[idx - 1]);
        }
        if (i < this.nx - 1) {
            cell.halfNeighbors.push(this.cells[idx + 1]);
        }

        // Add neighbors in y direction
        if (j > 0) {
            cell.halfNeighbors.push(this.cells[idx - this.nx]);
        }
        if (j < this.ny - 1) {
            cell.halfNeighbors.push(this.cells[idx + this.nx]);
        }

        // Add neighbors in z direction
        if (k > 0) {
            cell.halfNeighbors.push(this.cells[idx - this.nx * this.ny]);
        }
        if (k < this.nz - 1) {
            cell.halfNeighbors.push(this.cells[idx + this.nx * this.ny]);
        }
    }

    getCellFromLocation(x, y, z) {
        let i = Math.floor(this.nx * x / this.w);
        let j = Math.floor(this.ny * y / this.h);
        let k = Math.floor(this.nz * z / this.d);
        return this.cells[i + j * this.nx + k * this.nx * this.ny];
    }

    addParticleToCell(particle) {
        let cell = this.getCellFromLocation(particle.x, particle.y, particle.z);
        if (cell !== null) {
            cell.particles[cell.numParticles++] = particle;
        } else {
            console.log("Undefined grid cell!");
        }
    }

    reset() {
        this.cells.forEach(cell => {
            cell.numParticles = 0;
        });
    }

    hardReset() {
        this.cells.forEach(cell => {
            cell.numParticles = 0;
            cell.particles = new Array(INIT_MAX_PARTICLES_IN_CELL);
        });
    }
}

class Cell {
    constructor() {
        this.particles = new Array(INIT_MAX_PARTICLES_IN_CELL);
        this.halfNeighbors = [];
        this.numParticles = 0;
    }
}

