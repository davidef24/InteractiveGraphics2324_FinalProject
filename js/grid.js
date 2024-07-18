import Cell from './cell.js'

class Grid3D {
    constructor(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, domainScale, maxCellParticles) {
        this.rows = nx; // Number of cells in x direction
        this.columns = ny; // Number of cells in y direction
        this.layers = nz; // Number of cells in z direction
        //console.log(`Grid is nx= ${nx} ny= ${ny} nz= ${nz}`);
        this.xmin = xmin; 
        this.xmax = xmax;
        this.ymin = ymin;
        this.ymax = ymax;
        this.zmin = zmin;
        this.zmax = zmax;
        this.cells = new Array(nx * ny * nz);
        this.domainScale = domainScale;
        for(let i=0; i< nx *  ny * nz; i++){
            this.cells[i] = new Cell(maxCellParticles);
        }
        this.setNeighbors();
    }

    setParticlesToZero() {
        for (let c of this.cells) {
            c.numParticles = 0;
        }
    }

    setNeighbors() {
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.columns; j++) {
                for (let k = 0; k < this.layers; k++) {
                    let idx = i + j * this.rows + k * this.rows * this.columns;
                    let c = this.cells[idx];
                    this.setCellNeighbors(i, j, k, c);
                }
            }
        }
    }

    setCellNeighbors(i, j, k, cell) {        
        for (let di = -1; di <= 1; di++) {
            for (let dj = -1; dj <= 1; dj++) {
                for (let dk = -1; dk <= 1; dk++) {
                    // Skip the current cell itself
                    if (di === 0 && dj === 0 && dk === 0) continue;
                    
                    let ni = i + di;
                    let nj = j + dj;
                    let nk = k + dk;
                    
                    // Check if the neighbor is within bounds
                    if (ni >= 0 && ni < this.rows &&
                        nj >= 0 && nj < this.columns &&
                        nk >= 0 && nk < this.layers) {
                        
                        let neighborIdx = ni + nj * this.rows + nk * this.rows * this.columns;
                        cell.neighbors.push(this.cells[neighborIdx]);
                    }
                }
            }
        }
    }
    

    getCellIndexFromLocation(x, y, z) {
        // Normalize coordinates to range [0, 1]
        
        let normX = (x - this.xmin) / (this.xmax - this.xmin);
        let normY = (y - this.ymin) / (this.ymax - this.ymin);
        let normZ = (z - this.zmin) / (this.zmax - this.zmin);
        // Scale normalized coordinates to grid indices
        let i = Math.floor(normX * this.rows);
        let j = Math.floor(normY * this.columns);
        let k = Math.floor(normZ * this.layers);
        let index = i + j * this.rows + k * this.rows * this.columns;
        //console.log(`Location x = ${x} y = ${y} z = ${z}`);
        //console.log(`Cell xmin = ${this.xmin} ymin = ${this.ymin} zmin = ${this.zmin}`);
        //console.log(`Cell xmax = ${this.xmax} ymax = ${this.ymax} zmax = ${this.zmax}`);
        //console.log(`Cell normX = ${normX} normY = ${normY} normZ = ${normZ}`);
        //console.log(`Cell i = ${i} j = ${j} k = ${k}`);
        return index;
    }

    //here, positions are not from the grid coordinate system, so need to scale down the boundaries points
    // I mean that in the gridFromLocation function positions are all been modifed in the engine class that operates in positions all scaled by the
    // domainScale factor.
    // this function contains object given by the THREE.RayCaster().intersectObjects() that are in the scene domain
    getCellFromLocation2(x, y, z) {
        // Normalize coordinates to range [0, 1]
        let ds = this.domainScale;

        let normX = (x - (this.xmin * ds)) / ((this.xmax - this.xmin) * ds);
        let normY = (y - (this.ymin * ds)) / ((this.ymax - this.ymin) * ds);
        let normZ = (z - (this.zmin * ds)) / ((this.zmax - this.zmin) * ds);
        // Scale normalized coordinates to grid indices
        let i = Math.floor(normX * this.rows);
        let j = Math.floor(normY * this.columns);
        let k = Math.floor(normZ * this.layers);
        let index = i + j * this.rows + k * this.rows * this.columns;
        //console.log(`Location x = ${x} y = ${y} z = ${z}`);
        //console.log(`Cell xmin = ${this.xmin} ymin = ${this.ymin} zmin = ${this.zmin}`);
        //console.log(`Cell xmax = ${this.xmax} ymax = ${this.ymax} zmax = ${this.zmax}`);
        //console.log(`Cell normX = ${normX} normY = ${normY} normZ = ${normZ}`);
        //console.log(`Cell i = ${i} j = ${j} k = ${k}`);
        let c = this.cells[index];
        return c;
    }

    insertParticleIntoCell(particle) {
        const cellIndex = this.getCellIndexFromLocation(particle.x, particle.y, particle.z);
        const cell = this.cells[cellIndex];
      
        if (cell !== undefined && cell !== null && cell.numParticles < cell.max_particles) {
          cell.particles[cell.numParticles++] = particle;
        } else {
          // Try to find a neighboring cell with available space
          const nb = cell.neighbors
          for (const n of nb) {
            const neighborCell = n
            if (neighborCell !== undefined && neighborCell.numParticles < cell.max_particles) {
              neighborCell.particles[neighborCell.numParticles++] = particle;
              return; // Particle successfully inserted into neighbor cell
            }
          }
      
          // If no neighbor cell has space, handle the particle rejection
          console.error("[grid] No available space for particle!");
        }
      }
      

    // also clears references to particles
    clearParticleReferences() {
        this.cells.forEach(cell => {
            cell.particles = Array(this.maxParticles).fill(null);
        });
    }
}

export default Grid3D;