import Grid3D from './grid.js';
import Particle from './particles.js';
import * as THREE from './three.module.js'

let h = 1;     // smoothing length
const Wpoly6_coeff = (h) => 315.0 / (64 * Math.PI * Math.pow(h, 9));
const Wspiky_grad_coeff = (h) => -45.0 / (Math.PI * Math.pow(h, 6));
const Wvisc_lapl_coeff = (h) => 45.0 / (Math.PI * Math.pow(h, 5));
let m = 1.0;	    // Particle mass
let k = 120;				// Gas constant
let rho0 = 0;			// Rest density
let mu = 10;				// Viscosity
let gx = 0;				// Gravity-x
let gy = 0;			
let gz = 0;
//let clock= 0;

let velocityMultiplier = 1.5;

let numGridCellsX;
let numGridCellsY;
let numGridCellsZ; 

const domainScale = 0.04;

class Engine{

    //density
    kernel_wpoly6(r2) {
        let temp = Math.pow(h, 2) - r2;
        return Wpoly6_coeff(h) * temp * temp * temp;
    }
    
    //pressure
    // assumption: r is less than h
    kernel_wspiky_grad2(r) {
        let temp = h - r;
        return Wspiky_grad_coeff(h) * temp * temp / r;
    }
    
    //viscosity
    // assumption: r is less than h
    kernel_wvisc_lapl(r) {
        return Wvisc_lapl_coeff(h) * (1 - r / h);
    }
    
    dist2(p1, p2) {
        let dx = p2.x - p1.x;
        let dy = p2.y - p1.y;
        let dz = p2.z - p1.z;
        return dx * dx + dy * dy + dz * dz;
    }

    addParticleToCell(p) {
        let c = this.grid.getCellFromLocation(p.x, p.y, p.z);
        if (c !== null) {
            c.particles[c.numParticles++] = p;
        } else {
            console.log("Undefined grid cell [engine]!");
        }
    }

    setFluidProperties(mass, gasConstant, restDensity, viscosity, smoothingLength) {
        h = smoothingLength;
        m = mass;
        k = gasConstant;
        rho0 = restDensity;
        mu = viscosity;
    }

    setGravity(gr_x, gr_y, gr_z){
        gx = gr_x;
        gy = gr_y;
        gz = gr_z;
    }

    reset(){
        this.grid.setParticlesToZero();
        this.grid.clearParticleReferences();
    }

    updateParticleCount(n, boxDimensions, particleMeshes, particleMass) {
        //console.log("[ENGINE] particleMeshes received is ", particleMeshes);
        let i0 = 0;
        if (this.particles == undefined) {
            this.particles = new Array(n);
        }
        else {
            if (n < this.particles.length) {
                this.particles.length = n;
            }
            i0 = this.particles.length;
        }
        for (let i = i0; i < n; i++) {
            this.particles[i] = new Particle(particleMeshes[i], domainScale, particleMass);
            this.particles[i].rho = m * this.kernel_wpoly6(0);
            //console.log("[ENGINE] I-TH PARTICLE IS",i, this.particles[i]);
        }
        //console.log("[ENGINE] ALL PARTICLES ARE",this.particles);
    }
    

    constructor(width, height, depth, left, right, bottom, top, front, back){

        while (left >= width) {
            // assume two identical monitor setup arranged horizontally
            width += width;
        }
        this.xlimit = width / domainScale;
        this.xmin = left / domainScale;
        this.xmax = right / domainScale;

        this.ylimit = height / domainScale;
        this.ymin = bottom / domainScale;
        this.ymax = top / domainScale;

        this.zlimit = depth / domainScale;
        this.zmin = front / domainScale;
        this.zmax = back / domainScale;

        numGridCellsX = Math.floor(this.xlimit);
        if(numGridCellsX == 0) console.log("nx == 0");
        numGridCellsY = Math.floor(this.ylimit);
        if(numGridCellsY == 0) console.log("ny == 0");
        numGridCellsZ = Math.floor(this.zlimit);
        if(numGridCellsX == 0 || numGridCellsY == 0 || numGridCellsZ == 0) console.error("cells number cannot be zero");
        this.grid = new Grid3D(numGridCellsX, numGridCellsY, numGridCellsZ, this.xmin, this.xmax, this.ymin, this.ymax, this.zmin, this.zmax, domainScale);

        this.particles = []
        
        this.lastTime = performance.now();
        this.forceVelocityCell = null;  // cell where velocity should be forced when mouse is over it
        this.forceVx = 0;
        this.forceVy = 0;
        this.forceVz = 0;
    }
    
    computeDensity() {
        // Loop through all cells in the grid
        for (const cell of this.grid.cells) {
          // Loop through all particles in the current cell
          for (let i = 0; i < cell.numParticles; i++) {
            const particle1 = cell.particles[i];
      
            // Calculate density contribution from neighboring particles (within cell and neighbors)
            let densityContribution = 0;
            for (let j = i + 1; j < cell.numParticles; j++) {
              const particle2 = cell.particles[j];
              densityContribution += this.getDensityContribution(particle1, particle2);
            }
            for (const neighbor of cell.neighbors) {
              for (let j = 0; j < neighbor.numParticles; j++) {
                const particle2 = neighbor.particles[j];
                densityContribution += this.getDensityContribution(particle1, particle2);
              }
            }
      
            // Update particle density and pressure
            particle1.rho += densityContribution;
            particle1.P = Math.max(k * (particle1.rho - rho0), 0);
          }
        }
      }
      
      getDensityContribution(particle1, particle2) {
        const r2 = this.dist2(particle1, particle2);
        if (r2 < Math.pow(h, 2)) {
          return m * this.kernel_wpoly6(r2);
        } else {
          return 0;
        }
      }
      
    
    
      updateParticleForces(p1, p2) {
        let r2 = this.dist2(p1, p2);
        
        if (r2 < Math.pow(h, 2)) {
            let r = Math.sqrt(r2) + 1e-6; // Add a tiny bit to avoid divide by zero
    
            // Compute the common term for pressure force
            let avgPressure = (p1.P + p2.P) / 2;
            let pressureFactor = m * avgPressure / p2.rho;
            let pressureGradient = this.kernel_wspiky_grad2(r);
            let temp1 = pressureFactor * pressureGradient;
    
            // Compute pressure force components
            let fx_pressure = temp1 * (p2.x - p1.x);
            let fy_pressure = temp1 * (p2.y - p1.y);
            let fz_pressure = temp1 * (p2.z - p1.z);
    
            // Compute the common term for viscosity force
            let viscosityFactor = mu * m * this.kernel_wvisc_lapl(r) / p2.rho;
    
            // Compute viscosity force components
            let fx_viscosity = viscosityFactor * (p2.Vx - p1.Vx);
            let fy_viscosity = viscosityFactor * (p2.Vy - p1.Vy);
            let fz_viscosity = viscosityFactor * (p2.Vz - p1.Vz);
    
            // Total forces
            let fx_total = fx_pressure + fx_viscosity;
            let fy_total = fy_pressure + fy_viscosity;
            let fz_total = fz_pressure + fz_viscosity;
    
            // Apply forces to particles
            p1.Fx += fx_total;
            p1.Fy += fy_total;
            p1.Fz += fz_total;
    
            p2.Fx -= fx_total;
            p2.Fy -= fy_total;
            p2.Fz -= fz_total;
        }
    }
    
    
    addWallForces(p1) {
        // Define the kernel weight function
        const kernelWeight = (r) => m * p1.P / p1.rho * this.kernel_wspiky_grad2(r) * r;
    
        // Check and apply forces for the x boundaries
        if (p1.x < this.xmin + h) {
            let r = p1.x - this.xmin;
            p1.Fx -= kernelWeight(r);
        } else if (p1.x > this.xmax - h) {
            let r = this.xmax - p1.x;
            p1.Fx += kernelWeight(r);
        }
    
        // Check and apply forces for the y boundaries
        if (p1.y < this.ymin + h) {
            let r = p1.y - this.ymin;
            p1.Fy -= kernelWeight(r);
        } else if (p1.y > this.ymax - h) {
            let r = this.ymax - p1.y;
            p1.Fy += kernelWeight(r);
        }
    
        // Check and apply forces for the z boundaries
        if (p1.z < this.zmin + h) {
            let r = p1.z - this.zmin;
            p1.Fz -= kernelWeight(r);
        } else if (p1.z > this.zmax - h) {
            let r = this.zmax - p1.z;
            p1.Fz += kernelWeight(r);
        }
    }
    
    
    computeForces() {
        for (const cell of this.grid.cells) {
            for (let i = 0; i < cell.numParticles; i++) {
                const p1 = cell.particles[i];
                for (let j = i + 1; j < cell.numParticles; j++) {
                    const p2 = cell.particles[j];
                    this.updateParticleForces(p1, p2);
                }
                for (let neighbor of cell.neighbors) {
                    for (let j = 0; j < neighbor.numParticles; j++) {
                        const p2 = neighbor.particles[j];
                        this.updateParticleForces(p1, p2);
                    }
                }
                this.addWallForces(p1);
            }
        }
    }

    updateMeshPosition(i, position) {
        //console.log("[GET PARTICLE POSITION] PARTICLES ARE ", this.particles);
        let p = this.particles[i];
        position.x = p.x * domainScale;
        position.y = p.y * domainScale; 
        position.z = p.z * domainScale; 
    }

    removeParticleFromCell(cell, particle) {
        let cell_particles = cell.particles;
        let index = cell_particles.indexOf(particle);
    
        if (index !== -1) {
            cell_particles.splice(index, 1);  // Remove the particle at the found index
            cell.numParticles--;  // Decrement the particle count
        } else {
            console.log("Particle not found in the specified cell.");
        }
    }
    
    updatePosition(dT) {
        for (let p of this.particles) {
            let idx = this.grid.getCellIndexFromLocation(p.x, p.y, p.z);
            let c= this.grid.cells[idx];
            this.removeParticleFromCell(c,p);

            //console.log("[UPDATE POSITION] forces are ", p.Fx, p.Fy, p.Fz);
            //console.log("[UPDATE POSITION] rho is ", p.rho);
            let Ax = p.Fx / p.rho + gx;
            let Ay = p.Fy / p.rho + gy;
            let Az = p.Fz / p.rho + gz; 

            //console.log("[UPDATE POSITION] accelerations are ", Ax, Ay, Az);
    
            p.Vx += Ax * dT;
            p.Vy += Ay * dT;
            p.Vz += Az * dT; 

            //console.log("[UPDATE POSITION] velocities are ", p.Vx, p.Vy, p.Vz);

            //console.log("[UPDATE POSITION] gonna change position of particle from ", p.x, p.y, p.z);
    
            p.x += (p.Vx + 0.5 * Ax * dT) * dT;
            p.y += (p.Vy + 0.5 * Ay * dT) * dT;
            p.z += (p.Vz + 0.5 * Az * dT) * dT;  

            //console.log("[UPDATE POSITION] changed position of particle to ", p.x, p.y, p.z);

    
            if (p.x < this.xmin) {
                p.x = this.xmin + 1e-6;
                p.Vx *= -0.8;
            } else if (p.x > this.xmax) {
                p.x = this.xmax - 1e-6;
                p.Vx *= -0.8;
            }
    
            if (p.y < this.ymin) {
                p.y = this.ymin + 1e-6;
                p.Vy *= -0.8;
            } else if (p.y > this.ymax) {
                p.y = this.ymax - 1e-6;
                p.Vy *= -0.8;
            }
    
            if (p.z < this.zmin) {
                p.z = this.zmin + 1e-6;
                p.Vz *= -0.8;
            } else if (p.z > this.zmax) {
                p.z = this.zmax - 1e-6;
                p.Vz *= -0.8;
            }
    
            this.grid.insertParticleIntoCell(p);

            let kW = this.kernel_wpoly6(0);
    
            p.initForceDensity(kW);
        }
    }


    simulationStep() {
        this.computeDensity();
        this.computeForces();
        this.updatePosition(1.0/60.0);
        //clock++;
    }

    
    forceVelocity(int_Objects, vx, vy) {
        let sum = 0;
        for (let obj of int_Objects) {
            if (obj.object.geometry instanceof THREE.SphereGeometry) {
                let mesh = obj.object;
                this.forceVelocityCell = this.grid.getCellFromLocation2(mesh.position.x, mesh.position.y, mesh.position.z);
                sum += this.forceVelocityCell.particles.length;
                this.forceVx = vx;
                this.forceVy = -vy;
                this.forceVz = 0;  //just force velocity on x and y axes
                for (let i = 0; i < this.forceVelocityCell.numParticles; i++) {
                    let p = this.forceVelocityCell.particles[i];
                    p.Vx = vx * velocityMultiplier;
                    p.Vy = vy * velocityMultiplier;
                    //p.Fx = 0;
                    //p.Fy = 0;
                    //p.Fz = 0;
                    //force velocity even to all neighbours
                    for (let neighbor of this.forceVelocityCell.neighbors) {
                        for (let j = 0; j < neighbor.numParticles; j++) {
                            const p2 = neighbor.particles[j];
                            p2.Vx = vx * velocityMultiplier;
                            p2.Vy = vy * velocityMultiplier;
                            //p2.Fx = 0;
                            //p2.Fy = 0;
                            //p2.Fz = 0;
                        }
                    }

                }
            }
        }
    
    }
}

export default Engine;