let h = 1;     // smoothing length
let h2 = Math.pow(h, 2);
let h9 = Math.pow(h, 9);
let Wpoly6_coeff = 315.0 / (64 * Math.PI * h9);

let initialVelocityMultiplier= 10;

class Particle {
    constructor(particleMesh, domainScale, mass) {
        //console.log("[PARTICLE] particleMesh received is ", particleMesh);
        this.x = particleMesh.position.x / domainScale;
        this.y = particleMesh.position.y /domainScale;
        this.z = particleMesh.position.z / domainScale;
        //console.log(`Location x = ${this.x} y = ${this.y} z = ${this.z}`);
        this.Vx = 0;
        this.Vy = 0;
        this.Vz = 0;
        this.Fx = 0;
        this.Fy = 0;
        this.Fz = 0;
        this.rho = 0; //density
        this.P = 0; //pressure    
        this.m = mass;    
    }

    

    Wpoly6(r2) {
        let temp = h2 - r2;
        return Wpoly6_coeff * temp * temp * temp;
    }

    initForceDensity(kernelWeight) {
        ({ Fx: this.Fx, Fy: this.Fy, Fz: this.Fz } = { Fx: 0, Fy: 0, Fz: 0 });
        this.rho = this.m * kernelWeight;
    }
}

export default Particle;