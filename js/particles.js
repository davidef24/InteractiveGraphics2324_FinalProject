import * as THREE from './three.module.js'

class Particle {
    constructor(x, y, z, width, height, depth) {
        this.position = new THREE.Vector3(x, y, z);
        this.x = x;
        this.y = y;
        this.z = z;
        this.Vx = 0;
        this.Vy = 0;
        this.Vz = 0;
        this.Fx = 0;
        this.Fy = 0;
        this.Fz = 0;
        this.rho = 0; //density
        this.P = 0; //pressure

        this.width = width;
        this.height = height;
        this.depth = depth;
    }

    reset() {
        this.Fx = 0;
        this.Fy = 0;
        this.Fz = 0;
        this.rho = m * Wpoly6(0);
    }
}

export default Particle;