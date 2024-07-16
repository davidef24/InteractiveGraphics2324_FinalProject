
class Cell {
    constructor() {
        this.max_particles = 100;
        this.particles = new Array(this.max_particles);
        this.neighbors = [];
        this.numParticles = 0;
    }
}
export default Cell;