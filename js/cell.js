const max_particles = 20;
class Cell {
    constructor() {
        this.particles = new Array(max_particles);
        this.neighbors = [];
        this.numParticles = 0;
    }
}
export default Cell;