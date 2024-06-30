// Constants
import * as THREE from './three.module.js'
import { GUI } from './dat.gui.module.js'

// Constants
const boxWidth = 3;
const boxDepth = 2;
const boxHeight = 2;

class Particle {
    constructor(x, y, z, width, height, depth) {
        this.position = new THREE.Vector3(x, y, z);
        this.velocity = new THREE.Vector3(0, 0, 0);
        this.acceleration = new THREE.Vector3(0.01, 0.01, 0.01);
        this.density = 0;
        this.pressure = 0;

        this.width = width;
        this.height = height;
        this.depth = depth;
    }
}

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);

const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setAnimationLoop(animate);
document.body.appendChild(renderer.domElement);

const geometry = new THREE.BoxGeometry(boxWidth, boxHeight, boxDepth);

// Adjust material properties for transparency and wireframe
const material = new THREE.MeshBasicMaterial({
    color: 0xFAEBD7,
    transparent: true,
    wireframe: true,
    wireframeLinewidth: 2,
});

const cube = new THREE.Mesh(geometry, material);
scene.add(cube);

const particleGeometry = new THREE.SphereGeometry(0.01, 16, 8);
const particleMaterial = new THREE.MeshBasicMaterial({ color: 0x1E90FF, transparent: true, wireframe: true });

let particles = [];
let particleMeshes = [];

// Function to create particles
function createParticles(count) {
    // Remove existing particles
    particleMeshes.forEach(mesh => scene.remove(mesh));
    particles = [];
    particleMeshes = [];

    // Create new particles
    for (let i = 0; i < count; i++) {
        const x = (Math.random() - 0.5) * boxWidth;
        const y = 0.5 * boxHeight;
        const z = (Math.random() - 0.5) * boxDepth;

        const particle = new Particle(x, y, z, boxWidth, boxHeight, boxDepth);
        particles.push(particle);

        const particleMesh = new THREE.Mesh(particleGeometry, particleMaterial);
        particleMesh.position.copy(particle.position);
        particle.mesh = particleMesh;
        scene.add(particleMesh);
        particleMeshes.push(particleMesh);
    }
}

// Create initial particles
createParticles(2000);

camera.position.z = 5;

// GUI setup
const gui = new GUI();
const guiOptions = {
    particleCount: 2000,
};

gui.add(guiOptions, 'particleCount', 0, 2000).step(1).onChange(value => {
    createParticles(value);
});

let isRotating = false;
let previousMousePosition = { x: 0, y: 0 };

function handleMouseDown(event) {
    isRotating = true;
    previousMousePosition = { x: event.clientX, y: event.clientY };
}

function handleResize(event) {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}

function handleMouseMove(event) {
    if (isRotating) {
        let deltaMove = { x: event.clientX - previousMousePosition.x, y: event.clientY - previousMousePosition.y };

        cube.rotation.y += deltaMove.x * 0.01;
        cube.rotation.x += deltaMove.y * 0.01;

        previousMousePosition = { x: event.clientX, y: event.clientY };
    }
}

const keyState = {
    ArrowUp: false,
    ArrowDown: false
};

document.addEventListener('keydown', (event) => {
    if (event.key === 'ArrowUp' || event.key === 'ArrowDown') {
        keyState[event.key] = true;
    }
});

document.addEventListener('keyup', (event) => {
    if (event.key === 'ArrowUp' || event.key === 'ArrowDown') {
        keyState[event.key] = false;
    }
});

function handleMouseWheel(event) {
    let delta = event.deltaY * 0.001;
    camera.position.z += delta;
}

function handleMouseUp(event) {
    isRotating = false;
}

//mouse event listeners
renderer.domElement.addEventListener('mousedown', handleMouseDown, false);
renderer.domElement.addEventListener('mousemove', handleMouseMove, false);
renderer.domElement.addEventListener('mouseup', handleMouseUp, false);
renderer.domElement.addEventListener('wheel', handleMouseWheel, false);
window.addEventListener('resize', handleResize);

function animate() {
    requestAnimationFrame(animate);

    if (keyState.ArrowUp) {
        camera.position.y += 0.002;
    }
    if (keyState.ArrowDown) {
        camera.position.y -= 0.002;
    }
    particleMeshes.forEach((mesh, index) => {
        mesh.position.copy(particles[index].position);
    });

    renderer.render(scene, camera);
}

animate();
