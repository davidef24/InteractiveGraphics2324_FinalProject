import * as THREE from './three.module.js'
import { GUI } from './dat.gui.module.js'

import Particle from './particles.js'
import Engine from './engine.js'

window.addEventListener('load', init, false);

// Constants
const boxWidth = 2;
const boxDepth = 1;
const boxHeight = 1;
// GUI setup
const gui = new GUI();
const initial_particles_n = 500;
let initial_particle_radius = 0.008
const guiOptions = {
    particleCount: initial_particles_n,
    particleRadius: initial_particle_radius,
};
//mouse events handling
let isRotating = false;
let previousMousePosition = { x: 0, y: 0 };

//particles
let particles = [];
let particleMeshes = [];

let engine = new Engine();

//scene constants
const renderer = new THREE.WebGLRenderer();
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
const geometry = new THREE.BoxGeometry(boxWidth, boxHeight, boxDepth);
const cubeMaterial = new THREE.MeshBasicMaterial({
    color: 0x000000,
    transparent: true,
    wireframe: true,
    wireframeLinewidth: 2,
});
const particleMaterial = new THREE.MeshBasicMaterial({ color: 0x1E90FF});
const cube = new THREE.Mesh(geometry, cubeMaterial);



function init(){
    initScene();
    addListeners();
    addGUI();
    animate();
}


function initScene(){
    camera.position.z = 2;
    renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.setClearColor(0xFAEBD7, 1);  // Set background color to antiquewhite
    renderer.setAnimationLoop(animate);
    document.body.appendChild(renderer.domElement);
    scene.add(cube);
    createParticles(guiOptions.particleCount);
}


// Function to create particles
function createParticles(count) {
    // Remove existing particles
    particleMeshes.forEach(mesh => scene.remove(mesh));
    particles = [];
    particleMeshes = [];

    let particleGeometry = new THREE.SphereGeometry(guiOptions.particleRadius, 8, 4);

    // Create new particles
    for (let i = 0; i < count; i++) {
        const x = (Math.random() - 0.5) * boxWidth;
        const y = 0.5 * boxHeight; //spawn on top of the box
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

function addGUI(){
    gui.add(guiOptions, 'particleCount', 0, 4000).step(1).onChange(value => {
        createParticles(value);
    });
    
    gui.add(guiOptions, 'particleRadius', 0.004, 0.016).step(0.001).onChange(value => {
        guiOptions.particleRadius = value;
        createParticles(guiOptions.particleCount);
    });
}

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

        // Rotate camera around the cube
        const rotationSpeed = 0.005;
        const spherical = new THREE.Spherical();

        spherical.setFromVector3(camera.position);

        //adjust spherical coordinates, where theta is the horizontal angle and phi is the vertical angle
        spherical.theta -= deltaMove.x * rotationSpeed;
        spherical.phi -= deltaMove.y * rotationSpeed;

        spherical.phi = Math.max(0.01, Math.min(Math.PI - 0.01, spherical.phi));

        //reconvert into cartesian coordinates
        camera.position.setFromSpherical(spherical);

        //ensure camera always looks at the center of the cube
        camera.lookAt(cube.position);

        previousMousePosition = { x: event.clientX, y: event.clientY };
    }
}

function handleMouseWheel(event) {
    let delta = event.deltaY * 0.001;
    camera.position.z += delta;
}

function handleMouseUp(event) {
    isRotating = false;
}

function addListeners(){
    //mouse event listeners
    renderer.domElement.addEventListener('mousedown', handleMouseDown, false);
    renderer.domElement.addEventListener('mousemove', handleMouseMove, false);
    renderer.domElement.addEventListener('mouseup', handleMouseUp, false);
    renderer.domElement.addEventListener('wheel', handleMouseWheel, false);
    window.addEventListener('resize', handleResize);
}

function animate() {
    requestAnimationFrame(animate);
    engine.doPhysics();
    particleMeshes.forEach((mesh, index) => {
        mesh.position.copy(particles[index].position);
    });

    renderer.render(scene, camera);
}