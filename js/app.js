import * as THREE from './three.module.js'
import { GUI } from './dat.gui.module.js'
import Engine from './engine.js'

//window.addEventListener('load', init, false);

// Constants
const boxWidth = 0.80;
const boxDepth = 0.35;
const boxHeight = 1;
// GUI setup
const gui = new GUI();
const initial_particles_n = 1000;
let initial_particle_radius = 0.005;
let initial_mass = 1.0;	    
let k = 120; //gas constant		
let rho0 = 0; //rest density
let mu = 3; //viscosity		
let gx = 0;	//gravity		
let gy = -10;			
let gz = 0;
let h = 1; //smoothing length
//particles
let particleMeshes = [];
const panelOptions = {
    ShowWireFrame: true,
    ParticleCount: initial_particles_n,
    ParticleRadius: initial_particle_radius,
    Mass: initial_mass,
    RestDensity: rho0,
    Viscosity: mu,
    GasConstant: k,
    SmoothingLength: h,
    GravityX: gx,
    GravityY: gy,
    GravityZ: gz,
    startSim: function() {
        startSimulation();
    }
};
//mouse events handling
let isRotating = false;
let previousMousePosition = { x: 0, y: 0 };




const halfWidth = boxWidth/2;
const halfHeight = boxHeight/2;
const halfDepth = boxDepth/2;

let engine = new Engine(boxWidth, boxHeight, boxDepth, -halfWidth, halfWidth, -halfHeight, halfHeight, -halfDepth, halfDepth);

//scene constants
const renderer = new THREE.WebGLRenderer();
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
const geometry = new THREE.BoxGeometry(boxWidth, boxHeight, boxDepth);
const cubeMaterial = new THREE.MeshBasicMaterial({
    color: 0x000000,
    transparent: true,
    wireframe: true,
    wireframeLinewidth: 1,
});

const particleMaterial = new THREE.MeshPhongMaterial({ color: 0x1E90FF, transparent: true });
const cube = new THREE.Mesh(geometry, cubeMaterial);

// Lights
const ambientLight = new THREE.AmbientLight(0x404040, 1.5); // soft white light
scene.add(ambientLight);

const pointLight = new THREE.PointLight(0xffffff, 1, 100);
pointLight.position.set(1, 2, 2);
scene.add(pointLight);



function updateParticleCount(n) {
    particleMeshes.forEach(mesh => scene.remove(mesh));
    // Reset the engine if the new count is less than the current count
    if (n < particleMeshes.length) {
        engine.reset();
    }
    particleMeshes = new Array(n);

    let particleGeometry = new THREE.SphereGeometry(panelOptions.ParticleRadius, 16, 8);
    for (var i = 0; i < n; i++) {
        var m = new THREE.Mesh(particleGeometry, particleMaterial);
        m.position.x = (Math.random() - 0.5) * boxWidth;
        m.position.y = (Math.random() * 0.5) * boxHeight; 
        m.position.z = (Math.random() - 0.5) * boxDepth;
        //console.log(" [SET NUM PARTICLES APP] Particle MESH after definition ", m);
        particleMeshes[i] = m;
        //console.log("[SET NUM PARTICLES APP] Particle meshes after setting i-th array entry ", i, particleMeshes);
        scene.add(m);
    }
    //console.log("[SET NUM PARTICLES APP] Particle meshes after definition ", particleMeshes);
    engine.updateParticleCount(n, [boxWidth, boxHeight, boxDepth], particleMeshes, panelOptions.Mass);
    renderer.render(scene, camera);
}


function init(){
    updateFluidProperties();
    updateGravity();
    initScene();
    addListeners();
    updateParticleCount(panelOptions.ParticleCount);
    setGuiPanel();
}


function initScene(){
    camera.position.z = 1;
    renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.setClearColor(0xFAEBD7, 1);  // Set background color to antiquewhite
    document.body.appendChild(renderer.domElement); //renderer creates a canvas element that is used to draw and render the scene
    scene.add(cube);
}

function updateFluidProperties() {
    let h= panelOptions.SmoothingLength;
    let mass = panelOptions.Mass;
    let gasConstant = panelOptions.GasConstant;
    let restDensity = panelOptions.RestDensity;
    let viscosity = panelOptions.Viscosity;
    engine.setFluidProperties(mass, gasConstant, restDensity, viscosity, h);
}

function updateGravity() {
    let gr_x = panelOptions.GravityX;
    let gr_y = panelOptions.GravityY;
    let gr_z = panelOptions.GravityZ;
    engine.setGravity(gr_x, gr_y, gr_z);
}


function setGuiPanel(){
    gui.add(panelOptions, 'ParticleCount', 0, 10000).step(1).onChange(value => {
        updateParticleCount(panelOptions.ParticleCount);
    });
    
    gui.add(panelOptions, 'ParticleRadius', 0.001, 0.020).step(0.001).onChange(value => {
        panelOptions.particleRadius = value;
        updateParticleCount(panelOptions.ParticleCount);
    });
    gui.add(panelOptions, 'Mass', 1, 1.10).step(0.01).onChange(updateFluidProperties);
    gui.add(panelOptions, 'GasConstant', 1, 1000).onChange(updateFluidProperties);
    gui.add(panelOptions, 'RestDensity', 0, 5).step(0.5).onChange(updateFluidProperties);
    gui.add(panelOptions, 'Viscosity', 0, 11).onChange(updateFluidProperties);
    gui.add(panelOptions, 'SmoothingLength', 1, 1.2).step(0.001).onChange(updateFluidProperties);
    gui.add(panelOptions, 'GravityX', -100, 100).step(1).onChange(updateGravity);
    gui.add(panelOptions, 'GravityY', -100, 100).step(1).onChange(updateGravity);
    gui.add(panelOptions, 'GravityZ', -100, 100).step(1).onChange(updateGravity);

    gui.add(panelOptions, 'startSim').name('Start Simulation');
}

function getMousePosition(event, element){
    var boundingRect = element.getBoundingClientRect();
    return {
        //event.clientX provides the horizontal coordinate of the mouse cursor relative to the viewport.
        //boundingRect.left is the distance from the left edge of the viewport to the left edge of the element.
        x: event.clientX - boundingRect.left,
        y: event.clientY - boundingRect.top
    };
}

function mouseDownHandler(event) {
    isRotating = true;
    previousMousePosition = { x: event.clientX, y: event.clientY };
}

//to force velocity
function mapMouseToWorld(x, y, canvas, camera) {
    // get the normalized device coordinates (NDC) of the mouse, required for setFromCamera API of raycaster
    let ndcX = (x / canvas.clientWidth) * 2 - 1;
    let ndcY = - (y / canvas.clientHeight) * 2 + 1;  // "-"" sign is due to the fact that canvas y-coordinate increases downward while 
    //the NDC y-coordinate increases upward, so the conversion needs to invert the direction
    
    // Create a ray from the camera through the mouse position in NDC
    let raycaster = new THREE.Raycaster();
    let mouseVector = new THREE.Vector2(ndcX, ndcY);
    raycaster.setFromCamera(mouseVector, camera);

    const intersects = raycaster.intersectObjects( scene.children );
    if(intersects.length == 0) return null;
    else return intersects;
}



function mouseMoveHandler(event) {
    if (isRotating) {
        let deltaMove = { x: event.clientX - previousMousePosition.x, y: event.clientY - previousMousePosition.y };
        console.log(deltaMove);
        // rotate camera around the cube
        const rotationSpeed = 0.005;
        const spherical = new THREE.Spherical();  // a point spherical coordinates
        spherical.setFromVector3(camera.position);
        //adjust spherical coordinates, where theta is the horizontal angle and phi is the vertical angle
        spherical.theta -= deltaMove.x * rotationSpeed;  //azimuthal angle
        spherical.phi -= deltaMove.y * rotationSpeed;  //polar angle
        spherical.phi = Math.max(0.01, Math.min(Math.PI - 0.01, spherical.phi));
        //reconvert into cartesian coordinates
        camera.position.setFromSpherical(spherical);
        //ensure camera always looks at the center of the cube
        camera.lookAt(cube.position);
        previousMousePosition = { x: event.clientX, y: event.clientY };
        renderer.render(scene, camera);
    }

    let mousePos = getMousePosition(event, renderer.domElement);
    let int_Objects = mapMouseToWorld(mousePos.x, mousePos.y, renderer.domElement, camera);
    //console.log("[FORCE VELOCITY] intersected objects ", int_Objects);
    if(int_Objects != null) {
        engine.forceVelocity(int_Objects, event.movementX, event.movementY);
    }
}


function mouseWheelHandler(event) {
    let zoomFactor = 1.05; // Adjust this factor to control the zoom speed
    if (event.deltaY < 0) {
       camera.fov /= zoomFactor; // Zoom in
    } else {
        camera.fov *= zoomFactor; // Zoom out
    }
    camera.updateProjectionMatrix();
    renderer.render(scene, camera);
}
    
  

function mouseUpHandler(event) {
    isRotating = false;
}

function addListeners(){
    //mouse event listeners
    renderer.domElement.addEventListener('mousedown', mouseDownHandler, false);
    renderer.domElement.addEventListener('mousemove', mouseMoveHandler, false);
    renderer.domElement.addEventListener('mouseup', mouseUpHandler, false);
    renderer.domElement.addEventListener('wheel', mouseWheelHandler, false);
}

let clock = 0;

function startSimulation() {
    requestAnimationFrame(startSimulation);
    //if(clock < 7){
        engine.simulationStep();
        for (var i = 0; i < particleMeshes.length; i++){
            //console.log("[ANIMATE] particleMeshes before update is ", particleMeshes);
            engine.updateMeshPosition(i, particleMeshes[i].position);
            //console.log("[ANIMATE] particleMeshes after update is ", particleMeshes);
        }
        renderer.render(scene, camera);
    //}
    //clock++

}
init();